"""
EGA Patch Analysis — add files to an existing analysis.

Usage:
    python ega_patch_analysis.py <analysis_id> <folder>

Example:
    python ega_patch_analysis.py 128330 eso05_batch1
"""

import requests
import sys
import os

# ── Credentials ────────────────────────────────────────────────────────────────

def load_credentials(path=os.path.expanduser("~/.ega_credentials")):
    creds = {}
    with open(path) as f:
        for line in f:
            line = line.strip().removeprefix("export ")
            if "=" in line and not line.startswith("#"):
                key, _, value = line.partition("=")
                creds[key.strip()] = value.strip().strip('"').strip("'")
    username = creds.get("EGA_BOX")
    password = creds.get("ASPERA_SCP_PASS")
    if not username or not password:
        sys.exit("ERROR: EGA_BOX or ASPERA_SCP_PASS not found in ~/.ega_credentials")
    return username, password

TOKEN_URL = "https://idp.ega-archive.org/realms/EGA/protocol/openid-connect/token"
API_BASE  = "https://submission.ega-archive.org/api"

def get_token(username, password):
    resp = requests.post(TOKEN_URL, data={
        "grant_type": "password", "client_id": "sp-api",
        "username": username, "password": password,
    })
    resp.raise_for_status()
    return resp.json()["access_token"]

def headers(token):
    return {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python ega_patch_analysis.py <analysis_id> <folder>")

    analysis_id = int(sys.argv[1])
    folder      = sys.argv[2]

    username, password = load_credentials()
    token = get_token(username, password)
    print("  ✓ Token obtained")

    # Fetch inbox files for this folder
    resp = requests.get(f"{API_BASE}/files", headers=headers(token),
                        params={"prefix": f"/{folder}"})
    resp.raise_for_status()
    all_files = resp.json()
    print(f"  Found {len(all_files)} total files for prefix /{folder}")

    # Filter: .bam and .bam.bai only (exclude .md5 and other files)
    filtered = [f for f in all_files
                if f["relative_path"].endswith(".bam")
                or f["relative_path"].endswith(".bam.bai")]
    print(f"  After filter (.bam / .bam.bai): {len(filtered)} files")

    if not filtered:
        # Show sample paths to help diagnose
        print(f"  Sample relative_paths: {[f['relative_path'] for f in all_files[:5]]}")
        sys.exit("ERROR: No matching files found.")

    # Fetch existing analysis to build full payload for PUT
    resp = requests.get(f"{API_BASE}/analyses/{analysis_id}", headers=headers(token))
    resp.raise_for_status()
    existing = resp.json()

    # Merge new files with any already present
    existing_ids = {f["provisional_id"] for f in existing.get("files", [])}
    new_ids = [f["provisional_id"] for f in filtered]
    all_ids = list(existing_ids | set(new_ids))

    payload = {
        "alias":               existing.get("alias") or existing["title"].lower().replace(" ", "_"),
        "title":               existing["title"],
        "description":         existing["description"],
        "analysis_type":       existing["analysis_type"],
        "genome_id":           existing["genome_id"],
        "platform":            existing["platform"],
        "experiment_types":    existing["experiment_types"],
        "study_accession_id":  existing["study_accession_id"],
        "sample_accession_ids": [s["accession_id"] for s in existing["samples"]],
        "files":               all_ids,
    }
    resp = requests.put(f"{API_BASE}/analyses/{analysis_id}",
                        headers=headers(token), json=payload)
    if not resp.ok:
        sys.exit(f"  ✗ Failed: {resp.status_code} {resp.text}")

    print(f"  ✓ Analysis {analysis_id} updated with {len(filtered)} files")

if __name__ == "__main__":
    main()
