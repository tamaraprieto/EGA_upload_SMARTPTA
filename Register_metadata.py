"""
EGA Bulk Analysis Submission Script
------------------------------------
Creates REFERENCE_ALIGNMENT analyses via the EGA Submitter Portal API.

Usage:
    python ega_bulk_analyses.py <folder>

Example:
    python ega_bulk_analyses.py eso02_batch2

The sample name is extracted from the folder name (everything before '_batch'),
and looked up in SAMPLE_MAP to find the corresponding EGAN accession.
"""

import requests
import sys
import os

# ── Credentials (loaded from ~/.ega_credentials) ──────────────────────────────

def load_credentials(path=os.path.expanduser("~/.ega_credentials")):
    creds = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            line = line.removeprefix("export ")
            if "=" in line and not line.startswith("#"):
                key, _, value = line.partition("=")
                creds[key.strip()] = value.strip().strip('"').strip("'")
    username = creds.get("EGA_BOX")
    password = creds.get("ASPERA_SCP_PASS")
    if not username or not password:
        sys.exit("ERROR: EGA_BOX or ASPERA_SCP_PASS not found in ~/.ega_credentials")
    return username, password

# ── Configuration ─────────────────────────────────────────────────────────────

SUBMISSION_ID    = "EGA50000001666"
STUDY_ACCESSION  = "EGAS50000001793"
DATASET_ACCESSION = "EGAD50000002573"
ANALYSIS_TYPE    = "REFERENCE ALIGNMENT"
GENOME_ID = 15  # GRCh38 (GCA_000001405.15). See /api/enums/genomes for other IDs.

# Sample name → EGAN accession
SAMPLE_MAP = {
    "eso01": "EGAN50000419612",
    "eso02": "EGAN50000419609",
    "eso03": "EGAN50000419610",
    "eso04": "EGAN50000419611",
    "eso05": "EGAN50000419613",
}

# ── API endpoints ──────────────────────────────────────────────────────────────

TOKEN_URL = "https://idp.ega-archive.org/realms/EGA/protocol/openid-connect/token"
API_BASE  = "https://submission.ega-archive.org/api"

# ── Authentication ─────────────────────────────────────────────────────────────

def get_access_token(username, password):
    print("Authenticating...")
    resp = requests.post(TOKEN_URL, data={
        "grant_type": "password",
        "client_id":  "sp-api",
        "username":   username,
        "password":   password,
    })
    resp.raise_for_status()
    print("  ✓ Token obtained")
    return resp.json()["access_token"]

def auth_headers(token):
    return {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}

# ── File listing ───────────────────────────────────────────────────────────────

def list_inbox_files(token, prefix):
    for status in ["inbox", "submitted", "available"]:
        resp = requests.get(
            f"{API_BASE}/files",
            headers=auth_headers(token),
            params={"status": status, "prefix": f"/{prefix}"},
        )
        resp.raise_for_status()
        files = resp.json()
        if files:
            print(f"  Found files with status={status!r}")
            return files
    print(f"  ⚠ No files found for prefix: /{prefix}")
    return []

# ── Analysis creation ──────────────────────────────────────────────────────────

def create_analysis(token, folder, sample_name, batch_name, sample_accession):
    print(f"\nProcessing folder: {folder}")

    files = list_inbox_files(token, folder)
    if not files:
        return None
    # Only include encrypted BAM and BAI files
    files = [f for f in files if f["relative_path"].endswith(".bam.gpg") or f["relative_path"].endswith(".bam.bai.gpg")]
    print(f"  Found {len(files)} file(s) (excluding .md5)")

    payload = {
        "alias":               folder,
        "title":               f"BAM files for {sample_name} ({batch_name})",
        "description":         "Bam files resulting from merging all cram files per cell provided by Ultima Genomics",
        "analysis_type":       ANALYSIS_TYPE,
        "genome_id":           GENOME_ID,
        "platform":            "UG100",
        "experiment_types":    ["Whole genome sequencing"],
        "study_accession_id":  STUDY_ACCESSION,
        "sample_accession_ids": [sample_accession],
        "files": [f["provisional_id"] for f in files],
    }

    resp = requests.post(
        f"{API_BASE}/submissions/{SUBMISSION_ID}/analyses",
        headers=auth_headers(token),
        json=payload,
    )
    if not resp.ok:
        print(f"  ✗ Failed to create analysis: {resp.status_code} {resp.text}")
        return None

    analysis = resp.json()
    if isinstance(analysis, list):
        analysis = analysis[0]
    analysis_id = analysis.get("provisional_id") or analysis.get("id")
    print(f"  ✓ Analysis created: {analysis_id}")
    return analysis_id

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python ega_bulk_analyses.py <folder>\nExample: python ega_bulk_analyses.py eso02_batch2")

    folder = sys.argv[1]

    # Extract sample and batch from folder name: "eso02_batch2" → "eso02", "batch2"
    parts = folder.split("_batch")
    if len(parts) != 2:
        sys.exit("ERROR: Folder name must follow the pattern <sample>_batch<N>, e.g. eso02_batch2")
    sample_name = parts[0]
    batch_name  = f"batch{parts[1]}"

    sample_accession = SAMPLE_MAP.get(sample_name)
    if not sample_accession:
        sys.exit(f"ERROR: Sample '{sample_name}' not found in SAMPLE_MAP. Add it and retry.")

    print(f"Folder:  {folder}")
    print(f"Sample:  {sample_name} → {sample_accession}")
    print(f"Batch:   {batch_name}")

    username, password = load_credentials()
    token = get_access_token(username, password)

    create_analysis(token, folder, sample_name, batch_name, sample_accession)

if __name__ == "__main__":
    main()
