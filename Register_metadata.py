"""
EGA Bulk Analysis Submission Script
------------------------------------
Creates REFERENCE_ALIGNMENT analyses via the EGA Submitter Portal API.

Usage:
    python Register_metadata.py <folder>

Example:
    python Register_metadata.py eso02_batch2

The sample name is extracted from the folder name (everything before '_batch'),
and looked up in SAMPLE_MAP to find the corresponding EGAN accession.

Must be run from the directory containing <folder>/ (e.g.
/gpfs/commons/groups/landau_lab/ResolveOME/EGA_upload), since it reads the
local <folder>/ListFiles.txt manifest that was uploaded alongside the data.
"""

import requests
import sys
import os
import time

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

def _read_manifest_targets(folder):
    """
    Read the local <folder>/ListFiles.txt manifest (uploaded alongside the
    data) and return the set of expected post-upload filenames.

    EGA strips the .gpg suffix from relative_path once a file is in the
    inbox, so ".bam.gpg" / ".bam.bai.gpg" on disk become ".bam" / ".bam.bai"
    in the API. We only keep those two extensions here — checksum-only
    ".md5" entries are dropped since create_analysis() only wants
    .bam/.bam.bai anyway.
    """
    # Accept being run either from the parent directory (folder/ListFiles.txt)
    # or from inside the batch folder itself (./ListFiles.txt) — try both.
    candidates = [os.path.join(folder, "ListFiles.txt"), "ListFiles.txt"]
    manifest_path = next((p for p in candidates if os.path.isfile(p)), None)
    if manifest_path is None:
        raise FileNotFoundError(
            f"No manifest found (tried: {', '.join(candidates)}). This script "
            f"expects a ListFiles.txt uploaded alongside the data — run it "
            f"either from the parent directory containing {folder}/, or from "
            f"inside {folder}/ itself."
        )
    targets = set()
    with open(manifest_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name = os.path.basename(line)
            if name.endswith(".gpg"):
                name = name[:-4]
            if name.endswith(".bam") or name.endswith(".bam.bai"):
                targets.add(name)
    return sorted(targets)

def _get_with_retry(session, url, params, max_retries=5, backoff=1.5):
    """
    GET with retry/backoff on transient connection failures. EGA's API can
    reset the connection (RemoteDisconnected) under the volume of requests
    a large batch generates, so a single dropped connection shouldn't kill
    an entire 300+ file run.
    """
    for attempt in range(max_retries):
        try:
            resp = session.get(url, params=params, timeout=30)
            resp.raise_for_status()
            return resp
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
            if attempt == max_retries - 1:
                raise
            wait = backoff ** attempt
            print(f"    ⚠ request failed ({e.__class__.__name__}), retrying in {wait:.1f}s...")
            time.sleep(wait)

def list_inbox_files(token, folder):
    """
    Look up each file expected by the local ListFiles.txt manifest, by its
    exact filename. A full exact filename can only ever match 0 or 1 files
    at EGA, so this can never trip the "too much data" response-size cap
    that a broad folder-level prefix query hits on large batches — no
    guessing/splitting required.

    Each file is checked against all three possible statuses, since a file
    can move out of "inbox" (e.g. once EGA finishes checksum verification)
    before this script runs.
    """
    targets = _read_manifest_targets(folder)
    print(f"  {len(targets)} expected file(s) from ListFiles.txt")

    session = requests.Session()
    session.headers.update(auth_headers(token))

    seen_ids = set()
    files = []
    missing = []
    for i, name in enumerate(targets):
        matches = []
        for status in ["inbox", "submitted", "available"]:
            resp = _get_with_retry(
                session,
                f"{API_BASE}/files",
                params={"status": status, "prefix": f"/{folder}/{name}"},
            )
            matches = resp.json()
            if matches:
                break

        # A prefix match on "X.bam" will also match "X.bam.bai" (it's a
        # literal string-prefix of it), so filter to only the file whose
        # relative_path is exactly this target — not files that merely start
        # with it — and de-dup by provisional_id as a second safety net.
        found = False
        for f in matches:
            p = f["relative_path"]
            if p.endswith(".gpg"):
                p = p[:-4]
            if p != name:
                continue
            found = True
            fid = f.get("provisional_id")
            if fid in seen_ids:
                continue
            seen_ids.add(fid)
            files.append(f)
        if not found:
            missing.append(name)

        if (i + 1) % 50 == 0:
            print(f"    ...checked {i + 1}/{len(targets)}")
        time.sleep(0.1)  # small pacing to avoid tripping rate limits

    if missing:
        preview = ", ".join(missing[:5])
        more = f" (+{len(missing) - 5} more)" if len(missing) > 5 else ""
        print(f"  ⚠ {len(missing)} expected file(s) not found in inbox: {preview}{more}")

    return files

# ── Analysis creation ──────────────────────────────────────────────────────────

def create_analysis(token, folder, sample_name, batch_name, sample_accession):
    print(f"\nProcessing folder: {folder}")

    files = list_inbox_files(token, folder)
    if not files:
        return None

    # The EGA API strips the .gpg suffix from relative_path once files are in the
    # inbox, so paths come back as .bam / .bam.bai (not .bam.gpg / .bam.bai.gpg).
    # Strip a trailing .gpg before matching so this works either way.
    def is_bam(f):
        p = f["relative_path"]
        if p.endswith(".gpg"):
            p = p[:-4]
        return p.endswith(".bam") or p.endswith(".bam.bai")

    files = [f for f in files if is_bam(f)]
    print(f"  Found {len(files)} file(s) (excluding .md5)")
    if not files:
        return None

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
        sys.exit("Usage: python Register_metadata.py <folder>\nExample: python Register_metadata.py eso02_batch2")

    # Strip any trailing slash (e.g. from shell tab-completion, "eso02_batch2/"
    # -> "eso02_batch2") — a trailing slash would otherwise produce a
    # double-slash in API prefix queries later and silently match nothing.
    folder = sys.argv[1].rstrip("/")

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
