"""
EGA Bulk Analysis Submission Script
------------------------------------
Creates REFERENCE_ALIGNMENT analyses via the EGA Submitter Portal API.

Usage:
    python Register_metadata.py <folder> [--debug]

Example:
    python Register_metadata.py eso02_batch2
    python Register_metadata.py eso02_batch2 --debug

The sample name is extracted from the folder name (everything before '_batch'),
and looked up in SAMPLE_MAP to find the corresponding EGAN accession.

Must be run from the directory containing <folder>/ (e.g.
/gpfs/commons/groups/landau_lab/ResolveOME/EGA_upload), since it reads the
local <folder>/ListFiles.txt manifest that was uploaded alongside the data.

--debug prints every file EGA has actually indexed under /<folder>/ (all
statuses, no filename filter) before the per-file manifest check runs. If
list_inbox_files() reports files "not found", this shows whether that's a
processing delay (nothing under the prefix yet) or a mismatch (files ARE
there, just under different names than ListFiles.txt expects — e.g. the
wrong manifest got uploaded alongside the data).
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

_bad_statuses = set()  # statuses the API has already rejected this run — don't retry them

def _query_files(session, prefix, status, warn_context=""):
    """
    GET /files for a given status+prefix, tolerating a status value the API
    rejects outright (e.g. "invalid input value for enum fs.status_type").
    That's a 400, not transient, so retrying won't help — log it once per
    status and treat it as "no matches" so the caller can move on to the
    next status instead of crashing the whole run.
    """
    if status in _bad_statuses:
        return []
    try:
        resp = _get_with_retry(session, f"{API_BASE}/files", params={"status": status, "prefix": prefix})
        return resp.json()
    except requests.exceptions.HTTPError as e:
        code = e.response.status_code if e.response is not None else "?"
        if code == 400:
            print(f"  ⚠ API rejected status={status!r} ({code}){' ' + warn_context if warn_context else ''} — skipping this status for the rest of the run")
            _bad_statuses.add(status)
            return []
        raise

# A single response over roughly this many files is suspicious — could be a
# response-size cap silently truncating the listing rather than a genuine
# "that's everything" answer. Just a heads-up threshold, not a hard limit.
SUSPICIOUS_RESULT_SIZE = 500

def _fetch_folder_files(session, folder):
    """
    Fetch every file EGA has indexed under /{folder}/, across all three
    statuses, in one query per status (not one query per expected filename).

    This replaced a per-file exact-path prefix query (prefix=f"/{folder}/{name}")
    that turned out to reliably return zero matches even for files confirmed
    to exist at exactly that path — --debug proved this: a folder-level
    prefix (f"/{folder}/") found all 326 files, while the same files queried
    individually by full path all came back empty. Whatever EGA's "prefix"
    param is actually matching on, it isn't a plain string-prefix over the
    full relative_path. So: fetch the folder listing once, match locally.
    """
    seen_ids = set()
    files = []
    for status in ["inbox", "submitted", "available"]:
        matches = _query_files(session, f"/{folder}/", status, warn_context="while listing all files")
        if len(matches) >= SUSPICIOUS_RESULT_SIZE:
            print(f"  ⚠ status={status!r} returned {len(matches)} file(s) — at/above the "
                  f"{SUSPICIOUS_RESULT_SIZE} sanity threshold; verify this isn't a truncated page")
        for f in matches:
            fid = f.get("provisional_id")
            if fid in seen_ids:
                continue
            seen_ids.add(fid)
            files.append(f)
    return files

def list_inbox_files(token, folder):
    """
    Fetch the actual folder listing from EGA (see _fetch_folder_files) and
    match it locally against the local ListFiles.txt manifest, by exact
    filename (after stripping .gpg — EGA drops that suffix once a file is
    in the inbox).
    """
    targets = _read_manifest_targets(folder)
    print(f"  {len(targets)} expected file(s) from ListFiles.txt")

    session = requests.Session()
    session.headers.update(auth_headers(token))

    actual = _fetch_folder_files(session, folder)
    by_path = {}
    for f in actual:
        # relative_path comes back as the FULL path (e.g.
        # "/eso05_batch2/cu08_p1r4K8.bam.bai"), not just the filename —
        # confirmed by --debug's own printout. ListFiles.txt targets are
        # bare filenames (_read_manifest_targets already strips directories
        # via os.path.basename), so normalize this side to match: strip
        # .gpg, then take the basename too.
        p = f["relative_path"]
        if p.endswith(".gpg"):
            p = p[:-4]
        p = os.path.basename(p)
        by_path.setdefault(p, []).append(f)

    seen_ids = set()
    files = []
    missing = []
    for name in targets:
        matches = by_path.get(name, [])
        if not matches:
            missing.append(name)
            continue
        for f in matches:
            fid = f.get("provisional_id")
            if fid in seen_ids:
                continue
            seen_ids.add(fid)
            files.append(f)

    if missing:
        preview = ", ".join(missing[:5])
        more = f" (+{len(missing) - 5} more)" if len(missing) > 5 else ""
        print(f"  ⚠ {len(missing)} expected file(s) not found in inbox: {preview}{more}")

    return files

def debug_list_actual_files(token, folder):
    """
    Print every file EGA has actually indexed under /{folder}/, across all
    three statuses, with no filename filtering at all.

    Useful when list_inbox_files() reports files "not found" — this prints
    ground truth: what EGA thinks is under that folder path right now,
    regardless of what ListFiles.txt says should be there. Compare the two
    lists by eye to spot a mismatch fast.
    """
    print(f"\n[DEBUG] Actual files under /{folder}/ in EGA inbox (any status):")
    session = requests.Session()
    session.headers.update(auth_headers(token))

    actual = _fetch_folder_files(session, folder)
    for f in actual:
        print(f"  {f['relative_path']}")

    if not actual:
        print(f"  (nothing indexed under /{folder}/ in any status — either "
              f"the upload hasn't landed there yet, or it went to a "
              f"different path)")
    print(f"[DEBUG] {len(actual)} total file(s) found under /{folder}/\n")

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
    args = sys.argv[1:]
    debug = "--debug" in args
    args = [a for a in args if a != "--debug"]

    if len(args) != 1:
        sys.exit("Usage: python Register_metadata.py <folder> [--debug]\nExample: python Register_metadata.py eso02_batch2")

    # Strip any trailing slash (e.g. from shell tab-completion, "eso02_batch2/"
    # -> "eso02_batch2") — a trailing slash would otherwise produce a
    # double-slash in API prefix queries later and silently match nothing.
    folder = args[0].rstrip("/")

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

    if debug:
        debug_list_actual_files(token, folder)

    create_analysis(token, folder, sample_name, batch_name, sample_accession)

if __name__ == "__main__":
    main()
