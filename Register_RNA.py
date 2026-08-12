"""
EGA Run/Experiment Registration Script (fastq single-cell RNA-seq)
--------------------------------------------------------------------
Registers RUN objects (with their required parent EXPERIMENT) for
per-cell paired fastq files, via the EGA Submitter Portal API.

Unlike Register_metadata.py (which creates one Analysis per donor+batch
directly on the donor Sample), Runs cannot attach to a Sample directly:
each Run must link to an Experiment, and each Experiment links to exactly
one Sample. Since our samples are donor-level but library prep/sequencing
happens per cell, this script creates ONE EXPERIMENT + ONE RUN PER CELL,
all pointing back to the same donor Sample (EGAN) accession.

After all Runs in the folder are created, their provisional IDs are merged
into DATASET_ACCESSION's run_provisional_ids automatically (one PUT per
folder, not one per run) — see link_runs_to_dataset(). Pass
--no-link-dataset to skip this and link manually in the portal instead,
the way analyses are currently linked per EGA_submission_guide.md section 6.3.

Instrument model is NOT uniform across all cells — most were sequenced on
Illumina NovaSeq X Plus, at least one donor on Illumina NovaSeq 6000. Split
the LOCAL ListFiles.txt manifests into one folder per instrument (e.g.
rna_6000_encrypted, rna_xplus_encrypted) and pass --instrument-model to
match; there's no per-cell instrument lookup, it's one value per script run
(defaults to DEFAULT_INSTRUMENT_MODEL if not passed).

If the fastqs were actually uploaded to EGA under a single shared inbox
folder (e.g. all_rna_encrypted) rather than being re-uploaded per
instrument-split folder, pass --inbox-folder to point file lookups at the
real EGA inbox path while still reading the local, instrument-specific
ListFiles.txt from <folder>. Without --inbox-folder, the local folder name
and the EGA inbox prefix are assumed to be the same (the original,
single-folder-per-batch setup).

Usage:
    python Register_RNA.py <folder> [--dry-run] [--debug] [--no-link-dataset] \
        [--limit N] [--instrument-model "..."] [--inbox-folder NAME]

Example:
    python Register_RNA.py rna_6000_encrypted --instrument-model "Illumina NovaSeq 6000" \
        --inbox-folder all_rna_encrypted --dry-run
    python Register_RNA.py rna_xplus_encrypted --instrument-model "Illumina NovaSeq X Plus" \
        --inbox-folder all_rna_encrypted --dry-run

<folder> is just the inbox folder to register — pick whatever folder you
uploaded (its ListFiles.txt manifest must live alongside it locally). The
folder NAME is not parsed for donor/batch info; which donor each cell
belongs to is determined entirely from the fastq FILENAMES inside it (see
below), so a folder can contain cells from one donor or several.

Fastq filenames are expected as <donor_filename_prefix>_<cellID>_R1/R2.fastq.gz,
where <donor_filename_prefix> has no underscores (e.g. "cu01", "eso29") and
is looked up in DONOR_FILENAME_MAP to find the donor/sample (and from there
SAMPLE_MAP to find the EGAN accession). Files are grouped into cells by
stripping the donor prefix and the trailing _R1/_R2.fastq.gz suffix; each
group becomes one Experiment + one Run.

--dry-run prints every Experiment/Run payload that would be created,
without calling the API. IMPORTANT: EGA's exact JSON field names for
Experiments/Runs have not been empirically confirmed against this API the
way the Analysis payload in Register_metadata.py has (that one was tested
against the live API). Follow the guide's own advice (EGA_submission_guide.md
section 5.4): register ONE cell manually through the portal first, compare
the fields it produces against CELL_METADATA/build_experiment_payload below,
and adjust field names here if the API rejects something. Failed POSTs
print the full response body to make this fast to diagnose.

The Dataset-linking field names (run_provisional_ids, analysis_provisional_ids,
dataset_types, policy_accession_id, etc.) are similarly unconfirmed against
this exact API — sourced from a third-party EGA submission tool's template
rather than a tested response. link_runs_to_dataset() prints the fetched
Dataset object's actual keys before building the PUT payload, specifically
so you can eyeball whether the assumed field names match reality before
anything is written — check that output on the first real run.
"""

import re
import json
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

SUBMISSION_ID    = "EGA50000001666"   # same submission as the BAM analyses
STUDY_ACCESSION  = "EGAS50000001793"  # same study as the BAM analyses
DATASET_ACCESSION = "EGAD50000002573" # same dataset the BAM analyses link to

# Sample (donor) name → EGAN accession. Same donors as Register_metadata.py's
# SAMPLE_MAP — this is the Sample every cell's Experiment will link to.
SAMPLE_MAP = {
    "eso01": "EGAN50000419612",
    "eso02": "EGAN50000419609",
    "eso03": "EGAN50000419610",
    "eso04": "EGAN50000419611",
    "eso05": "EGAN50000419613",
}

# Donor name (as used in SAMPLE_MAP) → donor code actually embedded in the
# fastq filenames, e.g. files named "cu01_<cellID>_R1.fastq.gz" belong to
# donor "eso02". This is the ONLY way donor is determined — filenames are
# read straight out of the EGA inbox listing, not the folder name.
DONOR_FILENAME_MAP = {
    "eso01": "eso29",
    "eso02": "cu01",
    "eso03": "cu02",
    "eso04": "cu03",
    "eso05": "cu08",
}
REVERSE_DONOR_MAP = {v: k for k, v in DONOR_FILENAME_MAP.items()}

# Library/platform metadata — same for every cell in this batch.
LIBRARY_STRATEGY  = "RNA-Seq"
LIBRARY_SOURCE    = "TRANSCRIPTOMIC SINGLE CELL"
LIBRARY_SELECTION = "cDNA_oligo_dT"
LIBRARY_LAYOUT    = "PAIRED"
PLATFORM          = "ILLUMINA"
DESIGN_DESCRIPTION = "scRNA-seq of SMART-PTA"  # shows as "Design Name" in the portal
# Default instrument — confirmed against the portal's own dropdown for one
# donor, but most cells were actually sequenced on a different machine (see
# --instrument-model). Split fastqs into separate folders per instrument and
# run this script once per folder, passing the matching --instrument-model.
DEFAULT_INSTRUMENT_MODEL = "Illumina NovaSeq 6000"

# instrument model name -> id, as required by the `instrument_model_id`
# column (confirmed via a live POST 400: "null value in column
# instrument_model_id ... violates not-null constraint" -- the API wants an
# int here, not the model name string we were originally sending). IDs
# confirmed via `python Register_RNA.py --dump-enum platform_models`.
# Add more entries here as needed for other instruments/donors.
INSTRUMENT_MODEL_IDS = {
    "Illumina NovaSeq 6000":   25,
    "Illumina NovaSeq X Plus": 81,
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

# ── Enum sanity check ────────────────────────────────────────────────────────

def dump_enum(token, enum_name):
    """
    GET /api/enums/{enum_name} and pretty-print the raw JSON as-is, no
    parsing/assumptions about shape. Added after 'platform_models' turned
    out to be a real, working endpoint that our value-extraction logic
    couldn't parse (it errored with "unhashable type: 'dict'" instead of
    404ing) -- this bypasses that logic entirely so we can just see the
    actual structure and figure out the right way to read it, e.g. to find
    the numeric ID EGA wants for instrument_model_id (confirmed via a real
    experiment POST that this needs an int, not the instrument name string).
    """
    print(f"GET {API_BASE}/enums/{enum_name}")
    resp = requests.get(f"{API_BASE}/enums/{enum_name}", headers=auth_headers(token), timeout=15)
    print(f"  status: {resp.status_code}")
    try:
        data = resp.json()
        print(json.dumps(data, indent=2))
    except ValueError:
        print(resp.text)

def dump_endpoint(token, path):
    """
    GET {API_BASE}/{path} and pretty-print the raw JSON, no assumptions
    about shape. General-purpose version of dump_enum() for inspecting any
    object -- e.g. `runs/1133228` (a Run we just created) or
    `analyses/<id>` (an existing Analysis already linked to the dataset via
    the portal), to see how dataset-linking is actually represented on a
    real object. Needed because the Dataset object itself only exposes
    num_runs/num_analyses (counts), not an id list -- confirmed by
    link_runs_to_dataset()'s PUT failing outright.
    """
    path = path.strip("/")
    url = f"{API_BASE}/{path}"
    print(f"GET {url}")
    resp = requests.get(url, headers=auth_headers(token), timeout=15)
    print(f"  status: {resp.status_code}")
    try:
        data = resp.json()
        print(json.dumps(data, indent=2))
    except ValueError:
        print(resp.text)

def _try_enum(token, enum_name):
    """GET /api/enums/{enum_name}. Returns (values_set, None) on success,
    (None, error_description) on failure — never raises."""
    try:
        resp = requests.get(f"{API_BASE}/enums/{enum_name}", headers=auth_headers(token), timeout=15)
        resp.raise_for_status()
        data = resp.json()
        values = {item.get("value", item) if isinstance(item, dict) else item for item in data}
        return values, None
    except requests.exceptions.HTTPError as e:
        code = e.response.status_code if e.response is not None else "?"
        body = e.response.text[:200] if e.response is not None else ""
        return None, f"HTTP {code} — {body}"
    except Exception as e:
        return None, f"{e.__class__.__name__}: {e}"

def check_platform_model(token, instrument_model):
    """
    Verify INSTRUMENT_MODEL_IDS[instrument_model] against the live
    /api/enums/platform_models list — confirmed shape is
    [{"id": int, "platform": str, "model": str}, ...], NOT a plain value
    list like the other enums (that's why the old candidate-guessing loop
    never found it: it was looking for a set of strings, not this shape).
    Checks both that the model name exists AND that our hardcoded id still
    matches, so a stale id gets caught before a real POST uses it.
    """
    configured_id = INSTRUMENT_MODEL_IDS.get(instrument_model)
    if configured_id is None:
        print(f"  ⚠ instrument model {instrument_model!r} not in INSTRUMENT_MODEL_IDS — "
              f"registration will fail before any API calls (see the check in main())")
        return
    try:
        resp = requests.get(f"{API_BASE}/enums/platform_models", headers=auth_headers(token), timeout=15)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"  ⚠ Could not check platform_models ({e.__class__.__name__}: {e}) — skipping")
        return
    by_model = {row.get("model"): row.get("id") for row in data if isinstance(row, dict)}
    live_id = by_model.get(instrument_model)
    if live_id is None:
        close = [m for m in by_model if instrument_model.lower() in (m or "").lower()
                 or (m or "").lower() in instrument_model.lower()]
        print(f"  ⚠ platform_models: {instrument_model!r} not found live. Closest matches: {close[:5]}")
    elif live_id != configured_id:
        print(f"  ⚠ platform_models: {instrument_model!r} has id={live_id} live, but "
              f"INSTRUMENT_MODEL_IDS has {configured_id} on file — update the hardcoded map")
    else:
        print(f"  ✓ platform_models: {instrument_model!r} → id {live_id} (matches INSTRUMENT_MODEL_IDS)")

def check_enums(token, instrument_model):
    """
    Query EGA's public enum lists and warn if any configured library value
    isn't present verbatim, plus check_platform_model() for the instrument.
    Doesn't block execution — just surfaces likely typos/casing mismatches
    before you burn API calls on the real registration.
    """
    print("Checking configured values against /api/enums/...")
    check_platform_model(token, instrument_model)

    checks = [
        ("library_strategies", LIBRARY_STRATEGY, ["library_strategies"]),
        ("library_sources", LIBRARY_SOURCE, ["library_sources"]),
        ("library_selections", LIBRARY_SELECTION, ["library_selections"]),
    ]
    for label, configured_value, candidates in checks:
        values = None
        errors = []
        working_name = None
        for enum_name in candidates:
            values, err = _try_enum(token, enum_name)
            if values is not None:
                working_name = enum_name
                break
            errors.append(f"{enum_name}: {err}")

        if values is None:
            print(f"  ⚠ Could not check {label} — tried {candidates}, all failed:")
            for e in errors:
                print(f"      {e}")
            continue

        if working_name != label:
            print(f"  ℹ '{label}' 404'd, but '/api/enums/{working_name}' worked instead")

        if configured_value in values:
            print(f"  ✓ {working_name}: '{configured_value}' found")
        else:
            print(f"  ⚠ {working_name}: '{configured_value}' NOT found in enum list. "
                  f"Closest values: {sorted(values)[:10]}")

# ── File listing (mirrors Register_metadata.py's approach) ────────────────────

def _read_manifest_targets(folder):
    """
    Read <folder>/ListFiles.txt and return expected post-upload filenames.
    EGA strips .gpg once a file is in the inbox, so ".fastq.gz.gpg" on disk
    becomes ".fastq.gz" in the API.
    """
    candidates = [os.path.join(folder, "ListFiles.txt"), "ListFiles.txt"]
    manifest_path = next((p for p in candidates if os.path.isfile(p)), None)
    if manifest_path is None:
        raise FileNotFoundError(
            f"No manifest found (tried: {', '.join(candidates)}). Run this either from "
            f"the parent directory containing {folder}/, or from inside {folder}/ itself."
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
            if name.endswith(".fastq.gz"):
                targets.add(name)
    return sorted(targets)

def _get_with_retry(session, url, params, max_retries=5, backoff=1.5):
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

_bad_statuses = set()

class _TooMuchDataError(Exception):
    """Raised when a query is rejected for volume, not for an invalid status
    value. Prefix-specific — must NOT be cached as a permanently bad status,
    unlike an invalid-enum rejection (see _query_files)."""
    pass

def _query_files(session, prefix, status, warn_context=""):
    if status in _bad_statuses:
        return []
    try:
        resp = _get_with_retry(session, f"{API_BASE}/files", params={"status": status, "prefix": prefix})
        return resp.json()
    except requests.exceptions.HTTPError as e:
        code = e.response.status_code if e.response is not None else "?"
        body = e.response.text[:300] if e.response is not None else ""
        if code == 400:
            # Two different 400 causes need different handling:
            #  - "invalid input value for enum fs.status_type" — this STATUS
            #    VALUE is invalid, period, for every prefix. Safe to cache
            #    and skip for the rest of the run (true for submitted/available).
            #  - "too much data" — this PREFIX's result set is too large for
            #    this status. Prefix-specific: caching `status` as bad here
            #    would incorrectly skip it for every OTHER (possibly fine)
            #    prefix too. This bit us for real: cu01's volume error
            #    silently suppressed eso29's 'inbox' query right after it.
            if "invalid input value for enum" in body:
                print(f"  ⚠ /files?status={status!r}&prefix={prefix!r} returned 400{' ' + warn_context if warn_context else ''}: "
                      f"{body} — invalid status value, skipping it for the rest of the run")
                _bad_statuses.add(status)
                return []
            if "too much data" in body.lower():
                raise _TooMuchDataError(f"status={status!r} prefix={prefix!r}: {body}")
            print(f"  ⚠ /files?status={status!r}&prefix={prefix!r} returned 400{' ' + warn_context if warn_context else ''}: {body} — skipping this query only")
            return []
        raise

SUSPICIOUS_RESULT_SIZE = 500
# How many extra characters to extend a prefix by each time a query comes
# back "too much data". Small steps = more API calls but finer-grained
# splitting; kept small since donor-level alone wasn't enough (cu01 still
# too big), so we know we need to descend into the cell ID.
PREFIX_STEP = 3

def _fetch_prefix_recursive(session, folder, prefix_suffix, remaining_suffixes, warn_context=""):
    """
    Query /{folder}/{prefix_suffix}/... and, on a "too much data" rejection,
    adaptively narrow by extending prefix_suffix with more characters from
    remaining_suffixes (the not-yet-matched tail of each target filename
    past prefix_suffix) and recursing — rather than assuming one fixed
    split (e.g. donor-only) is narrow enough.

    Deliberately refuses to narrow into a filename's trailing _R1/_R2.fastq.gz
    suffix: a fully-specific, single-file prefix query has separately been
    observed on this API to reliably return zero matches even for files
    confirmed to exist (see Register_metadata.py's _fetch_folder_files
    comment) — narrowing that far would trade "too much data" for "no data
    at all", not actually fix anything.
    """
    full_prefix = f"/{folder}/{prefix_suffix}"
    results = []
    for status in ["inbox", "submitted", "available"]:
        try:
            matches = _query_files(session, full_prefix, status,
                                    warn_context=warn_context or f"while listing files under {full_prefix!r}")
        except _TooMuchDataError:
            safe_suffixes = [s for s in remaining_suffixes if s and not s.startswith("_R")]
            if not safe_suffixes:
                print(f"  ⚠ {full_prefix!r} ({status}) still too much data at the safe narrowing limit "
                      f"(won't narrow into the _R1/_R2 read suffix — exact-filename prefixes have "
                      f"reliably returned zero matches on this API before) — skipping")
                continue
            buckets = {}
            for s in safe_suffixes:
                key = s[:PREFIX_STEP]
                buckets.setdefault(key, []).append(s[PREFIX_STEP:])
            for key, sub_suffixes in buckets.items():
                results.extend(_fetch_prefix_recursive(session, folder, prefix_suffix + key, sub_suffixes, warn_context))
            continue
        if len(matches) >= SUSPICIOUS_RESULT_SIZE:
            print(f"  ⚠ prefix={full_prefix!r} status={status!r} returned {len(matches)} file(s) — at/above the "
                  f"{SUSPICIOUS_RESULT_SIZE} sanity threshold; verify this isn't a truncated page")
        results.extend(matches)
    return results

def _fetch_folder_files(session, folder, sub_prefixes=None, wanted_names=None):
    """
    Query EGA for files under /{folder}/, across all three statuses.

    If sub_prefixes + wanted_names are given, queries per sub-prefix and
    adaptively narrows further within each one on "too much data" (see
    _fetch_prefix_recursive) — needed for large shared inbox folders (e.g.
    all_rna_encrypted, ~10k files) where even a single donor's worth of
    files can still be too large for one query. Without wanted_names there's
    nothing to split on, so a volume rejection is just reported as-is (used
    by the --debug listing, which has no manifest to narrow against).
    Results are deduplicated by provisional_id across all sub-queries.
    """
    seen_ids = set()
    files = []
    sub_prefixes = sub_prefixes or [""]
    for sub in sub_prefixes:
        if wanted_names is not None:
            remaining = [n[len(sub):] for n in wanted_names if n.startswith(sub)]
        else:
            remaining = []
        matches = _fetch_prefix_recursive(session, folder, sub, remaining,
                                           warn_context=f"while listing files under '/{folder}/{sub}'")
        for f in matches:
            fid = f.get("provisional_id")
            if fid in seen_ids:
                continue
            seen_ids.add(fid)
            files.append(f)
    return files

def list_inbox_files(token, folder, inbox_prefix=None):
    """
    inbox_prefix lets the EGA inbox path being queried differ from the local
    folder holding ListFiles.txt — needed when fastqs were uploaded under one
    shared inbox folder (e.g. all_rna_encrypted) but you've since split the
    local manifests into per-instrument folders (e.g. rna_6000_encrypted)
    that were never actually re-uploaded under their own inbox path.
    Defaults to folder when not given, matching the old single-folder setup.
    """
    inbox_prefix = inbox_prefix or folder
    targets = _read_manifest_targets(folder)
    print(f"  {len(targets)} expected file(s) from ListFiles.txt")

    # Query per-donor-prefix rather than the whole inbox_prefix at once —
    # necessary when inbox_prefix is a large shared folder (see
    # _fetch_folder_files docstring). Harmless and slightly slower on small
    # folders too, so just always do it.
    donor_prefixes = sorted({t.split("_", 1)[0] for t in targets})
    print(f"  Querying {len(donor_prefixes)} donor prefix(es) under /{inbox_prefix}/: {donor_prefixes}")

    session = requests.Session()
    session.headers.update(auth_headers(token))

    actual = _fetch_folder_files(session, inbox_prefix, sub_prefixes=donor_prefixes, wanted_names=targets)
    by_path = {}
    for f in actual:
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
        print(f"  ⚠ {len(missing)} expected file(s) not found under inbox prefix /{inbox_prefix}/: {preview}{more}")

    return files

def debug_list_actual_files(token, folder, inbox_prefix=None):
    inbox_prefix = inbox_prefix or folder
    print(f"\n[DEBUG] Actual files under /{inbox_prefix}/ in EGA inbox (any status):")
    session = requests.Session()
    session.headers.update(auth_headers(token))
    actual = _fetch_folder_files(session, inbox_prefix)
    for f in actual:
        print(f"  {f['relative_path']}")
    if not actual:
        print(f"  (nothing indexed under /{inbox_prefix}/ in any status)")
    print(f"[DEBUG] {len(actual)} total file(s) found under /{inbox_prefix}/\n")

# ── Group fastqs by cell ────────────────────────────────────────────────────

# <donor_prefix>_<cellID>_R1/R2.fastq.gz — donor_prefix has no underscores
# (matches every value in DONOR_FILENAME_MAP: eso29, cu01, cu02, cu03, cu08).
READ_SUFFIX_RE = re.compile(r"^(?P<donor_prefix>[^_]+)_(?P<cell_id>.+)_R(?P<read_num>[12])\.fastq\.gz$")

def group_by_cell(files):
    """
    Group fastq files by cell, resolving each cell's donor straight from its
    filename (via REVERSE_DONOR_MAP) — the folder they live in is irrelevant
    to which donor/sample they belong to.

    Returns {cell_key: {"R1": file, "R2": file, "sample_name": str,
    "sample_accession": str}} for cells with both reads present and a
    recognized donor prefix.
    """
    cells = {}
    unmatched = []
    unknown_donor = []
    for f in files:
        name = os.path.basename(f["relative_path"].rstrip("/"))
        if name.endswith(".gpg"):
            name = name[:-4]
        m = READ_SUFFIX_RE.match(name)
        if not m:
            unmatched.append(name)
            continue
        donor_prefix = m.group("donor_prefix")
        cell_key     = f"{donor_prefix}_{m.group('cell_id')}"
        read_num     = m.group("read_num")

        sample_name = REVERSE_DONOR_MAP.get(donor_prefix)
        if not sample_name:
            unknown_donor.append(name)
            continue
        sample_accession = SAMPLE_MAP.get(sample_name)
        if not sample_accession:
            print(f"  ⚠ '{name}': donor '{sample_name}' has no EGAN accession in SAMPLE_MAP — skipping")
            continue

        entry = cells.setdefault(cell_key, {"sample_name": sample_name, "sample_accession": sample_accession})
        entry[f"R{read_num}"] = f

    if unmatched:
        preview = ", ".join(unmatched[:5])
        more = f" (+{len(unmatched) - 5} more)" if len(unmatched) > 5 else ""
        print(f"  ⚠ {len(unmatched)} file(s) didn't match <donor>_<cell>_R1/R2.fastq.gz pattern: {preview}{more}")

    if unknown_donor:
        preview = ", ".join(unknown_donor[:5])
        more = f" (+{len(unknown_donor) - 5} more)" if len(unknown_donor) > 5 else ""
        known = sorted(REVERSE_DONOR_MAP)
        print(f"  ⚠ {len(unknown_donor)} file(s) had an unrecognized donor prefix (known: {known}): {preview}{more}")

    incomplete = {k: v for k, v in cells.items() if "R1" not in v or "R2" not in v}
    for cell_key, entry in incomplete.items():
        reads_present = [k for k in ("R1", "R2") if k in entry]
        print(f"  ⚠ Cell '{cell_key}' has {len(reads_present)}/2 reads ({reads_present}) — skipping, needs both R1 and R2")

    complete = {k: v for k, v in cells.items() if "R1" in v and "R2" in v}
    return complete

# ── Experiment + Run creation ──────────────────────────────────────────────────

def build_experiment_payload(cell_key, sample_accession, instrument_model):
    instrument_model_id = INSTRUMENT_MODEL_IDS.get(instrument_model)
    if instrument_model_id is None:
        raise ValueError(
            f"Unknown instrument model {instrument_model!r} — not in INSTRUMENT_MODEL_IDS "
            f"({sorted(INSTRUMENT_MODEL_IDS)}). Run `python Register_RNA.py --dump-enum "
            f"platform_models` to find its id and add it there."
        )
    return {
        "alias":               f"{cell_key}_experiment",
        "title":               f"scRNA-seq experiment for cell {cell_key}",
        "description":         f"Single-cell RNA-seq library for cell {cell_key}, donor sample {sample_accession}",
        "study_accession_id":  STUDY_ACCESSION,
        "sample_accession_id": sample_accession,
        "design_description":  DESIGN_DESCRIPTION,
        "library_name":        cell_key,
        "library_strategy":    LIBRARY_STRATEGY,
        "library_source":      LIBRARY_SOURCE,
        "library_selection":   LIBRARY_SELECTION,
        "library_layout":      LIBRARY_LAYOUT,
        "platform":            PLATFORM,
        "instrument_model_id": instrument_model_id,
    }

def create_experiment(token, cell_key, sample_accession, instrument_model, dry_run=False):
    payload = build_experiment_payload(cell_key, sample_accession, instrument_model)
    if dry_run:
        print(f"  [dry-run] would POST /submissions/{SUBMISSION_ID}/experiments:\n    {payload}")
        return "DRY-RUN-EXPERIMENT-ID"

    resp = requests.post(
        f"{API_BASE}/submissions/{SUBMISSION_ID}/experiments",
        headers=auth_headers(token),
        json=payload,
    )
    if not resp.ok:
        print(f"  ✗ Failed to create experiment for {cell_key}: {resp.status_code} {resp.text}")
        return None

    experiment = resp.json()
    if isinstance(experiment, list):
        experiment = experiment[0]
    experiment_id = experiment.get("provisional_id") or experiment.get("id")
    print(f"  ✓ Experiment created for {cell_key}: {experiment_id}")
    return experiment_id

def build_run_payload(cell_key, experiment_id, entry):
    return {
        "alias":                     f"{cell_key}_run",
        "experiment_provisional_id": experiment_id,
        # Confirmed via a live 400 ("Sample (either provisional or
        # accession) is required") that Run needs its own direct sample
        # link -- not just inherited implicitly through the Experiment.
        "sample_accession_id":       entry["sample_accession"],
        "run_file_type":             "fastq",
        "files":                     [entry["R1"]["provisional_id"], entry["R2"]["provisional_id"]],
    }

def create_run(token, cell_key, experiment_id, entry, dry_run=False):
    payload = build_run_payload(cell_key, experiment_id, entry)
    if dry_run:
        print(f"  [dry-run] would POST /submissions/{SUBMISSION_ID}/runs:\n    {payload}")
        return "DRY-RUN-RUN-ID"

    resp = requests.post(
        f"{API_BASE}/submissions/{SUBMISSION_ID}/runs",
        headers=auth_headers(token),
        json=payload,
    )
    if not resp.ok:
        print(f"  ✗ Failed to create run for {cell_key}: {resp.status_code} {resp.text}")
        return None

    run = resp.json()
    if isinstance(run, list):
        run = run[0]
    run_id = run.get("provisional_id") or run.get("id")
    print(f"  ✓ Run created for {cell_key}: {run_id}")
    return run_id

# ── Dataset linking ─────────────────────────────────────────────────────────

def get_dataset(token):
    """
    Find DATASET_ACCESSION's object within this submission. Mirrors the
    /submissions/{id}/analyses listing pattern already used elsewhere in
    this codebase, applied to the "datasets" catalog object.
    """
    resp = requests.get(f"{API_BASE}/submissions/{SUBMISSION_ID}/datasets", headers=auth_headers(token))
    resp.raise_for_status()
    datasets = resp.json()
    for d in datasets:
        if d.get("accession_id") == DATASET_ACCESSION:
            return d
    found = [d.get("accession_id") for d in datasets]
    raise RuntimeError(f"Dataset {DATASET_ACCESSION} not found under submission {SUBMISSION_ID}. Found: {found}")

def link_runs_to_dataset(token, new_run_ids, dry_run=False):
    """
    Merge new_run_ids into DATASET_ACCESSION's run_provisional_ids (one PUT
    for the whole batch, not one per run) so you don't have to click through
    hundreds of runs in the portal by hand.
    """
    new_run_ids = [r for r in new_run_ids if r]
    if not new_run_ids:
        print("\nNo successfully created runs to link to the dataset.")
        return

    print(f"\nLinking {len(new_run_ids)} run(s) to dataset {DATASET_ACCESSION}...")

    if dry_run:
        print(f"  [dry-run] would fetch dataset {DATASET_ACCESSION} and PUT it back with "
              f"{len(new_run_ids)} new run_provisional_id(s) merged in")
        return

    dataset = get_dataset(token)
    print(f"  Fetched dataset object — keys present: {sorted(dataset.keys())}")
    dataset_id = dataset.get("provisional_id") or dataset.get("id")

    existing_run_ids = set(dataset.get("run_provisional_ids") or [])
    merged_run_ids = sorted(existing_run_ids | set(new_run_ids))

    payload = {
        "title":                     dataset.get("title"),
        "description":               dataset.get("description"),
        "dataset_types":             dataset.get("dataset_types"),
        "policy_accession_id":       dataset.get("policy_accession_id"),
        "run_provisional_ids":       merged_run_ids,
        "run_accession_ids":         dataset.get("run_accession_ids") or [],
        "analysis_provisional_ids":  dataset.get("analysis_provisional_ids") or [],
        "analysis_accession_ids":    dataset.get("analysis_accession_ids") or [],
    }

    resp = requests.put(f"{API_BASE}/datasets/{dataset_id}", headers=auth_headers(token), json=payload)
    if not resp.ok:
        print(f"  ✗ Failed to link runs to dataset: {resp.status_code} {resp.text}")
        print(f"  → You'll need to link these {len(new_run_ids)} run(s) manually in the portal instead.")
        return

    print(f"  ✓ Linked {len(new_run_ids)} new run(s) to dataset {DATASET_ACCESSION} "
          f"({len(merged_run_ids)} total run(s) now linked)")

def register_cells(token, folder, instrument_model, dry_run=False, link_dataset=True, limit=None, inbox_prefix=None):
    print(f"\nProcessing folder: {folder}")
    if inbox_prefix and inbox_prefix != folder:
        print(f"  EGA inbox prefix: /{inbox_prefix}/ (differs from local folder — manifest read locally, files looked up under this prefix)")
    print(f"  Instrument model for this run: {instrument_model}")
    files = list_inbox_files(token, folder, inbox_prefix=inbox_prefix)
    if not files:
        return

    # list_inbox_files() already restricted these to manifest-matched
    # .fastq.gz names (see _read_manifest_targets), so .md5/.gpg.md5
    # siblings never reach this point — no extra filtering needed here.
    print(f"  Found {len(files)} fastq file(s) in inbox")

    cells = group_by_cell(files)
    print(f"  {len(cells)} complete cell(s) (R1+R2 both present, donor recognized)")
    by_donor = {}
    for cell_key, entry in cells.items():
        by_donor.setdefault(entry["sample_name"], []).append(cell_key)
    for sample_name, cell_keys in sorted(by_donor.items()):
        print(f"    {sample_name}: {len(cell_keys)} cell(s)")

    items = sorted(cells.items())
    if limit is not None:
        print(f"  --limit {limit} set — only registering the first {limit} cell(s): "
              f"{[k for k, _ in items[:limit]]}")
        items = items[:limit]

    created_run_ids = []
    for cell_key, entry in items:
        experiment_id = create_experiment(token, cell_key, entry["sample_accession"], instrument_model, dry_run=dry_run)
        if experiment_id is None:
            continue
        run_id = create_run(token, cell_key, experiment_id, entry, dry_run=dry_run)
        created_run_ids.append(run_id)

    if link_dataset:
        link_runs_to_dataset(token, created_run_ids, dry_run=dry_run)
    else:
        print(f"\n--no-link-dataset set — skipping dataset link for {len(created_run_ids)} run(s). "
              f"Link them manually to {DATASET_ACCESSION} in the portal.")

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = sys.argv[1:]

    if "--dump-enum" in args:
        i = args.index("--dump-enum")
        try:
            enum_name = args[i + 1]
        except IndexError:
            sys.exit("ERROR: --dump-enum requires a name, e.g. --dump-enum platform_models")
        username, password = load_credentials()
        token = get_access_token(username, password)
        dump_enum(token, enum_name)
        return

    if "--dump-endpoint" in args:
        i = args.index("--dump-endpoint")
        try:
            path = args[i + 1]
        except IndexError:
            sys.exit("ERROR: --dump-endpoint requires a path, e.g. --dump-endpoint runs/1133228")
        username, password = load_credentials()
        token = get_access_token(username, password)
        dump_endpoint(token, path)
        return

    debug = "--debug" in args
    dry_run = "--dry-run" in args
    link_dataset = "--no-link-dataset" not in args

    limit = None
    if "--limit" in args:
        i = args.index("--limit")
        try:
            limit = int(args[i + 1])
        except (IndexError, ValueError):
            sys.exit("ERROR: --limit requires an integer argument, e.g. --limit 2")
        args = args[:i] + args[i + 2:]

    instrument_model = DEFAULT_INSTRUMENT_MODEL
    if "--instrument-model" in args:
        i = args.index("--instrument-model")
        try:
            instrument_model = args[i + 1]
        except IndexError:
            sys.exit("ERROR: --instrument-model requires a value, e.g. --instrument-model \"Illumina NovaSeq X Plus\"")
        args = args[:i] + args[i + 2:]

    if instrument_model not in INSTRUMENT_MODEL_IDS:
        sys.exit(f"ERROR: instrument model {instrument_model!r} not in INSTRUMENT_MODEL_IDS "
                  f"({sorted(INSTRUMENT_MODEL_IDS)}). Run `python Register_RNA.py --dump-enum "
                  f"platform_models` to find its id and add it there before retrying.")

    inbox_prefix = None
    if "--inbox-folder" in args:
        i = args.index("--inbox-folder")
        try:
            inbox_prefix = args[i + 1]
        except IndexError:
            sys.exit("ERROR: --inbox-folder requires a value, e.g. --inbox-folder all_rna_encrypted")
        args = args[:i] + args[i + 2:]

    args = [a for a in args if a not in ("--debug", "--dry-run", "--no-link-dataset")]

    if len(args) != 1:
        sys.exit("Usage: python Register_RNA.py <folder> [--dry-run] [--debug] [--no-link-dataset] "
                  "[--limit N] [--instrument-model \"...\"] [--inbox-folder NAME]\n"
                  "Example: python Register_RNA.py rna_6000_encrypted --instrument-model \"Illumina NovaSeq 6000\" "
                  "--inbox-folder all_rna_encrypted\n"
                  f"(default instrument model if not passed: \"{DEFAULT_INSTRUMENT_MODEL}\")\n"
                  "(--inbox-folder: use when the local ListFiles.txt folder differs from the folder\n"
                  " the fastqs were actually uploaded under in the EGA inbox — defaults to <folder>)")

    folder = args[0].rstrip("/")
    if inbox_prefix:
        inbox_prefix = inbox_prefix.rstrip("/")

    print(f"Folder: {folder}")
    print(f"Instrument model: {instrument_model}")
    if inbox_prefix:
        print(f"Inbox folder override: {inbox_prefix}")
    if dry_run:
        print("DRY RUN — no API writes will be made")
    if not link_dataset:
        print("--no-link-dataset set — runs will NOT be linked to the dataset automatically")
    if limit is not None:
        print(f"--limit {limit} set — only the first {limit} cell(s) will be registered")

    username, password = load_credentials()
    token = get_access_token(username, password)

    if debug:
        debug_list_actual_files(token, folder, inbox_prefix=inbox_prefix)

    check_enums(token, instrument_model)

    register_cells(token, folder, instrument_model, dry_run=dry_run, link_dataset=link_dataset,
                    limit=limit, inbox_prefix=inbox_prefix)

if __name__ == "__main__":
    main()
