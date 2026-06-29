# EGA Submission Guide — Landau Lab

This guide documents the full process for submitting genomic data to the [European Genome-phenome Archive (EGA)](https://ega-archive.org/).

It is specifically aimed at **SMART-PTA or single-cell whole-genome sequencing datasets** where each sample generates hundreds to thousands of high-depth BAM files. Manually loading files one by one through the EGA web portal is not feasible at this scale, so this guide covers both the portal workflow (for study and sample registration) and a programmatic approach using the EGA Submitter Portal API (for batch analysis registration).

---

## 1. Account Setup

### 1.1 Create a Personal EGA Account
[Register here](https://ega-archive.org/register/) and wait for account verification.

### 1.2 Obtain a Data Processing Agreement (DPA)
The DPA requires a signature from a legal representative at Cornell. You will need to provide:
- A title for the dataset
- Data elements description
- Number of subjects
- An IRB number

### 1.3 Set Up a DAC (Data Access Committee)


>Ask **Jane and Catherine** if you need help with DPA or DAC.

### 1.4 Request a Submitter Role
Once your account is validated, [request a submitter role](https://profile.ega-archive.org/submitter-request).

> **Note:** This step may not be required if you are granted access to the lab shared box instead (see below).

### 1.5 Request Access to the Lab Shared Box
Open a [Helpdesk ticket](https://ega-archive.org/need-help/) to request access to the lab common box: **ega-box-2193**.

EGA will send a PDF by email requesting:
- Authorization from **Dan**
- A list of all current users to be added

---

## 2. Lab EGA Box Credentials

Credentials can be stored on the cluster at:

```bash
cat ~/.ega_credentials
```

The file contains:
```bash
export EGA_BOX=ega-box-2193
export ASPERA_SCP_PASS=<password>
```

---

## 3. Directory Structure and Scripts

All EGA upload work lives at:

```
/gpfs/commons/groups/landau_lab/ResolveOME/EGA_upload/
```

Key files and folders:

```
EGA_upload/
├── data_paths/                               # Per-donor file lists (input to encryption)
│   └── <donor>_batches/
│       └── bam.list.<donor>.batch<N>.txt     # du-format: <size_kb>  <file_path>
├── EGA-Cryptor-2.0.0/                        # EGA encryption tool (Java)
├── <donor>_batch<N>/                         # Encrypted output folders (one per batch)
│   ├── *.bam.gpg                             # Encrypted BAM files
│   ├── *.bam.bai.gpg                         # Encrypted BAM index files
│   ├── *.md5                                 # Checksums (not submitted to EGA)
│   └── ListFiles.txt                         # List of files in this batch
├── split_by_size.py                          # Split file lists into ≤N TB batches
├── Upload_batch.sh                           # SLURM job: encrypt one batch
├── Upload_batch_missing.sh                   # SLURM job: re-encrypt failed files
├── Register_metadata.py                      # API script: create analysis for a batch
└── Register_metadata_manual.py              # API script: add files to existing analysis
```

**`split_by_size.py`** splits a full `du`-format BAM list into batches of up to 10 TB each, keeping each `.bam` and `.bam.bai` pair together. Run this first if your donor has more data than fits in a single transfer.

```bash
python3 split_by_size.py data_paths/bam.list.eso01.txt --max-tb 10
# Outputs: bam.list.eso01.batch1.txt, bam.list.eso01.batch2.txt, ...
```

**`Upload_batch.sh`** is a SLURM job that encrypts all BAM files for a given donor and batch using EGA-Cryptor. Submit with `--array` to specify the batch number:

```bash
sbatch --array=1 Upload_batch.sh eso01     # encrypts batch 1 for eso01
sbatch --array=1-5 Upload_batch.sh eso01   # encrypts batches 1 through 5
```

**`Upload_batch_missing.sh`** does the same but only re-encrypts files that failed or are missing from a previous run.

---

## 4. Encrypt the Data (~12 hours for 10TB)

Navigate to the upload directory and submit the encryption job:

```bash
cd /gpfs/commons/groups/landau_lab/ResolveOME/EGA_upload
sbatch --array=1 Upload_batch.sh eso01
```

Files should be encrypted with [Crypt4GH](https://ega-archive.org/submission/data/file-preparation/crypt4gh/) before upload. Encrypted files will have a `.gpg` extension.

---

## 4. Upload the Data

### Option A: Globus (recommended for large datasets)

Aspera was confirmed not working by the EGA team at the time of writing. Use Globus instead (see RESCOMP ticket [RESCOMP-20030](https://jira.nygenome.org/browse/RESCOMP-20030)).

**Setup:**
1. Log in to Globus using the Google Authenticator app
2. Set up the endpoint using the PDT: [NYGC Globus login instructions](https://wiki.nygenome.org/spaces/rescomp/pages/165609778/How+to+login+to+Globus+with+SSO)
3. Go to [Globus File Manager](https://app.globus.org/file-manager/gcp) and click **FILE MANAGER**

**Left panel (source — NYGC):**

- Data under `gpfs/commons/groups`: use the [groups endpoint](https://app.globus.org/file-manager/collections/529df0b0-5a48-4266-ac76-9feba72449a6/overview). Remove everything before `landau_lab` to navigate to your folder.
- Data under your home directory: use the [home endpoint](https://app.globus.org/file-manager/collections/c50f3e55-0997-4c13-92a4-ba7a4b8d74ad)

**Right panel (destination):** EMBL-EBI Private Collection

You will receive an email once the Globus transfer completes. Transfer speed is approximately 1 GB/hour, so plan accordingly for large datasets.

### Option B: FTP (small datasets only)

```bash
ftp ftp.ega.ebi.ac.uk
```

See [EGA FTP instructions](https://ega-archive.org/submission/data/uploading-files/ftp/).

---

## 5. Register Metadata in the Submitter Portal

Before registering analyses programmatically, the study and samples must be created manually through the portal. Go to [https://submission.ega-archive.org/submissions](https://submission.ega-archive.org/submissions).

### 5.1 Create a Submission

Click the **"Create a submission"** green button at the top right. This creates a submission object that will hold all your metadata (study, samples, analyses, dataset).

### 5.2 Register a Study

Inside the submission, go to **Register Study** and fill in the required fields (title, description, study type). This will generate an `EGAS...` accession.

### 5.3 Register Samples

Go to **Register Samples** and fill in the sample metadata interactively through the portal UI. Each sample will receive an `EGAN...` accession. Keep a record of the sample alias → EGAN mapping (you will need it for analysis registration).

### 5.4 Register Analyses via the Portal (first batch — recommended)

For your first batch, submit the analysis manually through the portal to verify all fields are correct before automating. Inside the submission:

1. Go to **Register Analysis**
2. Set type to **Reference Alignment**
3. Fill in title, description, genome, platform, and experiment type
4. Link the study and sample
5. Select the files from the inbox

> **Tip:** Use the first batch as a sanity check — verify that the analysis looks correct in the portal before running the script on all remaining batches. Once you are happy with the metadata, it is worth finalising the submission early to obtain an accession number (e.g. for a manuscript). You can continue adding analyses to the submission after finalisation while the dataset remains private pending data access approval. I recommend finalizing the submission so you can obtain an accession number. Then you can keep adding more analyses while the project remains public

---

## 6. Register Metadata Programmatically

Once files are uploaded, register metadata using the EGA Submitter Portal API. The scripts below automate analysis creation for the current study.

**Current study:** `EGAS50000001793`  
**Submission ID:** `EGA50000001666`  
**Dataset:** `EGAD50000002573`

### 6.1 Create a new analysis for a batch folder

```bash
python Register_metadata.py eso02_batch3
```

The script:
- Extracts the sample name from the folder (e.g. `eso02_batch3` → `eso02`)
- Looks up the corresponding EGAN accession
- Finds all `.bam` and `.bam.bai` files in the inbox under that folder
- Creates a `REFERENCE ALIGNMENT` analysis linked to the study and sample

> **Tip:** Submit the first batch manually via the portal (see section 5.4) to verify everything looks correct before running this script.

**Sample → EGAN accession map:**

| Sample | EGAN accession     |
|--------|--------------------|
| eso01  | EGAN50000419612    |
| eso02  | EGAN50000419609    |
| eso03  | EGAN50000419610    |
| eso04  | EGAN50000419611    |
| eso05  | EGAN50000419613    |

### 6.2 Add files to an existing analysis

If an analysis was already created and you need to update its file list:

```bash
python Register_metadata_manual.py <analysis_provisional_id> <folder>
```

Example:
```bash
python Register_metadata_manual.py 128330 eso05_batch1
```

### 6.3 Link analyses to the dataset

After creating analyses, link them to `EGAD50000002573` manually in the [Submitter Portal](https://submission.ega-archive.org/).

---

## 7. Finalise Submission

Once all analyses are registered and linked, open a [Helpdesk ticket](https://ega-archive.org/need-help/) to notify the EGA team for review and finalisation.

> **Note on ingestion:** File upload and file ingestion are separate steps. Ingestion — where EGA processes and validates your files into the archive — only begins after the submission is finalised. This can take up to 48 hours. You can monitor file status at [https://submission.ega-archive.org/files/archive](https://submission.ega-archive.org/files/archive). The submission will remain in a pending/review state until the helpdesk approves it.
