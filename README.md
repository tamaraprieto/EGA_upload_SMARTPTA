# EGA Submission Guide — Landau Lab

This guide documents the full process for submitting genomic data to the [European Genome-phenome Archive (EGA)](https://ega-archive.org/).

It is specifically aimed at **SMART-PTA or single-cell whole-genome sequencing datasets** where each sample generates hundreds to thousands of high-depth BAM files. Manually loading files one by one through the EGA web portal is not feasible at this scale, so this guide covers both the portal workflow (for study and sample registration) and a programmatic approach using the EGA Submitter Portal API (for batch analysis registration). Ultima Genomic delivers mapped alignments (BAM files) rather than raw reads. For that reason, the raw data will be deposited as ANALYSIS (Reference Alignment type) rather than RUNS.

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
Open a [Helpdesk ticket](https://ega-archive.org/need-help/) to request access to the lab common box: ega-box-2193. The EGA inbox has a default quota of 10 TB, so large submissions must be uploaded in batches. In practice there may be no hard limit, but if you anticipate exceeding 10 TB it is safer to contact the EGA helpdesk (once the access to be box is granted) to avoid being blocked.

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


## 3. Workflow

![EGA Submission Workflow](ega_workflow.svg)

Previous EGA upload work lives at:
```
/gpfs/commons/groups/landau_lab/ResolveOME/EGA_upload/
```

### 3.1 List size and full path of the bam and bai files

Create a list of files to be uploaded per donor (sample) annotated with their size
```
mkdir data_paths
du /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/*ba* > data_paths/bam.list.eso02.txt 
```
Manually check the lists as you don't want to upload incorrect files

File content example
```
78072160        /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/cu01_p1r1A10.bam
9400    /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/cu01_p1r1A10.bam.bai
77651120        /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/cu01_p1r1A11.bam
9384    /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/cu01_p1r1A11.bam.bai
73766488        /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/cu01_p1r1A12.bam
9376    /gpfs/commons/groups/landau_lab/ResolveOME/StartDir/cu01/dna/merged/cu01_p1r1A12.bam.bai
```


### 3.2 Create lists of files adding 10TB or less

**`split_by_size.py`** splits BAM list into batches of up to 10 TB each, keeping each `.bam` and `.bam.bai` pair together.

```bash
python3 split_by_size.py data_paths/bam.list.eso01.txt --max-tb 10
# Outputs: bam.list.eso01.batch1.txt, bam.list.eso01.batch2.txt, ...
```

### 3.3 Encrypt a batch (~12 hours for 10TB)
Files should be encrypted with [Crypt4GH](https://ega-archive.org/submission/data/file-preparation/crypt4gh/) before upload. Encrypted files will have a `.gpg` extension.

```bash
sbatch --array=1 Encrypt_batch.sh eso01     # encrypts batch 1 for eso01
sbatch --array=1-5 Encrypt_batch.sh eso01   # encrypts batches 1 through 5
```

### 3.4 Upload the Data

Once your files have been encrypted you can upload them to the box. 

If your dataset is small, you can upload them using [FTP](https://ega-archive.org/submission/data/uploading-files/ftp/).

```bash
ftp ftp.ega.ebi.ac.uk
```
However, if you have many batches, uploading through FTP could take months. Aspera is an alternative, but the EGA team confirmed it wasn't working at the time of writing. Use **Globus** instead.

The EGA team sent me the information on how to use Globus via email The information is not in their website. See RESCOMP ticket [RESCOMP-20030](https://jira.nygenome.org/browse/RESCOMP-20030) for details on how to set it up. You will need to set up the endpoint using the [NYGC Globus login instructions](https://wiki.nygenome.org/spaces/rescomp/pages/165609778/How+to+login+to+Globus+with+SSO) and the **Google Authenticator app** or equivalent.    

Then for each batch upload:
1. **Log in into [Globus](https://app.globus.org/file-manager/gcp)**. Your username is not your complete nygc email, only the name initial and last name. E.g. tprieto. Open the **Google Authenticator app** to obtain the one-time-code.  
3.  Click on **FILE MANAGER** at the left side bar.
4. Select a **Collection** under the left upper panel. This will be the folder(s) your want to upload.

- Data under `gpfs/commons/groups`: use the [groups endpoint](https://app.globus.org/file-manager/collections/529df0b0-5a48-4266-ac76-9feba72449a6/overview). Remove everything before `landau_lab` to navigate to your folder under **Path**.
- Data under your home directory: use the [home endpoint](https://app.globus.org/file-manager/collections/c50f3e55-0997-4c13-92a4-ba7a4b8d74ad)

5. Then select the **EMBL-EBI Private Collection** as **Collection** on the right panel. This is where you want to upload the data to (the EGA box). You might need to authenticate using your EGA box credentials. 

6. Press **`Start`** under the left-side panel. 
  
You will receive an email once the Globus transfer completes. Transfer speed is approximately 1 GB/hour, so plan accordingly for large datasets.

Once the upload has been completed, this **triggers ingestion on EGA's end (up to 48 hours)**.  

### 3.5a Register Metadata in the Submitter Portal

Before registering analyses programmatically, the study, samples, experiments, and dataset must be created manually through the portal. Go to [https://submission.ega-archive.org/submissions](https://submission.ega-archive.org/submissions).

I also recommend adding the first analysis (covering the first batch, up to 10TB) manually to understand how the programmatic registration of subsequent analyses works under the hood. If you have fewer than a few hundred files, you might choose to create analyses manually only. 

Click the **"Create a submission"** green button at the top right. Follow the instructions to complete the **info**, **study**, **samples**, **experiments**, first **analysis** and the **dataset**. Don't forget linking the analysis to the study under the study tab. You will not need to fill in **runs** if your data are mapped BAMs.


> **Tip:** Finalise the submission early to obtain an accession number, which may be needed for a manuscript. You can continue adding analyses after finalisation. The dataset will remain private until the EGA team contacts you to confirm you are ready to make it public. You can continue adding datasets after the submission has been made public. 
 


### 3.5b Register Metadata Programmatically

Once the a batch is uploaded using globus and ingested, register metadata using the EGA Submitter Portal API. The script Register_metadata.py automates analysis registration for the esophagus study. Replace with your IDs accordingly to register your own analysis. 

```
SUBMISSION_ID    = "EGA50000001666"       # Submission containing all analyses
STUDY_ACCESSION  = "EGAS50000001793"      # Study linked to each analysis
ANALYSIS_TYPE    = "REFERENCE ALIGNMENT"  # EGA analysis type
GENOME_ID        = 15                     # GRCh38 (see /api/enums/genomes)
PLATFORM         = "UG100"                # Ultima Genomics platform
EXPERIMENT_TYPES = ["Whole genome sequencing"]

# Title pattern:    "BAM files for {sample} ({batch})"
# Description:      "Bam files resulting from merging all cram files per cell provided by Ultima Genomics"

SAMPLE_MAP = {
    "eso01": "EGAN50000419612",
    "eso02": "EGAN50000419609",
    "eso03": "EGAN50000419610",
    "eso04": "EGAN50000419611",
    "eso05": "EGAN50000419613",
}  
```

```bash
python Register_metadata.py eso02_batch3
```

The script:
- Extracts the sample name from the folder (e.g. `eso02_batch3` → `eso02`)
- Looks up the corresponding EGAN accession
- Finds all `.bam` and `.bam.bai` files in the inbox under that folder
- Creates a `REFERENCE ALIGNMENT` analysis linked to the study and sample


After creating analyses, link them to the study manually in the [Submitter Portal](https://submission.ega-archive.org/).

### 3.6 Finalise Submission 

You do not need to wait until all batches are uploaded before finalizing the submission. You can finalize as many times as needed, especially if cluster space is limited. The submission will remain in a pending/review state until the helpdesk approves it (it can take up to one week). The files will be archived and automatically deleted from the EGA box. Once you have **verified that the files have been [archived](https://submission.ega-archive.org/files/archive)**, you can safely **delete** them from the cluster to free up space before preparing the next batch.  

