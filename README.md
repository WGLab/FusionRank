# FusionRank
FusionRank is a machine learning framework for prioritizing gene fusions detected from RNA sequencing data. FusionRank annotates the candidate fusions with biological, genomic, and structural features and ranks them according to their predicted oncogenic relevance using a Learning-to-Rank (LambdaRank) model trained on curated fusion databases. FusionRank supports fusion calls with or without breakpoint coordinate information and generates a ranked summary report for downstream review. 

<img width="1060" height="460" alt="Image" src="https://github.com/user-attachments/assets/9cb1e4d3-7e03-40bf-8eec-8274036af10b" />

## FusionRank supports two modes for prediction:
### 1. Breakpoint-Level Model
- **FR_gene_breakpoint_model**: Uses both fusion partner gene names, chromosomal breakpoint coordinates, and reported strandness.

### 2. Gene Name-Level Model
- **FR_gene_name_model**: Uses only fusion partner gene names for when breakpoint coordinates are unavailable. 

---

# Quick Start

## 1. Recommended Requirements

Create a conda environment: 
```bash
conda create -n fusionrank python=3.10.14
conda activate fusionrank
conda install -c bioconda bedtools
```

> **Note:** `bedtools` must be installed separately via conda before running FusionRank. It is required for breakpoint-based annotation. If you are using the gene name-based model (`FR_gene_name_model`) only, `bedtools` is not required.
 

## 2. FusionRank Installation (single pip command)
Clone the repository and install, from the repository root:
```bash
git clone https://github.com/WGLab/FusionRank.git
cd FusionRank
pip install .
```

Verify installation:
```bash
fusionrank --help
```
---

## 3. Download FusionRank Reference Files
```bash
wget https://github.com/WGLab/FusionRank/releases/download/v0.1.0/FR_refs.tar.gz
tar -xzf FR_refs.tar.gz
```
---

## 4. Download Required GENCODE Files
FusionRank requires the following GENCODE v49 files. Download them into the same reference directory.

### Basic Gene Annotation
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
```
### Transcript Sequences
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
```
### Protein-Coding Translation Sequences
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.pc_translations.fa.gz
```
---

# 5. Running FusionRank
```bash
fusionrank <input_file> \
    --sample-prefix <sample_prefix> \
    --output-dir <output_dir> \
    --refs-dir <reference_dir> \
    --model-dir <model_dir>
```

## Required Parameters

| Parameter | Description |
|:------------:|:------------:|
| `input_file` | Fusion calls as a TSV file (see below for input file format) |
| `--sample-prefix` | Prefix used for output files |
| `--output-dir` | Directory for Output folder with sample-prefix |
| `--refs-dir` | Reference directory |
| `--model-dir` | Directory containing trained model files (see below) |

## Model Selection
FusionRank provides two pretrained models. Select the appropriate model by specifying the path with `--model-dir`.

### Breakpoint-Level Model
Use when chromosomal breakpoint coordinates are available:
```text
--model-dir /path/to/FR_gene_breakpoint_model
```

### Gene-Level Model
Use when only fusion partner gene names are available:
```text
--model-dir /path/to/FR_gene_name_model
```

> **Note:** Replace the example paths above with the location of the downloaded model directory on your system. FusionRank automatically detects whether breakpoint columns are present in the input file. If breakpoint information is unavailable, breakpoint-dependent annotations are skipped and the gene-level model should be specified.


---

## Input Format
Input files must be tab-separated (TSV) files containing one fusion candidate per row.

### Required Columns

| sample_id | cancer_type | gene1 | gene2 |
|:------------:|:------------:|:--------:|:--------:|
| K562_rep1 | LAML | BCR | ABL1 |
| Sample_02 | LUAD | EML4 | ALK |

**Column descriptions**

| Column | Description |
|:----------:|:----------:|
| `sample_id` | Unique sample identifier |
| `cancer_type` | TCGA or pediatric cancer acronym |
| `gene1` | 5' fusion partner gene |
| `gene2` | 3' fusion partner gene |

> **Note:** FusionRank supports processing multiple samples within a single input file. Fusion candidates are automatically grouped by `sample_id`, and ranking is performed independently within each sample.

---

### Optional Additional Breakpoint (hg38) Columns
The following columns are required only when using the breakpoint-level model (`FR_gene_breakpoint_model`).

| sample_id | cancer_type | gene1 | gene2 | gene1_chr | gene1_bkp | gene1_strand | gene2_chr | gene2_bkp | gene2_strand |
|:------------:|:------------:|:--------:|:--------:|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|
| K562_rep1 | LAML | BCR | ABL1 | chr22 | 23290410 | + | chr9 | 130854063 | + |

**Breakpoint column descriptions**

| Column | Description |
|:----------:|:----------:|
| `sample_id` | Unique sample identifier |
| `cancer_type` | TCGA or pediatric cancer acronym |
| `gene1_chr` | Chromosome of the 5' fusion breakpoint |
| `gene1_bkp` | Genomic coordinate of the 5' breakpoint |
| `gene1_strand` | Strand of the 5' breakpoint (`+` or `-`) |
| `gene2_chr` | Chromosome of the 3' fusion breakpoint |
| `gene2_bkp` | Genomic coordinate of the 3' breakpoint |
| `gene2_strand` | Strand of the 3' breakpoint (`+` or `-`) |

> **Note:** If breakpoint columns are omitted, FusionRank automatically skips breakpoint-dependent annotations and should be run using the gene-level model (`FR_gene_name_model`).

---

# Supported Cancer Types

The `cancer_type` column must contain one of the following cancer acronyms:

| Cancer Acronym | Full Name |
|:---------:|:-----------:|
| ACC | Adrenocortical Carcinoma |
| ALL | Acute Lymphoblastic Leukemia |
| BLCA | Bladder Urothelial Carcinoma |
| BRCA | Breast Invasive Carcinoma |
| CCSK | Clear Cell Sarcoma of the Kidney |
| CESC | Cervical Squamous Cell Carcinoma |
| CHOL | Cholangiocarcinoma |
| COAD | Colon Adenocarcinoma |
| DLBC | Diffuse Large B-Cell Lymphoma |
| ESCA | Esophageal Carcinoma |
| GBM | Glioblastoma Multiforme |
| HNSC | Head and Neck Squamous Cell Carcinoma |
| KICH | Kidney Chromophobe |
| KIRC | Kidney Renal Clear Cell Carcinoma |
| KIRP | Kidney Renal Papillary Cell Carcinoma |
| LAML | Acute Myeloid Leukemia |
| LGG | Brain Lower Grade Glioma |
| LIHC | Liver Hepatocellular Carcinoma |
| LUAD | Lung Adenocarcinoma |
| LUSC | Lung Squamous Cell Carcinoma |
| MESO | Mesothelioma |
| NBL | Neuroblastoma |
| OS | Osteosarcoma |
| OV | Ovarian Serous Cystadenocarcinoma |
| PAAD | Pancreatic Adenocarcinoma |
| PCPG | Pheochromocytoma and Paraganglioma |
| PRAD | Prostate Adenocarcinoma |
| READ | Rectum Adenocarcinoma |
| RT | Rhabdoid Tumor |
| SARC | Sarcoma |
| SKCM | Skin Cutaneous Melanoma |
| STAD | Stomach Adenocarcinoma |
| TGCT | Testicular Germ Cell Tumors |
| THCA | Thyroid Carcinoma |
| THYM | Thymoma |
| UCS | Uterine Carcinosarcoma |
| UCEC | Uterine Corpus Endometrial Carcinoma |
| UVM | Uveal Melanoma |
| WT | Wilms Tumor |

---

# Output Files
Results are written to:
```text
<output_dir>/<sample_prefix>/
```

**Generated files include:**
| File | Description |
|:-----:|:------------:|
| `<sample_prefix>_annotated.tsv` | Full feature matrix with all annotation columns |
| `<sample_prefix>_annotated_ranked.tsv` | Annotated fusions with annoations, predicted scores, and per-sample ranks |
| `<sample_prefix>_summary.tsv` | Evidence-oriented summary with database, drug, biology, and panel annotations |
| `<sample_prefix>_summary.html` | Browser-viewable ranked fusion report with top 5 candidates highlighted |

---

# Processing Notes:

## Optional one-time PyEnsembl pre-index (recommended)
If the pipeline appears to stall during cDNA annotation (PyEnsembl indexing), run this once first:
```bash
python prepare_pyensembl_cache.py --refs-dir <reference_root>
```
Then run FusionRank normally.

--- 
Refs now use explicit OmniPath file keys: `ppi_adjacency_file`, `ppi_centrality_file`, and `ppi_edges_file` (auto-discovered when using `--refs-dir`).

Reference directory behavior: `--refs-dir` is the only reference input. The pipeline auto-discovers required files under that folder, creates `gencode_v49_downloads_cache/` inside it if missing, and expects these GENCODE files at the refs root: `gencode.v49.basic.annotation.gtf.gz`, `gencode.v49.transcripts.fa.gz`, `gencode.v49.pc_translations.fa.gz`.

If HPA cache files are not present yet, FusionRank can fall back to a raw HPA table if one of these files exists under `--refs-dir`: `hpa_expression.tsv`, `hpa_cancer_expression.tsv`, `HPA_cancer_expression.tsv`, `hpa_tcga_median_ptpm.tsv`, or `hpa_median_ptpm.tsv`.

HPA expression cache files can be placed directly in `--refs-dir`: `hpa_exp_dict.pkl`, `hpa_thresholds.tsv`, `hpa_low_everywhere.txt`.
