# FusionRank



## Installation (single pip command)

From the repo root, install locally with dependencies and CLI entrypoint:

```bash
pip install .
```

This installs the `fusionrank` command in your environment, so you can run:

```bash
fusionrank --help
```

No separate prerequisite-install step is required.


### Download References 
To download the reference files used for feature annotation, run the following command:
```
wget https://github.com/WGLab/FusionRank/releases/download/v0.1.0/FR_refs.tar.gz
```


# Additional References: 

Three main reference files are necessary and should be downloaded to the same folder as the previous reference files: https://www.gencodegenes.org/human/

## Basic gene annotation (CHR)
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
```

## Transcript Sequences (ALL)
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
```

## Protein-coding transcript translation sequences (ALL)

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.pc_translations.fa.gz
```

---


## FusionRank Run Command

```bash
fusionrank <input_file> \
  --sample-prefix <sample_prefix> \
  --output-dir <output_base_dir> \
  --refs-dir <reference_root> \
  --model-dir <model_dir>
```

## Optional one-time PyEnsembl pre-index (recommended)

If the pipeline appears to hang during cDNA annotation (PyEnsembl indexing), run this once first:

```bash
python prepare_pyensembl_cache.py --refs-dir <reference_root>
```

Then run FusionRank normally.


## Inputs
- `input_file` (positional): fusion input TSV.
- `--sample-prefix`: prefix used for output files and run directory naming.
- `--refs-dir`: path to the reference root folder; paths are auto-discovered by expected filenames/subfolders (no manual refs JSON needed).
- `--model-dir`: directory containing the trained model files (`fusion_lambdarank.txt`, `feature_cols.pkl`, `categorical_cols.pkl`). Two models are available:
  - `FR_gene_breakpoint_model` — requires breakpoint columns (`gene1_chr`, `gene1_bkp`, `gene1_strand`, `gene2_chr`, `gene2_bkp`, `gene2_strand`)
  - `FR_gene_name_model` — works without breakpoint columns (gene name level only)

> **Note:** Breakpoint columns are optional. If absent, breakpoint annotation and cDNA features are skipped automatically and `FR_gene_name_model` should be used.


## Outputs

Created under `<output_dir>/<sample_prefix>/`:

- `<sample_prefix>_annotated.tsv`
- `<sample_prefix>_annotated_ranked.tsv`
- `<sample_prefix>_summary.tsv`
- `<sample_prefix>_top5.tsv`

Summary files are now evidence-oriented and include: sample/cancer/rank/fusion, normalized score, best breakpoints, database evidence, drug evidence, gene biology, panel evidence, and all breakpoint candidates with normalized scores.


Refs now use explicit OmniPath file keys: `ppi_adjacency_file`, `ppi_centrality_file`, and `ppi_edges_file` (auto-discovered when using `--refs-dir`).


Reference directory behavior: `--refs-dir` is the only reference input. The pipeline auto-discovers required files under that folder, creates `gencode_v49_downloads_cache/` inside it if missing, and expects these GENCODE files at the refs root: `gencode.v49.basic.annotation.gtf.gz`, `gencode.v49.transcripts.fa.gz`, `gencode.v49.pc_translations.fa.gz`.

If HPA cache files are not present yet, FusionRank can fall back to a raw HPA table if one of these files exists under `--refs-dir`: `hpa_expression.tsv`, `hpa_cancer_expression.tsv`, `HPA_cancer_expression.tsv`, `hpa_tcga_median_ptpm.tsv`, or `hpa_median_ptpm.tsv`.

HPA expression cache files can be placed directly in `--refs-dir`: `hpa_exp_dict.pkl`, `hpa_thresholds.tsv`, `hpa_low_everywhere.txt`.
