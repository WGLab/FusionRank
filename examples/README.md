# Example Dataset

This directory contains example input files and expected outputs for FusionRank using K562 chronic myelogenous leukemia (CML) cell line fusions.


## Example Run Commands


### Run with breakpoint coordinates
```bash
fusionrank /path/to/examples/input/K562_input_fusions_breakpoints.tsv \
  --sample-prefix K562_breakpoints \
  --output-dir /path/to/examples/expected_output \
  --refs-dir /path/to/FR_refs \
  --model-dir /path/to/FR_gene_breakpoint_model
```

### Run without breakpoint coordinates
```bash
fusionrank /path/to/examples/input/K562_input_fusions_gene_names.tsv \
  --sample-prefix K562_gene_names \
  --output-dir /path/to/examples/output \
  --refs-dir /path/to/FR_refs \
  --model-dir /path/to/FR_gene_name_model
```

Expected outputs are provided in `expected_output/` for reference.



## Output: `_summary.html` Overview





An interactive browser-viewable report with the following columns:

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| Cancer | TCGA cancer type |
| Rank | Within-sample rank |
| Fusion | Gene fusion name |
| Score (0-1) | Normalized relevance score |
| Gene1 Breakpoint | Gene 1 chromosomal breakpoint coordinate |
| Gene2 Breakpoint | Gene 2 chromosomal breakpoint coordinate |
| Database Evidence | CIViC, COSMIC, Mitelman Cancer-specific matches |
| Drug Evidence | Actionable therapeutic target status from ChiTaRs |
| Gene Biology | Oncogene/TSG/Kinase roles |
| Panel Evidence | Gene-level clinical fusion panel membership |
| All Breakpoints (Norm. Score) | All fusion breakpoints for the given fusion and its normalized relevance score |

