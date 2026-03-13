#!/usr/bin/env python3

from pathlib import Path

REQUIRED_REF_KEYS = [
    'gtex_file',
    'babiceanu_file',
    'ariba_file',
    'gene_to_type_pkl',
    'civic_file',
    'cosmic_file',
    'mitelman_file',
    'cancermine_file',
    'hugo_gene_sets',
    'fusion_panels',
    'tad_pickle',
    'ppi_adjacency_file',
    'ppi_centrality_file',
    'ppi_edges_file',
    'druggable_file',
    'gencode_loci_bed',
    'gencode_bkp_bed',
    'gencode_gtf',
    'gencode_tx_fasta',
    'gencode_prot_fasta',
    'gencode_cache_dir',
]

REF_FILE_BASENAMES = {
    'gtex_file': 'GTEx_gencode_names_updated.tsv',
    'babiceanu_file': 'fusions_Babs_normal.tsv',
    'ariba_file': 'blacklist_genes_processed.tsv',
    'gene_to_type_pkl': 'gencode_gene_to_type.pkl',
    'civic_file': 'civicDB_fusions_processed.tsv',
    'cosmic_file': 'cosmic_fusions_to_TCGA_cancers.tsv',
    'mitelman_file': 'mitelman_db_cleaned_TCGA_cancer_types.tsv',
    'cancermine_file': 'processed_cancermine_genes.tsv',
    'hugo_gene_sets': 'hgnc_gene_type_sets.json',
    'fusion_panels': 'clinical_fusion_panels_genes.tsv',
    'tad_pickle': 'TAD_gene_lists.pkl',
    'druggable_file': 'druggable_fusions_all.tsv',
    'ppi_adjacency_file': 'ppi_adjacency.pkl',
    'ppi_centrality_file': 'ppi_centrality.tsv',
    'ppi_edges_file': 'ppi_edges.tsv',
    'gencode_loci_bed': 'gencode_v49_gene_loci.sorted.bed',
    'gencode_bkp_bed': 'gencode_v49_BED_format.bed',
    'gencode_gtf': 'gencode.v49.basic.annotation.gtf.gz',
    'gencode_tx_fasta': 'gencode.v49.transcripts.fa.gz',
    'gencode_prot_fasta': 'gencode.v49.pc_translations.fa.gz',
}

GTF_ROOT_KEYS = {'gencode_gtf', 'gencode_tx_fasta', 'gencode_prot_fasta'}


OPTIONAL_REF_FILE_BASENAMES = {
    # Optional precomputed HPA cache files (can live directly under --refs-dir).
    'hpa_exp_dict_file': ['hpa_exp_dict.pkl'],
    'hpa_thresholds_file': ['hpa_thresholds.tsv'],
    'hpa_low_everywhere_file': ['hpa_low_everywhere.txt'],

    # Used only if HPA cache files are missing.
    'hpa_file': [
        'hpa_expression.tsv',
        'hpa_cancer_expression.tsv',
        'HPA_cancer_expression.tsv',
        'hpa_tcga_median_ptpm.tsv',
        'hpa_median_ptpm.tsv',
    ],
}


def _require_ref_keys(refs):
    missing = [k for k in REQUIRED_REF_KEYS if k not in refs or not str(refs[k]).strip()]
    if missing:
        raise ValueError(f'Missing refs keys: {", ".join(missing)}')


def _find_unique(root, name, kind='file'):
    if kind == 'file':
        matches = [p for p in root.rglob(name) if p.is_file()]
    else:
        matches = [p for p in root.rglob(name) if p.is_dir()]

    if not matches:
        raise FileNotFoundError(f"Could not find {kind} '{name}' under {root}")
    if len(matches) > 1:
        raise ValueError(
            f"Found multiple {kind}s named '{name}' under {root}; please make filenames unique under --refs-dir. "
            f"Matches: {', '.join(str(m) for m in matches[:5])}"
        )
    return str(matches[0])



def _find_optional_by_basenames(root, basenames):
    hits = []
    for basename in basenames:
        hits.extend([p for p in root.rglob(basename) if p.is_file()])
    # de-duplicate paths while preserving order
    unique = []
    seen = set()
    for p in hits:
        sp = str(p)
        if sp not in seen:
            seen.add(sp)
            unique.append(sp)
    if len(unique) > 1:
        raise ValueError(
            'Found multiple possible HPA raw expression files under refs dir; '
            'please keep one of the supported names. Matches: ' + ', '.join(unique[:5])
        )
    return unique[0] if unique else None

def load_refs_from_dir(refs_dir):
    root = Path(refs_dir)
    if not root.exists() or not root.is_dir():
        raise NotADirectoryError(f'Reference directory not found: {refs_dir}')

    refs = {}
    for key, basename in REF_FILE_BASENAMES.items():
        if key in GTF_ROOT_KEYS:
            path = root / basename
            if not path.exists() or not path.is_file():
                raise FileNotFoundError(
                    f"Expected {basename} at reference root: {path}. "
                    "Place GENCODE v49 files directly in --refs-dir."
                )
            refs[key] = str(path)
        else:
            refs[key] = _find_unique(root, basename, kind='file')

    # Create required gencode cache directory inside the provided refs root.
    gencode_cache_dir = root / 'gencode_v49_downloads_cache'
    gencode_cache_dir.mkdir(parents=True, exist_ok=True)
    refs['gencode_cache_dir'] = str(gencode_cache_dir)

    hpa_exp_dict_file = _find_optional_by_basenames(root, OPTIONAL_REF_FILE_BASENAMES['hpa_exp_dict_file'])
    hpa_thresholds_file = _find_optional_by_basenames(root, OPTIONAL_REF_FILE_BASENAMES['hpa_thresholds_file'])
    hpa_low_everywhere_file = _find_optional_by_basenames(root, OPTIONAL_REF_FILE_BASENAMES['hpa_low_everywhere_file'])
    if hpa_exp_dict_file:
        refs['hpa_exp_dict_file'] = hpa_exp_dict_file
    if hpa_thresholds_file:
        refs['hpa_thresholds_file'] = hpa_thresholds_file
    if hpa_low_everywhere_file:
        refs['hpa_low_everywhere_file'] = hpa_low_everywhere_file

    hpa_file = _find_optional_by_basenames(root, OPTIONAL_REF_FILE_BASENAMES['hpa_file'])
    if hpa_file:
        refs['hpa_file'] = hpa_file



    _require_ref_keys(refs)
    return refs
