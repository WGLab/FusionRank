#!/usr/bin/env python3

import argparse
from pathlib import Path

from fusionrank_config import load_refs_from_dir


def parse_args():
    parser = argparse.ArgumentParser(
        description='Pre-build PyEnsembl (GENCODE v49) cache used by FusionRank cDNA annotation.'
    )
    parser.add_argument(
        '--refs-dir',
        required=True,
        help='Reference root folder (same value used for fusionrank.py --refs-dir).',
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print('Resolving reference paths...')
    refs = load_refs_from_dir(args.refs_dir)
    refs_dir = Path(args.refs_dir)

    cache_dir = refs_dir / 'gencode_v49_downloads_cache'
    cache_dir.mkdir(parents=True, exist_ok=True)
    print(f'Cache directory: {cache_dir}')

    print('\nBuilding PyEnsembl cache — this may take several minutes...')
    import pyensembl

    gd = pyensembl.Genome(
        reference_name='GRCh38',
        annotation_name='gencode',
        annotation_version='v49',
        gtf_path_or_url=str(refs['gencode_gtf']),
        transcript_fasta_paths_or_urls=str(refs['gencode_tx_fasta']),
        protein_fasta_paths_or_urls=str(refs['gencode_prot_fasta']),
        cache_directory_path=str(cache_dir),
        copy_local_files_to_cache=False,
    )
    gd.index()

    print('✅ PyEnsembl cache built successfully!')
    print(f'   Cache location: {cache_dir}')
    print('\nYou can now run FusionRank without the indexing delay.')


if __name__ == '__main__':
    main()
