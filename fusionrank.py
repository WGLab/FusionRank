#!/usr/bin/env python3

import argparse
from pathlib import Path

from fusionrank_config import load_refs_from_dir


def parse_args():
    parser = argparse.ArgumentParser(description='FusionRank annotation + ranking pipeline')
    parser.add_argument('input_file', help='Input fusion TSV file')
    parser.add_argument('--sample-prefix', required=True, help='Sample prefix used for output files')
    parser.add_argument('--output-dir', required=True, help='Base output directory')
    parser.add_argument('--refs-dir', required=True, help='Path to reference root folder (auto-resolve required files)')
    parser.add_argument('--model-dir', required=True, help='Directory containing trained ranking model files')
    return parser.parse_args()


def main():
    args = parse_args()

    from fusionrank_core import load_model, make_summary_tsv, rank_fusions, run_full_pipeline

    refs = load_refs_from_dir(args.refs_dir)

    run_output_dir = Path(args.output_dir) / args.sample_prefix
    run_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"{'='*60}")
    print('FUSION ANNOTATION PIPELINE')
    print(f"{'='*60}")
    print(f'Sample:     {args.sample_prefix}')
    print(f'Input:      {args.input_file}')
    print(f'Output dir: {run_output_dir}')
    print(f"{'='*60}\n")

    df = run_full_pipeline(
        input_file=args.input_file,
        output_dir=str(run_output_dir),
        project_name=args.sample_prefix,
        refs=refs,
    )

    print('\n[16] Loading model and ranking...')
    model, feature_cols, categorical_cols = load_model(args.model_dir)
    df = rank_fusions(df, model, feature_cols, categorical_cols)

    print('\n[17] Saving outputs...')
    full_out = run_output_dir / f'{args.sample_prefix}_annotated_ranked.tsv'
    df.to_csv(full_out, sep='\t', index=False)
    print(f'    Full:    {full_out}')

    summary_df = make_summary_tsv(df, refs)
    summary_out = run_output_dir / f'{args.sample_prefix}_summary.tsv'
    summary_df.to_csv(summary_out, sep='\t', index=False)
    print(f'    Summary: {summary_out}')

    top5_out = run_output_dir / f'{args.sample_prefix}_top5.tsv'
    rank_col = 'Rank' if 'Rank' in summary_df.columns else 'rank'
    summary_df[summary_df[rank_col] <= 5].to_csv(top5_out, sep='\t', index=False)
    print(f'    Top 5:   {top5_out}')

    print(f'\n{"="*60}')
    print(f'✅  Done. Final shape: {df.shape}')
    print(f'    Samples:  {df["sample_id"].nunique():,}')
    print(f'    Fusions:  {df["fusion_name"].nunique():,}')
    print(f'    Outputs:  {run_output_dir}')
    print(f'{"="*60}')


if __name__ == '__main__':
    main()
