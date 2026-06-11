#!/usr/bin/env python3

import argparse
from pathlib import Path

from fusionrank_config import load_refs_from_dir


def move_bed_outputs(output_dir):
    """Move BED artifacts into an output subdirectory after pipeline completion."""
    output_dir = Path(output_dir)
    bed_dir = output_dir / 'bed_files'
    moved = []

    for bed_file in sorted(output_dir.glob('*.bed')):
        if not bed_dir.exists():
            bed_dir.mkdir(parents=True, exist_ok=True)

        destination = bed_dir / bed_file.name
        if destination.exists():
            destination.unlink()
        bed_file.replace(destination)
        moved.append(destination)

    return moved


def parse_args():
    parser = argparse.ArgumentParser(description='FusionRank annotation + ranking pipeline')
    parser.add_argument('input_file', help='Input fusion TSV file')
    parser.add_argument('--sample-prefix', required=True, help='Sample prefix used for output files')
    parser.add_argument('--output-dir', required=True, help='Base output directory')
    parser.add_argument('--refs-dir', required=True, help='Path to reference root folder (auto-resolve required files)')
    parser.add_argument('--model-dir', required=True, help='Directory containing trained model files. Use FR_gene_breakpoint_model (requires breakpoint columns) or FR_gene_name_model (gene name level, no breakpoints required)')
    return parser.parse_args()


def main():
    args = parse_args()

    from fusionrank_core import (
        load_model,
        make_ranked_input_tsv,
        make_summary_tsv,
        rank_fusions,
        run_full_pipeline,
        write_summary_html,
    )

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

    print('\nLoading model and ranking...')
    model, feature_cols, categorical_cols = load_model(args.model_dir)
    df = rank_fusions(df, model, feature_cols, categorical_cols)

    print('\nSaving outputs...')
    ranked_out = run_output_dir / f'{args.sample_prefix}_ranked.tsv'
    ranked_df = make_ranked_input_tsv(df, args.input_file)
    ranked_df.to_csv(ranked_out, sep='\t', index=False)
    
    print(f'Ranked Fusions: {ranked_out}')
    legacy_annotated_ranked_out = run_output_dir / f'{args.sample_prefix}_annotated_ranked.tsv'
    if legacy_annotated_ranked_out.exists():
        legacy_annotated_ranked_out.unlink()
    
    summary_df = make_summary_tsv(df, refs)
    summary_out = run_output_dir / f'{args.sample_prefix}_summary.tsv'
    summary_df.to_csv(summary_out, sep='\t', index=False)

    html_out = run_output_dir / f'{args.sample_prefix}_summary.html'
    write_summary_html(summary_df, html_out, args.sample_prefix)

    if not html_out.exists():
        raise RuntimeError(f'HTML summary was not created: {html_out}')
    print(f'Annotation Summary: {html_out}')

    legacy_top5_out = run_output_dir / f'{args.sample_prefix}_top5.tsv'
    if legacy_top5_out.exists():
        legacy_top5_out.unlink()

    moved_bed_files = move_bed_outputs(run_output_dir)
    if moved_bed_files:
        print(f"BED Files: {run_output_dir / 'bed_files'}")

    print(f'\n{"="*60}')
    print(f'✅  Done. Final shape: {df.shape}')
    print(f'    Samples:  {df["sample_id"].nunique():,}')
    print(f'    Fusions:  {df["fusion_name"].nunique():,}')
    print(f'    Outputs:  {run_output_dir}')
    print(f'{"="*60}')


if __name__ == '__main__':
    main()
