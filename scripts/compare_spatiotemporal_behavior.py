#!/usr/bin/env python3
"""Compare spatiotemporal propagation summaries across groups.

This script builds on ``spatiotemporal_signal_propagation.py``.
Instead of analysing one experiment-site block, it loops over many blocks,
computes the same local-neighbour / near-future jump summary, and compares
those summaries across higher-level groups.

In principle, a natural default would be experimental conditions. In this
specific dataset, however, ``Conditions`` is constant across all sites, so the
script automatically falls back to ``Mutation`` as the default grouping.

Useful questions supported by this workflow include:

- Do PI3K-pathway mutants differ from WT in local spatiotemporal coordination?
- Is neighbour-linked propagation stronger for ERK than for FoxO?
- Are some experiment replicates systematically more coordinated than others?
- Does the result change when we vary the spatial radius or future window?
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

try:
    from spatiotemporal_signal_propagation import (
        add_track_deltas,
        annotate_spatial_exposure,
        assign_jump_events,
        build_spatial_edges,
        build_temporal_edges,
        compute_future_jump_flags,
        load_metadata,
        load_site_block,
        resolve_path,
        summarise_propagation,
    )
except ModuleNotFoundError:
    from scripts.spatiotemporal_signal_propagation import (
        add_track_deltas,
        annotate_spatial_exposure,
        assign_jump_events,
        build_spatial_edges,
        build_temporal_edges,
        compute_future_jump_flags,
        load_metadata,
        load_site_block,
        resolve_path,
        summarise_propagation,
    )


GROUP_CHOICES = ['auto', 'conditions', 'mutation', 'exp_id', 'site_id']


def parse_args() -> argparse.Namespace:
    """Read command-line options for the comparison workflow.

    This function lets us choose the signal, grouping variable, and
    analysis settings used when comparing many experiment-site blocks.
    """
    parser = argparse.ArgumentParser(
        description='Compare spatiotemporal propagation summaries across groups.'
    )
    parser.add_argument('--data-path', type=Path, default=Path('single-cell-tracks_exp1-6_noErbB2.csv.gz'))
    parser.add_argument('--meta-path', type=Path, default=Path('01-readme-experiment-description_2022-04-05.csv'))
    parser.add_argument(
        '--signal-col',
        type=str,
        default='ERKKTR_ratio',
        choices=['ERKKTR_ratio', 'FoxO3A_ratio'],
        help='Signal used to define jump events and propagation summaries.',
    )
    parser.add_argument(
        '--group-by',
        type=str,
        default='auto',
        choices=GROUP_CHOICES,
        help='How to group experiment-site blocks for comparison. Auto prefers Conditions when available, else Mutation.',
    )
    parser.add_argument(
        '--spatial-radius',
        type=float,
        default=60.0,
        help='Maximum Euclidean distance in image coordinates to connect spatial neighbours.',
    )
    parser.add_argument(
        '--future-window-frames',
        type=int,
        default=3,
        help='How many future frames count as near future when evaluating self jumps.',
    )
    parser.add_argument(
        '--jump-threshold',
        type=float,
        default=None,
        help='Absolute threshold on positive signal difference to call a jump. If omitted, a per-block quantile threshold is used.',
    )
    parser.add_argument(
        '--jump-quantile',
        type=float,
        default=0.9,
        help='Quantile of positive signal differences used when --jump-threshold is not provided.',
    )
    parser.add_argument('--chunksize', type=int, default=1_000_000)
    parser.add_argument(
        '--exclude-mutations',
        nargs='*',
        default=['ErbB2'],
        help='Mutations to exclude from the comparison table.',
    )
    parser.add_argument(
        '--exp-ids',
        nargs='*',
        type=int,
        default=None,
        help='Optional subset of experiments to analyse.',
    )
    parser.add_argument(
        '--site-ids',
        nargs='*',
        type=int,
        default=None,
        help='Optional subset of site ids to analyse.',
    )
    parser.add_argument(
        '--max-blocks',
        type=int,
        default=None,
        help='Optional cap on the number of experiment-site blocks, useful for testing.',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('analysis_outputs'),
        help='Directory where comparison tables and summaries are written.',
    )
    return parser.parse_args()


def find_available_blocks(data_path: Path, chunksize: int) -> pd.DataFrame:
    """Find which experiment-site blocks are actually present in the CSV file.

    This is useful because later comparisons should only iterate over blocks
    that really exist in the exported dataset.
    """
    pairs: set[tuple[int, int]] = set()
    usecols = ['Exp_ID', 'Image_Metadata_Site']
    for chunk in pd.read_csv(data_path, usecols=usecols, chunksize=chunksize):
        pairs.update((int(exp_id), int(site_id)) for exp_id, site_id in chunk.drop_duplicates().itertuples(index=False))

    pairs_df = pd.DataFrame(sorted(pairs), columns=['Exp_ID', 'Image_Metadata_Site'])
    return pairs_df


def choose_grouping(meta: pd.DataFrame, requested: str) -> tuple[str, str]:
    """Decide which metadata field will define the comparison groups.

    In this dataset, `Conditions` does not vary, so the automatic choice falls
    back to `Mutation`, which is the most informative default.
    """
    if requested == 'auto':
        if meta['Conditions'].nunique(dropna=False) > 1:
            return 'Conditions', 'conditions'
        return 'Mutation', 'mutation'

    mapping = {
        'conditions': ('Conditions', 'conditions'),
        'mutation': ('Mutation', 'mutation'),
        'exp_id': ('Exp_ID', 'exp_id'),
        'site_id': ('Image_Metadata_Site', 'site_id'),
    }
    return mapping[requested]


def prepare_block_table(meta: pd.DataFrame, available_blocks: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Build the table of blocks that will be analysed and compared.

    This step merges available blocks with metadata and applies any user
    filters, such as chosen experiments, sites, or mutations to exclude.
    """
    meta = meta.copy()
    block_table = available_blocks.merge(meta, on='Image_Metadata_Site', how='left')

    if args.exclude_mutations:
        block_table = block_table.loc[~block_table['Mutation'].isin(args.exclude_mutations)].copy()
    if args.exp_ids:
        block_table = block_table.loc[block_table['Exp_ID'].isin(args.exp_ids)].copy()
    if args.site_ids:
        block_table = block_table.loc[block_table['Image_Metadata_Site'].isin(args.site_ids)].copy()

    block_table = block_table.sort_values(['Exp_ID', 'Image_Metadata_Site']).reset_index(drop=True)
    if args.max_blocks is not None:
        block_table = block_table.head(args.max_blocks).copy()
    return block_table


def run_block_analysis(block_info: pd.Series, data_path: Path, frame_to_minutes: float, args: argparse.Namespace) -> dict:
    """Run the single-block propagation analysis for one row of the block table.

    This is the bridge between the first script and the comparison workflow:
    one block is analysed in detail, and only its summary numbers are returned.
    """
    exp_id = int(block_info['Exp_ID'])
    site_id = int(block_info['Image_Metadata_Site'])
    mutation = block_info.get('Mutation')

    block = load_site_block(data_path, exp_id, site_id, args.chunksize)
    block, default_threshold = add_track_deltas(block, args.signal_col, frame_to_minutes, args.jump_quantile)
    threshold = float(args.jump_threshold) if args.jump_threshold is not None else default_threshold
    block = assign_jump_events(block, threshold)

    spatial_edges = build_spatial_edges(block, args.spatial_radius)
    temporal_edges = build_temporal_edges(block, frame_to_minutes)
    block = compute_future_jump_flags(block, args.future_window_frames)
    block = annotate_spatial_exposure(block, spatial_edges)

    single_block_args = argparse.Namespace(
        exp_id=exp_id,
        site_id=site_id,
        signal_col=args.signal_col,
        spatial_radius=args.spatial_radius,
        future_window_frames=args.future_window_frames,
    )
    summary = summarise_propagation(block, spatial_edges, temporal_edges, threshold, single_block_args, mutation)
    summary['Conditions'] = block_info.get('Conditions')
    summary['group_site_id'] = site_id
    return summary


def aggregate_comparison(block_summaries: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """Aggregate many block-level summaries into one table per comparison group.

    For example, if the grouping is `Mutation`, this function combines all
    experiment-site blocks belonging to each mutation and computes average
    propagation metrics for that mutation.

    In the output table:

    - ``mean_exposed_jump_rate`` is the average probability of a cell making
      its own large upward step soon after at least one nearby cell does so.
    - ``mean_unexposed_jump_rate`` is the same probability when no nearby cell
      is making such a step.
    - ``mean_risk_difference`` is the absolute gap between those probabilities.
    - ``mean_relative_risk`` is the fold change between those probabilities.
    """
    aggregated = (
        block_summaries
        .groupby(group_col, dropna=False)
        .agg(
            n_blocks=('site_id', 'size'),
            n_unique_sites=('site_id', 'nunique'),
            n_unique_experiments=('exp_id', 'nunique'),
            total_nodes=('n_nodes', 'sum'),
            mean_relative_risk=('relative_risk', 'mean'),
            median_relative_risk=('relative_risk', 'median'),
            mean_risk_difference=('risk_difference', 'mean'),
            mean_exposed_jump_rate=('future_jump_rate_if_neighbor_jumps_now', 'mean'),
            mean_unexposed_jump_rate=('future_jump_rate_if_no_neighbor_jumps_now', 'mean'),
            mean_spatial_edges=('n_spatial_edges', 'mean'),
            mean_temporal_edges=('n_temporal_edges', 'mean'),
        )
        .reset_index()
        .sort_values('mean_relative_risk', ascending=False)
    )
    return aggregated


def main() -> None:
    """Run the full multi-block comparison workflow.

    This function finds blocks, chooses the grouping strategy, analyses each
    block one by one, aggregates the results, and saves comparison tables that
    can be explored later in notebooks.
    """
    args = parse_args()
    data_path = resolve_path(args.data_path)
    meta_path = resolve_path(args.meta_path)
    output_root = resolve_path(args.output_dir)

    meta = load_metadata(meta_path)
    available_blocks = find_available_blocks(data_path, args.chunksize)
    block_table = prepare_block_table(meta, available_blocks, args)
    if block_table.empty:
        raise ValueError('No experiment-site blocks remain after filtering.')

    group_col, group_label = choose_grouping(meta, args.group_by)
    frame_to_minutes = float(meta['Acquisition_frequency_min'].iloc[0])

    block_summaries = []
    for _, row in block_table.iterrows():
        summary = run_block_analysis(row, data_path, frame_to_minutes, args)
        summary['comparison_group'] = row[group_col] if group_col in row else row.get(group_col)
        block_summaries.append(summary)

    block_df = pd.DataFrame(block_summaries)
    comparison_df = aggregate_comparison(block_df, 'comparison_group')

    output_dir = output_root / f'comparison_{group_label}_{args.signal_col}'
    output_dir.mkdir(parents=True, exist_ok=True)
    block_df.to_csv(output_dir / 'block_level_summary.csv', index=False)
    comparison_df.to_csv(output_dir / 'group_level_summary.csv', index=False)

    task_description = {
        'comparison_question': (
            'Are local neighbour-linked signal jumps associated with stronger near-future jumps '
            'in some groups than in others?'
        ),
        'default_grouping_reason': (
            'Conditions does not vary in this dataset, so Mutation is the most informative default grouping.'
            if group_label == 'mutation' and meta['Conditions'].nunique(dropna=False) == 1
            else f'Using {group_col} as the comparison grouping.'
        ),
        'group_by': group_col,
        'signal_col': args.signal_col,
        'spatial_radius': args.spatial_radius,
        'future_window_frames': args.future_window_frames,
        'n_blocks_analysed': int(len(block_df)),
        'groups_found': comparison_df['comparison_group'].astype(str).tolist(),
    }
    with (output_dir / 'task_description.json').open('w', encoding='utf-8') as handle:
        json.dump(task_description, handle, indent=2)

    print('\nComparison task complete.')
    print(f'Grouping variable: {group_col}')
    print(f'Blocks analysed: {len(block_df)}')
    print(f'Output directory: {output_dir}')
    print('Interpretation of "jump": a large upward change since the previous frame.')
    print('mean_risk_difference = absolute probability gap between neighbour-exposed and unexposed cell-time points.')
    print('mean_relative_risk = how many times more likely a future jump is when a nearby cell jumps now.')
    print('\nTop grouped summary rows:')
    print(comparison_df.head().to_string(index=False))


if __name__ == '__main__':
    main()
