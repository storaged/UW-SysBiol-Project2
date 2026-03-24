#!/usr/bin/env python3
"""Build a spatiotemporal graph for one experiment-site block and quantify
whether local neighbour activity predicts near-future jumps in a target cell.

The script treats each row in the trajectory table as one node in a graph.
It creates two edge types:

1. spatial edges: cells close to one another at the same frame,
2. temporal edges: consecutive observations of the same track over time.

This allows a first-pass propagation analysis addressing questions like:
"If my neighbours jump in ERK activity now, am I more likely to jump soon after?"

The graph is exported as node and edge tables so it can be reused later in
network analysis, signal propagation modelling, or custom downstream scripts.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

PROJECT_ROOT = Path(__file__).resolve().parent.parent


REQUIRED_COLUMNS = [
    'Exp_ID',
    'Image_Metadata_Site',
    'track_id',
    'Image_Metadata_T',
    'ERKKTR_ratio',
    'FoxO3A_ratio',
    'Nuclear_size',
    'objNuclei_Location_Center_X',
    'objNuclei_Location_Center_Y',
]


def parse_args() -> argparse.Namespace:
    """Read command-line options.

    This is the function that lets us choose which experiment-site block,
    signal, radius, and time window should be analysed.
    """
    parser = argparse.ArgumentParser(
        description='Build a spatiotemporal neighbour graph and run a simple propagation analysis.'
    )
    parser.add_argument('--data-path', type=Path, default=Path('single-cell-tracks_exp1-6_noErbB2.csv.gz'))
    parser.add_argument('--meta-path', type=Path, default=Path('01-readme-experiment-description_2022-04-05.csv'))
    parser.add_argument('--exp-id', type=int, required=True, help='Experiment identifier to analyse.')
    parser.add_argument('--site-id', type=int, required=True, help='Site identifier to analyse.')
    parser.add_argument(
        '--signal-col',
        type=str,
        default='ERKKTR_ratio',
        choices=['ERKKTR_ratio', 'FoxO3A_ratio'],
        help='Signal used to define jump events and propagation summaries.',
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
        help='How many future frames count as the near future for the target-cell jump analysis.',
    )
    parser.add_argument(
        '--jump-threshold',
        type=float,
        default=None,
        help='Absolute threshold on positive signal difference to call a jump. If omitted, a data-driven threshold is used.',
    )
    parser.add_argument(
        '--jump-quantile',
        type=float,
        default=0.9,
        help='Quantile of positive signal differences used when --jump-threshold is not provided.',
    )
    parser.add_argument(
        '--chunksize',
        type=int,
        default=1_000_000,
        help='Chunk size used when scanning the large CSV file.',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('analysis_outputs'),
        help='Directory where node/edge tables and the summary JSON are written.',
    )
    return parser.parse_args()


def resolve_path(path: Path) -> Path:
    """Convert a relative path into an absolute path inside the project.

    Prefer the current working directory when the relative path exists there,
    but fall back to the project root so the scripts still work when launched
    from inside ``scripts/``.
    """
    if path.is_absolute():
        return path

    cwd_path = Path.cwd() / path
    if cwd_path.exists():
        return cwd_path

    return PROJECT_ROOT / path


def load_metadata(meta_path: Path) -> pd.DataFrame:
    """Load the experiment-description table and standardize the site column.

    The metadata table tells us which mutation belongs to which imaging site,
    so it provides the biological context for the graph analysis.
    """
    meta = pd.read_csv(meta_path, encoding='utf-8-sig')
    meta = meta.rename(columns={'Site': 'Image_Metadata_Site'})
    meta['Image_Metadata_Site'] = meta['Image_Metadata_Site'].astype(int)
    return meta


def load_site_block(data_path: Path, exp_id: int, site_id: int, chunksize: int) -> pd.DataFrame:
    """Extract one experiment-site block from the full large CSV file.

    We read the file in chunks because the full dataset is too large to load at
    once. The result is the subset of rows belonging to one chosen block.
    """
    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(data_path, usecols=REQUIRED_COLUMNS, chunksize=chunksize):
        mask = (chunk['Exp_ID'] == exp_id) & (chunk['Image_Metadata_Site'] == site_id)
        if mask.any():
            parts.append(chunk.loc[mask].copy())

    if not parts:
        raise ValueError(f'No rows found for Exp_ID={exp_id}, Site={site_id}.')

    block = pd.concat(parts, ignore_index=True)
    block = block.sort_values(['track_id', 'Image_Metadata_T']).reset_index(drop=True)
    return block


def add_track_deltas(
    block: pd.DataFrame,
    signal_col: str,
    frame_to_minutes: float,
    jump_quantile: float,
) -> tuple[pd.DataFrame, float]:
    """Add node ids, time in hours, and signal changes along each track.

    This is where the table starts to become a graph-ready object. We also use
    the distribution of positive signal changes to estimate a default jump
    threshold.
    """
    block = block.copy()
    block['node_id'] = np.arange(len(block), dtype=int)
    block['time_h'] = block['Image_Metadata_T'] * frame_to_minutes / 60.0
    block['signal_value'] = block[signal_col].astype(float)
    block['signal_delta'] = block.groupby('track_id')['signal_value'].diff()

    positive_deltas = block['signal_delta'].dropna()
    positive_deltas = positive_deltas[positive_deltas > 0]
    if positive_deltas.empty:
        threshold = 0.0
    else:
        threshold = float(positive_deltas.quantile(jump_quantile))

    return block, threshold


def assign_jump_events(block: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Label rows where the signal increase is large enough to count as a jump.

    Here, "jump" is just short for a relatively large upward step in the
    signal. Operationally, it means that the signal increase from the previous
    time point is at least as large as the chosen threshold.
    """
    block = block.copy()
    block['jump_event'] = block['signal_delta'].fillna(0.0) >= threshold
    return block


def build_spatial_edges(block: pd.DataFrame, spatial_radius: float) -> pd.DataFrame:
    """Create edges between cells that are close in space at the same frame.

    These edges capture local neighbourhood structure inside one image at one
    moment in time.
    """
    rows: list[dict[str, float | int | str]] = []

    for frame, frame_df in block.groupby('Image_Metadata_T', sort=True):
        if len(frame_df) < 2:
            continue

        coords = frame_df[['objNuclei_Location_Center_X', 'objNuclei_Location_Center_Y']].to_numpy(dtype=float)
        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=spatial_radius)
        node_ids = frame_df['node_id'].to_numpy()

        for left_idx, right_idx in pairs:
            left_node = int(node_ids[left_idx])
            right_node = int(node_ids[right_idx])
            distance = float(np.linalg.norm(coords[left_idx] - coords[right_idx]))
            rows.append(
                {
                    'source_node_id': left_node,
                    'target_node_id': right_node,
                    'edge_type': 'spatial',
                    'frame': int(frame),
                    'distance': distance,
                }
            )

    return pd.DataFrame(rows)


def build_temporal_edges(block: pd.DataFrame, frame_to_minutes: float) -> pd.DataFrame:
    """Create edges that connect one cell to itself across consecutive frames.

    These edges describe how a tracked cell moves through time, which is the
    temporal part of the spatiotemporal graph.
    """
    rows: list[dict[str, float | int | str]] = []

    for _, track_df in block.groupby('track_id', sort=False):
        track_df = track_df.sort_values('Image_Metadata_T')
        if len(track_df) < 2:
            continue

        node_ids = track_df['node_id'].to_numpy()
        frames = track_df['Image_Metadata_T'].to_numpy()

        for source_node, target_node, source_frame, target_frame in zip(node_ids[:-1], node_ids[1:], frames[:-1], frames[1:]):
            frame_gap = int(target_frame - source_frame)
            rows.append(
                {
                    'source_node_id': int(source_node),
                    'target_node_id': int(target_node),
                    'edge_type': 'temporal',
                    'frame_gap': frame_gap,
                    'time_gap_h': frame_gap * frame_to_minutes / 60.0,
                }
            )

    return pd.DataFrame(rows)


def compute_future_jump_flags(block: pd.DataFrame, future_window_frames: int) -> pd.DataFrame:
    """Mark whether each node is followed by a jump later in the same track.

    This turns the analysis question into a measurable outcome: after the
    current moment, does this cell itself show another large upward step
    within the next few frames?
    """
    block = block.copy()
    block['future_self_jump'] = False

    future_values = np.zeros(len(block), dtype=bool)
    for _, track_df in block.groupby('track_id', sort=False):
        track_df = track_df.sort_values('Image_Metadata_T')
        row_ids = track_df.index.to_numpy()
        frames = track_df['Image_Metadata_T'].to_numpy()
        jumps = track_df['jump_event'].to_numpy(dtype=bool)

        for i in range(len(track_df)):
            current_frame = frames[i]
            found = False
            j = i + 1
            while j < len(track_df) and frames[j] - current_frame <= future_window_frames:
                if jumps[j]:
                    found = True
                    break
                j += 1
            future_values[row_ids[i]] = found

    block['future_self_jump'] = future_values
    return block


def annotate_spatial_exposure(block: pd.DataFrame, spatial_edges: pd.DataFrame) -> pd.DataFrame:
    """Summarize what each node sees in its local spatial neighbourhood.

    For every cell-time node, we count neighbours, count how many neighbours
    are jumping now, and compute the mean neighbour signal level.
    """
    block = block.copy()
    if spatial_edges.empty:
        block['neighbor_count'] = 0
        block['neighbor_jump_count'] = 0
        block['neighbor_mean_signal'] = np.nan
        block['neighbor_jump_now'] = False
        return block

    signal_lookup = block.set_index('node_id')['signal_value'].to_dict()
    jump_lookup = block.set_index('node_id')['jump_event'].to_dict()
    jump_counts = defaultdict(int)
    signal_sums = defaultdict(float)
    neighbor_counts = defaultdict(int)

    for edge in spatial_edges.itertuples(index=False):
        source = int(edge.source_node_id)
        target = int(edge.target_node_id)

        neighbor_counts[source] += 1
        neighbor_counts[target] += 1
        signal_sums[source] += signal_lookup[target]
        signal_sums[target] += signal_lookup[source]

        if jump_lookup[target]:
            jump_counts[source] += 1
        if jump_lookup[source]:
            jump_counts[target] += 1

    block['neighbor_jump_count'] = block['node_id'].map(lambda n: int(jump_counts.get(int(n), 0)))
    block['neighbor_count'] = block['node_id'].map(lambda n: int(neighbor_counts.get(int(n), 0)))
    block['neighbor_mean_signal'] = block['node_id'].map(
        lambda n: (signal_sums[int(n)] / neighbor_counts[int(n)]) if neighbor_counts.get(int(n), 0) else np.nan
    )
    block['neighbor_jump_now'] = block['neighbor_jump_count'] > 0
    return block


def summarise_propagation(block: pd.DataFrame, spatial_edges: pd.DataFrame, temporal_edges: pd.DataFrame, threshold: float, args: argparse.Namespace, mutation: str | None) -> dict:
    """Reduce the full graph analysis to a compact set of summary numbers.

    The key comparison is between two kinds of cell-time points:

    - exposed: at least one nearby cell is making a large upward step now,
    - unexposed: no nearby cell is making such a step now.

    We then ask whether the focal cell is more likely to make its own large
    upward step soon after.

    ``risk_difference`` is the absolute gap between those two probabilities.
    For example, 0.05 means "5 percentage points more likely".

    ``relative_risk`` is the ratio of those two probabilities. For example,
    1.8 means "about 1.8 times as likely". A value near 1 means little or no
    difference between exposed and unexposed cell-time points.
    """
    # "Exposed" means that at least one nearby cell is jumping right now.
    exposure_mask = block['neighbor_jump_now'].fillna(False)
    exposed = block.loc[exposure_mask]
    unexposed = block.loc[~exposure_mask]

    exposed_rate = float(exposed['future_self_jump'].mean()) if len(exposed) else np.nan
    unexposed_rate = float(unexposed['future_self_jump'].mean()) if len(unexposed) else np.nan
    risk_difference = exposed_rate - unexposed_rate if pd.notna(exposed_rate) and pd.notna(unexposed_rate) else np.nan
    relative_risk = exposed_rate / unexposed_rate if pd.notna(exposed_rate) and pd.notna(unexposed_rate) and unexposed_rate > 0 else np.nan

    summary = {
        'exp_id': args.exp_id,
        'site_id': args.site_id,
        'mutation': mutation,
        'signal_col': args.signal_col,
        'spatial_radius': args.spatial_radius,
        'future_window_frames': args.future_window_frames,
        'jump_threshold': threshold,
        'n_nodes': int(len(block)),
        'n_spatial_edges': int(len(spatial_edges)),
        'n_temporal_edges': int(len(temporal_edges)),
        'n_tracks': int(block['track_id'].nunique()),
        'n_frames': int(block['Image_Metadata_T'].nunique()),
        'n_exposed_nodes': int(exposure_mask.sum()),
        'n_unexposed_nodes': int((~exposure_mask).sum()),
        'future_jump_rate_if_neighbor_jumps_now': exposed_rate,
        'future_jump_rate_if_no_neighbor_jumps_now': unexposed_rate,
        'risk_difference': risk_difference,
        'relative_risk': relative_risk,
    }
    return summary


def save_outputs(block: pd.DataFrame, spatial_edges: pd.DataFrame, temporal_edges: pd.DataFrame, summary: dict, output_dir: Path) -> None:
    """Write the graph tables and summary metrics to disk.

    Saving nodes, edges, and summary files separately makes it easier to reuse the same analysis outputs in later notebooks or scripts.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    node_cols = [
        'node_id', 'Exp_ID', 'Image_Metadata_Site', 'track_id', 'Image_Metadata_T', 'time_h',
        'objNuclei_Location_Center_X', 'objNuclei_Location_Center_Y', 'Nuclear_size',
        'ERKKTR_ratio', 'FoxO3A_ratio', 'signal_value', 'signal_delta', 'jump_event',
        'neighbor_count', 'neighbor_jump_count', 'neighbor_jump_now', 'neighbor_mean_signal', 'future_self_jump',
    ]
    available_node_cols = [col for col in node_cols if col in block.columns]
    block[available_node_cols].to_csv(output_dir / 'nodes.csv.gz', index=False)

    if not spatial_edges.empty:
        spatial_edges.to_csv(output_dir / 'spatial_edges.csv.gz', index=False)
    else:
        pd.DataFrame(columns=['source_node_id', 'target_node_id', 'edge_type', 'frame', 'distance']).to_csv(
            output_dir / 'spatial_edges.csv.gz', index=False
        )

    if not temporal_edges.empty:
        temporal_edges.to_csv(output_dir / 'temporal_edges.csv.gz', index=False)
    else:
        pd.DataFrame(columns=['source_node_id', 'target_node_id', 'edge_type', 'frame_gap', 'time_gap_h']).to_csv(
            output_dir / 'temporal_edges.csv.gz', index=False
        )

    with (output_dir / 'summary.json').open('w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2)


def main() -> None:
    """Run the full single-block spatiotemporal propagation workflow.

    In order, this function loads the data, builds graph components, computes
    jump-based propagation summaries, saves the results, and prints a short
    report to the terminal.
    """
    args = parse_args()

    data_path = resolve_path(args.data_path)
    meta_path = resolve_path(args.meta_path)
    output_dir = resolve_path(args.output_dir) / f'exp_{args.exp_id}_site_{args.site_id}_{args.signal_col}'

    meta = load_metadata(meta_path)
    site_row = meta.loc[meta['Image_Metadata_Site'] == args.site_id]
    mutation = str(site_row['Mutation'].iloc[0]) if not site_row.empty else None
    frame_to_minutes = float(meta['Acquisition_frequency_min'].iloc[0])

    block = load_site_block(data_path, args.exp_id, args.site_id, args.chunksize)
    block, default_threshold = add_track_deltas(block, args.signal_col, frame_to_minutes, args.jump_quantile)
    threshold = float(args.jump_threshold) if args.jump_threshold is not None else default_threshold
    block = assign_jump_events(block, threshold)

    spatial_edges = build_spatial_edges(block, args.spatial_radius)
    temporal_edges = build_temporal_edges(block, frame_to_minutes)
    block = compute_future_jump_flags(block, args.future_window_frames)
    block = annotate_spatial_exposure(block, spatial_edges)

    summary = summarise_propagation(block, spatial_edges, temporal_edges, threshold, args, mutation)
    save_outputs(block, spatial_edges, temporal_edges, summary, output_dir)

    print('\nSpatiotemporal graph analysis complete.')
    print(f"Output directory: {output_dir}")
    print(f"Mutation: {mutation}")
    print(f"Nodes: {summary['n_nodes']:,}")
    print(f"Spatial edges: {summary['n_spatial_edges']:,}")
    print(f"Temporal edges: {summary['n_temporal_edges']:,}")
    print(f"Jump threshold on {args.signal_col}: {summary['jump_threshold']:.4f}")
    print('Interpretation of "jump": a large upward change since the previous frame.')
    print(
        'Neighbour jump now means: at least one nearby cell shows such a large upward change at the current frame.'
    )
    print(
        f'Future self jump means: the same cell shows such a large upward change within the next '
        f'{args.future_window_frames} frames.'
    )
    print(
        'Future jump rate with neighbour jump now: '
        f"{summary['future_jump_rate_if_neighbor_jumps_now']}"
    )
    print(
        'Future jump rate without neighbour jump now: '
        f"{summary['future_jump_rate_if_no_neighbor_jumps_now']}"
    )
    print(
        'Risk difference: '
        f"{summary['risk_difference']} "
        '(absolute probability gap; 0.05 means 5 percentage points higher).'
    )
    print(
        'Relative risk: '
        f"{summary['relative_risk']} "
        '(times more likely; 1 means no difference, >1 means more likely).'
    )


if __name__ == '__main__':
    main()
