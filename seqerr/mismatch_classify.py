"""Classify mismatches, filter sequencing errors by QC and plot distributions.

Depth is computed only at positions containing non-VCF mismatches (errors).
Variants are retained, but depth QC is deliberately not applied to them.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Iterable
import re

import pysam

from .bam_reader import Mismatch, _min_mapq_ok


@lru_cache(maxsize=100_000)
def calculate_alignment_fraction(cigar: str) -> float:
    if not cigar:
        return 0.0

    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    aligned = 0
    read_length = 0

    for length_text, operation in operations:
        length = int(length_text)
        if operation in "MIS=X":
            read_length += length
        if operation in "M=X":
            aligned += length

    return aligned / read_length if read_length else 0.0


def classify_mismatches(
    mismatches: Iterable[Mismatch],
    variant_positions: set[tuple[str, int, str, str]],
) -> tuple[list[Mismatch], list[Mismatch]]:
    errors: list[Mismatch] = []
    variants: list[Mismatch] = []
    for m in mismatches:
        key = (m.chrom, m.ref_pos, m.ref_base, m.read_base)
        if key in variant_positions:
            variants.append(m)
        else:
            errors.append(m)
    return errors, variants


def calculate_depth(
    bam_path: str,
    interesting_positions: set[tuple[str, int]],
    min_alignment_fraction: float = 0.8,
    min_mapq: int = 0,
) -> dict[tuple[str, int], int]:
    """Count eligible read observations only at requested 1-based positions.

    Per chromosome, use one bounded, indexed pileup from minimum to maximum
    interesting position. This avoids rescanning the WHOLE BAM per worker.
    The numerator and denominator use the same read flag/MAPQ rules.
    """
    depth: dict[tuple[str, int], int] = {position: 0 for position in interesting_positions}
    if not interesting_positions:
        return depth

    positions_by_chrom: dict[str, set[int]] = defaultdict(set)
    for chrom, position in interesting_positions:
        positions_by_chrom[chrom].add(position)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for chrom, positions in positions_by_chrom.items():
            start = min(positions) - 1  # 0-based inclusive
            stop = max(positions)       # 0-based exclusive

            for column in bam.pileup(
                chrom, start, stop,
                truncate=True,
                stepper="nofilter",  # explicit consistent filtering below
                min_base_quality=0,
                min_mapping_quality=0,
                ignore_overlaps=False,  # match per-read extraction
                max_depth=100_000,
            ):
                pos = column.reference_pos + 1
                if pos not in positions:
                    continue

                coverage = 0
                for pileup_read in column.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    read = pileup_read.alignment
                    if not _min_mapq_ok(read, min_mapq):
                        continue
                    if calculate_alignment_fraction(read.cigarstring or "") < min_alignment_fraction:
                        continue
                    coverage += 1
                depth[(chrom, pos)] = coverage

    return depth


def filter_mismatches_by_qc(
    mismatches: Iterable[Mismatch],
    depth: dict[tuple[str, int], int],
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
) -> list[Mismatch]:
    filtered: list[Mismatch] = []
    for mismatch in mismatches:
        if calculate_alignment_fraction(mismatch.cigar) < min_alignment_fraction:
            continue
        if depth.get((mismatch.chrom, mismatch.ref_pos), 0) < min_depth:
            continue
        filtered.append(mismatch)
    return filtered


def _save_count_plot(
    counts: dict[int, int | float] | Counter,
    title: str,
    xlabel: str,
    ylabel: str,
    filepath: Path,
) -> None:
    if not counts:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = sorted(counts)
    y = [counts[k] for k in x]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x, y, width=0.9, edgecolor="black", linewidth=0.3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(filepath, dpi=150)
    plt.close(fig)


def plot_error_histograms_from_counts(
    base_qual_counts: Counter,
    depth_counts: Counter,
    mapq_counts: Counter,
    weighted_quality_counts: dict[int, float],
    technology: str,
    output_dir: str = "histograms",
) -> None:
    """Plots from aggregated counts: safe to use after multiprocessing."""
    output_path = Path(output_dir) / technology
    output_path.mkdir(parents=True, exist_ok=True)
    datasets = [
        (depth_counts, "Depth", "Number of errors", "depth_errors_histogram.png"),
        (base_qual_counts, "Base quality (Phred)", "Number of errors", "base_qual_errors_histogram.png"),
        (mapq_counts, "Mapping quality (MAPQ)", "Number of errors", "mapping_qual_errors_histogram.png"),
        (
            weighted_quality_counts,
            "Base quality (Phred)",
            "Sum of 1 / depth (coverage-weighted errors)",
            "base_qual_depth_weighted_histogram.png",
        ),
    ]
    for counts, xlabel, ylabel, filename in datasets:
        _save_count_plot(counts, f"{technology.upper()} errors: {xlabel}", xlabel, ylabel, output_path / filename)

"""
def plot_error_histograms(
    errors: list[Mismatch],
    technology: str,
    output_dir: str = "histograms",
) -> None:
    
    qual = Counter(m.base_qual for m in errors)
    depths = Counter(m.depth for m in errors)
    mapq = Counter(m.mapping_qual for m in errors)
    weighted: dict[int, float] = defaultdict(float)
    
    for m in errors:
        if m.depth > 0:
            weighted[m.base_qual] += 1.0 / m.depth
    plot_error_histograms_from_counts(qual, depths, mapq, dict(weighted), technology, output_dir)
"""
def plot_quality_depth_heatmap_from_counts(
    quality_by_depth_counts: Counter,
    total_by_depth_counts: Counter,
    technology: str,
    output_dir: str = "histograms",
) -> None:
    if not quality_by_depth_counts or not total_by_depth_counts:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    output_path = Path(output_dir) / technology
    output_path.mkdir(parents=True, exist_ok=True)

    depth_bins = sorted(total_by_depth_counts.keys())
    quality_bins = sorted({quality_bin for (_, quality_bin) in quality_by_depth_counts.keys()})

    matrix = np.zeros((len(quality_bins), len(depth_bins)))

    for row, quality_bin in enumerate(quality_bins):
        for col, depth_bin in enumerate(depth_bins):
            count = quality_by_depth_counts.get((depth_bin, quality_bin), 0)
            total = total_by_depth_counts.get(depth_bin, 0)

            if total > 0:
                matrix[row, col] = count / total
            else:
                matrix[row, col] = 0.0

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(matrix, aspect="auto", origin="lower")

    ax.set_xticks(range(len(depth_bins)))
    ax.set_xticklabels(depth_bins, rotation=45, ha="right")

    ax.set_yticks(range(len(quality_bins)))
    ax.set_yticklabels(quality_bins)

    ax.set_xlabel("Depth bin")
    ax.set_ylabel("Base quality bin")
    ax.set_title(f"{technology.upper()} error distribution: base quality vs depth")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Proportion of errors")

    fig.tight_layout()
    fig.savefig(output_path / "base_quality_vs_depth_heatmap.png", dpi=150)
    plt.close(fig)



def process_technology(
    mismatches: Iterable[Mismatch],
    variant_positions: set[tuple[str, int, str, str]],
    bam_path: str,
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
    *,
    min_mapq: int = 0,
    plot_histograms: bool = True,
    technology: str | None = None,
    output_dir: str = "histograms",
) -> tuple[list[Mismatch], list[Mismatch]]:
    """Classify & QC errors; keep VCF variants without minimum-depth filtering.

    Call on a bounded window rather than a whole 40+ GB collection.
    """
    errors, variants = classify_mismatches(mismatches, variant_positions)

    interesting_positions = {(m.chrom, m.ref_pos) for m in errors}
    depth = calculate_depth(
        bam_path,
        interesting_positions,
        min_alignment_fraction=min_alignment_fraction,
        min_mapq=min_mapq,
    )

    for mismatch in errors:
        mismatch.depth = depth.get((mismatch.chrom, mismatch.ref_pos), 0)

    errors = filter_mismatches_by_qc(errors,depth,min_alignment_fraction=min_alignment_fraction,min_depth=min_depth,)

    return errors, variants
