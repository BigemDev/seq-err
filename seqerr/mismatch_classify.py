"""
mismatch_classify.py

Classifies detected mismatches as either known variants (present in VCF) or sequencing errors.

"""
from __future__ import annotations

from typing import Iterable

from .bam_reader import Mismatch

from collections import defaultdict

import pandas as pd
import sys
import os
import re
import argparse
import pysam
import matplotlib.pyplot as plt

def calculate_alignment_fraction(cigar):

    if not cigar:
        return 0

    operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    aligned = 0
    read_length = 0

    for length, operation in operations:

        length = int(length)

        if operation in "MIS=X":
            read_length += length

        if operation in "M=X":
            aligned += length

    if read_length == 0:
        return 0

    return aligned / read_length

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
) -> dict[tuple[str, int], int]:
    """Calculate filtered read depth at selected 1-based positions."""
    depth: dict[tuple[str, int], int] = defaultdict(int)

    for position in interesting_positions:
        depth[position] = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for column in bam.pileup(
            truncate=False,
            stepper="all",
            min_base_quality=0,
            min_mapping_quality=0,
        ):
            chrom = column.reference_name
            position = column.reference_pos + 1
            key = (chrom, position)

            if key not in interesting_positions:
                continue

            coverage = 0

            for pileup_read in column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue

                alignment = pileup_read.alignment

                if alignment.is_duplicate or alignment.is_unmapped:
                    continue

                fraction = calculate_alignment_fraction(
                    alignment.cigarstring
                )

                if fraction < min_alignment_fraction:
                    continue

                coverage += 1

            depth[key] = coverage

    return dict(depth)

def filter_mismatches_by_qc(
    mismatches: Iterable[Mismatch],
    depth: dict[tuple[str, int], int],
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
) -> list[Mismatch]:
    """Keep mismatches meeting alignment-fraction and depth thresholds."""
    filtered: list[Mismatch] = []

    for mismatch in mismatches:
        fraction = calculate_alignment_fraction(mismatch.cigar)

        position_depth = depth.get(
            (mismatch.chrom, mismatch.ref_pos),
            0,
        )

        if fraction < min_alignment_fraction:
            continue

        if position_depth < min_depth:
            continue

        filtered.append(mismatch)

    return filtered

def plot_error_histograms(
    errors: list[Mismatch],
    technology: str,
    output_dir: str = "histograms",
):
    output_path = os.path.join(output_dir, technology)
    os.makedirs(output_path, exist_ok=True)

    datasets = [
        ("depth", "Depth", [m.depth for m in errors]),
        ("base_qual", "Base quality (Phred)", [m.base_qual for m in errors]),
        ("mapping_qual", "Mapping quality (MAPQ)", [m.mapping_qual for m in errors]),
    ]

    for filename, xlabel, values in datasets:
        if not values:
            continue

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(values, bins=30, edgecolor="black", alpha=0.75)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Number of errors")
        ax.set_title(f"{technology.upper()} error distribution: {xlabel}")

        fig.tight_layout()
        fig.savefig(os.path.join(output_path, f"{filename}_errors_histogram.png"),dpi=150,)
        plt.close(fig)

def process_technology(
    mismatches: Iterable[Mismatch],
    variant_positions: set[tuple[str, int, str, str]],
    bam_path: str,
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
) -> tuple[list[Mismatch], list[Mismatch]]:

    errors, variants = classify_mismatches(
        mismatches,
        variant_positions,
    )

    all_mismatches = errors + variants

    interesting_positions = {
        (m.chrom, m.ref_pos)
        for m in all_mismatches
    }

    depth = calculate_depth(
        bam_path,
        interesting_positions,
        min_alignment_fraction=min_alignment_fraction,
    )

    for mismatch in all_mismatches:
        mismatch.depth = depth.get(
            (mismatch.chrom, mismatch.ref_pos),
            0,
        )

    errors = filter_mismatches_by_qc(
        errors,
        depth,
        min_alignment_fraction=min_alignment_fraction,
        min_depth=min_depth,
    )

    variants = filter_mismatches_by_qc(
        variants,
        depth,
        min_alignment_fraction=min_alignment_fraction,
        min_depth=min_depth,
    )

    plot_error_histograms(
    errors,
    technology=all_mismatches[0].technology if all_mismatches else "unknown",
    output_dir="histograms",
)

    return errors, variants

