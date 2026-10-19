"""Classify mismatches, filter sequencing errors by QC and plot distributions.

Depth is computed only at positions containing non-VCF mismatches (errors).
Variants are retained, but depth QC is deliberately not applied to them.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from typing import Iterable
import re

import pysam

from .bam_reader import Mismatch, _min_mapq_ok
from .plots import plot_error_histograms


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
    min_mapq: int = 20,
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
                    if read.query_sequence is None or read.query_qualities is None:
                        continue
                    if not _min_mapq_ok(read, min_mapq):
                        continue
                    if calculate_alignment_fraction(read.cigarstring or "") < min_alignment_fraction:
                        continue
                    coverage += 1
                depth[(chrom, pos)] = coverage

    return depth


def count_rate_denominators(
    bam: pysam.AlignmentFile,
    fasta: pysam.FastaFile | None,
    chrom: str,
    start: int,
    end: int,
    variant_loci: set[tuple[str, int]],
    min_mapq: int = 20,
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
    min_base_qual: int = 0,
    depth_bin_size: int = 5,
    quality_bin_size: int = 5,
) -> tuple[Counter, Counter, Counter]:
    """Count eligible bases for real error rates across ALL reference positions.

    A region [start,end) is owned by exactly one worker, even when reads
    overlap window boundaries. VCF loci and non-ACGT reference loci are masked.

    Returns (bases_by_phred, bases_by_filtered_depth, bases_by_depth_and_phred_bin).
    All three are denominator counts, NOT error counts.
    Depth is evaluated before the min_base_qual threshold, consistent with
    calculate_depth(), while denominators count only quality-eligible bases.
    """
    if fasta is None:
        raise ValueError("count_rate_denominators requires an indexed FASTA")
    if depth_bin_size <= 0 or quality_bin_size <= 0:
        raise ValueError("Bin sizes must be positive")
    if start >= end:
        return Counter(), Counter(), Counter()
    reference = fasta.fetch(chrom, start, end).upper()
    quality_totals: Counter = Counter()
    depth_base_totals: Counter = Counter()
    quality_depth_base_totals: Counter = Counter()

    for column in bam.pileup(
        chrom, start, end,
        truncate=True,
        stepper="nofilter",
        min_base_quality=0,
        min_mapping_quality=0,
        ignore_overlaps=False,
        max_depth=100_000,
    ):
        ref0 = column.reference_pos
        pos1 = ref0 + 1
        if (chrom, pos1) in variant_loci:
            continue
        if reference[ref0 - start] not in "ACGT":
            continue

        depth = 0
        eligible_qualities: list[int] = []
        for pileup_read in column.pileups:
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            read = pileup_read.alignment
            if read.query_sequence is None or read.query_qualities is None:
                continue
            if not _min_mapq_ok(read, min_mapq):
                continue
            if calculate_alignment_fraction(read.cigarstring or "") < min_alignment_fraction:
                continue
            qpos = pileup_read.query_position
            if qpos is None:
                continue
            depth += 1
            qual = int(read.query_qualities[qpos])
            if qual >= min_base_qual:
                eligible_qualities.append(qual)

        if depth < min_depth:
            continue
        depth_bin = (depth // depth_bin_size) * depth_bin_size
        for qual in eligible_qualities:
            quality_totals[qual] += 1
            qual_bin = (qual // quality_bin_size) * quality_bin_size
            quality_depth_base_totals[(depth_bin, qual_bin)] += 1
        depth_base_totals[depth] += len(eligible_qualities)

    return quality_totals, depth_base_totals, quality_depth_base_totals


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



def process_technology(
    mismatches: Iterable[Mismatch],
    variant_positions: set[tuple[str, int, str, str]],
    bam_path: str,
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
    *,
    min_mapq: int = 20,
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

    errors = filter_mismatches_by_qc(
        errors,
        depth,
        min_alignment_fraction=min_alignment_fraction,
        min_depth=min_depth,
    )

    # Variants outside the error positions have no calculated depth.
    # Filtering them with min_depth would incorrectly remove valid variants.
    
    return errors, variants




