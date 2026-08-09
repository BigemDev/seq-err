"""
mismatch_classify.py

Classifies detected mismatches as either known variants (present in VCF) or sequencing errors.

"""
from __future__ import annotations

from typing import Iterable, Iterator

from .bam_reader import Mismatch


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


def iter_classify_mismatches(
    mismatches: Iterable[Mismatch],
    variant_positions: set[tuple[str, int, str, str]],
) -> Iterator[tuple[Mismatch, str]]:
    for m in mismatches:
        key = (m.chrom, m.ref_pos, m.ref_base, m.read_base)
        label = "variant" if key in variant_positions else "sequencing_error"
        yield m, label
