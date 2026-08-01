"""
vcf_reader.py

"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pysam


def load_variant_positions(
    vcf_path: str | Path,
    min_qual: Optional[float] = None,
    pass_only: bool = True,
) -> set[tuple[str, int, str, str]]:

    positions: set[tuple[str, int, str, str]] = set()
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            if pass_only and rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue
            if min_qual is not None and rec.qual is not None and rec.qual < min_qual:
                continue
            for alt in rec.alts or ():
                positions.add((rec.chrom, rec.pos, rec.ref, alt))
    return positions


def variant_position_only_index(
    variants: set[tuple[str, int, str, str]]
) -> set[tuple[str, int]]:

    return {(chrom, pos) for chrom, pos, _, _ in variants}
