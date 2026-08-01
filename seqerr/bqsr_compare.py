"""
bqsr_compare.py

"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .bam_reader import Mismatch


@dataclass
class BqsrDelta:
    read_name: str
    chrom: str
    ref_pos: int
    ref_base: str
    read_base: str
    qual_before: int
    qual_after: int
    delta: int  # qual_after - qual_before
    mapq_before: int
    mapq_after: int


def _index_by_key(mismatches: Iterable[Mismatch]) -> dict[tuple, Mismatch]:
    return {
        (m.read_name, m.chrom, m.ref_pos, m.is_read1): m
        for m in mismatches
    }


def compare_bqsr(
    mismatches_before: Iterable[Mismatch],
    mismatches_after: Iterable[Mismatch],
) -> list[BqsrDelta]:
    """
    Pair up mismatches found in a pre-BQSR BAM and a post-BQSR BAM at the
    same read/position and report how the base quality changed.

    Only positions present in BOTH sets are reported (a base that stopped
    being a mismatch after BQSR realignment/soft-clipping would simply drop
    out -- that is worth noting separately, see `dropped_after_bqsr`).
    """
    before_idx = _index_by_key(mismatches_before)
    after_idx = _index_by_key(mismatches_after)

    deltas = []
    for key, before in before_idx.items():
        after = after_idx.get(key)
        if after is None:
            continue
        deltas.append(
            BqsrDelta(
                read_name=before.read_name,
                chrom=before.chrom,
                ref_pos=before.ref_pos,
                ref_base=before.ref_base,
                read_base=before.read_base,
                qual_before=before.base_qual,
                qual_after=after.base_qual,
                delta=after.base_qual - before.base_qual,
                mapq_before=before.mapping_qual,
                mapq_after=after.mapping_qual,
            )
        )
    return deltas


def dropped_after_bqsr(
    mismatches_before: Iterable[Mismatch],
    mismatches_after: Iterable[Mismatch],
) -> list[Mismatch]:
    """Mismatches present before BQSR but no longer present after (e.g. the
    read was filtered, realigned, or the base was masked)."""
    after_keys = {
        (m.read_name, m.chrom, m.ref_pos, m.is_read1) for m in mismatches_after
    }
    return [
        m
        for m in mismatches_before
        if (m.read_name, m.chrom, m.ref_pos, m.is_read1) not in after_keys
    ]


def write_bqsr_deltas_csv(deltas: list[BqsrDelta], out_path: str | Path) -> None:
    fields = list(BqsrDelta.__dataclass_fields__.keys())
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for d in deltas:
            writer.writerow(d.__dict__)
