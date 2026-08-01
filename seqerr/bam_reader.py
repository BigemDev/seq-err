"""
bam_reader.py

"""
from __future__ import annotations

import csv
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator, Optional

import pysam


@dataclass
class Mismatch:
    read_name: str
    technology: str          # free-form label, e.g. "illumina" / "bgi"
    chrom: str
    ref_pos: int              # 1-based reference position
    ref_base: str
    read_base: str
    base_qual: int            # Phred base quality at this position
    mapping_qual: int         # read.mapping_quality
    cigar: str                 # full CIGAR string of the read
    is_read1: bool
    is_reverse: bool

    def as_row(self) -> dict:
        return asdict(self)


CSV_FIELDS = list(Mismatch.__dataclass_fields__.keys())


def _min_mapq_ok(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return (
        not read.is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and not read.is_duplicate
        and read.mapping_quality >= min_mapq
    )


def iter_mismatches(
    bam_path: str | Path,
    technology: str,
    reference_fasta: Optional[str | Path] = None,
    regions: Optional[list[tuple[str, int, int]]] = None,
    min_mapq: int = 1,
    min_base_qual: int = 0,
) -> Iterator[Mismatch]:
    fasta = pysam.FastaFile(str(reference_fasta)) if reference_fasta else None

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        fetch_iter = (
            (bam.fetch(c, s, e) for c, s, e in regions)
            if regions
            else [bam.fetch()]
        )
        for region_reads in fetch_iter:
            for read in region_reads:
                if not _min_mapq_ok(read, min_mapq):
                    continue
                if read.query_sequence is None or read.query_qualities is None:
                    continue

                chrom = read.reference_name
                aligned_pairs = read.get_aligned_pairs(
                    with_seq=(fasta is None),
                    matches_only=False,
                )

                for query_pos, ref_pos, *seq in aligned_pairs:
                    if query_pos is None or ref_pos is None:
                        continue  # insertion or deletion, not a substitution
                    read_base = read.query_sequence[query_pos].upper()
                    if fasta is not None:
                        ref_base = fasta.fetch(chrom, ref_pos, ref_pos + 1).upper()
                    else:
                        ref_base = seq[0].upper()

                    if ref_base == read_base or ref_base not in "ACGT":
                        continue  # match, or ref is N/ambiguous

                    base_qual = read.query_qualities[query_pos]
                    if base_qual < min_base_qual:
                        continue

                    yield Mismatch(
                        read_name=read.query_name,
                        technology=technology,
                        chrom=chrom,
                        ref_pos=ref_pos + 1,  # report 1-based
                        ref_base=ref_base,
                        read_base=read_base,
                        base_qual=base_qual,
                        mapping_qual=read.mapping_quality,
                        cigar=read.cigarstring or "",
                        is_read1=read.is_read1,
                        is_reverse=read.is_reverse,
                    )

    if fasta is not None:
        fasta.close()


def write_mismatches_csv(mismatches: Iterator[Mismatch], out_path: str | Path) -> int:
    """Stream Mismatch records to a CSV file. Returns the number of rows written."""
    n = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for m in mismatches:
            writer.writerow(m.as_row())
            n += 1
    return n
