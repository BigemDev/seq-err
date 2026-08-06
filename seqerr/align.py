"""
align.py

"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import mappy
import pysam


def _sniff_fastx_format(path: str | Path) -> str:
    path = Path(path)
    if not path.exists():
        raise ValueError(f"reads file not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"reads file is empty: {path}")

    is_gz = path.suffix == ".gz" or _looks_gzipped(path)
    opener = gzip.open if is_gz else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("@"):
                    return "fastq"
                raise ValueError(
                    f"{path} does not look like FASTA or FASTQ "
                    f"(first non-blank line starts with {line[:20]!r}, "
                    f"expected '>' for FASTA or '@' for FASTQ)"
                )
    except OSError as e:
        raise ValueError(f"could not read {path}: {e}") from e

    raise ValueError(f"reads file has no content lines: {path}")


def _looks_gzipped(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _cigar_tuples_from_mappy(hit, query_len) -> list[tuple[int, int]]:
    cigar = []

    if hit.q_st > 0:
        cigar.append((4, hit.q_st))

    cigar.extend((op, length) for length, op in hit.cigar)

    tail = query_len - hit.q_en
    if tail > 0:
        cigar.append((4, tail))

    return cigar

def _build_md_nm_mismatches(
    ref_seq: str, read_seq: str, ref_start: int, cigar: list[tuple[int, int]]
) -> tuple[str, int, list[dict]]:
    md_parts: list[str] = []
    run = 0
    nm = 0
    rpos = ref_start
    qpos = 0
    mismatches: list[dict] = []

    for op, length in cigar:
        if op in (0, 7, 8):  # M, =, X
            for _ in range(length):
                rb = ref_seq[rpos].upper()
                qb = read_seq[qpos].upper()
                if rb == qb:
                    run += 1
                else:
                    md_parts.append(str(run))
                    md_parts.append(rb)
                    run = 0
                    nm += 1
                    mismatches.append({
                        "ref_offset": rpos,   # offset within the fetched ref region
                        "query_offset": qpos,  # offset within the aligned read seq
                        "ref_base": rb,
                        "read_base": qb,
                    })
                rpos += 1
                qpos += 1
        elif op == 2:  # D
            md_parts.append(str(run))
            md_parts.append("^" + ref_seq[rpos:rpos + length].upper())
            run = 0
            nm += length
            rpos += length
        elif op == 1:  # I
            nm += length
            qpos += length
        elif op == 4:  # S
            qpos += length
        elif op == 3:  # N
            rpos += length

    md_parts.append(str(run))
    return "".join(md_parts), nm, mismatches


def _build_md_and_nm(
    ref_seq: str, read_seq: str, ref_start: int, cigar: list[tuple[int, int]]
) -> tuple[str, int]:
    md, nm, _ = _build_md_nm_mismatches(ref_seq, read_seq, ref_start, cigar)
    return md, nm


def map_reads_to_bam(
    reads_path: str | Path,
    reference_fasta: str | Path,
    out_bam: str | Path,
    preset: str = "sr",
    technology_tag: Optional[str] = None,
    min_mapq: int = 0,
) -> int:
    reference_fasta = str(reference_fasta)
    fmt = _sniff_fastx_format(reads_path)  # error on bad/empty/unrecognised input
    print(f"[align] detected {fmt.upper()} input: {reads_path}")

    aligner = mappy.Aligner(reference_fasta, preset=preset)
    if not aligner:
        raise RuntimeError(f"failed to load/index reference: {reference_fasta}")

    ref_fasta = pysam.FastaFile(reference_fasta)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": ref_fasta.get_reference_length(name)}
               for name in ref_fasta.references],
    }
    if technology_tag:
        header["RG"] = [{"ID": technology_tag, "SM": technology_tag}]

    tmp_bam = str(out_bam) + ".unsorted.tmp.bam"
    n_written = 0

    with pysam.AlignmentFile(tmp_bam, "wb", header=header) as bam_out:
        for name, seq, qual in mappy.fastx_read(str(reads_path)):
            hit = None
            for h in aligner.map(seq):
                if h.is_primary:
                    hit = h
                    break
            if hit is None:
                continue
            if hit.mapq < min_mapq:
                continue

            a = pysam.AlignedSegment(header=bam_out.header)
            a.query_name = name
            a.flag = 16 if hit.strand == -1 else 0
            a.reference_id = bam_out.get_tid(hit.ctg)
            a.reference_start = hit.r_st
            a.mapping_quality = hit.mapq

            query_seq = seq if hit.strand == 1 else mappy.revcomp(seq)
            a.cigar = _cigar_tuples_from_mappy(hit, len(query_seq))

            a.query_sequence = query_seq
            if qual:
                q = qual if hit.strand == 1 else qual[::-1]
                a.query_qualities = pysam.qualitystring_to_array(q)
            else:
                a.query_qualities = [30] * len(query_seq)  # no qualities

            ref_region = ref_fasta.fetch(hit.ctg, hit.r_st, hit.r_en)
           
            md, nm = _build_md_and_nm(ref_region, query_seq, 0, a.cigar)

            a.set_tag("MD", md)
            a.set_tag("NM", nm)
            if technology_tag:
                a.set_tag("RG", technology_tag)

            bam_out.write(a)
            n_written += 1

    ref_fasta.close()
    pysam.sort("-o", str(out_bam), tmp_bam)
    pysam.index(str(out_bam))
    Path(tmp_bam).unlink(missing_ok=True)
    return n_written


def _map_one(aligner, seq, min_mapq):
    hit = None
    for h in aligner.map(seq):
        if h.is_primary:
            hit = h
            break
    if hit is None or hit.mapq < min_mapq:
        return None
    return hit


def _build_segment(bam_out, aligner, ref_fasta, name, seq, qual, hit, technology_tag):
    a = pysam.AlignedSegment(header=bam_out.header)
    a.query_name = name
    a.flag = 16 if hit.strand == -1 else 0
    a.reference_id = bam_out.get_tid(hit.ctg)
    a.reference_start = hit.r_st
    a.mapping_quality = hit.mapq

    query_seq = seq if hit.strand == 1 else mappy.revcomp(seq)
    a.cigar = _cigar_tuples_from_mappy(hit, len(query_seq))

    a.query_sequence = query_seq
    if qual:
        q = qual if hit.strand == 1 else qual[::-1]
        a.query_qualities = pysam.qualitystring_to_array(q)
    else:
        a.query_qualities = [30] * len(query_seq)

    ref_region = ref_fasta.fetch(hit.ctg, hit.r_st, hit.r_en)
    md, nm = _build_md_and_nm(ref_region, query_seq, 0, a.cigar)
    a.set_tag("MD", md)
    a.set_tag("NM", nm)
    if technology_tag:
        a.set_tag("RG", technology_tag)
    return a


def map_paired_reads_to_bam(
    reads_r1: str | Path,
    reads_r2: str | Path,
    reference_fasta: str | Path,
    out_bam: str | Path,
    preset: str = "sr",
    technology_tag: Optional[str] = None,
    min_mapq: int = 0,
) -> int:
    reference_fasta = str(reference_fasta)
    _sniff_fastx_format(reads_r1)
    _sniff_fastx_format(reads_r2)

    aligner = mappy.Aligner(reference_fasta, preset=preset)
    if not aligner:
        raise RuntimeError(f"failed to load/index reference: {reference_fasta}")

    ref_fasta = pysam.FastaFile(reference_fasta)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": ref_fasta.get_reference_length(name)}
               for name in ref_fasta.references],
    }
    if technology_tag:
        header["RG"] = [{"ID": technology_tag, "SM": technology_tag}]

    tmp_bam = str(out_bam) + ".unsorted.tmp.bam"
    n_written = 0

    with pysam.AlignmentFile(tmp_bam, "wb", header=header) as bam_out:
        for (n1, s1, q1), (n2, s2, q2) in zip(
            mappy.fastx_read(str(reads_r1)), mappy.fastx_read(str(reads_r2))
        ):
            name = n1.split()[0]
            if name.endswith(("/1", "/2")):
                name = name[:-2]

            h1 = _map_one(aligner, s1, min_mapq)
            h2 = _map_one(aligner, s2, min_mapq)

            if h1 is None and h2 is None:
                continue

            a1 = _build_segment(bam_out, aligner, ref_fasta, name, s1, q1, h1, technology_tag) if h1 else None
            a2 = _build_segment(bam_out, aligner, ref_fasta, name, s2, q2, h2, technology_tag) if h2 else None

            for a, mate, is_read1 in ((a1, a2, True), (a2, a1, False)):
                if a is None:
                    continue
                a.flag |= 0x1
                a.is_read1 = is_read1
                a.is_read2 = not is_read1
                if mate is None:
                    a.flag |= 0x8  # mate unmapped
                else:
                    a.next_reference_id = mate.reference_id
                    a.next_reference_start = mate.reference_start
                    if mate.is_reverse:
                        a.flag |= 0x20
                    if a.reference_id == mate.reference_id:
                        span_end = max(a.reference_end, mate.reference_end)
                        span_start = min(a.reference_start, mate.reference_start)
                        tlen = span_end - span_start
                        a.template_length = tlen if a.reference_start <= mate.reference_start else -tlen
                    a.flag |= 0x2  # proper pair
                bam_out.write(a)
                n_written += 1

    ref_fasta.close()
    pysam.sort("-o", str(out_bam), tmp_bam)
    pysam.index(str(out_bam))
    Path(tmp_bam).unlink(missing_ok=True)
    return n_written


@dataclass
class SingleReadAlignment:
    read_name: str
    mapped: bool
    chrom: Optional[str] = None
    ref_start: Optional[int] = None      # 0-based, SAM convention
    ref_end: Optional[int] = None        # 0-based, exclusive
    strand: Optional[str] = None         # "+" or "-"
    mapping_qual: Optional[int] = None
    cigar_string: Optional[str] = None
    cigar_tuples: Optional[list[tuple[int, int]]] = None
    md_tag: Optional[str] = None
    nm: Optional[int] = None
    mismatches: Optional[list[dict]] = None  # [{ref_pos (1-based), ref_base, read_base}, ...]

    def pretty(self) -> str:
        if not self.mapped:
            return f"{self.read_name}: unmapped"
        lines = [
            f"{self.read_name} -> {self.chrom}:{self.ref_start + 1}-{self.ref_end} "
            f"({self.strand}) MAPQ={self.mapping_qual}",
            f"  CIGAR: {self.cigar_string}",
            f"  MD:    {self.md_tag}   NM={self.nm}",
        ]
        if self.mismatches:
            lines.append(f"  {len(self.mismatches)} mismatch(es):")
            for m in self.mismatches:
                lines.append(
                    f"    ref_pos={m['ref_pos']} {m['ref_base']}->{m['read_base']} "
                    f"base_qual={m['base_qual']}"
                )
        else:
            lines.append("  0 mismatches")
        return "\n".join(lines)


def map_single_read(
    sequence: str,
    reference_fasta: str | Path,
    read_name: str = "read1",
    quality: Optional[str] = None,
    preset: str = "sr",
) -> SingleReadAlignment:
    aligner = mappy.Aligner(str(reference_fasta), preset=preset)
    if not aligner:
        raise RuntimeError(f"failed to load/index reference: {reference_fasta}")

    hit = None
    for h in aligner.map(sequence):
        if h.is_primary:
            hit = h
            break
    if hit is None:
        return SingleReadAlignment(read_name=read_name, mapped=False)

    strand = "+" if hit.strand == 1 else "-"
    query_seq = sequence if hit.strand == 1 else mappy.revcomp(sequence)
    query_qual = quality if quality is None else (quality if hit.strand == 1 else quality[::-1])
    cigar_tuples = cigar_tuples = _cigar_tuples_from_mappy(hit, len(query_seq))
    cigar_string = "".join(f"{length}{'MIDNSHP=X'[op]}" for op, length in cigar_tuples)

    ref_fasta = pysam.FastaFile(str(reference_fasta))
    ref_region = ref_fasta.fetch(hit.ctg, hit.r_st, hit.r_en)
    ref_fasta.close()

    aligned_query = query_seq[hit.q_st:hit.q_en]

    #md, nm, raw_mismatches = _build_md_nm_mismatches(ref_region, aligned_query, 0, cigar_tuples)
    
    md, nm, raw_mismatches = _build_md_nm_mismatches(ref_region, query_seq, 0, cigar_tuples)

    mismatches = []
    for m in raw_mismatches:
        base_qual = None
        if query_qual is not None:
            base_qual = ord(query_qual[hit.q_st + m["query_offset"]]) - 33
        mismatches.append({
            "ref_pos": hit.r_st + m["ref_offset"] + 1,  # 1-based
            "ref_base": m["ref_base"],
            "read_base": m["read_base"],
            "base_qual": base_qual,
        })

    return SingleReadAlignment(
        read_name=read_name,
        mapped=True,
        chrom=hit.ctg,
        ref_start=hit.r_st,
        ref_end=hit.r_en,
        strand=strand,
        mapping_qual=hit.mapq,
        cigar_string=cigar_string,
        cigar_tuples=cigar_tuples,
        md_tag=md,
        nm=nm,
        mismatches=mismatches,
    )