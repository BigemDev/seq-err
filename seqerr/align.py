"""
align.py

Maps raw FASTA/FASTQ reads to a reference genome via bwa mem (bwapy) and
writes/returns CIGAR, MD, and mismatch data.

"""
from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import sys
import types
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pysam

if "imp" not in sys.modules:
    try:
        import imp  # noqa: F401
    except ModuleNotFoundError:
        # bwapy unconditionally does `import imp`, which was removed in
        # Python 3.12, and actually relies on imp.find_module(name) to
        # locate its compiled `bwalib` extension by file path so it can
        # dlopen() it directly (bwalib's module-init symbol doesn't match
        # its file name, so a normal `import bwalib` fails). Provide a
        # minimal stand-in backed by importlib so that lookup still works.
        import importlib.util as _importlib_util

        def _find_module(name, path=None):
            spec = _importlib_util.find_spec(name, path)
            if spec is None or spec.origin is None:
                raise ImportError(f"No module named {name!r}")
            return None, spec.origin, None

        _imp_shim = types.ModuleType("imp")
        _imp_shim.find_module = _find_module
        sys.modules["imp"] = _imp_shim

from bwapy import BwaAligner


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


def _fastx_read(path: str | Path):
    with pysam.FastxFile(str(path)) as fh:
        for rec in fh:
            yield rec.name, rec.sequence, rec.quality


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]



_BWA_CIGAR_OP = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _parse_bwa_cigar(cigar_string: str) -> list[tuple[int, int]]:
    return [(_BWA_CIGAR_OP[op], int(length)) for length, op in _CIGAR_RE.findall(cigar_string)]


def _ref_span(cigar_tuples: list[tuple[int, int]]) -> int:
    return sum(length for op, length in cigar_tuples if op in (0, 2, 3, 7, 8))


_BWA_INDEX_EXTS = (".bwt", ".pac", ".ann", ".amb", ".sa")


def _ensure_bwa_index(reference_fasta: str) -> None:
    if all(Path(reference_fasta + ext).exists() for ext in _BWA_INDEX_EXTS):
        return
    if shutil.which("bwa") is None:
        raise RuntimeError(
            "the 'bwa' command-line tool is required to build a bwa index "
            "(bwapy only loads pre-built indexes) but was not found on PATH"
        )
    result = subprocess.run(
        ["bwa", "index", reference_fasta],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bwa index failed for {reference_fasta}:\n{result.stderr}")




_PRESET_TO_BWA_OPTS = {
    "sr": "",
    "map-ont": "-x ont2d",
    "map-pb": "-x pacbio",
    "map-hifi": "-x pacbio",
    "intractg": "-x intractg",
}


def _bwa_options_for_preset(preset: str) -> str:
    try:
        return _PRESET_TO_BWA_OPTS[preset]
    except KeyError:
        raise ValueError(
            f"unknown preset {preset!r}; supported presets: "
            f"{', '.join(sorted(_PRESET_TO_BWA_OPTS))}"
        ) from None


def _get_bwa_aligner(reference_fasta: str, preset: str) -> BwaAligner:
    _ensure_bwa_index(reference_fasta)
    options = _bwa_options_for_preset(preset)
    try:
        return BwaAligner(reference_fasta, options=options)
    except ValueError as e:
        raise RuntimeError(f"failed to load/index reference: {reference_fasta}") from e


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

    aligner = _get_bwa_aligner(reference_fasta, preset)

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
        for name, seq, qual in _fastx_read(str(reads_path)):
            hits = aligner.align_seq(seq)
            if not hits:
                continue
            hit = hits[0]  # bwa mem returns its best-scoring hit first
            if hit.mapq < min_mapq:
                continue

            cigar_tuples = _parse_bwa_cigar(hit.cigar)

            a = pysam.AlignedSegment(header=bam_out.header)
            a.query_name = name
            a.flag = 16 if hit.orient == "-" else 0
            a.reference_id = bam_out.get_tid(hit.rname)
            a.reference_start = hit.pos
            a.mapping_quality = hit.mapq
            a.cigar = cigar_tuples

            query_seq = seq if hit.orient == "+" else _revcomp(seq)
            a.query_sequence = query_seq
            if qual:
                q = qual if hit.orient == "+" else qual[::-1]
                a.query_qualities = pysam.qualitystring_to_array(q)
            else:
                a.query_qualities = [30] * len(query_seq)  # no qualities

            ref_end = hit.pos + _ref_span(cigar_tuples)
            ref_region = ref_fasta.fetch(hit.rname, hit.pos, ref_end)

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
    hits = aligner.align_seq(seq)
    if not hits:
        return None
    hit = hits[0]  # bwa mem returns its best-scoring hit first
    if hit.mapq < min_mapq:
        return None
    return hit


def _build_segment(bam_out, ref_fasta, name, seq, qual, hit, technology_tag):
    cigar_tuples = _parse_bwa_cigar(hit.cigar)

    a = pysam.AlignedSegment(header=bam_out.header)
    a.query_name = name
    a.flag = 16 if hit.orient == "-" else 0
    a.reference_id = bam_out.get_tid(hit.rname)
    a.reference_start = hit.pos
    a.mapping_quality = hit.mapq
    a.cigar = cigar_tuples

    query_seq = seq if hit.orient == "+" else _revcomp(seq)
    a.query_sequence = query_seq
    if qual:
        q = qual if hit.orient == "+" else qual[::-1]
        a.query_qualities = pysam.qualitystring_to_array(q)
    else:
        a.query_qualities = [30] * len(query_seq)

    ref_end = hit.pos + _ref_span(cigar_tuples)
    ref_region = ref_fasta.fetch(hit.rname, hit.pos, ref_end)
    md, nm = _build_md_and_nm(ref_region, query_seq, 0, cigar_tuples)
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

    aligner = _get_bwa_aligner(reference_fasta, preset)

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
            _fastx_read(str(reads_r1)), _fastx_read(str(reads_r2))
        ):
            name = n1.split()[0]
            if name.endswith(("/1", "/2")):
                name = name[:-2]

            h1 = _map_one(aligner, s1, min_mapq)
            h2 = _map_one(aligner, s2, min_mapq)

            if h1 is None and h2 is None:
                continue

            a1 = _build_segment(bam_out, ref_fasta, name, s1, q1, h1, technology_tag) if h1 else None
            a2 = _build_segment(bam_out, ref_fasta, name, s2, q2, h2, technology_tag) if h2 else None

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
    reference_fasta = str(reference_fasta)
    aligner = _get_bwa_aligner(reference_fasta, preset)

    hits = aligner.align_seq(sequence)
    if not hits:
        return SingleReadAlignment(read_name=read_name, mapped=False)
    hit = hits[0]
    
    strand = hit.orient
    query_seq = sequence if hit.orient == "+" else _revcomp(sequence)
    query_qual = quality if quality is None else (quality if hit.orient == "+" else quality[::-1])
    cigar_tuples = _parse_bwa_cigar(hit.cigar)
    cigar_string = hit.cigar

    ref_end = hit.pos + _ref_span(cigar_tuples)
    ref_fasta = pysam.FastaFile(reference_fasta)
    ref_region = ref_fasta.fetch(hit.rname, hit.pos, ref_end)
    ref_fasta.close()

    md, nm, raw_mismatches = _build_md_nm_mismatches(ref_region, query_seq, 0, cigar_tuples)

    mismatches = []
    for m in raw_mismatches:
        base_qual = None
        if query_qual is not None:
            base_qual = ord(query_qual[m["query_offset"]]) - 33
        mismatches.append({
            "ref_pos": hit.pos + m["ref_offset"] + 1,  # 1-based
            "ref_base": m["ref_base"],
            "read_base": m["read_base"],
            "base_qual": base_qual,
        })

    return SingleReadAlignment(
        read_name=read_name,
        mapped=True,
        chrom=hit.rname,
        ref_start=hit.pos,
        ref_end=ref_end,
        strand=strand,
        mapping_qual=hit.mapq,
        cigar_string=cigar_string,
        cigar_tuples=cigar_tuples,
        md_tag=md,
        nm=nm,
        mismatches=mismatches,
    )