"""
bam_reader.py

Reads a BAM file and extracts per-base mismatches (position, bases, quality) 
between aligned reads and the reference.

"""
from __future__ import annotations

import csv
import multiprocessing as mp
import shutil
from collections import Counter, defaultdict
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
    depth: int = 0                # read depth at this position (filled later)

    def as_row(self) -> dict:
        return asdict(self)


CSV_FIELDS = list(Mismatch.__dataclass_fields__.keys())


def _min_mapq_ok(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return (
        not read.is_unmapped
        and not read.is_secondary
        and not read.is_duplicate
        and read.mapping_quality >= min_mapq
    )


_CIGAR_ALIGNED = (0, 7, 8)   # M, =, X : query and reference both advance
_CIGAR_QUERY_ONLY = (1, 4)    # I, S
_CIGAR_REF_ONLY = (2, 3)      # D, N


def _read_mismatches(
    read: pysam.AlignedSegment,
    ref_seq: Optional[str],
    ref_offset: int,
    technology: str,
    min_base_qual: int,
) -> Iterator[Mismatch]:
    seq = read.query_sequence
    quals = read.query_qualities
    chrom = read.reference_name

    if ref_seq is None:
        pairs = ((q, r, b) for q, r, b in read.get_aligned_pairs(with_seq=True, matches_only=False)
                 if q is not None and r is not None)
    else:
        pairs = None

    def make(qpos, rpos, ref_base):
        read_base = seq[qpos].upper()
        if ref_base == read_base or ref_base not in "ACGT" or not ref_base:
            return None
        bq = quals[qpos]
        if bq < min_base_qual:
            return None
        return Mismatch(
            read_name=read.query_name,
            technology=technology,
            chrom=chrom,
            ref_pos=rpos + 1,  # report 1-based
            ref_base=ref_base,
            read_base=read_base,
            base_qual=bq,
            mapping_qual=read.mapping_quality,
            cigar=read.cigarstring or "",
            is_read1=read.is_read1,
            is_reverse=read.is_reverse,
        )

    if pairs is not None:
        for qpos, rpos, rb in pairs:
            m = make(qpos, rpos, rb.upper())
            if m is not None:
                yield m
        return

    qpos = 0
    rpos = read.reference_start
    useq = seq.upper()
    for op, length in read.cigartuples or ():
        if op in _CIGAR_ALIGNED:
            rslice = ref_seq[rpos - ref_offset: rpos - ref_offset + length].upper()
            if rslice != useq[qpos:qpos + length]:
                for i in range(length):
                    m = make(qpos + i, rpos + i, rslice[i])
                    if m is not None:
                        yield m
            qpos += length
            rpos += length
        elif op in _CIGAR_QUERY_ONLY:
            qpos += length
        elif op in _CIGAR_REF_ONLY:
            rpos += length


def _scan_region(
    bam: pysam.AlignmentFile,
    fasta: Optional[pysam.FastaFile],
    technology: str,
    chrom: Optional[str],
    start: Optional[int],
    end: Optional[int],
    min_mapq: int,
    min_base_qual: int,
    own_reads_only: bool,
) -> Iterator[Mismatch]:
    reads = bam.fetch(chrom, start, end) if chrom is not None else bam.fetch()
    for read in reads:
        if own_reads_only and not (start <= read.reference_start < end):
            continue
        if not _min_mapq_ok(read, min_mapq):
            continue
        if read.query_sequence is None or read.query_qualities is None:
            continue
        ref_seq = None
        if fasta is not None:
            ref_seq = fasta.fetch(read.reference_name, read.reference_start,
                                  read.reference_end)
        yield from _read_mismatches(read, ref_seq, read.reference_start,
                                    technology, min_base_qual)


def iter_mismatches(
    bam_path: str | Path,
    technology: str,
    reference_fasta: Optional[str | Path] = None,
    regions: Optional[list[tuple[str, int, int]]] = None,
    min_mapq: int = 0,
    min_base_qual: int = 0,
) -> Iterator[Mismatch]:
    fasta = pysam.FastaFile(str(reference_fasta)) if reference_fasta else None
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            if regions:
                for c, s, e in regions:
                    yield from _scan_region(bam, fasta, technology, c, s, e,
                                            min_mapq, min_base_qual, False)
            else:
                yield from _scan_region(bam, fasta, technology, None, None, None,
                                        min_mapq, min_base_qual, False)
    finally:
        if fasta is not None:
            fasta.close()


def _row(m: Mismatch) -> list:
    return [m.read_name, m.technology, m.chrom, m.ref_pos, m.ref_base, m.read_base,
            m.base_qual, m.mapping_qual, m.cigar, m.is_read1, m.is_reverse,m.depth]


@dataclass
class ExtractionResult:
    n_mismatches: int
    n_errors: int
    n_variants: int
    error_hist: Counter
    errors_csv: Path
    variants_csv: Path
    weighted_error_hist: dict[int, float] | None = None

    quality_depth_distribution: dict[
        tuple[int, int], float
    ] | None = None


_W: dict = {}


def _init_worker(cfg: dict) -> None:
    _W.clear()
    _W.update(cfg)


def _worker_handles():
    if "bam" not in _W:
        _W["bam"] = pysam.AlignmentFile(str(_W["bam_path"]), "rb")
        ref = _W["reference"]
        _W["fasta"] = pysam.FastaFile(str(ref)) if ref else None
    return _W["bam"], _W["fasta"]


def _process_window(task: tuple) -> tuple:
    idx, chrom, start, end = task
    bam, fasta = _worker_handles()
    variants = _W["variant_positions"]
    tmp = Path(_W["tmp_dir"])
    err_path = tmp / f"errors_{idx:06d}.part"
    var_path = tmp / f"variants_{idx:06d}.part"

    hist: Counter = Counter()
    depth_hist: Counter = Counter()
    mapq_hist: Counter = Counter()
    quality_by_depth = Counter()
    total_by_depth = Counter()
    weighted_quality: dict[int, float] = defaultdict(float)
    n_all = n_err = n_var = 0

    with open(err_path, "w", newline="") as fe, open(var_path, "w", newline="") as fv:
        we, wv = csv.writer(fe), csv.writer(fv)
        mismatches = _scan_region(
            bam, fasta, _W["technology"], chrom, start, end,
            _W["min_mapq"], _W["min_base_qual"], True,
        )

        if _W["apply_qc"]:
            # Local import avoids bam_reader <-> mismatch_classify circular imports.
            from .mismatch_classify import process_technology

            # Only one bounded window is materialized in a worker.
            window_mismatches = list(mismatches)
            n_all = len(window_mismatches)
            errors, variant_records = process_technology(
                window_mismatches, variants, str(_W["bam_path"]),
                min_alignment_fraction=_W["min_alignment_fraction"],
                min_depth=_W["min_depth"],
                min_mapq=_W["min_mapq"],
                plot_histograms=False,  # plot ONCE after all windows merge
            )

            for m in errors:
                we.writerow(_row(m))
                n_err += 1
                hist[m.base_qual] += 1
                depth_hist[m.depth] += 1
                mapq_hist[m.mapping_qual] += 1
                if m.depth > 0:
                    weighted_quality[m.base_qual] += 1.0 / m.depth

                depth_bin = (m.depth // 10) * 10
                quality_bin = (m.base_qual // 10) * 10
                quality_by_depth[(depth_bin, quality_bin)] += 1
                total_by_depth[depth_bin] += 1

            for m in variant_records:
                wv.writerow(_row(m))
                n_var += 1

        else:
            # Original streaming extraction path, unchanged unless QC enabled.
            for m in mismatches:
                n_all += 1
                if variants and (m.chrom, m.ref_pos, m.ref_base, m.read_base) in variants:
                    wv.writerow(_row(m))
                    n_var += 1
                else:
                    we.writerow(_row(m))
                    hist[m.base_qual] += 1
                    n_err += 1

    return (idx, n_all, n_err, n_var, hist,
            str(err_path), str(var_path),
            depth_hist, mapq_hist,
            dict(weighted_quality),
            quality_by_depth,
            total_by_depth,
        )


def plan_windows(bam_path: str | Path, window_size: int) -> list[tuple[int, str, int, int]]:
    if window_size < 1:
        raise ValueError("window_size must be >= 1")
    windows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        mapped = {s.contig: s.mapped for s in bam.get_index_statistics()}
        for name, length in zip(bam.references, bam.lengths):
            if not mapped.get(name):
                continue
            for start in range(0, length, window_size):
                windows.append((len(windows), name, start, min(start + window_size, length)))
    return windows


def extract_errors_parallel(
    bam_path: str | Path,
    technology: str,
    errors_csv: str | Path,
    variants_csv: str | Path,
    reference_fasta: Optional[str | Path] = None,
    variant_positions: Optional[set] = None,
    min_mapq: int = 0,
    min_base_qual: int = 0,
    threads: int = 1,
    window_size: int = 5_000_000,
    apply_qc: bool = False,
    min_alignment_fraction: float = 0.8,
    min_depth: int = 10,
    histograms_dir: str = "histograms",
) -> ExtractionResult:
    errors_csv, variants_csv = Path(errors_csv), Path(variants_csv)
    errors_csv.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = errors_csv.parent / f".{errors_csv.stem}_parts"
    shutil.rmtree(tmp_dir, ignore_errors=True)
    tmp_dir.mkdir(parents=True)

    windows = plan_windows(bam_path, window_size)
    cfg = dict(bam_path=str(bam_path), reference=str(reference_fasta) if reference_fasta else None,
               technology=technology, variant_positions=variant_positions or set(),
               min_mapq=min_mapq, min_base_qual=min_base_qual, tmp_dir=str(tmp_dir),
               apply_qc=apply_qc, min_alignment_fraction=min_alignment_fraction,
               min_depth=min_depth)

    results: dict[int, tuple] = {}
    try:
        if threads <= 1 or len(windows) <= 1:
            _init_worker(cfg)
            for w in windows:
                r = _process_window(w)
                results[r[0]] = r
            _W.clear()
        else:
            try:
                ctx = mp.get_context("fork")
            except ValueError:
                ctx = mp.get_context()
            with ctx.Pool(min(threads, len(windows)), initializer=_init_worker,
                          initargs=(cfg,)) as pool:
                done = 0
                for r in pool.imap_unordered(_process_window, windows):
                    results[r[0]] = r
                    done += 1
                    if done % 25 == 0 or done == len(windows):
                        print(f"[extract] {technology}: {done}/{len(windows)} windows done",
                              flush=True)

        hist: Counter = Counter()
        depth_hist: Counter = Counter()
        mapq_hist: Counter = Counter()
        quality_by_depth_all = Counter()
        total_by_depth_all = Counter()
        weighted_quality: dict[int, float] = defaultdict(float)
        n_all = n_err = n_var = 0
        for out_path, kind in ((errors_csv, 5), (variants_csv, 6)):
            with open(out_path, "w", newline="") as out:
                csv.writer(out).writerow(CSV_FIELDS)
                for idx in sorted(results):
                    part = Path(results[idx][kind])
                    with open(part, "r", newline="") as src:
                        shutil.copyfileobj(src, out, 1 << 20)
                    part.unlink()
        for idx in sorted(results):
            (_, a, e, v, h, _, _, depth_counts, mapq_counts, quality_weights, window_quality_by_depth, window_total_by_depth) = results[idx]
            n_all += a; n_err += e; n_var += v
            hist.update(h)
            depth_hist.update(depth_counts)
            mapq_hist.update(mapq_counts)
            quality_by_depth_all.update(window_quality_by_depth)
            total_by_depth_all.update(window_total_by_depth)
            
            for q, weight in quality_weights.items():
                weighted_quality[q] += weight
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    distribution = {}

    for (depth_bin, quality_bin), count in quality_by_depth_all.items():

        distribution[(depth_bin, quality_bin)] = (count / total_by_depth_all[depth_bin])

    #print(distribution)
       
    if apply_qc:
        from .mismatch_classify import (
            #plot_error_histograms_from_counts,
            plot_quality_depth_heatmap_from_counts,
        )
        """
        plot_error_histograms_from_counts(
            hist, depth_hist, mapq_hist, dict(weighted_quality),
            technology, histograms_dir,
        )
        """
        plot_quality_depth_heatmap_from_counts(
            quality_by_depth_all,
            total_by_depth_all,
            technology,
            histograms_dir,
        )

    return ExtractionResult(n_all, n_err, n_var, hist, errors_csv, variants_csv, weighted_error_hist=(dict(weighted_quality) if apply_qc else None),
        quality_depth_distribution=(distribution if apply_qc else None),
    )


def _collect_window(task: tuple) -> list[Mismatch]:
    idx, chrom, start, end = task
    bam, fasta = _worker_handles()
    return list(_scan_region(bam, fasta, _W["technology"], chrom, start, end,
                             _W["min_mapq"], _W["min_base_qual"], True))


def iter_mismatches_parallel(
    bam_path: str | Path,
    technology: str,
    reference_fasta: Optional[str | Path] = None,
    min_mapq: int = 0,
    min_base_qual: int = 0,
    threads: int = 1,
    window_size: int = 5_000_000,
) -> Iterator[Mismatch]:
    if threads <= 1:
        yield from iter_mismatches(bam_path, technology, reference_fasta,
                                   min_mapq=min_mapq, min_base_qual=min_base_qual)
        return
    windows = plan_windows(bam_path, window_size)
    cfg = dict(bam_path=str(bam_path), reference=str(reference_fasta) if reference_fasta else None,
               technology=technology, variant_positions=set(),
               min_mapq=min_mapq, min_base_qual=min_base_qual, tmp_dir="")
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context()
    with ctx.Pool(min(threads, max(1, len(windows))), initializer=_init_worker,
                  initargs=(cfg,)) as pool:
        for chunk in pool.imap(_collect_window, windows):   # ordered
            yield from chunk
