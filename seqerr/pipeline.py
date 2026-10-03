"""
pipeline.py

CLI entry point exposing individual pipeline steps: map, extract, compare, bqsr.

"""
from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
from collections import Counter
from pathlib import Path

from .run_all import map_paired_bwa_mem
from .bam_reader import (iter_mismatches_parallel, extract_errors_parallel,
                         write_mismatches_csv, Mismatch, CSV_FIELDS)
from .vcf_reader import load_variant_positions
from .mismatch_classify import classify_mismatches
from .bqsr_compare import compare_bqsr, dropped_after_bqsr, write_bqsr_deltas_csv
from .quality_stats import QualityDistribution, summarize_comparison


def _read_mismatches_csv(path: str | Path) -> list[Mismatch]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        out = []
        for row in reader:
            row["ref_pos"] = int(row["ref_pos"])
            row["base_qual"] = int(row["base_qual"])
            row["mapping_qual"] = int(row["mapping_qual"])
            row["is_read1"] = row["is_read1"] == "True"
            row["is_reverse"] = row["is_reverse"] == "True"
            out.append(Mismatch(**row))
        return out


def cmd_map(args: argparse.Namespace) -> None:
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = map_paired_bwa_mem(
        reads_r1=args.reads,
        reads_r2=args.reads2,
        reference_fasta=args.reference,
        out_bam=out_path,
        label=args.technology,
        threads=args.threads,
        preset=args.preset,
    )
    

def cmd_extract(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    variant_positions = (
        load_variant_positions(args.vcf, threads=args.threads) if args.vcf else set()
    )
    res = extract_errors_parallel(
        bam_path=args.bam,
        technology=args.technology,
        errors_csv=out_dir / "errors.csv",
        variants_csv=out_dir / "variants.csv",
        reference_fasta=args.reference,
        variant_positions=variant_positions,
        min_mapq=args.min_mapq,
        min_base_qual=args.min_base_qual,
        threads=args.threads,
        window_size=args.window_mb * 1_000_000,
    )
    print(f"[extract] {res.n_mismatches} raw mismatches found")
    print(f"[extract] {res.n_errors} sequencing errors -> {out_dir/'errors.csv'}")
    print(f"[extract] {res.n_variants} known variants   -> {out_dir/'variants.csv'}")


def _qual_histogram_from_csv(path: str | Path) -> Counter:
    hist: Counter = Counter()
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        col = header.index("base_qual")
        for row in reader:
            hist[int(row[col])] += 1
    return hist


def cmd_compare(args: argparse.Namespace) -> None:
    paths = [args.errors_a, args.errors_b]
    if args.threads > 1:
        with mp.get_context().Pool(2) as pool:   # the two files are independent
            hist_a, hist_b = pool.map(_qual_histogram_from_csv, paths)
    else:
        hist_a, hist_b = (_qual_histogram_from_csv(p) for p in paths)

    dists = {
        args.label_a: QualityDistribution.from_histogram(args.label_a, hist_a),
        args.label_b: QualityDistribution.from_histogram(args.label_b, hist_b),
    }
    summary = summarize_comparison(dists)
    print(summary)
    if args.out:
        Path(args.out).write_text(summary + "\n")
        for tech, dist in dists.items():
            dist.write_histogram_csv(Path(args.out).with_name(f"{tech}_qual_histogram.csv"))


def cmd_bqsr(args: argparse.Namespace) -> None:
    variant_positions = (
        load_variant_positions(args.vcf, threads=args.threads) if args.vcf else set()
    )

    before_all = list(
        iter_mismatches_parallel(args.bam_before, args.technology, args.reference,
                                 min_mapq=args.min_mapq, min_base_qual=0,
                                 threads=args.threads)
    )
    after_all = list(
        iter_mismatches_parallel(args.bam_after, args.technology, args.reference,
                                 min_mapq=args.min_mapq, min_base_qual=0,
                                 threads=args.threads)
    )

    if variant_positions:
        before_errors, _ = classify_mismatches(before_all, variant_positions)
        after_errors, _ = classify_mismatches(after_all, variant_positions)
    else:
        before_errors, after_errors = before_all, after_all

    deltas = compare_bqsr(before_errors, after_errors)
    dropped = dropped_after_bqsr(before_errors, after_errors)

    print(f"[bqsr] {len(deltas)} error positions present before AND after BQSR")
    print(f"[bqsr] {len(dropped)} error positions present before but not after BQSR")
    if deltas:
        mean_delta = sum(d.delta for d in deltas) / len(deltas)
        print(f"[bqsr] mean quality delta (after - before): {mean_delta:+.3f}")

    write_bqsr_deltas_csv(deltas, args.out)
    print(f"[bqsr] wrote {args.out}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="seqerr", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    pm = sub.add_parser("map", help="Map raw FASTA/FASTQ reads to a reference -> sorted+indexed BAM")
    pm.add_argument("--reads", required=True, help="FASTA or FASTQ (.gz OK) reads file")
    pm.add_argument("--reference", required=True, help="reference FASTA")
    pm.add_argument("--out", required=True, help="output BAM path (index written alongside)")
    pm.add_argument("--technology", required=True, help='e.g. "illumina" or "bgi", stored as read-group')
    pm.add_argument("--preset", default="sr", help="mapping preset: sr (short reads, default), map-ont, map-pb, map-hifi, intractg")
    pm.add_argument("--reads2", help="R2 reads for paired-end mapping (--reads is then R1)")
    pm.add_argument("--min-mapq", type=int, default=0, help="kept for compatibility; filter at extract time")
    pm.add_argument("--threads", type=int, default=4, help="bwa mem / samtools threads")
    pm.set_defaults(func=cmd_map)

    pe = sub.add_parser("extract", help="Extract mismatches from a BAM, split errors vs variants")
    pe.add_argument("--bam", required=True)
    pe.add_argument("--technology", required=True, help='e.g. "illumina" or "bgi"')
    pe.add_argument("--reference", help="reference FASTA, needed only if BAM lacks MD tags")
    pe.add_argument("--vcf", help="VCF/BCF of called variants at these samples")
    pe.add_argument("--out-dir", required=True)
    pe.add_argument("--min-mapq", type=int, default=1)
    pe.add_argument("--min-base-qual", type=int, default=0)
    pe.add_argument("--threads", type=int, default=4, help="processes scanning the BAM in parallel")
    pe.add_argument("--window-mb", type=int, default=5, help="genomic window size in Mb (default: 5)")
    pe.set_defaults(func=cmd_extract)

    pc = sub.add_parser("compare", help="Compare error quality distributions of two technologies")
    pc.add_argument("--errors-a", required=True, help="errors.csv from `extract` for technology A")
    pc.add_argument("--errors-b", required=True, help="errors.csv from `extract` for technology B")
    pc.add_argument("--label-a", default="A")
    pc.add_argument("--label-b", default="B")
    pc.add_argument("--out", help="write summary text + per-tech histogram CSVs")
    pc.add_argument("--threads", type=int, default=2, help="the two CSVs are read in parallel when > 1")
    pc.set_defaults(func=cmd_compare)

    pb = sub.add_parser("bqsr", help="Compare base qualities at error sites before vs after BQSR")
    pb.add_argument("--bam-before", required=True)
    pb.add_argument("--bam-after", required=True)
    pb.add_argument("--technology", required=True)
    pb.add_argument("--reference", help="reference FASTA, needed only if BAMs lack MD tags")
    pb.add_argument("--vcf", help="VCF/BCF of called variants, to exclude real variants")
    pb.add_argument("--min-mapq", type=int, default=1)
    pb.add_argument("--out", required=True)
    pb.add_argument("--threads", type=int, default=4, help="processes scanning each BAM in parallel")
    pb.set_defaults(func=cmd_bqsr)

    return p


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()