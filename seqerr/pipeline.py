"""
pipeline.py

CLI entry point exposing individual pipeline steps: map, extract, compare, bqsr.

"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .align import map_reads_to_bam
from .bam_reader import iter_mismatches, write_mismatches_csv, Mismatch, CSV_FIELDS
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
    n = map_reads_to_bam(
        reads_path=args.reads,
        reference_fasta=args.reference,
        out_bam=out_path,
        preset=args.preset,
        technology_tag=args.technology,
        min_mapq=args.min_mapq,
    )
    print(f"[map] {n} reads mapped -> {out_path} (+ .bai index)")


def cmd_extract(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mismatches = list(
        iter_mismatches(
            bam_path=args.bam,
            technology=args.technology,
            reference_fasta=args.reference,
            min_mapq=args.min_mapq,
            min_base_qual=args.min_base_qual,
        )
    )
    print(f"[extract] {len(mismatches)} raw mismatches found")

    if args.vcf:
        variant_positions = load_variant_positions(args.vcf)
        errors, variants = classify_mismatches(mismatches, variant_positions)
    else:
        errors, variants = mismatches, []

    write_mismatches_csv(iter(errors), out_dir / "errors.csv")
    write_mismatches_csv(iter(variants), out_dir / "variants.csv")
    print(f"[extract] {len(errors)} sequencing errors -> {out_dir/'errors.csv'}")
    print(f"[extract] {len(variants)} known variants   -> {out_dir/'variants.csv'}")


def cmd_compare(args: argparse.Namespace) -> None:
    errors_a = _read_mismatches_csv(args.errors_a)
    errors_b = _read_mismatches_csv(args.errors_b)
    for m in errors_a:
        m.technology = args.label_a
    for m in errors_b:
        m.technology = args.label_b

    dists = {
        args.label_a: QualityDistribution.from_mismatches(args.label_a, errors_a),
        args.label_b: QualityDistribution.from_mismatches(args.label_b, errors_b),
    }
    summary = summarize_comparison(dists)
    print(summary)
    if args.out:
        Path(args.out).write_text(summary + "\n")
        for tech, dist in dists.items():
            dist.write_histogram_csv(Path(args.out).with_name(f"{tech}_qual_histogram.csv"))


def cmd_bqsr(args: argparse.Namespace) -> None:
    variant_positions = (
        load_variant_positions(args.vcf) if args.vcf else set()
    )

    before_all = list(
        iter_mismatches(args.bam_before, args.technology, args.reference,
                         min_mapq=args.min_mapq, min_base_qual=0)
    )
    after_all = list(
        iter_mismatches(args.bam_after, args.technology, args.reference,
                         min_mapq=args.min_mapq, min_base_qual=0)
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
    pm.add_argument("--preset", default="sr", help="minimap2 preset: sr (short reads, default), map-ont, map-pb, ...")
    pm.add_argument("--min-mapq", type=int, default=0)
    pm.set_defaults(func=cmd_map)

    pe = sub.add_parser("extract", help="Extract mismatches from a BAM, split errors vs variants")
    pe.add_argument("--bam", required=True)
    pe.add_argument("--technology", required=True, help='e.g. "illumina" or "bgi"')
    pe.add_argument("--reference", help="reference FASTA, needed only if BAM lacks MD tags")
    pe.add_argument("--vcf", help="VCF/BCF of called variants at these samples")
    pe.add_argument("--out-dir", required=True)
    pe.add_argument("--min-mapq", type=int, default=1)
    pe.add_argument("--min-base-qual", type=int, default=0)
    pe.set_defaults(func=cmd_extract)

    pc = sub.add_parser("compare", help="Compare error quality distributions of two technologies")
    pc.add_argument("--errors-a", required=True, help="errors.csv from `extract` for technology A")
    pc.add_argument("--errors-b", required=True, help="errors.csv from `extract` for technology B")
    pc.add_argument("--label-a", default="A")
    pc.add_argument("--label-b", default="B")
    pc.add_argument("--out", help="write summary text + per-tech histogram CSVs")
    pc.set_defaults(func=cmd_compare)

    pb = sub.add_parser("bqsr", help="Compare base qualities at error sites before vs after BQSR")
    pb.add_argument("--bam-before", required=True)
    pb.add_argument("--bam-after", required=True)
    pb.add_argument("--technology", required=True)
    pb.add_argument("--reference", help="reference FASTA, needed only if BAMs lack MD tags")
    pb.add_argument("--vcf", help="VCF/BCF of called variants, to exclude real variants")
    pb.add_argument("--min-mapq", type=int, default=1)
    pb.add_argument("--out", required=True)
    pb.set_defaults(func=cmd_bqsr)

    return p


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
