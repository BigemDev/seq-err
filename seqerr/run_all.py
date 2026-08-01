"""
run_all.py

"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .align import map_reads_to_bam
from .bam_reader import iter_mismatches, write_mismatches_csv
from .vcf_reader import load_variant_positions
from .mismatch_classify import classify_mismatches
from .quality_stats import QualityDistribution, summarize_comparison


@dataclass
class TechnologyResult:
    label: str
    bam_path: Path
    n_mismatches: int
    n_errors: int
    n_variants: int
    errors_csv: Path
    variants_csv: Path
    distribution: QualityDistribution


def _get_or_map_bam(
    label: str,
    out_dir: Path,
    reads: Optional[str],
    bam: Optional[str],
    reference: Optional[str],
    preset: str,
    min_mapq: int,
) -> Path:
    if bam:
        return Path(bam)
    if not reads:
        raise ValueError(f"[{label}] need either --reads-{label[0]} or --bam-{label[0]}")
    if not reference:
        raise ValueError(f"[{label}] --reference is required when mapping from reads")

    bam_path = out_dir / label / f"{label}.sorted.bam"
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    n = map_reads_to_bam(
        reads_path=reads,
        reference_fasta=reference,
        out_bam=bam_path,
        preset=preset,
        technology_tag=label,
        min_mapq=min_mapq,
    )
    print(f"[{label}] mapped {n} reads -> {bam_path}")
    return bam_path


def _extract_and_classify(
    label: str,
    bam_path: Path,
    out_dir: Path,
    reference: Optional[str],
    variant_positions: set,
    min_mapq: int,
    min_base_qual: int,
) -> TechnologyResult:
    mismatches = list(
        iter_mismatches(
            bam_path=bam_path,
            technology=label,
            reference_fasta=reference,
            min_mapq=min_mapq,
            min_base_qual=min_base_qual,
        )
    )
    if variant_positions:
        errors, variants = classify_mismatches(mismatches, variant_positions)
    else:
        errors, variants = mismatches, []

    tech_dir = out_dir / label
    tech_dir.mkdir(parents=True, exist_ok=True)
    errors_csv = tech_dir / "errors.csv"
    variants_csv = tech_dir / "variants.csv"
    write_mismatches_csv(iter(errors), errors_csv)
    write_mismatches_csv(iter(variants), variants_csv)

    dist = QualityDistribution.from_mismatches(label, errors)
    dist.write_histogram_csv(tech_dir / "error_qual_histogram.csv")

    print(f"[{label}] {len(mismatches)} raw mismatches "
          f"-> {len(errors)} errors, {len(variants)} known variants")

    return TechnologyResult(
        label=label,
        bam_path=bam_path,
        n_mismatches=len(mismatches),
        n_errors=len(errors),
        n_variants=len(variants),
        errors_csv=errors_csv,
        variants_csv=variants_csv,
        distribution=dist,
    )


def run_all(
    out_dir: str | Path,
    label_a: str,
    label_b: str,
    reads_a: Optional[str] = None,
    reads_b: Optional[str] = None,
    bam_a: Optional[str] = None,
    bam_b: Optional[str] = None,
    reference: Optional[str] = None,
    vcf: Optional[str] = None,
    preset: str = "sr",
    min_mapq: int = 1,
    min_base_qual: int = 0,
) -> tuple[TechnologyResult, TechnologyResult]:

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    variant_positions = load_variant_positions(vcf) if vcf else set()
    if vcf:
        print(f"[vcf] loaded {len(variant_positions)} known variant alleles from {vcf}")

    bam_path_a = _get_or_map_bam(label_a, out_dir, reads_a, bam_a, reference, preset, min_mapq)
    bam_path_b = _get_or_map_bam(label_b, out_dir, reads_b, bam_b, reference, preset, min_mapq)

    result_a = _extract_and_classify(label_a, bam_path_a, out_dir, reference,
                                      variant_positions, min_mapq, min_base_qual)
    result_b = _extract_and_classify(label_b, bam_path_b, out_dir, reference,
                                      variant_positions, min_mapq, min_base_qual)

    summary = summarize_comparison({label_a: result_a.distribution, label_b: result_b.distribution})
    summary_path = out_dir / "comparison_summary.txt"
    summary_path.write_text(summary + "\n")
    print()
    print(summary)
    print(f"\n[done] full report written under {out_dir}/")

    return result_a, result_b


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="seqerr.run_all", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reads-a", help="technology A: raw FASTA/FASTQ reads")
    p.add_argument("--reads-b", help="technology B: raw FASTA/FASTQ reads")
    p.add_argument("--bam-a", help="technology A: already-mapped, indexed BAM (skip mapping)")
    p.add_argument("--bam-b", help="technology B: already-mapped, indexed BAM (skip mapping)")
    p.add_argument("--label-a", default="illumina")
    p.add_argument("--label-b", default="bgi")
    p.add_argument("--reference", help="reference FASTA (required if mapping from reads)")
    p.add_argument("--vcf", help="VCF/BCF of called variants, to separate variants from errors")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--preset", default="sr", help="minimap2 preset for mapping (default: sr)")
    p.add_argument("--min-mapq", type=int, default=1)
    p.add_argument("--min-base-qual", type=int, default=0)
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    run_all(
        out_dir=args.out_dir,
        label_a=args.label_a,
        label_b=args.label_b,
        reads_a=args.reads_a,
        reads_b=args.reads_b,
        bam_a=args.bam_a,
        bam_b=args.bam_b,
        reference=args.reference,
        vcf=args.vcf,
        preset=args.preset,
        min_mapq=args.min_mapq,
        min_base_qual=args.min_base_qual,
    )


if __name__ == "__main__":
    main()
