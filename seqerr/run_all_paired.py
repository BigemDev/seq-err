"""
run_all_paired.py

"""
from __future__ import annotations
 
import argparse
from pathlib import Path
from typing import Optional

from .vcf_reader import load_variant_positions
from .run_all import (  # noqa: F401
    TechnologyResult,
    add_threading_args,
    _extract_and_classify,
    _intersect_vcfs,
    _maybe_call_variants,
    _prepare_reference,
    _write_summary,
    ensure_bwa_index,
    map_paired_bwa_mem,
    plan_intervals,
    read_fai,
    run_haplotypecaller_parallel,
)


def _get_or_map_bam_paired(
    label: str,
    out_dir: Path,
    reads_1: Optional[str],
    reads_2: Optional[str],
    bam: Optional[str],
    reference: Optional[str],
    threads: int,
) -> Path:
    if bam:
        return Path(bam)
    if not reads_1 or not reads_2:
        raise ValueError(
            f"Got R1={reads_1!r}, R2={reads_2!r}"
        )
    if not reference:
        raise ValueError(f"[{label}] --reference is required when mapping from reads")
 
    bam_path = out_dir / label / f"{label}.sorted.bam"
    n = map_paired_bwa_mem(
        reads_r1=reads_1,
        reads_r2=reads_2,
        reference_fasta=reference,
        out_bam=bam_path,
        label=label,
        threads=threads,
    )
    return bam_path
 
 
def run_all_paired(
    out_dir: str | Path,
    label_a: str,
    label_b: str,
    reads_a1: Optional[str] = None,
    reads_a2: Optional[str] = None,
    reads_b1: Optional[str] = None,
    reads_b2: Optional[str] = None,
    bam_a: Optional[str] = None,
    bam_b: Optional[str] = None,
    reference: Optional[str] = None,
    vcf: Optional[str] = None,
    min_mapq: int = 1,
    min_base_qual: int = 0,
    auto_gatk: bool = False,
    threads: int = 4,
    hc_chunk_mb: int = 25,
    gatk_java_options: Optional[str] = None,
    extract_window_mb: int = 5,
) -> tuple[TechnologyResult, TechnologyResult]:

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    variant_positions = load_variant_positions(vcf, threads=threads) if vcf else set()

    bam_path_a = _get_or_map_bam_paired(
        label_a, out_dir, reads_a1, reads_a2, bam_a, reference, threads)
    bam_path_b = _get_or_map_bam_paired(
        label_b, out_dir, reads_b1, reads_b2, bam_b, reference, threads)

    variant_positions = _maybe_call_variants(
        auto_gatk, reference, out_dir, label_a, label_b, bam_path_a, bam_path_b,
        variant_positions, threads, hc_chunk_mb, gatk_java_options)

    result_a = _extract_and_classify(label_a, bam_path_a, out_dir, reference, variant_positions,
                                     min_mapq, min_base_qual, threads, extract_window_mb)
    result_b = _extract_and_classify(label_b, bam_path_b, out_dir, reference, variant_positions,
                                     min_mapq, min_base_qual, threads, extract_window_mb)

    _write_summary(out_dir, label_a, label_b, result_a, result_b)
    return result_a, result_b


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="seqerr.run_all_paired", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reads-a1", help="technology A: R1 reads (FASTA/FASTQ, .gz OK)")
    p.add_argument("--reads-a2", help="technology A: R2 reads (FASTA/FASTQ, .gz OK)")
    p.add_argument("--reads-b1", help="technology B: R1 reads (FASTA/FASTQ, .gz OK)")
    p.add_argument("--reads-b2", help="technology B: R2 reads (FASTA/FASTQ, .gz OK)")
    p.add_argument("--bam-a", help="technology A: already-mapped, indexed BAM (skip mapping)")
    p.add_argument("--bam-b", help="technology B: already-mapped, indexed BAM (skip mapping)")
    p.add_argument("--label-a", default="illumina")
    p.add_argument("--label-b", default="bgi")
    p.add_argument("--reference", help="reference FASTA (required if mapping from reads)")
    p.add_argument("--vcf", help="VCF/BCF of called variants, to separate variants from errors")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-mapq", type=int, default=1)
    p.add_argument("--min-base-qual", type=int, default=0)
    p.add_argument("--auto-gatk", action="store_true", help="Auto gen VCF with GATK")
    add_threading_args(p)
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    run_all_paired(
        out_dir=args.out_dir,
        label_a=args.label_a,
        label_b=args.label_b,
        reads_a1=args.reads_a1,
        reads_a2=args.reads_a2,
        reads_b1=args.reads_b1,
        reads_b2=args.reads_b2,
        bam_a=args.bam_a,
        bam_b=args.bam_b,
        reference=args.reference,
        vcf=args.vcf,
        min_mapq=args.min_mapq,
        min_base_qual=args.min_base_qual,
        auto_gatk=args.auto_gatk,
        threads=args.threads,
        hc_chunk_mb=args.hc_chunk_mb,
        gatk_java_options=args.gatk_java_options,
        extract_window_mb=args.extract_window_mb,
    )


if __name__ == "__main__":
    main()