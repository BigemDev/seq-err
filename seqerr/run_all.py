"""
run_all.py

End-to-end pipeline: maps raw reads for two technologies,
extracts and classifies mismatches, and compares error-quality distributions.

"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pysam

from .bam_reader import extract_errors_parallel
from .vcf_reader import load_variant_positions
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


_BWA_INDEX_EXT = (".amb", ".ann", ".bwt", ".pac", ".sa")


def ensure_bwa_index(reference: str | Path) -> None:
    ref = str(reference)
    if all(Path(ref + ext).exists() for ext in _BWA_INDEX_EXT):
        return
    subprocess.run(["bwa", "index", ref], check=True)


_PRESET_TO_BWA_OPTS = {
    "sr": [],
    "map-ont": ["-x", "ont2d"],
    "map-pb": ["-x", "pacbio"],
    "map-hifi": ["-x", "pacbio"],
    "intractg": ["-x", "intractg"],
}


def map_paired_bwa_mem(
    reads_r1: str | Path,
    reads_r2: Optional[str | Path],
    reference_fasta: str | Path,
    out_bam: str | Path,
    label: str,
    threads: int = 4,
    sort_threads: Optional[int] = None,
    sort_mem: str = "1G",
    platform: str = "ILLUMINA",
    preset: str = "sr",
) -> int:
    """bwa mem | samtools sort | samtools index. reads_r2=None maps single-end."""
    if preset not in _PRESET_TO_BWA_OPTS:
        raise ValueError(f"unknown preset {preset!r}; supported presets: "
                         f"{', '.join(sorted(_PRESET_TO_BWA_OPTS))}")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    ensure_bwa_index(reference_fasta)

    if sort_threads is None:
        sort_threads = max(1, min(8, threads // 2))

    rg = f"@RG\\tID:{label}\\tSM:{label}\\tPL:{platform}"

    bwa_cmd = ["bwa", "mem", "-t", str(threads), "-R", rg,
               *_PRESET_TO_BWA_OPTS[preset], str(reference_fasta), str(reads_r1)]
    if reads_r2:
        bwa_cmd.append(str(reads_r2))
    sort_cmd = ["samtools", "sort", "-@", str(sort_threads), "-m", sort_mem,
                "-T", str(out_bam.with_suffix("")) + ".sorttmp",
                "-o", str(out_bam), "-"]

    bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
    srt = subprocess.Popen(sort_cmd, stdin=bwa.stdout)
    bwa.stdout.close()
    srt_rc = srt.wait()
    bwa_rc = bwa.wait()
    if bwa_rc != 0:
        raise RuntimeError(f"bwa mem failed with exit code {bwa_rc}")
    if srt_rc != 0:
        raise RuntimeError(f"samtools sort failed with exit code {srt_rc}")

    subprocess.run(["samtools", "index", "-@", str(sort_threads), str(out_bam)], check=True)

    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        return sum(s.mapped for s in bam.get_index_statistics())


def read_fai(reference: str | Path) -> list[tuple[str, int]]:
    fai = Path(f"{reference}.fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(reference)], check=True)
    contigs = []
    with open(fai) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            contigs.append((parts[0], int(parts[1])))
    return contigs


def plan_intervals(contigs: list[tuple[str, int]], chunk_size: int) -> list[list[str]]:
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    tasks: list[list[str]] = []
    cur: list[str] = []
    cur_size = 0
    for name, length in contigs:
        if length > chunk_size:
            if cur:
                tasks.append(cur)
                cur, cur_size = [], 0
            for start in range(1, length + 1, chunk_size):
                end = min(start + chunk_size - 1, length)
                tasks.append([f"{name}:{start}-{end}"])
        else:
            cur.append(name)
            cur_size += length
            if cur_size >= chunk_size:
                tasks.append(cur)
                cur, cur_size = [], 0
    if cur:
        tasks.append(cur)
    return tasks


def _run_task(i, intervals, bam, reference, work_dir, java_options):
    shard_vcf = work_dir / f"shard_{i:05d}.vcf.gz"
    log = work_dir / f"shard_{i:05d}.log"
    cmd = ["gatk"]
    if java_options:
        cmd += ["--java-options", java_options]
    cmd += ["HaplotypeCaller", "-R", str(reference), "-I", str(bam),
            "-O", str(shard_vcf), "--native-pair-hmm-threads", "1"]
    for iv in intervals:
        cmd += ["-L", iv]
    with open(log, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"HaplotypeCaller failed on task {i} ({intervals[:3]}...), see {log}")
    return shard_vcf


def run_haplotypecaller_parallel(
    bam_path: str | Path,
    reference: str | Path,
    out_vcf: str | Path,
    threads: int = 4,
    chunk_size: int = 25_000_000,
    java_options: Optional[str] = None,
    keep_tmp: bool = False,
) -> Path:
    out_vcf = Path(out_vcf)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    work_dir = out_vcf.parent / f"{out_vcf.stem}_hc_tmp"
    work_dir.mkdir(parents=True, exist_ok=True)

    tasks = plan_intervals(read_fai(reference), chunk_size)
    print(f"[gatk] {out_vcf.name}: {len(tasks)} interval tasks, {threads} parallel")

    shards: dict[int, Path] = {}
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        futures = {pool.submit(_run_task, i, t, bam_path, reference, work_dir, java_options): i
                   for i, t in enumerate(tasks)}
        done = 0
        try:
            for fut in as_completed(futures):
                shards[futures[fut]] = fut.result()
                done += 1
                print(f"[gatk] {out_vcf.name}: {done}/{len(tasks)} tasks done", flush=True)
        except BaseException:
            for f in futures:
                f.cancel()
            raise

    ordered = [str(shards[i]) for i in sorted(shards)]
    if len(ordered) == 1:
        subprocess.run(["bcftools", "view", "-O", "v", "-o", str(out_vcf), ordered[0]], check=True)
    else:
        subprocess.run(["bcftools", "concat", "-a", "-O", "v", "-o", str(out_vcf)] + ordered,
                       check=True)

    if not keep_tmp:
        shutil.rmtree(work_dir, ignore_errors=True)
    return out_vcf


def _get_or_map_bam(
    label: str,
    out_dir: Path,
    reads: Optional[str],
    bam: Optional[str],
    reference: Optional[str],
    preset: str,
    min_mapq: int,
    threads: int = 1,
) -> Path:
    if bam:
        return Path(bam)
    if not reads:
        raise ValueError(f"[{label}] need either --reads-{label[0]} or --bam-{label[0]}")
    if not reference:
        raise ValueError(f"[{label}] --reference is required when mapping from reads")

    bam_path = out_dir / label / f"{label}.sorted.bam"
    n = map_paired_bwa_mem(
        reads_r1=reads, reads_r2=None, reference_fasta=reference,
        out_bam=bam_path, label=label, threads=threads, preset=preset,
    )
    print(f"[{label}] mapped {n} alignments with bwa mem -> {bam_path}")
    return bam_path


def _extract_and_classify(
    label: str,
    bam_path: Path,
    out_dir: Path,
    reference: Optional[str],
    variant_positions: set,
    min_mapq: int,
    min_base_qual: int,
    threads: int = 1,
    window_mb: int = 5,
) -> TechnologyResult:
    tech_dir = out_dir / label
    tech_dir.mkdir(parents=True, exist_ok=True)
    errors_csv = tech_dir / "errors.csv"
    variants_csv = tech_dir / "variants.csv"

    res = extract_errors_parallel(
        bam_path=bam_path,
        technology=label,
        errors_csv=errors_csv,
        variants_csv=variants_csv,
        reference_fasta=reference,
        variant_positions=variant_positions,
        min_mapq=min_mapq,
        min_base_qual=min_base_qual,
        threads=threads,
        window_size=window_mb * 1_000_000,
    )

    dist = QualityDistribution.from_histogram(label, res.error_hist)
    dist.write_histogram_csv(tech_dir / "error_qual_histogram.csv")

    print(f"[{label}] {res.n_mismatches} raw mismatches "
          f"-> {res.n_errors} errors, {res.n_variants} known variants")

    return TechnologyResult(
        label=label,
        bam_path=bam_path,
        n_mismatches=res.n_mismatches,
        n_errors=res.n_errors,
        n_variants=res.n_variants,
        errors_csv=errors_csv,
        variants_csv=variants_csv,
        distribution=dist,
    )


def _prepare_reference(reference_path: Path):
    # samtools faidx
    if not Path(f"{reference_path}.fai").exists():
        subprocess.run(["samtools", "faidx", str(reference_path)], check=True)

    # gatk CreateSequenceDictionary
    dict_path = reference_path.with_suffix(".dict")
    if not dict_path.exists():
        subprocess.run([
            "gatk", "CreateSequenceDictionary",
            "-R", str(reference_path),
            "-O", str(dict_path)
        ], check=True)


def _run_variant_calling(bam_path: Path, reference: Path, out_vcf: Path, threads: int = 1,
                         chunk_size: int = 25_000_000, java_options: Optional[str] = None):
    run_haplotypecaller_parallel(bam_path, reference, out_vcf, threads=threads,
                                 chunk_size=chunk_size, java_options=java_options)


def _intersect_vcfs(vcf_a: Path, vcf_b: Path, out_vcf: Path, threads: int = 1):
    t = str(max(1, threads))
    for v in [vcf_a, vcf_b]:
        with open(f"{v}.gz", "wb") as f_out:
            subprocess.run(["bgzip", "-@", t, "-c", str(v)], stdout=f_out, check=True)
        subprocess.run(["bcftools", "index", "--threads", t, f"{v}.gz"], check=True)

    subprocess.run([
        "bcftools", "isec", "--threads", t, "-n=2", "-w", "1",
        f"{vcf_a}.gz", f"{vcf_b}.gz",
        "-O", "v", "-o", str(out_vcf)
    ], check=True)


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
    auto_gatk: bool = False,
    threads: int = 1,
    hc_chunk_mb: int = 25,
    gatk_java_options: Optional[str] = None,
    extract_window_mb: int = 5,
) -> tuple[TechnologyResult, TechnologyResult]:

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    variant_positions = load_variant_positions(vcf, threads=threads) if vcf else set()

    bam_path_a = _get_or_map_bam(label_a, out_dir, reads_a, bam_a, reference, preset, min_mapq, threads)
    bam_path_b = _get_or_map_bam(label_b, out_dir, reads_b, bam_b, reference, preset, min_mapq, threads)

    variant_positions = _maybe_call_variants(
        auto_gatk, reference, out_dir, label_a, label_b, bam_path_a, bam_path_b,
        variant_positions, threads, hc_chunk_mb, gatk_java_options)

    result_a = _extract_and_classify(label_a, bam_path_a, out_dir, reference, variant_positions,
                                     min_mapq, min_base_qual, threads, extract_window_mb)
    result_b = _extract_and_classify(label_b, bam_path_b, out_dir, reference, variant_positions,
                                     min_mapq, min_base_qual, threads, extract_window_mb)

    _write_summary(out_dir, label_a, label_b, result_a, result_b)
    return result_a, result_b


def _maybe_call_variants(auto_gatk, reference, out_dir, label_a, label_b, bam_a, bam_b,
                         variant_positions, threads, hc_chunk_mb, gatk_java_options):
    if not auto_gatk:
        return variant_positions
    if not reference:
        raise ValueError("--reference is required for --auto-gatk")
    _prepare_reference(Path(reference))

    vcf_a = out_dir / label_a / f"{label_a}.vcf"
    vcf_b = out_dir / label_b / f"{label_b}.vcf"
    consensus_vcf = out_dir / "consensus.vcf"

    for bam_p, vcf_p in ((bam_a, vcf_a), (bam_b, vcf_b)):
        run_haplotypecaller_parallel(
            bam_p, reference, vcf_p, threads=threads,
            chunk_size=hc_chunk_mb * 1_000_000, java_options=gatk_java_options)
    _intersect_vcfs(vcf_a, vcf_b, consensus_vcf, threads=threads)
    return load_variant_positions(str(consensus_vcf), threads=threads)


def _write_summary(out_dir, label_a, label_b, result_a, result_b):
    summary = summarize_comparison({label_a: result_a.distribution,
                                    label_b: result_b.distribution})
    summary_path = out_dir / "comparison_summary.txt"
    summary_path.write_text(summary + "\n")
    print()
    print(summary)
    print(f"\n[done] full report written under {out_dir}/")


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
    p.add_argument("--preset", default="sr", help="mapping preset for mapping (default: sr)")
    p.add_argument("--min-mapq", type=int, default=1)
    p.add_argument("--min-base-qual", type=int, default=0)
    p.add_argument("--auto-gatk", action="store_true", help="Auto gen VCF with GATK")
    add_threading_args(p)
    return p


def add_threading_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--threads", type=int, default=4,
                   help="threads/processes used by mapping, variant calling and mismatch extraction")
    p.add_argument("--hc-chunk-mb", type=int, default=25,
                   help="default: 25")
    p.add_argument("--gatk-java-options", default=None,
                   help='default: no limit')
    p.add_argument("--extract-window-mb", type=int, default=5,
                   help="default: 5")


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
        auto_gatk=args.auto_gatk,
        threads=args.threads,
        hc_chunk_mb=args.hc_chunk_mb,
        gatk_java_options=args.gatk_java_options,
        extract_window_mb=args.extract_window_mb,
    )


if __name__ == "__main__":
    main()