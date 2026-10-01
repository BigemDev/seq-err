"""
run_all_paired.py

"""
from __future__ import annotations
 
import argparse
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional
 
import pysam
 
from .quality_stats import summarize_comparison
from .vcf_reader import load_variant_positions
from .run_all import (
    TechnologyResult,
    _extract_and_classify,
    _intersect_vcfs,
    _prepare_reference,
)

_BWA_INDEX_EXT = (".amb", ".ann", ".bwt", ".pac", ".sa")
 
 
def ensure_bwa_index(reference: str | Path) -> None:
    ref = str(reference)
    if all(Path(ref + ext).exists() for ext in _BWA_INDEX_EXT):
        return
    subprocess.run(["bwa", "index", ref], check=True)
 
 
def map_paired_bwa_mem(
    reads_r1: str | Path,
    reads_r2: str | Path,
    reference_fasta: str | Path,
    out_bam: str | Path,
    label: str,
    threads: int = 4,
    sort_threads: Optional[int] = None,
    sort_mem: str = "2G",
    platform: str = "ILLUMINA",
) -> int:
    if threads < 1:
        raise ValueError("threads must be >= 1")
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    ensure_bwa_index(reference_fasta)
 
    if sort_threads is None:
        sort_threads = max(1, min(8, threads // 2))

    rg = f"@RG\\tID:{label}\\tSM:{label}\\tPL:{platform}"
 
    bwa_cmd = ["bwa", "mem", "-t", str(threads), "-R", rg,
               str(reference_fasta), str(reads_r1), str(reads_r2)]
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
) -> tuple[TechnologyResult, TechnologyResult]:
 
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
 
    variant_positions = load_variant_positions(vcf) if vcf else set()
 
    bam_path_a = _get_or_map_bam_paired(
        label_a, out_dir, reads_a1, reads_a2, bam_a, reference, threads)
    bam_path_b = _get_or_map_bam_paired(
        label_b, out_dir, reads_b1, reads_b2, bam_b, reference, threads)
 
    if auto_gatk:
        if not reference:
            raise ValueError("--reference is required for --auto-gatk")
        _prepare_reference(Path(reference))
 
        vcf_a = out_dir / label_a / f"{label_a}.vcf"
        vcf_b = out_dir / label_b / f"{label_b}.vcf"
        consensus_vcf = out_dir / "consensus.vcf"
 
        for bam_p, vcf_p in ((bam_path_a, vcf_a), (bam_path_b, vcf_b)):
            run_haplotypecaller_parallel(
                bam_p, reference, vcf_p, threads=threads,
                chunk_size=hc_chunk_mb * 1_000_000,
                java_options=gatk_java_options,
            )
        _intersect_vcfs(vcf_a, vcf_b, consensus_vcf)
 
        vcf = str(consensus_vcf)
        variant_positions = load_variant_positions(vcf)
 
    result_a = _extract_and_classify(label_a, bam_path_a, out_dir, reference,
                                     variant_positions, min_mapq, min_base_qual)
    result_b = _extract_and_classify(label_b, bam_path_b, out_dir, reference,
                                     variant_positions, min_mapq, min_base_qual)
 
    summary = summarize_comparison({label_a: result_a.distribution,
                                    label_b: result_b.distribution})
    summary_path = out_dir / "comparison_summary.txt"
    summary_path.write_text(summary + "\n")
    print()
    print(summary)
    print(f"\n[done] full report written under {out_dir}/")
 
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
    p.add_argument("--threads", type=int, default=4,
                   help="threads")
    p.add_argument("--hc-chunk-mb", type=int, default=25,
                   help="default: 25")
    p.add_argument("--gatk-java-options", default=None,
                   help='default: no limit')
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
    )
 
 
if __name__ == "__main__":
    main()
