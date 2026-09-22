"""
synthetic_variants.py
 
Generates a fully reproducible synthetic diploid dataset with a KNOWN set of
planted germline variants (the "truth set"), for validating that variant
calling.

A random reference is generated, then two haplotypes are derived from it:
 
  * haplotype A ("ref-like"): carries the homozygous truth variants only
  * haplotype B ("alt-like"): carries the homozygous *and* heterozygous
    truth variants
 
Reads are simulated by repeatedly picking a genome window and a haplotype
(50/50 A or B) to draw the read from. Because homozygous variants sit on
both haplotypes, ~100% of reads covering them show the alt allele; because
heterozygous variants sit only on haplotype B, ~50% of reads covering them
show the alt allele and ~50% show the reference allele -- i.e. realistic
zygosity. A separate, independent low-rate/low-quality per-base error model
is layered on top of every read to mimic sequencer noise; this noise is
*not* a real variant anywhere, so a correct caller should not call it.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

import pysam

@dataclass(frozen=True)
class TruthVariant:
    pos: int
    ref: str
    alt: str
    zygosity: str # "hom" or "het"

    @property
    def vcf_pos(self) -> int:
        return self.pos + 1

    @property
    def genotype(self) -> str:
        return "1/1" if self.zygosity == "hom" else "0/1"




CHROM = "chr1"
REF_LENGTH = 6000
READ_LENGTH = 150
TARGET_DEPTH = 60

TRUTH_VARIANTS: tuple[TruthVariant, ...] = (
    TruthVariant(pos=1000, ref="G", alt="T", zygosity="hom"), # homo SNP
    TruthVariant(pos=2000, ref="G", alt="A", zygosity="het"), # hete SNP
    TruthVariant(pos=3000, ref="GTC", alt="G", zygosity="hom"), # homo 2bp deletion
    TruthVariant(pos=4000, ref="G", alt="GAC", zygosity="het"), # hete 2bp deletion
)

def build_reference(length: int = REF_LENGTH, seed: int = 67) -> str:
    rng = random.Random(seed)
    bases = [rng.choice("ACGT") for _ in range(length)]
    for v in TRUTH_VARIANTS:
        bases[v.pos:v.pos + len(v.ref)] = list(v.ref)
    return "".join(bases)


def _build_haplotype(ref_seq: str, variants: list[TruthVariant]) -> str:
    variants = sorted(variants, key=lambda v: v.pos)
    out = []
    cursor = 0
    for v in variants:
        assert ref_seq[v.pos:v.pos + len(v.ref)] == v.ref, (
            f"reference base mismatch {v.pos}:"
            f"expected {v.ref!r}, instead {ref_seq[v.pos:v.pos + len(v.ref)]!r}"
        )
        out.append(ref_seq[cursor:v.pos])
        out.append(v.alt)
        cursor = v.pos + len(v.ref)
    out.append(ref_seq[cursor:])
    return "".join(out)


def _haplotype_offset(pos: int, variants: list[TruthVariant]) -> int:
    offset = 0
    for v in sorted(variants, key=lambda v: v.pos):
        if v.pos >= pos:
            break
        offset += len(v.alt) - len(v.ref)
    return offset


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def generate_reads_fastq(
    ref_seq: str,
    out_fastq: str | Path,
    variants: tuple[TruthVariant, ...] = TRUTH_VARIANTS,
    read_length: int = READ_LENGTH,
    depth: int = TARGET_DEPTH,
    error_rate: float = 0.003,
    error_qual: int = 8,
    base_qual: int = 40,
    seed: int = 1,
) -> int:
    hom_variants = [v for v in variants if v.zygosity == "hom"]
    het_variants = [v for v in variants if v.zygosity == "het"]

    hap_a = _build_haplotype(ref_seq, hom_variants)
    hap_b = _build_haplotype(ref_seq, hom_variants + het_variants)
    hap_a_variants = hom_variants
    hap_b_variants = hom_variants + het_variants

    rng = random.Random(seed)
    margin = read_length
    n_reads = (depth * (len(ref_seq) - 2 * margin)) // read_length
    lines = []
    n_written = 0
    for i in range(n_reads):
        g_start = rng.randint(margin, len(ref_seq) - margin - read_length)
        if rng.random() < 0.5:
            hap_seq, hap_variants = hap_a, hap_a_variants
        else:
            hap_seq, hap_variants = hap_b, hap_b_variants
 
        local_start = g_start + _haplotype_offset(g_start, hap_variants)
        read = list(hap_seq[local_start:local_start + read_length])
        if len(read) != read_length:
            continue

        quals = [base_qual] * read_length
        for pos in range(read_length):
            if rng.random() < error_rate:
                read[pos] = rng.choice([b for b in "ACGT" if b != read[pos]])
                quals[pos] = error_qual
 
        seq = "".join(read)
        if rng.random() < 0.5:
            seq = _revcomp(seq)
            quals = quals[::-1]
 
        qual_str = "".join(chr(33 + q) for q in quals)
        lines.append(f"@read{i}\n{seq}\n+\n{qual_str}\n")
        n_written += 1
 
    Path(out_fastq).write_text("".join(lines))
    return n_written


def write_truth_vcf(
    variants: tuple[TruthVariant, ...],
    ref_length: int,
    out_vcf: str | Path,
    chrom: str = CHROM,
) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom},length={ref_length}>",
        '##INFO=<ID=SYNTH,Number=0,Type=Flag,Description="Planted synthetic truth variant">',
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    for v in sorted(variants, key=lambda v: v.pos):
        lines.append(
            f"{chrom}\t{v.vcf_pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\tSYNTH\tGT\t{v.genotype}"
        )
    Path(out_vcf).write_text("\n".join(lines) + "\n")


@dataclass
class SyntheticVariantDataset:
    reference_fasta: Path
    reads_fastq: Path
    truth_vcf: Path
    variants: tuple[TruthVariant, ...]
    n_reads: int



def generate_synthetic_variant_dataset(
    out_dir: str | Path,
    seed: int = 67,
    depth: int = TARGET_DEPTH,
) -> SyntheticVariantDataset:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
 
    ref_seq = build_reference(REF_LENGTH, seed=seed)
    reference_fasta = out_dir / "synthetic_variant_ref.fa"
    reference_fasta.write_text(f">{CHROM}\n{ref_seq}\n")
    pysam.faidx(str(reference_fasta))
 
    reads_fastq = out_dir / "synthetic_variant_reads.fastq"
    n_reads = generate_reads_fastq(ref_seq, reads_fastq, depth=depth, seed=seed + 1000)
 
    truth_vcf = out_dir / "synthetic_variant_truth.vcf"
    write_truth_vcf(TRUTH_VARIANTS, REF_LENGTH, truth_vcf)
 
    return SyntheticVariantDataset(
        reference_fasta=reference_fasta,
        reads_fastq=reads_fastq,
        truth_vcf=truth_vcf,
        variants=TRUTH_VARIANTS,
        n_reads=n_reads,
    )
 
 
if __name__ == "__main__":
    ds = generate_synthetic_variant_dataset("synthetic_variant_data")
    print(f"Generated {ds.n_reads} reads {ds.reads_fastq}")
    print(f"Reference {ds.reference_fasta}")
    print(f"Truth VCF ({len(ds.variants)} variants) {ds.truth_vcf}")
