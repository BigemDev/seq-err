"""
test_haplotypecaller.py

Uses syntetic_variants to build a reproducible synthetic genome with a KNOWN set
of planted variants (2 SNP + 2 indels, one homozygous and one heterozygous of each).
"""

from __future__ import annotations
 
import shutil
 
import pysam
import pytest
 
from seqerr.align import map_reads_to_bam
from seqerr.run_all import run_all, _prepare_reference, _run_variant_calling
from seqerr.synthetic_variants import (
    CHROM,
    TRUTH_VARIANTS,
    build_reference,
    generate_reads_fastq,
    generate_synthetic_variant_dataset,
)
 
pytestmark = pytest.mark.skipif(
    not all(shutil.which(tool) for tool in ("bwa", "gatk", "samtools", "bcftools")),
    reason="see Dockerfile",
)
 
 
def _read_called_variants(vcf_path) -> set[tuple[str, int, str, str, str]]:
    calls = set()
    with pysam.VariantFile(str(vcf_path)) as vf:
        sample = next(iter(vf.header.samples))
        for rec in vf:
            gt = rec.samples[sample]["GT"]
            gt_str = "/".join(str(a) for a in gt)
            calls.add((rec.chrom, rec.pos, rec.ref, rec.alts[0], gt_str))
    return calls
 
 
def _expected_truth_calls(variants=TRUTH_VARIANTS) -> set[tuple[str, int, str, str, str]]:
    return {(CHROM, v.vcf_pos, v.ref, v.alt, v.genotype) for v in variants}
 
 
def test_haplotypecaller_recovers_known_variants(tmp_path):
    ds = generate_synthetic_variant_dataset(tmp_path)
 
    bam_path = tmp_path / "mapped.bam"
    n_mapped = map_reads_to_bam(
        ds.reads_fastq, ds.reference_fasta, bam_path, technology_tag="synth",
    )
    assert n_mapped == ds.n_reads
 
    _prepare_reference(ds.reference_fasta)
    called_vcf = tmp_path / "called.vcf"
    _run_variant_calling(bam_path, ds.reference_fasta, called_vcf)
 
    actual = _read_called_variants(called_vcf)
    expected = _expected_truth_calls(ds.variants)
 
    assert actual == expected
 
 
def test_run_all_auto_gatk_recovers_known_variants(tmp_path):
    ref_seq = build_reference(seed=67)
    reference_fasta = tmp_path / "ref.fa"
    reference_fasta.write_text(f">{CHROM}\n{ref_seq}\n")
    pysam.faidx(str(reference_fasta))

    fastq_a = tmp_path / "illumina.fastq"
    fastq_b = tmp_path / "bgi.fastq"
    generate_reads_fastq(ref_seq, fastq_a, seed=101)
    generate_reads_fastq(ref_seq, fastq_b, seed=202)
 
    out_dir = tmp_path / "results"
    result_a, result_b = run_all(
        out_dir=out_dir,
        label_a="illumina",
        label_b="bgi",
        reads_a=str(fastq_a),
        reads_b=str(fastq_b),
        reference=str(reference_fasta),
        min_mapq=1,
        auto_gatk=True,
    )
 
    consensus_vcf = out_dir / "consensus.vcf"
    assert consensus_vcf.exists()
    assert _read_called_variants(consensus_vcf) == _expected_truth_calls()
 
    assert result_a.n_variants > 0
    assert result_b.n_variants > 0
