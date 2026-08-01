import pysam
import pytest

from seqerr.bam_reader import iter_mismatches
from seqerr.mismatch_classify import classify_mismatches
from seqerr.quality_stats import QualityDistribution


@pytest.fixture
def synthetic_bam(tmp_path):
    ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(f">chr1\n{ref_seq}\n")
    pysam.faidx(str(ref_fa))

    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": len(ref_seq)}]}
    bam_path = tmp_path / "test.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        read_seq = list(ref_seq[0:10])
        mismatch_base = "T" if ref_seq[4] != "T" else "A"
        read_seq[4] = mismatch_base
        a.query_sequence = "".join(read_seq)
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigarstring = "10M"
        a.query_qualities = pysam.qualitystring_to_array("IIII#IIIII")
        bam.write(a)

    sorted_bam = tmp_path / "test.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return sorted_bam, ref_fa, ref_seq, mismatch_base


def test_iter_mismatches_finds_forced_mismatch(synthetic_bam):
    bam_path, ref_fa, ref_seq, mismatch_base = synthetic_bam
    mismatches = list(iter_mismatches(bam_path, "illumina", reference_fasta=ref_fa))

    assert len(mismatches) == 1
    m = mismatches[0]
    assert m.chrom == "chr1"
    assert m.ref_pos == 5  # 1-based
    assert m.ref_base == ref_seq[4]
    assert m.read_base == mismatch_base
    assert m.base_qual == 2  # '#' == Phred 2
    assert m.mapping_qual == 60
    assert m.technology == "illumina"


def test_classify_separates_variant_from_error(synthetic_bam):
    bam_path, ref_fa, ref_seq, mismatch_base = synthetic_bam
    mismatches = list(iter_mismatches(bam_path, "illumina", reference_fasta=ref_fa))

    # Case 1: no known variants -> everything is a sequencing error
    errors, variants = classify_mismatches(mismatches, set())
    assert len(errors) == 1 and len(variants) == 0

    # Case 2: the exact mismatch is a known variant -> reclassified
    known = {("chr1", 5, ref_seq[4], mismatch_base)}
    errors, variants = classify_mismatches(mismatches, known)
    assert len(errors) == 0 and len(variants) == 1


def test_quality_distribution_stats(synthetic_bam):
    bam_path, ref_fa, *_ = synthetic_bam
    mismatches = list(iter_mismatches(bam_path, "illumina", reference_fasta=ref_fa))
    dist = QualityDistribution.from_mismatches("illumina", mismatches)

    assert dist.n == 1
    assert dist.mean_qual == 2.0
    assert dist.histogram[2] == 1


def test_min_mapq_filters_read(synthetic_bam):
    bam_path, ref_fa, *_ = synthetic_bam
    mismatches = list(
        iter_mismatches(bam_path, "illumina", reference_fasta=ref_fa, min_mapq=61)
    )
    assert mismatches == []
