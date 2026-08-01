import gzip
import random

import pysam
import pytest

from seqerr.align import map_reads_to_bam, _sniff_fastx_format
from seqerr.bam_reader import iter_mismatches


def _make_reference(tmp_path, length=2000, seed=1):
    random.seed(seed)
    seq = "".join(random.choice("ACGT") for _ in range(length))
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(fa))
    return fa, seq


def test_map_fastq_reads_to_bam(tmp_path):
    ref_fa, ref_seq = _make_reference(tmp_path)

    # 100bp read from ref[500:600] with one forced mismatch at offset 50
    core = list(ref_seq[500:600])
    mismatch_base = "A" if core[50] != "A" else "C"
    core[50] = mismatch_base
    read_seq = "".join(core)
    qual = "I" * len(read_seq)

    fastq = tmp_path / "reads.fastq"
    fastq.write_text(f"@read1\n{read_seq}\n+\n{qual}\n")

    out_bam = tmp_path / "out.bam"
    n = map_reads_to_bam(fastq, ref_fa, out_bam, preset="sr", technology_tag="illumina")
    assert n == 1
    assert out_bam.exists()
    assert (tmp_path / "out.bam.bai").exists()

    mismatches = list(iter_mismatches(out_bam, "illumina"))  # MD-tag driven, no reference needed
    assert len(mismatches) == 1
    m = mismatches[0]
    assert m.ref_pos == 500 + 50 + 1  # 1-based
    assert m.read_base == mismatch_base


def test_map_fasta_reads_no_qualities(tmp_path):
    ref_fa, ref_seq = _make_reference(tmp_path)
    read_seq = ref_seq[100:180]

    fasta = tmp_path / "reads.fa"
    fasta.write_text(f">read1\n{read_seq}\n")

    out_bam = tmp_path / "out.bam"
    n = map_reads_to_bam(fasta, ref_fa, out_bam, preset="sr", technology_tag="bgi")
    assert n == 1

    mismatches = list(iter_mismatches(out_bam, "bgi"))
    assert mismatches == []  # perfect match, no mismatches expected

    with pysam.AlignmentFile(str(out_bam)) as bam:
        read = next(bam.fetch())
        assert all(q == 30 for q in read.query_qualities)


def test_sniff_detects_fasta_regardless_of_extension(tmp_path):
    f = tmp_path / "reads.weird_ext"
    f.write_text(">r1\nACGTACGT\n")
    assert _sniff_fastx_format(f) == "fasta"


def test_sniff_detects_fastq_gz(tmp_path):
    f = tmp_path / "reads.fastq.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("@r1\nACGTACGT\n+\nIIIIIIII\n")
    assert _sniff_fastx_format(f) == "fastq"


def test_sniff_rejects_garbage(tmp_path):
    f = tmp_path / "notes.txt"
    f.write_text("this is not a sequence file\n")
    with pytest.raises(ValueError, match="does not look like FASTA or FASTQ"):
        _sniff_fastx_format(f)


def test_sniff_rejects_empty_file(tmp_path):
    f = tmp_path / "empty.fasta"
    f.write_text("")
    with pytest.raises(ValueError, match="empty"):
        _sniff_fastx_format(f)


def test_sniff_rejects_missing_file(tmp_path):
    f = tmp_path / "does_not_exist.fastq"
    with pytest.raises(ValueError, match="not found"):
        _sniff_fastx_format(f)


def test_map_reads_to_bam_raises_on_garbage_input(tmp_path):
    ref_fa, _ = _make_reference(tmp_path)
    garbage = tmp_path / "garbage.txt"
    garbage.write_text("not a fastx file\n")
    with pytest.raises(ValueError, match="does not look like FASTA or FASTQ"):
        map_reads_to_bam(garbage, ref_fa, tmp_path / "out.bam")


def test_map_single_read_reports_cigar_and_mismatch(tmp_path):
    from seqerr.align import map_single_read

    ref_fa, ref_seq = _make_reference(tmp_path)
    core = list(ref_seq[500:600])
    core[42] = "A" if core[42] != "A" else "C"
    read_seq = "".join(core)
    qual = "I" * 42 + "#" + "I" * 57

    result = map_single_read(read_seq, ref_fa, read_name="single1", quality=qual)

    assert result.mapped is True
    assert result.chrom == "chr1"
    assert result.ref_start == 500
    assert result.cigar_string == "100M"
    assert result.nm == 1
    assert len(result.mismatches) == 1
    m = result.mismatches[0]
    assert m["ref_pos"] == 500 + 42 + 1
    assert m["base_qual"] == 2


def test_map_single_read_unmapped_for_garbage_sequence(tmp_path):
    from seqerr.align import map_single_read

    ref_fa, _ = _make_reference(tmp_path)
    result = map_single_read("N" * 50, ref_fa, read_name="junk")
    assert result.mapped is False
    assert result.mismatches is None


def test_map_paired_reads_flags_and_tlen(tmp_path):
    from seqerr.align import map_paired_reads_to_bam

    ref_fa, ref_seq = _make_reference(tmp_path, length=3000)

    def rc(s):
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(comp[b] for b in reversed(s))

    r1 = ref_seq[1000:1100]
    r2 = rc(ref_seq[1300:1400])
    fq1 = tmp_path / "r1.fastq"
    fq2 = tmp_path / "r2.fastq"
    fq1.write_text(f"@pair1/1\n{r1}\n+\n{'I'*100}\n")
    fq2.write_text(f"@pair1/2\n{r2}\n+\n{'I'*100}\n")

    out_bam = tmp_path / "paired.bam"
    n = map_paired_reads_to_bam(fq1, fq2, ref_fa, out_bam, technology_tag="illumina", min_mapq=0)
    assert n == 2

    with pysam.AlignmentFile(str(out_bam)) as bam:
        reads = list(bam.fetch())
    assert len(reads) == 2
    r1_seg = next(r for r in reads if r.is_read1)
    r2_seg = next(r for r in reads if r.is_read2)

    assert r1_seg.is_proper_pair and r2_seg.is_proper_pair
    assert r1_seg.next_reference_start == r2_seg.reference_start
    assert r2_seg.next_reference_start == r1_seg.reference_start
    assert r1_seg.template_length == -r2_seg.template_length
    assert r1_seg.mate_is_reverse == r2_seg.is_reverse
