import random

import pysam

from seqerr.run_all import run_all


def _make_reference(tmp_path, length=3000, seed=42):
    random.seed(seed)
    seq = "".join(random.choice("ACGT") for _ in range(length))
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(fa))
    return fa, seq


def _make_fastq(path, ref_seq, n_reads, mismatch_offset, mismatch_qual, base_qual=35):
    lines = []
    for i in range(n_reads):
        start = 200 + i * 100
        core = list(ref_seq[start:start + 100])
        core[mismatch_offset] = "A" if core[mismatch_offset] != "A" else "C"
        seq = "".join(core)
        quals = "".join(
            chr(33 + (mismatch_qual if j == mismatch_offset else base_qual))
            for j in range(len(seq))
        )
        lines.append(f"@read{i}\n{seq}\n+\n{quals}\n")
    path.write_text("".join(lines))


def test_run_all_end_to_end(tmp_path):
    ref_fa, ref_seq = _make_reference(tmp_path)

    illumina_fastq = tmp_path / "illumina.fastq"
    bgi_fastq = tmp_path / "bgi.fastq"
    # Both technologies make the same kind of error at the same read offset,
    # but Illumina reports it honestly with a low quality, BGI with a high one.
    _make_fastq(illumina_fastq, ref_seq, n_reads=5, mismatch_offset=40, mismatch_qual=2)
    _make_fastq(bgi_fastq, ref_seq, n_reads=5, mismatch_offset=40, mismatch_qual=38)

    out_dir = tmp_path / "results"
    result_a, result_b = run_all(
        out_dir=out_dir,
        label_a="illumina",
        label_b="bgi",
        reads_a=str(illumina_fastq),
        reads_b=str(bgi_fastq),
        reference=str(ref_fa),
        min_mapq=0,
    )

    assert result_a.n_errors == 5
    assert result_b.n_errors == 5
    assert result_a.distribution.mean_qual == 2.0
    assert result_b.distribution.mean_qual == 38.0

    # every expected output file should exist on disk
    assert (out_dir / "comparison_summary.txt").exists()
    assert (out_dir / "illumina" / "errors.csv").exists()
    assert (out_dir / "bgi" / "errors.csv").exists()
    assert (out_dir / "illumina" / "illumina.sorted.bam").exists()
    assert (out_dir / "bgi" / "bgi.sorted.bam.bai").exists()

    summary_text = (out_dir / "comparison_summary.txt").read_text()
    assert "illumina" in summary_text and "bgi" in summary_text


def test_run_all_accepts_pre_mapped_bam(tmp_path):
    """If BAMs are already mapped, run_all should skip alignment entirely."""
    ref_fa, ref_seq = _make_reference(tmp_path)
    from seqerr.align import map_reads_to_bam

    fastq = tmp_path / "reads.fastq"
    _make_fastq(fastq, ref_seq, n_reads=3, mismatch_offset=20, mismatch_qual=10)
    bam_a = tmp_path / "a.bam"
    bam_b = tmp_path / "b.bam"
    map_reads_to_bam(fastq, ref_fa, bam_a, technology_tag="illumina")
    map_reads_to_bam(fastq, ref_fa, bam_b, technology_tag="bgi")

    out_dir = tmp_path / "results2"
    result_a, result_b = run_all(
        out_dir=out_dir,
        label_a="illumina",
        label_b="bgi",
        bam_a=str(bam_a),
        bam_b=str(bam_b),
        min_mapq=0,
    )
    assert result_a.n_errors == 3
    assert result_b.n_errors == 3
