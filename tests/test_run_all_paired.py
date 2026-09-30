import random

import pysam
import pytest

from seqerr.run_all_paired import run_all_paired, _get_or_map_bam_paired


def _make_reference(tmp_path, length=3000, seed=42):
    random.seed(seed)
    seq = "".join(random.choice("ACGT") for _ in range(length))
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(fa))
    return fa, seq


def _revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _make_paired_fastq(path1, path2, ref_seq, n_pairs, mismatch_offset,
                       mismatch_qual, base_qual=35, read_len=100, insert=300):
    l1, l2 = [], []
    for i in range(n_pairs):
        start = 200 + i * 400
        frag = ref_seq[start:start + insert]
        r1 = list(frag[:read_len])
        r2 = list(_revcomp(frag)[:read_len])
        r1[mismatch_offset] = "A" if r1[mismatch_offset] != "A" else "C"
        r2[mismatch_offset] = "A" if r2[mismatch_offset] != "A" else "C"
        q = "".join(chr(33 + (mismatch_qual if j == mismatch_offset else base_qual))
                    for j in range(read_len))
        l1.append(f"@pair{i}/1\n{''.join(r1)}\n+\n{q}\n")
        l2.append(f"@pair{i}/2\n{''.join(r2)}\n+\n{q}\n")
    path1.write_text("".join(l1))
    path2.write_text("".join(l2))


def test_run_all_paired_end_to_end(tmp_path):
    ref_fa, ref_seq = _make_reference(tmp_path)

    ill1, ill2 = tmp_path / "ill_1.fq", tmp_path / "ill_2.fq"
    bgi1, bgi2 = tmp_path / "bgi_1.fq", tmp_path / "bgi_2.fq"
    _make_paired_fastq(ill1, ill2, ref_seq, 5, mismatch_offset=40, mismatch_qual=2)
    _make_paired_fastq(bgi1, bgi2, ref_seq, 5, mismatch_offset=40, mismatch_qual=38)

    out_dir = tmp_path / "results"
    a, b = run_all_paired(
        out_dir=out_dir,
        label_a="illumina", label_b="bgi",
        reads_a1=str(ill1), reads_a2=str(ill2),
        reads_b1=str(bgi1), reads_b2=str(bgi2),
        reference=str(ref_fa),
        min_mapq=0,
    )

    assert a.n_errors == 10
    assert b.n_errors == 10
    assert a.distribution.mean_qual == 2.0
    assert b.distribution.mean_qual == 38.0

    assert (out_dir / "comparison_summary.txt").exists()
    assert (out_dir / "illumina" / "errors.csv").exists()

    with pysam.AlignmentFile(str(a.bam_path)) as bam:
        recs = list(bam)
    assert any(r.is_read1 for r in recs) and any(r.is_read2 for r in recs)
    assert all(r.is_paired for r in recs)

