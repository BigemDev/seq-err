"""
Counting and analytics tests for base-quality calibration.

"""
import bisect
import math
import random
from collections import Counter

import pysam
import pytest

from seqerr.bam_reader import _min_mapq_ok, extract_errors_parallel
from seqerr.calibration import (
    Calibration, plot_calibration, summarize_calibration,
)
from seqerr.vcf_reader import variant_sites


# ---------------------------------------------------------------- oracle ---

def oracle_counts(bam_path, fa, seqs, sites, min_mapq=1, min_base_qual=0,
                  use_fasta=True):
    """Slow, obviously correct per-base reference implementation.

    Counts every aligned base of every read passing ``_min_mapq_ok`` using
    ``get_aligned_pairs``. Mirrors the counting definition exactly:
    base in an M/=/X block, base_qual >= min_base_qual, ref base in ACGT,
    read base not N, and position not a variant site.
    """
    total = Counter()
    errors = Counter()
    fasta = pysam.FastaFile(str(fa)) if (use_fasta and fa is not None) else None
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch():
                if not _min_mapq_ok(read, min_mapq):
                    continue
                if read.query_sequence is None or read.query_qualities is None:
                    continue
                chrom = read.reference_name
                seq = read.query_sequence.upper()
                quals = read.query_qualities
                site_list = sites.get(chrom)
                if use_fasta:
                    aligned = read.get_aligned_pairs(matches_only=False)
                    for qpos, rpos in aligned:
                        if qpos is None or rpos is None:
                            continue
                        ref_base = fasta.fetch(chrom, rpos, rpos + 1).upper()
                        _consider(qpos, rpos, seq[qpos], ref_base, quals[qpos],
                                  site_list, min_base_qual, total, errors)
                else:
                    aligned = read.get_aligned_pairs(with_seq=True, matches_only=False)
                    for qpos, rpos, rb in aligned:
                        if qpos is None or rpos is None:
                            continue
                        _consider(qpos, rpos, seq[qpos], rb.upper(), quals[qpos],
                                  site_list, min_base_qual, total, errors)
    finally:
        if fasta is not None:
            fasta.close()
    return total, errors


def _consider(qpos, rpos, read_base, ref_base, bq, site_list, min_base_qual,
              total, errors):
    pos1 = rpos + 1
    if ref_base not in "ACGT" or read_base == "N":
        return
    if bq < min_base_qual:
        return
    if site_list is not None:
        i = bisect.bisect_left(site_list, pos1)
        if i < len(site_list) and site_list[i] == pos1:
            return
    total[bq] += 1
    if read_base != ref_base:
        errors[bq] += 1


CONTIGS = [("chr1", 6000), ("chr2", 3500)]


def _make_reference(tmp_path, rng):
    seqs = {}
    lines = []
    for name, length in CONTIGS:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        seqs[name] = seq
        lines.append(f">{name}\n{seq}\n")
    fa = tmp_path / "ref.fa"
    fa.write_text("".join(lines))
    pysam.faidx(str(fa))
    return fa, seqs


def make_md(ref, start, seq, cig):
    """MD tag consistent with ``seq`` aligned to ``ref`` at ``start`` per ``cig``."""
    parts = []
    run = 0
    qpos = 0
    rpos = start
    for op, n in cig:
        if op == 0:
            for i in range(n):
                rb = ref[rpos + i]
                sb = seq[qpos + i]
                if rb == sb:
                    run += 1
                else:
                    if run:
                        parts.append(str(run))
                        run = 0
                    parts.append(rb.upper())
            qpos += n
            rpos += n
        elif op == 1:
            qpos += n
        elif op == 2:
            if run:
                parts.append(str(run))
                run = 0
            parts.append("^" + ref[rpos:rpos + n].upper())
            rpos += n
        elif op == 4:
            qpos += n
        elif op == 3:
            rpos += n
    if run:
        parts.append(str(run))
    return "".join(parts)


def _make_read(ref, start, rng):
    kind = rng.choice(["M", "M", "SM", "MIM", "MDM", "MS", "MNM"])
    q, cig, rpos = [], [], start

    def take(n):
        nonlocal rpos
        s = ref[rpos:rpos + n]
        rpos += n
        return s

    def junk(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def m(n):
        q.append(take(n).upper())
        cig.append((0, n))

    if kind == "M":
        m(rng.randint(60, 100))
    elif kind == "SM":
        c = rng.randint(1, 8)
        q.append(junk(c))
        cig.append((4, c))
        m(rng.randint(60, 90))
    elif kind == "MIM":
        m(rng.randint(20, 40))
        c = rng.randint(1, 4)
        q.append(junk(c))
        cig.append((1, c))
        m(rng.randint(20, 40))
    elif kind == "MDM":
        m(rng.randint(20, 40))
        c = rng.randint(1, 5)
        take(c)
        cig.append((2, c))
        m(rng.randint(20, 40))
    elif kind == "MS":
        m(rng.randint(60, 90))
        c = rng.randint(1, 8)
        q.append(junk(c))
        cig.append((4, c))
    else:
        m(rng.randint(20, 40))
        c = rng.randint(10, 60)
        take(c)
        cig.append((3, c))
        m(rng.randint(20, 40))
    seq = list("".join(q))
    aligned_q = set()
    qp = 0
    for op, n in cig:
        if op == 0:
            aligned_q.update(range(qp, qp + n))
        if op in (0, 1, 4):
            qp += n
    for i in aligned_q:
        r = rng.random()
        if r < 0.03:
            seq[i] = rng.choice([b for b in "ACGT" if b != seq[i]])
        elif r < 0.033:
            seq[i] = "N"
    return "".join(seq), cig


def _make_bam(tmp_path, seqs, rng, n_reads=400):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in CONTIGS]}
    recs = []
    for i in range(n_reads):
        tid = rng.choices([0, 1], weights=[6, 3])[0]
        name, length = CONTIGS[tid]
        start = rng.randrange(0, length - 260)
        seq, cig = _make_read(seqs[name], start, rng)
        a = pysam.AlignedSegment()
        a.query_name = f"r{i}"
        a.query_sequence = seq
        flag = rng.choice([0, 16, 1 | 64, 1 | 128 | 16, 1 | 64 | 16])
        x = rng.random()
        if x < 0.03:
            flag |= 1024
        elif x < 0.05:
            flag |= 256
        elif x < 0.07:
            flag |= 2048
        elif x < 0.09:
            flag |= 4
        a.flag = flag
        a.reference_id = tid
        a.reference_start = start
        a.mapping_quality = rng.choice([0, 1, 20, 60])
        a.cigartuples = cig
        a.query_qualities = pysam.qualitystring_to_array(
            "".join(chr(33 + rng.randint(2, 41)) for _ in seq))
        a.set_tag("MD", make_md(seqs[name], start, seq, cig))
        recs.append((tid, start, i, a))
    recs.sort(key=lambda t: t[:3])
    bam_path = tmp_path / "t.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for *_, a in recs:
            out.write(a)
    pysam.index(str(bam_path))
    return bam_path


def _make_variants(seqs, rng):
    s = seqs["chr1"]
    positions = set()
    for p in (100, 500, 1000):
        positions.add(("chr1", p + 1, s[p],
                       rng.choice([b for b in "ACGT" if b != s[p]])))
    positions.add(("chr1", 2001, s[2000:2003], "A"))
    return positions


# ------------------------------------------------- counting correctness ---

@pytest.mark.parametrize("threads,window", [(1, 5_000_000), (3, 97), (4, 500),
                                            (2, 1_000_000)])
@pytest.mark.parametrize("use_fasta", [True, False])
def test_counting_matches_oracle(tmp_path, threads, window, use_fasta):
    rng = random.Random(5)
    fa, seqs = _make_reference(tmp_path, rng)
    bam = _make_bam(tmp_path, seqs, rng)
    variant_positions = _make_variants(seqs, rng)
    sites = variant_sites(variant_positions)
    out = tmp_path / "out"
    out.mkdir()
    res = extract_errors_parallel(
        bam, "tech", out / "errors.csv", out / "variants.csv",
        reference_fasta=fa if use_fasta else None,
        variant_positions=variant_positions,
        variant_sites=sites,
        threads=threads, window_size=window,
    )
    ototal, oerrors = oracle_counts(bam, fa, seqs, sites, use_fasta=use_fasta)
    assert res.calib_total_hist == ototal
    assert res.calib_error_hist == oerrors


def test_md_tags_match_reference(tmp_path):
    """The MD tags written for the no-fasta path reproduce the reference."""
    rng = random.Random(5)
    fa, seqs = _make_reference(tmp_path, rng)
    bam = _make_bam(tmp_path, seqs, rng)
    fasta = pysam.FastaFile(str(fa))
    try:
        with pysam.AlignmentFile(str(bam), "rb") as bam:
            for read in bam.fetch():
                for qpos, rpos, rb in read.get_aligned_pairs(
                        with_seq=True, matches_only=False):
                    if qpos is None or rpos is None:
                        continue
                    ref_base = fasta.fetch(read.reference_name, rpos, rpos + 1).upper()
                    assert rb.upper() == ref_base
    finally:
        fasta.close()


def test_min_base_qual_filters_counting(tmp_path):
    rng = random.Random(3)
    fa, seqs = _make_reference(tmp_path, rng)
    bam = _make_bam(tmp_path, seqs, rng)
    sites = variant_sites()
    out = tmp_path / "out"; out.mkdir()
    res = extract_errors_parallel(
        bam, "tech", out / "e.csv", out / "v.csv", reference_fasta=fa,
        variant_positions=set(), variant_sites=sites,
        threads=1, window_size=5_000_000, min_base_qual=30)
    ototal, _ = oracle_counts(bam, fa, seqs, {}, min_base_qual=30,
                              use_fasta=True)
    assert res.calib_total_hist == ototal
    assert all(q >= 30 for q in res.calib_total_hist)


# ----------------------------------------------------------------- sites ---

def test_variant_sites_exclude_positions(tmp_path):
    ref = "".join(random.Random(1).choice("ACGT") for _ in range(300))
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chr1\n{ref}\n")
    pysam.faidx(str(fa))
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 300}]}

    def seg(name, start, seq, cig, mq, bq):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.query_sequence = seq
        a.flag = 0
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = mq
        a.cigartuples = cig
        a.query_qualities = pysam.qualitystring_to_array(bq * len(seq))
        return a

    # SNP site at 1-based 101 (0-based 100): mismatch there.
    # Deletion site at 1-based 201..203 (0-based 200..202): mismatch at 201.
    snp_pos = 101          # 1-based
    del_pos = 201          # 1-based, first base of 3-bp REF
    read1_seq = (ref[90:100] +
                 ("G" if ref[100] != "G" else "C") +
                 ref[101:110])
    read1 = seg("r1", 90, read1_seq, [(0, 20)], 60, "I")
    read2_seq = ref[195:200] + ("T" if ref[200] != "T" else "A") + ref[201:205]
    read2 = seg("r2", 195, read2_seq, [(0, 10)], 60, "I")

    bam_path = tmp_path / "t.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        out.write(read1)
        out.write(read2)
    pysam.index(str(bam_path))

    variant_positions = {("chr1", snp_pos, ref[100], "G"),
                         ("chr1", del_pos, ref[200:203], "A")}
    sites = variant_sites(variant_positions)

    out_no = tmp_path / "no"; out_no.mkdir()
    res_no = extract_errors_parallel(
        bam_path, "tech", out_no / "e.csv", out_no / "v.csv",
        reference_fasta=fa, variant_positions=variant_positions,
        variant_sites={}, threads=1, window_size=10)
    out_yes = tmp_path / "yes"; out_yes.mkdir()
    res_yes = extract_errors_parallel(
        bam_path, "tech", out_yes / "e.csv", out_yes / "v.csv",
        reference_fasta=fa, variant_positions=variant_positions,
        variant_sites=sites, threads=1, window_size=10)

    # errors.csv / variants.csv unchanged by passing variant_sites
    assert (out_no / "e.csv").read_bytes() == (out_yes / "e.csv").read_bytes()
    assert (out_no / "v.csv").read_bytes() == (out_yes / "v.csv").read_bytes()

    # oracle confirms site positions are removed from both counters
    seqs = {"chr1": ref}
    otot_no, oerr_no = oracle_counts(bam_path, fa, seqs, {}, use_fasta=True)
    otot_yes, oerr_yes = oracle_counts(bam_path, fa, seqs, sites, use_fasta=True)
    assert res_no.calib_total_hist == otot_no
    assert res_no.calib_error_hist == oerr_no
    assert res_yes.calib_total_hist == otot_yes
    assert res_yes.calib_error_hist == oerr_yes

    # the two mutated site bases are counted without sites, excluded with sites
    # every base has quality "I" == Phred 40
    assert res_no.calib_total_hist[40] == 30
    assert res_no.calib_error_hist[40] == 2
    assert res_yes.calib_total_hist[40] == 26
    assert res_yes.calib_error_hist[40] == 0


# -------------------------------------------------------------------- math ---

def test_calibration_rows_zero_errors(tmp_path):
    nb = Counter({20: 100, 30: 50, 40: 10})
    ne = Counter({20: 5, 30: 0, 40: 1})
    cal = Calibration("tech", nb, ne)
    rows = cal.rows()
    assert [r[:3] for r in rows] == [(20, 100, 5), (30, 50, 0), (40, 10, 1)]
    assert rows[0][3] == 0.05
    assert rows[0][4] == -10 * math.log10(0.05)
    assert rows[1][3] == 0.0
    assert rows[1][4] is None  # zero errors -> empty empirical_q
    assert rows[2][3] == 0.1
    assert rows[2][4] == -10 * math.log10(0.1)

    path = tmp_path / "cal.csv"
    cal.write_csv(path)
    lines = path.read_text().splitlines()
    assert lines[0] == ("technology,base_qual,n_bases,n_errors,"
                        "error_rate,empirical_q")
    q30 = [l for l in lines if l.startswith("tech,30,")][0]
    assert q30 == "tech,30,50,0,0.0,"


def test_summarize_calibration_table():
    cal = Calibration("t", Counter({20: 100}), Counter({20: 5}))
    s = summarize_calibration({"t": cal})
    assert "t" in s
    assert "observed/expected" in s
    # observed 0.05, expected 100*10^-2 / 100 = 0.01, ratio 5
    assert "5.000000" in s


# ------------------------------------------------- statistical sanity ---

def test_statistical_sanity(tmp_path):
    rng = random.Random(123)
    ref = "".join(rng.choice("ACGT") for _ in range(20000))
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(ref)}]}

    def make_reads(n, read_len, qual_char, p):
        reads = []
        for i in range(n):
            start = rng.randrange(0, len(ref) - read_len)
            bases = list(ref[start:start + read_len])
            for j in range(read_len):
                if rng.random() < p:
                    bases[j] = rng.choice([b for b in "ACGT" if b != bases[j]])
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "".join(bases)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, read_len)]
            a.query_qualities = pysam.qualitystring_to_array(qual_char * read_len)
            reads.append(a)
        return reads

    read_len = 150
    n = 2000  # 300k bases per quality
    q20_char = chr(33 + 20)
    q30_char = chr(33 + 30)
    reads = make_reads(n, read_len, q20_char, 0.01)   # Q20, p=0.01
    reads += make_reads(n, read_len, q30_char, 0.001) # Q30, p=0.001
    reads.sort(key=lambda a: a.reference_start)
    bam_path = tmp_path / "t.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(str(bam_path))

    out = tmp_path / "out"; out.mkdir()
    res = extract_errors_parallel(
        bam_path, "tech", out / "e.csv", out / "v.csv",
        reference_fasta=fa, variant_positions=set(), variant_sites={},
        threads=1, window_size=5_000_000)

    cal = Calibration("tech", res.calib_total_hist, res.calib_error_hist)
    q20 = [r for r in cal.rows() if r[0] == 20][0]
    q30 = [r for r in cal.rows() if r[0] == 30][0]
    # expected errors: p=0.01 -> 3000; p=0.001 -> 300 out of 300k bases
    assert q20[1] == n * read_len     # n_bases at Q20
    assert q30[1] == n * read_len     # n_bases at Q30
    # empirical Q within ~2 binomial standard errors of the declared Q
    assert abs(q20[4] - 20) < 2.0
    assert abs(q30[4] - 30) < 2.0


def test_plot_calibration_writes_png(tmp_path):
    cal = Calibration("tech", Counter({20: 100, 30: 200}), Counter({20: 50, 30: 2}))
    out_png = tmp_path / "cal.png"
    plot_calibration({"tech": cal}, out_png)
    assert out_png.exists()
    assert out_png.stat().st_size > 0
