"""
Correctness tests for the multi-process mismatch extraction.

The ORIGINAL sequential implementation (per-base reference fetch, list-based
classify, DictWriter CSV, statistics-based distribution) is kept below as an
oracle. The new code must reproduce it exactly: same mismatches, same order,
same CSV bytes, same statistics, for any thread count and window size.
"""
import csv
import random
from collections import Counter
from dataclasses import asdict

import pysam
import pytest

from seqerr.bam_reader import (
    CSV_FIELDS, Mismatch, extract_errors_parallel, iter_mismatches,
    iter_mismatches_parallel, plan_windows, _min_mapq_ok,
)
from seqerr.mismatch_classify import classify_mismatches
from seqerr.quality_stats import QualityDistribution


# ----------------------------------------------------------------- oracle ---

def oracle_iter_mismatches(bam_path, technology, reference_fasta=None,
                           min_mapq=1, min_base_qual=0):
    fasta = pysam.FastaFile(str(reference_fasta)) if reference_fasta else None
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if not _min_mapq_ok(read, min_mapq):
                continue
            if read.query_sequence is None or read.query_qualities is None:
                continue
            chrom = read.reference_name
            aligned_pairs = read.get_aligned_pairs(with_seq=(fasta is None),
                                                   matches_only=False)
            for query_pos, ref_pos, *seq in aligned_pairs:
                if query_pos is None or ref_pos is None:
                    continue
                read_base = read.query_sequence[query_pos].upper()
                if fasta is not None:
                    ref_base = fasta.fetch(chrom, ref_pos, ref_pos + 1).upper()
                else:
                    ref_base = seq[0].upper()
                if ref_base == read_base or ref_base not in "ACGT":
                    continue
                base_qual = read.query_qualities[query_pos]
                if base_qual < min_base_qual:
                    continue
                yield Mismatch(
                    read_name=read.query_name, technology=technology, chrom=chrom,
                    ref_pos=ref_pos + 1, ref_base=ref_base, read_base=read_base,
                    base_qual=base_qual, mapping_qual=read.mapping_quality,
                    cigar=read.cigarstring or "", is_read1=read.is_read1,
                    is_reverse=read.is_reverse)
    if fasta is not None:
        fasta.close()


def oracle_write_csv(mismatches, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for m in mismatches:
            w.writerow(asdict(m))


CONTIGS = [("chr1", 6000), ("chr2", 3500), ("chr3", 800), ("chrEmpty", 500)]


def _make_reference(tmp_path, rng):
    seqs = {}
    lines = []
    for name, length in CONTIGS:
        s = [rng.choice("ACGT") for _ in range(length)]
        for _ in range(length // 100):
            s[rng.randrange(length)] = "N"
        seq = "".join(s)
        a = rng.randrange(length - 200)
        seq = seq[:a] + seq[a:a + 150].lower() + seq[a + 150:]
        seqs[name] = seq
        lines.append(f">{name}\n{seq}\n")
    fa = tmp_path / "ref.fa"
    fa.write_text("".join(lines))
    pysam.faidx(str(fa))
    return fa, seqs


def _build_read(ref, start, rng):
    kind = rng.choice(["M", "M", "M", "SM", "MIM", "MDM", "MS", "MNM", "SMS"])
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
        c = rng.randint(1, 8); q.append(junk(c)); cig.append((4, c)); m(rng.randint(60, 90))
    elif kind == "MIM":
        m(rng.randint(20, 40)); c = rng.randint(1, 4); q.append(junk(c)); cig.append((1, c)); m(rng.randint(20, 40))
    elif kind == "MDM":
        m(rng.randint(20, 40)); c = rng.randint(1, 5); take(c); cig.append((2, c)); m(rng.randint(20, 40))
    elif kind == "MS":
        m(rng.randint(60, 90)); c = rng.randint(1, 8); q.append(junk(c)); cig.append((4, c))
    elif kind == "MNM":
        m(rng.randint(20, 40)); c = rng.randint(10, 60); take(c); cig.append((3, c)); m(rng.randint(20, 40))
    else:
        c = rng.randint(1, 6); q.append(junk(c)); cig.append((4, c)); m(rng.randint(50, 80))
        c2 = rng.randint(1, 6); q.append(junk(c2)); cig.append((4, c2))
    seq = list("".join(q))
    aligned_q = set()
    qp = 0
    for op, n in cig:
        if op in (0,):
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


def _make_bam(tmp_path, seqs, rng, n_reads=900):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in CONTIGS]}
    recs = []
    for i in range(n_reads):
        tid = rng.choices([0, 1, 2], weights=[6, 3, 1])[0]
        name, length = CONTIGS[tid]
        start = rng.randrange(0, length - 260)
        seq, cig = _build_read(seqs[name], start, rng)
        a = pysam.AlignedSegment()
        a.query_name = f"r{i}"
        a.query_sequence = seq
        flag = rng.choice([0, 16, 1 | 64, 1 | 128 | 16, 1 | 64 | 16])
        x = rng.random()
        if x < 0.03: flag |= 1024
        elif x < 0.05: flag |= 256
        elif x < 0.07: flag |= 2048
        elif x < 0.09: flag |= 4
        a.flag = flag
        a.reference_id = tid
        a.reference_start = start
        a.mapping_quality = rng.choice([0, 1, 20, 60])
        a.cigartuples = cig
        a.query_qualities = pysam.qualitystring_to_array(
            "".join(chr(33 + rng.randint(2, 41)) for _ in seq))
        recs.append((tid, start, i, a))
    recs.sort(key=lambda t: t[:3])
    bam_path = tmp_path / "t.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for *_, a in recs:
            out.write(a)
    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("par")
    rng = random.Random(7)
    fa, seqs = _make_reference(tmp, rng)
    bam = _make_bam(tmp, seqs, rng)
    return fa, bam, tmp



def test_new_sequential_scan_equals_oracle(data):
    fa, bam, _ = data
    expect = list(oracle_iter_mismatches(bam, "tech", fa))
    assert len(expect) > 500
    assert list(iter_mismatches(bam, "tech", fa)) == expect


@pytest.mark.parametrize("min_mapq,min_bq", [(0, 0), (1, 0), (20, 0), (1, 20), (60, 35)])
def test_filters_equal_oracle(data, min_mapq, min_bq):
    fa, bam, _ = data
    expect = list(oracle_iter_mismatches(bam, "t", fa, min_mapq, min_bq))
    assert list(iter_mismatches(bam, "t", fa, min_mapq=min_mapq, min_base_qual=min_bq)) == expect


@pytest.mark.parametrize("threads,window", [(1, 5_000_000), (3, 97), (4, 500), (2, 1_000_000), (5, 1)])
def test_iter_mismatches_parallel_equals_oracle(data, threads, window):
    fa, bam, _ = data
    expect = list(oracle_iter_mismatches(bam, "t", fa))
    got = list(iter_mismatches_parallel(bam, "t", fa, threads=threads, window_size=window))
    assert got == expect


@pytest.mark.parametrize("threads,window", [(1, 5_000_000), (3, 97), (4, 500), (6, 250)])
@pytest.mark.parametrize("with_variants", [False, True])
def test_extract_errors_parallel_equals_oracle(data, threads, window, with_variants):
    fa, bam, tmp = data
    allm = list(oracle_iter_mismatches(bam, "tech", fa))
    variants = set()
    if with_variants:
        rng = random.Random(3)
        variants = {(m.chrom, m.ref_pos, m.ref_base, m.read_base)
                    for m in rng.sample(allm, len(allm) // 3)}
        errors, vars_ = classify_mismatches(allm, variants)
    else:
        errors, vars_ = allm, []

    out = tmp / f"o_{threads}_{window}_{with_variants}"
    out.mkdir(exist_ok=True)
    res = extract_errors_parallel(bam, "tech", out / "errors.csv", out / "variants.csv",
                                  reference_fasta=fa, variant_positions=variants,
                                  threads=threads, window_size=window)

    oracle_write_csv(errors, out / "o_errors.csv")
    oracle_write_csv(vars_, out / "o_variants.csv")
    assert (out / "errors.csv").read_bytes() == (out / "o_errors.csv").read_bytes()
    assert (out / "variants.csv").read_bytes() == (out / "o_variants.csv").read_bytes()

    assert (res.n_mismatches, res.n_errors, res.n_variants) == (len(allm), len(errors), len(vars_))
    assert res.error_hist == Counter(m.base_qual for m in errors)
    assert not list(out.glob(".*_parts"))


def test_histogram_stats_equal_from_mismatches(data):
    fa, bam, _ = data
    errors = list(oracle_iter_mismatches(bam, "tech", fa))
    a = QualityDistribution.from_mismatches("tech", errors)
    b = QualityDistribution.from_histogram("tech", Counter(m.base_qual for m in errors))
    assert (a.n, a.histogram) == (b.n, b.histogram)
    assert b.mean_qual == pytest.approx(a.mean_qual, abs=1e-9)
    assert b.median_qual == a.median_qual
    assert b.stdev_qual == pytest.approx(a.stdev_qual, abs=1e-9)


@pytest.mark.parametrize("quals", [[], [7], [5, 9], [1, 2, 3], [4, 4, 4, 4], [2, 40, 40, 2, 17], [30] * 99 + [31]])
def test_from_histogram_small_cases(quals):
    ms = [Mismatch("r", "t", "c", 1, "A", "C", q, 60, "1M", True, False) for q in quals]
    a = QualityDistribution.from_mismatches("t", ms)
    b = QualityDistribution.from_histogram("t", Counter(quals))
    assert (a.n, a.histogram) == (b.n, b.histogram)
    assert b.mean_qual == pytest.approx(a.mean_qual)
    assert b.median_qual == a.median_qual
    assert b.stdev_qual == pytest.approx(a.stdev_qual)


def test_plan_windows_tiles_contigs_without_gaps_or_overlap(data):
    _, bam, _ = data
    w = plan_windows(bam, 700)
    assert [x[0] for x in w] == list(range(len(w)))
    assert "chrEmpty" not in {x[1] for x in w}
    for name, length in CONTIGS[:3]:
        spans = [(s, e) for _, c, s, e in w if c == name]
        assert spans[0][0] == 0 and spans[-1][1] == length
        assert all(spans[i][1] == spans[i + 1][0] for i in range(len(spans) - 1))


def test_md_tag_path_without_reference(tmp_path):
    rng = random.Random(11)
    ref = "".join(rng.choice("ACGT") for _ in range(2000))
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "c", "LN": 2000}]}
    path = tmp_path / "md.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(60):
            start = i * 25
            seq = list(ref[start:start + 80])
            md, run = "", 0
            for j, b in enumerate(seq):
                if rng.random() < 0.05:
                    alt = rng.choice([x for x in "ACGT" if x != b])
                    md += f"{run}{b}"; run = 0; seq[j] = alt
                else:
                    run += 1
            md += str(run)
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"; a.query_sequence = "".join(seq); a.flag = 0
            a.reference_id = 0; a.reference_start = start; a.mapping_quality = 60
            a.cigartuples = [(0, 80)]; a.set_tag("MD", md)
            a.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + rng.randint(2, 40)) for _ in seq))
            out.write(a)
    pysam.index(str(path))
    expect = list(oracle_iter_mismatches(path, "t", None))
    assert expect
    assert list(iter_mismatches(path, "t")) == expect


def test_pipeline_compare_threads_match(data, tmp_path):
    from seqerr.pipeline import _qual_histogram_from_csv
    fa, bam, tmp = data
    out = tmp_path / "x"; out.mkdir()
    res = extract_errors_parallel(bam, "t", out / "e.csv", out / "v.csv", reference_fasta=fa, threads=2)
    assert _qual_histogram_from_csv(out / "e.csv") == res.error_hist
