#!/usr/bin/env python3
import random
import pysam
from pathlib import Path


def _make_read(sequence, rng, mismatch_qual):
    bases = list(sequence)
    qualities = ["I"] * len(bases)
    for pos in rng.sample(range(len(bases)), rng.choice([1, 1, 2])):
        bases[pos] = rng.choice([base for base in "ACGT" if base != bases[pos]])
        qualities[pos] = chr(33 + mismatch_qual)
    return "".join(bases), "".join(qualities)


def make_paired_fastq(ref_seq, path1, path2, n_pairs, mismatch_qual, seed):
    rng = random.Random(seed)
    r1_lines, r2_lines = [], []
    complement = str.maketrans("ACGT", "TGCA")
    for i in range(n_pairs):
        start = rng.randint(0, len(ref_seq) - 300)
        fragment = ref_seq[start:start + 300]
        r1, q1 = _make_read(fragment[:150], rng, mismatch_qual)
        r2, q2 = _make_read(fragment[150:][::-1].translate(complement), rng, mismatch_qual)
        r1_lines.append(f"@read{i}/1\n{r1}\n+\n{q1}\n")
        r2_lines.append(f"@read{i}/2\n{r2}\n+\n{q2}\n")
    Path(path1).write_text("".join(r1_lines))
    Path(path2).write_text("".join(r2_lines))

def main():
    data_dir = Path(".")
    random.seed(67)
    ref_seq = "".join(random.choice("ACGT") for _ in range(5000))
    ref_path = data_dir / "synthetic_ref.fa"
    ref_path.write_text(f">chr1\n{ref_seq}\n")

    pysam.faidx(str(ref_path))
    fastq_illumina1 = data_dir / "synthetic_illumina_R1.fastq"
    fastq_illumina2 = data_dir / "synthetic_illumina_R2.fastq"
    fastq_bgi1 = data_dir / "synthetic_bgi_R1.fastq"
    fastq_bgi2 = data_dir / "synthetic_bgi_R2.fastq"

    make_paired_fastq(ref_seq, fastq_illumina1, fastq_illumina2,
                      n_pairs=200, mismatch_qual=3, seed=1)
    make_paired_fastq(ref_seq, fastq_bgi1, fastq_bgi2,
                      n_pairs=200, mismatch_qual=35, seed=2)
    print(f"Generated paired reads and reference under {data_dir.resolve()}")

if __name__ == "__main__":
    main()
