import random
from pathlib import Path

random.seed(42)

ref_seq = "".join(random.choice("ACGT") for _ in range(5000))
Path("synthetic_ref.fa").write_text(f">chr1\n{ref_seq}\n")


def make_fastq(path, n_reads, mismatch_qual, seed):
    random.seed(seed)
    lines = []
    for i in range(n_reads):
        start = random.randint(0, len(ref_seq) - 150)
        core = list(ref_seq[start:start + 150])
        n_errors = random.choice([1, 1, 2])
        err_positions = random.sample(range(150), n_errors)
        for pos in err_positions:
            core[pos] = random.choice([b for b in "ACGT" if b != core[pos]])
        seq = "".join(core)
        quals = ["I"] * 150
        for pos in err_positions:
            quals[pos] = chr(33 + mismatch_qual)
        lines.append(f"@read{i}\n{seq}\n+\n{''.join(quals)}\n")
    Path(path).write_text("".join(lines))

make_fastq("synthetic_illumina.fastq", n_reads=200, mismatch_qual=3, seed=1)
make_fastq("synthetic_bgi.fastq", n_reads=200, mismatch_qual=35, seed=2)

print("Generated: synthetic_ref.fa, synthetic_illumina.fastq, synthetic_bgi.fastq")