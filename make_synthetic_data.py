#!/usr/bin/env python3
import random
import shutil
import pysam
from pathlib import Path
from seqerr.run_all import run_all

def make_fastq(ref_seq, path, n_reads, mismatch_qual, seed):
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

def main():
    base_dir = Path("synthetic_Tests")
    if base_dir.exists():
        shutil.rmtree(base_dir)
    
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    data_dir.mkdir(parents=True)
    results_dir.mkdir(parents=True)
    random.seed(67)
    ref_seq = "".join(random.choice("ACGT") for _ in range(5000))
    ref_path = data_dir / "synthetic_ref.fa"
    ref_path.write_text(f">chr1\n{ref_seq}\n")

    pysam.faidx(str(ref_path))
    fastq_illumina = data_dir / "synthetic_illumina.fastq"
    fastq_bgi = data_dir / "synthetic_bgi.fastq"

    make_fastq(ref_seq, fastq_illumina, n_reads=200, mismatch_qual=3, seed=1)
    make_fastq(ref_seq, fastq_bgi, n_reads=200, mismatch_qual=35, seed=2)

    result_a, result_b = run_all(
        out_dir=results_dir,
        label_a="illumina",
        label_b="bgi",
        reads_a=str(fastq_illumina),
        reads_b=str(fastq_bgi),
        reference=str(ref_path),
        min_mapq=1,
        auto_gatk=True, 
    )
    print("\n" + "="*40)

    print(f"Raw data: {data_dir}")
    print(f"Results: {results_dir}\n")

if __name__ == "__main__":
    main()