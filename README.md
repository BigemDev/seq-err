# seq-err

Compare sequencing error profiles between two technologies (e.g. Illumina vs
BGI). Aligned mismatches are split into **known variants** (present in a VCF)
and **sequencing errors**, then error base qualities are summarized and checked
for calibration (declared Q vs empirically observed error rate).

## What it does

- Map paired reads with `bwa mem`, or use existing sorted/indexed BAMs.
- Extract per-base mismatches from CIGAR + reference (or MD tags).
- Separate variants from errors using a VCF, per-technology VCFs, or GATK
  HaplotypeCaller called in parallel.
- Compare error-quality distributions and produce a calibration plot/CSV.

## Requirements

The Docker image bundles `bwa`, `samtools`, `bcftools`, `tabix`, GATK, pysam and
matplotlib. Building the image downloads GATK and may take a few minutes.

```bash
docker build -t seq-err .
```

## Run the tests

The Dockerfile installs dependencies only, so mount the source when running:

```bash
docker run --rm -v "$PWD":/app seq-err python -m pytest -q
```

Tests that need `bwa`/`gatk`/`samtools`/`bcftools` are skipped automatically
when those tools are not on `PATH`.

## Quick start

Generate a small synthetic paired-end dataset (Illumina-like low error quality
and BGI-like high error quality):

```bash
docker run --rm -v "$PWD":/app seq-err python make_synthetic_data.py
```

Run the full comparison:

```bash
docker run --rm -v "$PWD":/app seq-err \
  python -m seqerr.run_all \
  --reads-a1 synthetic_illumina_R1.fastq --reads-a2 synthetic_illumina_R2.fastq \
  --reads-b1 synthetic_bgi_R1.fastq      --reads-b2 synthetic_bgi_R2.fastq \
  --label-a illumina --label-b bgi \
  --reference synthetic_ref.fa \
  --out-dir synthetic_results
```

Outputs under `synthetic_results/`:

```
comparison_summary.txt           error-quality + calibration summary
calibration.png                  declared vs empirical Q per technology
<label>/errors.csv               mismatches not explained by variants
<label>/variants.csv             mismatches matching known variants
<label>/error_qual_histogram.csv
<label>/calibration.csv
<label>/<label>.sorted.bam*      when mapping from reads
```

## Inputs

Provide either `--bam-a`/`--bam-b` (sorted, indexed) or both read files per
technology. When mapping from reads, `--reference` is required.

Variant sources (choose one):

- `--vcf FILE` — one VCF applied to both technologies.
- `--vcf-a FILE --vcf-b FILE` — per-technology VCFs; writes `consensus.vcf`
  via `bcftools isec` unless `--no-consensus`.
- `--auto-gatk` — call variants per technology with GATK HaplotypeCaller, then
  optionally intersect them.

Other useful flags: `--min-mapq`, `--min-base-qual`, `--preset`
(`sr`, `map-ont`, `map-pb`, `map-hifi`, `intractg`), `--threads`,
`--hc-chunk-mb`, `--extract-window-mb`.

## Individual steps

`seqerr.pipeline` exposes each stage separately:

```bash
python -m seqerr.pipeline map     --reads1 R1.fq --reads2 R2.fq --reference ref.fa \
                                  --out sample.bam --technology illumina
python -m seqerr.pipeline extract --bam sample.bam --technology illumina \
                                  --reference ref.fa --vcf calls.vcf --out-dir out
python -m seqerr.pipeline compare --errors-a out/a/errors.csv --errors-b out/b/errors.csv \
                                  --label-a illumina --label-b bgi --out summary.txt
python -m seqerr.pipeline bqsr    --bam-before raw.bam --bam-after bqsr.bam \
                                  --technology illumina --out deltas.csv
```

## Synthetic variant dataset

`seqerr.synthetic_variants` builds a reproducible reference with four planted
variants (2 SNPs + 2 indels, one homozygous and one heterozygous each) plus a
truth VCF, for validating variant calling:

```bash
python -m seqerr.synthetic_variants
```

Use `--auto-gatk` in `run_all` against this dataset to verify HaplotypeCaller
recovers the planted variants.
