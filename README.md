
# Seqerr



## Deployment

Run using synthetic data

```bash
docker run --rm -v "$(pwd)":/app -w /app python:3.11-slim bash -c "
apt-get update -qq && apt-get install -y -qq build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev bwa samtools >/dev/null &&
pip install -r requirements.txt --break-system-packages -q &&
python make_synthetic_data.py &&
python -m seqerr.run_all \
    --reads-a1 synthetic_illumina_R1.fastq --reads-a2 synthetic_illumina_R2.fastq --label-a illumina \
    --reads-b1 synthetic_bgi_R1.fastq      --reads-b2 synthetic_bgi_R2.fastq      --label-b bgi \
    --reference synthetic_ref.fa \
    --out-dir synthetic_results \
    --min-mapq 0
"
```
