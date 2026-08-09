
# Seqerr



## Deployment

Run using synthetic data

```bash
docker run --rm -v "$(pwd)":/app -w /app python:3.11-slim bash -c "
apt-get update -qq && apt-get install -y -qq build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev >/dev/null &&
pip install -r requirements.txt --break-system-packages -q &&
python make_synthetic_data.py &&
python -m seqerr.run_all \
    --reads-a synthetic_illumina.fastq --label-a illumina \
    --reads-b synthetic_bgi.fastq      --label-b bgi \
    --reference synthetic_ref.fa \
    --out-dir synthetic_results \
    --min-mapq 0
"
```

