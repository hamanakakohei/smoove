#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift ;;
    --cram) CRAM="$2"; shift ;;
    --site_vcf) SITE_VCF="$2"; shift ;;
    --sample) SAMPLE="$2"; shift ;;
    --threads) THREADS="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
CRAM_DIR=$(dirname "$(readlink -f "$CRAM")")
CRAM_BASE=$(basename "$CRAM")
REF_DIR=$(dirname "$(readlink -f "$REF")")
REF_BASE=$(basename "$REF")

OUT_DIR="results/03/${SAMPLE}"

mkdir -p "$OUT_DIR"


# 実行 genotype each sample at all merged sites
echo "[`date`] Starting smoove genotype for $SAMPLE"

docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd):/work \
  -v "$REF_DIR":/ref \
  -v $CRAM_DIR:/cram \
  brentp/smoove smoove genotype \
    -d -x \
    -p $THREADS \
    --name $SAMPLE \
    --vcf $SITE_VCF \
    --fasta /ref/$REF_BASE \
    --outdir $OUT_DIR \
    /cram/$CRAM_BASE

echo "[`date`] Finished smoove genotype for $SAMPLE"
