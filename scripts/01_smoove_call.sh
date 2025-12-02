#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift ;;
    --exclude) EXCLUDE="$2"; shift ;;
    --cram) CRAM="$2"; shift ;;
    --sample) SAMPLE="$2"; shift ;;
    --threads) THREADS="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
REF_DIR=$(dirname "$(readlink -f "$REF")")
EXCLUDE_DIR=$(dirname "$(readlink -f "$EXCLUDE")")
CRAM_DIR=$(dirname "$(readlink -f "$CRAM")")
OUT_DIR="results/01/${SAMPLE}"

mkdir -p "$OUT_DIR"


# 実行
echo "[`date`] Starting smoove call for $SAMPLE"

docker run --rm \
  --user $(id -u):$(id -g) \
  -v "$REF_DIR":/ref \
  -v "$EXCLUDE_DIR":/exclude \
  -v "$CRAM_DIR":/cram \
  -v "$(pwd)/results":/results \
  brentp/smoove smoove call \
    --name "$SAMPLE" \
    --fasta /ref/$(basename "$REF") \
    --exclude /exclude/$(basename "$EXCLUDE") \
    --genotype /cram/$(basename "$CRAM") \
    --outdir /results/01/"$SAMPLE"/ \
    -p "$THREADS"

echo "[`date`] Finished smoove call for $SAMPLE"


## 中間ファイル削除
: "${OUT_DIR:?OUT_DIR is not set or empty}"
rm "$OUT_DIR/$SAMPLE".{split.bam,split.bam.bai,split.bam.orig.bam,disc.bam,disc.bam.bai,disc.bam.orig.bam} #,histo}
#rm "$OUT_DIR/$SAMPLE"-lumpy-cmd.sh
