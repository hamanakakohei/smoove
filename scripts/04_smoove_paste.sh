#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf_list) VCF_LIST="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
OUT_DIR="results/04"
mkdir -p "$OUT_DIR"

mapfile -t VCF_PATHS < $VCF_LIST


# 実行
echo "[`date`] Starting smoove paste"

docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd):/work \
  brentp/smoove smoove paste \
  --name merged \
  --outdir $OUT_DIR \
  ${VCF_PATHS[@]}
  
echo "[`date`] Finished smoove paste"
