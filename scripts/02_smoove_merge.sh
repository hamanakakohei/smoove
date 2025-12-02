#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift ;;
    --vcf_list) VCF_LIST="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# いろいろ準備
OUT_DIR=results/02
mkdir -p $OUT_DIR

mapfile -t VCF_PATHS < $VCF_LIST


# 実行 Get the union of sites across all samples
docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd):/work \
  brentp/smoove smoove merge \
  --name merged \
  -f $REF \
  --outdir $OUT_DIR \
  ${VCF_PATHS[@]}
  
