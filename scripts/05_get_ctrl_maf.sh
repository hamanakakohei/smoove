#!/bin/bash
set -euo pipefail


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in_vcf) COHORT_VCF="$2"; shift ;;
    --out_vcf) CTRL_VCF="$2"; shift ;;
    --out_txt) CTRL_TXT="$2"; shift ;;
    --controls) CTRL_LIST="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# vcfをコントロールサンプルに限定して、AC, ANを抜き出してくる
export PATH=/path/to/bcftools-1.22/bin:$PATH

bcftools view \
  -Oz -G \
  --force-samples \
  -S $CTRL_LIST \
  -o $CTRL_VCF \
  $COHORT_VCF

(
  echo -e "ID\tAC\tAN"
  bcftools query \
    -f '%ID\t%INFO/AC\t%INFO/AN\n' \
    $CTRL_VCF \
) > $CTRL_TXT
