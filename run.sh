#!/bin/bash
set -euo pipefail

export PATH=/path/to/samtools-1.17/bin/:$PATH


# dockerコマンドのために相対パスで指定する
SAMPLE_CRAM_LIST=inputs/sample_cram.txt
REF=ref/Homo_sapiens_assembly38.fasta
EXCLUDE_BED=ref/exclude.cnvnator_100bp.GRCh38.20170403.bed

N_SAMPLE=$(wc -l < "$SAMPLE_CRAM_LIST")


# 01 smoove call
qsub -t 1-"$N_SAMPLE":1 -tc 10 qsubs/01_smoove_call.qsub \
  "$SAMPLE_CRAM_LIST" \
  "$REF" \
  "$EXCLUDE_BED"

bash scripts/01_qc.sh
find results/01/ -name "*.vcf.gz" > results/01/call_vcf.list


# 02 smoove merge
bash scripts/02_smoove_merge.sh \
  --ref $REF \
  --vcf_list results/01/call_vcf.list \
  > logs/02/log 2>&1


# 03 smoove genotype
qsub -t 1-"$N_SAMPLE":1 -tc 10 qsubs/03_smoove_genotype.qsub \
  "$SAMPLE_CRAM_LIST" \
  "$REF" \
  results/02/merged.sites.vcf.gz

bash scripts/03_qc.sh
find results/03/ -name "*.vcf.gz" > results/03/genotype_vcf.list


# 04 smoove paste (SNAMEが明らかに一部のサンプルしか表示されていない、バグ？)
bash scripts/04_smoove_paste.sh \
  --vcf_list results/03/genotype_vcf.list \
  > logs/04/log 2>&1

awk -v OFS="\t" \
  '{gsub(/;/, "and", $3); print}' \
  <(zcat results/04/merged.smoove.square.vcf.gz) \
  | bgzip -c \
  > results/04/merged.smoove.square.reID.vcf.gz


# 05 get ctrl MAF
CTRLS=inputs/controls_2175.txt
scripts/05_get_ctrl_maf.sh \
  --in_vcf results/04/merged.smoove.square.reID.vcf.gz \
  --controls $CTRLS \
  --out_vcf results/05/merged.smoove.square.reID.ctrl.site.vcf.gz \
  --out_txt results/05/svid_af.txt \
  > logs/05/log 2>&1


# 06 vepして、VCF colやctrl MAFつけて、フィルターして
PED=inputs/sr-wgs.ped
PROBAND_LIST=inputs/probands773.txt
N_PROBAND=$(wc -l < $PROBAND_LIST)

qsub -t 1-"$N_PROBAND":1 -tc 2 qsubs/06_downstream.qsub \
  results/04/merged.smoove.square.reID.vcf.gz \
  "$PROBAND_LIST" \
  "$PED" \
  results/05/svid_af.txt 

bash scripts/06_qc.sh


# 07 IGV撮影用インプットファイルを作る
source /path/to/python3/env.sh
scripts/07_make_input_for_igv.py \
  --probands $PROBAND_LIST \
  --ped $PED \
  --sample-cram $SAMPLE_CRAM_LIST \
  --vcf results/04/merged.smoove.square.reID.vcf.gz \
  --call-dir results/06/ \
  --prefix results/07/probands773


# 08 IGVで自動撮影する（cf. github.com/hamanakakohei/igv）
