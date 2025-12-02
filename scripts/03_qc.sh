#!/bin/bash
set -euo pipefail

OUT_DIR=results/03_qc/
mkdir -p $OUT_DIR

ls logs/03/*.log | wc -l > ${OUT_DIR}1.log
ls -l logs/03/   | awk '{print $5,$9}' | sort -n > ${OUT_DIR}2.log
ls logs/03/      | xargs -I{} bash -c 'echo -n {}" "; grep -e RROR -e rror -e found logs/03/{} || true; echo -e ""' > ${OUT_DIR}3.log
ls results/03/   | wc -l > ${OUT_DIR}4.log
ls results/03/   | xargs -I{} bash -c 'echo -n {}" "; ls results/03/{} | wc -l' | sort -k2,2n > ${OUT_DIR}5.log
ls results/03/*/*-smoove.genotyped.vcf.gz     | wc -l > ${OUT_DIR}6.log
ls results/03/*/*-smoove.genotyped.vcf.gz.csi | wc -l > ${OUT_DIR}7.log
ls -l results/03/*/*-smoove.genotyped.vcf.gz | awk '{print $5,$9}' | sort -n > ${OUT_DIR}8.log
ls results/03/*/*-smoove.genotyped.vcf.gz | xargs -P 40 -I{} bash -c 'echo {} $(zgrep  -c "^#" {})' | sort -k2,2n > ${OUT_DIR}9.log
ls results/03/*/*-smoove.genotyped.vcf.gz | xargs -P 40 -I{} bash -c 'echo {} $(zgrep -vc "^#" {})' | sort -k2,2n > ${OUT_DIR}10.log
ls results/03/*/*-smoove.genotyped.vcf.gz | xargs -P 40 -I{} bash -c 'echo {} $(zgrep -vc "^#" {} | cut -f3 | md5sum)' | sort -k2,2 > ${OUT_DIR}11.log
ls results/03/*/*-smoove.genotyped.vcf.gz | xargs -P 40 -I{} bash -c 'echo {} $(zcat {} | grep PRPOS)' | sort > ${OUT_DIR}12.log
ls results/03/*/*-smoove.genotyped.vcf.gz | xargs -P 40 -I{} bash -c 'echo {} $(zcat {} | grep PREND)' | sort > ${OUT_DIR}13.log
