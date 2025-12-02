#!/bin/bash
set -euo pipefail

OUT_DIR=results/06_qc/
mkdir -p $OUT_DIR

ls logs/06/*.log | wc -l > ${OUT_DIR}1.log
ls -l logs/06/   | awk '{print $5,$9}' | sort -n > ${OUT_DIR}2.log
ls logs/06/      | xargs -I{} bash -c 'echo -n {}" "; grep -e RROR -e rror -e found logs/06/{} || true; echo -e ""' > ${OUT_DIR}3.log
ls results/06/   | wc -l > ${OUT_DIR}4.log
ls results/06/   | xargs -I{} bash -c 'echo -n {}" "; ls results/06/{} | wc -l' | sort -k2,2n > ${OUT_DIR}5.log
ls results/06/*/*.fam.vep.annot.txt            | wc -l > ${OUT_DIR}6.log
ls results/06/*/*.fam.vep.annot.AD_XL.txt      | wc -l > ${OUT_DIR}7.log
ls results/06/*/*.fam.vep.annot.AD_XL.rare.txt | wc -l > ${OUT_DIR}8.log
ls results/06/*/*.fam.vep.annot.XR.txt         | wc -l > ${OUT_DIR}9.log
ls results/06/*/*.fam.vep.annot.XR.rare.txt    | wc -l > ${OUT_DIR}10.log
ls results/06/*/*.fam.vep.annot.txt            | xargs -I{} bash -c 'echo {} $(tail -n+2 {} | cut -f1 | sort -u | wc -l)' | sort -k2,2n > ${OUT_DIR}11.log
ls results/06/*/*.fam.vep.annot.AD_XL.txt      | xargs -I{} bash -c 'echo {} $(tail -n+2 {} | cut -f1 | sort -u | wc -l)' | sort -k2,2n > ${OUT_DIR}12.log
ls results/06/*/*.fam.vep.annot.AD_XL.rare.txt | xargs -I{} bash -c 'echo {} $(tail -n+2 {} | cut -f1 | sort -u | wc -l)' | sort -k2,2n > ${OUT_DIR}13.log
ls results/06/*/*.fam.vep.annot.XR.txt         | xargs -I{} bash -c 'echo {} $(tail -n+2 {} | cut -f1 | sort -u | wc -l)' | sort -k2,2n > ${OUT_DIR}14.log
ls results/06/*/*.fam.vep.annot.XR.rare.txt    | xargs -I{} bash -c 'echo {} $(tail -n+2 {} | cut -f1 | sort -u | wc -l)' | sort -k2,2n > ${OUT_DIR}15.log
