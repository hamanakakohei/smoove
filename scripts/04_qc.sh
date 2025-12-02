#!/bin/bash
set -euo pipefail

zgrep  -c "^#" results/04/merged.smoove.square.vcf.gz 
zgrep -vc "^#" results/04/merged.smoove.square.vcf.gz 
zgrep  -v "^#" results/04/merged.smoove.square.vcf.gz | awk -F"\t" '{print NF}' | sort -u
