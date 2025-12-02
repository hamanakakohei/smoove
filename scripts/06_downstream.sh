#!/usr/bin/env bash
set -euo pipefail

export PATH=/path/to/bcftools-1.22/bin:$PATH


# 引数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort_vcf) COHORT_VCF="$2"; shift ;;
    --proband) PROBAND="$2"; shift ;;
    --ped) PED="$2"; shift ;;
    --ctrl_af) CTRL_AF="$2"; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


# 0 家系のサンプルリストを作り
mkdir -p results/06/$PROBAND

SAMPLE_LIST=results/06/$PROBAND/samples.txt
samples=$(awk -v PROBAND=$PROBAND '$1==PROBAND {print $2}' $PED)
echo "$samples" | tr ' ' '\n' > $SAMPLE_LIST


# 1 家系ごとにvcfを分けて 
FAM_VCF=results/06/$PROBAND/$PROBAND.fam.vcf.gz
FAM_VCF_COLS=results/06/$PROBAND/svid_vcf_cols.txt

bcftools view \
  --force-samples \
  -S $SAMPLE_LIST \
  $COHORT_VCF \
  | bcftools view \
      --min-ac 1 \
  | bcftools annotate \
      -x INFO/SNAME \
      -Oz \
      -o $FAM_VCF

zgrep -v "^##" $FAM_VCF \
  > $FAM_VCF_COLS
  #| cut -f3,9- \


# 2 vepしてvcfで出力する（--tabだとBNDが無視されるようだ、、、
CACHE_DIR=/betelgeuse07/analysis/ogata/vep115_data/
VEP_VCF=${FAM_VCF/.vcf.gz/.vep.vcf}

docker run --rm \
  -u $(id -u):$(id -g) \
  -v ${CACHE_DIR}:/data \
  -v `pwd`:/work \
  ensemblorg/ensembl-vep:release_115.1 \
  vep \
    --cache \
    --offline \
    --force_overwrite \
    --verbose \
    --fork 4 \
    --assembly GRCh38 \
    --dir_cache /data \
    --per_gene \
    --ccds \
    --mane \
    --symbol \
    --max_sv_size 1000000000 \
    --vcf \
    --input_file /work/$FAM_VCF \
    --output_file /work/$VEP_VCF \
    #--tab \
    #--compress_output gzip \
    #--custom file=/work/results/05/merged.smoove.square.control.site.vcf.gz,short_name=controls1251,format=vcf,type=exact,coords=0,fields=AN%AC \
    #--most_severe \
    #--species homo_sapieens 
    #--everything \
    #--fasta /data/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    #${PLUGIN_OPTS} \


# 3 タブ表にする
VEP_TXT=${VEP_VCF/.vcf/.txt}

cat \
  <(echo -e "ID\t$(bcftools +split-vep -l $VEP_VCF | cut -f2 | tr '\n' '\t' | sed 's/\t$//')") \
  <(bcftools +split-vep $VEP_VCF -d -A tab -f '%ID\t%CSQ\n') \
  > $VEP_TXT


# 3 VCF colやctrl MAFつけつつ、フィルターする
(
  set +u
  source /path/to/python3/env.sh
  scripts/06_downstream.py \
    --vep $VEP_TXT \
    --ac_an $CTRL_AF \
    --vcf_cols $FAM_VCF_COLS \
    --ped $PED \
    --proband $PROBAND \
    --prefix results/06/$PROBAND/$PROBAND.fam.vep.annot
)
