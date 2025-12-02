#!/usr/bin/env python3
import pandas as pd
from cyvcf2 import VCF
from pathlib import Path
import re
import argparse

def get_sv_end_and_type(variant):
    """DEL/DUP/INV/BNDのend位置とタイプを返す"""
    svtype = variant.INFO.get("SVTYPE")
    chrom = variant.CHROM
    start = variant.POS
    end = None

    if svtype in {"DEL", "DUP", "INV"}:
        end = variant.INFO.get("END")
        return chrom, start, chrom, end, svtype
    elif svtype == "BND":
        alt = str(variant.ALT[0])
        # ALT例: N]chr2:12345], or [chr3:7890[N
        clean_alt = alt.replace("N", "").replace("[", "").replace("]", "")
        end_chrom, end = clean_alt.split(":")
        end = int(end)
        return chrom, start, end_chrom, end, "BND"

def determine_regions(chrom, start, end_chrom, end, svtype,
                      small_size, medium_size, extend_ratio, margin):
    """SVタイプとサイズに応じてchr, start, end, nameのリストを返す"""
    regions = []
    if svtype == "BND":
        regions.append((chrom, start - margin, start + margin, "BND_start"))
        regions.append((end_chrom, end - margin, end + margin, "BND_end"))
        return regions

    size = abs(end - start)
    if size <= small_size:
        regions.append((chrom, int(start - size * extend_ratio), int(end + size * extend_ratio), "small"))
    elif size <= medium_size:
        regions.append((chrom, start - margin, start + margin, "start"))
        regions.append((chrom, end - margin, end + margin, "end"))
        regions.append((chrom, int(start - size * extend_ratio), int(end + size * extend_ratio), "mid"))
    else:
        regions.append((chrom, start - margin, start + margin, "start"))
        regions.append((chrom, end - margin, end + margin, "end"))
    return regions

def main():
    parser = argparse.ArgumentParser(description="Prepare IGV snapshot regions for each family.")
    parser.add_argument("--probands", required=True)
    parser.add_argument("--ped", required=True)
    parser.add_argument("--sample-cram", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--call-dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--small_size", type=int, default=5000)
    parser.add_argument("--medium_size", type=int, default=50000)
    parser.add_argument("--extend_ratio", type=float, default=0.25)
    parser.add_argument("--margin", type=int, default=5000)
    args = parser.parse_args()

    # 各ファイル読み込み
    probands = pd.read_csv(args.probands, header=None, names=["proband"])
    ped = pd.read_csv(args.ped, sep="\t", header=None,
                      names=["FID", "IID", "father", "mother", "sex", "pheno"],
                      dtype=str)
    cram_df = pd.read_csv(args.sample_cram, sep=" ", header=None, names=["IID", "cram"])

    records = []
    vcf = VCF(args.vcf)
    vcf_dict = {v.ID: v for v in vcf if v.ID}  # ID→Variant辞書

    for _, row in probands.iterrows():
        proband_id = row["proband"]
        result_path = Path(args.call_dir) / proband_id / f"{proband_id}.fam.vep.annot.AD_XL.rare.txt"

        sv_df = pd.read_csv(result_path, sep="\t", usecols=["ID"], dtype={"ID": str}).drop_duplicates()
        if sv_df.empty:
            continue
        family_id = ped.loc[ped["IID"] == proband_id, "FID"].values[0]
        family_members = ped.loc[ped["FID"] == family_id, "IID"].tolist()
        cram_paths = cram_df.loc[cram_df["IID"].isin(family_members), "cram"].tolist()
        cram_str = ",".join(cram_paths)

        for sv_id in sv_df["ID"]:
            variant = vcf_dict.get(sv_id)
            chrom, start, end_chrom, end, svtype = get_sv_end_and_type(variant)
            regions = determine_regions(chrom, start, end_chrom, end, svtype,
                                        args.small_size, args.medium_size, args.extend_ratio, args.margin)
            for rchrom, rstart, rend, rname in regions:
                records.append([rchrom, rstart, rend, f"{proband_id}_{sv_id}_{rname}", cram_str])

    out_df = pd.DataFrame(records, columns=["chr", "start", "end", "name", "cram_paths"])
    out_df.to_csv(f"{args.prefix}.igv_regions.txt", sep="\t", index=False)

if __name__ == "__main__":
    main()
