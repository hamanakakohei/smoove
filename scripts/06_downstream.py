#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path
from typing import List


def genotype_match(geno: str, patterns: List[str]) -> bool:
    """Return True if genotype string starts with any of the given patterns."""
    return any(geno.startswith(p) for p in patterns)


def ad_xl_filter(row, affected_ids, unaffected_ids):
    affected_all_have = all(
        genotype_match(row[iid], ["0/1:", "1/1:"]) for iid in affected_ids if iid in row
    )
    unaffected_none = all(
        genotype_match(row[iid], ["0/0:", "./."]) for iid in unaffected_ids if iid in row
    )
    return affected_all_have and unaffected_none


def xr_filter(row, family_df):
    for _, member in family_df.iterrows():
        iid = member["IID"]
        sex = member["SEX"]
        pheno = member["PHENOTYPE"]

        if (iid not in row) or (pheno == -9):
            continue # ignore this member

        geno = row[iid]
        # Male
        if sex == 1:
            if pheno == 2 and not genotype_match(geno, ["1/1:"]):
                return False
            if pheno == 1 and not genotype_match(geno, ["0/0:"]):
                return False
        # Female
        elif sex == 2:
            if pheno == 2 and not genotype_match(geno, ["1/1:"]):
                return False
            if pheno == 1 and not genotype_match(geno, ["0/1:", "0/0:", "./."]):
                return False
        # Unknown sex
        else:
            if pheno == 2 and not genotype_match(geno, ["1/1:"]):
                return False
            if pheno == 1 and not genotype_match(geno, ["0/1:", "0/0:", "./."]):
                return False

    return True


def main():
    parser = argparse.ArgumentParser(description="Filter VEP")
    parser.add_argument("--vep", required=True, help="VEP output VCF")
    parser.add_argument("--ac_an", required=True, help="TSV with ID, AC, AN")
    parser.add_argument("--vcf_cols", required=True, help="TSV with ID and other vcf cols including genotypes")
    parser.add_argument("--ped", required=True)
    parser.add_argument("--proband", required=True, type=str)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--ac_threshold_ad_xl", default=1, type=int)
    parser.add_argument("--ac_threshold_xr", default=3, type=int)
    args = parser.parse_args()

    # load vep, MAF, genotype, etc.
    df_vep = pd.read_csv(Path(args.vep), sep='\t', low_memory=False)
    df_ac_an = pd.read_csv(Path(args.ac_an), sep='\t', low_memory=False)
    df_geno = pd.read_csv(Path(args.vcf_cols), sep='\t', low_memory=False)
    df_merged = df_vep.merge(df_ac_an, on="ID", how="left").merge(df_geno, on="ID", how="left")
    df_merged.to_csv(f"{args.prefix}.txt", sep="\t", index=False)

    # load PED
    ped_cols = ["FID", "IID", "PID", "MID", "SEX", "PHENOTYPE"]
    df_ped = pd.read_csv(args.ped, sep="\t", names=ped_cols, dtype={"IID": str})
    family_df = df_ped[df_ped["FID"] == df_ped.loc[df_ped["IID"] == args.proband, "FID"].values[0]]

    # AD/XL
    affected_ids = family_df.query("PHENOTYPE == 2")["IID"].tolist()
    unaffected_ids = family_df.query("PHENOTYPE == 1")["IID"].tolist()

    df_ad_xl = df_merged[df_merged.apply(lambda r: ad_xl_filter(r, affected_ids, unaffected_ids), axis=1)]
    df_ad_xl.to_csv(f"{args.prefix}.AD_XL.txt", sep="\t", index=False)

    df_ad_xl_rare = df_ad_xl[df_ad_xl["AC"] < args.ac_threshold_ad_xl]
    df_ad_xl_rare.to_csv(f"{args.prefix}.AD_XL.rare.txt", sep="\t", index=False)

    # XR
    df_xr = df_merged[
        (df_merged["#CHROM"] == "chrX") &
        (df_merged.apply(lambda r: xr_filter(r, family_df), axis=1))
    ]
    df_xr.to_csv(f"{args.prefix}.XR.txt", sep="\t", index=False)

    df_xr_rare = df_xr[df_xr["AC"] < args.ac_threshold_xr]
    df_xr_rare.to_csv(f"{args.prefix}.XR.rare.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()

