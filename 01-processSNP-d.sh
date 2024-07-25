#!/bin/bash

set -e
source ./config
module load plink/1.90-beta

for chr in {1..22}; do
    plink \
        --bfile ${rs_tabfile}_all \
        --chr ${chr} \
        --make-bed \
        --out ${rs_tabfile}${chr}
done

# Generate a List of SNPs from the 1000G Dataset:
for i in {1..22}; do
    cut -f2 ${oneKG}${i}.bim >> snps_1000G.txt
done

# Filter the AFA Dataset Based on the 1000G SNPs
for i in {1..22}; do
    plink --bfile ${rs_tabfile}${i} --extract snps_1000G.txt --make-bed --out ${rs_1000G}${i}
done

# Merge all the chrs
plink --bfile ${rs_1000G}1 --bmerge ${rs_1000G}2 --make-bed --out ${rs_1000G}_all
for i in {3..22}; do
plink --bfile ${rs_1000G}_all --bmerge ${rs_1000G}$i --make-bed --out ${rs_1000G}_all
done

rm -f ${rs_1000G}*.*~
# get MAF
plink --bfile ${rs_1000G}_all --freq --out ${res_dir}/rs_1000G_all


