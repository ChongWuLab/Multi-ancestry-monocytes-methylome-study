#!/bin/bash

set -e
source ./config
module load plink/1.90-beta
module load R

#SNP filtering
for i in {1..22}; do
  plink \
        --bfile ${bfile_raw}${i} \
        --geno ${snp_miss} \
        --make-bed \
        --out ${bfile}${i} 
done     

for i in {1..22}; do
  plink \
        --bfile ${bfile}${i} \
        --mind ${snp_imiss} \
        --make-bed \
        --out ${bfile}${i} 
done  

for i in {1..22}; do
  plink \
        --bfile ${bfile}${i} \
        --maf ${snp_maf} \
        --hwe ${snp_hwe} \
        --make-bed \
        --out ${bfile}${i}
done     

#merge
plink --bfile ${bfile}1 --bmerge ${bfile}2 --make-bed --out ${bfile_merge}
for i in {3..22}; do
plink --bfile ${bfile_merge} --bmerge ${bfile}$i --make-bed --out ${bfile_merge}
done

echo "Successfully split SNP data"


#extract SNPs from hapmap3
#plink --bfile ${bfile_merge} --extract ${hapmap3} --make-bed --out ${bfile_merge_hm3}