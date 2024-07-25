#!/bin/bash

set -e
source ./config
module load plink/1.90-beta
module load R

#cp ${bfile_merge}.bim ${bfileAll}.bim
#cp ${bfile_merge}.bed ${bfileAll}.bed
#cp ${bfile_merge}.fam ${bfileAll}.fam
cp ${bfile_merge}.bim ${bfileAll}.bim
cp ${bfile_merge}.bed ${bfileAll}.bed
cp ${bfile_merge}.fam ${bfileAll}.fam

# Change SNP ids to chr:position:{SNP/INDEL}
echo "Updating SNP ID coding"
cp ${bfileAll}.bim ${bfileAll}.bim.original
awk '{if (($5 == "A" || $5 == "T" || $5 == "C" || $5=="G") &&  ($6 == "A" || $6 == "T" || $6 == "C" || $6=="G")) print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;else print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;}' ${bfileAll}.bim.original > ${bfileAll}.bim


#Recode alleles to uniform format eg. I/D for INDELs
cp ${bfileAll}.bim ${bfileAll}.bim.original2
touch ${SNPfail1}
touch ${bfileAll}.duplicates.txt

Rscript /rsrch5/home/biostatistics/wzhang24/mQTL_project/resources/harmonization.R \
	${bfileAll}.bim \
	${SNPfail1}


# Checking for any duplicate SNPs
cp ${bfileAll}.bim ${bfileAll}.bim.original3
awk '{
	if (++dup[$2] > 1) { 
		print $1, $2".duplicate."dup[$2], $3, $4, $5, $6 
	} else {
		print $0 
	}}' ${bfileAll}.bim.original3 > ${bfileAll}.bim

grep "duplicate" ${bfileAll}.bim | awk '{ print $2 }' > ${bfileAll}.duplicates.txt




# Make GRMs
echo "Creating kinship matrix"
#gunzip -c ${hm3_snps} > temp_hm3snps.txt
plink \
	--bfile ${bfileAll} \
	--maf ${grm_maf_cutoff} \
	--make-grm-bin \
	--out ${grmfile_all} \
  --autosome \
  --threads ${nthreads}

	
#rm temp_hm3snps.txt

# Create pedigree matrix if family data, otherwise remove related individuals from existing kinship and data file
if [ "${related}" = "yes" ]
then
	echo "Creating pedigree GRM"
	Rscript resources/relateds/grm_relateds.R ${grmfile_all} ${grmfile_relateds} ${rel_cutoff}
elif [ "${related}" = "no" ]
then
	echo "Removing any cryptic relateds"
	plink \
		--grm-bin ${grmfile_all} \
		--rel-cutoff ${rel_cutoff} \
		--make-grm-bin \
		--out ${grmfile_unrelateds}

	plink  \
		--bfile ${bfileAll} \
		--keep ${grmfile_unrelateds}.grm.id \
		--make-bed \
		--out ${bfileAll}

else 
	echo "Error: Set related flag in config to yes or no"
	exit 1
fi

# Remove SNPs with low MAF, failing HWE again
plink \
	--bfile ${bfileAll} \
	--maf ${snp_maf} \
	--hwe ${snp_hwe} \
	--make-bed \
	--out ${bfileAll} \
	--chr 1-22 \
	--threads ${nthreads}

# PCA using hm3 list
      
gunzip -c ${hm3_snps_no_ld} > temp_hm3snpsnold.txt
plink \
	--bfile ${bfileAll} \
  --extract temp_hm3snpsnold.txt \
  --indep-pairwise 10000 5 0.1 \
  --maf 0.2 \
	--out ${pca} \
	--autosome \
  --threads ${nthreads}

if [ "${related}" = "no" ]
then
plink \
		--bfile ${bfileAll} \
		--extract ${pca}.prune.in \
		--pca 20 \
		--out ${pca} \
		--threads ${nthreads}
else

	${plink} \
		--bfile ${bfile} \
		--extract ${pca}.prune.in \
		--make-bed \
		--out ${bfile}_ldpruned \
		--threads ${nthreads}

	Rscript resources/genetics/pcs_relateds.R \
		${bfile}_ldpruned \
		${pca} \
		${n_pcs} \
		${nthreads}
fi

# Get genetic outliers
echo "Detecting genetic outliers"

Rscript ${home_dir}/resources/genetic_outliers.R \
	${pcs_all} \
	${pca_sd} \
	${n_pcs} \
	${genetic_outlier_ids} \
	${pcaplot}
 
# If there are any genetic outliers then remove them and recalculate PCs
# Otherwise don't do anything

n_outliers=`wc -l ${genetic_outlier_ids} | awk '{ print $1 }'`
if [ "${n_outliers}" = "0" ]
then
	echo "No genetic outliers detected"
else
	# Remove genetic outliers from data
	echo "Removing ${n_outliers} genetic outliers from data"
	plink \
		--bfile ${bfileAll} \
		--remove ${genetic_outlier_ids} \
		--make-bed \
		--out ${bfileAll} \
		--threads ${nthreads}

	${gcta} \
		--grm ${grmfile_all} \
		--remove ${genetic_outlier_ids} \
		--make-grm-bin \
		--out ${grmfile_all} \
		--thread-num ${nthreads}

fi


# From here on, we have clean data

if [ ! "${n_outliers}" = "0" ]
then
	echo "Recalculating PCs with outliers removed"

	if [ "${related}" = "no" ]
	then
    plink \
	    --bfile ${bfileAll} \
      --extract ${hapmap3} \
      --pca 20 \
	    --out ${pca} \
	    --autosome \
      --threads ${nthreads}

	else

		plink \
			--bfile ${bfileAll} \
			--extract ${pca}.prune.in \
			--make-bed \
			--out ${bfileAll}_ldpruned \
			--autosome \
			--threads ${nthreads}

		Rscript resources/genetics/pcs_relateds.R \
			${bfileAll}_ldpruned \
			${pca} \
			${n_pcs} \
			${nthreads}
	fi

fi

# Get frequencies, missingness, hwe, info scores
plink \
	--bfile ${bfileAll} \
	--freq gz \
	--hardy gz \
	--missing gz \
	--out ${res_dir}/data

# Check missingness
missingness=`zcat ${home_dir}/Results/01-SNP-res/data.imiss | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }'`

echo "Average missingness: ${missingness}"

# Update ids
awk '{print $1,$2}' < ${bfileAll}.fam > ${intersect_ids_plink}
awk '{print $2}' < ${bfileAll}.fam > ${intersect_ids}


rm -f ${bfileAll}.*~
#rm temp_hm3snpsnold.txt
cp ${bfileAll}.bed ${rs_tabfile}_all.bed
cp ${bfileAll}.bim ${rs_tabfile}_all.bim
cp ${bfileAll}.fam ${rs_tabfile}_all.fam

echo "Successfully formatted SNP data"

