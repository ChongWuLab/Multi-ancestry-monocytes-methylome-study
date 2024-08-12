# Multi-ancestry-monocytes-methylome-study
## Introduction
This is a comprehensive multi-ancestry methylome-wide association study (MWAS) conducted on purified monocytes from European American (EA) and African American (AA) populations.
Key Features:  
*	Whole-genome bisulfite sequencing (WGBS) data from 298 EA and 160 AA individuals  
*	Analysis of over 25 million methylation sites  
*	Identification of cis- and trans-methylation quantitative trait loci (meQTLs)  
*	Development of population-specific DNA methylation imputation models  
*	MWAS analysis of 41 complex traits using Million Veteran Program (MVP) data  

Our study provides:  
*	meQTL mappings results: Direct links between DNA methylation and genetic variants.  
*	CpG-trait associations: Direct links between specific methylation sites and complex traits  
*	Gene-trait associations: Aggregated effects of methylation on genes associated with various phenotypes  

This resource bridges the gap between genomics and the monocyte methylome, offering insights into:  
*	Genetic regulation of DNA methylation  
*	Novel methylation-phenotype associations  
*	Transferability of findings across diverse ancestries  

This work is partially supported by grants from the NIH (U19AG055373, R01AG061917, R01AR069055, P20GM109036, R01CA263494). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.    

Please cite the following manuscript for using DNAm models built and association results by our work:  

> Zhang, W., Zhang, X., Qiu, C., Zhang, Z., Su, K., Luo, Z., Liu, M., Zhao, B., Wu, L., Tian, Q., Shen, H., Wu, C. and Deng, H., 2024. An atlas of genetic effects on the monocyte methylome across European and African populations. Under Review.

The MWAS results are provided in our website: [gcbhub.org](https://www.gcbhub.org/)  
## Workflow
![Alt text](https://github.com/ChongWuLab/Multi-ancestry-monocytes-methylome-study/blob/main/Fig1.png)

## Tables of Contents

## Genetic data processing
`01-processSNP-a.sh`,`01-processSNP-b.sh`,`01-processSNP-c.sh`,`01-processSNP-d.sh`: We followed [GoDMC](https://github.com/MRCIEU/godmc)[1] pipeline to process genetic data.
## DNAm processing
`02a-normalize_DNAm.R`: DNAm data normalization  
`02b-calculate_nongenetic_pc.R`: Principal component analysis (PCA) was conducted on the DNAm data, focusing on the 20,000 most variable CpG sites to obtain the top ten nongenetic PCs.  
`02c-processCov.R`: We used age, BMI, smoking, alcohol consumption, proportion of blood cells with most variation (B cells, monocytes, neutrophils), and 10 genetic PCs as well as 10 nongenetic PCs to adjust for possible confounding.  
`02d-cpg_annotation.R`: Functions to get gene annotation for CpG sites.  

## Heritability
`03-heritability.R`: Calculate the cis-heritability for each CpG site using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#GREML).

## meQTL analysis
`04-mQTL.R`: Perform meQTL analysis using MatrixEQTL package in R.

## Fine mapping
`04b-fine-mapping.R`: Conduct fine-mapping using [susieR](https://github.com/stephenslab/susieR).

## Methylome-wide association studies
`05-construct_weights.R`: Build imputation models using elastic net regression for CpG sites on genetic data. Adjusting for potential confounders: age, BMI, smoking, alcohol consumption, proportion of blood cells with most variation (B cells, monocytes, neutrophils), and 10 genetic PCs as well as 10 nongenetic PCs. Use cross-validation for model evaluation.  
`06-association_MVP.R`: Perform MWAS on 41 traits in Million Veteran Program (MVP) dataset. We used the pipeline by [Haoran Xue, 2023](https://www.tandfonline.com/doi/full/10.1080/01621459.2023.2183127).[2]
`07-gene_mapping.R`: We mapped CpG sites from all associations to the genes and used Aggregated Cauchy Association Test (ACAT)[3] to detect significant gene-phenotype associations.  

## Post-MWAS analysis
`08a-coloc.R`, `08b-MR.R`: We aligned our MWAS findings with insights from Mendelian Randomization (MR) and Bayesian colocalization analyses for AA and EA ancestries separately.  

## Codes to generate figures in the manuscript
The codes in the `/ana_Res` folder were used to generate figures in the manuscript.  

## Supporting files
`support_MWAS.R`, `BLISSAssociation_Support.R` and `support.R` are the supporting functions to run the pipeline.  
## References
1. Min, J. L., Hemani, G., Hannon, E., Dekkers, K. F., Castillo-Fernandez, J., Luijk, R., ... & Visscher, P. M. (2021). Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation. Nature genetics, 53(9), 1311-1321.
2. Xue, H., Shen, X., & Pan, W. (2023). Causal inference in transcriptome-wide association studies with invalid instruments and GWAS summary data. Journal of the American Statistical Association, 118(543), 1525-1537.
3. Liu, Y., Chen, S., Li, Z., Morrison, A. C., Boerwinkle, E., & Lin, X. (2019). ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. The American Journal of Human Genetics, 104(3), 410-421.
