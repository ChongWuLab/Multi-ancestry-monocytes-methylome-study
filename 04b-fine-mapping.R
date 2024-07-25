library(susieR)
library(data.table)
library(BEDMatrix)
library(dplyr)
library(stringr)
setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
source("support_MWAS.R")
mqtl_dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/03-mQTL-V2/"
data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"
dat.wgbs.dir = "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"
save.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/07-finemapping/"

args = commandArgs(TRUE)
indx1 = as.numeric(args[[1]])

mqtl.file = list.files(mqtl_dir, pattern = paste0("_", indx1, ".RDS"), full.names = TRUE)
mqtl.res = readRDS(mqtl.file)
cat("Finish reading mQTL results\n")

cis.res = mqtl.res$cis_eqtls
cpg.list = unique(cis.res$gene)

# get chr id
file.name = basename(mqtl.file)
chr.id <- as.numeric(str_extract(file.name, "(?<=chr)\\d+"))

##### Process SNP data
geno = BEDMatrix(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id),simple_names=TRUE) %>% as.matrix() %>% PatchUp() %>% scale()
geno = geno[,colSums(is.na(geno)) == 0]
bim = fread(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id,".bim")) %>% as.data.frame()
colnames(bim) = c("CHR","SNP","V3","P0","A0","A1")
# bim should have the same snps with geno
bim = bim[bim$SNP %in% colnames(geno),]
cat("Finish reading SNP data\n")

##### Read methylation data
wgbs.file = list.files(dat.wgbs.dir,pattern=paste0("chunk_",indx1,".RDS"))
allpheno = readRDS(paste0(dat.wgbs.dir, wgbs.file))
rownames(allpheno) = paste0("CpG", allpheno$Index)
extract = allpheno[, 6:ncol(allpheno)]
annot = allpheno[,c(1,2,4,5)]
extract = t(extract)

##### read cov
cov = read.table("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/cov.txt") %>% as.matrix()

##### Run SuSIE
res.list = list()
for (cpg in cpg.list){
    pheno.used = extract[, cpg]
    cpg.pos = annot[cpg, "Pos38"] %>% as.numeric()
    # extract the geno
    bim.used = bim[bim[,"P0"] < cpg.pos + 1e6 &bim[,"P0"] > cpg.pos - 1e6,]
    if(dim(bim.used)[1] <10) {
        next
    }
    snps.used = bim.used[,"SNP"]
    geno.used = geno[,colnames(geno) %in% snps.used]

    # adjust covariates
    adjusted.pheno <- lm(pheno.used ~ cov)$residuals
    # Build SuSIE model
    susie_obj = suppressMessages(susieR::susie(
        X = geno.used,
        y = adjusted.pheno,
        L = 10,
        scaled_prior_variance = 0.1,
        estimate_residual_variance = TRUE,
        estimate_prior_variance = TRUE,
        verbose = FALSE,
        min_abs_corr = 0
    ))
    res.list[[cpg]] = susie_obj
    cat("Finish running SuSIE for CpG", cpg, "\n")
}
saveRDS(res.list, paste0(save.dir, "susie_res_", indx1, ".RDS"))
cat("Finish All")