#load packages
library(data.table)
library(dplyr)
library(stringr)
library(glmnet)
library(BEDMatrix)
setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA")
source("support_MWAS.R")

# obtain the data
args <- commandArgs(TRUE)
indx1 <- as.numeric(args[[1]])

#read the data
dat.wgbs.dir = "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"
#data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_WGS/"
data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"
#mQTL.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/03-mQTL-V2/"
#clumped.data.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/04-clumping/"
save.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/04-modelling-V2-nc/"
res.cov.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/"
h2.dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-AFA/"


wgbs.file = list.files(dat.wgbs.dir,pattern=paste0("chunk_",indx1,".RDS"))
file_name <- wgbs.file


# Define the regular expression pattern to match the number following "chr"
chrid <- as.numeric(str_extract(file_name, "(?<=chr)\\d+"))

#check if the file already exists
#if (file.exists(paste0(save.dir,"chr",chrid,"_",indx1,".RData"))) {
#    cat("File already exists, skip.\n")
#    quit()
#}

# read h2 file
h2_table = readRDS(file.path(h2.dir,"all_h2_filter.RDS"))
cpg_list = h2_table$CpG

#load reference genotype
ldref = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"
reference.bim <- paste0(ldref, chrid, ".bim") %>% fread(., data.table = FALSE) %>% as.data.frame()
reference.bed <- paste0(ldref, chrid) %>% BEDMatrix(., simple_names = TRUE) %>% as.matrix()
#check if there is duplicate
reference.dup <- reference.bim$V2 %>% duplicated()
if (reference.dup %>% sum() != 0) {
    reference.bim <- reference.bim[!reference.dup, ]
    reference.bed <- reference.bed[, !reference.dup]
}

rownames(reference.bim) <- reference.bim$V2
colnames(reference.bim) <- c("CHR","SNP","V3","Position","A1","A2")


# load geno
geno = BEDMatrix(paste0(data.wgs.dir,"AFA_1000G_chr",chrid),simple_names=TRUE) %>% as.matrix() %>% PatchUp() %>% scale()
geno = geno[,colSums(is.na(geno)) == 0]
bim = fread(paste0(data.wgs.dir,"AFA_1000G_chr",chrid,".bim"),data.table=FALSE)
colnames(bim) = c("CHR","SNP","V3","P0","A0","A1")
# bim should have the same snps with geno
bim = bim[bim$SNP %in% colnames(geno),]

#remove ambiguous snps
a1 = bim$A0
a2 = bim$A1
keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
bim = bim[keep,]
geno = geno[,keep]


#####process pheno
allpheno = readRDS(paste0(dat.wgbs.dir, wgbs.file))
allpheno = allpheno[rownames(allpheno) %in% cpg_list,]

pheno = allpheno[, 6:ncol(allpheno)] %>% t() %>% as.matrix()
annot = allpheno[,1:5]


### read cov
cov = read.table("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/cov.txt")

CpG_mat = data.frame(CpG = character(), chr = character(), p.0 = character(), p.1 = character(), Corr = numeric(),R2=numeric(),pval=numeric(),stringsAsFactors = FALSE)
MatrixLD = list()
SNPs = list()
idx=0 
cat("There are ",dim(pheno)[2]," CpGs to build models.\n")

# inidividal model:
for(j in 1:dim(pheno)[2]) {
  cpg.pos = annot[j,"Pos38"]
  
  # extract the snps 
  cpg.name = colnames(pheno)[j]

  # extract the geno
  bim.used = bim[bim[,"P0"] < cpg.pos + 5e5 &bim[,"P0"] > cpg.pos - 5e5,]
  if(dim(bim.used)[1] <10) {
      next
  }
  snps.used = bim.used[,"SNP"]
  geno.used = geno[,colnames(geno) %in% snps.used]

  # extract the pheno
  pheno.used = as.matrix(pheno[,j])
  cov$pheno.used = pheno[,j]
  model <- lm(pheno.used ~ ., data = cov)
  pheno.used = model$residuals
  
   # check if have erros when running the model
  error_occurred <- FALSE
  tryCatch({
      stdPWAS.weight = weights.stdPWAS(geno.used, as.matrix(pheno.used))
  }, error=function(e){
      cat("Warning: Cannot run ", j, cpg.name,"Error detail:",conditionMessage(e),"\n")
      # print the system error
  
      error_occurred <- TRUE
  })

  if (error_occurred) {
    next  # Skip the rest of this iteration and move to the next
  }
  
  cv.performance = stdPWAS.weight[[1]]
  weights = stdPWAS.weight[[2]]
  #End with tab
  R2 = cv.performance[1,1]
  corr = cv.performance[3,1]
  pval = cv.performance[2,1]
  cat("Finish ",j," ", cpg.name ,"R2:",round(R2,4)," nonZeros: ",sum(weights!=0)," ")
  # Check if cv.performance[1,1] is NA or < 0.01
  # check if number of weights != 0 > 1
  if (!is.na(cv.performance[1,1]) & cv.performance[1,1] >= 0.01 & sum(weights!=0) >= 2) {
      idx = idx + 1

      rownames(bim.used) = bim.used$SNP
      snps = bim.used
      #colnames(snps) = c("CHR","SNP","V3","P0","A0","A1")
      
      snps$weight = weights
      #keep only weights != 0
      keep = which(snps$weight!=0)
      used.snps = snps[keep,]

      R2 = cv.performance[1,1]
      corr = cv.performance[3,1]
      pval = cv.performance[2,1]
      CpG_mat[idx,] = c(cpg.name,chrid,cpg.pos,annot[cpg.name,"Pos38.right"],corr,R2,pval)

      snps.common <- intersect(used.snps$SNP, reference.bim$SNP)
      #reference.bim.keep <- reference.bim[snps.common, ]
      reference.bed.keep <- reference.bed[, snps.common] %>% PatchUp() %>% scale()
      matrix.LD <- t(reference.bed.keep) %*% reference.bed.keep / (nrow(reference.bed.keep) - 1)

      MatrixLD[[cpg.name]] <- matrix.LD
      SNPs[[cpg.name]] <- used.snps

      cat("Saved.\n")
  }
  else{
      cat("Not saved.\n")
  }
}


cat("CpG_mat contains ",dim(CpG_mat)[1]," CpGs in this chunk.\n")
if (dim(CpG_mat)[1] > 0) {
    save(CpG_mat, MatrixLD, SNPs, file = paste0(save.dir,"chr",chrid,"_",indx1,".RData"))
    cat("Finish saving chunk",indx1,"\n")
}else{
    cat("No CpGs saved in chunk",indx1,"\n")
}



