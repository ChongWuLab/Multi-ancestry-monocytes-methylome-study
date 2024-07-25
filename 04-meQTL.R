
suppressWarnings(library(data.table))
suppressWarnings(library(BEDMatrix))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressWarnings(library(stringr))
library(dplyr)

setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
#print(paste("Current working directory:", getwd()))

source("support_MWAS.R")
#read the data
dat.wgbs.dir = "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"
#data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_WGS/"
data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"
res.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/03-mQTL-allcis/"
res.cov.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/"



args = commandArgs(TRUE) 
indx1 = as.numeric(args[[1]]) # 1:12884

#allfiles = list.files(dat.wgbs.dir)

wgbs.file = list.files(dat.wgbs.dir,pattern=paste0("chunk_",indx1,".RDS"))
file_name = wgbs.file
message("Working on the methlytion file: ",file_name)

chr.id <- as.numeric(str_extract(file_name, "(?<=chr)\\d+"))

### Process SNP data
#geno = BEDMatrix(paste0(data.wgs.dir,"U19_trainSet_AFA_allSNP.QC.CHR",chrid),simple_names=TRUE)
geno = BEDMatrix(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id),simple_names=TRUE) %>% as.matrix() %>% PatchUp() %>% t()


bim = fread(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id,".bim")) %>% as.data.frame()
colnames(bim) = c("CHR","SNP","V3","P0","A0","A1")


a1 = bim$A0
a2 = bim$A1


##### Read methylation data
allpheno = readRDS(paste0(dat.wgbs.dir, wgbs.file))
#rownames(allpheno) = paste0("CpG", allpheno$Index)
extract = allpheno[, 5:ncol(allpheno)]
annot = allpheno[,1:4]

#extract = t(extract)
#extract = as.matrix(extract)


### read cov
cov = read.table("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/cov.txt") %>% as.matrix() %>% t()


# get position matrix
snpspos = bim[,c(2,1,4)]
colnames(snpspos) = c("snpid","chr", "pos")
genepos = annot
colnames(genepos) = c("geneid","chr", "left", "right")
genepos$geneid = rownames(genepos)
genepos$chr = as.numeric(sub("chr","",genepos$chr))

#remove NA
idx = which(apply(is.na(genepos),1,any))
if (length(idx) > 0) {
    genepos = genepos[-idx,]
    extract = extract[-idx,]}
idx = which(apply(is.na(snpspos),1,any))
if (length(idx) > 0) {
    snpspos = snpspos[-idx,]
    geno = geno[-idx,]}



#mQTL analysis
library(MatrixEQTL)
#set some parameters
slicesize = 1000
useModel = modelLINEAR
errorCovariance = numeric()
cis_threshold = 1
trans_threshold = 0
cisDist = 1e6


#Prepare data
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = slicesize
snps$CreateFromMatrix(geno)

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
#cvrt$fileOmitCharacter = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = slicesize
cvrt$CreateFromMatrix(cov)

gene = SlicedData$new()
gene$fileDelimiter = "\t"
#gene$fileOmitCharacter = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = slicesize
gene$CreateFromMatrix(as.matrix(extract))

# Ensure column names are characters for the matching operation
if (nrow(cvrt)>0){
  ids = Reduce(intersect, list(snps$columnNames, gene$columnNames, cvrt$columnNames))
  snps$ColumnSubsample(match(ids, snps$columnNames))
  gene$ColumnSubsample(match(ids, gene$columnNames))
  cvrt$ColumnSubsample(match(ids, cvrt$columnNames))
  stopifnot(all(snps$columnNames==gene$columnNames))
  stopifnot(all(snps$columnNames==cvrt$columnNames))
} else{
        ids = Reduce(intersect, list(snps$columnNames, gene$columnNames, cvrt$columnNames))
        snps$ColumnSubsample(match(ids, snps$columnNames))
        gene$ColumnSubsample(match(ids, gene$columnNames))
        stopifnot(all(snps$columnNames==gene$columnNames))}

message("Performing matrixeqtl analysis on ", length(snps$columnNames), " samples")

me <- Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = NULL,
        pvOutputThreshold = trans_threshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = NULL,
        pvOutputThreshold.cis = cis_threshold,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)

cat("Finish mQTL analysis.\n")


outfile = paste0(res.dir,"AFA_chr",chr.id,"_",indx1,".RDS")

cis_eqtls = me$cis$eqtls
trans_eqtls = me$trans$eqtls
saveRDS(list(cis_eqtls = cis_eqtls, trans_eqtls = trans_eqtls, cpg_list = rownames(extract)), outfile)


cat("The number of cis-meQTL is: ", nrow(cis_eqtls), "\n")
cat("The number of trans-meQTL is: ", nrow(trans_eqtls))














