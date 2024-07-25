library(data.table)
library(coloc)
#library(TwoSampleMR)
library(BEDMatrix)
#library(magrittr)
setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
source("support.R")
# Set data path
gwassum.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/processed/"
ldref = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"
mqtl.dir = "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/meQTL_chunk_results/AFR/"
data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"
out.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/09-coloc/"

args <- commandArgs(trailingOnly = TRUE)

# The first argument is stored in args[1]
file.id <- as.numeric(args[1])
indx1 <- as.numeric(args[2])

##### Process meQTL data
mqtl.file <- list.files(mqtl.dir, pattern = paste0("_", indx1, ".RDS"))
mqtl.data <- readRDS(paste0(mqtl.dir, mqtl.file))
allcis <- mqtl.data$cis_eqtls
colnames(allcis) <- c("SNP", "CpG", "Z", "P", "FDR", "BETA")
cpg.list <- mqtl.data$cpg_list

# check if the cpg is in the MWAS and MR cpg list
load("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/MWAS_and_MR_cpglist.RData")
cpg.list <- intersect(cpg.list, MWAS_and_MR_cpglist)

if (length(cpg.list) == 0) {
    cat("No cpgs for this file!")
    quit()
}
cat("Number of cpgs for this file: ", length(cpg.list), "\n")

allcis = allcis[allcis$CpG %in% cpg.list,]
# get chromosome
chr.id <- sub(".*chr(\\d+).*", "\\1", mqtl.file) %>% as.numeric()


# Read bim for snp info
bim = fread(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id,".bim"), data.table = FALSE)
rownames(bim) <- bim$V2
geno = BEDMatrix(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id),simple_names=TRUE) 


##### Prepare sumstats
info = fread("/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/dbGAP_Upload_Tracking_JUNE2023_selected_v2.csv",skip = 1,data.table=F)
info = info[info$Ethnicity=="AFR",]
# remove duplicate
info = info[!duplicated(info$`Association Filename`),]
row.names(info) = 1:nrow(info)
sumstatfile = info[file.id,"Association Filename"]
sumstatfile = gsub(".txt.gz",".sumstats",sumstatfile)
sumstatfile = gsub(".gz",".sumstats",sumstatfile)
sumstatfile = gsub(".txt",".sumstats",sumstatfile)

cat("sumstats file id: ",file.id,"\n")
cat("Chromosome: ",chr.id,"\n")
cat("Current sumstats file:",sumstatfile,"\n")

sumstats.org = fread(paste0(gwassum.dir,sumstatfile)) %>% as.data.frame()
if(sumstatfile=="MDD.EUR.MVP_UKBB.NatNeuro2021.sumstats") {
    sumstats.org$SAMPLESIZE = 286821
} else if (sumstatfile=="MDD.EUR.MVP.NatNeuro2021.sumstats") {
    sumstats.org$SAMPLESIZE = 250215
} else if (sumstatfile=="MDD.AFR.MVP.NatNeuro2021.sumstats") {
    sumstats.org$SAMPLESIZE = 59600
} else if (sumstatfile == "dbGAP_GAD2eur.sumstats") {
    sumstats.org$SAMPLESIZE = 175163
} else if (sumstatfile == "dbGAP_GAD2afr.sumstats") {
    sumstats.org$SAMPLESIZE = 24448
}

tmp.colname = colnames(sumstats.org)
tmp.colname[tmp.colname=="SAMPLESIZE"] = "N"
colnames(sumstats.org) = tmp.colname
n.sumstats = floor(mean(sumstats.org$N))
if(sum(grepl("chr",sumstats.org[,"CHROMOSOME"]))>100) {
    
} else {
    sumstats.org[,"CHROMOSOME"] = paste0("chr", sumstats.org[,"CHROMOSOME"])
}

sumstats <- sumstats.org[sumstats.org$CHROMOSOME == paste0("chr", chr.id), ]
sumstats = sumstats[rowSums(is.na(sumstats))==0, ]


# Drop duplicated SNPs if any
sumstats.dup <- sumstats$SNP %>% duplicated()
if (sumstats.dup %>% sum() != 0) {
    sumstats <- sumstats[!sumstats.dup, ]
}


#read LD reference panel; align the data and calculate the LD matrix
reference.bim <- paste0(ldref, chr.id, ".bim") %>% fread(., data.table = FALSE)
reference.bed <- paste0(ldref, chr.id) %>% BEDMatrix(., simple_names = TRUE)

reference.dup <- reference.bim$V2 %>% duplicated()
if (reference.dup %>% sum() != 0) {
    reference.bim <- reference.bim[!reference.dup, ]
    reference.bed <- reference.bed[, !reference.dup]
}

# Re-assign row names
rownames(sumstats) <- sumstats$SNP
#rownames(allcis)   <- allcis$snps
rownames(reference.bim) <- reference.bim$V2

# Get the unique snps and find the intersection with sumstats
unique_snps <- intersect(unique(allcis$SNP), reference.bim$V2) %>% intersect(., rownames(bim))
allcis <- allcis[allcis$SNP %in% unique_snps,]
cur.sumstats = sumstats[sumstats$SNP %in% unique_snps,]
cur.reference.bim = reference.bim[unique_snps,]
cur.reference.bed <- reference.bed[, unique_snps] %>% PatchUp()
cur.bim <- bim[unique_snps,]
cur.bed <- geno[,unique_snps] %>% PatchUp()

# calculate MAF for current bim
cur.bim$MAF = colSums(cur.bed)/(dim(cur.bed)[1]*2)

# initiate ressult table
res_table <- data.frame(
    CpG = character(),
    coloc.abf.H3 = numeric(),
    coloc.abf.H4 = numeric(),
    stringsAsFactors = FALSE
)

trait.name = gsub(".sumstats", "", sumstatfile)
save.dir <- paste0(out.dir, trait.name, "/")
#temp.outdir = paste0(save.dir,"tmp_",indx1,"/")
# create temp.outdir if not exist
if (!dir.exists(save.dir)) {
    dir.create(save.dir, recursive = TRUE)
}

j = 1
cat("==================================================\n")
for (i in 1:length(cpg.list)) {
    cur.cpg <- cpg.list[i]
    used.SS <- allcis[allcis$CpG == cur.cpg,]
    snps.common <- intersect(used.SS$SNP, cur.sumstats$SNP)
    if (length(snps.common) == 0) {
        next
    }
    # Drop duplicated SNPs if any
    SS.dup <- used.SS$SNP %>% duplicated()
    if (SS.dup %>% sum() != 0) {
        used.SS <- used.SS[!SS.dup, ]
    }
    rownames(used.SS) <- used.SS$SNP

    # Extract the common SNPs
    used.SS <- used.SS[snps.common,]
    used.sumstats <- cur.sumstats[snps.common,]
    used.reference.bed <- cur.reference.bed[, snps.common]
    used.reference.bim <- cur.reference.bim[snps.common,]
    used.bim <- cur.bim[snps.common,]
    # Allele-flip the phenotype ss and weight w.r.t. reference panel
    qc.1 <- allele.qc(
        used.sumstats$A2,
        used.sumstats$A1,
        used.reference.bim$V5,
        used.reference.bim$V6
    )

    if(!"MAF" %in% colnames(used.sumstats)) {
        used.sumstats$MAF = colSums(used.reference.bed)/(dim(used.reference.bed)[1]*2)
    } else {
        used.sumstats$MAF[qc.1$flip] = 1 - used.sumstats$MAF[qc.1$flip]
    }

    used.sumstats$Z[qc.1$flip] <- -1 * used.sumstats$Z[qc.1$flip]
    used.sumstats$BETA[qc.1$flip] <- -1 * used.sumstats$BETA[qc.1$flip]
    if (!"SE" %in% colnames(used.sumstats)) {
        used.sumstats$SE <- used.sumstats$BETA / used.sumstats$Z
    }

    qc.2 <- allele.qc(
        used.bim$V5,
        used.bim$V6,
        used.reference.bim$V5,
        used.reference.bim$V6
    )

    used.SS$BETA[qc.2$flip] <- -1 * used.SS$BETA[qc.2$flip]
    #used.bim$MAF[qc.2$flip] = 1  - used.bim$MAF[qc.2$flip]
    
    used.SS$SE <- used.SS$BETA / used.SS$Z
    #SS.original["correlation"] <- SS.original$Z / sqrt(SS.original$Z ^ 2 + SS.original$N - 1)

    keep <- qc.1$keep & qc.2$keep
    if (sum(!keep) > 0) {
        used.sumstats <- used.sumstats[keep, ]
        used.SS   <- used.SS[keep, ]
        used.reference.bim <- used.reference.bim[keep, ]
        used.reference.bed <- used.reference.bed[, keep]
        used.bim <- used.bim[keep, ]
    }

    matrix.LD <- cor(used.reference.bed)

    ##########################################################################
    ## prepare the data for coloc
    ##########################################################################
    if("BETA" %in% colnames(used.sumstats)) {
        D1 = list(beta = used.sumstats$BETA, varbeta = used.sumstats$SE^2, snp = used.sumstats$SNP, type ="cc", position = used.bim$V4, N = n.sumstats, MAF =  used.sumstats$MAF, LD = matrix.LD)

    } else {
        D1 = list( snp = used.sumstats$SNP, type ="cc", position = used.bim$V4, N = n.sumstats, MAF =  used.sumstats$MAF, pvalues = used.sumstats$P,  LD = matrix.LD)

    }

    D2 = list(beta = used.SS$BETA, varbeta = used.SS$SE^2, snp= used.SS$SNP, type ="quant", position = used.bim$V4, N = 298, MAF =  used.bim$MAF, LD = matrix.LD)

    abf.res = suppressMessages(suppressWarnings(coloc.abf(dataset1=D1, dataset2=D2)))

    # Save to the result table
    res_table[j, "CpG"] = cur.cpg
    res_table[j, "coloc.abf.H3"] = abf.res[[1]]["PP.H3.abf"]
    res_table[j, "coloc.abf.H4"] = abf.res[[1]]["PP.H4.abf"]
    cat(j, "\n")
    j <- j + 1
}

# filter res_table with coloc.abf.H3 > 0.8 or coloc.abf.H4 > 0.8
#res_table = res_table[res_table$coloc.abf.H4 > 0.8,]

# Save the result table
if (nrow(res_table) != 0) {
    cat("Number of significant results: ", nrow(res_table), "\n")
    res.table.file <- paste0(save.dir, "res_coloc_", indx1, ".RDS")
    saveRDS(res_table, res.table.file)
    cat("Done!")
} else {
    cat("No significant results!")
}











