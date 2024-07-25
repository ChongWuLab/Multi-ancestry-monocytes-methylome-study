library(data.table)
#library(coloc)
library(TwoSampleMR)
#library(BEDMatrix)
#library(magrittr)
setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
source("support.R")
# Set data path
gwassum.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/processed/"
ldref = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"
mqtl.dir = "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/meQTL_chunk_results/AFR/"
data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"
out.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/08-MR/"

args <- commandArgs(trailingOnly = TRUE)

# The first argument is stored in args[1]
file.id <- as.numeric(args[1])
indx1 <- as.numeric(args[2])

##### Process meQTL data
mqtl.file <- list.files(mqtl.dir, pattern = paste0("_", indx1, ".RDS"))
mqtl.data <- readRDS(paste0(mqtl.dir, mqtl.file))
allcis <- mqtl.data$cis_eqtls
colnames(allcis) <- c("SNP", "CpG", "Z", "P", "FDR", "BETA")

# Filter out the cis-eQTLs with pvalue > 1e-8
allcis <- allcis[allcis$P < 1e-8, ]
# Get the cpg list
cpg.list <- unique(allcis$CpG)
# get chromosome
chr.id <- sub(".*chr(\\d+).*", "\\1", mqtl.file) %>% as.numeric()


# Read bim for snp info
bim = fread(paste0(data.wgs.dir,"AFA_1000G_chr",chr.id,".bim"), data.table = FALSE)
rownames(bim) <- bim$V2


##### Prepare sumstats
info = fread("/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/dbGAP_Upload_Tracking_JUNE2023_selected_v2.csv",skip = 1,data.table=F)
info = info[info$Ethnicity=="AFR",]
# remove duplicate
info = info[!duplicated(info$`Association Filename`),]
# reset the rownames of info
rownames(info) = 1:nrow(info)
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

reference.dup <- reference.bim$V2 %>% duplicated()
if (reference.dup %>% sum() != 0) {
    reference.bim <- reference.bim[!reference.dup, ]
}

# Re-assign row names
rownames(sumstats) <- sumstats$SNP
rownames(reference.bim) <- reference.bim$V2

# Get the unique snps and find the intersection with sumstats
unique_snps <- intersect(unique(allcis$SNP), reference.bim$V2) %>% intersect(., rownames(bim))
allcis <- allcis[allcis$SNP %in% unique_snps,]
cur.sumstats = sumstats[sumstats$SNP %in% unique_snps,]
cur.reference.bim = reference.bim[unique_snps,]
cur.bim <- bim[unique_snps,]

# initiate ressult table
res_table <- data.frame(
    CpG = character(),
    MR.beta = numeric(),
    MR.se = numeric(),
    MR.p = numeric(),
    MR.niv = numeric(),
    stringsAsFactors = FALSE
)

trait.name = gsub(".sumstats", "", sumstatfile)
save.dir <- paste0(out.dir, trait.name, "/")
temp.outdir = paste0(save.dir,"tmp_",indx1,"/")
# create temp.outdir if not exist
if (!dir.exists(temp.outdir)) {
    dir.create(temp.outdir, recursive = TRUE)
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
    rownames(used.SS) <- used.SS$SNP

    # Extract the common SNPs
    used.SS <- used.SS[snps.common,]
    used.sumstats <- cur.sumstats[snps.common,]
    used.reference.bim <- cur.reference.bim[snps.common,]
    used.bim <- cur.bim[snps.common,]
    # Allele-flip the phenotype ss and weight w.r.t. reference panel
    qc.1 <- allele.qc(
        used.sumstats$A2,
        used.sumstats$A1,
        used.reference.bim$V5,
        used.reference.bim$V6
    )

    used.sumstats$Z[qc.1$flip] <- -1 * used.sumstats$Z[qc.1$flip]
    used.sumstats$BETA[qc.1$flip] <- -1 * used.sumstats$BETA[qc.1$flip]

    qc.2 <- allele.qc(
        used.bim$V5,
        used.bim$V6,
        used.reference.bim$V5,
        used.reference.bim$V6
    )

    used.SS$BETA[qc.2$flip] <- -1 * used.SS$BETA[qc.2$flip]  
    used.SS$SE <- used.SS$BETA / used.SS$Z
    #SS.original["correlation"] <- SS.original$Z / sqrt(SS.original$Z ^ 2 + SS.original$N - 1)

    keep <- qc.1$keep & qc.2$keep
    if (sum(!keep) > 0) {
        used.sumstats <- used.sumstats[keep, ]
        used.SS   <- used.SS[keep, ]
        used.reference.bim <- used.reference.bim[keep, ]
        used.bim <- used.bim[keep, ]
    }

    ##############################################
    # MR: used significant SNPs in runsusie
    # then run multivariate MR
    ##############################################
    # two way for this:

    # 1. extract the independent IVs:
    # Write ss; write this as a function
    
    ivs = clumping.iv(used.SS,0.001)
    tmp = used.SS[ivs, , drop=FALSE]
    tmp = tmp[tmp$P < 1e-8, , drop = FALSE]
    ivused = rownames(tmp)
        

    if(length(ivused)==0) {
        #cat(i, "No IVs found for CpG: ", cur.cpg, "\n")
        next
    } else if(length(ivused)==1) {
        exp = used.SS[ivused, , drop=FALSE]
        out = used.sumstats[ivused, , drop=FALSE]
        b_exp = exp$BETA
        b_out = out$BETA
        se_exp = exp$SE
        se_out = out$SE
        stdMR = unlist(mr_wald_ratio(b_exp, b_out, se_exp, se_out))
    } else {
        exp = used.SS[ivused, , drop = FALSE]
        out = used.sumstats[ivused, , drop = FALSE]
        b_exp = exp$BETA
        b_out = out$BETA
        se_exp = exp$SE
        se_out = out$SE
        stdMR = unlist(TwoSampleMR::mr_ivw(b_exp, b_out, se_exp, se_out))
    }
    
    # Save to the result table
    res_table[j, "CpG"] = cur.cpg
    res_table[j, c("MR.beta","MR.se","MR.p","MR.niv")] = stdMR[1:4]
    #res_table[j, c("MR2.beta","MR2.se","MR2.p","MR2.niv")] = corMR[1:4]
    j <- j + 1
}
# Save the result table
if (nrow(res_table) == 0) {
    cat("No CpG in chunk", indx1, "\n")
} else {
        res.table.file <- paste0(save.dir, "res_postPWAS_", indx1, ".RDS")
        cat("There are ", dim(res_table)[1], " CpGs with results in this chunk.\n")
        saveRDS(res_table, res.table.file)
}

# Delete the temp directory
unlink(temp.outdir, recursive = TRUE)

cat("Done!")










