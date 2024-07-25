library(data.table)
library(BEDMatrix)
library(dplyr)
library(MASS)

#library(purrr)


# Set the working directory to the directory containing the R script
setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
  
source("BLISSAssociation_Support.R")

# Set data path
gwassum.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/processed/"
ldref = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"
out.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/MVP/"
#res.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/04-modelling-V2-nc/"
res.dir = "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/models/AFA/"
args <- commandArgs(trailingOnly = TRUE)

# The first argument is stored in args[1]
file.id <- args[1] # file id of sumstats 1-42
#file.id = 11
chromosome = args[2] # chr 
chromosome = as.numeric(chromosome)
file.id = as.numeric(file.id)
#file.id=33


cat("sumstats file id: ",file.id,"\n")
cat("Chromosome: ",chromosome,"\n")



# obtain sample size
info = fread("/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/dbGAP_Upload_Tracking_JUNE2023_selected_v2.csv",skip = 1,data.table=FALSE)
info = info[info$Ethnicity=="AFR",]

files = list.files(gwassum.dir)

tmp = gsub(".sumstats",".txt.gz",files)
indx = info[,"Association Filename"] %in% tmp

tmp = gsub(".sumstats",".gz",files)
indx2 = info[,"Association Filename"] %in% tmp

tmp = gsub(".sumstats",".txt",files)
indx3 = info[,"Association Filename"] %in% tmp

indx = indx | indx2 | indx3

info = info[indx,]
rownames(info) = 1:nrow(info)
sumstatfile = info[file.id,"Association Filename"]
sumstatfile = gsub(".txt.gz",".sumstats",sumstatfile)
sumstatfile = gsub(".gz",".sumstats",sumstatfile)
sumstatfile = gsub(".txt",".sumstats",sumstatfile)
cat("Current sumstats file:",sumstatfile,"\n")

# load GWAS summary data
sumstats.org = fread(paste0(gwassum.dir,sumstatfile),data.table=F)


# get sample size
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
    sumstats.org[,"CHROMOSOME"] = paste0("chr",sumstats.org[,"CHROMOSOME"])
}



# 
ref = ldref
outindx = 1

reference.bim <- paste0(ref, chromosome, ".bim") %>% fread(., data.table = FALSE)
reference.bed <- paste0(ref, chromosome) %>% BEDMatrix(., simple_names = TRUE) 
n_sample <- nrow(reference.bed)

#check if there is duplicate
reference.dup <- reference.bim$V2 %>% duplicated()
if (reference.dup %>% sum() != 0) {
    reference.bim <- reference.bim[!reference.dup, ]
    #reference.bed <- reference.bed[, !reference.dup]
}
colnames(reference.bim) <- c("CHR","SNP","V3","Position","A1","A2")
rownames(reference.bim) <- reference.bim$SNP



sumstats <- sumstats.org[sumstats.org$CHROMOSOME == paste0("chr", chromosome), ]
sumstats = sumstats[rowSums(is.na(sumstats))==0, ]
# Drop duplicated SNPs if any
sumstats.dup <- sumstats$SNP %>% duplicated()
if (sumstats.dup %>% sum() != 0) {
    sumstats <- sumstats[!sumstats.dup, ]
}
rownames(sumstats) <- sumstats$SNP

trait.name = gsub(".sumstats","",sumstatfile)
save.dir = paste0(out.dir,trait.name)
#create save.dir if not exist, otherwise, skip
if (!dir.exists(save.dir)) {
    dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
}



# obtain all the chunks
list_of_files <- list.files(res.dir, pattern = paste0("chr",chromosome,"_"))
outres_list <- list()
cat("There are in total",length(list_of_files)," chunks for chromosome ",chromosome,"\n")
cat("--------------------------------------------------\n")
for (i in 1:length(list_of_files)){
    res.file = list_of_files[i]
    res.file.name = gsub(".RData","",res.file)
    load(paste0(res.dir,res.file))
    if (nrow(CpG_mat)==0) {
        cat("No CpG in chunk", i,res.file.name, "\n")
        next
    }
    all_data <- bind_rows(SNPs)
    #colnames(all_data) <- c("chr","SNP","V3","p0","A0","A1","weight")
    # Get unique SNP names
    unique_snps <- all_data %>%
      pull(SNP) %>%
      unique()
    cur.sumstats = sumstats[sumstats$SNP %in% unique_snps,]
    cur.reference.bim = reference.bim[unique_snps,]

    outres = data.frame(chr=character(),p0=character(),p1=character(),gene=character(),R2=numeric(),Zscore.classic=numeric(),p.classic=numeric(),beta_alt=numeric(),se_alt=numeric(),p_alt=numeric(),n_used_snp=numeric(),n_snp=numeric(),stringsAsFactors = FALSE)
    cat(i," Begin processing",res.file.name," ",nrow(CpG_mat),"CpGs in this chunk.\n")
    for(j in 1:nrow(CpG_mat)) {
        cpg.name = CpG_mat[j,"CpG"]
        snps = as.data.frame(SNPs[[cpg.name]])
        #outres[outindx,1:4] = used.lookup[used.lookup[,"gene"]==protein,c("chr","p.0","p.1","gene")]
        outres[j,1:5] = CpG_mat[j, c("chr","p.0","p.1","CpG","R2")]
        #weight = snps$weight
        outres[j,12] = dim(snps)[1]
        
        # read LD matrix
        matrix.LD <- MatrixLD[[cpg.name]]
        #used snps
        rownames(snps) = snps$SNP

        error = FALSE
        tryCatch({
            outres[j,6:11] = TestAssociation(sumstats=cur.sumstats,  
                                                   n.sumstats=n.sumstats, 
                                                   ref=cur.reference.bim, 
                                                   n_sample=n_sample,
                                                   snps=snps, 
                                                   LD = matrix.LD)
        }, error=function(e){
            #cat("Indx",outindx, cpg.name, "Can't do TestAssociation. Error Details:", conditionMessage(e), "\n")
            error=TRUE
        })
        
        #cat("Finish Indx",outindx,"\n")
        if (error) {
            outres[j,6:11] = NA
        }  
    }

    # remove NA and reset index
    outres = outres[!is.na(outres[,8]),]
    cat("   There are ",nrow(outres)," CpGs with results in this chunk.\n")
    # Add the outres data frame to the list
    if (exists("outres") && nrow(outres) > 0) {
        # Add the outres data frame to the list
        outres_list[[length(outres_list) + 1]] <- outres
    }
    else {
        cat("No result in chunk", i, res.file.name,"\n")
    }
}

# combine all data frames in the list into one data.table
final_res <- rbindlist(outres_list, use.names = TRUE, fill = TRUE)

output.name = paste0("res_AFA_MVP_file_",file.id,"_chr",chromosome,".RDS")
saveRDS(final_res, file = paste0(save.dir,"/",output.name))
cat("-----------------------------------------------\n")
cat("Finish All for", output.name)