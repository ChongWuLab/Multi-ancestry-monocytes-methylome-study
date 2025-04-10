library(parallel)
library(data.table)
##### merge all the chunks for all traits
res.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/08-MR/"
save.dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/"
dirs = list.dirs(res.dir, recursive = FALSE)

read_file = function(file){
    res <- readRDS(file)
    return(res)
}

for (dir in dirs){
    files <- list.files(dir, pattern = "\\.RDS$", full.names = TRUE)
    no_cores <- detectCores() - 1
    res.tables <- mclapply(files, read_file, mc.cores = no_cores)
    res_table <- rbindlist(res.tables)
    basename <- basename(dir)
    saveRDS(res_table, file = paste0(res.dir, basename, ".RDS"))
}
read_file = function(file){
    res <- readRDS(file)
    basename = basename(file)
    basename = gsub(".RDS","",basename)
    res$file = basename
    return(res)
}
files = list.files(res.dir, pattern = "\\.RDS$", full.names = TRUE)
res.tables <- lapply(files, read_file)
res.table <- rbindlist(res.tables)

res.table$MR.p.fdr = p.adjust(res.table$MR.p, method = "fdr")
res.table = res.table[order(res.table$MR.p.fdr),]
res.table = res.table[res.table$MR.p.fdr < 0.05,] # threshold 0.000327
saveRDS(res.table, file = paste0(res.dir, "AFA-MR-MVP-FDR.RDS"))

res.table$MR.p.bf = p.adjust(res.table$MR.p, method = "bonferroni")
res.table = res.table[order(res.table$MR.p.bf),]
res.table = res.table[res.table$MR.p.bf < 0.05,] # threshold: 1.2e-09
saveRDS(res.table, file = paste0(save.dir, "AFA-MR-MVP-bonferroni.RDS"))

##### merge the result for coloc 
# EUR
res.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/08a-coloc/"
dirs = list.dirs(res.dir, recursive = FALSE)

read_file = function(file){
    res <- readRDS(file)
    # Check if the data frame has zero rows
    if (nrow(res) == 0) {
        return(NULL)  # Skip this file by returning NULL
    }
    return(res)
}

for (dir in dirs) {
    files <- list.files(dir, pattern = "\\.RDS$", full.names = TRUE)
    no_cores <- parallel::detectCores() - 1
    res.tables <- parallel::mclapply(files, read_file, mc.cores = no_cores)
    # Filter out NULL elements before combining
    res.tables <- Filter(Negate(is.null), res.tables)
    if (length(res.tables) > 0) {  # Combine only if there are non-NULL data frames
        res_table <- data.table::rbindlist(res.tables)
        basename <- basename(dir)
        res_table$file <- basename
        saveRDS(res_table, file = paste0(res.dir, basename, ".RDS"))
    }
}

files = list.files(res.dir, pattern = "\\.RDS$", full.names = TRUE)
res.tables <- lapply(files, read_file)
res.table <- rbindlist(res.tables)
res.table = res.table[res.table$coloc.abf.H4 > 0.8,]
#res.table = res.table[res.table$coloc.abf.H4 > 0.7,]
saveRDS(res.table, file = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/CAU/CAU-coloc-PPH4-07.RDS")


# AFR
res.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/09-coloc/"
dirs = list.dirs(res.dir, recursive = FALSE)

read_file <- function(file) {
    res <- readRDS(file)
    # Check if the data frame has zero rows
    if (nrow(res) == 0) {
        return(NULL)  # Skip this file by returning NULL
    }
    return(res)
}

for (dir in dirs) {
    files <- list.files(dir, pattern = "\\.RDS$", full.names = TRUE)
    no_cores <- parallel::detectCores() - 1
    res.tables <- parallel::mclapply(files, read_file, mc.cores = no_cores)
    # Filter out NULL elements before combining
    res.tables <- Filter(Negate(is.null), res.tables)
    if (length(res.tables) > 0) {  # Combine only if there are non-NULL data frames
        res_table <- data.table::rbindlist(res.tables)
        basename <- basename(dir)
        res_table$file <- basename
        saveRDS(res_table, file = paste0(res.dir, basename, ".RDS"))
    }
}

files = list.files(res.dir, pattern = "\\.RDS$", full.names = TRUE)
res.tables <- lapply(files, read_file)
res.tables <- Filter(Negate(is.null), res.tables)
if (length(res.tables) > 0) {
    res.table <- data.table::rbindlist(res.tables)
    res.table <- res.table[res.table$coloc.abf.H4 > 0.7,]
    saveRDS(res.table, file = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/AFA-coloc-PPH4-07.RDS")
}

##### Start from here
##### Analysis for one trait and find if they have common CpGs with the MWAS result
library(dplyr)
library(data.table)
library(ggplot2)
### Process MWAS result
res.A = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA-association-MVP-bf.RDS")
res.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU-association-MVP-bf.RDS")
res.A <- res.A$res.table
res.C <- res.C$res.table

# apply nsnp cutoff
res.A <- res.A[res.A$n_used_snp > 10, ] # 2998
res.C <- res.C[res.C$n_used_snp > 10, ] # 31336

# remove duplicate
res.A <- res.A[!duplicated(res.A),] # 2998
res.C <- res.C[!duplicated(res.C),]

#res.A = res.A[, c("trait", "gene")]
#res.C = res.C[, c("trait", "gene")]


# Get traits information
Fingene_id <- fread("/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/dbGAP_Upload_Tracking_JUNE2023_selected_v2.csv", skip = 1, data.table = FALSE)
dim(Fingene_id)
Fingene_id <- Fingene_id[,c("Trait", "Ethnicity", "Association Filename","Broad Category")]
Fingene_id$`Association Filename` <- gsub(".txt.gz","",Fingene_id$`Association Filename`)
Fingene_id$`Association Filename` <- gsub(".gz","",Fingene_id$`Association Filename`)

# Merge trait info with res tables
res.A <- merge(res.A, Fingene_id, by.x = "trait", by.y = "Association Filename") # 2998
res.C <- merge(res.C, Fingene_id, by.x = "trait", by.y = "Association Filename") # 31264


# remove duplicates
res.A <- res.A[!duplicated(res.A), ] # 2,998
res.C <- res.C[!duplicated(res.C), ] # 31,334

# Merge res.A and res.C using Trait and gene
combined_res <- merge(res.A,res.C,by = c("Trait","gene", "Broad Category"), suffixes = c(".AFR",".EUR")) # 559 associations
combined_res = combined_res[!duplicated(combined_res),] # 559

combined_res$Inconsistency <- ifelse((combined_res$Zscore.AFR > 0 & combined_res$Zscore.EUR < 0) | 
                                     (combined_res$Zscore.AFR < 0 & combined_res$Zscore.EUR > 0),
                                     "Inconsistent pairs", as.character(combined_res$`Broad Category`))

combined_res$Inconsistency <- factor(combined_res$Inconsistency, levels = c("Addiction", "Anthropometry", 
                                                                            "CVD", "Lipids", "Metabolic", "Renal",
                                                                             "Inconsistent pairs"))
combined_res = combined_res[combined_res$Inconsistency != "Inconsistent pairs",] # 537

### Process MR data
save.dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/"
MR.data.CAU = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/CAU/CAU-MR-MVP-bonferroni.RDS")
MR.data.AFA = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/AFA-MR-MVP-bonferroni.RDS")

### Merge MWAS and postMWAS data
MR.data.AFA = MR.data.AFA[,c("file","CpG","MR.p.bf")] #1804
MR.data2.AFA = inner_join(MR.data.AFA, Fingene_id, by = c("file" = "Association Filename"))
MR.data.CAU = MR.data.CAU[,c("file","CpG","MR.p.bf")] #1804
MR.data2.CAU = inner_join(MR.data.CAU, Fingene_id, by = c("file" = "Association Filename"))

combined_res = combined_res[,c("Trait","gene","Broad Category")]
colnames(combined_res) = c("Trait","CpG","Broad Category")
#merged.data = merge(combined_res, MR.data2, by.x = c("Trait","gene"), by.y = c("Trait","CpG"))
# use inner join
merged.data = inner_join(combined_res, MR.data2, by = c("Trait","CpG"))

# To check if there are replications of Trait and CpG
merged.data = merged.data[!duplicated(merged.data),] #464

### check all the MR result EA 
#MR.res.A <- merge(res.A, MR.data2, by.x = c("Trait", "gene"), by.y = c("Trait", "CpG"), all = FALSE) # 879
MR.res.C = merge(res.C, MR.data2.CAU, by.x = c("Trait", "gene"), by.y = c("Trait", "CpG"), all = FALSE) # 21,483
MR.res.A = merge(res.A, MR.data2.AFA, by.x = c("Trait", "gene"), by.y = c("Trait", "CpG"), all = FALSE) # 1,166

############################################################
##### begin analysis of coloc
############################################################
coloc.data.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/08a-coloc/CAU-coloc.RDS")
coloc.data.A = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/AFA-coloc.RDS")

merged.C = merge(res.C, coloc.data.C, by.x = c("trait","gene"), by.y = c("file","CpG")) # 5328
merged.A = merge(res.A, coloc.data.A, by.x = c("trait","gene"), by.y = c("file","CpG")) # 243

#check coloc for consistently identified MWAS across both ancestries
merged.data = merge(combined_res, coloc.data.C, by.x = c("trait.EUR","gene"), by.y = c("file","CpG")) # 104
merged.data2 = merge(combined_res, coloc.data.A, by.x = c("trait.AFR","gene"), by.y = c("file","CpG")) # 106

# using PPH4 = 0.7
coloc.data.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/08a-coloc/CAU-coloc-PPH4-07.RDS")
coloc.data.A = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/AFA-coloc-PPH4-07.RDS")
merged.C = merge(res.C, coloc.data.C, by.x = c("trait","gene"), by.y = c("file","CpG")) # 11517
merged.A = merge(res.A, coloc.data.A, by.x = c("trait","gene"), by.y = c("file","CpG")) # 281
merged.data = merge(combined_res, coloc.data.C, by.x = c("trait.EUR","gene"), by.y = c("file","CpG")) # 230


### get MR and MWAS cpg list
# For AFA
MR.res = readRDS(paste0(save.dir, "AFA-MR-MVP-bonferroni.RDS"))
MR.cpg = unique(MR.res$CpG)
MWAS.res = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA-association-MVP-bf.RDS")
MWAS.res = MWAS.res$res.table
MWAS.cpg = unique(MWAS.res$gene)
MWAS_and_MR_cpglist = union(MR.cpg, MWAS.cpg)
save(MWAS_and_MR_cpglist, file = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/MWAS_and_MR_cpglist.RData")

# for EUR
MR.res = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/CAU/CAU-MR-MVP-bonferroni.RDS")
MR.cpg = unique(MR.res$CpG)
MWAS.res = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU-association-MVP-bf.RDS")
MWAS.res = MWAS.res$res.table
MWAS.cpg = unique(MWAS.res$gene)
MWAS_and_MR_cpglist = union(MR.cpg, MWAS.cpg)


# Prepare supplemnary table
# For AFA
all.res.A = res.A
all.res.A <- left_join(all.res.A, MR.data2.AFA, by = c("trait" = "file", "gene" = "CpG"))
all.res.A = all.res.A[!duplicated(all.res.A),]
all.res.A = left_join(all.res.A, coloc.data.A, by = c("trait" = "file", "gene" = "CpG"))
all.res.A = all.res.A[!duplicated(all.res.A),]
write.csv(all.res.A, file = "/home/wzhang24/mQTL_project/Results/06-association/AFA-all-res.csv", row.names = FALSE)

# for CAU
all.res.C = res.C
all.res.C <- left_join(all.res.C, MR.data2.CAU, by = c("trait" = "file", "gene" = "CpG"))
all.res.C = all.res.C[!duplicated(all.res.C),]
all.res.C = left_join(all.res.C, coloc.data.C, by = c("trait" = "file", "gene" = "CpG"))
all.res.C = all.res.C[!duplicated(all.res.C),]
write.csv(all.res.C, file = "/home/wzhang24/mQTL_project/Results/06-association/CAU-all-res.csv", row.names = FALSE)