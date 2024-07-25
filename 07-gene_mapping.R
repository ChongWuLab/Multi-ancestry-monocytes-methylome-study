
#=====================================Repeat the process for each trait===========================================================
library(dplyr)
library(purrr)
library(ACAT)
library(data.table)
res.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/MVP/"
##### process reference of ensemble gene
ref <- readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/ENSEMBL_GRch37_gene_list.RDS")
ref <- ref[,c("ensembl_gene_id","gene_biotype","chromosome_name","start_position","end_position")]
# remove duplicate
ref <- ref[!duplicated(ref$ensembl_gene_id),]
ref2 = ref[ref$gene_biotype=="protein_coding",]
setDT(ref2)

##### process cpg position
combined_pos19 <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19/combined_data.RDS")
colnames(combined_pos19) = c("gene","chr","Pos19")


res.files <- list.files(res.dir, full.names = TRUE)
all_res <- data.frame()
for (file in res.files){
  data <- readRDS(file)
  # get gene list
  cpg.table <- combined_pos19[combined_pos19$gene %in% data$gene,]
  cpg.table$chr = as.numeric(gsub("chr","",cpg.table$chr))
  setDT(cpg.table)
  result <- cpg.table[ref2, 
                    on = .(chr = chromosome_name, 
                           Pos19 >= start_position, 
                           Pos19 <= end_position),
                    .(gene, chr, Pos19 = x.Pos19, ensembl_gene_id = i.ensembl_gene_id)]
  # merge result and data by CpG
  data = data[,c("gene","chr","p_alt")]
  final_res <- merge(result, data, by="gene")
  # ACAT
  acat_results <- final_res[, .(p_acat = ACAT(p_alt)), by = ensembl_gene_id]
  acat_results$p_bf = p.adjust(acat_results$p_acat, method = "bonferroni")
  threshold = 0.05/length(acat_results$p_bf)
  cat(basename(file),sum(acat_results$p_bf < 0.05), "threshold:", threshold, "\n")
  all_res <- rbind(all_res,acat_results[acat_results$p_bf < 0.05,])
  saveRDS(acat_results,file = paste0("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/ACAT/",basename(file)))
}

read_file <- function(file) {
    res <- readRDS(file)
    res = res[res$p_bf < 0.05,]
    if (nrow(res) == 0) {
        return(NULL)  # Skip this file by returning NULL
    }
    trait_name <- tools::file_path_sans_ext(basename(file))
    res$trait <- trait_name
    return(res)
}

# List all .RDS files in the directory
files <- list.files("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/ACAT/", pattern = "\\.RDS$", full.names = TRUE)

# Read the files, add the 'trait' column, and filter out the empty data frames
res.tables <- lapply(files, read_file)
res.tables <- Filter(Negate(is.null), res.tables)

# Combine the non-empty data frames into one
if (length(res.tables) > 0) {
    combined_table <- rbindlist(res.tables)
    # Save the combined data frame
    saveRDS(combined_table, file = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA_associations_ACAT.RDS")
}



#=====================================Repeat the process for each trait===========================================================
#=====================================If we only use CpG covered by 450K/900K===========================================================
library(dplyr)
library(purrr)
library(ACAT)
library(data.table)
res.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/MVP/"
##### process reference of ensemble gene
ref <- readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/ENSEMBL_GRch37_gene_list.RDS")
ref <- ref[,c("ensembl_gene_id","gene_biotype","chromosome_name","start_position","end_position")]
# remove duplicate
ref <- ref[!duplicated(ref$ensembl_gene_id),]
ref2 = ref[ref$gene_biotype=="protein_coding",]
setDT(ref2)

##### read 450k and 900k
HM450 <- read.csv("/rsrch5/home/biostatistics/wzhang24/data/WGBS/humanmethylation450_15017482_v1-2.csv")
epic850k <- fread("/rsrch5/home/biostatistics/wzhang24/data/WGBS/EPIC-8v2-0_A1.csv", skip=7)
HM450 <- HM450[, c("CHR", "MAPINFO")]
HM450$CHR <- paste0("chr", HM450$CHR)
epic850k <- epic850k[, c("CHR", "MAPINFO")]

hg19_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19/combined_data.RDS")
hg38_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg38_annotation.RDS")
#hg38_annot <- hg38_annot[, c("CpG", "pos38")]

##### process cpg position
combined_pos19 <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19/combined_data.RDS")
colnames(combined_pos19) = c("gene","chr","Pos19")

# 450K
combined_pos19_450k <- combined_pos19 %>% inner_join(HM450, by = c("chr" = "CHR", "Pos19" = "MAPINFO"))
res.files <- list.files(res.dir, full.names = TRUE)
all_res <- data.frame()
for (file in res.files){
  data <- readRDS(file)
  # get gene list
  cpg.table <- combined_pos19_450k[combined_pos19_450k$gene %in% data$gene,]
  cpg.table$chr = as.numeric(gsub("chr","",cpg.table$chr))
  setDT(cpg.table)
  result <- cpg.table[ref2, 
                    on = .(chr = chromosome_name, 
                           Pos19 >= start_position, 
                           Pos19 <= end_position),
                    .(gene, chr, Pos19 = x.Pos19, ensembl_gene_id = i.ensembl_gene_id)]
  # merge result and data by CpG
  data = data[,c("gene","chr","p_alt")]
  final_res <- merge(result, data, by="gene")
  # ACAT
  acat_results <- final_res[, .(p_acat = ACAT(p_alt)), by = ensembl_gene_id]
  acat_results$p_bf = p.adjust(acat_results$p_acat, method = "bonferroni")
  threshold = 0.05/length(acat_results$p_bf)
  cat(basename(file),sum(acat_results$p_bf < 0.05), "threshold:", threshold, "\n")
  all_res <- rbind(all_res,acat_results[acat_results$p_bf < 0.05,])
  saveRDS(acat_results,file = paste0("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/ACAT_450K/",basename(file)))
}

res.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/ACAT_450K/"
res.files <- list.files(res.dir, full.names = TRUE)
all_res = data.frame()
for (file in res.files){
  data = readRDS(file)
  basename = basename(file)
  basename = gsub(".RDS","",basename)
  data$trait = basename
  all_res = rbind(all_res,data[data$p_bf < 0.05,])
}

saveRDS(all_res,file = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA_associations_ACAT_450K.RDS")


##### 900K
res.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/MVP/"
combined_pos38 <- hg38_annot[, c("CpG", "chr", "pos38")]
combined_pos38_900k <- combined_pos38 %>% inner_join(epic850k, by = c("chr" = "CHR", "pos38" = "MAPINFO"))
combined_pos38_900k <- combined_pos38_900k[!duplicated(combined_pos38_900k$CpG),]
combined_pos19_900k <- combined_pos38_900k %>% inner_join(combined_pos19, by = c("chr" = "chr", "CpG" = "gene"))
colnames(combined_pos19_900k) <- c("gene","chr","pos38", "Pos19")
res.files <- list.files(res.dir, full.names = TRUE)
all_res <- data.frame()
for (file in res.files){
  data <- readRDS(file)
  # get gene list
  cpg.table <- combined_pos19_900k[combined_pos19_900k$gene %in% data$gene,]
  cpg.table$chr = as.numeric(gsub("chr","",cpg.table$chr))
  setDT(cpg.table)
  result <- cpg.table[ref2, 
                    on = .(chr = chromosome_name, 
                           Pos19 >= start_position, 
                           Pos19 <= end_position),
                    .(gene, chr, Pos19 = x.Pos19, ensembl_gene_id = i.ensembl_gene_id)]
  # merge result and data by CpG
  data = data[,c("gene","chr","p_alt")]
  final_res <- merge(result, data, by="gene")
  # ACAT
  acat_results <- final_res[, .(p_acat = ACAT(p_alt)), by = ensembl_gene_id]
  acat_results$p_bf = p.adjust(acat_results$p_acat, method = "bonferroni")
  threshold = 0.05/length(acat_results$p_bf)
  cat(basename(file),sum(acat_results$p_bf < 0.05), "threshold:", threshold, "\n")
  all_res <- rbind(all_res,acat_results[acat_results$p_bf < 0.05,])
  saveRDS(acat_results,file = paste0("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/ACAT_900K/",basename(file)))
}

res.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/ACAT_900K/"
res.files <- list.files(res.dir, full.names = TRUE)
all_res = data.frame()
for (file in res.files){
  data = readRDS(file)
  basename = basename(file)
  basename = gsub(".RDS","",basename)
  data$trait = basename
  all_res = rbind(all_res,data[data$p_bf < 0.05,])
}

saveRDS(all_res,file = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA_associations_ACAT_900K.RDS")