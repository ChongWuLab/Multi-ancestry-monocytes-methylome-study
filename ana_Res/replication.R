library(data.table)
# setwd("/rsrch5/home/biostatistics/wzhang24/data/godmc_assoc/")

# data = fread("assoc_meta_all_chunk_aj.csv")
# data = data[data$cistrans == TRUE & data$clumped == TRUE, c("cpg", "snp", "pval", "beta_a1")] # 26220
# data = data[data$pval < 1e-8, ] # 24,994

# #data[, c("snp_chr", "snp_pos37") := tstrsplit(snp, ":", keep = c(1, 2))]
# #data[, snp_pos37 := as.numeric(snp_pos37)]

# HM450 <- read.csv("/rsrch5/home/biostatistics/wzhang24/data/WGBS/humanmethylation450_15017482_v1-2.csv")
# HM450 = HM450[, c("Name", "CHR", "MAPINFO")]
# HM450$CHR = paste0("chr", HM450$CHR)

# # Perform an inner join
# merged_data <- merge(
#   x = data,
#   y = HM450,
#   by.x = "cpg",
#   by.y = "Name",
#   all.x = FALSE, 
#   all.y = FALSE
# )

# dim(merged_data) # 6678374

# setnames(merged_data, c("cpg", "snp", "pval", "CHR", "cpg_pos37"))

# read meQTL from our analysis
meQTL = readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/final_cis_fdr_0.01_CAU.RDS")
#meQTL = meQTL[meQTL$pvalue < 1e-8, ]
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# Suppose your SNP list is:
snp_list <- unique(meQTL$snps)

# Use snpsById() to retrieve positions
# This returns a GRanges object
gr_snps <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, snp_list, ifnotfound = "drop")

# Convert GRanges to a data.frame
df_snps <- as.data.frame(gr_snps)

head(df_snps)

df_snps_dt <- as.data.table(df_snps)

# Rename columns to match your meQTL 'snps' column
setnames(df_snps_dt,
  old = c("RefSNP_id", "seqnames", "pos"),
  new = c("snps", "snp_chr37", "snp_pos37")
)

# Merge
meQTL_anno <- merge(meQTL, df_snps_dt, by = "snps", all.x = TRUE)
head(meQTL_anno)

# read cpg_annot
hg19_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19/combined_data.RDS")
meQTL_final = merge(meQTL_anno, hg19_annot, by.x = "gene", by.y = "CpG")
saveRDS(meQTL_final, "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/final_cis_fdr_0.01_CAU_annot.RDS")