library(data.table)
library(dplyr)
# Read the data
res.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU-association-MVP-bf.RDS")
res.A = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA-association-MVP-bf.RDS")

res.A <- res.A$res.table
res.C <- res.C$res.table

res.A <- res.A[!duplicated(res.A),]
res.C <- res.C[!duplicated(res.C),]

res.A <- res.A[res.A$n_used_snp > 10, ] # 2998
res.C <- res.C[res.C$n_used_snp > 10, ] # 31336

Fingene_id <- fread("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/dbGAP_Upload_Tracking_JUNE2023.csv", skip = 1, data.table = FALSE)
dim(Fingene_id)
Fingene_id <- Fingene_id[,c("Trait", "Ethnicity", "Association Filename","Broad Category")]
Fingene_id$`Association Filename` <- gsub(".txt.gz","",Fingene_id$`Association Filename`)
Fingene_id$`Association Filename` <- gsub(".gz","",Fingene_id$`Association Filename`)

# Merge trait info with res tables
res.A <- merge(res.A, Fingene_id, by.x = "trait", by.y = "Association Filename") # 2998
res.C <- merge(res.C, Fingene_id, by.x = "trait", by.y = "Association Filename") # 31264

# Extract T2D
res.A.T2D <- res.A[res.A$Trait == "T2D",] 
length(unique(res.A.T2D$gene)) # 21
res.C.T2D <- res.C[res.C$Trait == "T2D",] 
length(unique(res.C.T2D$gene)) # 4552

# Get the overlap of genes between AFA and CAU
overlap_cpgs <- intersect(unique(res.A.T2D$gene), unique(res.C.T2D$gene))
length(overlap_cpgs) # 13

# MR 
MR.data.CAU = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/CAU/CAU-MR-MVP-bonferroni.RDS")
MR.data.AFA = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/AFA-MR-MVP-bonferroni.RDS")

MR.data.AFA = MR.data.AFA[,c("file","CpG","MR.p.bf")] #1804
MR.data.AFA = inner_join(MR.data.AFA, Fingene_id, by = c("file" = "Association Filename"))
MR.data.CAU = MR.data.CAU[,c("file","CpG","MR.p.bf")] #1804
MR.data.CAU = inner_join(MR.data.CAU, Fingene_id, by = c("file" = "Association Filename"))

MR.data.AFA.T2D = MR.data.AFA[MR.data.AFA$Trait == "T2D",]
MR.data.CAU.T2D = MR.data.CAU[MR.data.CAU$Trait == "T2D",]

MR.cpg.AFA = unique(MR.data.AFA.T2D$CpG)
MR.cpg.CAU = unique(MR.data.CAU.T2D$CpG)

# Get the number of CpGs in res.A.T2D and res.C.T2D that are in MR.cpg.AFA and MR.cpg.CAU
length(intersect(MR.cpg.AFA, res.A.T2D$gene)) # 15
length(intersect(MR.cpg.CAU, res.C.T2D$gene)) # 3351

# get the number of CpGs in overlap_cpgs that are in either MR.cpg.AFA or MR.cpg.CAU
length(intersect(overlap_cpgs, MR.cpg.AFA)) # 11
length(intersect(overlap_cpgs, MR.cpg.CAU)) # 11

# Coloc
coloc.data.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/08a-coloc/CAU-coloc.RDS")
coloc.data.A = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/07-postPWAS/AFA/AFA-coloc.RDS")

coloc.data.C = coloc.data.C[coloc.data.C$coloc.abf.H4 > 0.8,]
coloc.data.A = coloc.data.A[coloc.data.A$coloc.abf.H4 > 0.8,]

coloc.data.C = inner_join(coloc.data.C, Fingene_id, by = c("file" = "Association Filename"))
coloc.data.A = inner_join(coloc.data.A, Fingene_id, by = c("file" = "Association Filename"))

coloc.data.C.T2D = coloc.data.C[coloc.data.C$Trait == "T2D",]
coloc.data.A.T2D = coloc.data.A[coloc.data.A$Trait == "T2D",]


# Get the number of CpGs in overlap_cpgs that are in either coloc.data.A.T2D or coloc.data.C.T2D
length(intersect(res.A.T2D$gene, coloc.data.A.T2D$CpG)) # 7
length(intersect(res.C.T2D$gene, coloc.data.C.T2D$CpG)) # 1353



# Discover independent genetic risk factors


# SNP CHROMOSOME POSITION A2 A1    MAF   BETA    SE      P SAMPLESIZE
# 1 rs576404767          1   544584  T  C 0.0024 -0.130 0.120 0.2787     826006
# 2 rs538567606          1   565196  C  T 0.0200  0.062 0.130 0.6334     826006
# 3 rs554127336          1   565469  T  C 0.0018 -0.035 0.160 0.8268     826006
# 4 rs544876160          1   565470  A  G 0.0011 -0.079 0.160 0.6215     826006
# 5 rs565235853          1   567006  T  G 0.0027  0.082 0.100 0.4122     826006
# 6   rs9285835          1   569004  C  T 0.0200 -0.013 0.067 0.8462     826006
#           Z
# 1 -1.0833333
# 2  0.4769231
# 3 -0.2187500
# 4 -0.4937500
# 5  0.8200000
# 6 -0.1940299


# Load required libraries if not already loaded
# Load required libraries
library(data.table)
library(dplyr)

# Define GWAS directory and files
gwassum.dir <- "/rsrch5/scratch/biostatistics/wzhang24/GWAS/sumstats/MVP/processed/"
file.names <- c("MVP.T2D.EUR.MAF0.001.combined.dbGaP.sumstats", 
                "MVP.T2D.EUR.MAF001.dbGaP.checked.sumstats", 
                "T2D_MVP_only_EUR_dbGaP.sumstats")

# Define the GWAS significance threshold
gwas_significance_threshold <- 5e-8

# Combine significant SNPs from all GWAS files
all_significant_snps <- data.frame()

# Process each GWAS file
for (file.name in file.names) {
  cat("Processing GWAS file:", file.name, "\n")
  
  # Read the GWAS data
  gwas_data <- fread(paste0(gwassum.dir, file.name)) %>% as.data.frame()
  
  # Filter for significant SNPs only
  significant_snps <- gwas_data %>% 
    filter(P < gwas_significance_threshold) %>%
    select(SNP, CHROMOSOME, POSITION, P)
  
  # Add to the combined data
  all_significant_snps <- bind_rows(all_significant_snps, significant_snps)
}

# Remove duplicates
unique_significant_snps <- all_significant_snps %>%
  distinct(SNP, CHROMOSOME, POSITION, .keep_all = TRUE)

cat("Combined", nrow(unique_significant_snps), "unique significant SNPs from all files\n")

# Get unique CpG sites for T2D from res.C.T2D (European ancestry)
# First, filter for T2D and significant sites
eur_t2d_cpgs <- res.C.T2D %>% 
  filter(Trait == "T2D") %>%
  filter(p_bf < 0.05) %>%
  distinct(gene, chr, p0, p1, .keep_all = TRUE)  # Get unique CpG sites

cat("Found", nrow(eur_t2d_cpgs), "unique significant T2D-associated CpG sites in European ancestry\n")

# Get unique CpG sites for T2D from res.A.T2D (African ancestry)
afr_t2d_cpgs <- res.A.T2D %>% 
  filter(Trait == "T2D") %>%
  filter(p_bf < 0.05) %>%
  distinct(gene, chr, p0, p1, .keep_all = TRUE)  # Get unique CpG sites

cat("Found", nrow(afr_t2d_cpgs), "unique significant T2D-associated CpG sites in African ancestry\n")

# Initialize counters for European ancestry analysis
eur_not_near_gwas <- 0
eur_near_gwas <- 0

# Check each European CpG site
for (i in 1:nrow(eur_t2d_cpgs)) {
  cpg_chr <- eur_t2d_cpgs$chr[i]
  cpg_pos <- mean(c(as.numeric(eur_t2d_cpgs$p0[i]), as.numeric(eur_t2d_cpgs$p1[i])))
  
  # Get significant SNPs on the same chromosome
  chr_snps <- unique_significant_snps %>% filter(CHROMOSOME == cpg_chr)
  
  if(nrow(chr_snps) == 0) {
    # No significant SNPs on this chromosome
    eur_not_near_gwas <- eur_not_near_gwas + 1
    next
  }
  
  # Check if any significant SNP is within 500kb
  distances <- abs(as.numeric(chr_snps$POSITION) - cpg_pos)
  if(any(distances <= 500000)) {
    eur_near_gwas <- eur_near_gwas + 1
  } else {
    eur_not_near_gwas <- eur_not_near_gwas + 1
  }
}

# Calculate percentage for European ancestry
eur_total <- eur_not_near_gwas + eur_near_gwas
eur_percentage <- (eur_not_near_gwas / eur_total) * 100 #24.25

# Initialize counters for African ancestry analysis
afr_not_near_gwas <- 0
afr_near_gwas <- 0

# Check each African CpG site
if(nrow(afr_t2d_cpgs) > 0) {
  for (i in 1:nrow(afr_t2d_cpgs)) {
    cpg_chr <- afr_t2d_cpgs$chr[i]
    cpg_pos <- mean(c(as.numeric(afr_t2d_cpgs$p0[i]), as.numeric(afr_t2d_cpgs$p1[i])))
    
    # Get significant SNPs on the same chromosome
    chr_snps <- unique_significant_snps %>% filter(CHROMOSOME == cpg_chr)
    
    if(nrow(chr_snps) == 0) {
      # No significant SNPs on this chromosome
      afr_not_near_gwas <- afr_not_near_gwas + 1
      next
    }
    
    # Check if any significant SNP is within 500kb
    distances <- abs(as.numeric(chr_snps$POSITION) - cpg_pos)
    if(any(distances <= 500000)) {
      afr_near_gwas <- afr_near_gwas + 1
    } else {
      afr_not_near_gwas <- afr_not_near_gwas + 1
    }
  }
}

# Calculate percentage for African ancestry
afr_total <- afr_not_near_gwas + afr_near_gwas
afr_percentage <- (afr_not_near_gwas / afr_total) * 100 # 24.25

# Print results
cat("\nResults for European Ancestry (EUR):\n")
cat(sprintf("%.2f%% (%d/%d) of significant T2D-associated CpG sites in EA were located in genomic regions with no significant GWAS signals\n", 
          eur_percentage, eur_not_near_gwas, eur_total))

if(nrow(afr_t2d_cpgs) > 0) {
  cat("\nResults for African Ancestry:\n")
  cat(sprintf("%.2f%% (%d/%d) of significant T2D-associated CpG sites in AA were located in genomic regions with no significant GWAS signals\n", 
            afr_percentage, afr_not_near_gwas, afr_total))
}






