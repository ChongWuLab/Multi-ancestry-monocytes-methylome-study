##########################################################################
##### gene-level analysis##################################################
##########################################################################

library(dplyr)
library(data.table)
library(ggplot2)
#library(grid)
res.A = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA_associations_ACAT.RDS")
res.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU_associations_ACAT.RDS")
save.dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/"

# Get traits information
Fingene_id <- fread("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/dbGAP_Upload_Tracking_JUNE2023.csv", skip = 1, data.table = FALSE)
dim(Fingene_id)
Fingene_id <- Fingene_id[,c("Trait", "Ethnicity", "Association Filename","Broad Category")]
Fingene_id$`Association Filename` <- gsub(".txt.gz","",Fingene_id$`Association Filename`)
Fingene_id$`Association Filename` <- gsub(".gz","",Fingene_id$`Association Filename`)

# Merge trait info with res tables
res.A <- merge(res.A, Fingene_id, by.x = "trait", by.y = "Association Filename") # 1109
res.C <- merge(res.C, Fingene_id, by.x = "trait", by.y = "Association Filename") # 8670
length(unique(res.A$Trait)) #17
length(unique(res.C$Trait)) #30
length(unique(res.A$`Broad Category`)) #7
length(unique(res.C$`Broad Category`)) #7

# map gene to ensembl gene id
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# AA
ensembl_id <- unique(res.A$ensembl_gene_id)
gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensembl_id, mart = ensembl)
geno_info <- gene_info[!duplicated(gene_info$ensembl_gene_id),]
res.A <- left_join(res.A, gene_info, by = "ensembl_gene_id")
write.csv(res.A, "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/AFA-gene-MVP-bf.csv", row.names = FALSE)

# EA
ensembl_id <- unique(res.C$ensembl_gene_id)
gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensembl_id, mart = ensembl)
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id),]
res.C <- left_join(res.C, gene_info, by = "ensembl_gene_id")
write.csv(res.C, "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/CAU-gene-MVP-bf.csv", row.names = FALSE)

# combined res
combined_res <- merge(res.A,res.C,by = c("Trait","ensembl_gene_id", "Broad Category"), suffixes = c(".A",".C")) # 747 associations
length(unique(combined_res$Trait)) #14
length(unique(combined_res$`Broad Category`)) #6
write.csv(combined_res, "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/combined-gene-MVP-bf.csv", row.names = FALSE)

# associations in AA but not in EA
res.A.only <- anti_join(res.A, combined_res, by = c("Trait", "ensembl_gene_id", "Broad Category"))
res.C.only <- anti_join(res.C, combined_res, by = c("Trait", "ensembl_gene_id", "Broad Category"))
write.csv(res.A.only, "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/AFA-only-gene-MVP-bf.csv", row.names = FALSE)
write.csv(res.C.only, "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/CAU-only-gene-MVP-bf.csv", row.names = FALSE)

# Start from reading the results if we have them
combined_res <- read.csv("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/combined-gene-MVP-bf.csv")
res.A <- read.csv("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/AFA-gene-MVP-bf.csv")
res.C <- read.csv("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/CAU-gene-MVP-bf.csv")
# analysis T2D
T2D_res <- combined_res[combined_res$Trait == "T2D",]
length(unique(T2D_res$ensembl_gene_id))

# "ENSG00000148737" "ENSG00000029534" "ENSG00000151532" "ENSG00000140718"
# "ENSG00000165066" "ENSG00000180176" "ENSG00000140382" "ENSG00000106633"
# "ENSG00000149948" "ENSG00000160360" "ENSG00000173517" "ENSG00000124721"
# "ENSG00000153814"
# sort T2D_res by p_bf
T2D_res <- T2D_res[order(T2D_res$p_bf.A),]
unique(T2D_res$ensembl_gene_id)

# TCF7L2, ANK1, VTI1A (no), FTO, 
# NKX6-3, TH, HMG20A, GCK
# HMGA2, GPSM1, PEAK1, DNAH8
# JAZF1

# find genes that are specific in AA but not in EA
T2D_res_A <- res.A[res.A$Trait == "T2D",]
T2D_res_C <- res.C[res.C$Trait == "T2D",]
# find genes that are specific in AA but not in EA
T2D_gene_AA_only = unique(T2D_res_A$ensembl_gene_id[!(T2D_res_A$ensembl_gene_id %in% T2D_res_C$ensembl_gene_id)])
# "ENSG00000106631" "ENSG00000182379" "ENSG00000166006" "ENSG00000255730"
# "ENSG00000248098" "ENSG00000077348" "ENSG00000164512"
T2D_res_A[T2D_res_A$ensembl_gene_id %in% T2D_gene_AA_only,]

# Analysis of AUD
AUD_res <- combined_res[combined_res$Trait == "AUD",]
# ENSG00000138813 C4orf17
# ENSG00000138823 MTTP
# no result

# analysis of CIHS
CIHS_res <- combined_res[combined_res$Trait == "CIHS",]
# ENSG00000130203 APOE
# ENSG00000130204 TOMM40

# height
height_res <- combined_res[combined_res$Trait == "Height",]
height_res <- height_res[order(height_res$p_bf.A),]
unique(height_res$ensembl_gene_id)
length(unique(height_res$ensembl_gene_id))
# "ENSG00000157766" "ENSG00000112033" "ENSG00000065029" "ENSG00000023892"
# ACAN, PPARD, ZNF76, DEF6

# BMI
BMI_res <- combined_res[combined_res$Trait == "BMI",]
BMI_res <- BMI_res[order(BMI_res$p_bf.A),]
unique(BMI_res$ensembl_gene_id)
length(unique(BMI_res$ensembl_gene_id))
# "ENSG00000140718" "ENSG00000188322" "ENSG00000196296" "ENSG00000198156"
# "ENSG00000168488" "ENSG00000178188" "ENSG00000196993" "ENSG00000196502"
# "ENSG00000177548" "ENSG00000197165"

# FTO, SBK1, ATP2A1, NPIPB6
# ATXN2L, SH2B1, 

# HDL
HDL_res <- combined_res[combined_res$Trait == "HDL",] # 66
HDL_res <- HDL_res[order(HDL_res$p_bf.EUR),]

##### Analysis the results for 450k and 900k
# 450k
res.A.450k = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA_associations_ACAT_450K.RDS")
res.C.450k = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU_associations_ACAT_450K.RDS")
save.dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/"

# Get traits information
Fingene_id <- fread("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/dbGAP_Upload_Tracking_JUNE2023.csv", skip = 1, data.table = FALSE)
dim(Fingene_id)
Fingene_id <- Fingene_id[,c("Trait", "Ethnicity", "Association Filename","Broad Category")]
Fingene_id$`Association Filename` <- gsub(".txt.gz","",Fingene_id$`Association Filename`)
Fingene_id$`Association Filename` <- gsub(".gz","",Fingene_id$`Association Filename`)

# Merge trait info with res tables
res.A.450k <- merge(res.A.450k, Fingene_id, by.x = "trait", by.y = "Association Filename") # 137
res.C.450k <- merge(res.C.450k, Fingene_id, by.x = "trait", by.y = "Association Filename") # 1359
res.A.450k <- res.A.450k %>% distinct(trait, ensembl_gene_id, .keep_all = TRUE) # 136
res.C.450k <- res.C.450k %>% distinct(trait, ensembl_gene_id, .keep_all = TRUE) # 1359
res.C <- res.C %>% distinct(trait, ensembl_gene_id, .keep_all = TRUE) # 8647
res.A <- res.A %>% distinct(trait, ensembl_gene_id, .keep_all = TRUE) # 1109
length(unique(res.A.450k$Trait)) #13
length(unique(res.C.450k$Trait)) #31
length(unique(res.A.450k$`Broad Category`)) #6
length(unique(res.C.450k$`Broad Category`)) #7

# overlappings with WGBS
overlap_C <- inner_join(res.C, res.C.450k, by = c("trait", "ensembl_gene_id", "Broad Category")) # 1105
overlap_A <- inner_join(res.A, res.A.450k, by = c("trait", "ensembl_gene_id", "Broad Category")) # 102

# 900k
res.A.900k = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA_associations_ACAT_900K.RDS")
res.C.900k = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU_associations_ACAT_900K.RDS")

res.C.900k <- merge(res.C.900k, Fingene_id, by.x = "trait", by.y = "Association Filename")
res.C.900k <- res.C.900k %>% distinct(trait, ensembl_gene_id, .keep_all = TRUE) # 1929
overlap_C_900k <- inner_join(res.C, res.C.900k, by = c("trait", "ensembl_gene_id", "Broad Category")) # 1644

res.A.900k <- merge(res.A.900k, Fingene_id, by.x = "trait", by.y = "Association Filename")
res.A.900k <- res.A.900k %>% distinct(trait, ensembl_gene_id, .keep_all = TRUE) # 200
overlap_A_900k <- inner_join(res.A, res.A.900k, by = c("trait", "ensembl_gene_id", "Broad Category")) # 153

# combined res
combined_res_900k <- merge(res.A.900k,res.C.900k,by = c("Trait","ensembl_gene_id", "Broad Category"), suffixes = c(".A",".C"))
combined_res_900k <- combined_res_900k %>% distinct(Trait, ensembl_gene_id, .keep_all = TRUE) # 55
# T2D
T2D_res_900k <- combined_res_900k[combined_res_900k$Trait == "T2D",] #3
# KCNQ1, DGKB and TCF7L2

##########################################################################
##### CpG-level analysis##################################################
##########################################################################
library(dplyr)
library(data.table)
library(ggplot2)
res.A = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/05-association/AFA-association-MVP-bf.RDS")
res.C = readRDS("/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/05-association/CAU-association-MVP-bf.RDS")
res.A <- res.A$res.table
res.C <- res.C$res.table
save.dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/"
# remove duplicate
res.A <- res.A[!duplicated(res.A),]
res.C <- res.C[!duplicated(res.C),]

# apply nsnp cutoff
res.A <- res.A[res.A$n_used_snp > 10, ] # 2998
res.C <- res.C[res.C$n_used_snp > 10, ] # 31336

write.csv(res.A, "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/AFA-association-MVP-bf.csv")

# Get traits information
Fingene_id <- fread("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/06-association/dbGAP_Upload_Tracking_JUNE2023.csv", skip = 1, data.table = FALSE)
dim(Fingene_id)
Fingene_id <- Fingene_id[,c("Trait", "Ethnicity", "Association Filename","Broad Category")]
Fingene_id$`Association Filename` <- gsub(".txt.gz","",Fingene_id$`Association Filename`)
Fingene_id$`Association Filename` <- gsub(".gz","",Fingene_id$`Association Filename`)

# Merge trait info with res tables
res.A <- merge(res.A, Fingene_id, by.x = "trait", by.y = "Association Filename") # 2998
res.C <- merge(res.C, Fingene_id, by.x = "trait", by.y = "Association Filename") # 31264

# Merge res.A and res.C using Trait and gene
combined_res <- merge(res.A,res.C,by = c("Trait","gene", "Broad Category"), suffixes = c(".A",".C")) # 468 associations
combined_res$Zscore.AFR <- combined_res$beta_alt.A/combined_res$se_alt.A
combined_res$Zscore.EUR <- combined_res$beta_alt.C/combined_res$se_alt.C

combined_res$Inconsistency <- ifelse((combined_res$Zscore.AFR > 0 & combined_res$Zscore.EUR < 0) | 
                                     (combined_res$Zscore.AFR < 0 & combined_res$Zscore.EUR > 0),
                                     "Inconsistent pairs", as.character(combined_res$`Broad Category`))

combined_res$Inconsistency <- factor(combined_res$Inconsistency, levels = c("Addiction", "Anthropometry", 
                                                                            "CVD", "Lipids", "Metabolic", "Renal",
                                                                             "Inconsistent pairs"))
length(unique(combined_res$gene)) # 433
sum(combined_res$Inconsistency=="Inconsistent pairs") # 22 
unique(combined_res$Trait) # 12
dim(combined_res) # 621
##### Consistency plot
# Color mapping for Broad Categories
colors <- c("Addiction" = "#a6cee3", 
            #"Anthropometry" = "#E68B81", 
            "Anthropometry" = "#1f78b4",
            "CVD" = "hotpink2",
            "Lipids" = "#ff7f00",
            "Metabolic" = "#6a3d9a",
            "Renal" = "#33a02c",
            "Inconsistent pairs" = "grey"
           )

ggplot(combined_res, aes(x = Zscore.AFR, y = Zscore.EUR, color = Inconsistency)) +
  geom_point(size = 1) +
  scale_color_manual(values = colors) +
  labs(x = "Z-score (AA)", y = "Z-score (EA)", title = "") +
  theme_minimal() +
  theme(legend.position = c(0.001, 1.01), # Adjust these values to move the legend inside the plot
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(hjust=1, size = 13),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.height = unit(1.15, "lines"), # Adjust legend key height
        legend.key.width = unit(1.15, "lines"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + # Vertical line at x = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Horizontal line at y = 0
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
ggsave(paste0(save.dir,"3b-scatter.jpg"))


##### manhattan plot

res.A = res.A[, c("p_alt", "Broad Category")]
res.C = res.C[, c("p_alt", "Broad Category")]
res.A$Ethnicity <- "AFR"
res.C$Ethnicity <- "EUR"
res.A$Broad_Category <- factor(res.A$`Broad Category`, levels = c("Addiction", "Anthropometry", "CVD", "Lipids", "Metabolic", "Renal"))
res.C[res.C$`Broad Category` == "Mental Health", "Broad Category"] <- "Mental health"
res.C$Broad_Category <- factor(res.C$`Broad Category`, levels = c("Addiction", "Anthropometry", "CVD", "Lipids", "Metabolic", "Renal", "Mental health"))

# Create a sequence of x-axis positions for the AFR points
x_positions_A <- seq(1, length(res.A$Broad_Category), by = 1)

# Assign x-axis positions to each AFR point in the dataframe
res.A <- res.A %>%
  mutate(xpos = x_positions_A)

# Offset for EUR points
offset <- max(res.A$xpos) + 3000  # The '+ 10' creates a space between AFR and EUR

# Create a sequence of x-axis positions for the EUR points, with an offset

res.C <- res.C %>%
  arrange(Broad_Category)

x_positions_C <- seq(1, length(res.C$Broad_Category), by = 1) + offset


res.C <- res.C %>%
  mutate(xpos = x_positions_C)
res.combined <- bind_rows(res.A, res.C)

x_label_pos = median(res.A$xpos)
x_label_pos2 = median(res.C$xpos)

# Define the colors for the categories, including "Mental Health" for EUR
colors <- c("Addiction" = "#a6cee3", 
            "Anthropometry" = "#1f78b4",
            "CVD" = "hotpink2",
            "Lipids" = "#ff7f00",
            "Metabolic" = "#6a3d9a",
            "Renal" = "#33a02c",
            "Mental health" = "grey")

# Plotting with both AFR and EUR data using xpos as the numeric x-axis
p <- ggplot(res.combined, aes(x = xpos, y = -log10(p_alt), color = Broad_Category)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +  # Increase base font size
  labs(y = "-log10(P-value)", x = "") +
  scale_y_continuous(limits = c(0, 120), breaks = c(10, 30, 50, 70, 90, 110), expand = c(0, 0)) +
  theme(legend.position = c(0.4, 0.95), # Adjust these values to move the legend inside the plot
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),  # Increase font size of x axis title
        axis.title.y = element_text(size = 14),  # Increase font size of y axis title
        axis.text.x = element_blank(),           # Remove x-axis text
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),  # Add y-axis ticks
        axis.text.y = element_text(size = 14),         # Ensure y-axis text is visible and adjusted size
        panel.background = element_rect(fill = "white", colour = NA), # White background
        panel.grid.major = element_blank(),      # Remove major grid lines
        panel.grid.minor = element_blank(),      # Remove minor grid lines
        panel.border = element_blank(),          # Remove panel border
        # add y axis border
        axis.line.y = element_line(color = "black")) + 
  guides(color = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(size = 3))) +
  annotate("text", x = x_label_pos, y = min(-log10(res.combined$p_alt)) - 1.5, label = "AA", size = 5, hjust = 0.5, vjust = 1) +
  annotate("text", x = x_label_pos2, y = min(-log10(res.combined$p_alt)) - 1.5, label = "EA", size = 5, hjust = 0.5, vjust = 1)
  # y limit from 10 to 120
  

ggsave(paste0(save.dir, "3a-manhattan.jpg"), p, width = 11, height = 7, units = "in", dpi = 300)


##### CpG annotation
annotation = readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg38_annotation.RDS")

# get consistent table
consistent_res <- combined_res[combined_res$Inconsistency != "Inconsistent pairs",]
consistent_res <- consistent_res[,c("Trait", "gene", "Broad Category","p_alt.A","p_alt.C")]

# merge
merged_res <- merge(consistent_res, annotation, by.x = "gene", by.y = "CpG")
sum(merged_res$annotation == "open_sea") # 529
sum(merged_res$annotation == "shore") # 30
sum(merged_res$annotation == "shelf") # 24
sum(merged_res$annotation == "island") # 16