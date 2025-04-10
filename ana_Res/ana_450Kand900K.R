library(data.table)
library(dplyr)
# how many models are not in 450K and 900K
# read 450k and 900k
HM450 <- read.csv("/rsrch5/home/biostatistics/wzhang24/data/WGBS/humanmethylation450_15017482_v1-2.csv")
epic850k <- fread("/rsrch5/home/biostatistics/wzhang24/data/WGBS/EPIC-8v2-0_A1.csv", skip=7)
HM450 <- HM450[, c("CHR", "MAPINFO")]
HM450$CHR <- paste0("chr", HM450$CHR)
epic850k <- epic850k[, c("CHR", "MAPINFO")]

hg19_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19/combined_data.RDS")
hg38_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg38_annotation.RDS")
hg38_annot <- hg38_annot[, c("CpG", "pos38")]

# read the table of models
models_AFA = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/04-modelling/AFA/all_res.RDS")
models_CAU = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/04-modelling/CAU/all_res.RDS")

models_AFA = left_join(models_AFA, hg19_annot, by="CpG")
models_AFA = models_AFA[, c("CpG", "Chr", "R2", "p.0", "Pos19")]
colnames(models_AFA) = c("CpG", "Chr", "R2", "pos38", "Pos19")
models_AFA$pos38 = as.numeric(models_AFA$pos38)
models_AFA$Pos19 = as.numeric(models_AFA$Pos19)
models_CAU = left_join(models_CAU, hg19_annot, by="CpG")
models_CAU = models_CAU[, c("CpG", "Chr", "R2", "p.0", "Pos19")]
colnames(models_CAU) = c("CpG", "Chr", "R2", "pos38", "Pos19")
models_CAU$pos38 = as.numeric(models_CAU$pos38)
models_CAU$Pos19 = as.numeric(models_CAU$Pos19)

models_AFA_hm450 <- models_AFA %>% inner_join(HM450, by = c("Chr" = "CHR", "Pos19" = "MAPINFO"))
models_AFA_epic850k <- models_AFA %>% inner_join(epic850k, by = c("Chr" = "CHR", "pos38" = "MAPINFO"))
dim(models_AFA_hm450) # 32870
dim(models_AFA_epic850k) # 59275

models_CAU_hm450 <- models_CAU %>% inner_join(HM450, by = c("Chr" = "CHR", "Pos19" = "MAPINFO"))
models_CAU_epic850k <- models_CAU %>% inner_join(epic850k, by = c("Chr" = "CHR", "pos38" = "MAPINFO"))
dim(models_CAU_hm450) # 27648
dim(models_CAU_epic850k) # 50339

median(models_AFA_hm450$R2) # 0.04
mean(models_AFA_hm450$R2) # 0.11
median(models_AFA_epic850k$R2) # 0.04
mean(models_AFA_epic850k$R2) # 0.10
median(models_CAU_hm450$R2) # 0.04
mean(models_CAU_hm450$R2) # 0.10
median(models_CAU_epic850k$R2) # 0.04

median(models_AFA$R2) # 0.11
mean(models_AFA$R2) # 0.19
sd(models_AFA$R2) # 0.18
median(models_CAU$R2) # 0.12
mean(models_CAU$R2) # 0.17
sd(models_CAU$R2) # 0.16

# conduct hypergeometric test to test if the mean of R2 is different between models_AFA and models_AFA_hm450
# hypergeometric test
# H0: mean(models_AFA$R2) = mean(models_AFA_hm450$R2)
# H1: mean(models_AFA$R2) != mean(models_AFA_hm450$R2)
t.test(models_AFA$R2, models_AFA_hm450$R2, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
t.test(models_AFA$R2, models_AFA_epic850k$R2, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)


### How many association are not in 450k and 900ks
# read tables of associations
asso_AFA = readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/associations/AFA-association-MVP-bf.RDS")
asso_CAU = readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/associations/CAU-association-MVP-bf.RDS")


dim(asso_AFA)
dim(asso_CAU)

asso_AFA = asso_AFA[, c("gene", "p0")]
asso_CAU = asso_CAU[, c("gene", "p0")]
colnames(asso_AFA) = c("CpG", "Pos38")
colnames(asso_CAU) = c("CpG", "Pos38")

asso_AFA = left_join(asso_AFA, hg19_annot, by="CpG")
asso_CAU = left_join(asso_CAU, hg19_annot, by="CpG")

asso_AFA$Pos38 = as.numeric(asso_AFA$Pos38)
asso_AFA_hm450 = asso_AFA %>% inner_join(HM450, by = c("Chr" = "CHR", "Pos19" = "MAPINFO")) # 54
# how many are not in 450k
nrow(asso_AFA) - nrow(asso_AFA_hm450) # 2885
asso_AFA_epic850k = asso_AFA %>% inner_join(epic850k, by = c("Chr" = "CHR", "Pos38" = "MAPINFO")) # 72
nrow(asso_AFA) - nrow(asso_AFA_epic850k) # 2872


asso_CAU$Pos38 = as.numeric(asso_CAU$Pos38)
asso_CAU_hm450 = asso_CAU %>% inner_join(HM450, by = c("Chr" = "CHR", "Pos19" = "MAPINFO")) # 480
# how many are not in 450k
nrow(asso_CAU) - nrow(asso_CAU_hm450) # 30758
asso_CAU_epic850k = asso_CAU %>% inner_join(epic850k, by = c("Chr" = "CHR", "Pos38" = "MAPINFO")) # 942
nrow(asso_CAU) - nrow(asso_CAU_epic850k) # 30303

# For T2D, how many are in 450k and 900k
asso_AFA = readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/associations/AFA-association-MVP-bf.RDS")
asso_CAU = readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/associations/CAU-association-MVP-bf.RDS")

asso_CAU= asso_CAU[asso_CAU$Trait == "T2D",]
asso_AFA = asso_AFA[asso_AFA$Trait == "T2D",]

asso_AFA = asso_AFA[, c("gene", "p0")]
asso_CAU = asso_CAU[, c("gene", "p0")]
colnames(asso_AFA) = c("CpG", "Pos38")
colnames(asso_CAU) = c("CpG", "Pos38")

asso_AFA = left_join(asso_AFA, hg19_annot, by="CpG")
asso_CAU = left_join(asso_CAU, hg19_annot, by="CpG")

asso_AFA$Pos38 = as.numeric(asso_AFA$Pos38)
asso_AFA_hm450 = asso_AFA %>% inner_join(HM450, by = c("Chr" = "CHR", "Pos19" = "MAPINFO")) # 54
length(unique(asso_AFA_hm450$CpG)) # 54
# how many are not in 450k
nrow(asso_AFA) - nrow(asso_AFA_hm450) # 2885
asso_AFA_epic850k = asso_AFA %>% inner_join(epic850k, by = c("Chr" = "CHR", "Pos38" = "MAPINFO")) # 72
length(unique(asso_AFA_epic850k$CpG)) # 0
nrow(asso_AFA) - nrow(asso_AFA_epic850k) # 2872


asso_CAU$Pos38 = as.numeric(asso_CAU$Pos38)
asso_CAU_hm450 = asso_CAU %>% inner_join(HM450, by = c("Chr" = "CHR", "Pos19" = "MAPINFO")) # 59
length(unique(asso_CAU_hm450$CpG)) # 59
# how many are not in 450k
nrow(asso_CAU) - nrow(asso_CAU_hm450) # 30758
asso_CAU_epic850k = asso_CAU %>% inner_join(epic850k, by = c("Chr" = "CHR", "Pos38" = "MAPINFO")) # 942
length(unique(asso_CAU_epic850k$CpG)) # 942
nrow(asso_CAU) - nrow(asso_CAU_epic850k) # 30303
















































