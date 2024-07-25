library(data.table)
library(ggplot2)
library(dplyr)
library(BEDMatrix)
#library(scales)
#library(patchwork)
#library(cowplot)
ref_EUR = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.EUR.ALLSNP.QC.CHR"
ref_AFR = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"

res_dir_CAU = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/CAU/summary/"
res_dir_AFA = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/AFA/summary/"
save_dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/"

final_cis_CAU = readRDS(paste0(res_dir_CAU,"final_cis_CAU.RDS"))
final_cis_AFA <- readRDS(paste0(res_dir_AFA, "final_cis_AFA.RDS"))
final_trans_CAU <- readRDS(paste0(res_dir_CAU, "final_trans_CAU.RDS"))
final_trans_AFA <- readRDS(paste0(res_dir_AFA, "final_trans_AFA.RDS"))


# colors to be chosen from
#E76254 red 
#EF8A47 orange
#F7AA58 light orange
#FFD06F yellow
#FFE6B7 light yellow
#AADCE0 light blue
#72BCD5 blue 
#528FAD dark blue
#376795 dark dark blue
#1E466E dark dark dark blue
#93BE6C green


#frq_AFA <- fread("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/01-SNP-res/AFA/rs_1000G_all.frq")
#frq_CAU <- fread("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/01-SNP-res/CAU/rs_1000G_all.frq")
#frq_AFA <- frq_AFA[, c("SNP", "MAF")]
#frq_CAU <- frq_CAU[, c("SNP", "MAF")]

#final_cis_CAU <- left_join(final_cis_CAU, frq_CAU, by = c("SNP" = "SNP"))
#final_cis_AFA <- left_join(final_cis_AFA, frq_AFA, by = c("SNP" = "SNP"))
#final_trans_CAU <- left_join(final_trans_CAU, frq_CAU, by = c("SNP" = "SNP"))
#final_trans_AFA <- left_join(final_trans_AFA, frq_AFA, by = c("SNP" = "SNP"))
# Read the annotation file
#cpg_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg38_annotation.RDS")
# Merge meQTL data with CpG annotations
#final_cis_AFA <- merge(final_cis_AFA, cpg_annot, by = "CpG")
#final_cis_CAU <- merge(final_cis_CAU, cpg_annot, by = "CpG")
#final_trans_AFA <- merge(final_trans_AFA, cpg_annot, by = "CpG")
#final_trans_CAU <- merge(final_trans_CAU, cpg_annot, by = "CpG")
#saveRDS(final_cis_CAU, file = paste0(res_dir_CAU,"final_cis_CAU_160.RDS"))
allele.qc <- function(a1, a2, ref1, ref2) {
    ref <- ref1
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip1 <- flip
    ref <- ref2
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip2 <- flip
    snp <- list()
    snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
    snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}


## remove ambiguous SNPs
#a1 = final_cis_AFA$A1
#a2 = final_cis_AFA$A2
#keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
#final_cis_AFA = final_cis_AFA[keep,]

#a1 = final_cis_CAU$A1
#a2 = final_cis_CAU$A2
#keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
#final_cis_CAU = final_cis_CAU[keep,]

#a1 = final_trans_AFA$A1
#a2 = final_trans_AFA$A2
#keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
#final_trans_AFA = final_trans_AFA[keep,]

#a1 = final_trans_CAU$A1
#a2 = final_trans_CAU$A2
#keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
#final_trans_CAU = final_trans_CAU[keep,]


#saveRDS(final_cis_AFA, file = paste0(res_dir_AFA,"final_cis_AFA.RDS"))
#saveRDS(final_cis_CAU, file = paste0(res_dir_CAU,"final_cis_CAU.RDS"))
#saveRDS(final_trans_AFA, file = paste0(res_dir_AFA, "final_trans_AFA.RDS"))
#saveRDS(final_trans_CAU, file = paste0(res_dir_CAU, "final_trans_CAU.RDS"))

# flip allele



### merge the two data frames and find common cpgs
data1_AFA <- final_cis_AFA[, c("CpG", "SNP", "beta")] #45143297
data1_CAU <- final_cis_CAU[, c("CpG", "SNP", "beta")] #98259931
data2 <- inner_join(data1_AFA, data1_CAU, by = c("CpG" = "CpG", "SNP" = "SNP"), suffix = c(".AFR", ".EUR"))
common_cpgs <- unique(data2$CpG) 
length(common_cpgs) # 880108
#common_cpgs2 <- intersect(final_cis_AFA$CpG, final_cis_CAU$CpG)
final2_cis_CAU <- final_cis_CAU[final_cis_CAU$CpG %in% common_cpgs, ]
final2_cis_AFA <- final_cis_AFA[final_cis_AFA$CpG %in% common_cpgs, ]

set.seed(123)
sampled_df_CAU <- final2_cis_CAU[sample(nrow(final2_cis_CAU), size = 5000), ]
set.seed(123)
sampled_df_AFA <- final2_cis_AFA[sample(nrow(final2_cis_AFA), size = 5000), ]

set.seed(123)
sampled_df_all <- data2[sample(nrow(data2), size = 5000), ]
saveRDS(sampled_df_all, file = paste0(res_dir_CAU, "sampled_data_all.RDS"))

saveRDS(sampled_df_CAU, file = paste0(res_dir_CAU,"sampled_data_CAU.RDS"))
saveRDS(sampled_df_AFA, file = paste0(res_dir_AFA,"sampled_data_AFA.RDS"))
# Or if only want the sampled_data, it can be read directly
sampled_df_CAU <- readRDS(paste0(res_dir_CAU,"sampled_data_CAU.RDS"))
sampled_df_AFA <- readRDS(paste0(res_dir_AFA,"sampled_data_AFA.RDS"))
median(abs(sampled_df_CAU$beta)) # 0.65
median(abs(sampled_df_AFA$beta)) # 0.84
wilcox.test(abs(sampled_df_CAU$beta), abs(sampled_df_AFA$beta)) # p-value = 0.0001

sampled_df_all <- readRDS(paste0("sampled_data_all.RDS"))

# Combine the two data frames and add a column to indicate the group
sampled_df_CAU$group <- "EA"
#remove the first column of sampled_df_CAU
sampled_df_AFA$group <- "AA"
combined_df <- rbind(sampled_df_CAU, sampled_df_AFA)
combined_df$MAF_transformed <- combined_df$MAF * (1 - combined_df$MAF)
combined_df$effect_size <- abs(combined_df$beta)

##### a. MAF(1-MAF) density plot
p <- ggplot(combined_df, aes(x = MAF_transformed, y = effect_size, color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "MAF(1-MAF)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD"))

ggsave(paste0(save_dir,"Fig2a-MAF-beta.jpg"), width = 6, height = 5)

# get the coefficient of the linear model
fit1 = lm(effect_size ~ MAF_transformed, data = combined_df[combined_df$group == "EA",]) # coeff: -2.04 pval: <2e-16
fit2 = lm(effect_size ~ MAF_transformed, data = combined_df[combined_df$group == "AA",])

p <- ggplot(combined_df, aes(x = MAF_transformed, y = -log(pvalue), color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "MAF(1-MAF)", y = "-log10(P-value)") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD"))

ggsave(paste0(save_dir,"a-MAF-pval.jpg"))

##### b. Scatter plot distance 
# Create the plot
# calculate coefficient and p values
fit1 = lm(-log10(pvalue) ~ abs(distance), data = combined_df[combined_df$group == "EA",]) # coeff: 1.81, pval
fit2 = lm(-log10(pvalue) ~ abs(distance), data = combined_df[combined_df$group == "AA",])
p <- ggplot(combined_df, aes(x = distance, y = -log10(pvalue))) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "-log10(P-value)") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), # Adjusted space for x-axis label
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Adjusted space for y-axis label
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"), # Set background color to white
        panel.grid.major = element_line(color = "grey90"), # Set major grid lines to a light grey
        panel.grid.minor = element_line(color = "grey95"), # Set minor grid lines to a very light grey
        axis.line = element_line(color = "black", linewidth = 0.3)) + # Set axis lines to black
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD")) +
  scale_x_continuous(breaks = c(-500000, 500000), labels = c("-500,000", "500,000"))

ggsave(paste0(save_dir,"b-distance-pval.jpg"))

p <- ggplot(combined_df, aes(x = distance, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), # Adjusted space for x-axis label
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Adjusted space for y-axis label
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"), # Set background color to white
        panel.grid.major = element_line(color = "grey90"), # Set major grid lines to a light grey
        panel.grid.minor = element_line(color = "grey95"), # Set minor grid lines to a very light grey
        axis.line = element_line(color = "black", linewidth = 0.3)) + # Set axis lines to black
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD")) +
  scale_x_continuous(breaks = c(-500000, 500000), labels = c("-500,000", "500,000"))

ggsave(paste0(save_dir,"b-distance-beta.jpg"), width = 6, height = 5)

###################################################################
##### draw the plots for overlapping cpgs
###################################################################

##### b. Scatter plot distance 
p <- ggplot(combined_df, aes(x = distance, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), # Adjusted space for x-axis label
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Adjusted space for y-axis label
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"), # Set background color to white
        panel.grid.major = element_line(color = "grey90"), # Set major grid lines to a light grey
        panel.grid.minor = element_line(color = "grey95"), # Set minor grid lines to a very light grey
        axis.line = element_line(color = "black", linewidth = 0.3)) + # Set axis lines to black
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD")) +
  scale_x_continuous(breaks = c(-500000, 500000), labels = c("-500,000", "500,000"))

ggsave(paste0(save_dir,"sup-distance-beta.jpg"))


###########################################################################
# Draw the figure of effect size and MAF
p <- ggplot(combined_df, aes(x = MAF, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  labs(x = "MAF", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.2)) +
  scale_color_manual(values = c("EUR" = "#E76254", "AFR" = "#528FAD"))

ggsave(paste0(save_dir,"c-MAF.jpg"))

# Draw the plot of effect size and sd of methylation data
sd_CAU <- readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/01-DNAm-var/CAU_SD.RDS")
sd_AFA <- readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/01-DNAm-var/AFA_SD.RDS")

sampled_df_CAU <- left_join(sampled_df_CAU, sd_CAU, by = c("CpG" = "CpG"))
sampled_df_AFA <- left_join(sampled_df_AFA, sd_AFA, by = c("CpG" = "CpG"))

combined_df <- rbind(sampled_df_CAU, sampled_df_AFA)
combined_df$effect_size <- abs(combined_df$beta)

p <- ggplot(combined_df, aes(x = SD, y = effect_size, color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "SD of DNAm", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.2)) +
  scale_color_manual(values = c("EUR" = "#E76254", "AFR" = "#528FAD"))

ggsave(paste0(save_dir,"d-SD.jpg"))

#######################################################
##### Different numbers of mQTL based on CpG annotation
#######################################################

# some basic stats
nrow(final_cis_AFA[final_cis_AFA$annotation== "island",])/nrow(final_cis_AFA)
nrow(final_cis_AFA[final_cis_AFA$annotation %in% c("shore", "shelf"),])/nrow(final_cis_AFA)
nrow(final_cis_AFA[final_cis_AFA$annotation== "shelf",])/nrow(final_cis_AFA)
nrow(final_cis_AFA[final_cis_AFA$annotation== "open_sea",])/nrow(final_cis_AFA)

nrow(final_cis_CAU[final_cis_CAU$annotation== "island",])/nrow(final_cis_CAU)
nrow(final_cis_CAU[final_cis_CAU$annotation %in% c("shore", "shelf"),])/nrow(final_cis_CAU)
nrow(final_cis_CAU[final_cis_CAU$annotation== "shelf",])/nrow(final_cis_CAU)
nrow(final_cis_CAU[final_cis_CAU$annotation== "open_sea",])/nrow(final_cis_CAU)

# Among all the CpG sites located in islands, how many of them have at least one cis-meQTL
# read annotation data
cpg_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg38_annotation.RDS")
islands = cpg_annot[cpg_annot$annotation == "island", "CpG"] 
adjacent = cpg_annot[cpg_annot$annotation %in% c("shore", "shelf"), "CpG"] 
open_sea = cpg_annot[cpg_annot$annotation == "open_sea", "CpG"] 
# make them to list
islands = islands$CpG
adjacent = adjacent$CpG
open_sea = open_sea$CpG
# calculate the proportion of CpG sites with at least one cis-meQTL
meqtl_cpgs = unique(final_cis_AFA$CpG)
length(intersect(islands, meqtl_cpgs))/length(islands) # 2.9%
length(intersect(adjacent, meqtl_cpgs))/length(adjacent) # 2.5%
length(intersect(open_sea, meqtl_cpgs))/length(open_sea) # 6.9%

meqtls_cpgs = unique(final_cis_CAU$CpG)
length(intersect(islands, meqtls_cpgs))/length(islands) # 1.8%
length(intersect(adjacent, meqtls_cpgs))/length(adjacent) # 4.2%
length(intersect(open_sea, meqtls_cpgs))/length(open_sea) # 5.6 %


# start plotting
# Count meQTLs per CpG site
count_AFR <- final_cis_AFA %>% group_by(CpG, annotation) %>% summarise(meQTL_count = n())
count_EUR <- final_cis_CAU %>% group_by(CpG, annotation) %>% summarise(meQTL_count = n())

# Calculate proportions for each annotation
#total_AFR <- nrow(final_data_AFR)
#total_EUR <- nrow(final_data_EUR)
#count_AFR <- count_AFR %>% mutate(proportion = meQTL_count / total_AFR)
#count_EUR <- count_EUR %>% mutate(proportion = meQTL_count / total_EUR)



# Adding a population column to each
count_AFR$population <- 'AA'
count_EUR$population <- 'EA'
count_AFR = as.data.frame(count_AFR)
count_EUR = as.data.frame(count_EUR)
mean(count_AFR[count_AFR$annotation=="island", "meQTL_count"])
mean(count_AFR[count_AFR$annotation=="shore", "meQTL_count"])
mean(count_AFR[count_AFR$annotation=="shelf", "meQTL_count"])
mean(count_AFR[count_AFR$annotation=="open_sea", "meQTL_count"])
# Combine the data for both populations
combined_data <- rbind(count_AFR, count_EUR)
combined_data$annotation <- factor(combined_data$annotation,
                                   levels = c( "open_sea", "island", "shore", "shelf"),
                                   labels = c("Open Sea", "Island", "Shore", "Shelf"))

# Create histograms
p <- ggplot(combined_data, aes(x = meQTL_count, fill = annotation)) +
  geom_histogram(binwidth = 1, position = "dodge", boundary = 0) +
  facet_wrap(~ population, nrow = 1, ncol = 2, labeller = labeller(population = c(AA = "AA", EA = "EA"))) +
  scale_x_continuous(limits = c(1, 20), breaks = c(1, 5, 10, 15, 20)) +
  scale_y_continuous(labels = scales::comma) + 
  labs(x = "Number of meQTLs per CpG site", y = "Number of CpG sites", fill = "DNAm site annotation") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_fill_manual(values = c("Open Sea" = "#93BE6C", "Island" = "#E76254", "Shore" = "#F7AA58", "Shelf" = "#72BCD5"))

ggsave(paste0(save_dir,"e-hist.jpg"), width = 10, height = 5)




##### Compare to 850k and 450k
# Read manifest data
manifest <- read.table("/rsrch5/home/biostatistics/wzhang24/data/WGBS/HM450.hg38.manifest.tsv", sep = "\t", header = TRUE)
manifest <- manifest[, c("CpG_chrm", "CpG_beg")]
colnames(manifest) <- c("chr", "pos")

##### Cis + trans together
set.seed(123)
sampled_cis_CAU <- final_cis_CAU[sample(nrow(final_cis_CAU), size = 5000), ]
set.seed(123)
sampled_cis_AFA <- final_cis_AFA[sample(nrow(final_cis_AFA), size = 5000), ]
set.seed(123)
sampled_trans_CAU <- final_trans_CAU[sample(nrow(final_trans_CAU), size = 5000), ]
set.seed(123)
sampled_trans_AFA <- final_trans_AFA[sample(nrow(final_trans_AFA), size = 5000), ]

sampled_cis_CAU = sampled_cis_CAU[, c("CpG", "SNP", "beta", "pvalue", "MAF", "distance")]
sampled_cis_AFA = sampled_cis_AFA[, c("CpG", "SNP", "beta", "pvalue", "MAF", "distance")]
sampled_trans_CAU = sampled_trans_CAU[, c("CpG", "SNP", "beta", "pvalue", "MAF", "distance")]
sampled_trans_AFA = sampled_trans_AFA[, c("CpG", "SNP", "beta", "pvalue", "MAF", "distance")]


# label the data
sampled_cis_CAU$group <- "EUR cis"
sampled_cis_AFA$group <- "AFR cis"
sampled_trans_CAU$group <- "EUR trans"
sampled_trans_AFA$group <- "AFR trans"

combined_df <- rbind(sampled_cis_CAU, sampled_cis_AFA, sampled_trans_CAU, sampled_trans_AFA)
combined_df$MAF_transformed <- combined_df$MAF * (1 - combined_df$MAF)
combined_df$effect_size <- abs(combined_df$beta)

# MAF(1-MAF) plot
p <- ggplot(combined_df, aes(x = MAF_transformed, y = effect_size, color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "MAF(1-MAF)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EUR cis" = "#E76254", "EUR trans" = "#E76254", "AFR cis" = "#528FAD", "AFR trans" = "#528FAD"))
ggsave(paste0(save_dir,"a-MAF-beta-all.jpg"))

# Distance plot
p <- ggplot(combined_df, aes(x = distance, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "free", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EUR cis" = "#E76254", "EUR trans" = "#E76254", "AFR cis" = "#528FAD", "AFR trans" = "#528FAD"))


ggsave(paste0(save_dir, "b-distance-beta-all.jpg"), width = 10, height = 8)





trans_data_check <- combined_df[combined_df$group %in% c("EUR trans", "AFR trans"),]
trans_data_check$distance <- abs(trans_data_check$distance)
summary(trans_data_check$distance)

p <- ggplot(trans_data_check, aes(x = distance, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "Effect size") +
  theme_classic() +
  scale_color_manual(values = c("EUR trans" = "#E76254", "AFR trans" = "#528FAD")) +
  ylim(0, 2.5)

ggsave(paste0(save_dir, "test.jpg"), width = 10, height = 8)




########################################################################
# get sumstats
########################################################################
dim(final_cis_AFA)
dim(final_cis_CAU)
data1_AFA <- final_cis_AFA[, c("CpG", "SNP")] #56112500
data1_CAU <- final_cis_CAU[, c("CpG", "SNP")] #122777173
data2 <- inner_join(final_cis_AFA, final_cis_CAU, by = c("CpG" = "CpG", "SNP" = "SNP"))
common_cpgs <- unique(data2$CpG)
length(unique(data2$CpG)) #861445

# Effect size
final_cis_AFA$effect_size <- abs(final_cis_AFA$beta)
final_cis_CAU$effect_size <- abs(final_cis_CAU$beta)
median(final_cis_AFA$effect_size) #
median(final_cis_CAU$effect_size)


################################################################
# to detect nonexistence or rare cis-meQTLs in other populations
################################################################
# for EUR
# read 1000G for AFR
ref_AFR = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"
SNP_in_1000G_AFR <- c()
for (chr in 1:22){
  reference.bed =  paste0(ref_AFR, chr) %>% BEDMatrix(., simple_names = TRUE) %>% as.matrix()
  variant_present <- reference.bed > 0
  variant_counts <- colSums(variant_present)
  not_rare_snps <- names(variant_counts[variant_counts > 2])
  SNP_in_1000G_AFR <- c(SNP_in_1000G_AFR, not_rare_snps)
}

final_cis_CAU$in_1000G <- ifelse(final_cis_CAU$SNP %in% SNP_in_1000G_AFR, 1, 0)
saveRDS(final_cis_CAU, file = paste0(res_dir_CAU, "final_cis_CAU.RDS"))
1 - sum(final_cis_CAU$in_1000G)/dim(final_cis_CAU)[1] # 0.08 nonexistence or rare



# for AFR
# read 1000G for EUR
ref_EUR = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.EUR.ALLSNP.QC.CHR"
SNP_in_1000G_EUR <- c()
for (chr in 1:22){
  reference.bed =  paste0(ref_EUR, chr) %>% BEDMatrix(., simple_names = TRUE) %>% as.matrix()
  variant_present <- reference.bed > 0
  variant_counts <- colSums(variant_present)
  not_rare_snps <- names(variant_counts[variant_counts > 2])
  SNP_in_1000G_EUR <- c(SNP_in_1000G_EUR, not_rare_snps)
}

final_cis_AFA$in_1000G <- ifelse(final_cis_AFA$SNP %in% SNP_in_1000G_EUR, 1, 0)
saveRDS(final_cis_AFA, file = paste0(res_dir_AFA, "final_cis_AFA.RDS"))
1 - sum(final_cis_AFA$in_1000G)/dim(final_cis_AFA)[1] # 0.3 nonexistence or rare

# plot
data <- data.frame(
  Category = c("AA", "EA", 
               "Sentinel AA", "Sentinel EA"),
  Percentage = c(30, 8, 49, 19)
)

# Create the bar chart
ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", show.legend = FALSE, width = 0.5, color = "black") +  # Adjust width here
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = "Percentage (%) of nonexistent or rare cis-meQTL"
  ) +
  scale_fill_manual(values = c("AA" = "#72BCD5", "Sentinel AA" = "#72BCD5", 
                               "EA" = "#E76254", "Sentinel EA" = "#E76254")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
ggsave(paste0(save_dir,"1d-nonexistence.jpg"), width = 6, height = 6)


# Sentinel meQTL
# Group by CpG, then filter for the row with the minimum p-value in each group
sentinel_meQTL_EUR <- final_cis_CAU %>%
  group_by(CpG) %>%
  filter(pvalue == min(pvalue)) %>%
  ungroup()


sentinel_meQTL_EUR <- as.data.frame(sentinel_meQTL_EUR)
sentinel_meQTL_EUR$in_1000G <- ifelse(sentinel_meQTL_EUR$SNP %in% SNP_in_1000G_AFR, 1, 0)
1 - sum(sentinel_meQTL_EUR$in_1000G)/dim(sentinel_meQTL_EUR)[1] # 0.19

# AFR
sentinel_meQTL_AFR <- final_cis_AFA %>%
  group_by(CpG) %>%
  filter(pvalue == min(pvalue)) %>%
  ungroup()

sentinel_meQTL_AFR <- as.data.frame(sentinel_meQTL_AFR)
sentinel_meQTL_AFR$in_1000G <- ifelse(sentinel_meQTL_AFR$SNP %in% SNP_in_1000G_EUR, 1, 0)
1 - sum(sentinel_meQTL_AFR$in_1000G)/dim(sentinel_meQTL_AFR)[1] # 0.49




########### plot consistent meQTL#############################
# For all meQTL
allele.qc <- function(a1, a2, ref1, ref2) {
    ref <- ref1
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip1 <- flip
    ref <- ref2
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip2 <- flip
    snp <- list()
    snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
    snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}

data1_AFA <- final_cis_AFA[, c("CpG", "SNP", "beta", "A1", "A2")] #45143297
data1_CAU <- final_cis_CAU[, c("CpG", "SNP", "beta", "A1", "A2")] #98259931
data2 <- inner_join(data1_AFA, data1_CAU, by = c("CpG" = "CpG", "SNP" = "SNP"), suffix = c(".AFR", ".EUR"))
# flip alleles

qc.1 <- allele.qc(
        data2$A1.EUR,
        data2$A2.EUR,
        data2$A1.AFR,
        data2$A2.AFR
    )
data2$beta.EUR[qc.1$flip] <- -1 * data2$beta.EUR[qc.1$flip]
data2 = data2[qc.1$keep,]
sampled_df = data2[sample(nrow(data2), size = 5000), ]
saveRDS(sampled_df, file = paste0(save_dir, "sampled_consistent_df.RDS"))

# plot all meQTL
sampled_df = readRDS(paste0(save_dir, "sampled_consistent_df.RDS"))
# corr
cor(sampled_df$beta.AFR, sampled_df$beta.EUR) # 0.96


ggplot(sampled_df, aes(x = beta.AFR, y = beta.EUR)) +
  geom_point(size = 1, color = "#B36A6F") +
  scale_color_manual(values = colors) +
  labs(x = "Effect size (AA)", y = "Effect size (EA)", title = "") +
  theme_minimal() +
  theme(legend.position = c(0.001, 1.01), # Adjust these values to move the legend inside the plot
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(hjust=1, size = 14),
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
ggsave(paste0(save_dir,"1e-consistent_all.jpg"), width = 5, height = 5)

# plot the sentinel cis-meQTLs across AFR and EUR
plot_data_EUR <- sentinel_meQTL_EUR[, c("CpG", "SNP", "beta")]
plot_data_AFR <- sentinel_meQTL_AFR[, c("CpG", "SNP", "beta")]
merged_data <- merge(plot_data_EUR, plot_data_AFR, by = c("SNP", "CpG"), suffixes = c(".EUR", ".AFR"))
set.seed(123)
sampled_data <- merged_data[sample(nrow(merged_data), 5000), ]

ggplot(sampled_df_all, aes(x = beta.AFR, y = beta.EUR)) +
  geom_point() +  # Add points
  theme_minimal() +  # Optional: use a minimal theme
  labs(
    x = "Effect Size in AFR",
    y = "Effect Size in EUR"
  )
ggsave(paste0(save_dir, "sentinel.jpg"))





#################################################################################
####################trans-meQTLs#################################################
#################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
res_dir_CAU = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/CAU/summary/"
res_dir_AFA = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/AFA/summary/"
save_dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/"

#final_cis_CAU = readRDS(paste0(res_dir_CAU,"final_cis_CAU.RDS"))
#final_cis_AFA <- readRDS(paste0(res_dir_AFA, "final_cis_AFA.RDS"))
final_trans_CAU <- readRDS(paste0(res_dir_CAU, "final_trans_CAU.RDS"))
final_trans_AFA <- readRDS(paste0(res_dir_AFA, "final_trans_AFA.RDS"))

# unique CpG
length(unique(final_trans_CAU$CpG)) # 103768
length(unique(final_trans_AFA$CpG)) # 159053
### merge the two data frames and find common cpgs
data1_AFA <- final_trans_AFA[, c("CpG", "SNP", "beta")] #45143297
data1_CAU <- final_trans_CAU[, c("CpG", "SNP", "beta")] #98259931
data2 <- inner_join(data1_AFA, data1_CAU, by = c("CpG" = "CpG", "SNP" = "SNP"), suffix = c(".AFR", ".EUR"))
common_cpgs <- unique(data2$CpG) 
length(common_cpgs)

set.seed(123)
sampled_trans_CAU <- final_trans_CAU[sample(nrow(final_trans_CAU), size = 5000), ]
set.seed(123)
sampled_trans_AFA <- final_trans_AFA[sample(nrow(final_trans_AFA), size = 5000), ]


# label
sampled_trans_CAU$group <- "EA"
sampled_trans_AFA$group <- "AA"

combined_df <- rbind(sampled_trans_CAU, sampled_trans_AFA)
combined_df$MAF_transformed <- combined_df$MAF * (1 - combined_df$MAF)
combined_df$effect_size <- abs(combined_df$beta)

a1 = combined_df$A1
a2 = combined_df$A2
keep <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
combined_df = combined_df[keep,]
##### a. MAF(1-MAF) density plot
combined_df2 = combined_df[combined_df$MAF>0.05, ]
fit1 = lm(effect_size ~ MAF_transformed, data = combined_df2[combined_df2$group == "EA",])
fit2 = lm(effect_size ~ MAF_transformed, data = combined_df2[combined_df2$group == "AA",])

p <- ggplot(combined_df, aes(x = MAF_transformed, y = effect_size, color = group)) +
  geom_point(alpha = 0.5) +
  #geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "MAF(1-MAF)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD"))

ggsave(paste0(save_dir,"4a-MAF-beta-trans.png"), width = 6, height = 5, dpi = 300)

p <- ggplot(combined_df, aes(x = distance, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), # Adjusted space for x-axis label
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Adjusted space for y-axis label
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"), # Set background color to white
        panel.grid.major = element_line(color = "grey90"), # Set major grid lines to a light grey
        panel.grid.minor = element_line(color = "grey95"), # Set minor grid lines to a very light grey
        axis.line = element_line(color = "black", linewidth = 0.3)) + # Set axis lines to black
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#528FAD")) +
  scale_x_continuous(breaks = c(-2e8, 0, 2e8), labels = c("-2e+08", "0", "2e+08"))


ggsave(paste0(save_dir,"4b-distance-beta-trans.png"), width = 6, height = 5, dpi = 300)

########### plot consistent meQTL#############################
# For all meQTL
data1_AFA <- final_trans_AFA[, c("CpG", "SNP", "beta", "A1", "A2")] #45143297
data1_CAU <- final_trans_CAU[, c("CpG", "SNP", "beta", "A1", "A2")] #98259931
data2 <- inner_join(data1_AFA, data1_CAU, by = c("CpG" = "CpG", "SNP" = "SNP"), suffix = c(".AFR", ".EUR"))
# flip alleles

qc.1 <- allele.qc(
        data2$A1.EUR,
        data2$A2.EUR,
        data2$A1.AFR,
        data2$A2.AFR
    )

data2$beta.EUR[qc.1$flip] <- -1 * data2$beta.EUR[qc.1$flip]
keep = qc.1$keep
data2 = data2[keep,]
sampled_df = data2[sample(nrow(data2), size = 5000), ]

# plot all meQTL

ggplot(sampled_df, aes(x = beta.AFR, y = beta.EUR)) +
  geom_point(size = 1, color = "#B36A6F") +
  scale_color_manual(values = colors) +
  labs(x = "Effect size (AA)", y = "Effect size (EA)", title = "") +
  theme_minimal() +
  theme(legend.position = c(0.001, 1.01), # Adjust these values to move the legend inside the plot
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(hjust=1, size = 14),
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
ggsave(paste0(save_dir,"4c-consistent_all_trans.jpg"), width = 5, height = 5, dpi = 300)


#################################################################################
####################Compare AFR and EUR in the same sample size##################
#You can start from here to run the following code, no need to run the previous codes
#################################################################################
ref_EUR = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.EUR.ALLSNP.QC.CHR"
ref_AFR = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.AFR.ALLSNP.QC.CHR"

res_dir_CAU = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/CAU/summary/"
res_dir_AFA = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/AFA/summary/"
save_dir = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/"

final_cis_CAU = readRDS(paste0(res_dir_CAU,"final_cis_CAU_160.RDS"))
final_cis_AFA <- readRDS(paste0(res_dir_AFA, "final_cis_AFA.RDS"))

### merge the two data frames and find common cpgs
data1_AFA <- final_cis_AFA[, c("CpG", "SNP", "beta")] #56112500
data1_CAU <- final_cis_CAU[, c("CpG", "SNP", "beta")] #48086725
data2 <- inner_join(data1_AFA, data1_CAU, by = c("CpG" = "CpG", "SNP" = "SNP"), suffix = c(".AFR", ".EUR"))
common_cpgs <- unique(data2$CpG) 
length(common_cpgs) # 880108
#common_cpgs2 <- intersect(final_cis_AFA$CpG, final_cis_CAU$CpG)
final2_cis_CAU <- final_cis_CAU[final_cis_CAU$CpG %in% common_cpgs, ]
final2_cis_AFA <- final_cis_AFA[final_cis_AFA$CpG %in% common_cpgs, ]

set.seed(123)
sampled_df_CAU <- final2_cis_CAU[sample(nrow(final2_cis_CAU), size = 5000), ]
set.seed(123)
sampled_df_AFA <- final2_cis_AFA[sample(nrow(final2_cis_AFA), size = 5000), ]

sampled_df_CAU = sampled_df_CAU[, c("CpG", "SNP", "pvalue", "beta", "Z", "A1", "A2", "distance", "MAF", "chr", "annotation")]
sampled_df_AFA = sampled_df_AFA[, c("CpG", "SNP", "pvalue", "beta", "Z", "A1", "A2", "distance", "MAF", "chr", "annotation")]

set.seed(123)
#sampled_df_all <- data2[sample(nrow(data2), size = 5000), ]
#saveRDS(sampled_df_all, file = "/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/sampled_data_all.RDS")

saveRDS(sampled_df_CAU, file = paste0(res_dir_CAU,"sampled_data_CAU_160.RDS"))
saveRDS(sampled_df_AFA, file = paste0(res_dir_AFA,"sampled_data_AFA.RDS"))
# Or if only want the sampled_data, it can be read directly
sampled_df_CAU <- readRDS(paste0(res_dir_CAU,"sampled_data_CAU_160.RDS"))
sampled_df_AFA <- readRDS(paste0(res_dir_AFA,"sampled_data_AFA.RDS"))


# Combine the two data frames and add a column to indicate the group
sampled_df_CAU$group <- "EUR"
#remove the first column of sampled_df_CAU
sampled_df_AFA$group <- "AFR"
combined_df <- rbind(sampled_df_CAU, sampled_df_AFA)
combined_df$MAF_transformed <- combined_df$MAF * (1 - combined_df$MAF)
combined_df$effect_size <- abs(combined_df$beta)

##### a. MAF(1-MAF) density plot
p <- ggplot(combined_df, aes(x = MAF_transformed, y = effect_size, color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "MAF(1-MAF)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EUR" = "#E76254", "AFR" = "#528FAD"))

ggsave(paste0(save_dir,"a-MAF-beta-same-samplesize.jpg"))

# get the coefficient of the linear model
fit1 = lm(effect_size ~ MAF_transformed, data = combined_df[combined_df$group == "EUR",]) # coeff: -2.04 pval: <2e-16
fit2 = lm(effect_size ~ MAF_transformed, data = combined_df[combined_df$group == "AFR",])

p <- ggplot(combined_df, aes(x = MAF_transformed, y = -log(pvalue), color = group)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#F7AA58") +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "MAF(1-MAF)", y = "-log10(P-value)") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
  scale_color_manual(values = c("EUR" = "#E76254", "AFR" = "#528FAD"))

ggsave(paste0(save_dir,"a-MAF-pval.jpg"))

##### b. Scatter plot distance 
# Create the plot
p <- ggplot(combined_df, aes(x = distance, y = -log(pvalue))) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "-log10(P-value)") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), # Adjusted space for x-axis label
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Adjusted space for y-axis label
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"), # Set background color to white
        panel.grid.major = element_line(color = "grey90"), # Set major grid lines to a light grey
        panel.grid.minor = element_line(color = "grey95"), # Set minor grid lines to a very light grey
        axis.line = element_line(color = "black", linewidth = 0.3)) + # Set axis lines to black
  scale_color_manual(values = c("EUR" = "#E76254", "AFR" = "#528FAD")) +
  scale_x_continuous(breaks = c(-500000, 500000), labels = c("-500,000", "500,000"))

ggsave(paste0(save_dir,"b-distance-pval.jpg"))

p <- ggplot(combined_df, aes(x = distance, y = effect_size)) +
  geom_point(aes(color = group), alpha = 0.5) +
  facet_wrap(~ group, scales = "fixed", ncol = 2) +
  labs(x = "Distance of SNP from DNAm site (bp)", y = "Effect size") +
  theme(strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        axis.title.x = element_text(size = 14, margin = margin(t = 10)), # Adjusted space for x-axis label
        axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Adjusted space for y-axis label
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"), # Set background color to white
        panel.grid.major = element_line(color = "grey90"), # Set major grid lines to a light grey
        panel.grid.minor = element_line(color = "grey95"), # Set minor grid lines to a very light grey
        axis.line = element_line(color = "black", linewidth = 0.3)) + # Set axis lines to black
  scale_color_manual(values = c("EUR" = "#E76254", "AFR" = "#528FAD")) +
  scale_x_continuous(breaks = c(-500000, 500000), labels = c("-500,000", "500,000"))

ggsave(paste0(save_dir,"b-distance-beta-same-samplesize.jpg"))