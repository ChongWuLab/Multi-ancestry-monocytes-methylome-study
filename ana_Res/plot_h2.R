# The preparation part. Do not need to do this if alread done before.
library(ggplot2)
library(data.table)
library(dplyr)
#library(rtracklayer)
h2_df_CAU = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-CAU/all_h2.RDS")
h2_df_CAU$h2 = as.numeric(h2_df_CAU$h2)
h2_df_CAU$bin <- cut(h2_df_CAU$h2, breaks = seq(0, 1, by = 0.1),right=FALSE)

h2_df_AFA = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-AFA/all_h2.RDS")
h2_df_AFA$h2 = as.numeric(h2_df_AFA$h2)
h2_df_AFA$bin <- cut(h2_df_AFA$h2, breaks = seq(0, 1, by = 0.1),right=FALSE)

# Calculate the number of unique CpG sites in each bin
cpg_count_CAU <- aggregate(CpG ~ bin, data = h2_df_CAU, FUN = function(x) length(unique(x)))
cpg_count_AFA <- aggregate(CpG ~ bin, data = h2_df_AFA, FUN = function(x) length(unique(x)))

# Find out which CpGs have mQTLs by checking if they appear in the mQTL dataframe
mqtl_df_CAU <- readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/CAU/final_data_CAU.RDS")
mqtl_df_AFA <- readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/03-mQTL-res/AFA/final_data_AFA.RDS")
h2_df_CAU$has_mQTL <- h2_df_CAU$CpG %in% mqtl_df_CAU$CpG
h2_df_AFA$has_mQTL <- h2_df_AFA$CpG %in% mqtl_df_AFA$CpG

saveRDS(h2_df_CAU,file="/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-CAU/all_h2.RDS")
saveRDS(h2_df_AFA,file="/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-AFA/all_h2.RDS")

#================================================================================================
# Begin with this part
#################################################################################################
library(ggplot2)
library(data.table)
library(dplyr)
#h2_df_AFA = h2_df_AFA %>% filter(h2 > 0) #25680499
#h2_df_CAU = h2_df_CAU %>% filter(h2 > 0) #25678678 
#saveRDS(h2_df_CAU,file="/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-CAU/all_h2.RDS")
#saveRDS(h2_df_AFA,file="/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-AFA/all_h2.RDS")
h2_df_CAU = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-CAU/all_h2.RDS")
h2_df_AFA = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2-AFA/all_h2.RDS")

# t-test
t.test(h2_df_CAU$h2, h2_df_AFA$h2, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95) # p < 2.2e-16

# check the 25% quantile
quantile(h2_df_CAU$h2, probs = 0.25)
quantile(h2_df_AFA$h2, probs = 0.25)
# check the 75% quantile
quantile(h2_df_CAU$h2, probs = 0.75)
quantile(h2_df_AFA$h2, probs = 0.75)

h2_df_CAU$bin <- cut(h2_df_CAU$h2, breaks = c(seq(0, 0.6, by = 0.1), 1), right=FALSE)
h2_df_AFA$bin <- cut(h2_df_AFA$h2, breaks = c(seq(0, 0.6, by = 0.1), 1), right=FALSE)

##### draw the figure of the relationship between h2 and p-value
h2_df_CAU$significance = ifelse(h2_df_CAU$pvalue < 0.05, "TRUE", "FALSE")
h2_df_AFA$significance = ifelse(h2_df_AFA$pvalue < 0.05, "TRUE", "FALSE")

h2_thresholds <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5)

# Function to calculate proportions for a given dataset and thresholds
calculate_proportions <- function(df, population) {
  lapply(h2_thresholds, function(threshold) {
    data <- df %>% 
      filter(h2 > threshold) %>% 
      summarize(proportion = mean(pvalue < 0.05, na.rm = TRUE))
    data$threshold = threshold
    data$Population = population
    return(data)
  }) %>% bind_rows()
}

# Calculate proportions for each population
proportions_df_CAU <- calculate_proportions(h2_df_CAU, 'EA')
proportions_df_AFA <- calculate_proportions(h2_df_AFA, 'AA')

# Combine the results
combined_df <- bind_rows(proportions_df_CAU, proportions_df_AFA)

# Plotting
ggplot(combined_df, aes(x = factor(threshold), y = proportion, fill = Population)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_x_discrete(labels = paste(">", h2_thresholds)) +
  labs(x = expression(italic(cis) * "-" * h^2), y = "Proportion of P-values < 0.05") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  scale_fill_manual(values = c("EA" = "#E76254", "AA" = "#72BCD5"))

ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2/2b-h2_pvalue_combined.png", width=8, height=5, units = "in", dpi = 300)



##### Historgram with mQTL together
# Calculate the percentage of CpG sites that have mQTLs in each bin
cpg_with_mQTL_count_CAU <- aggregate(has_mQTL ~ bin, data = h2_df_CAU, FUN = function(x) sum(x)/length(x))
cpg_with_mQTL_count_AFA <- aggregate(has_mQTL ~ bin, data = h2_df_AFA, FUN = function(x) sum(x)/length(x))

# Calculate the number of unique CpG sites in each bin
cpg_count_CAU <- aggregate(CpG ~ bin, data = h2_df_CAU, FUN = function(x) length(unique(x)))
cpg_count_AFA <- aggregate(CpG ~ bin, data = h2_df_AFA, FUN = function(x) length(unique(x)))

# merge the counts with the total counts to get the percentage
percentage_df_CAU <- merge(cpg_count_CAU, cpg_with_mQTL_count_CAU, by = "bin")
percentage_df_AFA <- merge(cpg_count_AFA, cpg_with_mQTL_count_AFA, by = "bin")
percentage_df_CAU$percentage <- percentage_df_CAU$has_mQTL * 100
percentage_df_AFA$percentage <- percentage_df_AFA$has_mQTL * 100


# Add an identifier for the population to each dataset
percentage_df_CAU$Population <- 'EA'
percentage_df_AFA$Population <- 'AA'

# Combine the two datasets
combined_df <- rbind(percentage_df_CAU, percentage_df_AFA)

# Now plot with ggplot2
ggplot(combined_df, aes(x = bin, y = percentage, color = Population, group = Population)) +
  geom_line() +
  geom_point(aes(shape = Population), size = 3) +
  scale_shape_manual(values = c("EA" = 16, "AA" = 16)) +
  scale_color_manual(values = c("EA" = "#E76254", "AA" = "#72BCD5")) +  # Change these colors as needed
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = expression(italic(cis) * "-" * h^2), y = "% of CpGs having meQTLs") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle=30, hjust=1, size = 12),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2/2c-h2_cpg_combined.png",width=5,height=5, units = "in", dpi = 300)
 
##### Histogram of h2 together
# Add a distinguishing column to each dataframe
h2_df_AFA_filter <- h2_df_AFA[h2_df_AFA$h2 > 0.01,]
h2_df_CAU_filter <- h2_df_CAU[h2_df_CAU$h2 > 0.01,]
h2_df_CAU_filter$Population <- "EA"
h2_df_AFA_filter$Population <- "AA"

combined_h2_df <- rbind(h2_df_CAU_filter, h2_df_AFA_filter)

# Use ggplot2 to create the side by side histogram
ggplot(combined_h2_df, aes(x = h2,fill = Population)) +
  geom_histogram(bins = 20, position = 'dodge', boundary = 0, binwidth = 0.05,color="black") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = expression(italic(cis) * "-" * h^2), y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(2, "lines"),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  scale_fill_manual(values = c("EA" = "#E76254", "AA" = "#72BCD5"))
ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2/2a-h2_hist_combined.png", width=8, height=5, units = "in", dpi = 300)

##### Boxplot of h2 together
#h2_df_CAU = h2_df_CAU %>% filter(h2 > 0.01)
#h2_df_AFA = h2_df_AFA %>% filter(h2 > 0.01)
h2_df_CAU$Category <- ifelse(h2_df_CAU$has_mQTL, "cis-meQTL CpGs", "non cis-meQTL CpGs")
h2_df_AFA$Category <- ifelse(h2_df_AFA$has_mQTL, "cis-meQTL CpGs", "non cis-meQTL CpGs")
# Sample 1000 rows from each dataframe
h2_df_CAU_sample <- h2_df_CAU %>% sample_n(50000)
h2_df_AFA_sample <- h2_df_AFA %>% sample_n(50000)

# Add a distinguishing column to each dataframe
#h2_df_CAU_all <- h2_df_CAU_sample
#h2_df_CAU_all$Category <- "All CpGs"
#h2_df_AFA_all <- h2_df_AFA_sample
#h2_df_AFA_all$Category <- "All CpGs"

h2_df_AFA_sample$Population <- "AA"
h2_df_CAU_sample$Population <- "EA"
#h2_df_CAU_all$Population <- "EUR"
#h2_df_AFA_all$Population <- "AFR"
# Combine the two dataframes

#combined_df <- rbind(h2_df_CAU_sample, h2_df_AFA_sample, h2_df_CAU_all, h2_df_AFA_all)
combined_df <- rbind(h2_df_CAU_sample, h2_df_AFA_sample)

p <- ggplot(combined_df, aes(x = Category, y = h2, fill = Population)) +
  #geom_violin(trim = TRUE, adjust = 1, alpha = 0.3, position = position_dodge(width = 0.8)) + # Trim the violin plot
  geom_boxplot(width = 0.4, position = position_dodge(width = 0.5), outlier.shape = NA) +
  scale_fill_manual(values = c("AA" = "#528FAD", "EA" = "#E76254")) +
  theme_minimal() +
  labs(x = "", y = expression(italic(cis) * "-" * h^2), fill = "") +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = -14, b = -8), 
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 12, hjust = 0.5), # Center the x-axis labels
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2/2d-h2_boxplot_combined2.jpg", width=5, height=5, units = "in", dpi = 300)



#================================================================================================
##### Boxplot of h2 from EPIC and HM450
HM450 <- read.csv("/rsrch5/home/biostatistics/wzhang24/data/WGBS/humanmethylation450_15017482_v1-2.csv")
epic850k <- fread("/rsrch5/home/biostatistics/wzhang24/data/WGBS/EPIC-8v2-0_A1.csv", skip=7)
HM450 <- HM450[, c("CHR", "MAPINFO")]
epic850k <- epic850k[, c("CHR", "MAPINFO")]
hg19_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19/combined_data.RDS")
hg38_annot <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg38_annotation.RDS")
hg38_annot <- hg38_annot[, c("CpG", "pos38")]

h2_df_CAU <- left_join(h2_df_CAU, hg19_annot, by="CpG")
h2_df_CAU <- left_join(h2_df_CAU, hg38_annot, by="CpG")
h2_df_CAU$Chr2 <- gsub("chr", "", h2_df_CAU$Chr)
h2_df_CAU_hm450 <- h2_df_CAU %>% inner_join(HM450, by = c("Chr2" = "CHR", "Pos19" = "MAPINFO"))
h2_df_CAU_epic850k <- h2_df_CAU %>% inner_join(epic850k, by = c("Chr" = "CHR", "pos38" = "MAPINFO"))
dim(h2_df_CAU_epic850k)
dim(h2_df_CAU_hm450)
h2_df_CAU_hm450 = h2_df_CAU_hm450[, c("CpG", "h2")]
h2_df_CAU_epic850k = h2_df_CAU_epic850k[, c("CpG", "h2")]

h2_df_AFA <- left_join(h2_df_AFA, hg19_annot, by="CpG")
h2_df_AFA <- left_join(h2_df_AFA, hg38_annot, by="CpG")
h2_df_AFA$Chr2 <- gsub("chr", "", h2_df_AFA$Chr)
h2_df_AFA_hm450 <- h2_df_AFA %>% inner_join(HM450, by = c("Chr2" = "CHR", "Pos19" = "MAPINFO"))
h2_df_AFA_epic850k <- h2_df_AFA %>% inner_join(epic850k, by = c("Chr" = "CHR", "pos38" = "MAPINFO"))
dim(h2_df_AFA_epic850k)
dim(h2_df_AFA_hm450)
h2_df_AFA_hm450 = h2_df_AFA_hm450[, c("CpG", "h2")]
h2_df_AFA_epic850k = h2_df_AFA_epic850k[, c("CpG", "h2")]

# some stats
mean(h2_df_AFA_hm450$h2)
IQR(h2_df_AFA_hm450$h2)
mean(h2_df_CAU_hm450$h2)
mean(h2_df_AFA_epic850k$h2)
mean(h2_df_CAU_epic850k$h2)


h2_df_AFA$Population <- "AA"
h2_df_CAU$Population <- "EA"
h2_df_AFA$Category <- "All"
h2_df_CAU$Category <- "All"
h2_df_AFA_hm450$Population <- "AA"
h2_df_CAU_hm450$Population <- "EA"
h2_df_AFA_hm450$Category <- "HM450k"
h2_df_CAU_hm450$Category <- "HM450k"
h2_df_AFA_epic850k$Population <- "AA"
h2_df_CAU_epic850k$Population <- "EA"
h2_df_AFA_epic850k$Category <- "EPIC900k"
h2_df_CAU_epic850k$Category <- "EPIC900k"

h2_df_AFA_sample <- h2_df_AFA %>% sample_n(100000)
h2_df_CAU_sample <- h2_df_CAU %>% sample_n(100000)
h2_df_CAU_sample <- h2_df_CAU_sample[, c("CpG", "h2", "Population", "Category")]
h2_df_AFA_sample <- h2_df_AFA_sample[, c("CpG", "h2", "Population", "Category")]
h2_df_AFA_hm450_sample <- h2_df_AFA_hm450 %>% sample_n(100000)
h2_df_CAU_hm450_sample <- h2_df_CAU_hm450 %>% sample_n(100000)
h2_df_AFA_epic850k_sample <- h2_df_AFA_epic850k %>% sample_n(100000)
h2_df_CAU_epic850k_sample <- h2_df_CAU_epic850k %>% sample_n(100000)

combined_df <- rbind(h2_df_CAU_sample, h2_df_AFA_sample, h2_df_CAU_hm450_sample, h2_df_AFA_hm450_sample, h2_df_CAU_epic850k_sample, h2_df_AFA_epic850k_sample)

p <- ggplot(combined_df, aes(x = Category, y = h2, fill = Population)) +
  #geom_violin(trim = TRUE, adjust = 1, alpha = 0.3, position = position_dodge(width = 0.8)) + # Trim the violin plot
  geom_boxplot(width = 0.4, position = position_dodge(width = 0.5), outlier.shape = NA) +
  #scale_fill_manual(values = c("AFR" = "#528FAD", "EUR" = "#E76254")) +
  scale_fill_manual(values = c("AA" = "#528FAD", "EA" = "#E76254")) +
  theme_minimal() +
  labs(x = "", y = expression(italic(cis) * "-" * h^2), fill = "") +
  coord_cartesian(ylim = c(0, 0.50)) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = -14, b = -8), 
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, hjust = 0.5), # Center the x-axis labels
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/02-h2/2e-h2_boxplot_combined3.jpg", width=5, height=5, units = "in", dpi = 300)


#### Some summary
# nonzero h2
median(h2_df_AFA$h2)
median(h2_df_CAU$h2)
#overlap
overlap = intersect(h2_df_AFA$CpG, h2_df_CAU$CpG)
length(overlap)
# h2 < 0.01
h2_df_AFA_filter = h2_df_AFA[h2_df_AFA$h2 < 0.01,]
h2_df_CAU_filter = h2_df_CAU[h2_df_CAU$h2 < 0.01,]
overlap = intersect(h2_df_AFA_filter$CpG, h2_df_CAU_filter$CpG)
dim(h2_df_AFA_filter)[1]/dim(h2_df_AFA)[1] #50.5%
dim(h2_df_CAU_filter)[1]/dim(h2_df_CAU)[1] #57.0%

# 0.01h2<=0.1
h2_df_AFA_filter = h2_df_AFA[h2_df_AFA$h2 >= 0.01 & h2_df_AFA$h2 <= 0.1,]
h2_df_CAU_filter = h2_df_CAU[h2_df_CAU$h2 >= 0.01 & h2_df_CAU$h2 <= 0.1,]
overlap = intersect(h2_df_AFA_filter$CpG, h2_df_CAU_filter$CpG)
dim(h2_df_AFA_filter)[1]/dim(h2_df_AFA)[1] #21.1%
dim(h2_df_CAU_filter)[1]/dim(h2_df_CAU)[1] #33.8%
# h2 > 0.1
h2_df_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.1,]
h2_df_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.1,]
overlap = intersect(h2_df_AFA_filter$CpG, h2_df_CAU_filter$CpG)
dim(h2_df_AFA_filter)[1]/dim(h2_df_AFA)[1] #28.3%
dim(h2_df_CAU_filter)[1]/dim(h2_df_CAU)[1] #9.2%
# h2 > 0.5
h2_df_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.5,]
h2_df_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.5,]
overlap = intersect(h2_df_AFA_filter$CpG, h2_df_CAU_filter$CpG)
dim(h2_df_AFA_filter)[1]/dim(h2_df_AFA)[1] #3.6%
dim(h2_df_CAU_filter)[1]/dim(h2_df_CAU)[1] #1.1%

# Among the CpGs with h2 > 0.1, how many of them have at least 5trone mQTL with p < 1e-8
h2_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.1,]
sum(h2_AFA_filter$has_mQTL==TRUE)/length(h2_AFA_filter$has_mQTL) # 20%
h2_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.1,]
sum(h2_CAU_filter$has_mQTL==TRUE)/length(h2_CAU_filter$has_mQTL) # 42%


# Among the CpGs with h2 > 0.5, how many of them have at least one mQTL with p < 1e-8
h2_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.5,]
sum(h2_CAU_filter$has_mQTL==TRUE)/length(h2_CAU_filter$has_mQTL) # 99%
h2_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.5,]
sum(h2_AFA_filter$has_mQTL==TRUE)/length(h2_AFA_filter$has_mQTL) # 86%

# h2 > 0.01 and pvalue < 0.05
h2_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.01, ]
h2_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.01, ]
h2_AFA_filter2 = h2_AFA_filter[h2_AFA_filter$pvalue < 0.05,]
h2_CAU_filter2 = h2_CAU_filter[h2_CAU_filter$pvalue < 0.05,]
merged_h2<- inner_join(h2_AFA_filter, h2_CAU_filter, by="CpG")
dim(h2_AFA_filter2)[1]/(dim(h2_AFA_filter)[1]) # 16.4%
dim(h2_CAU_filter2)[1]/(dim(h2_CAU_filter)[1]) # 18.2%
sum(h2_AFA_filter$has_mQTL==TRUE)/length(h2_AFA_filter$has_mQTL) #0.52
sum(h2_CAU_filter$has_mQTL==TRUE)/length(h2_CAU_filter$has_mQTL) #0.54


# h2 > 0.1 and pvalue < 0.05
h2_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.1, ]
h2_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.1, ]
h2_AFA_filter2 = h2_AFA_filter[h2_AFA_filter$pvalue < 0.05,]
h2_CAU_filter2 = h2_CAU_filter[h2_CAU_filter$pvalue < 0.05,]
dim(h2_AFA_filter2)[1]/(dim(h2_AFA_filter)[1]) # 28.6%
dim(h2_CAU_filter2)[1]/(dim(h2_CAU_filter)[1]) # 61.4%

# h2 > 0.5 and pvalue < 0.05
h2_AFA_filter = h2_df_AFA[h2_df_AFA$h2 > 0.5, ]
h2_CAU_filter = h2_df_CAU[h2_df_CAU$h2 > 0.5, ]
h2_AFA_filter2 = h2_AFA_filter[h2_AFA_filter$pvalue < 0.05,]
h2_CAU_filter2 = h2_CAU_filter[h2_CAU_filter$pvalue < 0.05,]
dim(h2_AFA_filter2)[1]/(dim(h2_AFA_filter)[1]) # 36.8%
dim(h2_CAU_filter2)[1]/(dim(h2_CAU_filter)[1]) # 55.1%

# cpgs having meQTLs
h2_AFA_filter = h2_df_AFA[h2_df_AFA$has_mQTL==TRUE,]
h2_CAU_filter = h2_df_CAU[h2_df_CAU$has_mQTL==TRUE,]
median(h2_AFA_filter$h2) #0.49
median(h2_CAU_filter$h2) #0.33