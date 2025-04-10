##############################################################################
# Load libraries
##############################################################################
library(ggplot2)
library(cowplot)  # for plot_grid (install with install.packages("cowplot") if needed)
library(dplyr)
##############################################################################
# 1) Process EA (GoDMC) data
##############################################################################

# Read EA data
final_data_ea <- readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/replication/replication_data_GoDMC_noclump.RDS")

# Find FDR threshold for EA
fdr_threshold_ea <- max(
  final_data_ea$pvalue[p.adjust(final_data_ea$pvalue, method = "fdr") < 0.05],
  na.rm = TRUE
)
cat("FDR threshold (EA):", fdr_threshold_ea, "\n")

# Add indicator whether the EA association replicates
final_data_ea$replicate <- final_data_ea$pval < fdr_threshold_ea

# Indicate if the beta and beta_a1 have the same sign
final_data_ea$same_sign <- final_data_ea$beta * final_data_ea$beta_a1 > 0

plot_data_ea = sample_n(final_data_ea, 10000)
# Create EA plot
p_ea <- ggplot(plot_data_ea, aes(x = beta, y = beta_a1, color = replicate)) +
  geom_point(size = 1, alpha = 0.3) +  # <--- alpha added
  scale_color_manual(
    name   = "Replication", 
    values = c("TRUE" = "#B36A6F", "FALSE" = "#6F6A6F"),
    labels = c("TRUE" = "True", "FALSE" = "False")
  ) +
  labs(
    x = "Effect size of EA (WGBS study)", 
    y = "Effect size of EA (GoDMC study)"
  ) +
  theme_minimal() +
  theme(
    legend.position      = c(0.001, 1.01),
    legend.justification = c(0, 1),
    legend.title         = element_text(size = 14),
    legend.text          = element_text(size = 14),
    axis.text.x          = element_text(hjust = 1, size = 14),
    axis.text.y          = element_text(size = 14),
    axis.title.x         = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y         = element_text(size = 14, margin = margin(r = 10)),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.height    = unit(1.15, "lines"),
    legend.key.width     = unit(1.15, "lines"),
    axis.ticks           = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length    = unit(0.2, "cm")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")

# Print replicate rate & correlation for EA
replicate_data_ea <- final_data_ea[final_data_ea$pval < fdr_threshold_ea, ]
replicate_rate_ea <- nrow(replicate_data_ea) / nrow(final_data_ea)
cat("Replication rate (EA):", replicate_rate_ea, "\n")

cat("Correlation (EA):\n")
print(cor.test(final_data_ea$beta, final_data_ea$beta_a1))

cat("Percentage with same sign (EA):", 
    sum(final_data_ea$same_sign) / nrow(final_data_ea), "\n")

##############################################################################
# 2) Process AA (GENOA) data
##############################################################################

# Read AA data
final_data_aa <- readRDS("/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/replication/replication_data_GENOA_AFA.RDS")

# Find FDR threshold for AA
fdr_threshold_aa <- max(
  final_data_aa$pvalue[p.adjust(final_data_aa$pvalue, method = "fdr") < 0.05],
  na.rm = TRUE
)
cat("FDR threshold (AA):", fdr_threshold_aa, "\n")

# Add indicator whether the AA association replicates
final_data_aa$replicate <- final_data_aa$p_wald < fdr_threshold_aa

# Indicate if the beta.x and beta.y have the same sign
# (In the GENOA dataset: beta.x => WGBS study, beta.y => GENOA study)
final_data_aa$same_sign <- final_data_aa$beta.x * final_data_aa$beta.y > 0

plot_data_aa = sample_n(final_data_aa, 10000)
# Create AA plot
p_aa <- ggplot(plot_data_aa, aes(x = beta.x, y = beta.y, color = replicate)) +
  geom_point(size = 1, alpha = 0.3) +
  scale_color_manual(
    name   = "Replication", 
    values = c("TRUE" = "#B36A6F", "FALSE" = "#6F6A6F"),
    labels = c("TRUE" = "True", "FALSE" = "False")
  ) +
  labs(
    x = "Effect size of AA (WGBS study)", 
    y = "Effect size of AA (GENOA study)"
  ) +
  theme_minimal() +
  theme(
    legend.position      = c(0.001, 1.01),
    legend.justification = c(0, 1),
    legend.title         = element_text(size = 14),
    legend.text          = element_text(size = 14),
    axis.text.x          = element_text(hjust = 1, size = 14),
    axis.text.y          = element_text(size = 14),
    axis.title.x         = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y         = element_text(size = 14, margin = margin(r = 10)),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.height    = unit(1.15, "lines"),
    legend.key.width     = unit(1.15, "lines"),
    axis.ticks           = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length    = unit(0.2, "cm")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")

# Print replicate rate & correlation for AA
replicate_data_aa <- final_data_aa[final_data_aa$p_wald < fdr_threshold_aa, ]
replicate_rate_aa <- nrow(replicate_data_aa) / nrow(final_data_aa)
cat("Replication rate (AA):", replicate_rate_aa, "\n") # 78.1%

cat("Correlation (AA):\n")
print(cor.test(final_data_aa$beta.x, final_data_aa$beta.y))

cat("Percentage with same sign (AA):", 
    sum(final_data_aa$same_sign) / nrow(final_data_aa), "\n")

##############################################################################
# 3) Combine EA and AA plots into one figure
##############################################################################

# Combine the two plots side by side
combined_plot <- plot_grid(p_ea, p_aa, 
                           ncol = 2, 
                           align = "hv")

# Save combined plot
save_dir <- "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/figures/"
ggsave(
  filename = paste0(save_dir, "Fig2f-combined_EA_AA_replication.jpg"),
  plot     = combined_plot,
  width    = 10,   # Adjust as needed
  height   = 5,    # Adjust as needed
  dpi      = 300
)

