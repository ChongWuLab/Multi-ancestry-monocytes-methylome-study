library(ggplot2)

res_fm_CAU = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/08-finemapping/CAU/unique_variants_count.RDS")
res_fm_CAU = unlist(res_fm_CAU)
# remove those names = ""
res_fm_CAU = res_fm_CAU[names(res_fm_CAU) != ""]
median(res_fm_CAU)
IQR(res_fm_CAU)

res_fm_AFA = readRDS("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/08-finemapping/AFA/unique_variants_count.RDS")
res_fm_AFA = unlist(res_fm_AFA)
# remove those names = ""
res_fm_AFA = res_fm_AFA[names(res_fm_AFA) != ""]
median(res_fm_AFA)
IQR(res_fm_AFA)


wilcox.test(res_fm_AFA, res_fm_CAU)

res_df_AFA = data.frame(CpG = names(res_fm_AFA), num_variants = res_fm_AFA, Population = "AA")
res_df_CAU = data.frame(CpG = names(res_fm_CAU), num_variants = res_fm_CAU, Population = "EA")
combined_df = rbind(res_df_AFA, res_df_CAU)
# plot
p <- ggplot(combined_df, aes(x = Population, y = num_variants, fill = Population)) +
  geom_violin(adjust = 1, alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c("AA" = "#528FAD", "EA" = "#E76254")) +
  theme_minimal() +
  labs(x = "", y = "Number of variants in credible sets", fill = "") +
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
ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/08-finemapping/1f-num_variants.jpg", p, width = 5, height = 5)
