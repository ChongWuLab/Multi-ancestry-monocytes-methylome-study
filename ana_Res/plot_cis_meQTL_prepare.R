library(data.table)
library(dplyr)
library(BEDMatrix)

res_dir = "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/"
save_dir = "/rsrch5/home/biostatistics/chongwulab/wzhang24/MWAS/Results/meQTL/figures/"

final_cis_CAU = readRDS(paste0(res_dir,"final_cis_fdr_0.01_CAU.RDS"))
final_cis_AFA <- readRDS(paste0(res_dir, "final_cis_fdr_0.01_AFA.RDS"))

# # Basic stats (results paragraph 2)
# nrow(final_cis_AFA) # 45143297
# nrow(final_cis_CAU) # 98259931

# length(unique(final_cis_AFA$gene)) # 3,634,383
# length(unique(final_cis_CAU$gene)) # 3,261,407
# length(unique(final_cis_AFA$gene))/25721231 # 14.1%
# length(unique(final_cis_CAU$gene))/25721231 # 12.7%

data1_AFA <- final_cis_AFA[, c("gene", "snps", "beta")] #45143297
data1_CAU <- final_cis_CAU[, c("gene", "snps", "beta")] #98259931
data2 <- inner_join(data1_AFA, data1_CAU, by = c("gene" = "gene", "snps" = "snps"), suffix = c(".AFR", ".EUR"))
common_cpgs <- unique(data2$gene) 
# length(common_cpgs) # 1,046,098
# length(common_cpgs)/length(unique(final_cis_AFA$gene)) # 28.8%
# length(common_cpgs)/length(unique(final_cis_CAU$gene)) # 32.1%
# get the sample data
# get the shared sampledata
final2_cis_CAU <- final_cis_CAU[final_cis_CAU$gene %in% common_cpgs, ]
final2_cis_AFA <- final_cis_AFA[final_cis_AFA$gene %in% common_cpgs, ]

set.seed(123)
sampled_df_CAU <- final2_cis_CAU[sample(nrow(final2_cis_CAU), size = 5000), ]
set.seed(123)
sampled_df_AFA <- final2_cis_AFA[sample(nrow(final2_cis_AFA), size = 5000), ]



saveRDS(sampled_df_AFA, file = paste0(save_dir, "sampled_df_AFA.RDS"))
saveRDS(sampled_df_CAU, file = paste0(save_dir, "sampled_df_CAU.RDS"))
