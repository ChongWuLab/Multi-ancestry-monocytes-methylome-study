library(data.table)
library(dplyr)
library(ggplot2)
##### read data
AFA.c.path = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/04-modelling-V2/"
CAU.c.path = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/04-modelling-V2/"
CAU.nc.path = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/CAU/04-modelling-V2-nc/"

files.AFA.c = list.files(AFA.c.path,full.names = TRUE)
#files.CAU.c = list.files(CAU.c.path,full.names = TRUE)
#files.CAU.nc = list.files(CAU.nc.path,full.names = TRUE)
files.CAU.c = c()
for (i in 1:100){
  files.CAU.c = c(files.CAU.c, paste0(CAU.c.path,"chr1_",i,".RData"))
}
files.CAU.nc = c()
for (i in 1:100){
  files.CAU.nc = c(files.CAU.nc, paste0(CAU.nc.path,"chr1_",i,".RData"))
}


CpG_list <- lapply(files.AFA.c[1:100], function(f) {
  load(f)
  return(as.data.table(CpG_mat))
})
AFA.c <- rbindlist(CpG_list)

CpG_list <- lapply(files.CAU.c, function(f) {
  load(f)
  return(as.data.table(CpG_mat))
})
CAU.c <- rbindlist(CpG_list)

CpG_list <- lapply(files.CAU.nc, function(f) {
  load(f)
  return(as.data.table(CpG_mat))
})
CAU.nc <- rbindlist(CpG_list)

AFA.c = AFA.c[,c("CpG","R2")]
CAU.c = CAU.c[,c("CpG","R2")]

# merge the two data by inner join
merged_df <- inner_join(AFA.c, CAU.c, by = "CpG", suffix=c(".AFA",".CAU")) %>%
  filter(complete.cases(.))

##### AFA clumping vs CAU clumping
ggplot(merged_df, aes(x = R2.AFA, y = R2.CAU)) +
  geom_point(alpha = 0.6) +  # Adding points with some transparency
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  # y=x line
  labs(x = "R2 for AFA",
       y = "R2 for CAU",
       title = "Scatter Plot of R2 Values") +
  theme_minimal()  # Using a minimal theme for better visual
ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/04-modelling/R2_scatter.png", width = 6, height = 6, dpi = 300)

##### CAU clumping vs CAU noclumping
length(intersect(CAU.c$CpG, CAU.nc$CpG))
merged_df <- inner_join(CAU.c, CAU.nc, by = "CpG", suffix=c(".c",".nc")) 
# fill na with 0
merged_df[is.na(merged_df)] <- 0
ggplot(merged_df, aes(x = R2.c, y = R2.nc)) +
  geom_point(alpha = 0.6) +  # Adding points with some transparency
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  # y=x line
  labs(x = "R2 for EUR clumping",
       y = "R2 for EUR noclumping",
       title = "Scatter Plot of R2 Values") +
  theme_minimal()  # Using a minimal theme for better visual
ggsave("/rsrch5/home/biostatistics/wzhang24/mQTL_project/Results/04-modelling/R2_scatter_CAU_fill.png", width = 6, height = 6, dpi = 300)
