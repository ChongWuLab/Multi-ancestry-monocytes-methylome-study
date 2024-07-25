##### map to 1000G
library(data.table)
library(dplyr)
data.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"

bim = fread(paste0(data.wgs.dir,"AFA_chr_all.bim"))

bim_org = fread("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/dataAll.bim.original")
#merge bim and bim_org by V4(position) and V1(chr)
bim_all = bim %>% left_join(bim_org, by = c("V4" = "V4", "V1" = "V1"))

bim_rs = bim_all[,c("V1","V2.y","V3.x","V4","V5.x","V6.x")]


write.table(bim_rs, "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/AFA_chr_all.bim",row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

