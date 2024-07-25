suppressWarnings(suppressMessages(library(data.table)))
library(dplyr)
setwd("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/")
cov = as.data.frame(fread("/rsrch5/home/biostatistics/wzhang24/data/WGBS/AFA_185meta_w20Kwgbs.csv"))


##### process cov
rownames(cov) = cov$LOS_id
#cov = cov[,colnames(cov) %in% c("Age","BMI","Gender","Race","Smoke","Alcohol_Drinking","Exercise_Regular")]
cov = cov[,colnames(cov) %in% c("Age","BMI","Smoke","Alcohol_Drinking","Monocytes","Bcells", "Neutrophils")]

ids = read.table("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/ids/AFA/intersect_ids.txt") 
ids = ids$V1
cov = cov[ids,]

# add 5 genetic PCs
pc = as.data.frame(fread("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/pca.eigenvec"))
pc = pc[,2:12]
#colnames(pc) = c("Subj","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
colnames(pc) = c("Subj","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
rownames(pc) = pc$Subj
pc = pc[,-1]

#add ten non-genetic PCs
ng_pcs = readRDS("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/nongenetic_PCs/AFA/pcs_most_variated.RDS") %>% as.data.frame()
ng_pcs = ng_pcs[,1:10]
colnames(ng_pcs) = paste0("ng",colnames(ng_pcs))

#Combine cov and PCs
cov = cbind(cov,pc)
cov = cbind(cov,ng_pcs)

qcovar = cov[,c("Age","BMI","Bcells","Monocytes","Neutrophils","PC1","PC2","PC3","PC4","PC5","ngPC1","ngPC2","ngPC3","ngPC4","ngPC5","ngPC6","ngPC7","ngPC8","ngPC9","ngPC10")]
covar = cov[,c("Smoke","Alcohol_Drinking")]

covar$Smoke = as.factor(covar$Smoke)
covar$Alcohol_Drinking = as.factor(covar$Alcohol_Drinking)


qcovar = cbind(rownames(qcovar), qcovar)
covar = cbind(rownames(covar), covar)
write.table(qcovar,"qcovar.txt",row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(covar,"covar.txt", row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(cov,"cov.txt", row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
