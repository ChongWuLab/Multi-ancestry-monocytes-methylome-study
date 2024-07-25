suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyverse)))
library(dplyr)

setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
dat.wgbs.dir = "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"
dat.wgs.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/genetic_data/AFA/dataAll/rs_tabfile/"
res.dir = "/rsrch5/scratch/biostatistics/wzhang24/GWAS/res/AFA/02-heritability-V2/"
res.cov.dir = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA/"
gcta = "/rsrch5/home/biostatistics/wzhang24/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
cov_dir="/rsrch5/home/biostatistics/wzhang24/data/Processed_data/cov/AFA"
qcovar=file.path(cov_dir,"qcovar.txt")
covar=file.path(cov_dir,"covar.txt")


args = commandArgs(TRUE)
indx1 = as.numeric(args[[1]])

# Read methylation data
wgbs.file = list.files(dat.wgbs.dir,pattern=paste0("chunk_",indx1,".RDS"))



# get the chr
chr.id <- as.numeric(str_extract(wgbs.file, "(?<=chr)\\d+"))
bfile = paste0(dat.wgs.dir,"AFA_chr",chr.id)

#####process pheno
allpheno = readRDS(paste0(dat.wgbs.dir, wgbs.file))
#rownames(allpheno) = paste0("CpG",allpheno$Index)
extract = allpheno[, 6:ncol(allpheno)] %>% t() %>% as.matrix()
#annot = allpheno[,1:4]



##### process SNP
# read bim
bim_all = fread(paste0(dat.wgs.dir,"AFA_chr",chr.id,".bim")) %>% as.data.frame()
colnames(bim_all) = c("CHR","SNP","V3","P0","A0","A1")



#create directory to store the result 
filename = gsub(".RDS","",wgbs.file)
res_dir <- paste0(res.dir,filename)
if (!file.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
}
    
tmp_dir <- paste0("/rsrch5/home/biostatistics/wzhang24/mQTL_project/tmp/",filename) 
if (!file.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
}

##### estimate cis heritability
files = list.files(res_dir)
cpg_list = colnames(extract)
for (i in 1:ncol(extract)) {
    pheno_pos = allpheno[i,"Pos38"]
    cpg.name = cpg_list[i]
    #if (is.na(pheno_pos)) {
    #    cat("NA position for CpG ", i, cpg.name,"\n")
    #    next   
    #}
    # Check if the file already exists    
    name = paste0(cpg.name,".hsq")
    file_path = paste0(res_dir,"/",name)
    check=FALSE
    if (file.exists(file_path)) {
        # Try to read the file
        can_read <- tryCatch({
            data <- read.table(file_path, header = TRUE, fill = TRUE)
            h2_value <- data$Variance[data$Source == "V(G)/Vp"]
            # Check the length of h2_value
            if (length(h2_value) > 0) {
                check=TRUE  # File can be read and data is valid, return TRUE
            } else {
                check=FALSE  # h2_value has length 0, return FALSE
            }
        }, error = function(e) {
            check=FALSE  # There was an error, return FALSE
        })
        # If the file can be read and data is valid, skip the iteration
        if (can_read) {
            check=TRUE
        } else {
            check=FALSE
        }
    }
    if (check) {
        cat("File can be read and data is valid, skipping iteration: ", i, "\n")
        next
    }
    # Create the phenotype file
    pheno = as.data.frame(extract[,i])
    pheno = cbind(rownames(pheno),pheno)

    pheno_file = file.path(tmp_dir,"pheno.txt")
    cis_snp_file = file.path(tmp_dir,"cis_snps_list")
    out_prefix = file.path(tmp_dir,"AFA_cis")

    write.table(pheno,pheno_file, row.names=TRUE, col.names=FALSE, quote=FALSE,sep="\t")
    
    # get the cis SNPs
    cis_snps <- bim_all[(bim_all$P0 >= (pheno_pos - 1e6)) & (bim_all$P0 <= (pheno_pos + 1e6)), ]
    # save the cis SNPs ID
    write.table(cis_snps$SNP, cis_snp_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
    # Use PLINK to subset the data
    sys = sprintf("module load plink/1.90-beta && plink --bfile %s --extract %s --make-bed --out %s > /dev/null 2>&1", bfile, cis_snp_file, out_prefix)
    system(sys)
    
    
    #1. generate grm for cis only SNPs
    sys = sprintf("module load plink/1.90-beta && plink --allow-no-sex --bfile %s --make-grm-bin --out %s --threads 16 > /dev/null 2>&1", out_prefix, out_prefix)
    system(sys)
    #2. estimate heritability
    #sys = paste0("%s --grm AFA_cis"," --pheno pheno.txt --qcovar ${qcovar} --covar ${covar} --out %s --reml --threads 16")
    sys <- sprintf("%s --grm %s --pheno %s --qcovar %s --covar %s --out %s --reml > /dev/null 2>&1",
               gcta, out_prefix, pheno_file,
               qcovar, covar,
               file.path(res_dir, cpg.name))
    
    system(sys)

    cat("Finish Iteration:", i, "\n")
}

unlink(tmp_dir, recursive = TRUE)
cat("Finish processing chunk ",indx1,"\n")

