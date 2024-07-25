library(data.table)
library(dplyr)
library(stringr)
library(parallel)
library(RNOmni)

dat.wgbs.dir <- "/rsrch5/home/biostatistics/wzhang24/data/WGBS/processed_AFA"
out.dir <- "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"

files <- list.files(dat.wgbs.dir, full.names = TRUE)
ids = read.table("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/ids/AFA/intersect_ids.txt") 
ids = ids$V1
cov.all = as.data.frame(fread("/rsrch5/home/biostatistics/wzhang24/data/WGBS/AFA_185meta_w20Kwgbs.csv"))

id_mapping = setNames(cov.all$LOS_id, cov.all$Subj)

# Function to process each file
process_file <- function(file, out.dir, i, ids, id_mapping) {
  filename <- tools::file_path_sans_ext(basename(file))
  chromosome_match <- regmatches(filename, regexpr("chr(\\d+)|chr(X|Y)", filename))
  chr.id <- gsub("chr", "", chromosome_match)
  liftover = fread(paste0("/rsrch5/home/biostatistics/wzhang24/data/WGBS/hg19_to_hg38/hglft_hg38_genome_chr",chr.id,".bed"))
  liftover = liftover[-1,]
  colnames(liftover) = c("Chr","Pos38","V3","Index")
  allpheno <- fread(file, data.table = FALSE) 
  allpheno = left_join(allpheno, liftover, by="Index")
  rownames(allpheno) <- paste0("CpG", allpheno$Index)
  # remove the rows with NAs in allpheno$Pos38
  allpheno = allpheno[!is.na(allpheno$Pos38),]
  
  extract <- allpheno[, 6:(ncol(allpheno)-3)]
  #replace MM id with NO id
  colnames(extract) = id_mapping[colnames(extract)]
  extract = extract[,ids]

  annot <- allpheno[, c(1, 2, 3, 192, 193)]
  colnames(annot) <- c("Index", "Chr", "Pos19","Pos38", "Pos38.right")
  
  extract <- apply(extract, 2, RankNorm)
  dat = cbind(annot, extract)
  # Save the normalized data
  saveRDS(dat, file = paste0(out.dir, "AFA_chr", chr.id, "_chunk_", i, ".RDS"))
}


# Detect the number of available cores
no_cores <- detectCores() - 1
cat("Number of cores available: ", no_cores, "\n")
# Setup cluster
cl <- makeCluster(no_cores)

# Apply process_file function in parallel
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  library(stringr)
  library(RNOmni)
})

# Export functions and objects to cluster
clusterExport(cl, c("process_file", "files", "out.dir", "ids", "id_mapping"))

# Apply process_file function in parallel
parLapply(cl, seq_along(files), function(i) {
  process_file(files[i], out.dir, i, ids=ids, id_mapping=id_mapping)
})

# Stop the cluster
stopCluster(cl)

#data = readRDS("/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_CAU/CAU_chr9_chunk_12884.RDS") %>% as.data.frame()



