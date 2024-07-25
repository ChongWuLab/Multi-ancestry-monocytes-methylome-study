library(data.table)
library(parallel)
library(matrixStats)
library(BEDMatrix)
library(dplyr)
setwd("/rsrch5/home/biostatistics/wzhang24/mQTL_project/codes/AFA/")
#print(paste("Current working directory:", getwd()))
source("support_MWAS.R")


# Path to the directory where RDS files are stored
dir_path <- "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"
# List all RDS files
files <- list.files(dir_path, full.names = TRUE, pattern = "\\.RDS$")

# Function to read a single RDS file
read_rds_file <- function(file_path) {
  dt <- readRDS(file_path)
  dt[, 5:ncol(dt)]
}

# Set up a cluster using available cores (or specify the number you want to use)
# check if the PC file already exists
if (!file.exists("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/nongenetic_PCs/AFA/pcs_most_variated.RDS")) {
    cat("PCs file doesn't exists. Running from the beginning. \n")

    # Detect the number of available cores
    no_cores <- detectCores() - 1
    cat("Using", no_cores, "cores\n")
    cl <- makeCluster(no_cores)

    # Export the read_rds_file function to each node
    clusterExport(cl, varlist = c("read_rds_file"))

    # Read the files in parallel
    list_of_dfs <- parLapply(cl, files, read_rds_file)
    for (i in 1:length(list_of_dfs)) {
        data <- list_of_dfs[[i]]
    }
    cat("Finished reading", length(list_of_dfs), "files\n")
    # Stop the cluster
    stopCluster(cl)

    # Combine all data.tables into one
    combined_dt <- rbindlist(list_of_dfs)
    cat("Combined data.table of DNAm has dimensions", dim(combined_dt), "\n")

    #Identifying the most variable CpGs
    calculate_variances <- function(dt) {
    apply(dt, 1, var)  # Calculate variance for each row
    }

    variances <- rowVars(as.matrix(combined_dt)) # Calculate variance for each row
    cat("Finished calculating variances\n")

    # Select the top 2000 CpGs with the highest variance
    top_sites_indices <- order(variances, decreasing = TRUE)[1:20000]
    top_sites_matrix <- combined_dt[top_sites_indices, ] %>% t()



    # Run PCA
    pca_result <- prcomp(top_sites_matrix)

    # keep all PCs that cumulatively explained 80% of the variance.
    explained_variance <- summary(pca_result)$importance[2,]
    cum_explained_variance <- cumsum(explained_variance)
    num_pcs_to_keep <- which(cum_explained_variance >= 0.80)[1]
    pcs_to_keep <- pca_result$x[, 1:num_pcs_to_keep]
    cat("Keeping", num_pcs_to_keep, "PCs\n")
    # save the PCs to disk
    saveRDS(pcs_to_keep, file = "/rsrch5/home/biostatistics/wzhang24/data/Processed_data/nongenetic_PCs/AFA/pcs_most_variated.RDS")
    cat("Saved PCs with most variated 20000 CpGs and that cumulatively explained 80% of the variance.\n")
    cat("Saved PCs have dimensions", dim(pcs_to_keep), "\n")
} else{
    cat("PCs file already exists. Loading from disk. \n")
    pcs_to_keep <- readRDS("/rsrch5/home/biostatistics/wzhang24/data/Processed_data/nongenetic_PCs/AFA/pcs_most_variated.RDS") %>% as.matrix()
    cat("Loaded PCs with most variated 20000 CpGs and that cumulatively explained 80% of the variance.\n")
    cat("Loaded PCs have dimensions", dim(pcs_to_keep), "\n")
}
