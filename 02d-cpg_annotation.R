dat.wgbs.dir = "/rsrch5/home/biostatistics/wzhang24/data/WGBS/normalized_AFA/"

files = list.files(dat.wgbs.dir, full.names = TRUE)

library(GenomicRanges)
library(annotatr)
library(parallel)
library(data.table)

annotations <- build_annotations(genome= "hg38", annotations = c("hg38_cpg_islands"))
get_annotation = function(file){
    DNAm_data = readRDS(file)
    cpg_sites = DNAm_data[,c(1,2,4,5)]
    colnames(cpg_sites) = c("CpG","chr","pos38","pos38.right")
    cpg_sites$CpG = paste0("CpG",cpg_sites$CpG)
    gr_cpg_sites <- GRanges(seqnames = cpg_sites[,2], IRanges(start = cpg_sites[,3], end = cpg_sites[,4]))
    #Create Shores
    shores <- suppressWarnings(flank(annotations, width = 2000, both = TRUE))
    shores <- setdiff(shores, annotations)

    #Create shelves
    shelves <- suppressWarnings(flank(shores, width = 2000, both = TRUE))
    shelves <- setdiff(shelves, c(annotations, shores))

    #Create open seas
    open_seas <- GRanges(seqnames = cpg_sites[,1], IRanges(start = 1, end = 3000))
    open_seas <- setdiff(open_seas, c(annotations, shores, shelves))

    # Initialize the annotation column
    cpg_sites$annotation <- "open_sea"

    # Annotate based on overlaps
    overlaps_islands <- findOverlaps(gr_cpg_sites, annotations)
    cpg_sites$annotation[queryHits(overlaps_islands)] <- "island"

    overlap_shores <- findOverlaps(gr_cpg_sites, shores)
    cpg_sites$annotation[queryHits(overlap_shores)] <- "shore"

    overlap_shelves <- findOverlaps(gr_cpg_sites, shelves)
    cpg_sites$annotation[queryHits(overlap_shelves)] <- "shelf"
    return(cpg_sites)
}

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

tables <- mclapply(files, get_annotation, mc.cores = no_cores)
res_table <- rbindlist(tables)

saveRDS(res_table, file = "/rsrch5/home/biostatistics/wzhang24/data/WGBS/annotation.RDS")
