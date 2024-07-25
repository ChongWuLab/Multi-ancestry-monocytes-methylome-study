library(BEDMatrix)
library(data.table)
library(dplyr)
#library(lasso2)
library(MASS)
#library(TScML)

# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

# allele.qc
allele.qc <- function(a1, a2, ref1, ref2) {
    ref <- ref1
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip1 <- flip
    ref <- ref2
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip2 <- flip
    snp <- list()
    snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C") | (a1 == "I" | a2 == "D") | (a1 == "D" | a2 == "I"))
    snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}



clumping.iv <- function(SS.original, r2.threshold = 0.001) {
    ss = SS.original
    chromosome = chr.id
    # Create a temporary SNP list for Plink
    temp.snpfile <- paste0(temp.outdir, cur.cpg, ".txt")
    write.table(
        ss$SNP,
        file = temp.snpfile,
        quote = FALSE,
        col.names = FALSE,
        row.names = FALSE
    )

    # Put together the command and execute
    command <- "module load plink/1.90-beta && plink"
    command <- paste0(command, " --bfile ",ldref, chromosome)
    command <- paste0(command, " --chr ", chromosome)
    command <- paste0(command, " --extract ", temp.snpfile)
    command <- paste0(command, " --make-bed")
    command <- paste0(command, " --out ", temp.outdir, cur.cpg)
    
    system(command, intern = TRUE)
    

    # clumping
    # Write ss
    write.table(
        ss,
        file = temp.snpfile,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )
    
    command <- "module load plink/1.90-beta && plink"
    command <- paste0(command, " --bfile ", temp.outdir, cur.cpg)
    command <- paste0(command, " --clump ", temp.snpfile)
    command <- paste0(command, " --clump-p1 1")
    command <- paste0(command, " --clump-p2 1")
    command <- paste0(command, " --clump-r2 ",r2.threshold)
    command <- paste0(command, " --clump-kb 10000")
    command <- paste0(command, " --out ", temp.outdir, cur.cpg)

    
    system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    
    if (paste0(temp.outdir, cur.cpg, ".clumped") %>% file.exists()) {
        clump.keep <- paste0(temp.outdir, cur.cpg, ".clumped") %>% fread(., data.table = FALSE, header = TRUE)
    } else {
        clump.keep <- NULL
        return(clump.keep)
    }

    if (ncol(clump.keep) == 1) {
        return(clump.keep[, 1])
    } else {
       return(clump.keep$SNP)
    }
}