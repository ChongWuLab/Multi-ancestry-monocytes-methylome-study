library(BEDMatrix)
library(data.table)
library(dplyr)
library(MASS)

# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

# ComputeAlpha
ComputeAlpha <- function(w, Z, n, n_0, LD.matrix) {
    # BIG SIGMA
    #SIGMA <- (t(ref %*% w) %*% (ref %*% w)) / n_0
    SIGMA <- t(w) %*% LD.matrix %*% w * (n_0-1)/n_0
    SIGMA <- SIGMA[1, 1]

    # Alpha hat
    alpha <- ((w %*% Z) / sqrt(n)) / SIGMA
    alpha <- alpha[1, 1]

    # small sigma
    sigma <- sqrt(1 - 2 * ((w %*% Z) / sqrt(n)) * alpha + alpha ^ 2 * SIGMA)
    sigma <- sigma[1, 1]

    # se of alpha
    se <- sqrt((1 / n + 1 / n_0) * (SIGMA * alpha ^ 2) / SIGMA + sigma ^ 2 / (n * SIGMA))

    return(c(alpha, se))
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
    snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
    snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}

# TScMLVar_alt
TScMLVar_alt <- function (Z.ref.original, Stage1FittedModel, betaalpha.hat.stage2, 
          Est.Sigma1Square, Est.Sigma2Square, n1, n2, n.ref) 
{
    scale.Z.ref = scale(Z.ref.original)
    set.A = which(Stage1FittedModel != 0)
    set.B = which(betaalpha.hat.stage2 != 0)[-1] - 1
    hat.gamma.A = Stage1FittedModel[set.A]
    U.AA = cov(scale.Z.ref[, set.A, drop = FALSE])
    if (length(set.B) == 1) {
        U.BB = diag(1)
    }
    else {
        U.BB = cov(scale.Z.ref[, set.B])
    }
    U.BA = cov(scale.Z.ref[, set.B], scale.Z.ref[, set.A, drop = FALSE])
    U.AB = t(U.BA)
    Est.Theta = solve(cov(scale.Z.ref[, set.A, drop = FALSE]) + diag(1e-05, 
                                                       length(set.A)))
    Est.beta = betaalpha.hat.stage2[1]
    Est.alpha = (betaalpha.hat.stage2[-1])[set.B]
    Est.Sigma = rbind(cbind(hat.gamma.A %*% U.AA %*% hat.gamma.A, 
                            hat.gamma.A %*% U.AB), cbind(U.BA %*% hat.gamma.A, U.BB))
    Est.Sigma = Est.Sigma + diag(1e-05, nrow(Est.Sigma))
    inv.Est.Sigma = solve(Est.Sigma)
    Matrix.A = (t(t(c(Est.beta, Est.alpha))) %*% t(c(Est.beta, 
                                                     Est.alpha)))
    Matrix.B = n2 * Est.Sigma2Square * Est.Sigma + Est.beta^2 * 
        Est.Sigma1Square/n1 * ((n2^2 + n2) * rbind(cbind(hat.gamma.A %*% 
                                                             U.AA %*% Est.Theta %*% U.AA %*% hat.gamma.A, hat.gamma.A %*% 
                                                             U.AA %*% Est.Theta %*% U.AB), cbind(U.BA %*% Est.Theta %*% 
                                                                                                     U.AA %*% hat.gamma.A, U.BA %*% Est.Theta %*% U.AB)) + 
                                   n2 * sum(diag(U.AA %*% Est.Theta)) * rbind(cbind(hat.gamma.A %*% 
                                                                                        U.AA %*% hat.gamma.A, hat.gamma.A %*% U.AB), cbind(U.BA %*% 
                                                                                                                                               hat.gamma.A, U.BB)))
    const.p = length(set.B) + 1
    constant.d = (n.ref - const.p) * (n.ref - const.p - 1) * 
        (n.ref - const.p - 3)
    Est.Cov.Mat = ((n.ref - const.p - 1)/constant.d * (n2 + n2^2) * 
                       Matrix.A + n2 * (n2 + n.ref)/constant.d * sum(diag(Matrix.A %*% 
                                                                              Est.Sigma)) * inv.Est.Sigma) * n.ref^2/n2^2 + n.ref^2/n2^2 * 
        (1/constant.d * sum(diag(inv.Est.Sigma %*% Matrix.B)) * 
             inv.Est.Sigma + (n.ref - const.p - 1)/constant.d * 
             (inv.Est.Sigma %*% Matrix.B %*% inv.Est.Sigma)) - 
        n.ref^2/(n.ref - const.p - 1)^2 * Matrix.A
    return(Est.Cov.Mat[1, 1])
}

#################
# Main function #
#################

# TestAssociation
TestAssociation <- function(sumstats, 
                            n.sumstats,
                            ref, 
                            n_sample, 
                            snps, 
                            LD) {
    # Exception: SS.original contains SNPs from multiple chromosomes
    
    # Get data
    sumstats.temp <- sumstats
    snps.temp <- snps
    # Find common SNPs
    snps.common <- intersect(sumstats.temp$SNP, snps.temp$SNP)
    if(length(snps.common) <= 1) {
    stop("No overlapping SNPs found for the provided datasets.")
    } 
    # Keep only common SNPs
    sumstats.temp <- sumstats.temp[snps.common, ]
    snps.temp  <- snps.temp[snps.common, ]
    #weights        <- snps.temp$weight
    ref <- ref[snps.common, ]
    
    #reference.bed <- reference.bed[, snps.common] %>% as.matrix() %>% PatchUp() %>% scale()
    LD = LD[snps.common,snps.common]
    LD = as.matrix(LD)
    # Allele-flip the phenotype ss and weight w.r.t. reference panel
    qc.1 <- allele.qc(
        sumstats.temp$A2,
        sumstats.temp$A1,
        ref$A1,
        ref$A2
    )

    sumstats.temp$Z[qc.1$flip] <- -1 * sumstats.temp$Z[qc.1$flip]

    # Remove strand ambiguous SNPs (if any)
    if ( sum(!qc.1$keep) > 0 ) {
        sumstats.temp = sumstats.temp[qc.1$keep,]
        snps.temp = snps.temp[qc.1$keep,]
        ref = ref[qc.1$keep,]
    }
    
    
    qc.2 <- allele.qc(
        snps.temp$A0,
        snps.temp$A1,
        ref$A1,
        ref$A2
    )

    snps.temp[qc.2$flip,"weight"] <- -1 *  snps.temp[qc.2$flip,"weight"]

    if (!("N" %in% colnames(sumstats.temp))) {
        sumstats.temp["N"] <- n.sumstats
    } else {
        n.sumstats = floor(mean(sumstats.temp[,"N"]))
    }
    
    sumstats.temp["correlation"] <- sumstats.temp$Z / sqrt(sumstats.temp$Z ^ 2 + sumstats.temp$N - 1)

    # Remove strand ambiguous SNPs (if any)
    keep <- qc.1$keep & qc.2$keep
    if (sum(!keep) > 0) {
        sumstats.temp <- sumstats.temp[keep, ]
        snps.temp   <- snps.temp[keep, ]
        #weights        <- weights[keep]
        ref <- ref[keep, ]
        #reference.bed <- reference.bed[, keep]
        LD <- LD[keep,keep]
    }

    ################
    # Classic PWAS #
    ################

    # Keep the non-zero components of weights vector
    weights.used <- as.vector(snps.temp$weight)

    # make weights.used to a vector
    nusedsnp = length(weights.used)
    
    # Compute TWAS z-score, r2, and p-value
    z.twas  <- as.numeric(weights.used %*% sumstats.temp$Z)
    #matrix.LD <- t(reference.bed) %*% reference.bed / (nrow(reference.bed) - 1)
    r2.twas <- as.numeric(weights.used %*% LD %*% weights.used)

    p.classic <- 2 * (pnorm(abs(z.twas / sqrt(r2.twas)), lower.tail = FALSE))

    ####################
    # Alternative PWAS #
    ####################
    
    if (length(weights.used) > 1) {
        out <- ComputeAlpha(
            w   = weights.used,
            Z   = sumstats.temp$Z,
            n   = n.sumstats,
            n_0 = n_sample,
            LD.matrix = LD
        )

        out[3] <- 2 * (pnorm(abs(out[1] / out[2]), lower.tail = FALSE))
    } else {
        out <- c(NA, NA, NA)
    }

    output <- c(z.twas / sqrt(r2.twas), p.classic, out, nusedsnp)
    names(output) <- c(
        "Zscore.classic", "p.classic",
        "beta_alt", "se_alt", "p_alt","n_used_snp"
    )

    return(output)
}