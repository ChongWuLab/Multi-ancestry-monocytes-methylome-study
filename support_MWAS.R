
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

PatchUp <- function(M) {
  for (p in 1:ncol(M)) {
    # Get the current column vector
    vec <- M[, p]
    
    # Skip if no NA is found
    if (sum(is.na(vec)) == 0) {
      next()
    }

    # Impute by averaging
    M[is.na(vec), p] <- mean(vec, na.rm = TRUE)
  }
  
  return(M)
}

# Elastic Net
weights.enet = function( genos , pheno , alpha=0.5 ) {
    set.seed(1)
    eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
    # remove monomorphics
    sds = apply( genos  , 2 , sd )
    keep = sds != 0 & !is.na(sds)
    enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
    eff.wgt[ keep ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
    return( eff.wgt )
}

# Define ComputeCoef()
ComputeCoef <- function(Z, Y, libraryNames, verbose, obsWeights, ...) {
    cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x - Y) ^ 2))
    names(cvRisk) <- libraryNames
    fit.nnls <- nnls(sqrt(obsWeights) * Z, sqrt(obsWeights) * Y)
    if (verbose) {
        message(paste("Non-Negative least squares convergence:", fit.nnls$mode == 1))
    }
    initCoef <- coef(fit.nnls)
    initCoef[is.na(initCoef)] <- 0
    if (sum(initCoef) > 0) {
        coef <- initCoef / sum(initCoef)
    } else {
        warning("All algorithms have zero weight", call. = FALSE)
        coef <- initCoef
    }
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = fit.nnls)
    return(out)
}

weights.stdPWAS <- function(genos, pheno, alpha = 0.5,crossval = 5) {
    
    set.seed(1)
    
    M = 1
    models = "elastic_net"
    cv.performance = matrix(NA,nrow=4,ncol=M)
    rownames(cv.performance) = c("rsq","pval", "correlation", "cor_pval")
    colnames(cv.performance) = models
    
    cv.all = pheno
    #if we scale the phenotype, we need
    #cv.all = scale(pheno)
    N = nrow(cv.all)
    cv.sample = sample(N)
    cv.all = cv.all[ cv.sample , ,drop=FALSE]
    folds = cut(seq(1,N),breaks=crossval,labels=FALSE)

    cv.calls = matrix(NA,nrow=N,ncol=1)

    for (i in 1:crossval) {
        indx = which(folds==i,arr.ind=TRUE)
        cv.train = cv.all[-indx,]
        # store intercept
        #intercept = mean( cv.train[,3] )
        #cv.train[,3] = scale(cv.train[,3]) # not sure if we really need this?
        
        # hide current fold
        genostmp = genos[cv.sample[ -indx ],]
        
        pred.wgt = weights.enet( genostmp , as.matrix(cv.train) , alpha=0.5 )

        # predict from weights into sample
        pred.wgt[ is.na(pred.wgt) ] = 0
        cv.calls[ indx , 1 ] =  genos[ cv.sample[ indx ] , ] %*% pred.wgt
    }
    
    # compute rsq + P-value for each model
    mod = 1
    if ( !is.na(sd(cv.calls[,mod])) && sd(cv.calls[,mod]) != 0 ) {
            reg = summary(lm( cv.all ~ cv.calls[,mod] ))
            cor.res <- cor.test(cv.all, cv.calls[,mod])
            cv.performance[ 1, mod ] = reg$adj.r.sq
            cv.performance[ 2, mod ] = reg$coef[2,4]
            cv.performance[ 3, mod ] = cor.res$estimate
            cv.performance[ 4, mod ] = cor.res$p.value
    } else {
            cv.performance[ 1, mod ] = NA
            cv.performance[ 2, mod ] = NA
            cv.performance[ 3, mod ] = NA
            cv.performance[ 4, mod ] = NA
    }
    
    wgt.matrix = matrix(0,nrow=ncol(genos),ncol=M)
    colnames(wgt.matrix) = models
    rownames(wgt.matrix) = colnames(genos)

    wgt.matrix[,mod] = weights.enet( genos, as.matrix(pheno) , alpha=0.5 )

    return(list(cv.performance = cv.performance, wgt.matrix = wgt.matrix))
}
