# NB-GLM for association between gene expression and mode membership

# Estimate the size factors using DESeq2
library(DESeq2)
coldata <- data.frame(condition = as.factor(rep(1, ncol(rna_mat))))
rownames(coldata) <- colnames(rna_mat)

dds <- DESeqDataSetFromMatrix(countData = rna_mat,
                              colData = coldata,
                              design= ~ 1)
dds <- estimateSizeFactors(dds)
lib.size <- sizeFactors(dds)

##################
# VGLM-NB to test whether the ATAC-seq modes are associated with gene expression
compareModels <- function(full_models, reduced_models){
  
  stopifnot(length(full_models) == length(reduced_models))
  
  test_res <- mapply(function(x,y) { 
    if (is.null(x) == FALSE && is.null(y) == FALSE) {
      lrt <- VGAM::lrtest(x,y) 
      pval=lrt@Body["Pr(>Chisq)"][2,]
      family = x@family@vfamily
      if (length(family) > 1)
        family = family[1]
      data.frame(status = "OK", family=family, pval=pval)
    } else { data.frame(status = "FAIL", family=NA, pval=1.0) } 
  } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
  
  test_res <- do.call(rbind.data.frame, test_res)
  test_res$qval <- p.adjust(test_res$pval, method="BH")
  test_res
}



diff_test_NB <- function(gexpr, 
                         fullModelFormulaStr, 
                         reducedModelFormulaStr
){ 
  
  reducedModelFormulaStr <- paste("genexpr", reducedModelFormulaStr, sep="~")
  fullModelFormulaStr <- paste("genexpr", fullModelFormulaStr, sep="~")
  
  test_res <- tryCatch({
    
    full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), negbinomial(lmu = "loglink", imethod = 1), data = gexpr, offset = log(lib.size))
    reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), negbinomial(lmu = "loglink", imethod = 1), data = gexpr, offset = log(lib.size))
    
    compareModels(list(full_model_fit), list(reduced_model_fit))
  }, 
  
  #warning = function(w) { FM_fit },
  error = function(e) { 
    print (e);
    data.frame(status = "FAIL", family="negbinomial", pval=1.0, qval=1.0)
  }
  )
  test_res
}


reducedModelFormulaStr <- '1'
fullModelFormulaStr <- 'beta'

library(VGAM)
diff_lrt_nb <- NULL


for (i in 1:nrow(rna_2ormore)) {
  
  gexpr <- data.frame(genexpr = rna_2ormore[i,],
                      beta = as.factor(DP_cluster_rna[i,]),
                      lib.size = lib.size)
  
  diff_lrt_01 <- diff_test_NB(gexpr, fullModelFormulaStr, reducedModelFormulaStr)
  
  diff_lrt_nb <- rbind(diff_lrt_nb, diff_lrt_01) 
}

rownames(diff_lrt_nb) <- rownames(rna_2ormore)
save(diff_lrt_nb, file = "./RObjects/diff_lrt_nb.rds")

#select sites with q-val < 0.01
diff_nb_sig <- diff_lrt_nb[diff_lrt_nb$qval < 0.01,]

