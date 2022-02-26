# running logistic regression model to check for association between accessibility and mode membership 
# calculate parameters of the logistic regression model
# Functions for running logistic regression model
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



diff_test_helper <- function(gexpr, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily
){ 
  
  reducedModelFormulaStr <- paste("accessibility", reducedModelFormulaStr, sep="~")
  fullModelFormulaStr <- paste("accessibility", fullModelFormulaStr, sep="~")
  
  test_res <- tryCatch({
    
    full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily, data = gexpr)
    reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily, data = gexpr)
    
    compareModels(list(full_model_fit), list(reduced_model_fit))
  }, 
  
  #warning = function(w) { FM_fit },
  error = function(e) { 
    print (e);
    data.frame(status = "FAIL", family=expressionFamily@vfamily, pval=1.0, qval=1.0)
  }
  )
  test_res
}


reducedModelFormulaStr <- 'beta'
fullModelFormulaStr <- 'alpha + beta'
expressionFamily <- "binomialff"

diff_lrt <- NULL


for (i in 1:nrow(snare_binary_2orm)) {
  
  gexpr <- data.frame(accessibility = as.factor(snare_binary_2orm[i,]),
                      alpha = as.factor(DP_cluster_2orm[i,]),
                      beta = beta)
  
  diff_lrt_01 <- diff_test_helper(gexpr, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily)
  
  diff_lrt <- rbind(diff_lrt, diff_lrt_01) 
}

rownames(diff_lrt) <- rownames(snare_binary_2orm)
save(diff_lrt, file = "./RObjects/diff_lrt.rds")

#select sites with q-val < 0.01
diff_lrt_sig <- diff_lrt[diff_lrt$qval < 0.01,]


