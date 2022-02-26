# Code for identifying HVG
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7017299/?fbclid=IwAR2lSsdnsVOsj_vuRqqhWO4pPVpDK1gPVwLoBkzdHotmEEo-2biTWjqXgjM
# Code adopted from https://github.com/cailab-tamu/HVG/tree/master/R

# Find HVG
findHVG <- function(X, cutOff = 0.01){
  nCells <- ncol(X)
  cellBarcodes <- colnames(X)
  
  temp <- X[rowMeans(X!=0) > 0.05,]
  means <- rowMeans(temp, na.rm = T)
  vars <- apply(temp, 1, var, na.rm=T)
  cv2 <- vars/(means^2)
  recip.means <- 1/means
  recip.means[is.infinite(recip.means)] <- 0
  fit <- glm(cv2~recip.means, family = Gamma(link = 'identity'))
  pFit <- predict(fit)
  pVal <- pchisq((cv2/pFit)*(nCells-1),(nCells-1), lower.tail = FALSE)
  FC <- log2(cv2/pFit)
  pAdj <- p.adjust(pVal, method = 'fdr')
  hvgStat <- cbind(means, cv2, pFit,FC, pVal, pAdj)
  hvgStat <- as.data.frame(hvgStat)
  colnames(hvgStat) <- c('mean', 'cv2obs', 'cv2exp','log2FC', 'p.value', 'p.adj')
  HVG <- names(pAdj[pAdj < cutOff & FC > log2(1.5)])
  length(HVG)
  out <- list()
  out$HVG <- HVG
  out$stat <- hvgStat
  return(out)
}


