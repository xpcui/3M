require(glmnet)
require(dplyr)

#########################################################################################
# INTERACTION TREES
#########################################################################################

IT <- function(Y, Tr, G, ntree=50){
  it.data <- cbind.data.frame(Y, Tr, G)
  
  SNP.ID <- colnames(G); colnames(it.data) <- c("Response","Treatment",SNP.ID)
  
  ctg <- NA; col.y <- 1; col.Tr <- 2; split.var <- 3:ncol(it.data)
  
  fit <- Build.RF.IT(dat=it.data, col.y=col.y, col.trt=col.Tr, split.var=split.var,
                     ctg=ctg, N0=20, n0=5, max.depth=20, ntree=ntree, mtry = max(floor(length(split.var)/5), 1), 
                     avoid.nul.tree=T)
  
  # VARIABLE IMPORTANCE
  VI <- Variable.Importance(fit, n0=5, sort=F, details=F, truncate.zeros=T)
  
  VI <- VI %>% abs() %>% round(3)
  
  return(VI)
}

#########################################################################################
# LASSO
#########################################################################################

Lasso <- function(Y, Tr, G){
  fit <- lm(Y ~ Tr)
  resid = fit$residuals
  
  n <- dim(G)[1]; m <- dim(G)[2]; SNP.ID <- colnames(G)
  X <- cbind(G, Tr*G)
  
  crossval <-  cv.glmnet(x = X, y = resid)
  penalty <- crossval$lambda.min
  fit <- glmnet(x = X, y = resid, alpha = 1, lambda = penalty, intercept = FALSE)
  
  beta.vec <- coef(fit) %>% as.vector(); beta.vec <- beta.vec[-1]
  
  beta_G <- beta.vec[1:m]; beta_GT <- beta.vec[(m+1):(2*m)]
  
  names(beta_GT) <- SNP.ID
  
  beta_GT <- beta_GT %>% abs() %>% round(3)
  
  return(beta_GT)
}

#########################################################################################
# P-VALUE RANKING
#########################################################################################

pvalue.rank <- function(Y, Tr, G){
  n <- dim(G)[1]; m <- dim(G)[2]; SNP.ID <- colnames(G)
  
  pvalue.est <- double()
  
  for(i in 1:m){
    fit <- summary(lm(Y ~ Tr*G[,i]))
    pvalue.est[i] <- coef(fit)[4,4]
  }
  
  names(pvalue.est) <- SNP.ID
  
  # -LOG10 P-VALUE
  re <- -log10(pvalue.est)
  return(re)
}

