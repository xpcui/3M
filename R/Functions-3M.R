require(dplyr)

# ==============================================================================
# LSmeans (Estimation) + MCB (Selection)
# ==============================================================================

# ------------------------------------------------------------------------------
# Main functions
# ------------------------------------------------------------------------------

## LSmeans + MCB
LSmeans_MCB <- function(Y, Tr, G, alpha = 0.05, nboot = 100, length.out = 10){
  ID <- colnames(G); nsubj <- dim(G)[1]; nsnp <- dim(G)[2]
  
  alphabeta_mat <- cal_alphabeta_mat(Y, Tr, G)
  
  alphabeta_star <- apply(alphabeta_mat, 1, min, na.rm = TRUE)
  
  alphabeta_mcb <- sapply(1:nsnp, convert_MCB, alphabeta_star)
  
  alphabeta_boot <- sapply(1:nboot, generate_boot, Y, Tr, G) %>% t() # Time cost in this step
  
  width <- sapply(1:nsnp, calculate_delta, alphabeta_star, alphabeta_boot, nboot, length.out, alpha)
  
  CIs <- sapply(1:nsnp, construct_CI, alphabeta_star, width) %>% t()
  
  re0 <- cbind.data.frame(ID, CIs[,1], alphabeta_mcb, CIs[,2], width/2, alphabeta_star)
  colnames(re0) <- c("SNPID","lowerbound","mcb","upperbound","width","beta")
  rownames(re0) <- NULL
  return(re0)
}

## LSmeans
LSmeans_MCB_est <- function(Y, Tr, G){
  ID <- colnames(G); nsubj <- dim(G)[1]; nsnp <- dim(G)[2]
  
  alphabeta_mat <- cal_alphabeta_mat(Y, Tr, G)
  
  alphabeta_star <- apply(alphabeta_mat, 1, min, na.rm = TRUE)
  
  names(alphabeta_star) <- ID
  alphabeta_star
}

# ------------------------------------------------------------------------------
# Supplementary functions
# ------------------------------------------------------------------------------

## Calculate conditional estimates matrix
cal_alphabeta_mat <- function(Y, Tr, G){
  nsnp <- dim(G)[2]
  
  alphabeta_mat <- matrix(NA, ncol = nsnp, nrow = nsnp)
  for (i in 1:(nsnp-1)) {
    for(k in (i+1):nsnp){
      fit <- lm(Y ~ Tr + G[,i] + G[,k] + Tr:G[,i] + Tr:G[,k]) # Tr:G[,i]
      alphabeta_mat[i, k] <- ifelse(is.na(fit$coefficients[5]), 0, fit$coefficients[5]) 
      alphabeta_mat[k, i] <- ifelse(is.na(fit$coefficients[6]), 0, fit$coefficients[6]) 
    }
  }
  
  alphabeta_mat
}

## Bootstrap
generate_boot <- function(nitr, Y, Tr, G){
  nsnp <- ncol(G); nsubj <- nrow(G)
  
  index.boot <- sample(1:nsubj, nsubj, replace = TRUE) %>% sort()
  Y.boot <- Y[index.boot]
  Tr.boot <- Tr[index.boot]
  G.boot <- G[index.boot,]
  
  alphabeta_mat_boot <- cal_alphabeta_mat(Y.boot, Tr.boot, G.boot)
  
  alphabeta_star_boot <- apply(alphabeta_mat_boot, 1, min, na.rm = TRUE)
  return(alphabeta_star_boot)
}

## Calculate widths of confidence intervals based on bootstrap results
calculate_delta <- function(i, alphabeta_star, alphabeta_boot, nboot, length.out, alpha){
  obs <- alphabeta_star[i] - alphabeta_star[-i] # MCC
  Delta <- double()
  for (j in 1:nboot) {
    boot <- alphabeta_boot[j,i] - alphabeta_boot[j,-i] # MCC boot
    Delta <- rbind(Delta, boot - obs)
  }
  gap <- seq(0, max(Delta), length.out = length.out)
  Count <- double()
  for (p in 1:length(gap)) {
    gap_i <- gap[p]
    count <- sapply(1:nboot, in_and_out, Delta, gap_i)
    Count <- c(Count, sum(count))
  }
  index <- which(Count >= nboot*(1-alpha))[1]
  return(gap[index])
}

in_and_out <- function(k, Delta, gap_i){
  vec <- Delta[k,]
  re <- sum(vec <= gap_i) == length(vec)
  return(re)
}

## Calculate confidence intervals
construct_CI <- function(i, alphabeta_star, width){
  di <- width[i]/2
  mcb <- alphabeta_star[i] - max(alphabeta_star[-i])
  Di.minus <- min(mcb-di, 0)
  Di.plus <- max(mcb+di, 0)
  return(c(Di.minus, Di.plus))
}

## Calculate MCB parameters
convert_MCB <- function(i, alphabeta_star){
  re <- alphabeta_star[i] - max(alphabeta_star[-i])
  return(re)
}





