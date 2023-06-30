library(STAREG)
source("./simulation/data_generation.R")

n = 100 # number of replicates
m = 10000
xi00 = 0.9
xi01 = 0.025
xi10 = 0.025
xi11 = 0.05
mu1 = 2
mu2 = 2
sigma1 = 1
sigma2 = 1
alphas <- seq(0.01, 0.1, 0.01)

fdp_naiveBH  <- matrix(NA, n, length(alphas))
fdp_maxP     <- matrix(NA, n, length(alphas))
fdp_STAREG    <- matrix(NA, n, length(alphas))

pd_naiveBH  <- matrix(NA, n, length(alphas))
pd_maxP     <- matrix(NA, n, length(alphas))
pd_STAREG    <- matrix(NA, n, length(alphas))

xi00.hat <- rep(NA, n)
xi01.hat <- rep(NA, n)
xi10.hat <- rep(NA, n)
xi11.hat <- rep(NA, n)

for(i in 1:n){
  cat('i = ', i, '\n')

  data.obj <- data_generation(m, xi00, xi01, xi10, xi11, mu1, mu2, sigma1, sigma2)
  p1 = data.obj$pvals1
  p2 = data.obj$pvals2
  states1 = data.obj$states1
  states2 = data.obj$states2
  
  # BH
  pval1.bh <- p.adjust(p1, method = "BH")
  pval2.bh <- p.adjust(p2, method = "BH")
  
  # MaxP
  maxp <- apply(cbind(p1, p2), 1, max)
  pval.maxp <- p.adjust(maxp, method = "BH")
  
  # STAREG
  res.rep <- stareg(p1, p2)
  pval.rep <- res.rep$fdr
  xi00.hat[i] <- res.rep$xi00
  xi01.hat[i] <- res.rep$xi01
  xi10.hat[i] <- res.rep$xi10
  xi11.hat[i] <- res.rep$xi11
  
  for(j in 1:length(alphas)){
    alpha = alphas[j]
    
    # naive BH
    fdp_naiveBH[i,j] <- sum(pval1.bh <= alpha & pval2.bh <= alpha & !(states1 * states2))/ max(sum(pval1.bh <= alpha & pval2.bh <= alpha), 1)
    pd_naiveBH[i,j]  <- sum(pval1.bh <= alpha & pval2.bh <= alpha & (states1 * states2)) / sum(states1 * states2)
    
    # MaxP
    # thr <- maxP(p1, p2, alpha)
    fdp_maxP[i,j] <- sum(pval.maxp <= alpha & !(states1 * states2))/ max(sum(pval.maxp <= alpha), 1)
    pd_maxP[i,j]  <- sum(pval.maxp <= alpha & (states1 * states2)) / sum(states1 * states2)
    
    # STAREG
    fdp_STAREG[i,j] <- sum(pval.rep <= alpha & !(states1 * states2))/ max(sum(pval.rep <= alpha), 1)
    pd_STAREG[i,j] <- sum(pval.rep <= alpha & (states1 * states2)) / sum(states1 * states2)
  }
}

fdr_BH <- colMeans(fdp_naiveBH, na.rm = TRUE) 
fdr_MaxP <- colMeans(fdp_maxP, na.rm = TRUE) 
fdr_STAREG <- colMeans(fdp_STAREG, na.rm = TRUE)

power_BH <- colMeans(pd_naiveBH, na.rm = TRUE) 
power_MaxP <- colMeans(pd_maxP, na.rm = TRUE)
power_STAREG <- colMeans(pd_STAREG, na.rm = TRUE)  

