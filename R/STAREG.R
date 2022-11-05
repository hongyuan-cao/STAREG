############################################################################################
# Package: STAREG
# Version: 0.1.0
# Data: 2022-09-29
############################################################################################
#' @title Replicability analysis of two multiple testing studies providing FDR control
#'
#' @param pvals1 a numeric vector of p-values from study 1.
#' @param pvals2 a numeric vector of p-values from study 2.
#' @param minIter an integer value indicating the minimum number of iterations.
#' @param maxIter an integer value indicating the maximum number of iterations.
#' @param est.pi0 a logical value indicating whether to estimate the proportion of true nulls in each study for the initialization of xi's or just specify the value as 0.9. Default is FALSE.
#' @param pi0.const a numeric number between 0 and 1 to initialize xi's if est.pi0 is false. 
#' @param tol a numeric value giving the tolerance in the relative change in the log likelihood below which the EM algorithm is considered to be converged. Default is set to 1e-3.
#' @param trace a logical value indicating whether to print out the EM process.
#'
#' @return a list with the following elements:
#' \item{call}{the function call made.}
#' \item{lfdr.rep}{a numeric vector of the estimated values of joint local FDR.}
#' \item{fdr.rep}{a numeric vector of the adjusted p-values for FDR control based on the joint local FDR.}
#' \item{xi00}{a numeric value of the estimated probability for the joint hidden state of (0, 0)}
#' \item{xi01}{a numeric value of the estimated probability for the joint hidden state of (0, 1)}
#' \item{xi10}{a numeric value of the estimated probability for the joint hidden state of (1, 0)}
#' \item{xi11}{a numeric value of the estimated probability for the joint hidden state of (1, 1)}
#' \item{f1}{a vector of the estimated densities for the p-value under the alternatives in study 1}
#' \item{f2}{a vector of the estimated densities for the p-value under the alternatives in study 1}
#' \item{loglik}{a numeric vector of the log likelihood at each step}
#' \item{converge}{a list containing the convergence information. \code{code}: 1 - converged, 0 - not converged; \code{iter}: the number of iterations performed.}
#'
#' @authors Yan Li, Xianyang Zhang, Run Chen, Zhi Wei, Xianyang Zhang and Hongyuan Cao.
#' @reference Replicability analysis of spatial variable gene detection across multiple samples in spatially resolved transcriptomic data.
#' @importFrom qvalue pi0est
#' @importFrom Iso pava
#'
#' @example
#' ## obtained two sequence of p-values from two SVG detection studies: pvals1, pvals2
#' data.obj <- data_generation(m = 10000, xi00 = 0.9, xi01 = 0.025, xi10 = 0.025, xi11 = 0.05, mu1 = 2, mu2 = 2, sigma1 = 1, sigma2 = 1)
#' pvals1 = data.obj$pvals1
#' pvals1 = data.obj$pvals2
#'
#' ## replicability analysis
#' stareg.obj <- STAREG(pvals1, pvals2)
#' pval.rep <- stareg.obj$fdr.rep
#'
#' ## calculate the false discovery proportion and true positive proportion
#' states1 = data.obj$states1
#' states2 = data.obj$states2
#' fdp = sum(pval.rep <= 0.05 & !(states1 * states2))/ max(sum(pval.rep <= 0.05), 1)
#' tpp = sum(pval.rep <= 0.05 & (states1 * states2)) / sum(states1 * states2)
#'
#' @export


STAREG <- function(pvals1, pvals2, minIter = 3, maxIter = 100, est.pi0 = FALSE, pi0.const = 0.9, tol = 1e-3, trace = FALSE) {

  k.init = 0.75
  pvals.cutoff = 1e-15
  f1.cutoff = 1e-15

  if(est.pi0){
    pi0.pvals1 <- min(pi0est(pvals1, pi0.method = 'smoother')$pi0, 0.999)
    pi0.pvals2 <- min(pi0est(pvals2, pi0.method = 'smoother')$pi0, 0.999)
  }else{
    pi0.pvals1 = pi0.const
    pi0.pvals2 = pi0.const
  }

  m <- length(pvals1)
  m1 <- length(pvals2)
  if(m!=m1){
    stop("p-value vectors must be of equal lengths!")
  }
  pvals1[pvals1 == 0] <- min(min(pvals1[pvals1 != 0]), pvals.cutoff)
  pvals2[pvals2 == 0] <- min(min(pvals2[pvals2 != 0]), pvals.cutoff)

  out1 <- sort(pvals1, index.return = TRUE)
  p1 <- out1$x
  ix1 <- out1$ix
  p1.diff <- c(p1[1], diff(p1))

  out2 <- sort(pvals2, index.return = TRUE)
  p2 <- out2$x
  ix2 <- out2$ix
  p2.diff <- c(p2[1], diff(p2))

  # Initialization
  xi00 <- pi0.pvals1 * pi0.pvals2
  xi01 <- pi0.pvals1 * (1 - pi0.pvals2)
  xi10 <- (1 - pi0.pvals1) * pi0.pvals2
  xi11 <- (1 - pi0.pvals1) * (1 - pi0.pvals2)

  f0 <- 1
  f1 <- (1 - k.init) * pvals1 ^ (-k.init)
  f2 <- (1 - k.init) * pvals2 ^ (-k.init)

  iter <- 0
  loglik <- -Inf
  converge <- list()

  if (trace)
    cat('STAREG: Begin iteration ...\n')

  while (TRUE) {
    iter <- iter + 1

    # E-step
    f <- xi00 * f0 * f0 + xi10 * f1 * f0 + xi01 * f0 * f2 + xi11 * f1 * f2
    gamma00 <- xi00 * f0 * f0 / f
    gamma01 <- xi01 * f0 * f2 / f
    gamma10 <- xi10 * f1 * f0 / f
    gamma11 <- 1 - gamma00 - gamma01 - gamma10

    # M-step
    # Update the f1 distribution
    Q1 <- gamma01 + gamma00
    Q2 <- gamma10 + gamma00
    Q1 <- Q1[ix1]
    Q2 <- Q2[ix2]

    y1 <- -p1.diff * sum(1 - Q1, na.rm = TRUE) / (1 - Q1)
    y2 <- -p2.diff * sum(1 - Q2, na.rm = TRUE) / (1 - Q2)

    y1[is.na(y1)] <- -Inf
    y1[y1 == 0] <- max(y1[y1 != 0])
    y1[y1 == -Inf] <- min(y1[y1 != -Inf])

    y2[is.na(y2)] <- -Inf
    y2[y2 == 0] <- max(y2[y2 != 0])
    y2[y2 == -Inf] <- min(y2[y2 != -Inf])

    res1 <- pava(y1, 1 - Q1, decreasing = TRUE, long.out = FALSE, stepfun = FALSE)
    res2 <- pava(y2, 1 - Q2, decreasing = TRUE, long.out = FALSE, stepfun = FALSE)

    f1 <- -1 / res1
    f1 <-  f1 / sum(f1 * p1.diff, na.rm = TRUE)
    f1[ix1] <- f1
    f1[which(is.na(f1))] <- min(f1[!is.na(f1)])
    f1[f1 <= 0] <- min(min(f1[f1 > 0]), f1.cutoff)

    f2 <- -1 / res2
    f2 <- f2 / sum(f2 * p2.diff, na.rm = TRUE)
    f2[ix2] <- f2
    f2[which(is.na(f2))] <- min(f2[!is.na(f2)])
    f2[f2 <= 0] <- min(min(f2[f2 > 0]), f1.cutoff)

    # Update the xi's
    xi00 <- mean(gamma00)
    xi01 <- mean(gamma01)
    xi10 <- mean(gamma10)
    xi11 <- mean(gamma11)

    loglik0 <-
      sum((gamma10 + gamma11) * log(f1) + (gamma01 + gamma11) * log(f2)) +
      sum(gamma11 * log(xi11) + gamma10 * log(xi10) + gamma01 * log(xi01) + gamma00 *
            log(xi00))
    loglik <- c(loglik, loglik0)
    loglik.delta <- abs((loglik0 - loglik[iter]) / loglik[iter])

    if (trace) cat('Iter ', iter, ': loglik = ', loglik0, ', eps = ', loglik.delta, '\n', sep = '')
    if (iter >= minIter & loglik.delta <= tol & iter <= maxIter) {
      if (trace) cat('Converged!\n')
      converge$code <- 1
      converge$iter <- iter
      break
    }
    if (iter > maxIter) {
      if (trace) cat('Maximum iteration reached!\n')
      converge$code <- 0
      converge$iter <- iter
      break
    }
  }

  f <- xi00 * f0 * f0 + xi10 * f1 * f0 + xi01 * f0 * f2 + xi11 * f1 * f2
  gamma00 <- xi00 * f0 * f0 / f
  gamma01 <- xi01 * f0 * f2 / f
  gamma10 <- xi10 * f1 * f0 / f

  lfdr <- gamma00 + gamma01 + gamma10

  out <- sort(lfdr, index.return = TRUE)
  fdr <- cumsum(out$x) / (1:m)
  fdr <- fdr[order(out$ix)]

  return(
    list(
      call = match.call(),
      lfdr.rep = lfdr,
      fdr.rep = fdr,
      xi00 = xi00,
      xi01 = xi01,
      xi10 = xi10,
      xi11 = xi11,
      f1 = f1,
      f2 = f2,
      loglik = loglik,
      converge = converge
    )
  )
}

