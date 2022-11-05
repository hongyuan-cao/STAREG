## numeric data generated from normal distribution and z-test
#' @param m the number of hypotheses in each study
#' @param xi00 the prior probability of state (0, 0)
#' @param xi01 the prior probability of state (0, 1)
#' @param xi10 the prior probability of state (1, 0)
#' @param xi11 the prior probability of state (1, 1)
#' @param mu1 the mean value for the alternative distribution in study 1
#' @param mu2 the mean value for the alternative distribution in study 2
#' @param sigma1 the standard deviation for the alternative distribution in study 1
#' @param sigma2 the standard deviation for the alternative distribution in study 2
#'
#' @return a list with the following elements
#' {\item}{pvals1}{the simulated p-values obtained from study 1}
#' {\item}{pvals2}{the simulated p-values obtained from study 2}
#' {\item}{states1}{the simulated hidden states for study 1}
#' {\item}{states2}{the simulated hidden states for study 2}
#'
data_generation <- function(m = 10000, xi00 = 0.9, xi01 = 0.025, xi10 = 0.025, xi11 = 0.05,
                            mu1 = 2, mu2 = 2, sigma1 = 1, sigma2 = 1)
{
  h = sample(0:3, m, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
  states1 = rep(0, m)
  states1[c(which(h==2), which(h==3))] = 1
  states2 = rep(0, m)
  states2[c(which(h==1), which(h==3))] = 1

  stat1 = rnorm(m, states1*mu1, sigma1)
  stat2 = rnorm(m, states2*mu2, sigma2)

  p1 = 1 - pnorm(stat1, mean = 0, sd = sigma1)
  p2 = 1 - pnorm(stat2, mean = 0, sd = sigma2)

  return(list(pvals1 = p1, pvals2 = p2, states1 = states1, states2 = states2))
}
