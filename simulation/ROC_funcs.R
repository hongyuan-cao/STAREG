##------------------------------------------------------------------------------------------
## A revised version of the classical ROC function to make fair comparisons of the power
## We use the empirical FDR, the proportion of false discoveries in all discoveries,
## but not the false positive rate, as x-axis
##------------------------------------------------------------------------------------------
## The revised ROC function for replicability analysis based on joint FDR control procedures
#' @param states the joint states for the replicability null hypotheses across two studies.
#' @param padj the adjusted p-values for the replicability null hypotheses.
#'
#' @return a list with the following elements
#' \item{fdr}{a sequence of empirical FDR based on a range of FDR cutoffs.}
#' \item{tpr}{a sequence of power based on a range of FDR cutoffs.}
#'
revisedROC <- function(states, padj){
  fdr = 0
  tpr = 0
  thresholds = c(sort(unique(padj)),1)
  
  for (i in 1: length(thresholds)){
    thr = thresholds[i]
    fdr = c(fdr, sum(padj <= thr & !states)/max(sum(padj <= thr), 1))
    tpr = c(tpr, sum(padj <= thr & states)/sum(states))
  }

  dups <- duplicated(tpr) & duplicated(fdr)
  fdr <- fdr[!dups]
  tpr <- tpr[!dups]
  
  return(list(fdr = fdr, tpr = tpr))
}
## The revised ROC function for replicability analysis based on naive separate FDR control procedures
#' @param states the joint states for the replicability null hypotheses across two studies.
#' @param padj1 the adjusted p-values for study 1.
#' @param padj1 the adjusted p-values for study 2.
#'
#' @return a list with the following elements
#' \item{fdr}{a sequence of empirical FDR based on a range of FDR cutoffs.}
#' \item{tpr}{a sequence of power based on a range of FDR cutoffs.}
#'
revisedROC2 <- function(states, padj1, padj2){
  fdr = 0
  tpr = 0
  thresholds = c(sort(unique(c(padj1, padj2))),1)
  
  for (i in 1: length(thresholds)){
    thr = thresholds[i]
    fdr = c(fdr, sum(padj1 <= thr & padj2 <= thr & !states)/max(sum(padj1 <= thr & padj2 <= thr), 1))
    tpr = c(tpr, sum(padj1 <= thr & padj2 <= thr & states)/sum(states))
  }

  dups <- duplicated(tpr) & duplicated(fdr)
  fdr <- fdr[!dups]
  tpr <- tpr[!dups]
  
  return(list(fdr = fdr, tpr = tpr))
}
