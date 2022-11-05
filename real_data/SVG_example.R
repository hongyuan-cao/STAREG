## Replicability analysis of SVG detection: mouse olfactory bulb
# Load the count data for mouse olfactory bulb in Replicate 1 and Replicate 8 
counts1 <- read.table("./real_data/Rep1_MOB_count_matrix-1.tsv", check.names = F)
counts8 <- read.table("./real_data/Rep8_MOB_count_matrix-1.tsv", check.names = F)
counts1 <- t(counts1)
counts8 <- t(counts8)

# Analyze the two SRT datasets with SPARK (referring to the SPARK pacakge)
library(SPARK)
# Replicate 1
location1 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 2)))
rownames(location1) <- colnames(counts1)
spark1 <- CreateSPARKObject(counts = counts1, location = location1, 
                            percentage = 0.1, min_total_counts = 10)
spark1@lib_size <- apply(spark1@counts, 2, sum)
spark1 <- spark.vc(spark1, covariates = NULL, lib_size = spark1@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark1 <- spark.test(spark1, check_positive = T, verbose = T)

# Replicate 8
location8 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts8), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts8), split = "x"), "[", 2)))
rownames(location8) <- colnames(counts8)
spark8 <- CreateSPARKObject(counts = counts8, location = location8, 
                            percentage = 0.1, min_total_counts = 10)
spark8@lib_size <- apply(spark8@counts, 2, sum)
spark8 <- spark.vc(spark8, covariates = NULL, lib_size = spark8@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark8 <- spark.test(spark8, check_positive = T, verbose = T)

# Extract the well-calibrated p-values from the analysis results of the two studies and match them by gene.
p1 <- spark1@res_mtest$combined_pvalue
names(p1) <- rownames(spark1@counts)
p2 <- spark8@res_mtest$combined_pvalue
names(p2) <- rownames(spark8@counts)
overlap <- intersect(names(p1),names(p2))
pvals1 = p1[overlap]
pvals2 = p2[overlap]

# Replicability analysis of the two studies using RepEM
library(STAREG)
alpha <- 0.05
rep.obj <- STAREG(pvals1, pvals2)
# Output the replicable SVGs across the two SRT studies
rep.svgs <- overlap[which(rep.obj$fdr.rep <= alpha)]  