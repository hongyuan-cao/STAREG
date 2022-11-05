library(STAREG)
library(ggplot2)
source('./simulation/ROC_funcs.R')
source('./simulation/data_generation.R')

## Data generation
m = 10000
xi00 = 0.9
xi01 = 0.025
xi10 = 0.025
xi11 = 0.05
mu1 = 2
mu2 = 2.5
sigma1 = 1
sigma2 = 1

data.obj <- data_generation(m, xi00, xi01, xi10, xi11, mu1, mu2, sigma1, sigma2)
p1 = data.obj$pvals1
p2 = data.obj$pvals2
states1 = data.obj$states1
states2 = data.obj$states2

## Replicability analysis
# BH
padj1.bh <- p.adjust(p1, method = 'BH')
padj2.bh <- p.adjust(p2, method = 'BH')

# MaxP
maxp <- apply(cbind(p1, p2), 1, max)
padj.maxp <- p.adjust(maxp, method = "BH")

# STAREG
res.rep <- STAREG(p1, p2, est.pi0 = TRUE)
padj.rep <- res.rep$fdr.rep
f1.rep = res.rep$f1
f2.rep = res.rep$f2
xi00.hat = res.rep$xi00
xi01.hat = res.rep$xi01
xi10.hat = res.rep$xi10
xi11.hat = res.rep$xi11

################################################################################
## TPR~FDR plot
################################################################################

rroc.bh <- revisedROC2(states1*states2, padj1.bh, padj2.bh)
rroc.maxp <- revisedROC(states1*states2, padj.maxp)
rroc.rep <- revisedROC(states1*states2, padj.rep)

rroc.bh <- as.data.frame(rroc.bh)
rroc.maxp <- as.data.frame(rroc.maxp)
rroc.rep <- as.data.frame(rroc.rep)

rroc.bh$method <- "BH"
rroc.maxp$method <- "MaxP"
rroc.rep$method <- "STAREG"

res <- rbind(rroc.maxp, rroc.bh, rroc.rep)

ggplot(res, aes(fdr,tpr, group = method, color = method, shape = method)) + 
  # geom_point() + 
  scale_colour_manual(name="",  
                      values = c("MaxP"="#6496D2", "BH"="#F4B183", "STAREG"="#8FBC8F")) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) + 
  geom_line(size = 1.5) + xlab("FDR") + ylab("Power") +
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1, colour = "#44546A", size = 1, linetype = 2) +
  theme(legend.position = c(0.02,0.98),
        legend.justification = c(0.02,0.98),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.6, linetype = "solid"))




