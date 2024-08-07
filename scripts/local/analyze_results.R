library(tidyverse)

# load in estimation results 
estimates <- matrix(nrow = 45, ncol = 758)
for (i in 1:45) {
  emuEst <- readRDS(paste0("results/estimation_results/coef", i, ".rds"))
  estimates[i, ] <- emuEst$estimate
}

# calculate differences 
diffs <- sapply(2:45, function(x) {estimates[x, 1] - estimates[1, 1]})
summary(estimates[1, ])
mean(abs(estimates[1, ] < 0.15))

# load in inference results 
pvals <- matrix(nrow = 45, ncol = 758)
for (i in 1:45) {
  emuScore <- readRDS(paste0("results/score_test_results/score_res", i, ".rds"))
  pvals[i, ] <- emuEst$pval
}
