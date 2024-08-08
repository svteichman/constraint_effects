library(tidyverse)
library(radEmu)
library(fastEmu)

# get batch as command line argument 
args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  batch <- 1
} else {
  arg <- args[length(args)]
  batch <- abs(readr::parse_number(arg))
}

# load data 
data("wirbel_sample")
data("wirbel_otu")

# prep data
wirbel_sample$Group <- factor(wirbel_sample$Group, levels = c("CTR","CRC"))
ch_study_obs <- which(wirbel_sample$Country %in% c("CHI"))

wirbel_otu_ch <- wirbel_otu[ch_study_obs, ]
sum(rowSums(wirbel_otu_ch) == 0) # no samples have a count sum of 0 
sum(colSums(wirbel_otu_ch) == 0) # 87 categories have a count sum of 0

categories_to_rm <- which(colSums(wirbel_otu_ch) == 0)
wirbel_otu_ch <- wirbel_otu_ch[, -categories_to_rm]
sum(colSums(wirbel_otu_ch) == 0)

# get estimation results
fit <- readRDS("constraint_effects/results/full_model/fit4.rds")

# load constraint sets list 
constraint_sets <- readRDS("constraint_effects/results/constraint_sets.rds")
constraint_set <- constraint_sets[[4]]
constraint_fn <- (function(x) {
  radEmu:::pseudohuber_center(x[constraint_set], d = .1)
})
constraint_grad_fn <- (function(x) {
  grad <- rep(0, length(x))
  grad[constraint_set] <-
    radEmu:::dpseudohuber_center_dx(x[constraint_set], d = .1)
  return(grad)
})

score_res <- data.frame(j = batch,
                        pval = NA,
                        time = NA)

# run score tests 
start <- proc.time()
emu_res <- emuFit(formula = ~ Group, 
                  data = wirbel_sample[ch_study_obs, ],
                  Y = wirbel_otu_ch,
                  fitted_model = fit,
                  refit = FALSE,
                  test_kj = data.frame(k = 2, j = batch),
                  constraint_fn = constraint_fn, 
                  constraint_grad_fn = constraint_grad_fn,
                  constraint_param = NA)
end <- proc.time() - start
score_res$pval <- emu_res$coef[batch, "pval"]
score_res$time <- end[3]
write_rds(score_res, paste0("constraint_effects/results/score_test_results/s4_score_res",
                            batch, ".rds"))