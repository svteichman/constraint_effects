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

# load constraint sets list 
constraint_sets <- readRDS("constraint_effects/results/constraint_sets.rds")
constraint_set <- constraint_sets[[batch]]

# get estimation results
fit <- readRDS(paste0("constraint_effects/results/full_model/fit",
                      batch, ".rds"))

J <- ncol(wirbel_otu_ch)
score_res <- data.frame(pval = rep(NA, J),
                        time = NA)
# run score tests 
if (batch == 1) {
  for (j in 1:J) {
    start <- proc.time()
    emu_res <- emuFit(formula = ~ Group, 
                      data = wirbel_sample[ch_study_obs, ],
                      Y = wirbel_otu_ch,
                      fitted_model = fit,
                      refit = FALSE,
                      test_kj = data.frame(k = 2, j = j))
    end <- proc.time() - start
    score_res$pval[j] <- emu_res$coef[j, "pval"]
    score_res$time[j] <- end[3]
  }
} else {
  for (j in 1:J) {
    start <- proc.time()
    emu_res <- fastEmuTest(constraint_cats = constraint_set, 
                           estimate_full_model = FALSE,
                           formula = ~ Group, 
                           data = wirbel_sample[ch_study_obs, ],
                           Y = wirbel_otu_ch,
                           fitted_model = fit,
                           refit = FALSE,
                           test_kj = data.frame(k = 2, j = j))
    end <- proc.time() - start
    score_res$pval[j] <- emu_res$coef[1, "pval"]
    score_res$time[j] <- end[3]
  }
}
write_rds(score_res, paste0("constraint_effects/results/score_test_results/score_res",
                            batch, ".rds"))