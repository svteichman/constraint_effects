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
fit <- readRDS("constraint_effects/results/full_model/fit1.rds")
J <- ncol(wirbel_otu_ch)
est_df <- data.frame(ind = 1:J, est = fit$B[2, ])
constraint_cats <- est_df %>% arrange(abs(est)) %>% head(25) %>% pull(ind)
constraint_set <- constraint_cats
new_B <- fit$B - 
  matrix(c(radEmu:::pseudohuber_center(fit$B[1, constraint_cats], d = 0.1),
           radEmu:::pseudohuber_center(fit$B[2, constraint_cats], d = 0.1)),
         nrow = nrow(fit$B), ncol = ncol(fit$B))

score_res <- data.frame(j = batch,
                        pval = NA,
                        time = NA)

# run score tests 
start <- proc.time()
emu_res <- fastEmuTest(constraint_cats = constraint_set, 
                       estimate_full_model = FALSE,
                       formula = ~ Group, 
                       data = wirbel_sample[ch_study_obs, ],
                       Y = wirbel_otu_ch,
                       B = new_B,
                       refit = FALSE,
                       test_kj = data.frame(k = 2, j = batch))
end <- proc.time() - start
score_res$pval <- emu_res$coef["pval"]
score_res$time <- end[3]
write_rds(score_res, paste0("constraint_effects/results/score_test_results/constraint_score_res",
                            batch, ".rds"))