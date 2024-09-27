library(tidyverse)
library(radEmu)

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
wirbel_sample_ch <- wirbel_sample[ch_study_obs, ]

# load constraint sets list 
constraint_sets <- readRDS("constraint_effects/results/dd_constraint_sets.rds")
constraint_set <- constraint_sets[[batch]]

# get test set
if (batch %in% 1:4) {
  test_otu <- wirbel_otu_ch
  test_samp <- wirbel_sample_ch
} else if (batch %in% 5:44) {
  trial <- ceiling((batch - 4) / 4)
  test_inds <- readRDS("constraint_effects/results/ss_test_sets.rds")
  test_otu <- wirbel_otu_ch[test_inds[[trial]], ]
  test_samp <- wirbel_sample_ch[test_inds[[trial]], ]
} else {
  trial <- ceiling((batch - 44) / 4)
  test_sets <- readRDS("constraint_effects/results/thin_test_sets.rds")
  test_otu <- test_sets[[trial]]
  test_samp <- wirbel_sample_ch
}
# remove columns or rows of Y with only 0 observations in test set 
zero_cols <- which(colSums(test_otu) == 0)
if (length(zero_cols) > 0) {
  test_otu <- test_otu[, -zero_cols]
}
zero_rows <- which(rowSums(test_otu) == 0)
if (length(zero_rows) > 0) {
  test_samp <- test_samp[-zero_rows, ]
  test_otu <- test_otu[-zero_rows, ]
}

J <- ncol(wirbel_otu_ch)
if (length(zero_cols) > 0) {
  ind_df <- data.frame(original = 1:J,
                       new = NA)
  ind_df$new[-zero_cols] <- 1:ncol(test_otu)
  upd_constraint_set <- ind_df$new[ind_df$original %in% constraint_set]
  upd_constraint_set <- upd_constraint_set[!is.na(upd_constraint_set)]
} else {
  upd_constraint_set <- constraint_set
}

# constraint functions
constraint_fn <- (function(x) {
  radEmu:::pseudohuber_center(x[upd_constraint_set], d = .1)
})
constraint_grad_fn <- (function(x) {
  grad <- rep(0, length(x))
  grad[upd_constraint_set] <-
    radEmu:::dpseudohuber_center_dx(x[upd_constraint_set], d = .1)
  return(grad)
})

# fit model
ch_fit <- emuFit(formula = ~ Group, 
                 data = test_samp,
                 Y = test_otu,
                 run_score_tests = FALSE,
                 constraint_fn = constraint_fn, 
                 constraint_grad_fn = constraint_grad_fn,
                 constraint_param = NA) 
write_rds(ch_fit$coef, paste0("constraint_effects/results/estimation_results/coef_dd",
                              batch, ".rds"))
write_rds(ch_fit, paste0("constraint_effects/results/full_model/fit_dd",
                              batch, ".rds"))

