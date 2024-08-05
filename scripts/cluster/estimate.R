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

# load constraint sets list 
constraint_sets <- readRDS("constraint_effects/results/constraint_sets.rds")
constraint_set <- constraint_sets[[batch]]

# run estimation 
ch_fit <- emuFit(formula = ~ Group, 
                 data = wirbel_sample[ch_study_obs, ],
                 Y = wirbel_otu_ch,
                 run_score_tests = FALSE) 
write_rds(ch_fit$coef, paste0("constraint_effects/results/estimation_results/coef",
                              batch, ".rds"))
write_rds(ch_fit, paste0("constraint_effects/results/full_model/fit",
                              batch, ".rds"))

