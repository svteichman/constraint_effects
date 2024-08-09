library(tidyverse)
library(radEmu)

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
J <- ncol(wirbel_otu_ch)

# run estimation 
constraint_fn <- (function(x) {
  mean(x)
})
constraint_grad_fn <- (function(x) {
  grad <- rep(1/J, length(x))
  return(grad)
})
ch_fit <- emuFit(formula = ~ Group, 
                 data = wirbel_sample[ch_study_obs, ],
                 Y = wirbel_otu_ch,
                 run_score_tests = FALSE,
                 constraint_fn = constraint_fn, 
                 constraint_grad_fn = constraint_grad_fn,
                 constraint_param = NA) 
write_rds(ch_fit$coef, "constraint_effects/results/estimation_results/coef_mean.rds")
write_rds(ch_fit, "constraint_effects/results/full_model/fit_mean.rds")

