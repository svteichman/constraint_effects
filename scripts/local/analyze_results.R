library(tidyverse)
library(qvalue)

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
ggplot(data.frame(x = estimates[1,], y = estimates[3,]), aes(x = x, y = y)) + 
  geom_point() + geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  labs(x = "Estimates with full constraint set",
       y = "Estimates with subset of 30 taxon constraint set")

# load in inference results 
pvals <- matrix(nrow = 45, ncol = 758)
times <- matrix(nrow = 45, ncol = 758)
for (i in 2:45) {
  emuScore <- readRDS(paste0("results/score_test_results/score_res", i, ".rds"))
  pvals[i, ] <- emuScore$pval
  times[i, ] <- emuScore$time
}
# analyze times
rowMeans(times)
apply(times, 1, max)/60
rowSums(times)/3600
# analyze significant results 
qvals <- apply(pvals[2:45, ], 1, function(x) {qvalue(x)$qvalues})
sig_cats <- lapply(1:44, function(x) {which(qvals[, x] <= 0.1)})
sig_cats_comb <- unlist(sig_cats)
sig_cats_unique <- unique(sig_cats_comb)
plot_df <- data.frame(p10 = pvals[2, ], p30 = pvals[3, ], p50 = pvals[4, ],
                      p100 = pvals[5, ])
plot(plot_df)
plot_df <- plot_df %>%
  mutate(rank10 = rank(p10),
         rank30 = rank(p30),
         rank50 = rank(p50),
         rank100 = rank(p100))
cor(plot_df %>% dplyr::select(contains("rank")))
plot(plot_df %>% dplyr::select(contains("rank")))
plot_df %>% dplyr::select(contains("rank")) %>%
  arrange(rank100) %>% head(20)

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

