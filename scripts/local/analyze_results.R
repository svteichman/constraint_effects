library(tidyverse)
library(qvalue)

# load in estimation results 
estimates <- matrix(nrow = 46, ncol = 758)
for (i in 1:45) {
  emuEst <- readRDS(paste0("results/estimation_results/coef", i, ".rds"))
  estimates[i, ] <- emuEst$estimate
}
emuEst <- readRDS("results/estimation_results/coef_mean.rds")
estimates[46, ] <- emuEst$estimate

# calculate differences 
diffs <- sapply(2:46, function(x) {estimates[x, 1] - estimates[1, 1]})
summary(estimates[1, ])
mean(abs(estimates[1, ] < 0.15))
ggplot(data.frame(x = estimates[1,], y = estimates[4,]), aes(x = x, y = y)) + 
  geom_point() + geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  labs(x = "Estimates with full constraint set",
       y = "Estimates with subset of 50 taxon constraint set") + 
  theme_bw()
ggsave("figures/estimate_med_p50.png", height = 6, width = 6)

constraint_sets <- readRDS("results/constraint_sets.rds")
summary(estimates[1, ])
summary(estimates[1, constraint_sets[[2]]])
summary(estimates[1, constraint_sets[[3]]])
summary(estimates[1, constraint_sets[[4]]])
summary(estimates[1, constraint_sets[[5]]])

# load in inference results 
pvals <- matrix(nrow = 47, ncol = 758)
times <- matrix(nrow = 47, ncol = 758)
# full radEmu score tests
for (j in 1:758) {
  emuScore <- readRDS(paste0("results/score_test_results/full_score_res", j, ".rds"))
  pvals[1, j] <- emuScore$pval
  times[1, j] <- emuScore$time
}
# fastEmu score tests, various constraints
for (i in 2:45) {
  emuScore <- readRDS(paste0("results/score_test_results/score_res", i, ".rds"))
  pvals[i, ] <- emuScore$pval
  times[i, ] <- emuScore$time
}
# radEmu with constraint 4
for (j in 1:758) {
  emuScore <- readRDS(paste0("results/score_test_results/s4_score_res", j, ".rds"))
  pvals[46, j] <- emuScore$pval
  times[46, j] <- emuScore$time
}
# radEmu with mean constraint
for (j in 1:758) {
  emuScore <- readRDS(paste0("results/score_test_results/mean_score_res", j, ".rds"))
  pvals[47, j] <- emuScore$pval
  times[47, j] <- emuScore$time
}
# analyze times
rowMeans(times)
apply(times, 1, max)/60
rowSums(times)/3600
# analyze significant results 
qvals <- apply(pvals[1:47, ], 1, function(x) {qvalue(x)$qvalues})
sig_cats <- lapply(1:47, function(x) {which(qvals[, x] <= 0.1)})
plot_df <- data.frame(p_med = pvals[1, ], p10 = pvals[2, ], p30 = pvals[3, ], 
                      p50 = pvals[4, ], p100 = pvals[5, ])
plot(plot_df)
cor(plot_df)
plot_df <- plot_df %>%
  mutate(rank_med = rank(p_med),
         rank10 = rank(p10),
         rank30 = rank(p30),
         rank50 = rank(p50),
         rank100 = rank(p100))
cor(plot_df %>% dplyr::select(contains("rank")))
plot(plot_df %>% dplyr::select(contains("rank")))
plot_df %>% dplyr::select(contains("rank")) %>%
  arrange(rank_med) %>% head(20)
ggplot(plot_df, aes(x = p_med, y = p50)) + 
  geom_point() + geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  labs(x = "P-values with full constraint set",
       y = "P-values with subset of 50 taxon constraint set") + 
  theme_bw()
ggsave("figures/pval_med_p50.png", height = 6, width = 6)
plot_df$q_med <- qvals[, 1]
plot_df$q50 <- qvals[, 4]
ggplot(plot_df, aes(x = q_med, y = q50)) + 
  geom_point() + geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  labs(x = "Q-values with full constraint set",
       y = "Q-values with subset of 50 taxon constraint set") + 
  theme_bw() + 
  xlim(c(0, 0.5)) + 
  ylim(c(0, 0.5))
ggsave("figures/qval_med_p50.png", height = 6, width = 6)

pval_mse <- apply(pvals, 1, function(x) {sum((x - pvals[1, ])^2)})
qval_mse <- apply(qvals, 1, function(x) {sum((x - qvals[1, ])^2)})
sizes <- c(10, 30, 50, 100)
plot_df <- data.frame(p_mse = pval_mse[2:45], diff = abs(diffs)[2:45], 
                      q_mse = qval_mse[2:45],
                      type = c(paste0("presence ", sizes),
                               rep(paste0("random ", sizes), each = 10)),
                      agg_type = c(rep("presence", 4), rep("random", 40)),
                      full_time = apply(times[2:45, ], 1, sum))
plot_df$num_disc_agree <- sapply(2:45, function(x) {sum(sig_cats[[x]] %in% sig_cats[[1]])})
plot_df$num_disc_disagree <- sapply(2:45, function(x) {sum(!(sig_cats[[x]] %in% sig_cats[[1]]))})
plot_df$type <- factor(plot_df$type, levels = c(paste0("presence ", sizes),
                                                paste0("random ", sizes)))
ggplot(plot_df, aes(x = diff, y = p_mse, color = type, shape = agg_type)) + 
  geom_point() + 
  scale_shape_manual(values = c(8, 19)) + 
  labs(x = "Absolute difference in estimates",
       y = "MSE across p-values",
       color = "Constraint",
       shape = "Constraint Type") + 
  theme_bw()
ggsave("figures/p_mse_by_est_diff.png", height = 6, width = 6)
ggplot(plot_df, aes(x = diff, y = q_mse, color = type, shape = agg_type)) + 
  geom_point() + 
  scale_shape_manual(values = c(8, 19)) + 
  labs(x = "Absolute difference in estimates",
       y = "MSE across q-values",
       color = "Constraint",
       shape = "Constraint Type") + 
  theme_bw()
ggsave("figures/q_mse_by_est_diff.png", height = 6, width = 6)
ggplot(plot_df, aes(x = diff, y = full_time, color = type, shape = agg_type)) + 
  geom_point() + 
  scale_shape_manual(values = c(8, 19)) + 
  labs(x = "Absolute difference in estimates",
       y = "Time to run all tests (seconds)",
       color = "Constraint",
       shape = "Constraint Type") + 
  theme_bw()
ggsave("figures/full_time_by_est_diff.png", height = 6, width = 6)
ggplot(plot_df, aes(x = num_disc_agree, y = num_disc_disagree, color = type, shape = agg_type)) + 
  geom_point() + 
  scale_shape_manual(values = c(8, 19)) + 
  labs(x = "Number discoveries also in pseudo-median discovery set",
       y = "Number discoveries not in pseudo-median discovery set",
       color = "Constraint",
       shape = "Constraint Type") + 
  theme_bw()
ggsave("figures/discoveries.png", height = 6, width = 6)

# radEmu with constraint 4
median(times[46, ])
median(times[4, ])
median(times[46, ]) / median(times[4, ])
max(times[46, ])
max(times[4, ])
max(times[46, ]) / max(times[4, ])
sum(times[46, ])
sum(times[4, ])
sum(times[46, ]) / sum(times[4, ])
cor(pvals[4, ], pvals[46, ])

# data
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

