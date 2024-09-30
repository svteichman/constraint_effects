library(tidyverse)
library(qvalue)
library(ggpubr)
library(radEmu)
library(latex2exp)

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
wirbel_sample_ch <- wirbel_sample[ch_study_obs, ]

# load in estimation results 
estimates <- matrix(nrow = 46, ncol = 758)
for (i in 1:45) {
  emuEst <- readRDS(paste0("results/estimation_results/coef", i, ".rds"))
  estimates[i, ] <- emuEst$estimate
}
emuEst <- readRDS("results/estimation_results/coef_mean.rds")
estimates[46, ] <- emuEst$estimate
dd_estimates <- matrix(nrow = 84, ncol = 758)
colnames(dd_estimates) <- colnames(wirbel_otu_ch)
for (i in 1:84) {
  emuEst <- readRDS(paste0("results/estimation_results/coef_dd", i, ".rds"))
  dd_estimates[i, colnames(dd_estimates) %in% emuEst$category] <- emuEst$estimate
}

# calculate differences 
sq_diffs <- sapply(2:46, function(x) {(estimates[x, 1] - estimates[1, 1])^2})
sq_dd_diffs <- c(as.vector(sapply(1:4, function(x) {(dd_estimates[x, 1] - estimates[1, 1])^2})), 
              as.vector(sapply(5:84, function(x) {mean((dd_estimates[x, ] - estimates[1, ])^2, na.rm = TRUE)})))

# compare estimates to radEmu median, using constraint of size 50 
plot_df <- data.frame(med = estimates[1, ], prev = estimates[4, ], 
                      dd = dd_estimates[3, ], ss = dd_estimates[27, ],
                      thin = dd_estimates[47, ])
p1 <- ggplot(plot_df, aes(x = med, y = dd)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = TeX(r"(Estimates ($S_{all}$))"), 
       y = TeX(r"(Estimates ($S_{dd}$))")) + 
  theme_bw(base_size = 14)
p2 <- ggplot(plot_df, aes(x = med, y = ss)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = TeX(r"(Estimates ($S_{all}$))"), 
       y = TeX(r"(Estimates ($S_{ss}$))")) + 
  theme_bw(base_size = 14)
p3 <- ggplot(plot_df, aes(x = med, y = thin)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = TeX(r"(Estimates ($S_{all}$))"), 
       y = TeX(r"(Estimates ($S_{th}$))")) + 
  theme_bw(base_size = 14)
p_est <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
annotate_figure(p_est, top = text_grob("Comparing estimates across reference sets of size 50", size = 16))
p_est
ggsave("figures/compare_ests.pdf", height = 6, width = 6)

# see how many NA's for each run 
sapply(1:84, function(x) {sum(is.na(dd_estimates[x, ]))})[seq(1, 84, 4)]

# load in inference results 
pvals <- matrix(nrow = 46, ncol = 758)
times <- matrix(nrow = 46, ncol = 758)
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
# radEmu with mean constraint
for (j in 1:758) {
  emuScore <- readRDS(paste0("results/score_test_results/mean_score_res", j, ".rds"))
  pvals[46, j] <- emuScore$pval
  times[46, j] <- emuScore$time
}
# data driven constraint results 
pvals_dd <- matrix(nrow = 84, ncol = 758)
times_dd <- matrix(nrow = 84, ncol = 758)
for (i in 1:84) {
  emuScore <- readRDS(paste0("results/score_test_results/score_res_dd", i, ".rds"))
  pvals_dd[i, colnames(wirbel_otu_ch) %in% emuScore$cat] <- emuScore$pval
  times_dd[i, colnames(wirbel_otu_ch) %in% emuScore$cat] <- emuScore$time
}

# analyze times
rowMeans(times)
apply(times, 1, max)/60
rowSums(times)/3600

rowMeans(times_dd, na.rm = T)
apply(times_dd, 1, max, na.rm = T)/60
rowSums(times_dd, na.rm = T)/3600

# analyze significant results 
qvals <- apply(pvals[1:46, ], 1, function(x) {qvalue(x)$qvalues})
sig_cats <- lapply(1:46, function(x) {which(qvals[, x] <= 0.1)})

qvals_dd <- apply(pvals_dd[1:84, ], 1, function(x) {qvalue(x)$qvalues})
sig_cats_dd <- lapply(1:84, function(x) {which(qvals_dd[, x] <= 0.1)})

# compare subset of p-values and q-values 
plot_df <- data.frame(med = pvals[1, ], prev = pvals[4, ], 
                      dd = pvals_dd[3, ], ss = pvals_dd[27, ],
                      thin = pvals_dd[47, ])
p2 <- ggplot(plot_df, aes(x = med, y = dd)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = TeX(r"(P-values ($S_{all}$))"), 
       y = TeX(r"(P-values ($S_{dd}$))")) +
  theme_bw(base_size = 14)
p3 <- ggplot(plot_df, aes(x = med, y = ss)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = TeX(r"(P-values ($S_{all}$))"), 
       y = TeX(r"(P-values ($S_{ss}$))")) +
  theme_bw(base_size = 14)
p4 <- ggplot(plot_df, aes(x = med, y = thin)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = TeX(r"(P-values ($S_{all}$))"), 
       y = TeX(r"(P-values ($S_{th}$))")) +
  theme_bw(base_size = 14)
p_inf <- ggarrange(p2, p3, p4, ncol = 1, nrow = 3)
annotate_figure(p_inf, top = text_grob("Comparing p-values across reference sets of size 50", size = 16))
ggsave("figures/compare_pvals.pdf", height = 6, width = 6)
p_all <- ggarrange(p_est, p_inf, nrow = 1)
annotate_figure(p_all, top = text_grob("Comparing estimates and p-values across reference sets of size 50", size = 16))
ggsave("figures/compare_ests_and_pvals.pdf", height = 8, width = 8)


plot_df <- data.frame(med = qvals[, 1], prev = qvals[, 4], 
                      dd = qvals_dd[, 3], ss = qvals_dd[, 27],
                      thin = qvals_dd[, 47])
p1 <- ggplot(plot_df, aes(x = med, y = prev)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = "Pseudo median", y = "Prevalence subset") +
  xlim(c(0, 0.55)) + ylim(c(0, 0.55))
p2 <- ggplot(plot_df, aes(x = med, y = dd)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = "Pseudo median", y = "Data-driven subset") +
  xlim(c(0, 0.55)) + ylim(c(0, 0.55))
p3 <- ggplot(plot_df, aes(x = med, y = ss)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = "Pseudo median", y = "Data-driven subset (sample splitting)") +
  xlim(c(0, 0.55)) + ylim(c(0, 0.55))
p4 <- ggplot(plot_df, aes(x = med, y = thin)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
  labs(x = "Pseudo median", y = "Data-driven subset (thinning)") +
  xlim(c(0, 0.55)) + ylim(c(0, 0.55))
p <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
annotate_figure(p, top = "Comparing q-values across constraints over subsets of size 50")
ggsave("figures/compare_qvals.png", height = 6, width = 6)

# overall comparison plots
pval_mse <- c(apply(pvals[2:45, ], 1, function(x) {sqrt(mean((x - pvals[1, ])^2))}), 
              apply(pvals_dd, 1, function(x) {sqrt(mean((x - pvals[1, ])^2, na.rm = T))}))
#qval_mse <- c(apply(qvals[, 2:45], 2, function(x) {mean((x - qvals[, 1])^2)}), 
#              apply(qvals_dd, 2, function(x) {mean((x - qvals[, 1])^2, na.rm = T)}))
est_mse <- c(apply(estimates[2:45, ], 1, function(x) {sqrt(mean((x - estimates[1, ])^2))}), 
             apply(dd_estimates, 1, function(x) {sqrt(mean((x - estimates[1, ])^2, na.rm = T))}))
sizes <- c(10, 30, 50, 100)
plot_df <- data.frame(est_mse = est_mse,
                      p_mse = pval_mse, 
                      #q_mse = qval_mse,
                      type = c(paste0("presence ", sizes),
                               rep(paste0("random ", sizes), each = 10),
                               paste0("data-driven ", sizes),
                               rep(paste0("data-driven (sample split) ", sizes), 10),
                               rep(paste0("data-driven (thin) ", sizes), 10)),
                      agg_type = c(rep("presence", 4), rep("random", 40),
                                   rep("data-driven", 4), rep("data-driven (sample split)", 40),
                                   rep("data-driven (thin)", 40)),
                      full_time = c(apply(times[2:45, ], 1, sum),
                                    apply(times_dd, 1, sum, na.rm = T)))
#plot_df$num_disc_agree <- c(sapply(2:45, function(x) {sum(sig_cats[[x]] %in% sig_cats[[1]])}),
#                            sapply(1:84, function(x) {sum(sig_cats_dd[[x]] %in% sig_cats[[1]])}))
#plot_df$num_disc_disagree <- c(sapply(2:45, function(x) {sum(!(sig_cats[[x]] %in% sig_cats[[1]]))}),
#                               sapply(1:84, function(x) {sum(!(sig_cats_dd[[x]] %in% sig_cats[[1]]))}))
plot_df$type <- factor(plot_df$type, levels = c(paste0("data-driven ", sizes),
                                                paste0("data-driven (sample split) ", sizes),
                                                paste0("data-driven (thin) ", sizes),
                                                paste0("presence ", sizes),
                                                paste0("random ", sizes)))
plot_df <- plot_df %>%
  filter(!(type %in% c("presence 10", "presence 30", "presence 50", "presence 100")))
#ggplot(plot_df, aes(x = est_mse, y = p_mse, color = type, shape = agg_type)) + 
p1 <- ggplot(plot_df, aes(x = sqrt(est_mse), y = sqrt(p_mse), color = agg_type)) + 
  geom_point() + 
  #scale_shape_manual(values = c(17, 13, 12, 19, 8)) + 
  labs(x = "Root mean squared difference\nbetween estimates",
       y = "Root mean squared difference between p-values",
       color = "Reference set") + 
  theme_bw(base_size = 14) + 
  ylim(c(0,0.6)) + 
  guides(color = guide_legend(position = "bottom", nrow = 1)) + 
  scale_color_hue(labels = c("data-driven" = TeX(r"($S_{dd}$)"),
                             "data-driven (sample split)" = TeX(r"($S_{ss}$)"),
                             "data-driven (thin)" = TeX(r"($S_{th}$)"),
                             "random" = "random"))
p1
ggsave("figures/p_mse_by_est_mse.png", height = 6, width = 6)
# p2 <- ggplot(plot_df, aes(x = sqrt(est_mse), y = sqrt(q_mse), color = agg_type)) + 
#   geom_point() + 
#   #scale_shape_manual(values = c(8, 19)) + 
#   labs(x = "Square root of MSE of estimates",
#        y = "Square root of MSE across q-values",
#        color = "Constraint",
#        shape = "Constraint Type") + 
#   theme_bw()
# p2
# ggsave("figures/q_mse_by_est_mse.png", height = 6, width = 6)
p3 <- ggplot(plot_df, aes(x = sqrt(est_mse), y = full_time/60, color = agg_type)) + 
  geom_point() + 
  #scale_shape_manual(values = c(8, 19)) + 
  labs(x = "Root mean squared difference\nbetween estimates",
       y = "Time to run all tests serially (minutes)",
       color = "Constraint",
       shape = "Constraint Type") + 
  theme_bw(base_size = 14) + 
  guides(color = guide_legend(position = "bottom", nrow = 2)) + 
  scale_color_hue(labels = c("data-driven" = TeX(r"($S_{dd}$)"),
                             "data-driven (sample split)" = TeX(r"($S_{ss}$)"),
                             "data-driven (thin)" = TeX(r"($S_{th}$)"),
                             "random" = "random"))
p3
ggsave("figures/full_time_by_est_mse.png", height = 6, width = 6)
# p4 <- ggplot(plot_df, aes(x = num_disc_agree, y = num_disc_disagree, color = agg_type)) + 
#   geom_jitter() + 
#   #scale_shape_manual(values = c(8, 19)) + 
#   labs(x = "Number discoveries also in pseudo-median discovery set",
#        y = "Number discoveries not in pseudo-median discovery set",
#        color = "Constraint",
#        shape = "Constraint Type") + 
#   theme_bw()
# p4
# ggsave("figures/discoveries.png", height = 6, width = 6)

p <- ggarrange(p1, p3, nrow = 1, common.legend = T, legend = "bottom")
annotate_figure(p, 
                top = text_grob("Comparing results using reference sets to full taxon set",
                                size = 16))
ggsave("figures/est_mse_p_mse_time.pdf", height = 8, width = 8)


