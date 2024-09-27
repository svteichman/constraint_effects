library(tidyverse)
library(radEmu)
library(datathin)

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

# check on sparsity
mean(wirbel_otu_ch == 0)
# 75% zeros 

## make list of initial settings 

constraint_sets <- vector(mode = "list", length = 45)
# entire set of otus
J <- ncol(wirbel_otu_ch)
constraint_sets[[1]] <- 1:J
# subsets of otus that are most present 
pres_vec <- colSums(wirbel_otu_ch > 0)
pres_order <- data.frame(index = 1:J,
                         present_num = pres_vec) %>%
  arrange(desc(present_num))
ggplot(pres_order, aes(x = present_num)) + 
  geom_histogram(bins = 15) + 
  labs(x = "Presence in samples", y = "Count")
sizes <- c(10, 30, 50, 100)
for (s in 1:length(sizes)) {
  size <- sizes[s]
  constraint_sets[[s + 1]] <- pres_order$index[1:size]
}
# random subsets of otus
ind <- 6
for (size in sizes) {
  for (seed in 1:10) {
    set.seed(seed) 
    constraint_sets[[ind]] <- sample(1:J, size)
    ind <- ind + 1
  }
}

## make list of data-driven settings 

dd_constraint_sets <- vector(mode = "list", length = 84)
ss_test_sets <- vector(mode = "list", length = 10)
thin_test_sets <- vector(mode = "list", length = 10)

# fully data snooping 
fit <- readRDS("results/full_model/fit.rds")
J <- ncol(wirbel_otu_ch)
n <- nrow(wirbel_otu_ch)
n_samp <- floor(0.25 * n)
est_df <- data.frame(ind = 1:J, est = fit$B[2, ])
dd_constraint_sets[[1]] <- est_df %>% arrange(abs(est)) %>% head(10) %>% pull(ind)
dd_constraint_sets[[2]] <- est_df %>% arrange(abs(est)) %>% head(30) %>% pull(ind)
dd_constraint_sets[[3]] <- est_df %>% arrange(abs(est)) %>% head(50) %>% pull(ind)
dd_constraint_sets[[4]] <- est_df %>% arrange(abs(est)) %>% head(100) %>% pull(ind)

# sample splitting 
samp_split <- function(trial) {
  print(trial)
  set.seed(trial) 
  
  # make training and test sets
  train_ind <- sample(1:n, n_samp)
  test_ind <- (1:n)[-train_ind]
  
  # run estimation on training set 
  train_samp <- wirbel_sample_ch[train_ind, ]
  train_otu <- wirbel_otu_ch[train_ind, ]
  zero_cols <- which(colSums(train_otu) == 0)
  train_otu <- train_otu[, -zero_cols]
  zero_rows <- which(rowSums(train_otu) == 0)
  if (length(zero_rows) > 0) {
    train_samp <- train_samp[-zero_rows, ]
    train_otu <- train_otu[-zero_rows, ]
    train_ind <- train_ind[-zero_rows]
  }
  est <- emuFit(formula = ~ Group, 
                data = train_samp,
                Y = train_otu,
                run_score_tests = FALSE,
                verbose = TRUE) 
  est_df <- data.frame(ind = (1:J)[-zero_cols], est = est$coef$estimate)
  trial_inds <- (trial) * 4 + 1:4
  constraints <- vector(mode = "list", length = 4)
  constraints[[1]] <- est_df %>% arrange(abs(est)) %>% head(10) %>% pull(ind)
  constraints[[2]] <- est_df %>% arrange(abs(est)) %>% head(30) %>% pull(ind)
  constraints[[3]] <- est_df %>% arrange(abs(est)) %>% head(50) %>% pull(ind)
  constraints[[4]] <- est_df %>% arrange(abs(est)) %>% head(100) %>% pull(ind)
  
  return(list(constraint_list = constraints, test_ind = test_ind))
}
for (t in 1:10) {
  res <- samp_split(t) 
  dd_constraint_sets[t * 4 + 1:4] <- res$constraint_list
  ss_test_sets[[t]] <- res$test_ind
}

# thinning 
thin <- function(trial) {
  print(trial)
  set.seed(trial) 
  
  # make training and test sets
  otu_thin <- datathin(wirbel_otu_ch, family = "poisson", epsilon = c(0.25, 0.75))
  train_otu <- otu_thin[, , 1]
  test_otu <- otu_thin[, , 2]
  
  # run estimation on training set 
  zero_cols <- which(colSums(train_otu) == 0)
  train_otu <- train_otu[, -zero_cols]
  zero_rows <- which(rowSums(train_otu) == 0)
  if (length(zero_rows) > 0) {
    train_samp <- train_samp[-zero_rows, ]
    train_otu <- train_otu[-zero_rows, ]
    train_ind <- train_ind[-zero_rows]
  }
  est <- emuFit(formula = ~ Group, 
                data = wirbel_sample_ch,
                Y = train_otu,
                run_score_tests = FALSE,
                verbose = TRUE) 
  est_df <- data.frame(ind = (1:J)[-zero_cols], est = est$coef$estimate)
  constraints <- vector(mode = "list", length = 4)
  constraints[[1]] <- est_df %>% arrange(abs(est)) %>% head(10) %>% pull(ind)
  constraints[[2]] <- est_df %>% arrange(abs(est)) %>% head(30) %>% pull(ind)
  constraints[[3]] <- est_df %>% arrange(abs(est)) %>% head(50) %>% pull(ind)
  constraints[[4]] <- est_df %>% arrange(abs(est)) %>% head(100) %>% pull(ind)
  
  return(list(constraint_list = constraints, test = test_otu))
}
for (t in 2:10) {
  res <- thin(t) 
  dd_constraint_sets[40 + t * 4 + 1:4] <- res$constraint_list
  thin_test_sets[[t]] <- res$test
}
