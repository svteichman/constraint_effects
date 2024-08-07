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

# check on sparsity
mean(wirbel_otu_ch == 0)
# 75% zeros 

## make list of settings 

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

