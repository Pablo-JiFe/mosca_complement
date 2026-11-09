library(effectsize)
library(dplyr)
library(tibble)
library(tidyr)

effector_comp <- readRDS("./results/effector_loss_dom.RDS")

# Obtain median of each perturbed source

effector_comp_summ <- 
  effector_comp %>% 
  group_by(source, condition) %>% 
  summarise(median_loss = median(absolute_loss)) %>% 
  ungroup()

# Sep by condition


cond_numb_summ <- effector_comp %>%
  filter(condition == "numb") %>% 
  group_by(source, condition) %>% 
  summarise(median_loss = median(absolute_loss)) %>% 
  ungroup()

cond_rand_summ <- effector_comp %>%
  filter(condition == "rand") %>% 
  group_by(source, condition) %>% 
  summarise(median_loss = median(absolute_loss))  %>% 
  ungroup()



wil <- wilcox.test(cond_numb_summ$median_loss, cond_rand_summ$median_loss)

rank_biserial(cond_numb_summ$median_loss, cond_rand_summ$median_loss)



effector_comp %>% 
  ggplot(aes(x = condition, y = absolute_loss)) + 
  geom_boxplot() + 
  scale_y_log10()
