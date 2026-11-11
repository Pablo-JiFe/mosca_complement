library(effectsize)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

effector_comp <- readRDS("./results/effector_loss_dom.RDS")

# Obtain median of each perturbed source

effector_comp_summ <- 
  effector_comp %>% 
  group_by(source, condition) %>% 
  summarise(median_loss = median(absolute_loss)) %>% 
  ungroup() %>% 
  mutate(rank = rank((effector_comp %>%
                        group_by(source, condition) %>% 
                        summarise(median_loss = median(absolute_loss)) %>% 
                        ungroup())$median_loss))

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


effector_comp_summ %>% 
  mutate(median_loss = median_loss) %>% 
  ggplot(aes(x = condition, y = median_loss, fill = condition)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("#334499", "#994433"), labels = c(
    "numb" = "Selected units",
    "rand" = "Randomly selected units"
  )) +
  labs(x = "Condition", y = "Loss log10 transformed", fill = "Condition") +
  scale_x_discrete(labels = c(
    "numb" = "Selected units",
    "rand" = "Randomly selected units"
  )) + 
  theme_minimal(base_size = 25)
