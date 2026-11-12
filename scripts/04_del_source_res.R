library(effectsize)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)

effector_comp <- readRDS("./results/effector_loss_dom_40.RDS")
sources_to_del_data <- readRDS("results/source_to_del_data.RDS")
metrics <- readr::read_tsv("results/perturbome_metrics.tsv")
eff_metr <- readRDS("results/dom_source.RDS")

# Obtain median of each perturbed source

effector_comp_summ <- 
  effector_comp %>% 
  group_by(source, condition) %>% 
  summarise(median_loss = median(absolute_loss)) %>% 
  ungroup() %>% 
  mutate(rank = rank(-(effector_comp %>%
                        group_by(source, condition) %>% 
                        summarise(median_loss = median(absolute_loss)) %>% 
                        ungroup())$median_loss)) %>%
  left_join(
    metrics %>%
      select(source, breadth), 
    by  = "source"
  )

# Sep by condition


cond_numb_summ <- effector_comp_summ %>%
  filter(condition == "numb") 

cond_rand_summ <- effector_comp_summ %>%
  filter(condition == "rand") 


wil <- wilcox.test(cond_numb_summ$median_loss, cond_rand_summ$median_loss)

rank_biserial(cond_numb_summ$median_loss, cond_rand_summ$median_loss)

# Accumulated median

effector_comp_summ_sel <- effector_comp_summ %>% 
  filter(condition =="numb") %>%
  arrange((rank)) %>% 
  mutate(accumulated_median = cumsum(median_loss)) 


find_elbow <- function(x, y) {
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  p1 <- c(x_norm[1], y_norm[1])
  p2 <- c(x_norm[length(x_norm)], y_norm[length(y_norm)])
  
  line_vec <- p2 - p1
  line_len <- sqrt(sum(line_vec^2))
  
  distances <- sapply(seq_along(x_norm), function(i) {
    pt <- c(x_norm[i], y_norm[i])
    vec_from_p1 <- pt - p1
    abs(line_vec[1] * vec_from_p1[2] - line_vec[2] * vec_from_p1[1]) / line_len
  })

  dist_ratio <- distances / max(distances)
  best_idx <- which.max(distances)
  
  return(list(
    optimal = x[best_idx],
    distances = distances,
    distance_ratios = dist_ratio
  ))
}


res <- find_elbow(x = effector_comp_summ_sel$rank, y = effector_comp_summ_sel$accumulated_median)



effector_comp_summ_sel %>% 
  ggplot(aes(x = rank, y  = accumulated_median)) + 
  geom_line() + 
  geom_text_repel(
    data = subset(effector_comp_summ_sel, rank <= res$optimal ),
    aes(label = source, x = rank, y = accumulated_median),
    size = 3,
    vjust = -2,
    hjust = 1.5,
    segment.color = "gray50",     
    segment.size = 1.2,           
    segment.alpha = 0.7,          
    min.segment.length = 0        
  ) + 
  geom_point(aes(x = rank, y  = accumulated_median)) + 
  annotate("segment", x = min(effector_comp_summ_sel$rank), y = min(effector_comp_summ_sel$accumulated_median), xend = max(effector_comp_summ_sel$rank), yend = max(effector_comp_summ_sel$accumulated_median), color = "red", linewidth = 1) + 
  annotate("point", x = res$optimal, y = effector_comp_summ_sel$accumulated_median[effector_comp_summ_sel$rank == res$optimal], size = 4, colour = "red") +
  scale_x_continuous(breaks = seq(0, 140, by = 5)) + 
  theme_minimal(base_size = 25) + 
  labs( x = "Rank", y = "Accumulated median loss")




sources_to_del_data <- sources_to_del_data %>% 
  mutate(rank_imp  = rank(-n),
         rank_eff = rank(-median)) %>% 
  left_join(effector_comp_summ %>% 
              select(-breadth) %>% 
              rename("dominant_source" = "source"),
            by = "dominant_source") %>% 
  mutate(over_elb = ifelse(rank <= res$optimal, 1, 0)) 

elb_abs <- wilcox.test(formula = sources_to_del_data$median_loss ~ factor(sources_to_del_data$over_elb))

rank_biserial(x = sources_to_del_data$median_loss ~ factor(sources_to_del_data$over_elb))


saveRDS(effector_comp_summ, "results/effector_comp_summ.RDS")
saveRDS(sources_to_del_data, "results/source_to_del_data.RDS")
