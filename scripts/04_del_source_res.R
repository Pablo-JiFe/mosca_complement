library(effectsize)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)

effector_comp <- readRDS("./results/effector_loss_dom_40.RDS")
sources_to_del_data <- readRDS("results/source_to_del_data.RDS")
metrics <- read_tsv("results/perturbome_metrics.tsv")
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

rank_biserial(x = breadth_median$median_loss ~ factor(breadth_median$over_med))





#Boxplot

first_bxpl <- effector_comp_summ %>% 
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
  theme_minimal(base_size = 25) +
  coord_fixed(ratio = 0.15)

# Breadth plot


breadth_plt <- effector_comp_summ %>% 
  ggplot(aes(x = breadth, y = median_loss, colour = condition)) + 
  geom_point(position = "jitter") +
  theme_minimal(base_size = 20) + 
  labs(x = "Breadth", y = "Median loss", colour = "Condition") + 
  scale_color_manual(values = c("#334499", "#994433"), labels = c(
    "numb" = "Selected units",
    "rand" = "Randomly selected units"
  )) +
  geom_text_repel(
    data = subset(effector_comp_summ, rank <= 18 & condition == "numb"),
    aes(label = source, x = breadth, y = median_loss),
    size = 3,
    vjust = -2,
    segment.color = "gray",     
    segment.size = 1.1,           
    segment.alpha = 0.7,          
    min.segment.length = 0        
  ) 

second_bxplt <- sources_to_del_data %>% 
  ggplot(aes(x = factor(over_elb), y = median_loss, fill = factor(over_elb))) + 
  geom_boxplot() +
  scale_x_discrete(labels = c(
    "0" = "Under elbow threshold",
    "1" = "Over elbow threshold"
  )) + 
  scale_y_log10() +
  theme_minimal(base_size = 20) + 
  labs(x = "Sources with breadth over the median", y = "Median loss", colour = "Sources with breadth over the median", fill = "Over or under elbow threshold") + 
  scale_fill_manual(values = c("#334499", "#994433"), labels = c(
    "0" = "Under elbow threshold",
    "1" = "Over elbow threshold"
  )) +
  coord_fixed(ratio = 0.15)


fill_manual_pal <-   scale_fill_manual(
  values = c(
    "#4477AA", "#EE6677", "#228833", "#888888",
    "#66CCEE", "lightpink", "#324537", "#332288",
    "#44AA99", "slateblue1", "coral1", "#DDCC77",
    "#DDAA33", "#88CCEE", "#661100", "#343434",
    "#AA4466", "#004488"
  ))

alluvial <- 
  eff_metr %>% 
  mutate(over_elb = factor(ifelse(dominant_source %in% sources_to_del_data$dominant_source[sources_to_del_data$over_elb == 1], 1, 0))) %>% 
  filter(over_elb == 1) %>% 
  select(effector, dominant_source, value_dom_source, over_elb) %>% 
  ggplot(aes(axis1 = effector, axis2 = dominant_source)) + 
  geom_alluvium(aes(fill = dominant_source), width = 1/6) +
  geom_stratum(width = 1/6, color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
  scale_x_discrete(limits = c("Effector group", "Dominant source"), expand = c(.05, .05)) +
  fill_manual_pal + 
  theme_minimal(base_size = 15) +
    theme(legend.position = "none")

  
(((first_bxpl | second_bxplt | patchwork::plot_spacer()) + patchwork::plot_layout(widths = c(1,1,0.1))) / ((breadth_plt | alluvial) + patchwork::plot_layout(widths = c(3,4)))) + 
  patchwork::plot_layout(heights = c(2,3)) +
    patchwork::plot_annotation(tag_levels = "A")
  
