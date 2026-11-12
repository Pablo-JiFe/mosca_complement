library(ggplot2)
library(dplyr)
library(patchwork)
library(ggnewscale)
library(readr)
library(ggsci)
library(cowplot)
library(grid)
library(forcats) 


pca_var <- read_tsv("results/perturbome_pca_variance.tsv", show_col_types = FALSE)
parallel <- read_tsv("results/perturbome_parallel_analysis.tsv", show_col_types = FALSE)
metrics <- read_tsv("results/perturbome_metrics.tsv")
umap <- read_tsv("results/perturbome_umap.tsv", show_col_types = FALSE) |>
  left_join(metrics, by = "source")
deg <- read_tsv("results/degeneracy_structural.tsv", show_col_types = FALSE) |>
  mutate(neighbor_jaccard = pmax(out_neighbor_jaccard, in_neighbor_jaccard))
effector_comp_summ <- readRDS("results/effector_comp_summ.RDS")
eff_metr <- readRDS("results/dom_source.RDS")
sources_to_del_data <- readRDS("results/source_to_del_data.RDS")


fill_manual_pal <-   scale_fill_manual(
  values = c(
    "#4477AA", "#EE6677", "#228833", "#CCBB44",
    "#66CCEE", "lightpink", "#BBBBBB", "#332288",
    "#44AA99", "slateblue1", "coral1", "#DDCC77",
    "black", "#88CCEE", "#661100", "#6699CC",
    "#AA4466", "#004488", "#DDAA33"
  ),
  name = "Superclass"
) 

pca_var %>% 
  ggplot(aes(x = pc, y = cumulative)) +
  geom_point() +
  geom_vline(xintercept = 53, colour = "red") +
  theme_minimal() +
  labs(x = "PC", y = "Cumulative variance explained") +
  theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.5)) + 
  scale_x_continuous(breaks = c(1, seq(10, 90, by = 10), 99)) 


parallel %>% 
  ggplot(aes(x = pc)) + 
  geom_line(
    aes(y = observed_eigenvalue, colour = "Observed", linetype = "Observed"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = null_mean, colour = "Null mean", linetype = "Null mean"),
    linewidth = 1
  ) +
  labs(
    x = "PC",
    y = "Eigenvalue",
    colour = "Distribution",
    linetype = "Distribution"
  ) +
  scale_colour_manual(
    values = c("Observed" = "red", "Null mean" = "blue")
  ) +
  scale_linetype_manual(
    values = c("Observed" = "dashed", "Null mean" = "solid")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = c(1, seq(10, 90, by = 10), 99)) +
  geom_vline(xintercept = 53, colour = "black")

library(scales)

umap1 <- 
  umap %>%
  ggplot(aes(
    x = UMAP1,
    y = UMAP2,
    colour = breadth
  )) +
  paletteer::scale_colour_paletteer_c("grDevices::Blue-Red 3") + 
  geom_point(
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_minimal(base_size = 20) +
  coord_fixed(ratio = 1)


umap2 <- 
  umap %>%
  ggplot(aes(
    x = UMAP1,
    y = UMAP2,
    fill = superclass,
  )) +
  geom_point(
    shape = 21,
    size = 5,
    color = "black",
    stroke = 0.4,
    alpha = 0.8
  ) +
  fill_manual_pal +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 9, face = "bold"),
    legend.spacing.y = unit(0.05, "cm")
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE,
      override.aes = list(
        size = 10,
        alpha = 1,
        color = "black",
        stroke = 0.5
      )
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(
        alpha = 1,
        color = "black"
      )
    )
  )


umap_no_legend <- umap2 + theme(legend.position = "none")


umap_legend <- get_legend(umap2)



plot1_3 <- umap %>% 
  mutate(superclass = fct_reorder(superclass, breadth, .fun = mean, .desc = TRUE)) %>%
  ggplot(aes(y = breadth, x = superclass, fill = superclass)) + 
  geom_boxplot(size = 0.5) + 
  theme_minimal(base_size = 25) +
  theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.5)) + 
  fill_manual_pal +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) + 
  labs(
    x = "Breadth",
    y = "Superclass",
    fill = "Superclass"
  ) 

tag_theme <- theme(
  plot.tag = element_text(size = 16, face = "bold"),
  plot.tag.position = c(0, 1)
)

umap1_tag <- umap1 +
  labs(tag = "A") +
  tag_theme

umap2_tag <- umap_no_legend +
  labs(tag = "B") +
  tag_theme

plot1_3_tag <- plot1_3 +
  labs(tag = "C") +
  tag_theme


left_plots <- ((umap1_tag + umap2_tag) / plot1_3_tag)


final_layout <- (left_plots | wrap_elements(umap_legend)) +
  plot_layout(widths = c(4, 1))


final_layout

neighbor_plot_1 <- deg |>
  ggplot() +
  aes(y = neighbor_jaccard, x = cosine_similarity, colour = factor(ifelse(direct_edge_a_to_b | direct_edge_b_to_a, 1, 0))) +
  geom_point(size = 1, alpha = 0.3) +
  scale_color_manual(values = c("#53a4fc", "#333333"), labels = c("No direct connection", "Direct connection")) + 
  labs(
    y = "Max(in/out neighbor Jaccard)",
    x = "Perturbation cosine similarity",
    colour = "Direct connection"
  ) +
  theme_minimal(base_size = 25)



neighbor_plot_2 <- deg %>% 
  ggplot(aes(x = cosine_similarity)) + 
  geom_histogram(binwidth = 0.005, fill = "#53a4fc", colour = "white") +
  geom_vline(xintercept = median(deg$cosine_similarity), linetype = "dashed", colour = "#982222", linewidth = 1.5) + 
  annotate(
    geom = "text",
    x = median(deg$cosine_similarity),              
    y = 3000,              
    label = "Median pair cosine = 0.994",
    hjust = 1.05
  ) + 
  geom_vline(xintercept = 0.8354, linetype = "dashed", colour = "#333333", linewidth = 1.5) + 
  annotate(
    geom = "text",
    x = 0.8354,              
    y = 3000,              
    label = "Median pair cosine null = 0.8354",
    hjust = 1.05
  ) + 
  theme_minimal(base_size = 20) + 
  labs(x = "Cosine similarity between nearest neighbor pairs from different neuronal types", y = "Number of unique pairs")



neighbor_plot_2 / neighbor_plot_1 + 
  patchwork::plot_annotation(tag_levels = "A")



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

