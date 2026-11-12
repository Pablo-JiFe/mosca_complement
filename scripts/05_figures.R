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


umap1

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


