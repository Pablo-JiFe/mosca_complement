library(ggplot2)
library(dplyr)
library(patchwork)
library(ggnewscale)
library(readr)
library(ggsci)
library(cowplot)
library(gridExtra)
library(grid)
library(ggalluvial)
library(forcats) 
library(stringr)


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
sensitivity <- readr::read_tsv("results/sensitivity/model_sensitivity.tsv")



fill_manual_pal <-   scale_fill_manual(
  labels = superclass_names,
  values = c(
    "#4477AA", "#EE6677", "#228833", "#CCBB44",
    "#66CCEE", "lightpink", "#BBBBBB", "#332288",
    "#44AA99", "slateblue1", "coral1", "#DDCC77",
    "black", "#88CCEE", "#661100", "#6699CC",
    "#AA4466", "#004488", "#DDAA33"
  ),
  name = "Superclass"
) 

fill_manual_supp <-   scale_fill_manual(
  values = c(
    "lightpink", "#EE6677", 
    "#66CCEE", "#004488","#332288", 
    "slateblue1", "#228833", "#44AA99","#53d377", "#DDCC77",
    "#DDAA33", "#88CCEE", "#661100", "#343434",
    "#AA4466"
  ))



fill_manual_alluvial <-   scale_fill_manual(
  values = c(
    "#4477AA", "#EE6677", "#228833", "#661100",
    "#66CCEE", "lightpink", "#324537", "#332288",
    "#44AA99", "slateblue1", "coral1", "#DDCC77",
    "#DDAA33", "#88CCEE", "#888888", "#343434",
    "#AA4466", "#004488"
  ))

superclass_names <- 
  (
    metrics %>% 
  count(superclass) %>% 
  mutate(names = gsub("(.*)_(.*)", "\\1 \\2", str_to_title(superclass)))
  )$names

names(superclass_names) <-   (
  metrics %>% 
    count(superclass) %>% 
    mutate(names = gsub("(.*)_(.*)", "\\1 \\2", str_to_title(superclass)))
)$superclass

library(scales)

umap1 <- 
  umap %>%
  ggplot(aes(
    x = UMAP1,
    y = UMAP2,
    fill = breadth
  )) +
  paletteer::scale_fill_paletteer_c("ggthemes::Purple") + 
  geom_point(shape = 21,
             size = 5,
             color = "black",
             stroke = 0.4,
             alpha = 0.8) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    fill = "Breadth"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 9, face = "bold"),
    legend.spacing.y = unit(0.05, "cm")
  )


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
  theme_classic(base_size = 30) + 
  guides(color = guide_legend(override.aes = list(size = 8, alpha = 1)))



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
  theme_classic(base_size = 28) +
  labs(x = "Cosine similarity between nearest neighbor pairs from different neuronal types", y = "Number of unique pairs")



neighbor_plot_2 / neighbor_plot_1 + 
  patchwork::plot_annotation(tag_levels = "A")



#Boxplot

first_bxpl <- effector_comp_summ %>% 
  ggplot(aes(x = condition, y = median_loss, fill = condition)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("#334499", "#994433"), labels = c(
    "numb" = "High impact",
    "rand" = "Randomly High impact"
  )) +
  labs(x = "Condition", y = "Median loss (log10)", fill = "Condition") +
  scale_x_discrete(labels = c(
    "numb" = "High impact",
    "rand" = "Randomly High impact"
  )) + 
  theme_minimal(base_size = 20) +
  theme(legend.position = "none") 

# Breadth plot



table_subset <- effector_comp_summ %>%
  filter(rank <= 18 & condition == "numb") %>%
  select(source, breadth, median_loss) %>% 
  mutate(breadth = round(breadth, 3),
         median_loss = round(median_loss, 3)) %>% 
  rename(
    "Source" = "source",
    "Breadth" = "breadth",
    "Median loss (log10)" = "median_loss"
  )



my_table <- tableGrob(
  table_subset, 
  rows = NULL, 
  theme = ttheme_minimal(base_size = 15, base_colour = "black", base_family = "",
                 parse = FALSE, padding = unit(c(1, 1), "mm"))
  
)



breadth_plt <- effector_comp_summ %>% 
  ggplot(aes(x = breadth, y = median_loss, colour = condition)) + 
  geom_point(position = "jitter") +
  theme_minimal(base_size = 20) + 
  labs(x = "Breadth", y = "Median loss", colour = "Condition") + 
  scale_color_manual(values = c("#334499", "#994433"), labels = c(
    "numb" = "High impact",
    "rand" = "Randomly High impact"
  )) +
  ggrepel::geom_text_repel(
    data = subset(effector_comp_summ, rank <= 18 & condition == "numb"),
    aes(label = source, x = breadth, y = median_loss),
    size = 3,
    vjust = -4.2,
    hjust = 0.5,
    segment.color = "gray",     
    segment.size = 1.1,           
    segment.alpha = 0.7,          
    min.segment.length = 0        
  ) +
  annotation_custom(
    my_table, 
    xmin = 0.02,   
    xmax = 0.35,   
    ymin = 0.1,   
    ymax = 0.40    
  ) 


second_bxplt <- 
  sources_to_del_data %>% 
  mutate(over_elb = relevel(factor(over_elb), ref = 1)) %>% 
  ggplot(aes(x = fct_reorder(over_elb, as.numeric(over_elb), .desc = TRUE), y = median_loss, fill = over_elb)) + 
  geom_boxplot() +
  scale_x_discrete(labels = c(
    
    "1" = "Over elbow threshold",
    "0" = "Under elbow threshold"
  )) + 
  scale_y_log10() +
  theme_minimal(base_size = 20) + 
  labs(x = "Sources divided by maximmum distance", y = "Median loss (log10)", fill = "Elbow threshold") + 
  scale_fill_manual(values = c("#004488", "slateblue1"), labels = c(
    
    "1" = "Over elbow threshold",
    "0" = "Under elbow threshold"
  )) 

alluvial <- eff_metr %>% 
  mutate(over_elb = factor(ifelse(dominant_source %in% sources_to_del_data$dominant_source[sources_to_del_data$over_elb == 1], 1, 0))) %>% 
  filter(over_elb == 1, dominant_source %in% (sources_to_del_data %>% 
                                                arrange(rank))$dominant_source[1:5]) %>% 
  mutate(effector = gsub("(.*)_(.*)", "\\1 \\2", str_to_title(effector))) %>% 
  select(effector, dominant_source, value_dom_source, over_elb) %>% 
  ggplot(aes(axis1 = effector, axis2 = dominant_source)) + 
  geom_alluvium(aes(fill = dominant_source), width = 0.7) +
  geom_stratum(width = 0.4, color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Effector group", "Dominant source") ,
    expand = expansion(mult = c(0.2, 0.2))
    ) + 
  fill_manual_alluvial + 
  theme_minimal(base_size = 1) +
  theme(legend.position = "none")
  



(((first_bxpl | breadth_plt) + patchwork::plot_layout(widths = c(2,4))) / 
    ((second_bxplt | alluvial) + patchwork::plot_layout(widths = c(3,4)))) + 
  patchwork::plot_annotation(tag_levels = "A")

  pca_var_p <- pca_var %>% 
  ggplot(aes(x = pc, y = cumulative)) +
  geom_point() +
  geom_vline(xintercept = 53, colour = "red") +
  theme_minimal(base_size = 25) +
  labs(x = "PC", y = "Cumulative variance explained") +
  theme(panel.grid.major = element_line(color = "grey80", linewidth = 0.5)) + 
  scale_x_continuous(breaks = c(1, seq(10, 90, by = 10), 99)) + 
  coord_fixed(50)


null_pca <- 
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
  theme_minimal(base_size = 25) +
  theme(
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = c(1, seq(10, 90, by = 10), 99)) +
  geom_vline(xintercept = 53, colour = "black") + 
  coord_fixed(15000)



(pca_var_p | null_pca) + 
  patchwork::plot_annotation(tag_levels = "A") 

sensitivity_pl <- sensitivity %>%  
  mutate(
    effector_grouping = gsub("(.*)_(.*)_(.*)", "\\1 \\2 \\3", effector_grouping),
    scenario = factor(scenario, levels = c("primary", "rho", "synapse_threshold", "effector_taxonomy")),
    
    "type" = case_when(
      scenario == "primary" ~ "primary",
      scenario == "rho" ~ paste0("Rho: ", rho),
      scenario == "synapse_threshold" ~ as.character(min_synapses),
      scenario == "effector_taxonomy" ~ paste0("E", substr(start = 2, stop = stringr::str_length(effector_grouping), x = effector_grouping))
    ),
    
    type = factor(
      type,
      levels = c(
        "primary",
        paste0("Rho: ", sort(unique(rho[scenario == "rho"]))),
        as.character(sort(unique(min_synapses[scenario == "synapse_threshold"]))),
        unique(type[scenario == "effector_taxonomy"])
      )
    )
  ) %>%  
  tidyr::pivot_longer(
    cols = c(median_breadth, median_nearest_distinct_type_cosine, pc50),
    names_to = "test",
    values_to = "values"
  ) %>%  
  ggplot(aes(x = type, y = values, colour = scenario)) +  
  geom_point() +
  geom_line(aes(x = type, y = values, group = scenario)) + 
  facet_grid(
    test ~ scenario, 
    scales = "free", 
    space = "free_x",
    labeller = labeller(
      scenario = c(
        "primary" = "Primary",
        "rho" = "Rho Sensitivity",
        "synapse_threshold" = "Synapse Threshold",
        "effector_taxonomy" = "Effector Grouping"
      ),
      test = c(
        "median_breadth" = "Median Breadth",
        "median_nearest_distinct_type_cosine" = "Nearest Distinct Cosine",
        "pc50" = "PC50"
      )
    )
  ) +  
  scale_colour_manual(values = c("#df1a1a", "#0077af", "#5a32cc", "#1a02cc")) + 
  theme_minimal(base_size = 25) +  
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(x = "Test", y = "Values") +
  theme(legend.position = "none")


sensitivity_pl
