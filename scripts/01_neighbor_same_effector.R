library(readr)
library(dplyr)

metrics <- read_tsv( "results/perturbome_metrics.tsv")
degeneracy_pairs <- read_tsv( "results/degeneracy_pairs.tsv")

# Keep active

metrics_active <- metrics[metrics$active,]

# Join t o observe which dominant effector belongs to each source

degeneracy_pairs_fus <- 
  degeneracy_pairs %>% 
  left_join(
    metrics_active %>% 
      select(source, dominant_effector) %>% 
      dplyr::rename("source_a" = "source",
                    "dominant_effector_a" = "dominant_effector"),
    by = "source_a"
  ) %>% 
  left_join(
    metrics_active %>% 
      select(source, dominant_effector) %>% 
      dplyr::rename("source_b" = "source",
                    "dominant_effector_b" = "dominant_effector"),
    by = "source_b"
  ) %>% 
  mutate(same_dom_eff = ifelse(
    dominant_effector_a == dominant_effector_b,
    TRUE,
    FALSE
  ))




write_tsv(degeneracy_pairs_fus, "results/degeneracy_pairs.tsv")
