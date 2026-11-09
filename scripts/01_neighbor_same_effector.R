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


set.seed(37)

equal <- list()

# Null model

for (i in 1:20000) {
  
  # Random unit
  source_sel <- metrics_active[metrics_active$source %in% sample(metrics_active$source, 1), ]
  
  type_sel_sour <- source_sel$type
  
  # Rand neighbor from diff type
  rand_neighb <- sample(metrics_active$source[!(metrics_active$type %in% type_sel_sour)], 1)
  
  # Both superclasses
  
  superclass_sel <- metrics_active$superclass[metrics_active$source %in% source_sel$source]
  superclass_rand <- metrics_active$superclass[metrics_active$source %in% rand_neighb]
  
  # Comp
  
  equal[[i]] <- superclass_sel == superclass_rand
  print(i)
}

# Chi squared

chi <- 
  data.frame(
  row.names = c("FALSE", "TRUE"),
  deg = (degeneracy_pairs_fus %>% 
    count(same_dom_eff) %>% 
    select(n) %>% 
    rename("deg" = "n")),
  null = as.numeric(table(unlist(equal)))
) %>% 
  chisq.test()

print(chi)
chi$residuals
