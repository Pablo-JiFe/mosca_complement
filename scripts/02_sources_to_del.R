library(arrow)
library(dplyr)
library(FNN)
library(irlba)
library(purrr)
library(readr)
library(tibble)
library(tidyr)
library(uwot)

source("scripts/_config.R")
source("scripts/_lib.R")

set.seed(RNG_SEED)

p <- read_parquet("results/perturbome_indirect.parquet")
meta <- read_parquet("data/derived/unit_meta.parquet")
metrics <- read_tsv("results/perturbome_metrics.tsv")

# Singature 

sig <- signature_matrix_from_long(p, source_order = meta$unit)
X <- sig$X

# Transpose because then normalization will be applied to effectors 

X_t <- t(X)

a <- normalize_signatures(X_t)

# object with effector, its source that has hihghest influence over the effector, the scaled value and that sourece units dominannt effector

eff_metr <- tibble(
  effector = rownames(a$prob),
  dominant_source = colnames(a$prob)[max.col(a$prob)],
  value_dom_source = a$prob[cbind(seq_len(nrow(a$prob)), max.col(a$prob))]
) %>%
  left_join(
    metrics %>%
      select(source, dominant_effector) %>%
      rename("dominant_source" = source),
    by  = "dominant_source"
  ) 

#group by source to obtain over how many was it the highest evffector

eff_stat <- eff_metr %>%
  group_by(dominant_source) %>%
  summarize(n = n(), median = median(value_dom_source)) %>%
  left_join(
    metrics %>%
      select(source, superclass, dominant_effector) %>%
      rename("dominant_source" = source),
    by  = "dominant_source"
  ) %>%
  arrange(desc(n)) 


# names of top 5 highest effectros for perturbation 

saveRDS(eff_stat$dominant_source, "results/source_to_del.RDS")
