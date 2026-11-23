suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(glmnet)
  library(igraph)
  library(purrr)
  library(readr)
  library(tibble)
})

source("scripts/_config.R")
source("scripts/_lib.R")

set.seed(RNG_SEED)

edges <- read_parquet("data/derived/type_edges.parquet")
meta <- read_parquet("data/derived/unit_meta.parquet")
metrics <- read_tsv("results/perturbome_metrics.tsv", show_col_types = FALSE)
pairs <- read_tsv("results/degeneracy_pairs.tsv", show_col_types = FALSE)
p <- read_parquet("results/perturbome_indirect.parquet")
degeneracy_null <- readr::read_tsv("results/degeneracy_null_9999.tsv")
structural_summary <- readr::read_tsv("results/degeneracy_structural_summary.tsv")
degeneracy_pairs_fus <- readr::read_tsv("results/degeneracy_pairs.tsv")




stopifnot(
  all(c("source", "effector_group", "influence") %in% names(p)),
  !anyDuplicated(meta$unit),
  setequal(unique(p$source), meta$unit),
  all(is.finite(p$influence)),
  all(p$influence >= 0)
)

n_groups <- n_distinct(p$effector_group)
stopifnot(nrow(p) == nrow(meta) * n_groups)

sig <- signature_matrix_from_long(p, source_order = meta$unit)
X <- sig$X
sources <- sig$sources
normed <- normalize_signatures(X, ZERO_TOL)

row_total <- normed$total
active <- normed$active
prob <- normed$prob

prob_active <- prob[active, , drop = FALSE]
sources_active <- sources[active]




stopifnot(
  nrow(edges) > 0L,
  !anyDuplicated(meta$unit),
  setequal(metrics$source, meta$unit),
  all(edges$pre %in% meta$unit),
  all(edges$post %in% meta$unit)
)

g <- graph_from_data_frame(
  edges |>
    transmute(from = pre, to = post, weight = count),
  directed = TRUE,
  vertices = meta |> transmute(name = unit)
)

# H3 structural divergence: functionally nearest distinct-type pairs should not
# be interpreted as degeneracy if they are simply duplicate local wiring.
out_adj <- as_adj_list(g, mode = "out")
in_adj <- as_adj_list(g, mode = "in")
vertex_names <- V(g)$name
vertex_index <- setNames(seq_along(vertex_names), vertex_names)

jaccard_for <- function(a, b, adjacency) {
  ia <- unname(vertex_index[[a]])
  ib <- unname(vertex_index[[b]])
  
  # Check for missing lookups FIRST before subsetting adjacency
  if (is.null(ia) || is.null(ib) || is.na(ia) || is.na(ib)) {
    warning(sprintf("Vertex lookup failed for: %s or %s", a, b))
    return(NA_real_)
  }
  
  na <- as.integer(adjacency[[ia]])
  nb <- as.integer(adjacency[[ib]])
  
  union_n <- length(union(na, nb))
  if (union_n == 0L) return(0)
  
  length(intersect(na, nb)) / union_n
}

edge_keys <- paste(edges$pre, edges$post, sep = "\r")
edge_key_set <- unique(edge_keys)



# H3 null: preserve the marginal distribution of every effector channel, but
# independently permute it across perturbations and renormalize rows. We then
# repeat the SAME nearest-distinct-type search, so the null includes the
# nearest-neighbour selection effect rather than comparing a maximum to random
# single pairs. A fixed random subset keeps this paper-grade check tractable.


null_n <- min(DEGENERACY_NULL_SAMPLE, nrow(prob_active))
null_ids <- sample(seq_len(nrow(prob_active)), null_n, replace = FALSE)
prob_null_base <- prob_active[null_ids, , drop = FALSE]
sources_null <- sources_active[null_ids]

observed_null_subset <- nearest_distinct_type(
  prob_null_base, sources_null, meta
)$nearest_cosine
observed_null_median <- median(observed_null_subset, na.rm = TRUE)


degeneracy_null <- map_dfr(seq_len(DEGENERACY_PERMUTATIONS), function(b) {
  permuted <- vapply(
    seq_len(ncol(prob_null_base)),
    \(j) sample(prob_null_base[, j], replace = FALSE),
    numeric(nrow(prob_null_base))
  )
  permuted <- permuted / pmax(rowSums(permuted), .Machine$double.eps)
  
  nn_null <- nearest_distinct_type(permuted, sources_null, meta)
  
  same_dom_eff <-
    tibble(
      source = nn_null$source,
      nearest_source  = nn_null$nearest_source
    ) %>%
    left_join(metrics[metrics$active, ] %>%
                select(source, dominant_effector, superclass),
              by = "source") %>%
    left_join(metrics[metrics$active, ] %>%
                select(source, dominant_effector, superclass) %>%
                rename("nearest_source" = "source",
                       "dominant_effector_near" = "dominant_effector",
                       "superclass_near" = "superclass"),
              by = "nearest_source") %>%
    mutate(same_dom_eff = dominant_effector == dominant_effector_near,
           direct_edge_a_to_b = paste(source, nearest_source, sep = "\r") %in% edge_key_set,
           direct_edge_b_to_a = paste(nearest_source, source, sep = "\r") %in% edge_key_set,
           superclass_similarity = superclass == superclass_near,
           out_neighbor_jaccard = map2_dbl(
             source, nearest_source,
             \(a, b) jaccard_for(a, b, out_adj)
           ),
           in_neighbor_jaccard = map2_dbl(
             source, nearest_source,
             \(a, b) jaccard_for(a, b, in_adj)
           )
    )
  
  same_dom_eff_quant <- sum(same_dom_eff$same_dom_eff, na.rm = TRUE)
  direct_con <- sum(ifelse(same_dom_eff$direct_edge_a_to_b | same_dom_eff$direct_edge_b_to_a, TRUE, FALSE))
  superclass_similarity <- sum(same_dom_eff$superclass_similarity, na.rm = TRUE)
  n_valid <- sum(!is.na(same_dom_eff$same_dom_eff))
  
  tibble(
    replicate = b,
    median_nearest_distinct_type_cosine = median(
      nn_null$nearest_cosine, na.rm = TRUE
    ),
    q90_nearest_distinct_type_cosine = quantile(
      nn_null$nearest_cosine, 0.90, na.rm = TRUE
    ),
    same_dom_eff_quant = sum(same_dom_eff$same_dom_eff, na.rm = TRUE),
    direct_con = sum(ifelse(same_dom_eff$direct_edge_a_to_b |same_dom_eff$direct_edge_b_to_a, TRUE, FALSE)),
    superclass_similarity = sum(same_dom_eff$superclass_similarity, na.rm = TRUE),
    n_valid = sum(!is.na(same_dom_eff$same_dom_eff)),
    prop_same_dom = same_dom_eff_quant / n_valid,
    prop_sup_dom = superclass_similarity / n_valid,
    prop_direct_con = direct_con / n_valid,
    median_out = median(same_dom_eff$out_neighbor_jaccard, na.rm = TRUE),
    median_in = median(same_dom_eff$in_neighbor_jaccard, na.rm = TRUE),
    mean_out = mean(same_dom_eff$out_neighbor_jaccard, na.rm = TRUE),
    q95_out = quantile(same_dom_eff$out_neighbor_jaccard, 0.95, na.rm = TRUE),
    mean_in = mean(same_dom_eff$in_neighbor_jaccard, na.rm = TRUE),
    q95_in = quantile(same_dom_eff$in_neighbor_jaccard, 0.95, na.rm = TRUE)
    
    
    
  )
})




chi_table_same_dom <- 
  data.frame(
    true = c(sum(degeneracy_pairs_fus$same_dom_eff), sum(degeneracy_null$same_dom_eff_quant)),
    false = c(sum(!degeneracy_pairs_fus$same_dom_eff), sum(degeneracy_null$n_valid - degeneracy_null$same_dom_eff_quant)),
    row.names = c("Observed", "Null")
  )

chi_table_same_superclass <- 
  data.frame(
    true = c(sum(degeneracy_pairs_fus$same_superclass), sum(degeneracy_null$superclass_similarity)),
    false = c(sum(!degeneracy_pairs_fus$same_superclass), sum(degeneracy_null$n_valid - degeneracy_null$superclass_similarity)),
    row.names = c("Observed", "Null")
  )

degeneracy_null_summary <- tibble(
  n_sources = null_n,
  permutations = DEGENERACY_PERMUTATIONS,
  observed_median = observed_null_median,
  observed_dir_con = 1 - structural_summary$fraction_no_direct_edge,
  observed_jaccard_out = structural_summary$median_out_neighbor_jaccard,
  observed_jaccard_in = structural_summary$median_in_neighbor_jaccard,
  observed_super_sim = sum(degeneracy_pairs_fus$same_superclass) / length(degeneracy_pairs_fus$same_superclass),
  observed_domeff_sim = sum(degeneracy_pairs_fus$same_dom_eff) / length(degeneracy_pairs_fus$same_dom_eff),
  null_median_mean = mean(
    degeneracy_null$median_nearest_distinct_type_cosine
  ),
  null_median_q95 = quantile(
    degeneracy_null$median_nearest_distinct_type_cosine, 0.95
  ),
  null_dir_con_mean = mean(
    degeneracy_null$prop_direct_con
  ),
  null_jacc_out_mean = mean(degeneracy_null$median_out),
  null_jacc_in_mean = mean(degeneracy_null$median_in),
  null_prop_super = sum(degeneracy_null$superclass_similarity) / sum(degeneracy_null$n_valid),
  null_prop_domeff = sum(degeneracy_null$same_dom_eff_quant) / sum(degeneracy_null$n_valid),
  
  permutation_p_one_sided = (
    1 + sum(
      degeneracy_null$median_nearest_distinct_type_cosine >=
        observed_null_median
    )
  ) / (1 + DEGENERACY_PERMUTATIONS),
  
  permutation_p_one_sided_con = (
    1 + sum(
      degeneracy_null$prop_direct_con >=
        1 - structural_summary$fraction_no_direct_edge
    )
  ) / (1 + DEGENERACY_PERMUTATIONS),
  
  permutation_p_one_sided_jacc_in = (
    1 + sum(
      degeneracy_null$median_in >=
        structural_summary$median_in_neighbor_jaccard
    )
  ) / (1 + DEGENERACY_PERMUTATIONS),
  
  permutation_p_one_sided_jacc_out = (
    1 + sum(
      degeneracy_null$median_in >=
        structural_summary$median_out_neighbor_jaccard
    )
  ) / (1 + DEGENERACY_PERMUTATIONS),
  
  
  p_value_chi_same_dom <- chisq.test(chi_table_same_dom)$p.value,
  
  p_value_chi_same_super <- chisq.test(chi_table_same_superclass)$p.value
  
  
  
  
)

degeneracy_null_summary$permutation_p_one_sided_jacc_in

mean(degeneracy_null$median_nearest_distinct_type_cosine)

write_tsv(degeneracy_null, "./degeneracy_null_9999.tsv")
