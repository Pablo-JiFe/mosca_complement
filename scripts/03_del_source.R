library(arrow)
library(dplyr)
library(Matrix)
library(readr)
library(tibble)
library(tidyr)
library(foreach)
library(doParallel)



raw_long <- read_parquet("./results/perturbome_raw.parquet")
meta <- read_parquet("./data/derived/unit_meta.parquet")
edges <- read_parquet("./data/derived/type_edges.parquet")
effectors <- read_parquet("./data/derived/effector_membership.parquet")
metrics <- read_tsv("./results/perturbome_metrics.tsv")
degeneracy_pairs <- read_tsv("./results/degeneracy_pairs.tsv")
source_to_del <- readRDS("./results/source_to_del.rds")

total_runs <- 140

source("_lib.R")

# Same initial weight mat

wm <- build_weight_matrix(edges = edges,
                          units = meta$unit,
                          rho = 0.9)

W <- wm$W

# And effector mat

Bm <- build_effector_matrix(effectors = effectors,
                            units = meta$unit,
                            group_col = "effector_group_fine")

B <- Bm$B

# Raw to wide format

raw_wide <- raw_long %>%
  pivot_wider(values_from = influence, names_from = source) %>%
  column_to_rownames("effector_group") %>%
  as.matrix()

# Control to see if they where formoded correctly 

x <- (Matrix::Diagonal(dim(W)[1]) - t(W)) %*% t(raw_wide)
y <- (B %>% as.matrix())

all(colnames(x) == rownames(y))
all(colnames(y) == rownames(x))



set.seed(37)

# Obtian random units not belonging tho the soruces to delete

random_sources <- sample(
  metrics$source[
    metrics$active &
      !(metrics$source %in% source_to_del)
  ],
  size = total_runs - length(source_to_del),
  replace = FALSE
)

# Add all to same obj

all_sources <- c(source_to_del, random_sources)

# Id of type of analysis

type_of_source <- c(rep("numb", length(source_to_del)), rep("rand", total_runs - length(source_to_del)))

cores = 14
cl <- makeCluster(cores) 
clusterEvalQ(cl, {
  library(Matrix)
  library(tibble)
})
registerDoParallel(cl)

effector_loss_par_res <- 
  foreach(e=1:total_runs, .packages = c("Matrix", "tibble")) %dopar% {
    
    a <- all_sources[e]
    
    I <- Matrix::Diagonal(nrow(W))
    
    X_original <- raw_wide
    
    
    W_lesion <- W
    W_lesion[a, ] <- 0
    W_lesion[, a] <- 0
    
    B_lesion <- B
    B_lesion[, a] <- 0
    
    X_lesion <- t(solve(I - t(W_lesion), t(B_lesion)))
    
    effector_loss_par <- tibble(
      lesion_id = e,
      source = a,
      condition = type_of_source[e],
      effector_group = rownames(X_original),
      original = rowSums(X_original),
      post_lesion = rowSums(X_lesion),
      absolute_loss = rowSums(X_original) - rowSums(X_lesion)
    )
    
    print(paste(e))
    
    effector_loss_par
    
  }
stopCluster(cl)


effector_tibb <- effector_loss_par_res[[1]] %>%
  mutate(name = names(effector_loss_par_res)[1])

for (i in 2:total_runs) {
  effector_tibb <-
    rbind(effector_tibb, effector_loss_par_res[[i]] %>%
            mutate(name = names(effector_loss_par_res)[i]))
}


effector_comp <- effector_tibb %>% 
  left_join(
    metrics %>%
      select(source, type, dominant_effector, superclass),
    by = "source", multiple = "all"
  ) 

saveRDS(effector_comp, "effector_loss_dom.RDS")
