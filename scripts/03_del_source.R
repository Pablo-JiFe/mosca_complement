library(arrow)
library(dplyr)
library(Matrix)
library(readr)
library(tibble)
library(tidyr)
library(foreach)
library(doParallel)
library(parallel)



raw_long <- read_parquet("./results/perturbome_raw.parquet")
meta <- read_parquet("./data/derived/unit_meta.parquet")
edges <- read_parquet("./data/derived/type_edges.parquet")
effectors <- read_parquet("./data/derived/effector_membership.parquet")
metrics <- read_tsv("./results/perturbome_metrics.tsv")
degeneracy_pairs <- read_tsv("./results/degeneracy_pairs.tsv")
source_to_del <- readRDS("./results/source_to_del.rds")


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
  size = 100 - length(source_to_del),
  replace = FALSE
)

# Add all to same obj

all_sources <- c(source_to_del, random_sources)

# Id of type of analysis, first n correspond to selected and next untul 100 random

type_of_source <- c(rep("numb", length(source_to_del)), rep("rand", 100 - length(source_to_del)))

#cores

cores = detectCores(logical = TRUE) - 1
cl <- makeCluster(cores) 
clusterEvalQ(cl, {
  library(Matrix)
  library(tibble)
})
registerDoParallel(cl)

# parallel loop

effector_loss_par_res <- 
  foreach(e=1:100, .packages = c("Matrix", "tibble")) %dopar% {
    
    #source at loop e
    
    a <- all_sources[e]
    
    I <- Matrix::Diagonal(nrow(W))
    
    X_original <- raw_wide
    
    # change values of source at loop e to 0
    
    W_lesion <- W
    W_lesion[a, ] <- 0
    W_lesion[, a] <- 0
    
    B_lesion <- B
    B_lesion[, a] <- 0
    
    # loentief inverse
    
    X_lesion <- t(solve(I - t(W_lesion), t(B_lesion)))
    
    # result tibble 
    
    effector_loss_par <- tibble(
      lesion_id = e,
      source = a,
      condition = type_of_source[e],
      effector_group = rownames(X_original),
      original = rowSums(X_original),
      post_lesion = rowSums(X_lesion),
      absolute_loss = rowSums(X_original) - rowSums(X_lesion) # substraction
    )
    
    print(paste(e))
    
    effector_loss_par
    
  }
stopCluster(cl)

# First object to then join with others

effector_tibb <- effector_loss_par_res[[1]] %>%
  mutate(name = names(effector_loss_par_res)[1])

# Join with rest of result

for (i in 2:100) {
  effector_tibb <-
    rbind(effector_tibb, effector_loss_par_res[[i]] %>%
            mutate(name = names(effector_loss_par_res)[i]))
}

# Join with metrics for information

effector_comp <- effector_tibb %>% 
  left_join(
    metrics %>%
      select(source, type, dominant_effector, superclass),
    by = "source", multiple = "all"
  ) 


saveRDS(effector_comp, "effector_loss_dom.RDS")
