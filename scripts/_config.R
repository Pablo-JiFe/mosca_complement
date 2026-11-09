# Global analysis constants. Environment variables can override the expensive
# permutation counts without changing the committed analysis definition.

RNG_SEED <- 1L

PRIMARY_MIN_SYNAPSES <- 5L
PRIMARY_RHO <- 0.90
PRIMARY_EFFECTOR_GROUP <- "effector_group_fine"

RHO_GRID <- c(0.50, 0.70, 0.90, 0.95, 0.99)
SYNAPSE_GRID <- c(1L, 3L, 5L, 10L)
EFFECTOR_GROUP_GRID <- c(
  "effector_group_fine",
  "effector_group_coarse",
  "effector_group_superclass"
)

HORN_PERMUTATIONS <- as.integer(Sys.getenv("HORN_PERMUTATIONS", "200"))
DN_PERMUTATIONS <- as.integer(Sys.getenv("DN_PERMUTATIONS", "10000"))
DEGENERACY_PERMUTATIONS <- as.integer(Sys.getenv("DEGENERACY_PERMUTATIONS", "9999"))
DEGENERACY_NULL_SAMPLE <- as.integer(Sys.getenv("DEGENERACY_NULL_SAMPLE", "3000"))
MIN_GROUP_N <- 20L
PCA_MAX_RANK <- 50L
ZERO_TOL <- 1e-12
NEGATIVE_TOL <- 1e-8

UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.05

.detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
N_THREADS <- if (is.na(.detected_cores)) 1L else max(1L, .detected_cores - 1L)
rm(.detected_cores)

stopifnot(
  PRIMARY_MIN_SYNAPSES %in% SYNAPSE_GRID,
  PRIMARY_RHO > 0,
  PRIMARY_RHO < 1,
  all(RHO_GRID > 0 & RHO_GRID < 1),
  HORN_PERMUTATIONS >= 20L,
  DN_PERMUTATIONS >= 100L,
  DEGENERACY_PERMUTATIONS >= 19L,
  DEGENERACY_NULL_SAMPLE >= 100L
)
