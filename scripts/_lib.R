first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0L) NA_character_ else x[[1L]]
}

build_weight_matrix <- function(edges, units, rho) {
  stopifnot(
    is.numeric(rho), length(rho) == 1L, rho > 0, rho < 1,
    !anyDuplicated(units),
    all(c("pre", "post", "count") %in% names(edges))
  )

  unit_index <- stats::setNames(seq_along(units), units)

  e <- edges |>
    dplyr::filter(pre %in% units, post %in% units, count > 0) |>
    dplyr::mutate(
      i = unname(unit_index[post]),
      j = unname(unit_index[pre])
    )

  stopifnot(nrow(e) > 0L, !anyNA(e$i), !anyNA(e$j))

  # W[target, source]. Before spectral rescaling, each row is normalized by
  # total synaptic input to the postsynaptic unit.
  W_count <- Matrix::sparseMatrix(
    i = e$i,
    j = e$j,
    x = e$count,
    dims = c(length(units), length(units)),
    dimnames = list(units, units)
  )

  input_count <- rowSums(W_count)
  input_scale <- ifelse(input_count > 0, 1 / input_count, 0)
  W_base <- Matrix::Diagonal(x = input_scale) %*% W_count

  leading <- RSpectra::eigs(W_base, k = 1L, which = "LM")$values[[1L]]
  lambda <- Mod(leading)

  stopifnot(is.finite(lambda), lambda > 0)

  W <- W_base * (rho / lambda)
  dimnames(W) <- list(units, units)

  list(
    W = W,
    lambda_before_rescale = lambda,
    leading_eigenvalue = leading,
    n_zero_input_units = sum(input_count == 0),
    n_edges = nrow(e)
  )
}

build_effector_matrix <- function(effectors, units, group_col) {
  stopifnot(
    group_col %in% names(effectors),
    all(c("unit", "n_neurons") %in% names(effectors))
  )

  unit_index <- stats::setNames(seq_along(units), units)

  e <- effectors |>
    dplyr::filter(
      unit %in% units,
      !is.na(.data[[group_col]]),
      .data[[group_col]] != "",
      n_neurons > 0
    ) |>
    dplyr::transmute(
      unit,
      effector_group = .data[[group_col]],
      n_neurons
    ) |>
    dplyr::group_by(unit, effector_group) |>
    dplyr::summarise(n_neurons = sum(n_neurons), .groups = "drop")

  stopifnot(nrow(e) > 0L)

  groups <- sort(unique(e$effector_group))
  group_index <- stats::setNames(seq_along(groups), groups)

  weights <- e |>
    dplyr::group_by(effector_group) |>
    dplyr::mutate(weight = n_neurons / sum(n_neurons)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      i = unname(group_index[effector_group]),
      j = unname(unit_index[unit])
    )

  stopifnot(!anyNA(weights$i), !anyNA(weights$j))

  B <- Matrix::sparseMatrix(
    i = weights$i,
    j = weights$j,
    x = weights$weight,
    dims = c(length(groups), length(units)),
    dimnames = list(groups, units)
  )

  group_weight_sums <- rowSums(B)
  stopifnot(max(abs(group_weight_sums - 1)) < 1e-10)

  list(B = B, groups = groups)
}

compute_perturbome <- function(
  edges,
  meta,
  effectors,
  rho,
  group_col,
  negative_tol = 1e-8
) {
  units <- meta$unit
  stopifnot(length(units) > 1L, !anyDuplicated(units))

  wm <- build_weight_matrix(edges, units, rho)
  em <- build_effector_matrix(effectors, units, group_col)

  A <- Matrix::Diagonal(length(units)) - wm$W

  # P = B %*% solve(I - W). Solving the transposed system needs one sparse
  # factorization with g right-hand sides instead of n separate solves.
  X <- solve(t(A), t(em$B))
  P_raw <- t(as.matrix(X))

  # solve()/as.matrix() can discard one or both dimnames. Restore them
  # explicitly before any serialization.
  dimnames(P_raw) <- list(em$groups, units)

  P_indirect <- P_raw - as.matrix(em$B)
  dimnames(P_indirect) <- dimnames(P_raw)

  stopifnot(all(is.finite(P_raw)), all(is.finite(P_indirect)))

  min_indirect <- min(P_indirect)
  if (min_indirect < -negative_tol) {
    stop(
      "Indirect influence contains values below numerical tolerance: ",
      signif(min_indirect, 6)
    )
  }

  # Only clip numerical noise; preserve all genuine positive influence.
  P_indirect[P_indirect < 0] <- 0

  stopifnot(
    identical(rownames(P_raw), em$groups),
    identical(colnames(P_raw), units),
    identical(dimnames(P_indirect), dimnames(P_raw))
  )

  list(
    raw = P_raw,
    indirect = P_indirect,
    lambda_before_rescale = wm$lambda_before_rescale,
    leading_eigenvalue = wm$leading_eigenvalue,
    n_zero_input_units = wm$n_zero_input_units,
    n_edges = wm$n_edges,
    n_effector_groups = length(em$groups)
  )
}


compute_perturbome_iterative <- function(
  edges,
  meta,
  effectors,
  rho,
  group_col,
  negative_tol = 1e-8,
  tol = 1e-10,
  maxit = 5000L
) {
  units <- meta$unit
  stopifnot(
    length(units) > 1L,
    !anyDuplicated(units),
    is.numeric(tol), length(tol) == 1L, tol > 0,
    is.numeric(maxit), length(maxit) == 1L, maxit >= 1L
  )

  wm <- build_weight_matrix(edges, units, rho)
  em <- build_effector_matrix(effectors, units, group_col)

  B <- as.matrix(em$B)
  dimnames(B) <- list(em$groups, units)

  # Same equation as compute_perturbome():
  #   P = B %*% solve(I - W)
  # hence
  #   P = B + P %*% W.
  #
  # Because build_weight_matrix() rescales W to spectral radius rho < 1,
  # this fixed-point iteration converges to the same unique solution without
  # constructing a sparse LU factorization (and its potentially large fill-in).
  P_raw <- B
  residual <- Inf

  for (iter in seq_len(maxit)) {
    P_next <- B + as.matrix(P_raw %*% wm$W)
    dimnames(P_next) <- dimnames(B)

    # P_next - P_raw is exactly the fixed-point residual:
    # B + P_raw W - P_raw.
    residual <- max(abs(P_next - P_raw)) /
      max(1, max(abs(P_next)))

    P_raw <- P_next

    if (is.finite(residual) && residual <= tol) break
  }

  if (!is.finite(residual) || residual > tol) {
    stop(
      "Iterative perturbome solver did not converge; residual = ",
      signif(residual, 6),
      ", iterations = ",
      maxit
    )
  }

  P_indirect <- P_raw - B
  dimnames(P_indirect) <- dimnames(P_raw)

  stopifnot(all(is.finite(P_raw)), all(is.finite(P_indirect)))

  min_indirect <- min(P_indirect)
  if (min_indirect < -negative_tol) {
    stop(
      "Indirect influence contains values below numerical tolerance: ",
      signif(min_indirect, 6)
    )
  }

  P_indirect[P_indirect < 0] <- 0

  stopifnot(
    identical(rownames(P_raw), em$groups),
    identical(colnames(P_raw), units),
    identical(dimnames(P_indirect), dimnames(P_raw))
  )

  list(
    raw = P_raw,
    indirect = P_indirect,
    lambda_before_rescale = wm$lambda_before_rescale,
    leading_eigenvalue = wm$leading_eigenvalue,
    n_zero_input_units = wm$n_zero_input_units,
    n_edges = wm$n_edges,
    n_effector_groups = length(em$groups),
    iterations = iter,
    residual = residual
  )
}

matrix_to_long <- function(x) {
  stopifnot(!is.null(rownames(x)), !is.null(colnames(x)))

  as.data.frame(x, check.names = FALSE) |>
    tibble::rownames_to_column("effector_group") |>
    tidyr::pivot_longer(
      cols = -effector_group,
      names_to = "source",
      values_to = "influence"
    )
}

signature_matrix_from_long <- function(p, source_order = NULL) {
  stopifnot(all(c("source", "effector_group", "influence") %in% names(p)))

  if (!is.null(source_order)) {
    p <- p |>
      dplyr::mutate(source = factor(source, levels = source_order)) |>
      dplyr::arrange(source, effector_group)
  }

  wide <- p |>
    tidyr::pivot_wider(
      names_from = effector_group,
      values_from = influence,
      values_fill = 0
    )

  sources <- as.character(wide$source)
  X <- wide |>
    dplyr::select(-source) |>
    as.matrix()

  rownames(X) <- sources
  list(X = X, sources = sources)
}

normalize_signatures <- function(X, zero_tol = 1e-12) {
  totals <- rowSums(X)
  active <- totals > zero_tol

  prob <- matrix(0, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
  prob[active, ] <- X[active, , drop = FALSE] / totals[active]

  list(total = totals, active = active, prob = prob)
}

pca_eigenvalues <- function(X) {
  stopifnot(nrow(X) >= 2L, ncol(X) >= 2L)
  Xc <- sweep(X, 2L, colMeans(X), FUN = "-")
  vals <- eigen(
    crossprod(Xc) / (nrow(Xc) - 1),
    symmetric = TRUE,
    only.values = TRUE
  )$values
  pmax(vals, 0)
}

pc_count_for_fraction <- function(eigenvalues, fraction) {
  stopifnot(fraction > 0, fraction < 1)
  total <- sum(eigenvalues)
  if (total <= 0) return(NA_integer_)
  which(cumsum(eigenvalues / total) >= fraction)[[1L]]
}

nearest_distinct_type <- function(prob_active, sources_active, meta) {
  stopifnot(nrow(prob_active) == length(sources_active))

  if (nrow(prob_active) < 2L) {
    return(tibble::tibble(
      source = character(),
      nearest_source = character(),
      nearest_cosine = numeric()
    ))
  }

  type_active <- meta$type[match(sources_active, meta$unit)]
  stopifnot(!anyNA(type_active))

  l2 <- sqrt(rowSums(prob_active^2))
  X <- prob_active / l2

  # At most m-1 other rows can share the source type, where m is the largest
  # multiplicity of any type. We ask for m+1 candidates and explicitly reject
  # self matches, which is safe even if a backend changes self-neighbour rules.
  max_type_multiplicity <- max(table(type_active))
  k <- min(
    nrow(X) - 1L,
    max(1L, as.integer(max_type_multiplicity) + 1L)
  )

  if (k < 1L || length(unique(type_active)) < 2L) {
    return(tibble::tibble(
      source = sources_active,
      nearest_source = NA_character_,
      nearest_cosine = NA_real_
    ))
  }

  nn <- FNN::get.knn(X, k = k)

  pick_one <- function(row_id) {
    candidates <- nn$nn.index[row_id, ]
    keep <-
      candidates != row_id &
      type_active[candidates] != type_active[[row_id]]
    pos <- which(keep)[1L]

    if (is.na(pos)) {
      return(c(index = NA_real_, cosine = NA_real_))
    }

    dist <- nn$nn.dist[row_id, pos]
    c(
      index = candidates[[pos]],
      cosine = max(-1, min(1, 1 - dist^2 / 2))
    )
  }

  picked <- purrr::map(seq_len(nrow(X)), pick_one)

  tibble::tibble(
    source = sources_active,
    nearest_source = purrr::map_chr(
      picked,
      \(z) if (is.na(z[["index"]])) NA_character_ else sources_active[[as.integer(z[["index"]])]]
    ),
    nearest_cosine = purrr::map_dbl(picked, \(z) z[["cosine"]])
  )
}

summarize_perturbome_matrix <- function(P_indirect, meta, zero_tol = 1e-12) {
  # P is effector_group x source; signatures are source x effector_group.
  X <- t(P_indirect)
  normed <- normalize_signatures(X, zero_tol)

  active_prob <- normed$prob[normed$active, , drop = FALSE]
  active_sources <- rownames(active_prob)

  k <- ncol(active_prob)
  entropy <- -rowSums(ifelse(active_prob > 0, active_prob * log(active_prob), 0))
  breadth <- if (k > 1L) entropy / log(k) else rep(0, nrow(active_prob))
  hill_q2 <- 1 / rowSums(active_prob^2)

  eigenvalues <- pca_eigenvalues(active_prob)
  nn <- nearest_distinct_type(active_prob, active_sources, meta)

  tibble::tibble(
    n_sources = nrow(X),
    n_active = sum(normed$active),
    active_fraction = mean(normed$active),
    n_effector_groups = ncol(X),
    median_total_influence = stats::median(normed$total[normed$active]),
    median_breadth = stats::median(breadth),
    median_effective_channels = stats::median(hill_q2),
    median_nearest_distinct_type_cosine = stats::median(nn$nearest_cosine, na.rm = TRUE),
    pc50 = pc_count_for_fraction(eigenvalues, 0.50),
    pc75 = pc_count_for_fraction(eigenvalues, 0.75),
    pc90 = pc_count_for_fraction(eigenvalues, 0.90)
  )
}
