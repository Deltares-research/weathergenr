#' K-Nearest Neighbor (KNN) Sampling from Candidates
#'
#' @description
#' Given a set of candidate vectors and a target vector, selects indices of the
#' \code{k} nearest candidates (using weighted Euclidean distance) to the target,
#' then samples \code{n} indices from these neighbors, either uniformly or with
#' rank-based or distance-based probabilities. Optionally, a random seed can be
#' set for reproducibility.
#'
#' @param candidates Numeric matrix or data frame. Each row is a candidate vector.
#' @param target Numeric vector. The target vector for comparison.
#' @param k Integer. Number of nearest neighbors to consider (k <= number of candidates).
#' @param n Integer. Number of samples to draw (default = 1).
#' @param prob Logical. If TRUE, sampling probabilities favor closer neighbors;
#'   if FALSE (default), sampling is uniform among k nearest neighbors.
#' @param weights Optional numeric vector of length equal to ncol(candidates).
#'   Feature weights for distance calculation. Default is equal weights.
#' @param seed Optional integer. If provided, sets random seed for reproducible
#'   sampling. The seed state is restored on exit.
#' @param sampling Character. Sampling method when prob = TRUE. Either "rank"
#'   (default, probability decreases with neighbor rank: 1, 1/2, 1/3, ...) or
#'   "distance" (probability based on Gaussian kernel of distances).
#'   Ignored when prob = FALSE.
#' @param bandwidth Numeric. Bandwidth parameter for distance-based sampling
#'   kernel. If NULL (default), uses median nearest neighbor distance.
#'   Only used when sampling = "distance" and prob = TRUE.
#' @param epsilon Numeric. Small constant added to distance-based probabilities
#'   to prevent zero probabilities (default = 1e-8). Only used when
#'   sampling = "distance" and prob = TRUE.
#'
#' @return
#' Integer vector of length \code{n}, giving indices of sampled candidates
#' (rows of \code{candidates}).
#'
#' @details
#' The function computes the weighted Euclidean distance between each candidate
#' and the target. The \code{k} nearest neighbors (smallest distances) are
#' identified.
#'
#' Sampling modes:
#' \describe{
#'   \item{prob = FALSE}{All k neighbors sampled with equal probability (uniform)}
#'   \item{prob = TRUE, sampling = "rank"}{Probability = 1/rank, normalized.
#'         Closest neighbor has highest probability.}
#'   \item{prob = TRUE, sampling = "distance"}{Probability based on Gaussian
#'         kernel: exp(-distance^2 / (2 * bandwidth^2))}
#' }
#'
#' Sampling is with replacement. If \code{seed} is set, the random seed will be
#' temporarily changed and restored on exit.
#'
#' @examples
#' set.seed(42)
#' candidates <- matrix(rnorm(50), ncol = 2)
#' target <- c(0, 0)
#'
#' # Sample 1 index from 5 nearest neighbors, uniform probability
#' knn_sample(candidates, target, k = 5, n = 1)
#'
#' # Sample 3 indices from 5 nearest neighbors, rank-weighted probability
#' knn_sample(candidates, target, k = 5, n = 3, prob = TRUE, seed = 123)
#'
#' # Using feature weights (weight first dimension more heavily)
#' knn_sample(candidates, target, k = 5, n = 2, weights = c(2, 1), seed = 10)
#'
#' # Distance-based sampling with custom bandwidth
#' knn_sample(candidates, target, k = 10, n = 5, prob = TRUE,
#'            sampling = "distance", bandwidth = 1.5)
#' @export
knn_sample <- function(
    candidates,
    target,
    k,
    n = 1,
    prob = FALSE,
    weights = NULL,
    seed = NULL,
    sampling = c("rank", "distance"),
    bandwidth = NULL,
    epsilon = 1e-8
) {

  sampling <- match.arg(sampling)

  # -------------------------------------------------
  # RNG handling
  # -------------------------------------------------
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    set.seed(seed)
    on.exit(.Random.seed <<- old_seed, add = TRUE)
  }

  candidates <- as.matrix(candidates)
  nc <- nrow(candidates)
  p <- ncol(candidates)
  if (nc == 0) {
    stop("No candidates provided to knn_sample")
  }

  if (is.null(weights)) {
    weights <- rep(1, p)
  } else {
    if (length(weights) != p)
      stop("Length of weights must equal number of columns in candidates.")
  }

  # -------------------------------------------------
  # OPTIMIZED: Weighted squared Euclidean distances
  # Compute squared distances (defer sqrt until needed)
  # -------------------------------------------------
  if (p == 2) {
    # Special case for p=2 (common: precipitation + temperature)
    # Unrolled for maximum speed
    d1 <- candidates[, 1] - target[1]
    d2 <- candidates[, 2] - target[2]
    d2_sq <- weights[1] * d1 * d1 + weights[2] * d2 * d2

  } else if (p <= 5) {
    # Vectorized operations for small p
    diffs <- sweep(candidates, 2, target, "-")
    weighted_sq <- sweep(diffs^2, 2, weights, "*")
    d2_sq <- rowSums(weighted_sq)

  } else {
    # Loop for large p (memory efficient)
    d2_sq <- numeric(nc)
    for (j in seq_len(p)) {
      dj <- candidates[, j] - target[j]
      d2_sq <- d2_sq + weights[j] * dj * dj
    }
  }

  # -------------------------------------------------
  # OPTIMIZED: Partial sorting when k << n
  # -------------------------------------------------
  k_eff <- min(k, nc)

  if (k_eff < nc * 0.2) {
    # Partial sorting: faster when k is small relative to n
    # Uses sort.int with partial parameter: O(n) + O(k log k) instead of O(n log n)

    # Find the k-th smallest squared distance (partition threshold)
    threshold <- sort(d2_sq, partial = k_eff)[k_eff]

    # Get all candidate indices with squared distance <= threshold
    # (may include ties at the boundary)
    candidates_idx <- which(d2_sq <= threshold)

    # Fully sort only this subset to get exact k nearest
    sorted_subset_order <- order(d2_sq[candidates_idx])
    nn_indices <- candidates_idx[sorted_subset_order[seq_len(k_eff)]]

  } else {
    # k is large relative to n: full sort is fine
    nn_indices <- order(d2_sq)[seq_len(k_eff)]
  }

  # OPTIMIZED: Only compute sqrt for the k nearest neighbors
  # (avoided n - k sqrt operations)
  nn_dists <- sqrt(d2_sq[nn_indices])

  # -------------------------------------------------
  # Sampling probabilities
  # -------------------------------------------------
  if (!prob) {
    # Uniform sampling
    probs <- rep(1 / k_eff, k_eff)

  } else if (sampling == "rank") {
    # Rank-based probabilities (default for prob = TRUE)
    probs <- (1 / seq_len(k_eff))
    probs <- probs / sum(probs)

  } else if (sampling == "distance") {
    # Distance-based probabilities (Gaussian kernel)
    if (is.null(bandwidth)) {
      # Automatic bandwidth: median NN distance
      bandwidth <- median(nn_dists, na.rm = TRUE)
    }

    if (!is.finite(bandwidth) || bandwidth <= 0) {
      # Fallback to uniform if bandwidth invalid
      probs <- rep(1 / k_eff, k_eff)
    } else {
      probs <- exp(-(nn_dists^2) / (2 * bandwidth^2)) + epsilon
      probs <- probs / sum(probs)
    }
  }

  # -------------------------------------------------
  # Sample neighbors
  # -------------------------------------------------
  sampled_rel <- sample.int(k_eff, n, replace = TRUE, prob = probs)
  nn_indices[sampled_rel]
}

