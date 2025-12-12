#' K-Nearest Neighbor (KNN) Sampling from Candidates
#'
#' @description
#' Given a set of candidate vectors and a target vector, selects indices of the \code{k} nearest candidates
#' (using weighted Euclidean distance) to the target, then samples \code{n} indices from these neighbors,
#' either uniformly or with rank-based probabilities favoring closer neighbors. Optionally, a random seed can be set
#' for reproducibility.
#'
#' @param candidates Numeric matrix or data frame. Each row is a candidate vector.
#' @param target Numeric vector. The target vector for comparison.
#' @param k Integer. Number of nearest neighbors to consider (\eqn{k \leq} number of candidates).
#' @param n Integer. Number of samples to draw (default = 1).
#' @param prob Logical. If \code{TRUE}, sampling probabilities decrease with neighbor rank; if \code{FALSE}, sampling is uniform (default).
#' @param weights Optional numeric vector of length equal to \code{ncol(candidates)}. Feature weights for distance calculation.
#' @param seed Optional integer. If provided, sets random seed for reproducible sampling.
#'
#' @return
#' Integer vector of length \code{n}, giving indices of sampled candidates (rows of \code{candidates}).
#'
#' @details
#' The function computes the weighted Euclidean distance between each candidate and the target. The \code{k}
#' nearest neighbors (smallest distances) are identified. If \code{prob = TRUE}, the closest neighbor has the highest
#' probability (\code{1}), the second-closest half that, and so on (probabilities sum to 1). If \code{prob = FALSE},
#' all \code{k} neighbors are sampled with equal probability. Sampling is with replacement.
#' If \code{seed} is set, the random seed will be temporarily changed and restored on exit.
#'
#' @examples
#' set.seed(42)
#' candidates <- matrix(rnorm(50), ncol = 2)
#' target <- c(0, 0)
#'
#' # Sample 1 index from 5 nearest neighbors, uniform probability
#' knn_sample(candidates, target, k = 5, n = 1)
#'
#' # Sample 3 indices from 5 nearest neighbors, rank-weighted probability, with seed
#' knn_sample(candidates, target, k = 5, n = 3, prob = TRUE, seed = 123)
#'
#' # Using feature weights
#' knn_sample(candidates, target, k = 5, n = 2, weights = c(2, 1), seed = 10)
#'
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
    on.exit({
      if (exists("old_seed", inherits = FALSE)) {
        .Random.seed <<- old_seed
      }
    })
  }

  candidates <- as.matrix(candidates)

  p <- ncol(candidates)
  if (is.null(weights)) {
    weights <- rep(1, p)
  } else {
    if (length(weights) != p)
      stop("Length of weights must equal number of columns in candidates.")
  }

  # -------------------------------------------------
  # Weighted Euclidean distances
  # -------------------------------------------------
  diffs <- candidates - matrix(target, nrow(candidates), p, byrow = TRUE)
  weighted_sq_diffs <- diffs^2 * rep(weights, each = nrow(candidates))
  dists <- sqrt(rowSums(weighted_sq_diffs))

  # -------------------------------------------------
  # k nearest neighbors
  # -------------------------------------------------
  ord <- order(dists)
  nn_indices <- ord[seq_len(min(k, length(ord)))]
  nn_dists   <- dists[nn_indices]
  k_eff      <- length(nn_indices)

  # -------------------------------------------------
  # Sampling probabilities
  # -------------------------------------------------
  if (!prob) {

    probs <- rep(1 / k_eff, k_eff)

  } else if (sampling == "rank") {

    # Rank-based probabilities (default)
    probs <- (1 / seq_len(k_eff))
    probs <- probs / sum(probs)

  } else if (sampling == "distance") {

    # Distance-based probabilities
    if (is.null(bandwidth)) {
      # Automatic bandwidth: median NN distance (robust)
      bandwidth <- stats::median(nn_dists, na.rm = TRUE)
    }

    if (!is.finite(bandwidth) || bandwidth <= 0) {
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

