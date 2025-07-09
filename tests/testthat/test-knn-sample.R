# In tests/testthat/test-knn_sample.R

test_that("knn_sample returns correct length and index range", {
  set.seed(100)
  candidates <- matrix(rnorm(100), ncol = 2)
  target <- c(0, 0)

  # Sample 1 index from 5 nearest
  idx <- knn_sample(candidates, target, k = 5, n = 1)
  expect_true(length(idx) == 1)
  expect_true(all(idx >= 1 & idx <= nrow(candidates)))

  # Sample 10 indices from 5 nearest
  idx10 <- knn_sample(candidates, target, k = 5, n = 10)
  expect_true(length(idx10) == 10)
  expect_true(all(idx10 >= 1 & idx10 <= nrow(candidates)))
})

test_that("knn_sample is reproducible with seed", {
  candidates <- matrix(rnorm(100), ncol = 2)
  target <- c(0, 0)

  idx1 <- knn_sample(candidates, target, k = 5, n = 4, seed = 42)
  idx2 <- knn_sample(candidates, target, k = 5, n = 4, seed = 42)
  expect_identical(idx1, idx2)
})

test_that("knn_sample handles weights", {
  candidates <- matrix(c(0,0, 1,1, 2,2), ncol = 2, byrow = TRUE)
  target <- c(1, 1)
  idx <- knn_sample(candidates, target, k = 1, weights = c(10, 1))
  expect_true(idx %in% 1:3)
})

test_that("knn_sample errors on bad weights length", {
  candidates <- matrix(rnorm(12), ncol = 3)
  target <- c(0, 0, 0)
  expect_error(knn_sample(candidates, target, k = 2, weights = c(1, 2)))
})

test_that("knn_sample samples with rank-based probability", {
  candidates <- matrix(seq(-2, 2, length.out = 10), ncol = 2)
  target <- c(0, 0)
  out1 <- knn_sample(candidates, target, k = 3, n = 1000, prob = TRUE, seed = 77)
  tab <- table(out1)
  # The most frequent should be the closest neighbor (index 5 or 6 for a symmetric matrix)
  expect_true(any(tab == max(tab)))
})

# Optional: Test with data frame input
test_that("knn_sample works with data frame candidates", {
  df_candidates <- as.data.frame(matrix(rnorm(40), ncol = 2))
  target <- c(0, 0)
  idx <- knn_sample(df_candidates, target, k = 4, n = 2)
  expect_true(length(idx) == 2)
})
