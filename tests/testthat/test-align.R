test_that("Greedy Cosine Similarity Alignment Works", {
  fm <- list()
  fm$xcoef <- matrix(c(0.9, 0.1, 0, -1, 0, 0, .5, 0, .5), 3, 3)
  fm$ycoef <- diag(3)

  ref <- list()
  ref$xcoef <- diag(3)
  ref$ycoef <- diag(3)

  map <- c(2, 1, 3)
  signs <- c(-1, 1, 1)
  fm_aligned <- fm
  fm_aligned$xcoef <- fm$xcoef[, map] %*% diag(signs)
  fm_aligned$ycoef <- fm$ycoef[, map] %*% diag(signs)

  expect_equal(cca_align_greedy_cosx(fm, ref), fm_aligned)
})