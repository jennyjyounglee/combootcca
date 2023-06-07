test_that("Greedy Cosine Similarity Alignment Works", {
  fm <- list()
  fm$xcoef <- matrix(c(0.9, 0.1, 0, -1, 0, 0, .5, 0, .5), 3, 3)
  fm$ycoef <- fm$xcoef
  fm$cor <- c(0.9, 0.8, 0.7)

  ref <- list()
  ref$xcoef <- diag(3)
  ref$ycoef <- diag(3)
  ref$cor <- c(0.9, 0.8, 0.7)

  map <- c(2, 1, 3)
  signs <- c(-1, 1, 1)
  fm_aligned <- fm
  fm_aligned$xcoef <- fm$xcoef[, map] %*% diag(signs)
  fm_aligned$ycoef <- fm$ycoef[, map] %*% diag(signs)
  fm_aligned$cor <- fm$cor[map]

  expect_equal(cca_align_greedy(fm, ref), fm_aligned)
})


test_that("Generalized Hungarian Algorithm Works: Square", {
  C <- diag(c(-10, 5, 1))
  P <- diag(c(-1, 1, 1))

  expect_equal(hungarian_max_signflip(C), P)

  C <- matrix(c(0, 5, 0, -10, 0, 0, 0, 0, 1), 3, 3)
  P <- matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 1), 3, 3)

  expect_equal(hungarian_max_signflip(C), P)
})

test_that("Generalized Hungarian Algorithm Works: Wide", {

  C <- matrix(
    c(0, -10, 3,
      5, 0, 25),
    byrow = TRUE,
    2, 3)

  P <- matrix(
    c(
      0, 0,
      -1, 0,
      0, 1
    ),
    byrow = TRUE,
    3, 2
  )

  expect_equal(hungarian_max_signflip(C), P)
})

test_that("Signflip Works", {
  fm <- list()
  fm$xcoef <- matrix(c(0.9, 0.1, 0, -0.9, 0, 0, .5, 0, -1.5), 3, 3)
  fm$ycoef <- diag(3)
  fm$cor <- c(0.9, 0.5, 0.1)

  ref <- list()
  ref$xcoef <- diag(3)
  ref$ycoef <- diag(3)
  ref$cor <- c(0.9, 0.5, 0.1)

  signs <- c(1, 1, -1)

  fm_aligned <- fm
  fm_aligned$xcoef <- fm$xcoef %*% diag(signs)
  fm_aligned$ycoef <- fm$ycoef %*% diag(signs)

  expect_equal(cca_align_signflip(fm, ref), fm_aligned)


})
