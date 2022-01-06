p <- 2
q <- 3
k <- 2

fm <- list(
  xcoef = matrix(1:4, p, k),
  ycoef = matrix(5:10, q, k)
)

vec <- c(1, 2, 5, 6, 7, 3, 4, 8, 9, 10)
mat <- matrix(vec, p + q, k)

test_that("Parameter reshaping functions give expected output", {
  expect_equal(fm2mat(fm), mat)
  expect_equal(fm2vec(fm), vec)
  expect_equal(mat2fm(mat, p), fm)
  expect_equal(mat2vec(mat), vec)
  expect_equal(vec2fm(vec, p, q), fm)
  expect_equal(vec2mat(vec, p, q), mat)

})

test_that("Parameter reshaping functions have well-defined inverses", {
  expect_equal(fm2mat(mat2fm(mat, p)), mat)
  expect_equal(fm2vec(vec2fm(vec, p, q)), vec)
  expect_equal(mat2fm(fm2mat(fm), p), fm)
  expect_equal(mat2vec(vec2mat(vec, p, q)), vec)
  expect_equal(vec2fm(fm2vec(fm), p, q), fm)
  expect_equal(vec2mat(mat2vec(mat), p, q), mat)
})
