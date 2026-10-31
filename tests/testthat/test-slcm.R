test_that("slcm() stores every draw when burnin = 0", {
  set.seed(1)
  y <- matrix(rbinom(60 * 12, 1, 0.5), nrow = 60, ncol = 12)

  fit <- slcm(y, k = 2, burnin = 0, chain_length = 20)

  # Each stored draw of the class shares sums to one; unfilled columns are 0.
  expect_equal(unname(colSums(fit$chain$pi)), rep(1, 20))
})

test_that("slcm() rejects responses with missing values", {
  set.seed(1)
  y <- matrix(rbinom(60 * 12, 1, 0.5), nrow = 60, ncol = 12)
  y[1, 1] <- NA

  expect_error(
    slcm(y, k = 2, burnin = 5, chain_length = 10),
    "must not contain missing values"
  )
})
