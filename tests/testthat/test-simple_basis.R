context("simple_basis function")

test_that("square domain works", {
  # set a square domain
  dat <- expand.grid(list(x = seq(0, 100, by = 10), y = seq(0, 100, by = 10)))
  # get basis nodes
  bfs <- simple_basis(10, dat)

  # check number of nodes on horizontal axis equals 10
  expect_equal(length(unique(bfs[, 1L])), 10)

  # check number of nodes on vertical axis equals 10
  expect_equal(length(unique(bfs[, 2L])), 10)
})

test_that("tall domain works", {
  # set a tall domain
  dat <- expand.grid(list(x= seq(0, 50, by = 10), y= seq(0, 100, by = 10)))
  # get basis nodes
  bfs <- simple_basis(10, dat)

  # check number of nodes on horizontal axis equals 5
  expect_equal(length(unique(bfs[, 1L])), 5)

  # check number of nodes on vertical axis equals 10
  expect_equal(length(unique(bfs[, 2L])), 10)
})

test_that("wide domain works", {
  # set a tall domain
  dat <- expand.grid(list(x= seq(0, 100, by = 10), y= seq(0, 50, by = 10)))
  # get basis nodes
  bfs <- simple_basis(10, dat)

  # check number of nodes on horizontal axis equals 5
  expect_equal(length(unique(bfs[, 1L])), 10)

  # check number of nodes on vertical axis equals 10
  expect_equal(length(unique(bfs[, 2L])), 5)
})
