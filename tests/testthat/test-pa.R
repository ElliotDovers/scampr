context("Presence/Absence Model Function: pa()")
# usethis::use_test() when particular R function file is open to create these

test_that("ipp pa model creates scampr model", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  mod <- pa(Y ~ mnt, dat_pa, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (va) pa model doesn't work", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_pa)
  testthat::expect_error(pa(Y ~ mnt, dat_pa, model.type = "variational", simple.basis = bfs))
})

test_that("lgcp (laplace) pa model creates scampr model", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_pa)
  mod <- pa(Y ~ mnt, dat_pa, model.type = "laplace", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp pa model creates scampr model w. subset", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  mod <- pa(Y ~ mnt, dat_pa, model.type = "ipp", subset = 1:nrow(dat_pa) < (0.5 * nrow(dat_pa)))
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp pa model creates scampr model works with interaction terms in formula", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  dat_pa$mxt <- scale(dat_pa$TMP_MAX)
  mod <- pa(Y ~ mnt*mxt, dat_pa, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (laplace) pa model works with interaction terms in formula", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  dat_pa$mxt <- scale(dat_pa$TMP_MAX)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_pa)
  mod <- pa(Y ~ mnt*mxt, dat_pa, model.type = "laplace", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (laplace) pa model works w. subset", {
  dat_pa <- eucalypt[["pa"]]
  dat_pa$mnt <- scale(dat_pa$TMP_MIN)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_pa)
  mod <- pa(Y ~ mnt, dat_pa, model.type = "laplace", simple.basis = bfs, subset = 1:nrow(dat_pa) < (0.5 * nrow(dat_pa)))
  testthat::expect_equal(class(mod), "scampr")
})
