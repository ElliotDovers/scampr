context("Combined Data Model Function: popa()")

test_that("ipp popa model creates scampr model", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp popa model creates scampr model without either intercept", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  mod <- popa(pres ~ -1 + MNT + D.Main, sp1 ~ -1 + MNT, dat_po, dat_pa, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp popa model doesn't work with PO Intercept and without PA intercept", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  testthat::expect_error(mod <- popa(pres ~ MNT + D.Main, sp1 ~ -1 + MNT, dat_po, dat_pa, model.type = "ipp"))
})

test_that("ipp popa model doesn't work with PA Intercept and without PO intercept", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  testthat::expect_error(mod <- popa(pres ~ -1 + MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "ipp"))
})

test_that("lgcp (va) popa model doesn't work", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  testthat::expect_error(popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "variational", simple.basis = bfs))
})

test_that("lgcp (laplace) popa model creates a scampr model", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "laplace", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})
