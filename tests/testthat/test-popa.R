context("Combined Data Model Function: popa()")

test_that("ipp popa model creates scampr model", {
  dat_po <- eucalypt[["po"]]
  dat_pa <- eucalypt[["pa"]]
  std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
  dat_po$mnt <- std.tmp[1:nrow(dat_po)]
  dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
  dat_po$rd <- scale(dat_po$D_MAIN_RDS)
  mod <- popa(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (va) popa model doesn't work", {
  dat_po <- eucalypt[["po"]]
  dat_pa <- eucalypt[["pa"]]
  std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
  dat_po$mnt <- std.tmp[1:nrow(dat_po)]
  dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
  dat_po$rd <- scale(dat_po$D_MAIN_RDS)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  testthat::expect_error(popa(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, model.type = "variational", simple.basis = bfs))
})

test_that("lgcp (laplace) popa model creates a scampr model", {
  dat_po <- eucalypt[["po"]]
  dat_pa <- eucalypt[["pa"]]
  std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
  dat_po$mnt <- std.tmp[1:nrow(dat_po)]
  dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
  dat_po$rd <- scale(dat_po$D_MAIN_RDS)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  mod <- popa(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, model.type = "laplace", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})
