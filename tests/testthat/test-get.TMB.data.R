context("Main internal function for modelling framework get.TMB.data.input()")

test_that("gets correct fixed effect names", {
  dat_po <- eucalypt[["po"]]
  dat_pa <- eucalypt[["pa"]]
  std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
  dat_po$mnt <- std.tmp[1:nrow(dat_po)]
  dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
  dat_po$rd <- scale(dat_po$D_MAIN_RDS)
  res.list <- get.TMB.data.input(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, data.type = "popa")
  testthat::expect_equal(res.list$fixed.names, c("(Intercept)", "mnt", "(Bias Intercept)", "rd"))
})
