context("Main internal function for modelling framework get.TMB.data.input()")

test_that("gets correct fixed effect names", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  res.list <- get.TMB.data.input(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, data.type = "popa")
  testthat::expect_equal(res.list$fixed.names, c("(Intercept)", "MNT", "(Bias Intercept)", "D.Main"))
})
