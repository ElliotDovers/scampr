context("simulate_model_pp function")

# test_that("ipp po model works from fitted quad", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   mod <- po(pres ~ elev, dat, model.type = "ipp")
#   pp <- simulate(mod)
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("ipp popa model works from fitted quad", {
#   dat_po <- eucalypt[["po"]]
#   dat_pa <- eucalypt[["pa"]]
#   std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
#   dat_po$mnt <- std.tmp[1:nrow(dat_po)]
#   dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
#   dat_po$rd <- scale(dat_po$D_MAIN_RDS)
#   mod <- popa(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, model.type = "ipp")
#   pp <- simulate(mod)
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("ipp po model works from provided quad", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   mod <- po(pres ~ elev, dat, model.type = "ipp")
#   pp <- simulate(mod, domain.data = dat[dat$pres == 0, ])
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("ipp popa model works from provided quad", {
#   dat_po <- eucalypt[["po"]]
#   dat_pa <- eucalypt[["pa"]]
#   std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
#   dat_po$mnt <- std.tmp[1:nrow(dat_po)]
#   dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
#   dat_po$rd <- scale(dat_po$D_MAIN_RDS)
#   mod <- popa(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, model.type = "ipp")
#   pp <- simulate(mod, domain.data = dat_po[dat_po$pres == 0, ])
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("ipp pa model does not work", {
#   dat_pa <- eucalypt[["pa"]]
#   dat_pa$mnt <- scale(dat_pa$TMP_MIN)
#   mod <- pa(Y ~ mnt, dat_pa, model.type = "ipp")
#   testthat::expect_error(simulate(mod))
# })
#
# test_that("lgcp (va) po model works from fitted quad", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
#   mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
#   pp <- simulate(mod)
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("lgcp (laplace) po model works from fitted quad", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
#   mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
#   pp <- simulate(mod)
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("lgcp popa model works from fitted quad", {
#   dat_po <- eucalypt[["po"]]
#   dat_pa <- eucalypt[["pa"]]
#   std.tmp <- scale(c(dat_po$TMP_MIN, dat_pa$TMP_MIN))
#   dat_po$mnt <- std.tmp[1:nrow(dat_po)]
#   dat_pa$mnt <- std.tmp[(nrow(dat_po) + 1):length(std.tmp)]
#   dat_po$rd <- scale(dat_po$D_MAIN_RDS)
#   bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
#   mod <- popa(pres ~ mnt + rd, Y ~ mnt, dat_po, dat_pa, model.type = "laplace", simple.basis = bfs)
#   pp <- simulate(mod)
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("lgcp (laplace) po model works from fitted quad w random coeffs from posterior", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
#   mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
#   pp <- simulate(mod)
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("lgcp (va) po model works from fitted quad w random coeffs from posterior", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
#   mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
#   pp <- simulate(mod, rcoef.density = "posterior", which.intensity = "sample")
#   testthat::expect_equal(class(pp), "ppp")
# })
#
# test_that("lgcp (va) po model works from fitted quad w random coeffs from prior", {
#   dat <- scampr::gorillas
#   dat$elev <- scale(dat$elevation)
#   bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
#   mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
#   pp <- simulate(mod, rcoef.density = "prior", which.intensity = "sample")
#   testthat::expect_equal(class(pp), "ppp")
# })
