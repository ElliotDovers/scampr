context("simulate function")

test_that("ipp po model works from fitted quad", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp")
  pp <- simulate(mod, return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("ipp popa model works from fitted quad", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "ipp")
  pp <- simulate(mod, return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("ipp po model works from provided quad", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp")
  pp <- simulate(mod, domain.data = dat[dat$pres == 0, ], return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("ipp popa model can return data frame", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "ipp")
  pp <- simulate(mod, return.type = "data.frame")
  testthat::expect_equal(class(pp), "data.frame")
})

test_that("ipp pa model does not work", {
  dat_pa <- flora$pa
  mod <- pa(sp1 ~ MNT, dat_pa, model.type = "ipp")
  testthat::expect_error(simulate(mod))
})

test_that("lgcp (va) po model works with expected posterior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "posterior", which.intensity = "expected", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (laplace) po model works with expected posterior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "posterior", which.intensity = "expected", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (va) po model works with expected prior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "prior", which.intensity = "expected", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (laplace) po model works with expected prior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "prior", which.intensity = "expected", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (va) po model works with sampled posterior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "posterior", which.intensity = "sample", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (laplace) po model works with sampled posterior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "posterior", which.intensity = "sample", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (va) po model works with sampled prior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "prior", which.intensity = "sample", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp (laplace) po model works with sampled prior intensity", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, rcoef.density = "posterior", which.intensity = "sample", return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("ipp model works with multiple sims", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp")
  pp <- simulate(mod, nsim = 3, return.type = "ppp")
  testthat::expect_equal(class(pp),  c("ppplist", "solist", "anylist", "listof", "list"))
})

test_that("lgcp (va) po model works with multiple sims (sampled intensity)", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, which.intensity = "sample", nsim = 3, return.type = "ppp")
  testthat::expect_equal(class(pp), c("ppplist", "solist", "anylist", "listof", "list"))
})

test_that("ipp model can return data frame", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp")
  pp <- simulate(mod, return.type = "data.frame")
  testthat::expect_equal(class(pp), "data.frame")
})

test_that("ipp model can return data frame from multiple sims", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp")
  pp <- simulate(mod, return.type = "data.frame", nsim = 3)
  testthat::expect_equal(class(pp[[1L]]), "data.frame")
})

test_that("lgcp (va) po model can return data frame (sampled intensity)", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, which.intensity = "sample", return.type = "data.frame")
  testthat::expect_equal(class(pp), "data.frame")
})

test_that("lgcp (va) po model can return data frame with multiple sims (sampled intensity)", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 3, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  pp <- simulate(mod, which.intensity = "sample", nsim = 3, return.type = "data.frame")
  testthat::expect_equal(class(pp[[1L]]), "data.frame")
})

test_that("lgcp popa model works from fitted quad", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, return.type = "ppp")
  testthat::expect_equal(class(pp), "ppp")
})

test_that("lgcp popa model can return data frame from sampled intensity", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, which.intensity = "sample", return.type = "data.frame")
  testthat::expect_equal(class(pp), "data.frame")
})

test_that("lgcp popa model works with multiple simulations", {
  dat_po <- flora$po$sp1
  dat_pa <- flora$pa
  dat_po <- rbind.data.frame(dat_po, flora$quad)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat_po)
  mod <- popa(pres ~ MNT + D.Main, sp1 ~ MNT, dat_po, dat_pa, model.type = "laplace", simple.basis = bfs)
  pp <- simulate(mod, nsim = 3, return.type = "ppp")
  testthat::expect_equal(class(pp[[1L]]), "ppp")
})
