context("Presence-Only Model Function: po()")
# usethis::use_test() when particular R function file is open to create these

test_that("ipp po model creates scampr model", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (va) po model creates scampr model", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (laplace) po model creates scampr model", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "laplace", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp po model creates scampr model w. logical subset", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  row.sample <- c(sample(rownames(dat[dat$pres == 1, ]), 0.5 * sum(dat$pres)), rownames(dat[dat$pres == 0, ]))
  mod <- po(pres ~ elev, dat, model.type = "ipp", subset = rownames(dat) %in% row.sample)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp po model creates scampr model w. numbered subset", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  mod <- po(pres ~ elev, dat, model.type = "ipp", subset = 1:1000)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (va) model creates scampr model w. logical subset", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  row.sample <- c(sample(rownames(dat[dat$pres == 1, ]), 0.5 * sum(dat$pres)), rownames(dat[dat$pres == 0, ]))
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
  mod <- po(pres ~ elev, dat, model.type = "variational", simple.basis = bfs, subset = rownames(dat) %in% row.sample)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("ipp po model creates scampr model works with interaction terms in formula", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  dat$dist2water <- scale(dat$waterdist)
  mod <- po(pres ~ elev*dist2water, dat, model.type = "ipp")
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (va) po model works with interaction terms in formula", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  dat$dist2water <- scale(dat$waterdist)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
  mod <- po(pres ~ elev*dist2water, dat, model.type = "variational", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})

test_that("lgcp (laplace) po model works with interaction terms in formula", {
  dat <- scampr::gorillas
  dat$elev <- scale(dat$elevation)
  dat$dist2water <- scale(dat$waterdist)
  bfs <- simple_basis(nodes.on.long.edge = 5, data = dat)
  mod <- po(pres ~ elev*dist2water, dat, model.type = "laplace", simple.basis = bfs)
  testthat::expect_equal(class(mod), "scampr")
})
