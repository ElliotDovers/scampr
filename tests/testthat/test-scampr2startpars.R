context("INTERNAL: scampr model to starting pars list")

# Presence-Only Models:

test_that("ipp po model creates scampr model", {
  fake.mod <- list()
  fake.mod$par <- rep(0, 6)
  names(fake.mod$par) <- rep("fixed", 6)
  fake.mod$approx.type <- NA
  startpars <- scampr2startpars(fake.mod, "ipp")
})

test_that("lgcp (laplace) po model creates scampr model", {

})

test_that("lgcp (va) po model creates scampr model", {

})

# Presence/Absence Models:

test_that("ipp pa model creates scampr model", {

})

test_that("lgcp (va) pa model doesn't work", {
})

test_that("lgcp (laplace) pa model creates scampr model", {

  })

# Combined Data Models:

test_that("ipp popa model creates scampr model", {

})

test_that("lgcp (va) popa model doesn't work", {
})

test_that("lgcp (laplace) popa model creates a scampr model", {
})
