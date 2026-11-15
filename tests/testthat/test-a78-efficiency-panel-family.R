## efficiency() on the models whose score is stored under a name the fallback
## chain did not know (gap A78): psfm()'s GTRE/TRE family, which carries the
## transient score as $U and the persistent one as $H, and esfm(), which
## writes $Efficiency with a capital E.
##
## Most of this is asserted against fits BUILT BY HAND. What the gap was about
## is which slot the chain reads and what it says when there is nothing to
## read, and a hand-built object states that contract without depending on an
## optimizer reaching any particular point. One real GTRE fit then checks that
## the contract holds on an object psfm() actually produces.

## A minimal "sfareg" carrying only the slots the chain looks at.
fake_fit <- function(..., model_name = "GTRE") {
  z <- list(model_name = model_name, out = matrix(0, 1, 3,
            dimnames = list("beta", c("par", "st_err", "t-val"))))
  z <- utils::modifyList(z, list(...))
  class(z) <- "sfareg"
  z
}

test_that("the chain reads $U, so a GTRE-shaped fit is not refused", {
  f <- fake_fit(U = c(0.8, 0.6, 0.9), H = c(0.7, 0.75))
  expect_equal(suppressMessages(efficiency(f)), c(0.8, 0.6, 0.9))
  expect_equal(efficiency(f, component = "persistent"), c(0.7, 0.75))
})

test_that("the chain reads $Efficiency, so esfm() is not refused", {
  f <- fake_fit(Efficiency = c(0.5, 0.55), model_name = "ESFM")
  expect_equal(efficiency(f), c(0.5, 0.55))
})

test_that("the earlier slots still win, so nothing already working moves", {
  ## exp_u_hat before efficiency before Efficiency before U: a fit that has
  ## more than one must keep answering from the first.
  f <- fake_fit(exp_u_hat = c(0.1, 0.2), efficiency = c(0.3, 0.4),
                Efficiency = c(0.5, 0.6), U = c(0.7, 0.8))
  expect_equal(efficiency(f), c(0.1, 0.2))
  expect_equal(efficiency(fake_fit(efficiency = c(0.3, 0.4), U = c(0.7, 0.8))),
               c(0.3, 0.4))
})

test_that("a fit with no persistent component refuses instead of substituting", {
  ## The whole point: falling back to the transient score here would be gap
  ## A9's component mix-up arriving through a different door.
  f <- fake_fit(U = c(0.8, 0.6), model_name = "TRE")
  expect_error(efficiency(f, component = "persistent"),
               "stores no persistent efficiency")
  expect_error(efficiency(fake_fit(exp_u_hat = 0.5, model_name = "NHN"),
                          component = "persistent"),
               "stores no persistent efficiency")
})

test_that("the refusal names a type the fit can actually answer", {
  ## The old message advised type = "jlms" on a fit that stores no E[u | e],
  ## so following it failed again.
  f <- fake_fit(U = c(0.8, 0.6), H = 0.7)
  ## No message either: a call that is about to error should not also advise.
  expect_silent(try(efficiency(f, type = "jlms"), silent = TRUE))
  e <- tryCatch(efficiency(f, type = "jlms"), error = function(z) conditionMessage(z))
  expect_match(e, "can answer type = \"bc\"", fixed = TRUE)
  expect_false(grepl("use type = \"jlms\"", e, fixed = TRUE))

  ## And a fit with nothing at all says so rather than pointing anywhere.
  e2 <- tryCatch(efficiency(fake_fit(model_name = "NONE")),
                 error = function(z) conditionMessage(z))
  expect_match(e2, "stores no efficiency predictor at all")
})

test_that("jlms is never derived from the stored bc score", {
  ## -log(E[exp(-u)|e]) is not E[u|e]; they differ by Jensen. Reporting one
  ## under the other's name is the defect gap A68 recorded, so the fit must
  ## refuse rather than convert.
  f <- fake_fit(U = c(0.8, 0.6), H = 0.7)
  expect_error(efficiency(f, type = "jlms"), "neither E\\[u \\| e\\]")
})

test_that("a per-firm score refuses the level scale rather than recycling", {
  f <- fake_fit(U = rep(0.8, 6), H = c(0.7, 0.75))
  expect_error(efficiency(f, component = "persistent", logDepVar = FALSE),
               "one value per firm")
})

test_that("the two-efficiency message fires only when the caller did not choose", {
  f <- fake_fit(U = c(0.8, 0.6), H = 0.7)
  expect_message(efficiency(f), "reports two efficiencies")
  expect_silent(efficiency(f, component = "transient"))
  expect_silent(efficiency(f, component = "persistent"))
  ## A fit with one efficiency says nothing.
  expect_silent(efficiency(fake_fit(U = c(0.8, 0.6), model_name = "TRE")))
})

test_that("on a real GTRE fit, efficiency() and predict() finally agree", {
  skip_on_cran()
  d <- data_gen_p(N = 25, t = 4, rand = 2, sig_u = 1, sig_v = .5, sig_r = .5,
                  sig_h = .5, cons = .5, beta1 = .5, beta2 = .5)
  ## suppressWarnings deliberately: on a panel this small sigr often collapses
  ## to its boundary, which psfm() reports at length and which is documented
  ## behaviour rather than a failure. Nothing asserted below depends on where
  ## the optimizer lands -- only on the length, the range and the two generics
  ## agreeing -- which is what lets the fixture stay this cheap.
  f <- suppressWarnings(
    psfm(y_gtre ~ x1 + x2, data = d, model_name = "GTRE", individual = "name"))

  te <- suppressMessages(efficiency(f))
  expect_length(te, nrow(d))
  expect_true(all(te > 0 & te <= 1))
  ## The two generics disagreed about whether this model has an efficiency
  ## measure at all; predict() always carried the $U rung and efficiency() did
  ## not. They are now the same number.
  expect_equal(unname(te), unname(as.numeric(predict(f, type = "efficiency"))))

  hp <- efficiency(f, component = "persistent")
  expect_length(hp, length(unique(d$name)))
  expect_true(all(hp > 0 & hp <= 1))
  expect_equal(unname(hp), unname(as.numeric(f$H)))
})
