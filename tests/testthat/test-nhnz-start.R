## NHN_Z/NE_Z used to stop short of the maximum on an uncentred z, below the
## NHN model they nest (JSS referee item 17; notes/code_history/sfm.md).

test_that("NHN_Z and NE_Z reach the optimum of the same model fitted on the vhet path", {
  skip_on_cran()
  data(panel89, package = "sfa", envir = environment())
  d <- panel89[seq(1, nrow(panel89), by = 4), ]
  fz <- y ~ q1 + q2 + q3 + q4 + q5 + w1 + w2 + w3 + w4 | z
  for (mn in c("NHN_Z", "NE_Z")) {
    f <- suppressWarnings(sfm(fz, model_name = mn, inefdec = FALSE, data = d))
    h <- suppressWarnings(sfm(fz, model_name = mn, vhet = ~1, inefdec = FALSE, data = d))
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(h)), tolerance = 1e-4, info = mn)
    ## the returned layout is the homoskedastic-v one, not the vhet one
    expect_equal(rownames(f$out)[1], "sigma_v", info = mn)
  }
  n <- suppressWarnings(sfm(y ~ q1 + q2 + q3 + q4 + q5 + w1 + w2 + w3 + w4,
    model_name = "NHN", inefdec = FALSE, data = d
  ))
  f <- suppressWarnings(sfm(fz, model_name = "NHN_Z", inefdec = FALSE, data = d))
  expect_gte(as.numeric(logLik(f)), as.numeric(logLik(n)) - 1e-6)
})
