library("testthat");


source("tipping.R");


context("auxiliary stuff");
test_that("commaSep", {
  expect_identical(commaSep(c(1.2, 3.4), "%3.1f"), "1.2, 3.4");
});

context("cubic differential equation basics");
test_that("cubic analysis", {
  e <- tippingCubicExtrema(1, 1, 0);
  expect_identical(e[1], -e[2]);
  expect_identical(length(tippingCubicStableFixedPoints(1, 1, -tippingCubicCcrit - 1e-7)), 1L);
  expect_identical(length(tippingCubicStableFixedPoints(1, 1, -tippingCubicCcrit + 1e-7)), 2L);
  expect_identical(length(tippingCubicStableFixedPoints(1, 1, tippingCubicCcrit + 1e-7)), 1L);
  expect_identical(length(tippingCubicStableFixedPoints(1, 1, tippingCubicCcrit - 1e-7)), 2L);
});

context("coupled tipping elements");
test_that("coupled tipping fixed points", {
  sse <- function(fp, ts)
  {
    d <- fp - ts;
    return(sum(d * d));
  }
  
  set.seed(33);
  n <- 3L;
  p <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -0.2, 0.2);
  initialState <- makeCoupledTippingInitialState(rnorm(n), p);
  fp <- coupledTippingStableFixedPoint(p, initialState);
  ts1 <- coupledTippingTimeSeries(p, 20L, 0.1, initialState);
  ts2 <- coupledTippingTimeSeries(p, 50L, 0.1, initialState);
  ts3 <- coupledTippingTimeSeries(p, 100L, 0.1, initialState);
  sse1 <- sse(fp, ts1[nrow(ts1), 2:ncol(ts1)]);
  sse2 <- sse(fp, ts2[nrow(ts2), 2:ncol(ts2)]);
  sse3 <- sse(fp, ts3[nrow(ts3), 2:ncol(ts3)]);
  ## message(sprintf("sse1 = %f, sse2 = %f, sse3 = %f", sse1, sse2, sse3));
  expect_true(sse1 >= sse2);
  expect_true(sse2 >= sse3);
});


