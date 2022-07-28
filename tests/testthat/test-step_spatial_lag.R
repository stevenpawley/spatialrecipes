library(sp)

test_that("test step_spatial_lag", {
  data("meuse")

  rec_obj <- meuse %>%
    recipe(lead ~ dist + elev + x + y) %>%
    step_spatial_lag(x, y, outcome = "lead", neighbors = 3, weight_func = "gaussian")

  prepped <- prep(rec_obj)
  result <- juice(prepped)
  testthat::expect_named(result, c(
    "dist", "elev", "x", "y", "lead",
    "lead_lag_3_gaussian"
  ))
})
