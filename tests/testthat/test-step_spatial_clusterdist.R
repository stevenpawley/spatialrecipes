library(modeldata)
data("Sacramento")

test_that("test step_clusterdist", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_spatial_clusterdist(all_numeric_predictors(),
      ref_lon = "longitude",
      ref_lat = "latitude", num_comp = 5
    )

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  expect_named(
    object = result,
    expected = c(names(Sacramento), paste0("spatial_clusterdist", 1:5)),
    ignore.order = TRUE
  )
})
