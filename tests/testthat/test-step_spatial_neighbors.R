library(modeldata)
data("Sacramento")

test_that("test step_spatial_neighbors basic", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_spatial_neighbors(latitude, longitude, outcome = "type", neighbors = 3)

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  expect_true(all(
    c(
      "nn_type_1",
      "nn_type_2",
      "nn_type_3",
      "dist_type_1",
      "dist_type_2",
      "dist_type_3"
    ) %in% names(result)
  ))
})
