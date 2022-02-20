library(modeldata)
data("Sacramento")

test_that("test step_knn_classdist basic", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_knn_classdist(latitude, longitude, class = "type", neighbors = 1)

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  names(result)

  expect_true(all(
    c(
      "Condo_neighbor1",
      "Multi_Family_neighbor1",
      "Residential_neighbor1"
    ) %in% names(result)
  ))
})
