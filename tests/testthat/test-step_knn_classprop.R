library(modeldata)
data("Sacramento")

test_that("test step_knn_classdist basic", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_knn_classprop(latitude, longitude, class = "type", neighbors = 3)

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  expect_true(all(
    c(
      "prop_Condo",
      "prop_Multi_Family",
      "prop_Residential"
    ) %in% names(result)
  ))
})
