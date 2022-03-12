library(sp)

test_that("test step_knn", {
  data("meuse")

  rec_obj <- meuse %>%
    recipe(lead ~ dist + elev + x + y) %>%
    step_knn(x, y, outcome = "lead", neighbors = 3, weight_func = "gaussian")

  prepped <- prep(rec_obj)
  result <- juice(prepped)
})
