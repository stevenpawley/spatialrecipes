library(sp)

test_that("test step_knn", {
  data("meuse")

  rec_obj <- meuse %>%
    recipe(lead ~ dist + elev + x + y) %>%
    step_knn(x, y, outcome = "lead", neighbors = 3, weight_func = "inv")

  prepped <- prep(rec_obj)
  result <- juice(prepped)
})
