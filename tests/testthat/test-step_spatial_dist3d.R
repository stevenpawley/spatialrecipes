library(modeldata)
library(testthat)

data("Sacramento")

test_that("test step_geodist3d incorrect parameters", {
  # catch error if no reference locations are provided
  expect_error({
    rec_obj <- Sacramento %>%
      recipe(price ~ .) %>%
      step_geodist3d(latitude, longitude, sqft)
  })

  # catch error if reference locations have different lengths
  expect_error({
    rec_obj <- Sacramento %>%
      recipe(price ~ .) %>%
      step_spatial_dist3d(latitude, longitude, sqft,
        ref_lat = 40, ref_lon = -120,
        ref_height = c(500, 600)
      )
  })
})


test_that("test step_geodist3d single location", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_spatial_dist3d(latitude, longitude, sqft,
      ref_lat = 40, ref_lon = -120,
      ref_height = 500
    )

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  expect_length(names(result), ncol(Sacramento) + 1)
  expect_true("spatial_dist3d" %in% names(result))
  expect_true(inherits(result$spatial_dist3d, "numeric"))
  expect_equal(sum(is.na(result$spatial_dist3d)), 0)
})


test_that("test step_geodist3d multiple locations", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_spatial_dist3d(latitude, longitude, sqft,
      ref_lat = c(40, 50, 42),
      ref_lon = c(-120, -110, -105),
      ref_height = c(50, 55, 5)
    )

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  expect_length(names(result), ncol(Sacramento) + 3)

  expect_true(
    all(c("spatial_dist3d_1", "spatial_dist3d_2", "spatial_dist3d_3") %in% names(result))
  )

  expect_true(inherits(result$spatial_dist3d_1, "numeric"))
  expect_true(inherits(result$spatial_dist3d_2, "numeric"))
  expect_true(inherits(result$spatial_dist3d_3, "numeric"))

  expect_equal(sum(is.na(result$spatial_dist3d_1)), 0)
  expect_equal(sum(is.na(result$spatial_dist3d_2)), 0)
  expect_equal(sum(is.na(result$spatial_dist3d_3)), 0)
})


test_that("test step_geodist3d minimum distance to multiple locations", {
  rec_obj <- Sacramento %>%
    recipe(price ~ .) %>%
    step_spatial_dist3d(latitude, longitude, sqft,
      ref_lat = c(40, 50, 42),
      ref_lon = c(-120, -110, -105),
      ref_height = c(50, 55, 5), minimum = TRUE
    )

  prepped <- prep(rec_obj)
  result <- juice(prepped)

  expect_length(names(result), ncol(Sacramento) + 1)
  expect_true("spatial_dist3d" %in% names(result))
  expect_true(inherits(result$spatial_dist3d, "numeric"))
  expect_equal(sum(is.na(result$spatial_dist3d)), 0)
})
