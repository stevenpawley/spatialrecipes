#' Distance between locations in two dimensions
#'
#' `step_spatial_clusterdist` creates a *specification* of a recipe step that will
#' calculate the distance between points on a map to a reference location.
#'
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param ref A character with the name(s) of the columns in the training data
#'   that represents the spatial domain. For example,
#'   `ref = c("lat", "lon", "depth")` would calculate the centroids of clusters
#'   using the terms specified in the selector function, and then would measure
#'   the euclidean distances to these cluster centroids based on the `ref`
#'   variables. The variables specified in `ref` also need to be included
#'   as terms in the selector function.
#' @param num_comp Number of clusters to generate.
#' @param role or model term created by this step, what analysis role should be
#'   assigned?. By default, the function assumes that resulting distance will be
#'   used as a predictor in a model.
#' @param columns A character string of variable names that will be populated
#'   (eventually) by the `terms` argument.
#' @param centers The cluster geographic locations of the cluster centroids
#'   based on `ref,` and the clusters generated from the columns
#'   specified by `terms`.
#' @param name A single character value to use for the new predictor column. If
#'   a column exists with this name, an error is issued. If multiple clusters
#'   are used, then the name is appended with `1, 2, 3 ...` representing the
#'   index of each reference location.
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any).
#' @keywords datagen
#' @concept preprocessing
#' @export
#' @examples
#' library(tidyverse)
#' library(tidymodels)
#'
#' rec_obj <- Sacramento %>%
#' recipe(price ~ .) %>%
#'   step_spatial_clusterdist(
#'     price, latitude, longitude,
#'     ref = c("longitude", "latitude"),
#'     num_comp = 5
#'  )
#'
#'  prepped <- prep(rec_obj
#'  bake(prepped, new_data = NULL)
step_spatial_clusterdist <-
  function(recipe,
           ...,
           ref,
           num_comp = 5,
           role = "predictor",
           trained = FALSE,
           name = "spatial_clusterdist",
           columns = NULL,
           centers = NULL,
           skip = FALSE,
           id = recipes::rand_id("spatial_clusterdist")) {
    terms <- recipes::ellipse_check(...)

    if (length(name) != 1 || !is.character(name)) {
      rlang::abort("`name` should be a single character value.")
    }

    recipes::add_step(
      recipe,
      step_spatial_clusterdist_new(
        terms = terms,
        ref = ref,
        num_comp = num_comp,
        role = role,
        trained = trained,
        name = name,
        columns = columns,
        centers = centers,
        skip = skip,
        id = id
      )
    )
  }


step_spatial_clusterdist_new <-
  function(terms,
           ref,
           num_comp,
           role,
           trained,
           name,
           columns,
           centers,
           skip,
           id) {
    recipes::step(
      subclass = "spatial_clusterdist",
      terms = terms,
      ref = ref,
      num_comp = num_comp,
      role = role,
      trained = trained,
      name = name,
      columns = columns,
      centers = centers,
      skip = skip,
      id = id
    )
  }


#' @export
prep.step_spatial_clusterdist <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)

  km <- kmeans(
    training[col_names],
    centers = x$num_comp,
    algorithm = "Lloyd",
    iter.max = 1000
  )

  step_spatial_clusterdist_new(
    terms = x$terms,
    ref = x$ref,
    num_comp = x$num_comp,
    role = x$role,
    trained = TRUE,
    name = x$name,
    columns = col_names,
    centers = as_tibble(km$centers),
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_spatial_clusterdist <- function(object, new_data, ...) {
  ref_colnames <- object$ref
  ref_locations <- object$centers[ref_colnames]

  dist_vals <-
    apply(ref_locations, 1, function(ref, data) {
      euc_dist(x1 = data, x2 = ref)
    }, data = new_data[, ref_colnames])

  dist_vals <- as.data.frame(dist_vals)

  dist_vals <- setNames(
    dist_vals,
    paste0(object$name, seq_len(ncol(dist_vals)))
  )
  dist_vals <- as_tibble(dist_vals)

  bind_cols(new_data, dist_vals)
}

print.step_spatial_clusterdist <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat(
      "Geographical distances from",
      format(x$ref, digits = 10),
      "\n"
    )
    invisible(x)
  }

#' @rdname step_spatial_clusterdist
#' @param x A `step_spatial_clusterdist` object.
#' @export
tidy.step_spatial_clusterdist <- function(x, ...) {
  if (is_trained(x)) {
    centres <- setNames(as.data.frame(t(x$centers)), paste0(x$name, seq_len(x$num_comp)))

    res <- tibble(
      terms = x$columns,
      ref = paste(x$ref, collapse = ", "),
      name = x$name,
    )
    res <- bind_cols(res, centres)

  } else {
    res <- tibble(
      terms = recipes::sel2char(x$terms),
      ref = x$ref,
      name = x$name
    )
  }
  res$id <- x$id
  res
}

#' @export
tunable.step_spatial_clusterdist <- function(x, ...) {
  tibble::tibble(
    name = c("num_comp"),
    call_info = list(
      list(pkg = "dials", fun = "num_comp", range = c(1L, 10L))
    ),
    source = "recipe",
    component = "step_spatial_clusterdist",
    component_id = x$id
  )
}
