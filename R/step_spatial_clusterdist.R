#' Distance between locations in two dimensions
#'
#' `step_spatial_clusterdist` creates a *specification* of a recipe step that will
#' calculate the distance between points on a map to a reference location.
#'
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param ref_lon A character with the name of the column in the training data
#'   that represents that x-coordinate of the spatial domain.
#' @param ref_lat A character with the name of the column in the training data
#'   that represents that y-coordinate of the spatial domain.
#' @param num_comp Number of clusters to generate.
#' @param role or model term created by this step, what analysis role should be
#'   assigned?. By default, the function assumes that resulting distance will be
#'   used as a predictor in a model.
#' @param columns A character string of variable names that will be populated
#'   (eventually) by the `terms` argument.
#' @param centers The cluster geographic locations of the cluster centroids
#'   based on `ref_lon,` `ref_lat` and the clusters generated from the columns
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
step_spatial_clusterdist <-
  function(recipe,
           ...,
           ref_lon,
           ref_lat,
           num_comp = 5,
           role = "predictor",
           trained = FALSE,
           name = "clusterdist",
           columns = NULL,
           centers = NULL,
           skip = FALSE,
           id = recipes::rand_id("clusterdist")) {
    terms <- recipes::ellipse_check(...)

    if (length(name) != 1 || !is.character(name)) {
      rlang::abort("`name` should be a single character value.")
    }

    recipes::add_step(
      recipe,
      step_spatial_clusterdist_new(
        terms = terms,
        ref_lon = ref_lon,
        ref_lat = ref_lat,
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
           ref_lon,
           ref_lat,
           num_comp,
           role,
           trained,
           name,
           columns,
           centers,
           skip,
           id) {
    recipes::step(
      subclass = "clusterdist",
      terms = terms,
      ref_lon = ref_lon,
      ref_lat = ref_lat,
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
  col_names <- recipes::terms_select(terms = x$terms, info = info)
  km <- kmeans(training[col_names],
    centers = x$num_comp, algorithm = "Lloyd",
    iter.max = 1000
  )

  step_spatial_clusterdist_new(
    terms = x$terms,
    ref_lon = x$ref_lon,
    ref_lat = x$ref_lat,
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
  geo_dist_2d_calc <- function(df, a, b) {
    apply(df, 1, function(x, a, b, c) {
      sqrt((x[1] - a)^2 + (x[2] - b)^2)
    },
    a = a, b = b
    )
  }

  cols <- c(object$ref_lon, object$ref_lat)

  dist_vals <- geo_dist_2d_calc(
    df = new_data[, cols],
    a = object$centers[[object$ref_lon]],
    b = object$centers[[object$ref_lat]]
  )

  if (inherits(dist_vals, "numeric")) {
    dist_vals <- tibble(dist_vals)
    names(dist_vals) <- object$name
  } else if (inherits(dist_vals, "matrix")) {
    dist_vals <- as.data.frame(t(dist_vals))
    dist_vals <- setNames(
      dist_vals,
      paste0(object$name, seq_len(ncol(dist_vals)))
    )
  }

  bind_cols(new_data, dist_vals)
}

print.step_spatial_clusterdist <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat(
      "Geographical distances from",
      format(x$ref_lat, digits = 10),
      format(x$ref_lon, digits = 10),
      "\n"
    )
    invisible(x)
  }

#' @rdname step_spatial_clusterdist
#' @param x A `step_spatial_clusterdist` object.
#' @export
tidy.step_spatial_clusterdist <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble(
      terms = x$columns,
      ref_latitude = x$ref_lat,
      ref_longitude = x$ref_lon,
      name = x$name
    )
  } else {
    res <- tibble(
      terms = recipes::sel2char(x$terms),
      ref_latitude = x$ref_lat,
      ref_longitude = x$ref_lon,
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
