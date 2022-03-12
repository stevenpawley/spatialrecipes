#' Distance between locations in three dimensions
#'
#' `step_spatial_dist3d` creates a *specification* of a recipe step that will
#' calculate the distance between points on a map to a reference location.
#'
#' @param lon,lat,height Selector functions to choose which variables are
#'   affected by the step. See selections() for more details.
#' @param ref_lon,ref_lat,ref_height Numeric values for the location of the
#'   reference points. New features representing the euclidean distance to each
#'   `ref_lon`, `ref_lat` and `ref_height` location will be created, unless
#'   `minimum = TRUE` when the minimum distance to the combined features will be
#'   created as a single new feature.
#' @param minimum A logical for whether only the minimum distance to the
#'   reference locations should be calculated.
#' @param role or model term created by this step, what analysis role should be
#'   assigned?. By default, the function assumes that resulting distance will be
#'   used as a predictor in a model.
#' @param log A logical: should the distance be transformed by the natural log
#'   function?
#' @param columns A character string of variable names that will be populated
#'   (eventually) by the `terms` argument.
#' @param name A single character value to use for the new predictor column. If
#'   a column exists with this name, an error is issued. If multiple reference
#'   locations are used, then the name is appended with `1, 2, 3 ...`
#'   representing the index of each reference location.
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any). For the `tidy` method, a tibble with
#'   columns echoing the values of `lat`, `lon`, `height`, `ref_lat`, `ref_lon`,
#'   `ref_height`, `name`, and `id`.
#' @keywords datagen
#' @concept preprocessing
#' @export
step_spatial_dist3d <- function(recipe,
                           lat = NULL,
                           lon = NULL,
                           height = NULL,
                           ref_lat = NULL,
                           ref_lon = NULL,
                           ref_height = NULL,
                           minimum = FALSE,
                           log = FALSE,
                           role = "predictor",
                           trained = FALSE,
                           name = "geo_dist3d",
                           columns = NULL,
                           skip = FALSE,
                           id = recipes::rand_id("geodist3d")) {
  vect_lengths <- c(length(ref_lon), length(ref_lat), length(ref_height))

  if (!all(vect_lengths == vect_lengths[1])) {
    rlang::abort("`ref_lon`, `ref_lat` and `ref_height` should be the same length.")
  }

  if (all(vect_lengths == 0)) {
    rlang::abort("at least a single reference location must be provided.")
  }

  if (length(log) != 1 || !is.logical(log)) {
    rlang::abort("`log` should be a single logical value.")
  }

  if (length(name) != 1 || !is.character(name)) {
    rlang::abort("`name` should be a single character value.")
  }

  recipes::recipes_pkg_check("nabor")

  recipes::add_step(
    recipe,
    step_spatial_dist3d_new(
      lon = rlang::enquos(lon),
      lat = rlang::enquos(lat),
      height = rlang::enquos(height),
      role = role,
      trained = trained,
      ref_lon = ref_lon,
      ref_lat = ref_lat,
      ref_height = ref_height,
      minimum = minimum,
      log = log,
      name = name,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}


step_spatial_dist3d_new <-
  function(lon, lat, height, role, trained, ref_lon, ref_lat, ref_height,
           minimum, log, name, columns, skip, id) {
    recipes::step(
      subclass = "geodist3d",
      lon = lon,
      lat = lat,
      height = height,
      role = role,
      trained = trained,
      ref_lon = ref_lon,
      ref_lat = ref_lat,
      ref_height = ref_height,
      minimum = minimum,
      log = log,
      name = name,
      columns = columns,
      skip = skip,
      id = id
    )
  }


#' @export
prep.step_spatial_dist3d <- function(x, training, info = NULL, ...) {
  lat_name <- recipes::terms_select(terms = x$lat, info = info)
  lon_name <- recipes::terms_select(terms = x$lon, info = info)
  height_name <- recipes::terms_select(terms = x$height, info = info)

  step_spatial_dist3d_new(
    lon = x$lon,
    lat = x$lat,
    height = x$height,
    role = x$role,
    trained = TRUE,
    ref_lon = x$ref_lon,
    ref_lat = x$ref_lat,
    ref_height = x$ref_height,
    minimum = x$minimum,
    log = x$log,
    name = x$name,
    columns = c(lat_name, lon_name, height_name),
    skip = x$skip,
    id = x$id
  )
}


#' Euclidean distance calculation
#'
#' @param x A numeric vector of locations
#' @param a The reference longitude.
#' @param b The reference latitude.
#' @param c The reference height.
#'
#' @return Will return a numeric vector for a single reference location, or a
#'   matrix for multiple locations with (n_locations, n_distances).
#' @keywords internal
geo_dist_3d_calc <- function(x, a, b, c) {
  apply(x, 1, function(x, a, b, c) {
    sqrt((x[1] - a)^2 + (x[2] - b)^2 + (x[3] - c)^2)
  },
  a = a, b = b, c = c
  )
}


#' @export
bake.step_spatial_dist3d <- function(object, new_data, ...) {
  if (isFALSE(object$minimum)) {
    dist_vals <- geo_dist_3d_calc(
      x = new_data[, object$columns],
      a = object$ref_lat,
      b = object$ref_lon,
      c = object$ref_height
    )
  } else {
    refs <- tibble(
      y = object$ref_lat,
      x = object$ref_lon,
      z = object$ref_height
    )
    nn <- nabor::knn(
      data = as.matrix(refs),
      query = as.matrix(new_data[, object$columns]),
      k = 1
    )
    dist_vals <- as.numeric(nn$nn.dists)
  }

  if (object$log) {
    dist_vals <- log(dist_vals)
  }

  if (inherits(dist_vals, "numeric")) {
    new_data[[object$name]] <- dist_vals
  }

  if (inherits(dist_vals, "matrix")) {
    dist_vals <- as.data.frame(t(dist_vals))
    feature_names <- paste(object$name, seq_len(ncol(dist_vals)), sep = "_")
    names(dist_vals) <- feature_names
    new_data <- bind_cols(new_data, dist_vals)
  }

  new_data
}

#' @export
print.step_spatial_dist3d <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat(
      "Geographical distances from",
      format(x$ref_lat, digits = 10),
      format(x$ref_lon, digits = 10),
      format(x$ref_height, digits = 10),
      "\n"
    )
    invisible(x)
  }



#' @rdname step_spatial_dist3d
#' @param x A `step_spatial_dist3d` object.
#' @export
tidy.step_spatial_dist3d <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble(
      latitude = x$columns[1],
      longitude = x$columns[2],
      height = x$columns[3],
      ref_latitude = x$ref_lat,
      ref_longitude = x$ref_lon,
      ref_height = x$ref_height,
      name = x$name
    )
  } else {
    res <- tibble(
      latitude = recipes::sel2char(x$lat),
      longitude = recipes::sel2char(x$lon),
      height = recipes::sel2char(x$height),
      ref_latitude = x$ref_lat,
      ref_longitude = x$ref_lon,
      ref_height = x$height,
      name = x$name
    )
  }
  res$id <- x$id
  res
}
