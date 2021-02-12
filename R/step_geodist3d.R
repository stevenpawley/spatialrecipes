#' Distance between two locations
#'
#' `step_geodist3d` creates a *specification* of a
#'  recipe step that will calculate the distance between
#'  points on a map to a reference location.
#'
#' @param lon,lat,z Selector functions to choose which variables are
#'  affected by the step. See selections() for more details.
#' @param ref_lon,ref_lat,ref_z Single numeric values for the location
#'  of the reference point.
#' @param role or model term created by this step, what analysis
#'  role should be assigned?. By default, the function assumes
#'  that resulting distance will be used as a predictor in a model.
#' @param log A logical: should the distance be transformed by
#'  the natural log function?
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#' @param name A single character value to use for the new
#'  predictor column. If a column exists with this name, an error is
#'  issued.
#' @return An updated version of `recipe` with the new step added
#'  to the sequence of existing steps (if any). For the `tidy`
#'  method, a tibble with columns echoing the values of `lat`,
#'  `lon`, `z`, `ref_lat`, `ref_lon`, `ref_z`, `name`, and `id`.
#' @keywords datagen
#' @concept preprocessing
#' @export
step_geodist3d <- function(recipe,
                           lat = NULL,
                           lon = NULL,
                           z = NULL,
                           role = "predictor",
                           trained = FALSE,
                           ref_lat = NULL,
                           ref_lon = NULL,
                           ref_z = NULL,
                           log = FALSE,
                           name = "geo_dist3d",
                           columns = NULL,
                           skip = FALSE,
                           id = recipes::rand_id("geodist3d")) {
  if (length(ref_lon) != 1 || !is.numeric(ref_lon))
    rlang::abort("`ref_lon` should be a single numeric value.")

  if (length(ref_lat) != 1 || !is.numeric(ref_lat))
    rlang::abort("`ref_lat` should be a single numeric value.")

  if (length(ref_z) != 1 || !is.numeric(ref_z))
    rlang::abort("`ref_z` should be a single numeric value.")

  if (length(log) != 1 || !is.logical(log))
    rlang::abort("`log` should be a single logical value.")

  if (length(name) != 1 || !is.character(name))
    rlang::abort("`name` should be a single character value.")

  recipes::add_step(
    recipe,
    step_geodist3d_new(
      lon = rlang::enquos(lon),
      lat = rlang::enquos(lat),
      z = rlang::enquos(z),
      role = role,
      trained = trained,
      ref_lon = ref_lon,
      ref_lat = ref_lat,
      ref_z = ref_z,
      log = log,
      name = name,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

step_geodist3d_new <-
  function(lon, lat, z, role, trained, ref_lon, ref_lat,
           ref_z, log, name, columns, skip, id) {
    recipes::step(
      subclass = "geodist3d",
      lon = lon,
      lat = lat,
      z = z,
      role = role,
      trained = trained,
      ref_lon = ref_lon,
      ref_lat = ref_lat,
      ref_z = ref_z,
      log = log,
      name = name,
      columns = columns,
      skip = skip,
      id = id
    )
  }


#' @export
prep.step_geodist3d <- function(x, training, info = NULL, ...) {
  lat_name <- recipes::terms_select(terms = x$lat, info = info)
  lon_name <- recipes::terms_select(terms = x$lon, info = info)
  z_name <- recipes::terms_select(terms = x$z, info = info)

  step_geodist3d_new(
    lon = x$lon,
    lat = x$lat,
    z = x$z,
    role = x$role,
    trained = TRUE,
    ref_lon = x$ref_lon,
    ref_lat = x$ref_lat,
    ref_z = x$ref_z,
    log = x$log,
    name = x$name,
    columns = c(lat_name, lon_name, z_name),
    skip = x$skip,
    id = x$id
  )
}

geo_dist_3d_calc <- function(x, a, b, c)
  apply(x, 1, function(x, a, b, c)  {
    sqrt((x[1] - a) ^ 2 + (x[2] - b) ^ 2 + (x[3] - c) ^ 2)
  },
  a = a, b = b, c = c)

#' @export
bake.step_geodist3d <- function(object, new_data, ...) {
  dist_vals <-
    geo_dist_3d_calc(new_data[, object$columns], object$ref_lat, object$ref_lon,
                     object$ref_z)

  if (object$log) {
    new_data[, object$name] <- log(dist_vals)
  } else {
    new_data[, object$name] <- dist_vals
  }
  new_data
}

print.step_geodist3d <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Geographical distances from",
        format(x$ref_lat, digits = 10),
        format(x$ref_lon, digits = 10),
        format(x$ref_z, digits = 10),
        "\n")
    invisible(x)
  }



#' @rdname step_geodist3d
#' @param x A `step_geodist3d` object.
#' @export
tidy.step_geodist3d <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble(
      latitude = x$columns[1],
      longitude = x$columns[2],
      z = x$columns[3],
      ref_latitude = x$ref_lat,
      ref_longitude = x$ref_lon,
      ref_z = x$ref_z,
      name = x$name
    )
  } else {
    res <- tibble(
      latitude = recipes::sel2char(x$lat),
      longitude = recipes::sel2char(x$lon),
      z = recipes::sel2char(x$z),
      ref_latitude = x$ref_lat,
      ref_longitude = x$ref_lon,
      ref_z = x$z,
      name = x$name
    )
  }
  res$id <- x$id
  res
}
