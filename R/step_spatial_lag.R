gaussian_kernel <- function(X, std = 2) {
  exp(-X^2 / std^2)
}

lag_train <- function(formula, x, y, k, weight_func, dist_power = 2) {
  target_variable <-
    rlang::f_lhs(formula) %>%
    as.character()
  term_variables <- attr(terms(formula), "term.labels")

  # split data
  train_data <- x[term_variables]
  query_data <- y[term_variables]

  # for training use k[2:k+1] i.e. use only neighboring points but not the
  # point itself
  if (identical(train_data, query_data)) {
    neighbors <- nabor::knn(data = train_data, query = query_data, k = k + 1)
    neighbors$nn.idx <- neighbors$nn.idx[, 2:ncol(neighbors$nn.idx)]
    neighbors$nn.dists <- neighbors$nn.dists[, 2:ncol(neighbors$nn.dists)]

  } else {
    neighbors <- nabor::knn(data = train_data, query = query_data, k = k)
  }

  # get ids and distances to neighbors
  neighbor_ids <- neighbors$nn.idx
  D <- neighbors$nn.dists

  # create initial row-normalized weights from distances
  W <- D / rowSums(D)

  # get values of neighbors
  neighbor_vals <- x[as.integer(neighbor_ids), ][[target_variable]]
  neighbor_vals <- matrix(neighbor_vals, ncol = k)

  # calculate weights
  W <- switch(
    weight_func,
    inv = 1 / W^dist_power,
    rectangular = matrix(1, nrow = nrow(neighbor_vals), ncol = k),
    gaussian = gaussian_kernel(W, dist_power)
  )

  # calculate weighted mean/mode of neighbors
  if (inherits(x[[target_variable]], "numeric")) {
    # W <- W / rowSums(W)
    # fitted <- rowSums(neighbor_vals * W)

    fitted <- sapply(seq_len(nrow(W)), function(i)
      weighted.mean(x = neighbor_vals[i,], w = W[i,]))

  } else {
    fitted <- sapply(seq_len(nrow(W)), function(i) {
      collapse::fmode(x = neighbor_vals[i,], w = W[i,])
    })
  }

  new_feature_name <- paste(target_variable, "lag", k, weight_func, sep = "_")

  fitted <- tibble(fitted) %>%
    rlang::set_names(new_feature_name)

  fitted
}


#' Spatial lag step
#'
#' `step_spatial_lag` creates a *specification* of a recipe step that will add a
#' new 'lag' feature to a dataset based on the inverse distance-weighted mean of
#' surrounding observations.
#'
#' @param recipe A recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param outcome Selector function to choose which variable will be used to
#'   create a new feature based on the inverse distance-weighted mean of
#'   surrounding observations.
#' @param role role or model term created by this step, what analysis role
#'   should be assigned?. By default, the function assumes that resulting
#'   distance will be used as a predictor in a model.
#' @param trained A logical that will be updated once the step has been trained.
#' @param neighbors The number of closest neighbours to use in the distance
#'   weighting.
#' @param weight_func A single character for the kernel function used to weight
#'   the distances between samples.
#' @param dist_power Power function for "inv".
#' @param data Used internally to store the training data.
#' @param skip A logical to skip training.
#' @param id An identifier for the step. If omitted then this is generated
#'   automatically.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any).
#' @keywords datagen
#' @concept preprocessing
#' @export
#' @details
#' No details yet!
#'
#' @examples
#' library(modeldata)
#' library(recipes)
#' library(spatialrecipes)
#' data(ames)
#'
#' rec_obj <- ames %>%
#' recipe(Sale_Price ~ Latitude + Longitude) %>%
#' step_spatial_lag(Latitude, Longitude, outcome = "Sale_Price", neighbors = 3,
#'                  weight_func = "inv", dist_power = 2)
#'
#' prepped <- prep(rec_obj)
#' juice(prepped)
step_spatial_lag <- function(
  recipe, ...,
  outcome = NULL,
  role = "predictor",
  trained = FALSE,
  neighbors = NA,
  weight_func = NULL,
  dist_power = NA,
  data = NULL,
  columns = NULL,
  skip = FALSE,
  id = recipes::rand_id("spatial_lag")) {

  recipes::recipes_pkg_check("nabor")

  terms <- recipes::ellipse_check(...)

  if (neighbors <= 0)
    rlang::abort("`neighbors` should be greater than 0.")

  if (!weight_func %in% c("inv", "gaussian", "rectangular"))
    rlang::abort("`weight_func` should be either 'inv', 'gaussian', or 'rectangular'")

  recipes::add_step(
    recipe,
    step_spatial_lag_new(
      terms = terms,
      outcome = rlang::enquos(outcome),
      trained = trained,
      role = role,
      neighbors = neighbors,
      weight_func = weight_func,
      dist_power = dist_power,
      data = data,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

# wrapper around 'step' function that sets the class of new step objects
step_spatial_lag_new <- function(terms, role, trained, outcome, neighbors,
                                 weight_func, dist_power, data, columns, skip,
                                 id) {
  recipes::step(
    subclass = "spatial_lag",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    neighbors = neighbors,
    weight_func = weight_func,
    dist_power = dist_power,
    data = data,
    columns = columns,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_spatial_lag <- function(x, training, info = NULL, ...) {
  # First translate the terms argument into column name
  col_names <- recipes::terms_select(terms = x$terms, info = info)
  outcome_name <- recipes::terms_select(x$outcome, info = info)

  # Use the constructor function to return the updated object
  # Note that `trained` is set to TRUE
  step_spatial_lag_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    outcome = outcome_name,
    neighbors = x$neighbors,
    weight_func = x$weight_func,
    dist_power = x$dist_power,
    data = training,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_spatial_lag <- function(object, new_data, ...) {
  f <- as.formula(
    paste(object$outcome, paste(object$columns, collapse = " + "), sep = " ~ ")
  )

  lags <-
    lag_train(
      formula = f,
      x = object$data,
      y = new_data,
      k = object$neighbors,
      weight_func = object$weight_func,
      dist_power = object$dist_power
    )

  new_data <- dplyr::bind_cols(new_data, lags)
  new_data
}

#' @export
print.step_spatial_lag <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Spatial lags")

    if (recipes::is_trained(x)) {
      cat(paste0(" (", x$neighbors, " neighbors)"))
    }
    cat("\n")

    invisible(x)
  }

#' @rdname step_spatial_lag
#' @param x A `step_spatial_lag` object.
#' @export
tidy.step_spatial_lag <- function(x, ...) {
  term_names <- sel2char(x$terms)

  res <- tibble(
    id = rep(x$id, times = length(term_names)),
    terms = term_names,
    neighbors = rep(x$neighbors, times = length(term_names)),
    weight_func = rep(x$weight_func, times = length(term_names)),
    dist_power = rep(x$dist_power, times = length(term_names))
  )

  res
}

#' @rdname tunable.step
#' @export
tunable.step_spatial_lag <- function(x, ...) {
  tibble::tibble(
    name = c("neighbors", "weight_func", "dist_power"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1L, 10L)),
      list(pkg = "dials", fun = "weight_func", values = c("rectangular", "inv", "gaussian")),
      list(pkg = "dials", fun = "dist_power", range = c(1, 2))
    ),
    source = "recipe",
    component = "step_spatial_lag",
    component_id = x$id
  )
}

#' @rdname required_pkgs.step
#' @export
required_pkgs.step_spatial_lag <- function(x, ...) {
  c("nabor")
}
