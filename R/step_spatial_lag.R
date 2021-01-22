gaussian_kernel <- function(W, k) {
  alpha <- 1 / (2 * (k + 1))
  qua <- abs(qnorm(alpha))
  W <- W * qua
  W <- dnorm(x = W, sd = 1)
  W
}

optimal_kernel <- function(W, k, d) {
  optKernel <- function(k, d = 1)
    1 / k * (1 + d/2 - d / (2 * k^(2 / d)) * ((1:k)^(1 + 2 / d) - (0:(k - 1))^(1 + 2 / d)))
  p <- dim(W)[1]
  W <- rep(optKernel(k, d = d), each = p)
  W <- matrix(W, nrow = p, ncol = k)
  W
}

lag_train <- function(formula, x, y, k, weight_func, dist = FALSE) {
  # get response and term variables
  target_variable <- formula %>%
    rlang::f_lhs() %>%
    as.character()
  term_variables <- attr(terms(formula), "term.labels")

  # split data
  response_data <- x[[target_variable]]
  train_data <- x[term_variables]
  query_data <- y[term_variables]

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

  # create initial row-standardized weights from distances
  maxdist <- D[, k]
  maxdist[maxdist < 1.0e-6] <- 1.0e-6
  W <- D / maxdist
  W <- pmin(W, 1 - (1e-6))
  W <- pmax(W, 1e-6)
  d <- length(term_variables)

  # get values of neighbors
  neighbor_vals <- x[as.numeric(neighbor_ids), ][[target_variable]]
  neighbor_vals <- matrix(neighbor_vals, ncol = k)

  # calculate weights (functions from kknn)
  if (weight_func == "rank") W <- (k + 1) - t(apply(D, 1, rank))
  if (weight_func == "inv") W <- 1/W
  if (weight_func == "triangular") W <- 1 - W
  if (weight_func == "rectangular") W <- matrix(1, nrow = nrow(neighbor_vals), ncol = k)
  if (weight_func == "epanechnikov") W <- 0.75*(1 - W^2)
  if (weight_func == "biweight") W <- dbeta((W + 1) / 2, 3, 3)
  if (weight_func == "triweight") W <- dbeta((W + 1) / 2, 4, 4)
  if (weight_func == "cos") W <- cos(W * pi / 2)
  if (weight_func == "triweights") W <- 1
  if (weight_func == "gaussian") W <- gaussian_kernel(W, k)
  if (weight_func == "optimal") W <- optimal_kernel(W, k, d)

  # calculate weighted mean of neighbors
  fitted <- rowSums(W*neighbor_vals) / pmax(rowSums(W))

  new_feature_name <- paste(target_variable, "lag", k, weight_func, sep = "_")
  fitted <- tibble(fitted) %>%
    rlang::set_names(new_feature_name)

  if (dist) {
    W <- matrix(1, nrow = nrow(D), ncol = k)
    new_feature_name <- paste(target_variable, "dist", k, sep = "_")
    fitted[[new_feature_name]] <- rowSums(W*D) / pmax(rowSums(W))
  }

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
#'   the distances between samples. Valid choices are: "rectangular",
#'   "triangular", "epanechnikov", "biweight", "tri-weight", "cos", "inv",
#'   "gaussian", "rank", and "optimal".
#' @param dist Whether to also return the weighted mean of the distances to the
#'   neighbours.
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
#' No examples yet!
step_spatial_lag <- function(
  recipe, ...,
  outcome = NULL,
  role = "predictor",
  trained = FALSE,
  neighbors = NA,
  weight_func = NULL,
  dist = FALSE,
  data = NULL,
  columns = NULL,
  skip = FALSE,
  id = recipes::rand_id("spatial_lag")) {

  recipes::recipes_pkg_check("nabor")

  terms <- recipes::ellipse_check(...)

  if (neighbors <= 0)
    rlang::abort("`neighbors` should be greater than 0.")

  recipes::add_step(
    recipe,
    step_spatial_lag_new(
      terms = terms,
      outcome = rlang::enquos(outcome),
      trained = trained,
      role = role,
      neighbors = neighbors,
      weight_func = weight_func,
      dist = dist,
      data = data,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

# wrapper around 'step' function that sets the class of new step objects
step_spatial_lag_new <- function(terms, role, trained, outcome, neighbors,
                                 weight_func, dist, data, columns, skip, id) {
  recipes::step(
    subclass = "spatial_lag",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    neighbors = neighbors,
    weight_func = weight_func,
    dist = dist,
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
    dist = x$dist,
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
      dist = object$dist
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
  if (recipes::is_trained(x)) {
    res <- tibble(terms = x$terms)
  } else {
    term_names <- recipes::sel2char(x$terms)
    res <- tibble(terms = rlang::na_chr)
  }
  res$id <- x$id
  res$neighbors <- x$neighbors
  res$weight_func <- x$weight_func
  res
}


#' @export
tunable.step_spatial_lag <- function(x, ...) {
  tibble(
    name = c("neighbors", "weight_func"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1, 10)),
      list(pkg = "dials", fun = "weight_func")
    ),
    source = "recipe",
    component = "step_spatial_lag",
    component_id = x$id
  )
}
