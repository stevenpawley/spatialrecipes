#' Compute which points in `y` are the k-neighbors of `x`
#' 
#' If `x` and `y` represent identical datasets then the points where x[i] ==
#' y[i] are not considered neighbors, i.e. the query point is not considered to
#' be its own neighbor.
#'
#' @param formula A formula.
#' @param x A data.frame of training data.
#' @param y A data.frame of data to predict, can be the same as `x`.
#' @param k An integer specifying the number of neighbours to use.
#'
#' @return
return_neighbours <- function(formula, x, y, k) {
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
  
  # get values and distances of neighbours
  idx <- neighbors$nn.idx
  D <- neighbors$nn.dists
  W <- x[as.numeric(idx), ][[target_variable]]
  W <- matrix(W, ncol = k)
  
  # return as tibble
  prefix <- paste("nn", target_variable, sep = "_")
  colnames(W) <- paste(prefix, 1:ncol(W), sep = "_")
  colnames(D) <- paste("dist", 1:ncol(D), sep = "_")
  
  W <- as_tibble(W)
  D <- as_tibble(D)

  bind_cols(W, D)
}


#' Spatial si step
#' 
#' `step_spatial_si` creates a *specification* of a recipe step that will add a new
#' 'si' features to a dataset based on the inverse distance-weighted mean of
#' surrounding observations.
#'
#' @param recipe A recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param outcome Selector function to choose which variable will be used to
#'   create a new feature based on the inverse distance-weighted mean of
#'   surrounding observations.
#' @param role role or model term created by this step, what analysis
#'  role should be assigned?. By default, the function assumes
#'  that resulting distance will be used as a predictor in a model.
#' @param trained A logical that will be updated once the step has been trained.
#' @param neighbors The number of closest neighbours to use in the distance
#'   weighting.
#' @param data Used internally to store the training data.
#' @param skip A logical to skip training.
#' @param id An identifier for the step. If omitted then this is generated
#' automatically.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any).
#' @export
step_spatial_si <- function(
  recipe, ...,
  outcome = NULL,
  role = "predictor",
  trained = FALSE,
  neighbors = NA,
  data = NULL,
  columns = NULL,
  skip = FALSE,
  id = rand_id("spatial_si")) {
  
  if (!"nabor" %in% installed.packages()[, 1])
    stop("step_infgain requires the package `nabor` to be installed")
  
  terms <- ellipse_check(...)
  
  if (neighbors <= 0)
    rlang::abort("`neighbors` should be greater than 0.")
  
  recipes::add_step(
    recipe,
    step_spatial_si_new(
      terms = terms,
      outcome = rlang::enquos(outcome),
      trained = trained,
      role = role,
      neighbors = neighbors,
      data = data,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}


# wrapper around 'step' function that sets the class of new step objects
#' @importFrom recipes step
step_spatial_si_new <- function(terms, role, trained, outcome, neighbors, 
                                 data, columns, skip, id) {
  recipes::step(
    subclass = "spatial_si",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    neighbors = neighbors,
    data = data,
    columns = columns,
    skip = skip,
    id = id
  )
}


prep.step_spatial_si <- function(x, training, info = NULL, ...) {
  
  # First translate the terms argument into column name
  col_names <- terms_select(terms = x$terms, info = info)
  outcome_name <- terms_select(x$outcome, info = info)
  
  # Use the constructor function to return the updated object
  # Note that `trained` is set to TRUE
  step_spatial_si_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    outcome = outcome_name,
    neighbors = x$neighbors,
    data = training,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}


bake.step_spatial_si <- function(object, new_data, ...) {
  
  f <- as.formula(
    paste(object$outcome, paste(object$columns, collapse = " + "), sep = " ~ ")
  )
  
  lags <-
    return_neighbours(
      formula = f,
      x = object$data,
      y = new_data,
      k = object$neighbors
    )
  
  new_data <- dplyr::bind_cols(new_data, lags)
  new_data
}


tidy.step_spatial_si <- function(x, ...) {
  res <- tibble::tibble(
    neighbors = x$neighbors
  )
  res
}


tunable.step_spatial_si <- function(x, ...) {
  tibble::tibble(
    name = c("neighbors"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1, 10))
    ),
    source = "recipe",
    component = "step_spatial_si",
    component_id = x$id
  )
}
