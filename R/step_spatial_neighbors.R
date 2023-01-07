return_neighbours <- function(formula, x, y, k, prefix = "nn") {
  target_variable <-
    rlang::f_lhs(formula) %>%
    as.character()

  term_variables <- attr(terms(formula), "term.labels")

  # split data
  train_data <- x[term_variables]
  query_data <- y[term_variables]

  if (identical(train_data, query_data)) {
    neighbors <- nabor::knn(data = train_data, query = query_data, k = k + 1)
    neighbors$nn.idx <- neighbors$nn.idx[, 2:ncol(neighbors$nn.idx)]
    neighbors$nn.dists <- neighbors$nn.dists[, 2:ncol(neighbors$nn.dists)]
  } else {
    neighbors <- nabor::knn(data = train_data, query = query_data, k = k)
  }

  # get values and distances of neighbors
  idx <- neighbors$nn.idx
  D <- neighbors$nn.dists
  W <- x[as.numeric(idx), ][[target_variable]]

  W <- matrix(W, ncol = k)
  D <- matrix(D, ncol = k)

  # return as tibble
  prefix <- paste(prefix, target_variable, sep = "_")
  prefix_dist <- paste(prefix, "dist", target_variable, sep = "_")
  colnames(W) <- paste(prefix, 1:ncol(W), sep = "_")
  colnames(D) <- paste(prefix_dist, 1:ncol(D), sep = "_")

  W <- as_tibble(W)
  D <- as_tibble(D)

  bind_cols(W, D)
}


#' Spatial neighbors step
#'
#' `step_spatial_neighbors` creates a *specification* of a recipe step that will
#' add a new 'si' features to a dataset based on the inverse distance-weighted
#' mean of surrounding observations.
#'
#' @param recipe A recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param outcome Selector function to choose which variable will be used to
#'   create a new feature based on the inverse distance-weighted mean of
#'   surrounding observations.
#' @param neighbors The number of closest neighbors to use in the distance
#'   weighting. The default is 3.
#' @param prefix Prefix to use for the newly created variables. Default is "nn".
#' @param role role or model term created by this step, what analysis role
#'   should be assigned?. By default, the function assumes that resulting
#'   distance will be used as a predictor in a model.
#' @param trained A logical that will be updated once the step has been trained.
#' @param data Used internally to store the training data.
#' @param columns A character string that contains the names of columns used in
#'   the transformation. This is `NULL` until computed by `prep.recipe()`.
#' @param skip A logical to skip training.
#' @param id An identifier for the step. If omitted then this is generated
#' automatically.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any).
#' @export
step_spatial_neighbors <- function(recipe, ...,
                                   outcome = NULL,
                                   role = "predictor",
                                   trained = FALSE,
                                   neighbors = 3,
                                   prefix = "nn",
                                   data = NULL,
                                   columns = NULL,
                                   skip = FALSE,
                                   id = recipes::rand_id("spatial_neighbors")) {
  if (!"nabor" %in% installed.packages()[, 1]) {
    stop("step_spatial_neighbors requires the package `nabor` to be installed")
  }

  recipes::recipes_pkg_check("nabor")
  terms <- recipes::ellipse_check(...)

  if (neighbors <= 0) {
    rlang::abort("`neighbors` should be greater than 0.")
  }

  recipes::add_step(
    recipe,
    step_spatial_neighbors_new(
      terms = terms,
      outcome = rlang::enquos(outcome),
      neighbors = neighbors,
      prefix = prefix,
      trained = trained,
      role = role,
      data = data,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}


# wrapper around 'step' function that sets the class of new step objects
step_spatial_neighbors_new <- function(terms, role, trained, outcome, neighbors,
                                       prefix, data, columns, skip, id) {
  recipes::step(
    subclass = "spatial_neighbors",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    neighbors = neighbors,
    prefix = prefix,
    data = data,
    columns = columns,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_spatial_neighbors <- function(x, training, info = NULL, ...) {
  # First translate the terms argument into column name
  col_names <- recipes_eval_select(x$terms, training, info)
  outcome_name <- recipes_eval_select(x$outcome, training, info)

  # Use the constructor function to return the updated object
  # Note that `trained` is set to TRUE
  step_spatial_neighbors_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    outcome = outcome_name,
    neighbors = x$neighbors,
    prefix = x$prefix,
    data = training,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_spatial_neighbors <- function(object, new_data, ...) {
  f <- as.formula(
    paste(object$outcome, paste(object$columns, collapse = " + "), sep = " ~ ")
  )

  new_X <- return_neighbours(
    formula = f,
    x = object$data,
    y = new_data,
    k = object$neighbors,
    prefix = object$prefix
  )

  new_data <- dplyr::bind_cols(new_data, new_X)
  new_data
}

#' @export
tidy.step_spatial_neighbors <- function(x, ...) {
  res <- tibble::tibble(
    terms = recipes::sel2char(x$terms),
    outcome = x$outcome,
    neighbors = x$neighbors
  )

  res
}

#' @export
tunable.step_spatial_neighbors <- function(x, ...) {
  tibble::tibble(
    name = c("neighbors"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1, 10))
    ),
    source = "recipe",
    component = "step_spatial_neighbors",
    component_id = x$id
  )
}
