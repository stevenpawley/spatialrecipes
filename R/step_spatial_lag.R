knn_train <- function(formula, x, y, k, weight_func, prefix = "lag") {
  target_variable <-
    rlang::f_lhs(formula) %>%
    as.character()

  term_variables <- attr(terms(formula, data = x), "term.labels")

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
  D[D < 1e-6] <- 1e-6

  # get values of neighbors
  neighbor_vals <- x[as.integer(neighbor_ids), ][[target_variable]]
  neighbor_vals <- matrix(neighbor_vals, ncol = k)

  # calculate weights
  W <- switch(weight_func,
    inv = 1 / D,
    rectangular = matrix(1, nrow = nrow(neighbor_vals), ncol = k),
    gaussian = gaussian_kernel(D, k)
  )

  # calculate weighted mean/mode of neighbors
  if (inherits(x[[target_variable]], "numeric")) {
    denom <- rowSums(W)
    num <- rowSums(neighbor_vals * W)
    fitted <- num / denom
  } else {
    fitted <- sapply(seq_len(nrow(W)), function(i) {
      collapse::fmode(x = neighbor_vals[i, ], w = W[i, ])
    })

    target_type <- typeof(x[[target_variable]])
    fitted <- as(fitted, target_type)
    fitted <- factor(fitted, levels = levels(x[[target_variable]]))
  }

  new_feature_name <- paste(target_variable, prefix, k, weight_func, sep = "_")

  fitted <- tibble(fitted) %>%
    rlang::set_names(new_feature_name)

  fitted
}


#' Spatial lag step
#'
#' `step_spatial_lag` creates a *specification* of a recipe step that will add a
#' new 'lag' feature to a dataset based on the weighted or unweighted mean or
#' mode of neighbouring observations.
#'
#' @param recipe A recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param outcome Character, name of which variable will be used to create a new
#'   feature based on the unweighted or distance-weighted mean of surrounding
#'   observations.
#' @param role role or model term created by this step, what analysis role
#'   should be assigned?. By default, the function assumes that resulting
#'   distance will be used as a predictor in a model.
#' @param trained A logical that will be updated once the step has been trained.
#' @param neighbors The number of closest neighbours to use in the distance
#'   weighting.
#' @param weight_func A single character for the kernel function used to weight
#'   the distances between samples. The default is 'rectangular' and the
#'   available choices are 'rectangular', 'inv', 'gaussian'.
#' @param prefix Name to use as prefix for the created variable. Default is
#'   'lag.
#' @param data Used internally to store the training data.
#' @param columns A character string that contains the names of columns used in
#'   the transformation. This is `NULL` until computed by `prep.recipe()`.
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
#'   recipe(Sale_Price ~ Latitude + Longitude) %>%
#'   step_spatial_lag(Latitude, Longitude,
#'     outcome = "Sale_Price", neighbors = 3,
#'     weight_func = "inv"
#'   )
#'
#' prepped <- prep(rec_obj)
#' juice(prepped)
step_spatial_lag <- function(recipe, ...,
                             outcome = NULL,
                             role = "predictor",
                             trained = FALSE,
                             neighbors = 3,
                             weight_func = "rectangular",
                             prefix = "lag",
                             data = NULL,
                             columns = NULL,
                             skip = FALSE,
                             id = recipes::rand_id("spatial_lag")) {
  recipes::recipes_pkg_check("nabor")

  terms <- recipes::ellipse_check(...)

  if (neighbors <= 0) {
    rlang::abort("`neighbors` should be greater than 0.")
  }

  recipes::add_step(
    recipe,
    step_spatial_lag_new(
      terms = terms,
      outcome = outcome,
      trained = trained,
      role = role,
      neighbors = neighbors,
      weight_func = weight_func,
      prefix = prefix,
      data = data,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

# wrapper around 'step' function that sets the class of new step objects
step_spatial_lag_new <- function(terms, role, trained, outcome, neighbors,
                                 weight_func, prefix, data, columns, skip, id) {
  recipes::step(
    subclass = "spatial_lag",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    neighbors = neighbors,
    weight_func = weight_func,
    prefix = prefix,
    data = data,
    columns = columns,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_spatial_lag <- function(x, training, info = NULL, ...) {
  # First translate the terms argument into column name
  col_names <- recipes::recipes_eval_select(x$terms, training, info)

  # Use the constructor function to return the updated object
  # Note that `trained` is set to TRUE
  step_spatial_lag_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    outcome = x$outcome,
    neighbors = x$neighbors,
    weight_func = x$weight_func,
    prefix = x$prefix,
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

  lags <- knn_train(
    formula = f,
    x = object$data,
    y = new_data,
    k = object$neighbors,
    weight_func = object$weight_func,
    prefix = object$prefix
  )

  new_data <- dplyr::bind_cols(new_data, lags)
  new_data
}

#' @export
print.step_spatial_lag <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Spatial lags")
    cat(paste0(" (", x$neighbors, " neighbors)"))
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
    weight_func = rep(x$weight_func, times = length(term_names))
  )

  res
}

#' @export
tunable.step_spatial_lag <- function(x, ...) {
  tibble::tibble(
    name = c("neighbors", "weight_func"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1L, 10L)),
      list(pkg = "dials", fun = "weight_func", values = c("rectangular", "inv", "gaussian"))
    ),
    source = "recipe",
    component = "step_spatial_lag",
    component_id = x$id
  )
}
