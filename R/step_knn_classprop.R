knn_class_proportions <- function(formula, x, y, k) {
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

  # get values of neighbors
  neighbor_vals <- x[as.integer(neighbor_ids), ][[target_variable]]
  neighbor_vals <- matrix(neighbor_vals, ncol = k)

  # calculate proportions of each class
  class_variables <- levels(x[[target_variable]])

  df <- as.data.frame(matrix(nrow = 1, ncol = length(class_variables)))
  names(df) <- class_variables

  props <- apply(neighbor_vals, 1, function(x) {
    row_prop <- prop.table(table(x))
    df[1, names(row_prop)] <- row_prop
    df
  })

  props <- do.call(rbind, props)
  props[is.na(props)] <- 0
  names(props) <- paste("prop", names(props), sep = "_")

  props
}


#' Spatial class proportions
#'
#' `step_knn_classprop` creates a *specification* of a recipe step that will add
#' a new features to a dataset based on the proportion of each class in the
#' surrounding observations.
#'
#' @param recipe A recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param class Selector function to choose which variable will be used to
#'   create a new feature based proportion of occurrencs of each class within
#'   the spatial neighborhood.
#' @param role role or model term created by this step, what analysis role
#'   should be assigned?. By default, the function assumes that resulting
#'   distance will be used as a predictor in a model.
#' @param trained A logical that will be updated once the step has been trained.
#' @param neighbors The number of closest neighbours to use in the distance
#'   weighting.
#' @param data Used internally to store the training data.
#' @param columns A character string that contains the names of columns used in the
#' transformation. This is `NULL` until computed by `prep.recipe()`.
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
#'   step_knn_classprop(Latitude, Longitude, class = "Sale_Price", neighbors = 3)
#'
#' prepped <- prep(rec_obj)
#' juice(prepped)
step_knn_classprop <- function(recipe, ...,
                               class = NULL,
                               role = "predictor",
                               trained = FALSE,
                               neighbors = 3,
                               data = NULL,
                               columns = NULL,
                               skip = FALSE,
                               id = recipes::rand_id("knn_classprop")) {
  recipes::recipes_pkg_check("nabor")

  terms <- recipes::ellipse_check(...)

  if (neighbors <= 1) {
    rlang::abort("`neighbors` should be greater than 1.")
  }

  recipes::add_step(
    recipe,
    step_knn_classprop_new(
      terms = terms,
      class = rlang::enquos(class),
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
step_knn_classprop_new <- function(terms, role, trained, class, neighbors,
                                   data, columns, skip, id) {
  recipes::step(
    subclass = "knn_classprop",
    terms = terms,
    role = role,
    trained = trained,
    class = class,
    neighbors = neighbors,
    data = data,
    columns = columns,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_knn_classprop <- function(x, training, info = NULL, ...) {
  # First translate the terms argument into column name
  col_names <- recipes::terms_select(terms = x$terms, info = info)
  class_name <- recipes::terms_select(x$class, info = info)

  # Use the constructor function to return the updated object
  # Note that `trained` is set to TRUE
  step_knn_classprop_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    class = class_name,
    neighbors = x$neighbors,
    data = training,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_knn_classprop <- function(object, new_data, ...) {
  f <- as.formula(
    paste(object$class, paste(object$columns, collapse = " + "), sep = " ~ ")
  )

  lags <- knn_class_proportions(
    formula = f,
    x = object$data,
    y = new_data,
    k = object$neighbors
  )

  new_data <- dplyr::bind_cols(new_data, lags)
  new_data
}

#' @export
print.step_knn_classprop <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Spatial lags")

    if (recipes::is_trained(x)) {
      cat(paste0(" (", x$neighbors, " neighbors)"))
    }
    cat("\n")

    invisible(x)
  }

#' @rdname step_knn_classprop
#' @param x A `step_knn_classprop` object.
#' @export
tidy.step_knn_classprop <- function(x, ...) {
  term_names <- sel2char(x$terms)

  res <- tibble(
    id = rep(x$id, times = length(term_names)),
    terms = term_names,
    neighbors = rep(x$neighbors, times = length(term_names))
  )

  res
}

#' @export
tunable.step_knn_classprop <- function(x, ...) {
  tibble::tibble(
    name = c("neighbors", "weight_func", "dist_power"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1L, 10L))
    ),
    source = "recipe",
    component = "step_knn_classprop",
    component_id = x$id
  )
}
