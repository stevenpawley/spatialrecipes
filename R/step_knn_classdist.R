return_knn_classdists <- function(formula, x, y, k) {
  # get the classes of the target variable
  target_variable <-
    rlang::f_lhs(formula) %>%
    as.character()

  classes <- levels(x[[target_variable]])
  term_variables <- attr(terms(formula), "term.labels")

  # split data
  train_data <- x[term_variables]
  query_data <- y[term_variables]

  # store distances to neighbors in each class
  class_names <- rep(classes, each = k)
  neighbor_names <- paste0("neighbor", seq_len(k))
  feature_names <- paste(class_names, neighbor_names, sep = "_")

  class_dists <- matrix(nrow = nrow(query_data), ncol = length(feature_names))
  colnames(class_dists) <- feature_names

  for (clf in classes) {
    query_rows <- x[[target_variable]] == clf

    if (identical(train_data, query_data)) {
      neighbors <-
        nabor::knn(data = train_data[query_rows, ], query = query_data, k = k + 1)
      neighbors$nn.idx <- neighbors$nn.idx[, 2:ncol(neighbors$nn.idx)]
      neighbors$nn.dists <- neighbors$nn.dists[, 2:ncol(neighbors$nn.dists)]

    } else {
      neighbors <-
        nabor::knn(data = train_data[query_rows, ], query = query_data, k = k)
    }

    # get values and distances of neighbors
    D <- matrix(neighbors$nn.dists, ncol = k)
    columns <- paste(clf, paste(paste0("neighbor", seq_len(k))), sep = "_")
    class_dists[, columns] <- D
  }

  as_tibble(class_dists)
}


#' Spatial class distances step
#'
#' `step_knn_classdist` creates a *specification* of a recipe step that will add
#' new features to a dataset based on the distances to each neighbor, per class.
#'
#' For example, if a dataset's outcome variable contains two classes, c("sand",
#' "clay"), this step will append a new feature matrix with (n_samples,
#' n_classes * neighbors) to the dataset. For the previous example with
#' neighbors = 3, the step would create new features called "sand_neighbor1,
#' "sand_neighbor2, "sand_neighbor3", "clay_neighbor1", "clay_neighbor2",
#' "clay_neighbor3" with each feature representing the distance (in the feature
#' space) to the closest k points of each class.
#'
#' This step is based on the concept described in the 'fastknn' package
#' (https://github.com/davpinto/fastknn) to perform a non-linear mapping of the
#' similarity between classes and reproject it into a linear space.
#'
#' @param recipe A recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used.
#' @param class Selector function to choose which variable will be used to
#'   create features based on the distances to each neighbor per class.
#' @param neighbors The number of closest neighbors to use in the distance
#'   weighting. Default is 3.
#' @param role role or model term created by this step, what analysis role
#'   should be assigned?. By default, the function assumes that resulting
#'   distance will be used as a predictor in a model.
#' @param trained A logical that will be updated once the step has been trained.
#' @param data Used internally to store the training data.
#' @param columns A character string that contains the names of columns used in
#'   the transformation. This is `NULL` until computed by `prep.recipe()`.
#' @param skip A logical to skip training.
#' @param id An identifier for the step. If omitted then this is generated
#'   automatically.
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any).
#' @export
step_knn_classdist <- function(
  recipe, ...,
  class = NULL,
  role = "predictor",
  neighbors = 3,
  trained = FALSE,
  data = NULL,
  columns = NULL,
  skip = FALSE,
  id = recipes::rand_id("knn_classdist")) {

  recipes::recipes_pkg_check("nabor")
  terms <- recipes::ellipse_check(...)

  if (neighbors <= 0)
    rlang::abort("`neighbors` should be greater than 0.")

  recipes::add_step(
    recipe,
    step_knn_classdist_new(
      terms = terms,
      class = class,
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
step_knn_classdist_new <- function(terms, role, trained, class, neighbors,
                                   data, columns, skip, id) {
  recipes::step(
    subclass = "knn_classdist",
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
prep.step_knn_classdist <- function(x, training, info = NULL, ...) {

  # First translate the terms argument into column name
  col_names <- terms_select(terms = x$terms, info = info)
  class_name <- x$class

  # Check selections
  if (!inherits(training[[class_name]], "factor"))
    rlang::abort("The selected `class` must be a factor variable")

  col_types <- sapply(col_names, function(nm) inherits(training[[nm]], "numeric"))
  if (!all(col_types))
    rlang::abort("`step_knn_classdist` can only be applied to numeric features")

  n_per_class <- table(training[[class_name]])

  if (!all(n_per_class >= x$neighbors + 1))
    rlang::abort("There needs to be >= neighbors+1 cases")

  # Use the constructor function to return the updated object
  # Note that `trained` is set to TRUE
  step_knn_classdist_new(
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
bake.step_knn_classdist <- function(object, new_data, ...) {
  f <- as.formula(
    paste(object$class, paste(object$columns, collapse = " + "), sep = " ~ "))

  new_X <- return_knn_classdists(
    formula = f,
    x = object$data,
    y = new_data,
    k = object$neighbors
  )

  new_data <- dplyr::bind_cols(new_data, new_X)
  new_data
}

#' @export
tidy.step_knn_classdist <- function(x, ...) {
  res <- tibble::tibble(
    terms = recipes::sel2char(x$terms),
    class = x$class,
    neighbors = x$neighbors
  )
  res
}

#' @export
tunable.step_knn_classdist <- function(x, ...) {
  tibble::tibble(
    name = c("neighbors"),
    call_info = list(
      list(pkg = "dials", fun = "neighbors", range = c(1, 10))
    ),
    source = "recipe",
    component = "step_knn_classdist",
    component_id = x$id
  )
}
