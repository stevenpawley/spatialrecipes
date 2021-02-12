#' Parameter functions for spatial feature engineering recipes
#'
#' New features derived from spatial distances require normalization. `norm` is
#' for specifying the type of row-normalization.
#'
#' @param values A character string of possible values. Either "l1", "l2", or
#' "max".
#'
#' @return A function with classes "qual_param" and "param"
#' @export
#'
#' @examples
#' norm("l1")
norm <- function(values = values_norm_func) {
  dials::new_qual_param(
    type = "character",
    values = values,
    label = c(norm = "Normalization method"),
    finalize = NULL
  )
}


values_norm_func <- c("l1", "l2", "max")
