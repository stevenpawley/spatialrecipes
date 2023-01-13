#' Euclidean distance calculation
#'
#' @param x1 A numeric vector of locations
#' @param x2 A reference location with the same number of dimensions as x1
#'
#' @return Will return a numeric vector for a single reference location
#' @keywords internal
#'
#' @examples
#' x1 <- matrix(c(0,0,0,1,1,1,2,2,2), nrow = 3, byrow = TRUE)
#' x2 <- c(1, 1, 1)
#' euc_dist(x1, x2)
euc_dist <- function(x1, x2) {
  # check that dimensions of locations and reference are equal
  data_dim <- dim(x1)[2]
  reference_dim <- length(x2)

  if (reference_dim != data_dim) {
    rlang::abort("Dimensionality of reference locations does not match the data")
  }

  # euclidean distance calculation to a single reference location
  apply(x1, 1, function(x1, x2) {
    sqrt(sum((x1 - x2) ^ 2))
  }, x2 = x2)
}
