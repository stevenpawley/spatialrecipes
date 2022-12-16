gaussian_kernel <- function(distances, k) {
  # standardize distances [0, 1] relative to the maximum distance for each obs
  maxdist <- distances[, k]

  # add a small constant to avoid zero distances
  maxdist[maxdist < 1.0e-6] <- 1.0e-6
  W <- distances / maxdist
  W <- pmin(W, 1 - (1e-6))
  W <- pmax(W, 1e-6)

  # create adaptive window width
  # get probability based on the number of neighbours
  # idea is to set sigma equal to the inverse of twice the number of dimensions in the data
  # for example 1 / (2 * 2) = 1/4 for two dimensions
  # for three dimensions, 1 / (2 * 3) = 1/6
  # the idea that the distance between points in high-dimensional space is
  # typically smaller than the distance between points in low-dimensional space.
  # Therefore, a smaller sigma value is needed to capture the similarity between
  # points in high-dimensional space.

  # for KNN using a large number of neighbors, 1 / 2k creates a small kernel
  # bandwith that will be more sensitive to points that are close to the prediction
  # point. Conversely, KNN with a small number of neighbours, 1 / 2k makes a
  # large bandwidth that will smooth out the prediction and will reduce the
  # variance of the KNN model
  p <- 1 / (2 * k)

  # get percentile for probability - how many standard deviations from the mean
  # for a small number of neighbours, this yields a lower percentile that
  # represents a narrower width of the bell-shaped curve, for a large number of
  # neighbours a higher percentile will represent a wider area of the
  # bell-shaped curve so this adjusts the kernel width adaptively
  qua <- abs(qnorm(p))

  # weight normalized distances by the percentile
  Wq <- W * qua

  # use the normal density function to get the density corresponding to the
  # percentile - these are used as the weights
  Wd <- dnorm(Wq, sd = 1)

  return(Wd)
}
