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
