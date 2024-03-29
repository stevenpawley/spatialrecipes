% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step_spatial_classprop.R
\name{step_spatial_classprop}
\alias{step_spatial_classprop}
\alias{tidy.step_spatial_classprop}
\title{Spatial class proportions}
\usage{
step_spatial_classprop(
  recipe,
  ...,
  class = NULL,
  role = "predictor",
  trained = FALSE,
  neighbors = 3,
  data = NULL,
  columns = NULL,
  skip = FALSE,
  id = recipes::rand_id("spatial_classprop")
)

\method{tidy}{step_spatial_classprop}(x, ...)
}
\arguments{
\item{recipe}{A recipe.}

\item{...}{One or more selector functions to choose which variables are
affected by the step. See selections() for more details. For the tidy
method, these are not currently used.}

\item{class}{Selector function to choose which variable will be used to
create a new feature based proportion of occurrencs of each class within
the spatial neighborhood.}

\item{role}{role or model term created by this step, what analysis role
should be assigned?. By default, the function assumes that resulting
distance will be used as a predictor in a model.}

\item{trained}{A logical that will be updated once the step has been trained.}

\item{neighbors}{The number of closest neighbours to use in the distance
weighting.}

\item{data}{Used internally to store the training data.}

\item{columns}{A character string that contains the names of columns used in the
transformation. This is `NULL` until computed by `prep.recipe()`.}

\item{skip}{A logical to skip training.}

\item{id}{An identifier for the step. If omitted then this is generated
automatically.}

\item{x}{A `step_spatial_classprop` object.}
}
\value{
An updated version of `recipe` with the new step added to the
  sequence of existing steps (if any).
}
\description{
`step_spatial_classprop` creates a *specification* of a recipe step that will add
a new features to a dataset based on the proportion of each class in the
surrounding observations.
}
\details{
No details yet!
}
\examples{
library(modeldata)
library(recipes)
library(spatialrecipes)
data(ames)

rec_obj <- ames \%>\%
  recipe(Sale_Price ~ Latitude + Longitude) \%>\%
  step_spatial_classprop(Latitude, Longitude, class = "Sale_Price", neighbors = 3)

prepped <- prep(rec_obj)
juice(prepped)
}
\concept{preprocessing}
\keyword{datagen}
