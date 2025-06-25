#' Delay times according to function
#'
#' Function can be deterministic (e.g., a single number) or stochastic (e.g.,
#' sampling from a distribution like rgamma()).
#'
#' @param t vector of numerical times
#' @param delayfunction non-negative numeric or function that outputs
#' non-negative value to add to time
#' If a function, it should not have any input arguments
#' @param sortlogic optional logic whether to sort the new times. Sorting is
#' useful for outputting time series. Not sorting is useful to keep the order
#' according to the latent counts associated with the times
#' default is TRUE which designates sorting
#' @return updated times
#' @export
delaytime_byfunction <- function(t, delayfunction, sortlogic = TRUE) {
  if (is.numeric(delayfunction)) {
    addfun <- replicate(length(t), delayfunction)
  } else {
    addfun <- replicate(length(t), delayfunction())
  }

  if (min(addfun) < 0) {
    stop("Delayed values must be non-negative. Check input function")
  }
  delayt <- t + addfun
  if (sortlogic) {
    delayt <- sort(delayt)
  }

  delayt
}
