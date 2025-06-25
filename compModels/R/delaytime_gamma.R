#' Randomly delay times according to Gamma distribution
#'
#' Note: if an integer shape is input then the sampled gamma distribution
#' matches an erlang distribution, which is the distribution of times resulting
#' from chained transition processes.
#'
#' @param t vector of numerical times
#' @param meandelay average time to add, sets mean of erlang distribution to
#' sample from
#' @param shape number of steps in for the analogous chained transition
#' @param sortlogic optional logic whether to sort the new times. Sorting is
#' useful for outputting time series. Not sorting is useful to keep the order
#' according to the latent counts associated with the times
#' default is TRUE which designates sorting
#' @return updated times
#' @export
delaytime_gamma <- function(t, meandelay, shape, sortlogic = TRUE) {
  delayt <- t + stats::rgamma(length(t), shape, scale = meandelay / shape)
  if (sortlogic) {
    delayt <- sort(delayt)
  }

  delayt
}
