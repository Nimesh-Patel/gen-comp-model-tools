#' Randomly delay times according to exponential distribution
#'
#' @param t vector of numerical times
#' @param meandelay average time to add, sets exponential distribution to sample
#' from
#' @param sortlogic optional logic whether to sort the new times. Sorting is
#' useful for outputting time series. Not sorting is useful to keep the order
#' according to the latent counts associated with the times
#' default is TRUE which designates sorting
#' @return updated times
#' @export
delaytime_exp <- function(t, meandelay, sortlogic = TRUE) {
  delayt <- t + stats::rexp(length(t), 1 / meandelay)
  if (sortlogic) {
    delayt <- sort(delayt)
  }

  delayt
}
