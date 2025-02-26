#' Calculates rates along model in/out flows from compartments

#' Helper function to calculate modified change rates for in and out flows in
#' a stochastic generalized compartmental model using [adaptivetau::ssa.exact()]
#' @param rate_eqns a list of the transition rate equations between flows,
#'                  expressed with dplyr::expr() and listed in the same order
#'                  as the transitions supplied to the model runner function
#' @return the rate(s) at a specific time point
#' @export

generalized_rates <- function(rate_eqns) {
  function(state, params, t) {
    sapply(rate_eqns, function(eq) {
      eval(parse(text = eq), envir = as.list(c(state, params, t)))
    })
  }
}
