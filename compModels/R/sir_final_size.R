#' Calculates final size of deterministic SIR model
#'
#' Function containing the final size relation equation for a
#' basic SIR model for use with root finding tools. Equation from JC Miller
#' Bull Math Biol (2012) 74:2126. This works for the case when the total
#' population is N = 1.
#'
#' Note: Due to equation being for N = 1. the deprecated compModels function
#' sir_deriv() uses explicit frequency dependent transmission on the flow from
#' susceptibles to infecteds.
#'
#' @param r_infinity Proportion (between 0 and 1) of recovered individuals
#' as time goes toward infinity
#' @param r_0 Reproductive number
#' @param s0 Initial proportion of susceptibles
#' @param r0 Initial proportion of recovereds
#' @return single value vector (typeL double) containing the final outbreak
#' size, eg proportion of recovered individuals of the total population
#' @export
#' @examples
#' \dontrun{
#' sir_final_size(1, 3, 1, 0)
#' sir_final_size(1, 3, .5, .5)
#' sir_final_size(1, 3, .2, .8)
#' }
##############################################################################
# NOTE: This function is PLANNED TO BE MODIFIED / INCORPORATED DIFFERENTLY in
#       this repo as part of the merge of compModels with SIRmodelbuilder
#       functionality. It may be determined that part of the functionality is
#       useful and, if so, it will remain or be incorporated into other
#       functionality.
##############################################################################
sir_final_size <- function(r_infinity, r_0, s0, r0) {
  # A = 1 -s(0) e^-R0(A-r(0)) where A is the final proportion infected (eg
  # the final value of r(infinity))
  1 - s0 * exp(-r_0 * (r_infinity - r0))
}
