#' Wrapper to run [deSolve::ode()] from compiled model and formatted inputs
#'
#' @param y the initial (state) values for the ODE system, a named vector
#' @param times vector of sequence of times for which output is wanted; first
#' value must be the initial time
#' @param compiledmodel compiled model as a list, output from compilemodel()
#' @param parms named vector of parameters passed to the ode function created by
#' model2desolvefunction using the compiledmodel
#' @return data frame of deSolve output run on provided input
#' @export
#' @examples
#' \dontrun{
#' base_states <- c("S", "I", "R")
#' sir <- define_states(base_states) |>
#'   add_infection("I", "S", "I", "beta") |>
#'   add_transition("I", "R", "tau")
#' sircompiled <- compilemodel(sir)
#' wrap_ode(
#'   c("S" = 999, "I" = 1, "R" = 0),
#'   seq(0, 25, by = 0.01),
#'   sircompiled,
#'   c(beta = 2, tau = 1)
#' )
#' #'
#' }
wrap_ode <- function(y, times, compiledmodel, parms) {
  odefun <- model2desolvefunction(compiledmodel)
  x0 <- define_initialstate(compiledmodel, y) |>
    output_initialstate()
  parameters <- parms
  times <- times
  dyn <- deSolve::ode(y = x0, times = times, func = odefun, parms = parameters)
  data.frame(dyn)
}
