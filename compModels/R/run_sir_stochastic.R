#' Runs SIR stochastic compartmental model
#'
#' Runs basic stochastic SIR model using [GillespieSSA::ssa.d()] given
#' parameters, initial conditions, number of timesteps, and state change
#' matrix conditions. No seed is set in the function so results are expected
#' to vary each time the function is run.
#'
#' @param beta transmission rate
#' @param gamma recovery rate
#' @param s0 initial count of susceptibles
#' @param i0 initial count of infecteds
#' @param r0 initial count of recovereds
#' @param n_pop population size
#' @param n_timesteps number of time steps
#' @param n_sims number of simulations desired
#' @param propensity_fns provide transition equations as a vector of equations
#' @param change_matrix matrix that gives transition information between states
#' @return list of dataframes of the data part of GillespieSSA object
#' containing the time and states of the simulation. Number of elements should
#' match the number of simulations
#' @export
#' @examples
#' \dontrun{
#' run_sir_stochastic(
#'   0.001, 0.1, 499, 100, 0, 500, 100,
#'   c("beta*s*i", "gamma*i"),
#'   matrix(c(-1, 0, +1, -1, 0, +1),
#'     nrow = 3, byrow = TRUE
#'   ), 10
#' )
#' }
run_sir_stochastic <- function(beta, gamma, s0, i0, r0, n_pop, n_timesteps,
                               propensity_fns, change_matrix, n_sims) {
  validate_sir_stoch_input(
    beta, gamma, s0, i0, r0, n_pop, n_timesteps,
    method = "GillespieSSA",
    propensity_fns = propensity_fns,
    change_matrix = change_matrix,
    n_sims = n_sims
  )

  parms <- c(beta = beta, gamma = gamma)
  init <- c(s = s0, i = i0, r = r0)
  a <- propensity_fns # Propensity Functions
  nu <- change_matrix # State Change Matrix
  t <- n_timesteps

  sims <- lapply(1:n_sims, function(i) {
    GillespieSSA::ssa(
      x0 = init,
      a = a,
      nu = nu,
      parms = parms,
      tf = t,
      method = GillespieSSA::ssa.d(),
      simName = "Stochastic SIR"
    )$data
  })

  lapply(sims, function(x) data.frame(x))
}
