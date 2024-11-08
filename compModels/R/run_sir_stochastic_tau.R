#' Runs SIR stochastic compartmental model with tau leaping
#'
#' Runs basic stochastic SIR model using [adaptivetau::ssa.exact()] given
#' parameters, initial conditions, number of timesteps, and state change matrix.
#' No seed is set in the function so results are expected to vary each time
#' the function is run.
#'
#' @param beta transmission rate
#' @param gamma recovery rate
#' @param s0 initial count of susceptibles
#' @param i0 initial count of infecteds
#' @param r0 initial count of recovereds
#' @param n_pop population size
#' @param n_timesteps number of time steps
#' @param transitions list object with transition information between states
#' @param n_sims number of simulations to run
#' @return data frame object of the time and states of the simulation
#' @export
#' @examples
#' \dontrun{
#' run_sir_stochastic_tau(
#'   0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
#'   list(c(s = -1, i = 1), c(i = -1, r = 1)),
#'   10
#' )
#' }
run_sir_stochastic_tau <- function(beta, gamma, s0, i0, r0, n_pop,
                                   n_timesteps, transitions, n_sims) {
  validate_sir_stoch_input(
    beta, gamma, s0, i0, r0, n_pop, n_timesteps,
    method = "adaptivetau",
    transitions = transitions,
    n_sims = n_sims
  )

  parms <- list(beta = beta, gamma = gamma)
  init <- c(s = s0, i = i0, r = r0)
  sir_transitions <- transitions # Transitions function

  t <- n_timesteps

  # Rate Function: sir_rates() stand alone function in compModels package

  out_stoch2 <- lapply(1:n_sims, function(i) {
    adaptivetau::ssa.exact(
      init.values = init,
      transitions = sir_transitions,
      rateFunc = sir_rates,
      params = parms,
      tf = t
    )
  })

  lapply(out_stoch2, function(x) data.frame(x))
}
