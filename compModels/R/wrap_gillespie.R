#' Wrapper to run Gillespie algorithm with [GillespieSSA::ssa()] from compiled
#' model and formatted inputs
#'
#' @param init_vals the starting values for populations in each compartment
#' @param compiledmodel compiled model as a list, output from compilemodel()
#' @param parameters vector of parameter values
#' @param n_timesteps number of time steps
#' @param n_sims number of simulations desired
#' @return list of data frames object of the data part of GillespieSSA object
#' containing the time and states of the simulation. Number of elements should
#' match the number of simulations
#' @export
#' @examples
#' \dontrun{
#' # SIR
#' base_states <- c("S", "I", "R")
#' sir <- define_states(base_states) |>
#'   add_infection("I", "S", "I", "beta") |>
#'   add_transition("I", "R", "tau")
#' sircompiled <- compilemodel(sir)
#' wrap_gillespie(
#'   c("S" = 999, "I" = 1, "R" = 0),
#'   sircompiled,
#'   c(beta = 2, tau = 1),
#'   25,
#'   10
#' )
#'
#' # SEIR
#' base_states <- c("S", "E", "I", "R")
#' seir <- define_states(base_states) |>
#'   add_infection("I", "S", "E", "beta") |>
#'   add_transition("I", "R", "tau") |>
#'   add_transition("E", "I", "taue")
#' seircompiled <- compilemodel(seir)
#' test2 <- wrap_gillespie(
#'   c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
#'   seircompiled,
#'   c(beta = 2, tau = 1, taue = .5),
#'   25,
#'   10
#' )
#' }
wrap_gillespie <- function(init_vals, compiledmodel, parameters, n_timesteps,
                           n_sims) {
  x0 <- define_initialstate(
    compiledmodel,
    init_vals
  ) |>
    output_initialstate()
  model_rates <- compiledmodel$modeloutstructions$processrates # propensity
  model_peter <- compiledmodel$modeloutstructions$petermatrix # change matrix
  parameters <- parameters

  nu <- matrix(as.numeric(model_peter),
    ncol = dim(model_peter)[2],
    nrow = dim(model_peter)[1]
  )
  t <- n_timesteps

  sims <- lapply(1:n_sims, function(i) {
    current_parameters <- parameters

    sim_data <- GillespieSSA::ssa(
      x0 = x0,
      a = model_rates,
      nu = nu,
      parms = current_parameters,
      tf = t,
      method = GillespieSSA::ssa.d()
    )$data |> data.frame()
  })
}
