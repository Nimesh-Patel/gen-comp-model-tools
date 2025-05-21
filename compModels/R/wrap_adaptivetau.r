#' Wrapper to run adaptivetau functions [adaptivetau::ssa.exact()] or
#' [adaptivetau::ssa.adaptivetau()] from compiled model and formatted inputs
#'
#' @param init_vals the starting values for populations in each compartment
#' @param compiledmodel compiled model as a list, output from compilemodel()
#' @param rate_func function to calculate rates, defaults to NULL which uses
#' [compModels::generalized_rates()]
#' @param parameters vector of parameter values
#' @param n_timesteps number of time steps
#' @param n_sims number of simulations desired
#' @param method method to use in adaptivetau, defaults to "exact". Other option
#' is "adaptivetau"
#' @param usestates "initial" is default option to use the user supplied
#' init_vals vector as as initial values, "current" option allows init_vals to
#' be used as the current values  of the simulation for when the user needs to
#' start from a spot other than initial (e.g. stitching together simulations
#' with different time blocks of interventions)
#' @return list of data frames object of the data part of adaptivetau object
#' containing the time and states of the simulation. Number of elements should
#' match the number of simulations
#' @export
#' @examples
#' \dontrun{
#' # SIR
#' # Note this example does not name each argument so rate_func = NULL is needed
#' # for R to recognize the correct arguments in order
#' base_states <- c("S", "I", "R")
#' sir <- define_states(base_states) |>
#'   add_infection("I", "S", "I", "beta") |>
#'   add_transition("I", "R", "tau")
#' sircompiled <- compilemodel(sir)
#' wrap_adaptivetau(
#'   c("S" = 999, "I" = 1, "R" = 0),
#'   sircompiled,
#'   rate_func = NULL,
#'   c(beta = 2, tau = 1),
#'   25,
#'   10
#' )
#'
#' # SIR with tau leaping
#' base_states <- c("S", "I", "R")
#' sir <- define_states(base_states) |>
#'   add_infection("I", "S", "I", "beta") |>
#'   add_transition("I", "R", "tau")
#' sircompiled <- compilemodel(sir)
#' wrap_adaptivetau(
#'   c("S" = 999, "I" = 1, "R" = 0),
#'   sircompiled,
#'   rate_func = NULL,
#'   c(beta = 2, tau = 1),
#'   25,
#'   10,
#'   "adaptivetau"
#' )
#'
#' # SEIR
#' # Note: this example names each argument so rate_func can be left out
#' base_states <- c("S", "E", "I", "R")
#' seir <- define_states(base_states) |>
#'   add_infection("I", "S", "E", "beta") |>
#'   add_transition("I", "R", "tau") |>
#'   add_transition("E", "I", "taue")
#' seircompiled <- compilemodel(seir)
#' test2 <- wrap_adaptivetau(
#'   c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
#'   seircompiled,
#'   rate_func = NULL,
#'   c(beta = 2, tau = 1, taue = .5),
#'   25,
#'   10
#' )
#' }
wrap_adaptivetau <- function(init_vals, compiledmodel, rate_func = NULL,
                             parameters, n_timesteps, n_sims,
                             method = "exact", usestates = "initial") {
  if (usestates == "initial") {
    x0 <- define_initialstate(
      compiledmodel,
      init_vals
    ) |>
      output_initialstate()
  } else if (usestates == "current") {
    x0 <- define_currentstate(
      compiledmodel,
      init_vals # filled with current values by user
    ) |>
      output_initialstate()
  }

  # rateeqns, transitions, parameters
  model_rates <- as.list(compiledmodel$modeloutstructions$processrates)
  model_peter <- as.matrix(compiledmodel$modeloutstructions$petermatrix)
  parameters <- parameters

  t <- n_timesteps

  if (is.null(rate_func)) {
    rate_func <- generalized_rates
  }

  if (method == "exact") {
    sims <- lapply(1:n_sims, function(i) {
      adaptivetau::ssa.exact(
        init.values = x0,
        transitions = model_peter,
        rateFunc = rate_func(model_rates),
        params = parameters,
        tf = t
      ) |> data.frame()
    })
  } else if (method == "adaptivetau") {
    sims <- lapply(1:n_sims, function(i) {
      adaptivetau::ssa.adaptivetau(
        init.values = x0,
        transitions = model_peter,
        rateFunc = rate_func(model_rates),
        params = parameters,
        tf = t
      ) |> data.frame()
    })
  } else {
    stop("Method not recognized. Please use 'exact' or 'adaptivetau'")
  }
}
