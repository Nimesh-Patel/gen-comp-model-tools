#' Checks validity of input to stochastic SIR model
#'
#' Helper function that runs checks on basic SIR stochastic model input
#' parameters, times, and initial conditions
#'
#' @param beta transmission rate
#' @param gamma recovery rate
#' @param s0 initial count of susceptibles
#' @param i0 initial count of infecteds
#' @param r0 initial count of recovereds
#' @param n_pop population size
#' @param n_timesteps number of time steps
#' @param method character name of stochastic method used, default NULL. Example
#' options include "GillespieSSA" or "adaptivetau"
#' @param propensity_fns provide transition equations as a vector of equations,
#' default NULL
#' @param change_matrix matrix that gives transition information between states,
#' default NULL
#' @param transitions list object with transition information between states,
#' default NULL
#' @param n_sims number of simulations desired
#' @return stops with information or indicates Input Checks Passed
#' @export
#' @examples
#' \dontrun{
#' validate_sir_stoch_input(
#'   0.001, 0.1, 499, 100, 0, 500, 100,
#'   method = "GillespieSSA",
#'   c("beta*s*i", "gamma*i"),
#'   matrix(c(-1, 0, +1, -1, 0, +1),
#'     nrow = 3, byrow = TRUE
#'   ), 10
#' )
#'
#' validate_sir_stoch_input(
#'   0.001, 0.1, 499, 100, 0, 500, 100,
#'   method = "adaptivetau",
#'   list(c(s = -1, i = 1), c(i = -1, r = 1)),
#'   10
#' )
#' }
validate_sir_stoch_input <- function(beta, gamma, s0, i0, r0, n_pop,
                                     n_timesteps,
                                     method = NULL,
                                     propensity_fns = NULL,
                                     change_matrix = NULL,
                                     transitions = NULL,
                                     n_sims) {
  # Checks on validity of user input
  sapply(
    X = c("beta", "gamma", "s0", "i0", "r0", "n_pop", "n_timesteps", "n_sims"),
    FUN = function(x) {
      if (!is.numeric(get(x)) | get(x) < 0) {
        stop(x, " must be numeric and non-negative. It is ", get(x))
      }
    }
  )

  if (method == "GillespieSSA") {
    if (!is.character(propensity_fns)) stop("propensity_fns must be character")
    if (!is.matrix(change_matrix)) stop("change_matrix must be a matrix")
  } else if (method == "adaptivetau") {
    if (!is.list(transitions)) stop("transitions must be a list")
    sapply(
      X = transitions,
      FUN = function(x) {
        if (!is.numeric(x)) {
          stop("All elements of transitions list must be numeric")
        }
      }
    )
    sapply(
      X = transitions,
      FUN = function(x) {
        if (!all(!is.null(names(x)))) {
          stop("All elements of transitions list must be named numeric
                    vectors")
        }
      }
    )
  }

  sapply(
    X = c("n_timesteps", "n_sims"),
    FUN = function(x) {
      if (get(x) == 0) {
        stop(x, " must be > 0. It is ", get(x))
      }
    }
  )
  if (length(n_timesteps) > 1) {
    stop("times should contain the number of time steps, not a vector of times")
  }
  if (any(c(s0, i0, r0, n_pop, n_timesteps, n_sims) > .Machine$integer.max)) {
    stop("Initial conditions exceed the maximum integer size.")
  }
  print("Input checks passed")
}
