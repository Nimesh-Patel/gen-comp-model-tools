#' Calculates initial growth rate of outbreak
#'
#' Helper function to calculate the initial growth rate for the outbreak from
#' user supplied compartmental model output by regressing the log of the
#' cumulative incidence on time to estimate the rate of exponential increase r.
#' This r can then be used to estimate R0 with the relationship R0 = V*r + 1
#' where V is the serial interval.
#'
#' @param output data frame output from a compartmental model
#' @param time_var character string specifying the name of the time variable,
#' defaults to "time"
#' @param i_var character string specifying the name of the infected state
#' variable of interest, defaults to "i"
#' @return calculated initial growth rate (numeric)
#' @export
#' @examples
#' \dontrun{
#' # Deterministic Model
#' out <- run_sir(
#'   init = c(s = 1e05 - 1, i = 1, r = 0),
#'   time = seq(0.1, 100, by = 0.1),
#'   parms = c(beta = 0.00001, gamma = 0.1)
#' )
#' calculate_initial_growth_rate(out)
#'
#' # Stochastic Model with 10 simulations
#' init_vals <- c(s = 990, e = 100, i1 = 10, i2 = 0, r = 0)
#' change_matrix <- matrix(
#'   c(
#'     -1, 0, 0, 0, # S -> E
#'     +1, -1, 0, 0, # E -> I1
#'     0, +1, -1, 0, # I1 -> I2
#'     0, 0, +1, -1, # I2 -> R
#'     0, 0, 0, +1
#'   ),
#'   nrow = length(init_vals),
#'   byrow = TRUE,
#'   dimnames = list(names(init_vals), NULL)
#' )
#' modelout <- run_gen_stochastic(
#'   parms_vec = c(
#'     beta = 0.001, sigma = 0.1,
#'     alpha = 0.05, gamma = 0.1
#'   ),
#'   propensity_fns = c(
#'     "beta * s * i1",
#'     "sigma * e",
#'     "alpha * i1",
#'     "gamma * i2"
#'   ),
#'   init_vals = c(
#'     s = 990, e = 100,
#'     i1 = 10, i2 = 0, r = 0
#'   ),
#'   n_timesteps = 100,
#'   change_matrix = change_matrix,
#'   n_sims = 10
#' )
#' calculate_initial_growth_rate(modelout, time_var = "t", i_var = "e")
#' }
calculate_initial_growth_rate <- function(output, time_var = "time",
                                          i_var = "i") {
  if (is.data.frame(output)) {
    output$cum_inf <- cumsum(output[[i_var]])
    formula_string <- paste("log(cum_inf) ~", time_var)
    fit <- lm(formula_string, data = output)
    r <- fit$coef[[time_var]]
    r
  } else if (is.list(output)) {
    mod_output <- calculate_cuminf(output, i_var = i_var)
    lapply(output,
      calculate_initial_growth_rate, # recursion
      time_var = time_var,
      i_var = i_var
    )
  }
}
