#' Calculates initial growth rate of outbreak
#'
#' Helper function used by [compModels::calculate_initial_growth_rate()] to
#' calculate the initial growth rate for the outbreak from
#' user supplied compartmental model output by regressing the log of the
#' incidence on time to estimate the rate of exponential increase r.
#' This r can then be used to estimate R0 with the relationship R0 = V*r + 1
#' where V is the serial interval.
#'
#' Currently this calculation is for an SIR model and uses the equation
#' exp(i0) * exp((beta - gamma) * t)
#'
#' @param output data frame output from a compartmental model
#' @param time_var character string specifying the name of the time variable,
#' defaults to "time"
#' @param i_var character string specifying the name of the infected state
#' variable of interest, defaults to "i"
#' @param beta numeric value for transmission rate, default is NULL
#' @param gamma numeric value for recovery rate, default is NULL
#' @param lower_time numeric value for time step on which will serve as lower
#' bound in initial growth rate calculation, default is NULL
#' @param upper_time numeric value for time step on which will serve as upper
#' bound in initial growth rate calculation, default is NULL
#' @return prints initial growth rate(s) to console (numeric) and returns one
#' or more plots of the simulated observations, calculated theoretical growth
#' rate in the user-specified time range, and observed growth rate based on the
#' calculated growth rate r.
#' @import ggplot2
#' @importFrom rlang .data :=
#' @examples
#' \dontrun{
#' # Deterministic Model
#' out <- run_sir(
#'   init = c(s = 1e05 - 1, i = 1, r = 0),
#'   time = seq(0.1, 100, by = 0.1),
#'   parms = c(beta = 0.5, gamma = 0.1)
#' )
#' calculate_initial_growth_rate(out,
#'   time_var = "time",
#'   beta = 0.5, gamma = 0.1,
#'   lower_time = 10,
#'   upper_time = 25
#' )
#'
#' # Stochastic Model with 10 simulations
#' init_vals <- c(s = 499, i = 1, r = 0)
#' change_matrix <- matrix(
#'   c(
#'     -1, 0, # S -> I #nolint
#'     +1, -1, # I -> R #nolint
#'     0, +1
#'   ),
#'   nrow = length(init_vals),
#'   byrow = TRUE,
#'   dimnames = list(names(init_vals), NULL)
#' )
#'
#' modelout <- run_gen_stochastic(
#'   parms_vec = c(
#'     beta = 0.5, gamma = 0.1
#'   ),
#'   propensity_fns = c(
#'     "(beta * s * i)/n",
#'     "gamma * i"
#'   ),
#'   init_vals = init_vals,
#'   n_timesteps = 100,
#'   change_matrix = change_matrix,
#'   n_sims = 10
#' )
#' plot_stoch_model(modelout, compartments = c("s", "i", "r"))
#' calculate_initial_growth_rate(modelout,
#'   time_var = "t", i_var = "i",
#'   beta = 0.5, gamma = 0.1,
#'   lower_time = 0,
#'   upper_time = 25
#' )
#' }
apply_growth_rate_calc <- function(output, time_var = NULL,
                                   i_var = "i",
                                   beta = NULL, gamma = NULL,
                                   lower_time = NULL,
                                   upper_time = NULL) {
  output$inc <- output[[i_var]]
  output$log_inc <- log(output[[i_var]])
  max_log_inc <- max(output$log_inc)
  names(output)[names(output) == time_var] <- "time"


  formula_string <- paste("log_inc ~ time")
  # This line below needed for stochastic applications
  output <- output |>
    dplyr::filter(
      !is.na(.data$log_inc) &
        !is.nan(.data$log_inc) &
        is.finite(.data$log_inc)
    )
  output$time <- as.numeric(output$time)
  output2 <- output |>
    dplyr::filter(dplyr::between(.data$time, lower_time, upper_time))

  fit <- lm(stats::as.formula(formula_string), data = output2)
  r <- fit$coef[["time"]]
  print(paste("Growth rate is:", r))

  i_0 <- output$log_inc[1]
  print(paste("log(i_0) is:", i_0))
  stoch_t <- output$time
  stoch_t_f <- stoch_t[stoch_t >= lower_time & stoch_t <= upper_time]

  theorygrowth <-
    function(beta, gamma, t, i0) {
      outdyn <- exp(i0) * exp((beta - gamma) * t)
    }
  predicted_inc <- theorygrowth(beta, gamma, stoch_t_f, i_0)

  r_growth <-
    function(r, t, i0) {
      outdyn <- exp(i0) * exp(r * t)
    }
  second_r_line <- r_growth(r, stoch_t_f, i_0)

  output3 <- output |>
    dplyr::select("time", "log_inc") |>
    dplyr::filter(dplyr::between(.data$time, lower_time, upper_time))

  print(length(stoch_t_f))
  print(length(predicted_inc))
  print(length(second_r_line))
  print(length(log(predicted_inc)))
  print(length(output3$log_inc))

  exp_growth_df <- data.frame(
    time = stoch_t_f,
    pred = predicted_inc,
    rline = log(second_r_line), # theory using calculated r
    log_pred_inc = log(predicted_inc), # theory
    log_inc = output3$log_inc # observations
  )
  utils::str(exp_growth_df)

  exp_growth_df2 <- exp_growth_df |>
    tidyr::pivot_longer(!"time",
      names_to = "group",
      values_to = "value"
    )
  utils::str(exp_growth_df2)


  p <- ggplot(
    data = subset(exp_growth_df2, !"group" %in% c("log_inc", "pred")),
    aes(x = .data$time, y = .data$value, color = .data$group)
  ) +
    geom_point(data = output, aes(
      x = .data$time, y = .data$log_inc,
      color = "log_inc"
    )) +
    geom_line() +
    labs(
      title = "Log-Transformed Incidence Over Time with
                  Exponential Growth Curve",
      y = "Log(Incidence)",
      color = "Legend"
    ) +
    scale_color_manual(
      labels = c(
        "Simulated Observations",
        "Theoretical Expected Growth",
        "Expected Growth (Calculated r)"
      ),
      values = c(
        log_inc = "grey",
        log_pred_inc = "red",
        rline = "blue"
      )
    ) +
    theme_bw() +
    theme(legend.position = "top") +
    scale_x_continuous(limits = c(
      min(output$time),
      max(output$time)
    ))
}
