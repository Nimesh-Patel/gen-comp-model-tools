#' Plots stochastic compartmental models
#'
#' Plots stochastic compartmental model output using [ggplot2::ggplot()]
#'
#' @param output data frame output from a compartmental model
#' @param compartments character string specifying which compartment to plot,
#' defaults to "i" for infected
#' @param time_var character string specifying the name of the time variable,
#' default is NULL, internal logic checks vector of possible time variable names
#' (t, time)
#' @param show_intervention_period optional logical argument (default FALSE) to
#' add line(s) showing the start and/or end of the intervention
#' @param intervention_period optional argument (default NULL) to add a vector
#' (type: numeric) with the start and end time of an intervention
#' @return ggplot2 object of compartmental model output
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' modelout <- run_sir_stochastic(
#'   0.001, 0.1, 499, 100, 0, 500, 100,
#'   c("beta*s*i", "gamma*i"),
#'   matrix(c(-1, 0, +1, -1, 0, +1),
#'     nrow = 3, byrow = TRUE
#'   ), 10
#' )
#' plot_stoch_model(modelout)
#'
#' modelout2 <- run_sir_stochastic_tau(
#'   0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
#'   list(c(s = -1, i = 1), c(i = -1, r = 1)),
#'   10
#' )
#' plot_stoch_model(modelout2)
#' }
plot_stoch_model <- function(output,
                             compartments = c("i"),
                             time_var = NULL,
                             show_intervention_period = FALSE,
                             intervention_period = NULL) {
  possible_time_vars <- c("t", "time")

  if (is.null(time_var)) {
    found_time_var <- FALSE
    for (possible_time_var in possible_time_vars) {
      if (possible_time_var %in% colnames(output[[1]])) {
        time_var <- possible_time_var
        found_time_var <- TRUE
        print(paste0("Found time variable:", possible_time_var))
        break
      }
    }
    if (!found_time_var) {
      stop("Time variable not found in output data.")
    }
  }

  combined_sims_df <- dplyr::bind_rows(lapply(seq_along(output), function(i) {
    cbind(simulation = i, output[[i]])
  }), .id = "simulation_id")

  p <- ggplot() +
    theme_classic()

  for (compartment in compartments) {
    p <- p + geom_line(
      data = combined_sims_df,
      aes_string(
        x = time_var,
        y = compartment,
        group = "simulation_id",
        color = factor(compartment)
      ),
      alpha = 0.5
    )
  }

  if (show_intervention_period == TRUE && !is.null(intervention_period)) {
    p <- p + geom_vline(
      xintercept = intervention_period[1],
      linetype = "longdash"
    ) + geom_vline(
      xintercept = intervention_period[2],
      linetype = "longdash"
    )
  }

  # Add labels and colors
  p <- p + labs(
    title = "Stochastic SIR Model Simulations",
    x = "Time",
    y = "Number of Individuals"
  ) + scale_color_discrete(name = "Compartment")

  return(p)
}
