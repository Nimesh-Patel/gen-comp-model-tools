#' Plots deterministic compartmental models
#'
#' Plots deterministic compartmental model output in using [ggplot2::ggplot()]
#'
#' @param output data frame output from a compartmental model
#' @param stratify_by optional argument where user can name stratifications to
#' include. Each name will result in one plot per strata
#' @param show_intervention_period optional logical argument (default FALSE) to
#' add line(s) showing the start and/or end of the intervention
#' @param intervention_period optional argument (default NULL) to add a vector
#' (type: numeric) with the start and end time of an intervention
#' @param colors optional argument taking a vector of colors to use in plotting,
#' default is RColorBrewer::brewer.pal(8, "Dark2") and colors are used in
#' [ggplot2::scale_color_manual()] values argument
#' @return ggplot2 plot of compartmental model output
#' @import ggplot2
#' @importFrom rlang .data
#' @rdname plot_determin_model
#' @export
#' @examples
#' \dontrun{
#' modelout <- run_sir(
#'   init = c(s = 1e05 - 1, i = 1, r = 0),
#'   time = seq(0.1, 100, by = 0.1),
#'   parms = c(beta = 0.5, gamma = 0.1)
#' )
#' plot_determin_model(modelout)
#'
#' modelout2 <- run_sir_intervention(
#'   init = c(s = 1e05 - 1, i = 1, r = 0),
#'   time = seq(0.1, 100, by = 0.1),
#'   parms = list(
#'     beta = 0.5, gamma = 0.1,
#'     intervention_start_time = 20,
#'     intervention_end_time = 30,
#'     intervention_impact = 0.3
#'   )
#' )
#' plot_determin_model(modelout2,
#'   show_intervention_period = TRUE,
#'   intervention_period = c(20, 30)
#' )
#'
#' init_vals_base <- c(s = 99000, e = 0, i = 10000, r = 0)
#' comp_names_base <- c("s", "e", "i", "r")
#' beta <- .9 # transmission
#' sigma <- .1 # progression E->I
#' gamma <- 0.05 # recovery
#' trans_matrix_base <- matrix(
#'   c(
#'     -beta, 0, 0, 0,
#'     beta, -sigma, 0, 0,
#'     0, sigma, -gamma, 0,
#'     0, 0, gamma, 0
#'   ),
#'   nrow = length(comp_names_base), byrow = TRUE,
#'   dimnames = list(comp_names_base, comp_names_base)
#' )
#' times <- seq(0, 100, by = 1)
#' age_grps_vector <- c("child", "adult", "elderly")
#' vacc_status_vector <- c("vaccinated", "unvaccinated")
#' gender_cat_vector <- c("male", "female")
#'
#' modelout3 <- run_gen_deterministic(
#'   init_vals_base = init_vals_base,
#'   times = times,
#'   comp_names_base = comp_names_base,
#'   trans_matrix_base = trans_matrix_base,
#'   subgroups_list =
#'     list(
#'       age_groups = age_grps_vector,
#'       vaccination_statuses = vacc_status_vector,
#'       gender_categories = gender_cat_vector
#'     ),
#'   intervention_start_time = NULL,
#'   intervention_end_time = NULL,
#'   modifier_function = NULL
#' )
#' plot_determin_model(modelout3,
#'   stratify_by = c("child", "adult", "elderly")
#' )
#' }
plot_determin_model <- function(output, stratify_by = NULL,
                                show_intervention_period = FALSE,
                                intervention_period = NULL,
                                colors = RColorBrewer::brewer.pal(8, "Dark2")) {
  out_long <- as.data.frame(output) |>
    tidyr::pivot_longer(-"time", names_to = "variable", values_to = "value")
  plot_list <- list()

  unique_variables <- unique(out_long$variable)
  num_colors <- length(unique_variables)
  colors <- colors[1:num_colors]
  print(paste0(
    "There are ", num_colors, " unique compartments called ",
    paste(unique_variables, collapse = ", "),
    ", which are given these colors:",
    paste(colors, collapse = ", ")
  ))

  if (!is.null(stratify_by) && length(stratify_by) > 0) {
    for (stratum in stratify_by) {
      pattern <- paste0("(^|_)", stratum, "($|_)")
      filtered_data <- subset(out_long, grepl(pattern, out_long$variable))
      # it is 'variable' in the line above that causes a NOTE in check()
      p <- filtered_data |>
        ggplot(aes(x = .data$time, y = .data$value, color = .data$variable)) +
        geom_line() +
        scale_color_manual(values = colors) +
        xlab("Time") +
        ylab("Number") +
        labs(colour = "Compartment") +
        theme_classic()

      if (show_intervention_period && !is.null(intervention_period)) {
        p <- p +
          geom_vline(
            xintercept = intervention_period[1],
            linetype = "longdash",
            color = "blue"
          ) +
          geom_vline(
            xintercept = intervention_period[2],
            linetype = "longdash",
            color = "red"
          )
      }
      plot_list[[stratum]] <- p
    }
    do.call(gridExtra::grid.arrange, c(plot_list, ncol = 1))
  } else {
    p <- out_long |>
      ggplot(aes(x = .data$time, y = .data$value, color = .data$variable)) +
      geom_line() +
      scale_color_manual(values = colors) +
      xlab("Time") +
      ylab("Number") +
      labs(colour = "Compartment") +
      theme_classic()
    if (show_intervention_period && !is.null(intervention_period)) {
      p <- p + geom_vline(
        xintercept = intervention_period[1],
        linetype = "longdash",
        color = "blue"
      ) +
        geom_vline(
          xintercept = intervention_period[2],
          linetype = "longdash", color = "red"
        )
    }
    return(p)
  }
}
