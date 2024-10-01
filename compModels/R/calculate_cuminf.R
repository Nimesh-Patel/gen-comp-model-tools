#' Calculates cumulative incidence on state variable
#'
#' Helper function to calculate the cumulative sum of a state value,
#' specifically intended for summing the cumulative number of infections at a
#' time step and giving the cumulative sum (cumulative incidence) in a new
#' data frame column called cum_inf.
#'
#' @param output data frame output from a compartmental model
#' @param i_var character string specifying the name of the state
#' variable of interest, defaults to "i"
#' @return data frame with additional column cum_inf that gives the cumulative
#' incidence at that time (numeric)
#' @importFrom stats lm
#' @export
calculate_cuminf <- function(output, i_var = "i") {
  lapply(
    output,
    function(df) {
      df$cum_inf <- cumsum(df[[i_var]])
      df
    }
  )
}
