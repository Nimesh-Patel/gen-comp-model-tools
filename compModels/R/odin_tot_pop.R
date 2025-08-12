#' Create calculate total population for odin model
#'
#' Total value is called N and is the sum of provided states
#' Total N can then be used in model calculations
#'
#' @param model_pop model population information from
#' compiled model $modelinstructions$tblntotal
#' @return string formatted to calculate total population
#' @family calculations
#' @export
odin_tot_pop <- function(model_pop) {
  glue::glue("Ntotal <- {model_pop$totalpop}")
}
