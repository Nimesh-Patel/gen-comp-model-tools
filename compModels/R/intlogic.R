#' Determine if vector has integer values
#'
#' @param vec vector to check values
#' @return logical for containing integers
intlogic <- function(vec) {
  outlogic <- ((!(FALSE %in% (vec %% 1 == 0))) && (sum(vec) > 1))
  outlogic
}
