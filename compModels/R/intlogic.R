#' Determine if popsize has integer values
#'
#' @param tblpopsize table with popsize column
#' @return logical for containing integers
intlogic <- function(tblpopsize) {
  vec <- tblpopsize[["popsize"]]
  logic_decimal <- !(FALSE %in% (vec %% 1 == 0))

  if ("metapopulation" %in% colnames(tblpopsize)) {
    tblpopsize <- tblpopsize |> dplyr::group_by("metapopulation")
  }

  maxsizeinmetapop <-
    tblpopsize |>
    dplyr::summarize(sumvals = sum(.data$popsize)) |>
    dplyr::pull("sumvals") |>
    max()
  logic_max <- maxsizeinmetapop > 1
  outlogic <- logic_max && logic_decimal

  outlogic
}
