#' Makes a tibble row for a given interaction
#'
#' Seems unnecessary, consider refactoring
#'
#' @param tblprocesses appears unnecessary
#' @param beforestates vector of state names
#' @param afterstates vector of statenames
#' @param rate numeric or character for
#' interaction rate
#' @param metapopulation metapopulation it occurs,
#' default is "" which specifies all metapopulations
#' @param environment environments it occurs,
#' default is "" which specifies all environments
#' @param name name of process
#' default is ""
#' @param namegroup type of process, rename
#' default is ""
#' @param beforestategroup specifies which state groups this
#' interaction applies to
#' default is "" which specifies all groups
#' @param afterstategroup specifies which state groups this
#' interaction applies to
#' default is "" which specifies all groups
#' @return tibble row with input values
process_row <- function(
    tblprocesses, beforestates, afterstates,
    rate, metapopulation = "", environment = "", name = "",
    namegroup = "", beforestategroup = "", afterstategroup = "") {
  add2tbl <- tibble::tibble(
    beforestates = beforestates,
    afterstates = afterstates,
    rate = rate,
    metapopulation = metapopulation,
    environment = environment,
    name = name,
    namegroup = namegroup,
    beforestategroup = beforestategroup,
    afterstategroup = afterstategroup
  )
  return(add2tbl)
}
