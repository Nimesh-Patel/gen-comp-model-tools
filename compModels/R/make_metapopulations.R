#' Set metapopulations in model
#'
#' This function may be deprecated, check input
#' output. Inputs instructions specifying metapopulations
#' in model. Overwrites any prior
#' instructions to add_metapopulation. Input can be
#' nested, although this functionality is
#' unnecessarilly complex and should be simplified.
#'
#' @param metapop_names vector of metapopulation names
#' default is ""
#' @param environment_names vector/list of environment names
#' default is ""
#' @param interactionscale vector same length as
#' input metapopulations/environment_names that
#' scales all interactions in given metapop/environment
#' default is 1
#' @param transitionscale vector same length as
#' input metapopulations/environment_names that
#' scales all transitions in given metapop/environment
#' default is 1
#' @param basestates specifies which states can
#' visit metapopulation/environment.
#' default is "" which specifies all basestates
#' @return list combining all information.
#' @export
#' @importFrom rlang .data
make_metapopulations <- function(
    metapop_names = "",
    environment_names = "", interactionscale = 1,
    transitionscale = 1, basestates = "") {
  variablenames <- c(
    "environment_names", "metapop_names",
    "interactionscale", "transitionscale"
  )
  variables <- lapply(variablenames, function(x) {
    eval(parse(text = x))
  })
  names(variables) <- variablenames
  parsevariablelength <- function(x) {
    eval(parse(text = paste0("length(unlist(", x, "))")))
  }
  variablenameslength <- sapply(variables, function(x) {
    length(unlist(x))
  })
  # choose (first) variable with most elements to match to
  matchthisvariablestring <- names(
    variablenameslength
  )[variablenameslength == max(
    variablenameslength
  )][1]
  matchthisvariable <- variables[[matchthisvariablestring]]
  # coerce all inputs to have same length as matched variable
  coercelengths <- sapply(matchthisvariable, length)
  total_elements <- sum(coercelengths)
  coercedinputs <- lapply(
    variables,
    function(x) {
      coerceinput2matchshape(x, matchthisvariable)
    }
  )

  # there is a nested structure, update metapop
  # names if it's a single value
  vnamelengthlogic <- variablenameslength[["metapop_names"]] == 1
  if ((max(coercelengths) > 1) && (vnamelengthlogic)) {
    warning("An input variable has a nested
    structure, but input metapopulation names
    contain only one element. Renaming metapopulations
    by numerical indices.")
    newmetapopnames <- as.character(seq_along(matchthisvariable))
    coercedinputs[["metapop_names"]] <- coerceinput2matchshape(
      newmetapopnames, matchthisvariable
    )
  }
  # deal with base states, make sure each element
  # base_state is a character vector and then
  # repeat appropriately
  coercedinputs[["basestates"]] <-
    coerceinput2matchshape_basestate(basestates, matchthisvariable)
  return(coercedinputs)
}
