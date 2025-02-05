#' Generate deSolve function
#'
#' @param expressions expressions of instantaneous
#' population changes
#' @return deSolve input function
odegenerator <- function(expressions) {
  outode <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      evaluated_values <- sapply(
        expressions,
        function(expr) eval(parse(text = expr))
      )
      list(evaluated_values)
    }) # end with(as.list ...
  }
  return(outode)
}
