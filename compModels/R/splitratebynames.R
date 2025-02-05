#' Replaces substrings in string
#'
#' Useful when updating the basestate in a
#' rate with the actual state name
#'
#' @param currrate String to replace substrings
#' @param oldnames substrings to replace (commonly basestates)
#' @param newnames Replacing substrings in same order
#' as oldnames (commonly updated state names)
#' @return Full string with replaced substrings
splitratebynames <- function(currrate, oldnames, newnames) {
  currrate_tokened <- tokenizerate(currrate)
  splittokenlogic <- currrate_tokened %in% oldnames

  currpaste0 <- ""
  for (currlogicidx in seq(splittokenlogic)) {
    if (splittokenlogic[currlogicidx]) {
      curroldnamelogic <- currrate_tokened[currlogicidx] == oldnames
      currpaste0 <- paste0(currpaste0, newnames[curroldnamelogic, ])
    } else {
      currpaste0 <- paste0(currpaste0, currrate_tokened[currlogicidx])
    }
  }
  return(currpaste0)
}
