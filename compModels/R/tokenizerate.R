#' Used to isolate user input state and parameter
#' names within mathematical functions
#'
#' @param ratestr Character string
#' @return Character vector where each element
#' are the substrings separated by common
#' mathematical operations.
tokenizerate <- function(ratestr) {
  tokensplits <- c("+", "-", "*", "/", "(", ")", "^")
  whichtoken <-
    do.call(
      rbind,
      stringr::str_locate_all(ratestr, stringr::fixed(tokensplits))
    )[, 1]
  ratestrchars <- unlist(strsplit(ratestr, ""))
  ratestrnumchars <- length(ratestrchars)
  whichtoken_wends <- c(0, whichtoken, ratestrnumchars)
  allstart <- sort(unique(c(1, whichtoken, whichtoken + 1)))
  allend <- sort(unique(c(whichtoken - 1, whichtoken, ratestrnumchars)))
  allstart <- allstart[(allstart > 0) & (allstart <= ratestrnumchars)]
  allend <- allend[(allend > 0) & (allend <= ratestrnumchars)]
  if (length(allstart) != length(allend)) {
    print("different lengths")
  }
  updatedrates <- c()
  for (i in seq_along(allstart)) { # linter changed
    startidx <- allstart[i]
    endidx <- allend[i]
    currtoken <- stringr::str_flatten(ratestrchars[startidx:endidx])
    updatedrates <- c(updatedrates, stringr::str_flatten(currtoken))
  }
  return(updatedrates)
}
