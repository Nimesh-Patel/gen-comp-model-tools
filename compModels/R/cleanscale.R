#' Simplify auto-generated math expressions
#'
#' This function appends "*" to rates and converts "1*"
#' and "*" to empty strings "".
#' This facilitates and simplifies specifying rates
#' auto-generated during model compilation
#'
#'
#' @param currtbl A dataframe/tibble with columns
#' "interactionscale" and "transitionscale"
#' @return tibble with updated interactionscale
#' and transitionscale values
#' @importFrom rlang .data
cleanscale <- function(currtbl) {
  tblout <- currtbl |>
    dplyr::mutate(
      interactionscale = paste0(as.character(.data$interactionscale), "*"),
      transitionscale = paste0(as.character(.data$transitionscale), "*")
    ) |>
    dplyr::mutate(
      interactionscale =
        dplyr::case_when(.data$interactionscale == "*" ~ "",
          .data$interactionscale == "1*" ~ "",
          .default = .data$interactionscale
        )
    ) |>
    dplyr::mutate(
      transitionscale =
        dplyr::case_when(.data$transitionscale == "*" ~ "",
          .data$transitionscale == "1*" ~ "",
          .default = .data$transitionscale
        )
    )
  return(tblout)
}
