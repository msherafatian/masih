#' Inverted versions of in, is.null and is.na
#'
#' @noRd
#'
#' @examples
#' 1 %not_in% 1:10
#' not_null(NULL)
`%not_in%` <- Negate(`%in%`)

not_null <- Negate(is.null)

not_na <- Negate(is.na)

#' Remove NULL elements from list
#'
#' @param .x A list
#' @noRd
drop_nulls <- function(.x) {
  .x[!sapply(.x, is.null)]
}

#' Check if all elements are NULL
#'
#' @param ... Elements to check
#' @noRd
all_null <- function(...) {
  all(sapply(list(...), is.null))
}

#' Check if any elements are NULL  
#'
#' @param ... Elements to check
#' @noRd
any_null <- function(...) {
  any(sapply(list(...), is.null))
}

#' Safe data loading with error handling
#'
#' @param file_path Path to file
#' @param type Type of file (rds, csv, etc)
#' @noRd
safe_load_data <- function(file_path, type = "rds") {
  tryCatch({
    switch(type,
           "rds" = readRDS(file_path),
           "csv" = read.csv(file_path),
           stop("Unsupported file type")
    )
  }, error = function(e) {
    message("Error loading file: ", e$message)
    NULL
  })
}