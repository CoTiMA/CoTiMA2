.ctmaCheckRetiredArguments <- function(call, function_name, allow_other = FALSE) {
  dots <- as.list(call[["..."]])
  if (!length(dots)) return(invisible(NULL))
  retired <- intersect(names(dots), c("activateRPB", "experimental"))
  if (length(retired)) {
    hints <- vapply(retired, function(name) {
      value <- dots[[name]]
      if (identical(value, FALSE)) {
        paste0("Remove ", name, " = FALSE.")
      } else if (identical(value, TRUE)) {
        paste0("Remove ", name, " = TRUE. Omitting it does not reproduce the former option.")
      } else {
        paste0("Remove ", name, " and check whether the script relied on the former option.")
      }
    }, character(1))
    stop(function_name, "() no longer supports ", paste(retired, collapse = ", "),
         ". ", paste(hints, collapse = " "),
         " Use explicit names for the remaining arguments.", call. = FALSE)
  }
  if (!allow_other)
    stop(function_name, "() received unused arguments: ",
         paste(vapply(dots, deparse1, character(1)), collapse = ", "),
         ". Check the argument names in help().", call. = FALSE)
  invisible(NULL)
}
