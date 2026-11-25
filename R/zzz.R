.ctma_package_name <- "CoTiMA2"
.ctma_package_version <- "0.14.69.9003"
.ctma_call_state <- new.env(parent = emptyenv())
.ctma_call_state$depth <- 0L

.ctmaVersionNote <- function() {
  if (.ctma_call_state$depth == 0L) {
    .ctma_print_message(paste0(
      "This is the experimental CoTiMA 2.0 package (version ",
      utils::packageVersion(.ctma_package_name), ")."
    ))
  }
  invisible(NULL)
}

.ctmaRegisterExecutionTimer <- function(function_name) {
  started <- proc.time()[["elapsed"]]
  .ctma_call_state$depth <- .ctma_call_state$depth + 1L
  finished <- FALSE
  force(function_name)
  function() {
    if (!finished) {
      .ctma_call_state$depth <- max(0L, .ctma_call_state$depth - 1L)
      finished <<- TRUE
    }
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed < 60) {
      formatted <- sprintf("%.3f seconds", elapsed)
    } else if (elapsed < 3600) {
      formatted <- sprintf("%.2f minutes", elapsed / 60)
    } else {
      formatted <- sprintf("%.2f hours", elapsed / 3600)
    }
    message(
      function_name,
      " total wall-clock execution time: ",
      formatted,
      "."
    )
    invisible(elapsed)
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(.ctmaFormatNote(paste0(
    "This is the experimental CoTiMA 2.0 package (version ",
    utils::packageVersion(pkgname, lib.loc = libname), ")."
  )))
}
