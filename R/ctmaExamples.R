.ctma_example_registry <- new.env(parent = emptyenv())

.ctmaExampleDirectory <- function() {
  installed <- system.file("example-data", package = "CoTiMA2")
  if (nzchar(installed)) return(installed)
  source_directory <- file.path(getwd(), "inst", "example-data")
  if (dir.exists(source_directory)) return(source_directory)
  stop("The CoTiMA2 example-data directory could not be found.", call. = FALSE)
}

#' Load CoTiMA2 user-guide example objects explicitly
#'
#' @description Loads packaged user-guide objects into a chosen environment.
#' They are deliberately not attached with the package, preventing accidental
#' use by functions such as \code{ctmaPrep()}.
#'
#' @param names Character vector of object names. With \code{NULL}, all
#' packaged example objects are loaded.
#' @param envir Target environment; defaults to the global environment.
#' @param overwrite Whether existing objects may be overwritten. The default
#' is \code{FALSE}.
#' @return Invisibly, the names loaded.
#' @export
ctmaLoadExamples <- function(names = NULL, envir = .GlobalEnv,
    overwrite = FALSE) {
  .ctmaVersionNote()
  directory <- .ctmaExampleDirectory()
  files <- list.files(directory, pattern = "\\.[Rr][Dd][Aa]$", full.names = TRUE)
  available <- tools::file_path_sans_ext(basename(files))
  if (is.null(names)) names <- available
  names <- unique(as.character(names))
  missing <- setdiff(names, available)
  if (length(missing)) stop("Unknown CoTiMA2 example object(s): ",
    paste(missing, collapse = ", "), call. = FALSE)
  conflicts <- names[vapply(names, exists, logical(1L), envir = envir,
    inherits = FALSE)]
  if (length(conflicts) && !isTRUE(overwrite)) stop(
    "Object(s) already exist and were not overwritten: ",
    paste(conflicts, collapse = ", "), call. = FALSE)

  loaded <- character()
  for (name in names) {
    previous_exists <- exists(name, envir = envir, inherits = FALSE)
    previous_value <- if (previous_exists) {
      get(name, envir = envir, inherits = FALSE)
    } else NULL
    temporary <- new.env(parent = emptyenv())
    object_names <- load(files[match(name, available)], envir = temporary)
    if (!name %in% object_names) stop("Example file does not contain '", name,
      "'.", call. = FALSE)
    value <- get(name, envir = temporary, inherits = FALSE)
    assign(name, value, envir = envir)
    key <- paste0(format(envir), "\r", name)
    assign(key, list(
      loaded_value = value,
      previous_exists = previous_exists,
      previous_value = previous_value
    ), envir = .ctma_example_registry)
    loaded <- c(loaded, name)
  }
  message("Loaded ", length(loaded), " CoTiMA2 example object(s) into ",
    environmentName(envir), ".")
  invisible(loaded)
}

#' Unload unchanged CoTiMA2 user-guide example objects
#'
#' @param names Character vector of loaded example names. With \code{NULL},
#' all objects loaded by \code{ctmaLoadExamples()} into \code{envir} are
#' considered.
#' @param envir Environment from which unchanged objects are removed.
#' @return Invisibly, a list naming removed, modified, and missing objects.
#' @export
ctmaUnloadExamples <- function(names = NULL, envir = .GlobalEnv) {
  .ctmaVersionNote()
  prefix <- paste0(format(envir), "\r")
  tracked <- ls(.ctma_example_registry, all.names = TRUE)
  tracked <- tracked[startsWith(tracked, prefix)]
  tracked_names <- substring(tracked, nchar(prefix) + 1L)
  if (is.null(names)) names <- tracked_names
  names <- unique(as.character(names))
  removed <- modified <- absent <- untracked <- character()
  for (name in names) {
    key <- paste0(prefix, name)
    if (!exists(key, envir = .ctma_example_registry, inherits = FALSE)) {
      untracked <- c(untracked, name)
    } else if (!exists(name, envir = envir, inherits = FALSE)) {
      absent <- c(absent, name)
      rm(list = key, envir = .ctma_example_registry)
    } else {
      current <- get(name, envir = envir, inherits = FALSE)
      record <- get(key, envir = .ctma_example_registry, inherits = FALSE)
      if (identical(current, record$loaded_value)) {
        if (isTRUE(record$previous_exists)) {
          assign(name, record$previous_value, envir = envir)
        } else {
          rm(list = name, envir = envir)
        }
        rm(list = key, envir = .ctma_example_registry)
        removed <- c(removed, name)
      } else modified <- c(modified, name)
    }
  }
  if (length(modified)) message("Kept modified object(s): ",
    paste(modified, collapse = ", "), ".")
  message("Removed ", length(removed), " unchanged CoTiMA2 example object(s).")
  invisible(list(removed = removed, modified = modified, missing = absent,
    untracked = untracked))
}
