#' ctmaSaveFile
#'
#' @description Internal fcuntion to save files
#'
#' @param activeDirectory directory name
#' @param SaveObject object to save
#' @param FileName filename
#' @param Directory directory to save file in
#' @param silentOverwrite Explicitly allow replacing an existing file. In
#'   noninteractive sessions, an existing file otherwise raises an error; no
#'   input is requested. Interactive users may overwrite, skip, rename or cancel.
#'
#' @importFrom crayon red blue
#'
#' @return No return value. Just saves files
#'
ctmaSaveFile <- function(
                           activeDirectory=activeDirectory,
                           SaveObject,
                           FileName,
                           Directory,
                           silentOverwrite=FALSE)
{
  if (!is.logical(silentOverwrite) || length(silentOverwrite) != 1L ||
      is.na(silentOverwrite)) stop("silentOverwrite must be TRUE or FALSE.")
  directory <- file.path(activeDirectory, Directory)
  target <- file.path(directory, FileName)
  if (file.exists(target) && !silentOverwrite) {
    if (!.ctmaSaveFileInteractive()) stop(
      "Output file already exists: ", target,
      ". Non-interactive use requires silentOverwrite = TRUE or a different filename.",
      call. = FALSE)

    repeat {
      choice <- tolower(trimws(.ctmaSaveFileRead(
        paste0("File exists: ", target, ". Overwrite [y], skip [n], change name [c], quit [q]? "))))
      if (!nzchar(choice) || choice == "q") stop("Saving cancelled.", call. = FALSE)
      if (choice == "n") return(invisible(NULL))
      if (choice == "y") break
      if (choice == "c") {
        name <- trimws(.ctmaSaveFileRead("New filename: "))
        if (!nzchar(name)) stop("Saving cancelled.", call. = FALSE)
        if (!grepl("[.]rds$", name, ignore.case = TRUE)) name <- paste0(name, ".rds")
        target <- file.path(directory, name)
        # A renamed target must pass the same overwrite check.
        if (!file.exists(target)) break
      } else message("Please answer y, n, c or q.")
    }
  }
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
    stop("Could not create output directory: ", directory, call. = FALSE)
  saveRDS(SaveObject, file = target)
  invisible(NULL)
}

.ctmaSaveFileInteractive <- function() interactive()
.ctmaSaveFileRead <- function(prompt) readline(prompt)
