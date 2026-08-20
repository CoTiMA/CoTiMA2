.ctma_package_name <- "CoTiMA2"
.ctma_package_version <- "0.9.91"

.ctmaVersionNote <- function() {
  message(
    "Note: This is the experimental CoTiMA 2.0 package (version ",
    .ctma_package_version,
    ")."
  )
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Experimental CoTiMA 2.0 package (CoTiMA2 version ",
    .ctma_package_version,
    "). The original CoTiMA package can remain installed separately."
  )
}
