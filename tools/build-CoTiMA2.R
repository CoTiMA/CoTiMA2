## Build CoTiMA2 and verify that the resulting source package installs.
## Run from the package root with: source("tools/build-CoTiMA2.R")

package_dir <- normalizePath(getwd(), mustWork = TRUE)
description <- file.path(package_dir, "DESCRIPTION")
if (!file.exists(description)) {
  stop("Run this script with the CoTiMA package directory as the working directory.")
}

desc <- read.dcf(description)
package_name <- unname(desc[1L, "Package"])
package_version <- unname(desc[1L, "Version"])
parent_dir <- dirname(package_dir)
archive <- file.path(parent_dir, sprintf("%s_%s.tar.gz", package_name, package_version))
r_bin <- file.path(R.home("bin"), "R")

vignette_dir <- file.path(package_dir, "vignettes")
vignette_sources <- list.files(vignette_dir, pattern = "\\.[Rr]nw$",
  full.names = TRUE)
if (length(vignette_sources)) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required to compile the Rnw guides.")
  }
  doc_dir <- file.path(package_dir, "inst", "doc")
  dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
  build_dir <- tempfile("cotima2-vignettes-")
  dir.create(build_dir)
  on.exit(unlink(build_dir, recursive = TRUE, force = TRUE), add = TRUE)
  inputs <- list.files(vignette_dir, full.names = TRUE, all.files = TRUE,
    no.. = TRUE)
  file.copy(inputs, build_dir, recursive = TRUE)
  vignette_wd <- setwd(build_dir)
  on.exit(setwd(vignette_wd), add = TRUE)
  for (source in basename(vignette_sources)) {
    tex <- sub("\\.[Rr]nw$", ".tex", source)
    pdf <- sub("\\.[Rr]nw$", ".pdf", source)
    knitr::knit(source, output = tex, quiet = TRUE)
    tools::texi2pdf(tex, clean = TRUE, quiet = TRUE)
    if (!file.exists(pdf)) stop("Vignette compilation failed for ", source)
    file.copy(pdf, file.path(doc_dir, pdf), overwrite = TRUE)
  }
  setwd(vignette_wd)
}

old_wd <- setwd(parent_dir)
on.exit(setwd(old_wd), add = TRUE)
status <- system2(
  r_bin,
  c("CMD", "build", "--no-build-vignettes",
    shQuote(basename(package_dir))),
  stdout = "",
  stderr = ""
)
if (!identical(status, 0L) || !file.exists(archive)) {
  stop("R CMD build failed.")
}

verify_lib <- tempfile("cotima2-install-")
dir.create(verify_lib)
on.exit(unlink(verify_lib, recursive = TRUE, force = TRUE), add = TRUE)
status <- system2(
  r_bin,
  c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(verify_lib),
    shQuote(archive)),
  stdout = "",
  stderr = ""
)
if (!identical(status, 0L)) {
  stop("The tarball was built, but clean-library installation failed. Check dependencies above.")
}

installed <- as.character(utils::packageVersion(package_name, lib.loc = verify_lib))
stopifnot(identical(installed, package_version))
message("Verified tarball: ", archive)
