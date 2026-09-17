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
archive <- file.path(
  package_dir,
  sprintf("%s_%s.tar.gz", package_name, package_version)
)
r_bin <- file.path(R.home("bin"), "R")

vignette_dir <- file.path(package_dir, "vignettes")
expected_user_guide <- sprintf(
  "CoTiMA2_User_Guide_%s.Rnw",
  package_version
)
expected_user_guide_path <- file.path(vignette_dir, expected_user_guide)
if (!file.exists(expected_user_guide_path)) {
  available_guides <- list.files(
    vignette_dir,
    pattern = "^(CoTiMA2?|CoTiMA)_User_Guide_.*\\.[Rr]nw$"
  )
  stop(
    "Documentation version mismatch. Expected '", expected_user_guide,
    "' for package version ", package_version, ". Found: ",
    if (length(available_guides)) paste(available_guides, collapse = ", ")
    else "no versioned User Guide", "."
  )
}
user_guide_lines <- readLines(expected_user_guide_path, warn = FALSE)
required_guide_text <- c(
  "%\\VignetteIndexEntry{CoTiMA2 User Guide}",
  paste0("Version ", package_version)
)
missing_guide_text <- required_guide_text[
  !vapply(required_guide_text, function(x) {
    any(grepl(x, user_guide_lines, fixed = TRUE))
  }, logical(1L))
]
if (length(missing_guide_text)) {
  stop(
    "The User Guide does not match CoTiMA2 version ", package_version,
    ". Missing text: ", paste(missing_guide_text, collapse = "; "), "."
  )
}
legacy_sources <- list.files(vignette_dir, pattern = "^CoTiMA_.*\\.[Rr]nw$")
if (length(legacy_sources)) {
  stop(
    "Rename package documentation sources from 'CoTiMA_' to 'CoTiMA2_': ",
    paste(legacy_sources, collapse = ", "), "."
  )
}
vignette_sources <- list.files(vignette_dir, pattern = "\\.[Rr]nw$",
  full.names = TRUE)
if (length(vignette_sources)) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required to compile the Rnw guides.")
  }
  doc_dir <- file.path(package_dir, "inst", "doc")
  dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
  expected_user_guide_pdf <- sub("\\.[Rr]nw$", ".pdf", expected_user_guide)
  stale_user_guide_pdfs <- list.files(
    doc_dir,
    pattern = "^(CoTiMA2?|CoTiMA)_User_Guide_.*\\.pdf$",
    full.names = TRUE
  )
  stale_user_guide_pdfs <- stale_user_guide_pdfs[
    basename(stale_user_guide_pdfs) != expected_user_guide_pdf
  ]
  if (length(stale_user_guide_pdfs)) unlink(stale_user_guide_pdfs)
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

old_wd <- setwd(package_dir)
on.exit(setwd(old_wd), add = TRUE)
status <- system2(
  r_bin,
  c("CMD", "build", "--no-build-vignettes",
    "."),
  stdout = "",
  stderr = ""
)
if (!identical(status, 0L) || !file.exists(archive)) {
  stop("R CMD build failed.")
}

verify_lib <- tempfile("cotima2-install-")
dir.create(verify_lib)
status <- system2(
  r_bin,
  c(
    "CMD", "INSTALL", "--no-test-load",
    paste0("--library=", verify_lib),
    shQuote(basename(archive))
  ),
  stdout = "",
  stderr = ""
)
if (!identical(status, 0L)) {
  stop("The tarball was built, but clean-library installation failed. Check dependencies above.")
}

installed <- as.character(utils::packageVersion(package_name, lib.loc = verify_lib))
stopifnot(identical(installed, package_version))
message("Verified tarball: ", archive)
