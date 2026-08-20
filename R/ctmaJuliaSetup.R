#' Prepare Julia for the optional CoTiMA estimation backend
#'
#' @description
#' Locates a Julia executable, installs the Julia-side packages needed by
#' JuliaCall, instantiates the Julia project bundled with CoTiMA, and
#' precompiles that project. All work is performed in separate Julia
#' processes; this function does not initialize JuliaCall in the current R
#' session. Restart R after setup before requesting a different Julia thread
#' count or after a failed JuliaCall initialization. On Linux, Julia setup
#' subprocesses do not inherit \code{LD_LIBRARY_PATH} or \code{LD_PRELOAD};
#' this prevents R, Java, or RStudio Server libraries from shadowing Julia's
#' own native dependencies.
#'
#' @param julia Optional path to the Julia executable. When omitted, the
#'   function searches PATH, common Juliaup launcher locations, and
#'   Juliaup's versioned Linux installation directories.
#' @param julia_home Optional directory containing the Julia executable. It is
#'   used only when `julia` is omitted.
#' @param project Optional path to the Julia project. The default is the
#'   `julia` directory installed with CoTiMA.
#' @param bootstrap_juliacall JuliaCall dependency policy. The default
#'   `"auto"` installs `Suppressor` or `RCall` only when missing and builds
#'   RCall only when newly installed. `TRUE` also forces an RCall rebuild;
#'   `FALSE` skips this bootstrap. The function deliberately does not broadly
#'   precompile the default Julia environment because that can trigger
#'   unrelated packages.
#' @param instantiate Instantiate the CoTiMA Julia project without automatic
#'   precompilation.
#' @param precompile Precompile the CoTiMA Julia project.
#' @param precompile_tasks Positive number of concurrent Julia precompilation
#'   tasks. The conservative default of one avoids excessive worker creation
#'   on high-core-count machines.
#' @param r_home R installation directory used to build Julia's RCall package.
#' @param dry_run Return the commands without executing them.
#' @param verbose Print commands and subprocess output.
#'
#' @return Invisibly, a list containing the Julia executable, Julia home,
#'   Julia version, project path, commands, and execution results.
#' @export
#'
#' @examples
#' \dontrun{
#' setup <- ctmaJuliaSetup()
#' setup$julia_home
#'
#' fit <- ctmaJuliaFit(
#'   model,
#'   estimation = "ML",
#'   julia_home = setup$julia_home
#' )
#' }
ctmaJuliaSetup <- function(
    julia = NULL,
    julia_home = NULL,
    project = NULL,
    bootstrap_juliacall = "auto",
    instantiate = TRUE,
    precompile = TRUE,
    precompile_tasks = 1L,
    r_home = R.home(),
    dry_run = FALSE,
    verbose = TRUE) {
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaJuliaSetup")
  on.exit(.ctma_timer(), add = TRUE)

  scalar_flag <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  if (is.character(bootstrap_juliacall) &&
      length(bootstrap_juliacall) == 1L &&
      identical(tolower(bootstrap_juliacall), "auto")) {
    bootstrap_mode <- "auto"
  } else if (is.logical(bootstrap_juliacall) &&
             length(bootstrap_juliacall) == 1L &&
             !is.na(bootstrap_juliacall)) {
    bootstrap_mode <- if (bootstrap_juliacall) "force" else "skip"
  } else stop(
    "'bootstrap_juliacall' must be \"auto\", TRUE, or FALSE.",
    call. = FALSE
  )
  instantiate <- scalar_flag(instantiate, "instantiate")
  precompile <- scalar_flag(precompile, "precompile")
  dry_run <- scalar_flag(dry_run, "dry_run")
  verbose <- scalar_flag(verbose, "verbose")

  precompile_tasks <- suppressWarnings(as.integer(precompile_tasks)[1L])
  if (is.na(precompile_tasks) || precompile_tasks < 1L) stop(
    "'precompile_tasks' must be a positive integer.", call. = FALSE
  )

  julia <- .ctmaJuliaFindExecutable(julia = julia, julia_home = julia_home)
  julia_home <- dirname(julia)

  if (is.null(project)) {
    project <- system.file("julia", package = "CoTiMA2")
  }
  if (!is.character(project) || length(project) != 1L || !nzchar(project) ||
      !dir.exists(project)) stop(
    "The CoTiMA Julia project was not found. Reinstall CoTiMA with its ",
    "inst/julia directory, or supply 'project'.", call. = FALSE
  )
  project <- normalizePath(project, mustWork = TRUE)
  project_file <- file.path(project, "Project.toml")
  backend_file <- file.path(project, "CtsemJulia.jl")
  if (!file.exists(project_file) || !file.exists(backend_file)) stop(
    "'project' must contain Project.toml and CtsemJulia.jl.", call. = FALSE
  )
  if (!is.character(r_home) || length(r_home) != 1L ||
      !dir.exists(r_home)) stop("'r_home' is not a directory.", call. = FALSE)
  r_home <- normalizePath(r_home, mustWork = TRUE)

  commands <- list()
  add_command <- function(label, args, env = character()) {
    commands[[length(commands) + 1L]] <<- list(
      label = label, executable = julia, args = args, env = env
    )
  }

  if (bootstrap_mode != "skip") {
    bootstrap_code <- if (bootstrap_mode == "auto") paste(
      "import Pkg;",
      "required = [\"Suppressor\", \"RCall\"];",
      "missing = filter(package -> Base.find_package(package) === nothing, required);",
      "isempty(missing) || Pkg.add(missing);",
      "\"RCall\" in missing && Pkg.build(\"RCall\")"
    ) else paste(
      "import Pkg;",
      "required = [\"Suppressor\", \"RCall\"];",
      "missing = filter(package -> Base.find_package(package) === nothing, required);",
      "isempty(missing) || Pkg.add(missing);",
      "Pkg.build(\"RCall\")"
    )
    add_command(
      "Bootstrap JuliaCall dependencies",
      c("--startup-file=no", "-e", shQuote(bootstrap_code)),
      c(
        paste0("R_HOME=", shQuote(r_home)),
        paste0("JULIA_NUM_PRECOMPILE_TASKS=", precompile_tasks),
        "JULIA_PKG_PRECOMPILE_AUTO=0"
      )
    )
  }

  project_argument <- paste0("--project=", shQuote(project))
  if (instantiate) {
    instantiate_code <- paste(
      "import Pkg;",
      "Pkg.instantiate(; allow_autoprecomp=false);",
      "Pkg.status()"
    )
    add_command(
      "Instantiate the CoTiMA Julia project",
      c("--startup-file=no", project_argument,
        "-e", shQuote(instantiate_code)),
      "JULIA_PKG_PRECOMPILE_AUTO=0"
    )
  }
  if (precompile) {
    add_command(
      "Precompile the CoTiMA Julia project",
      c("--startup-file=no", project_argument,
        "-e", shQuote("import Pkg; Pkg.precompile()")),
      paste0("JULIA_NUM_PRECOMPILE_TASKS=", precompile_tasks)
    )
  }

  render_command <- function(command) {
    command <- .ctmaJuliaCommandInvocation(command)
    paste(
    c(command$env, shQuote(command$executable), command$args),
    collapse = " "
  )}
  if (dry_run) {
    result <- structure(list(
      julia = julia,
      julia_home = julia_home,
      julia_version = NA_character_,
      project = project,
      commands = lapply(commands, render_command),
      results = list(),
      restart_required = FALSE
    ), class = "ctmaJuliaSetup")
    if (verbose) cat(paste(unlist(result$commands), collapse = "\n"), "\n")
    return(invisible(result))
  }

  version_command <- .ctmaJuliaCommandInvocation(list(
    executable = julia, args = "--version", env = character()
  ))
  version_output <- system2(
    version_command$executable, version_command$args,
    env = version_command$env, stdout = TRUE, stderr = TRUE
  )
  version_status <- attr(version_output, "status")
  if (!is.null(version_status) && version_status != 0L) stop(
    "Julia was found but could not be executed:\n",
    paste(version_output, collapse = "\n"), call. = FALSE
  )
  julia_version <- paste(version_output, collapse = " ")
  if (verbose) cat("Using ", julia_version, " at ", julia, "\n", sep = "")

  results <- vector("list", length(commands))
  for (index in seq_along(commands)) {
    command <- .ctmaJuliaCommandInvocation(commands[[index]])
    if (verbose) cat("\n", command$label, "...\n", sep = "")
    output <- system2(
      command$executable,
      command$args,
      env = command$env,
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    results[[index]] <- list(
      label = command$label, status = status, output = output
    )
    if (verbose && length(output)) cat(paste(output, collapse = "\n"), "\n")
    if (status != 0L) stop(
      command$label, " failed with status ", status, ":\n",
      paste(utils::tail(output, 30L), collapse = "\n"),
      "\nRestart R before retrying if JuliaCall initialization has failed.",
      call. = FALSE
    )
  }

  result <- structure(list(
    julia = julia,
    julia_home = julia_home,
    julia_version = julia_version,
    project = project,
    commands = lapply(commands, render_command),
    results = results,
      restart_required = bootstrap_mode != "skip"
  ), class = "ctmaJuliaSetup")
  if (verbose) {
    cat("\nJulia setup completed successfully.\n")
    if (bootstrap_mode != "skip") cat(
      "Restart R before initializing JuliaCall, then use julia_home = ",
      shQuote(julia_home), ".\n", sep = ""
    )
  }
  invisible(result)
}

# Build a subprocess invocation without executing it. RStudio Server and some
# Linux R installations prepend R/Java library directories to LD_LIBRARY_PATH.
# Julia can then load ABI-incompatible libcurl, libstdc++, or related native
# libraries and terminate with status 139. `env -u` genuinely removes inherited
# variables; system2(env=) can add or replace variables but cannot unset them.
.ctmaJuliaCommandInvocation <- function(command, linux = identical(
    unname(Sys.info()[["sysname"]]), "Linux")) {
  if (!isTRUE(linux)) return(command)
  env_executable <- if (file.exists("/usr/bin/env")) "/usr/bin/env" else
    unname(Sys.which("env"))
  if (!nzchar(env_executable)) stop(
    "Linux Julia setup requires the 'env' executable to sanitize native ",
    "library paths.", call. = FALSE
  )
  list(
    label = command$label,
    executable = env_executable,
    args = c(
      "-u", "LD_LIBRARY_PATH", "-u", "LD_PRELOAD",
      command$env, shQuote(command$executable), command$args
    ),
    env = character()
  )
}

.ctmaJuliaFindExecutable <- function(julia = NULL, julia_home = NULL) {
  executable_name <- if (.Platform$OS.type == "windows") "julia.exe" else "julia"

  if (!is.null(julia)) {
    if (!is.character(julia) || length(julia) != 1L || !nzchar(julia)) stop(
      "'julia' must be one executable path.", call. = FALSE
    )
    candidates <- julia
  } else if (!is.null(julia_home)) {
    if (!is.character(julia_home) || length(julia_home) != 1L ||
        !nzchar(julia_home)) stop(
      "'julia_home' must be one directory path.", call. = FALSE
    )
    candidates <- file.path(julia_home, executable_name)
  } else {
    path_candidate <- unname(Sys.which(executable_name))
    candidates <- c(
      path_candidate[nzchar(path_candidate)],
      path.expand(file.path("~/.juliaup/bin", executable_name)),
      path.expand(file.path("~/.julia", "juliaup", "bin", executable_name)),
      Sys.glob(path.expand(
        file.path("~/.julia", "juliaup", "julia-*", "bin", executable_name)
      ))
    )
  }
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (!length(candidates)) stop(
    "Julia was not found. Install Julia with Juliaup from ",
    "https://julialang.org/install/, restart RStudio, or supply 'julia'.",
    call. = FALSE
  )
  normalizePath(candidates[[length(candidates)]], mustWork = TRUE)
}
