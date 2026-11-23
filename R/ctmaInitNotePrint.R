# Preserve every original message; only its default printing changes.
.ctmaInitNoteFormat <- function(notes) {
  raw <- as.character(notes)
  matched <- grepl("^Study [^:]+: ", raw)
  study <- sub("^Study ([^:]+): .*$", "\\1", raw[matched])
  text <- sub("^Study [^:]+: ", "", raw[matched])
  type <- rep("O", length(text))
  type[startsWith(text, "unavailable uncertainty for ")] <- "U"
  type[grepl("^[0-9]+ weak parameter direction", text)] <- "W"
  type[startsWith(text, "At most two observed waves per person")] <- "T"
  type[text == "saturated parameter transformation."] <- "S"
  type[text == "random-effect correlation at its boundary."] <- "B"
  type[startsWith(text, "unstable curvature-based interval widths for ")] <- "I"
  ids <- unique(study)
  numeric_ids <- suppressWarnings(as.numeric(ids))
  ids <- if (length(ids) && all(is.finite(numeric_ids))) ids[order(numeric_ids)] else sort(ids)
  table <- data.frame(Study = ids, Weak = rep(NA_character_, length(ids)),
    Unavailable = rep(NA_integer_, length(ids)), Flags = rep("-", length(ids)),
    stringsAsFactors = FALSE)
  for (i in seq_along(ids)) {
    at <- study == ids[i]
    w <- unique(sub("^([0-9]+) weak.*", "\\1", text[at & type == "W"]))
    if (length(w)) table$Weak[i] <- paste(w, collapse = "/")
    u <- text[at & type == "U"]
    if (length(u)) table$Unavailable[i] <- length(unique(trimws(unlist(strsplit(
      sub("^unavailable uncertainty for ", "", u), ",", fixed = TRUE)))))
    flags <- c("T", "S", "B", "I", "O")
    flags <- flags[flags %in% type[at]]
    if (length(flags)) table$Flags[i] <- paste(flags, collapse = ",")
  }
  counts <- vapply(c("U", "W", "T", "S", "B", "I", "O"), function(code)
    length(unique(study[type == code])), integer(1))
  list(table = table, counts = counts,
    wave_explanation = unique(text[type == "T"]),
    other = unique(c(raw[!matched], raw[matched][type == "O"])))
}

#' @export
print.ctmaIdentificationNotes <- function(x, ..., details = FALSE) {
  if (!is.logical(details) || length(details) != 1L || is.na(details))
    stop("details must be TRUE or FALSE.", call. = FALSE)
  if (details) {
    print(as.character(x), ...)
    return(invisible(x))
  }
  f <- .ctmaInitNoteFormat(x)
  cat("Identification and uncertainty notes: ", length(x), " entries for ",
    nrow(f$table), " studies.\n", sep = "")
  if (nrow(f$table)) print(f$table, row.names = FALSE, right = FALSE)
  width <- max(40L, min(120L, getOption("width", 80L)))
  say <- function(text) cat(paste(strwrap(text, width = width), collapse = "\n"), "\n", sep = "")
  studies <- function(n) paste(n, if (n == 1L) "study" else "studies")
  if (f$counts[["W"]]) say(paste0("Weak: number of weak parameter directions (", studies(f$counts[["W"]]),
    "). Optimizer convergence does not establish identification."))
  if (f$counts[["U"]]) say(paste0("Unavailable: number of parameter coordinates with unavailable uncertainty (",
    studies(f$counts[["U"]]), "), not the number of weak directions. Full parameter names are retained in the detailed notes."))
  if (f$counts[["T"]]) {
    say(paste0("[T] Two-wave design (", studies(f$counts[["T"]]), "):"))
    for (note in f$wave_explanation) say(note)
  }
  legend <- c(S = "Saturated parameter transformation.", B = "Random-effect correlation at its boundary.",
    I = "Unstable curvature-based interval widths.")
  for (code in names(legend)) if (f$counts[[code]]) say(paste0("[", code, "] ", legend[[code]],
    " Affected studies: ", f$counts[[code]], "."))
  if (length(f$other)) {
    say("[O] Other notes (shown in full):")
    for (note in f$other) say(note)
  }
  say("NA means not recorded in these notes; absence of a flag does not certify a clean fit. No study is excluded by this report.")
  say("Full messages: print(fit$summary$identification_notes, details = TRUE), or as.character(fit$summary$identification_notes).")
  invisible(x)
}

.ctmaInitWaveMessage <- function(note, first = TRUE) {
  if (first) paste0("[T] ", note) else
    "[T] At most two observed waves per person; see the full [T] explanation earlier in this ctmaInit() call."
}
