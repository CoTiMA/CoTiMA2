.ctmaPubSurnameKey <- function(name) {
  name <- trimws(gsub("^&\\s*", "", as.character(name)))
  surname <- vapply(strsplit(name, ";", fixed = TRUE), function(parts) {
    trimws(parts[length(parts)])
  }, character(1L))
  has_comma <- grepl(",", surname, fixed = TRUE)
  surname[has_comma] <- trimws(vapply(
    strsplit(surname[has_comma], ",", fixed = TRUE), `[`, character(1L), 1L))
  surname <- iconv(surname, from = "", to = "ASCII//TRANSLIT")
  tolower(gsub("[^[:alnum:]]", "", surname))
}

.ctmaPubStudyAuthors <- function(source) {
  source <- as.character(source)
  year_positions <- which(grepl("^[0-9]{4}$", trimws(source)))
  if (!length(year_positions)) {
    stop("The study source contains no four-digit publication year.",
      call. = FALSE)
  }
  year_position <- year_positions[length(year_positions)]
  authors <- if (year_position > 1L) source[seq_len(year_position - 1L)] else
    character()
  authors <- authors[!grepl("^study\\s*[0-9]*$", trimws(authors),
    ignore.case = TRUE)]
  list(authors = authors, year = as.integer(source[year_position]))
}

.ctmaPubAuthorCounts <- function(frequencies, target_year, recency,
    years_to_exclude, aggregate_function) {
  years <- suppressWarnings(as.integer(names(frequencies)))
  frequencies <- as.numeric(frequencies)
  valid <- is.finite(frequencies) & !is.na(years) &
    !(years %in% years_to_exclude)
  frequencies <- frequencies[valid]
  years <- years[valid]
  before <- years < target_year
  recent <- before & years >= target_year - recency
  list(
    total = if (any(before)) aggregate_function(
      frequencies[before], na.rm = TRUE) else NA_real_,
    recent = if (any(recent)) aggregate_function(
      frequencies[recent], na.rm = TRUE) else NA_real_
  )
}
