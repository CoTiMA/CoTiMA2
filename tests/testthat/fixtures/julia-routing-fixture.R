ctma_julia_routing_fixture <- function() {
  setup <- structure(
    list(
      marker = "small stored Julia routing fixture",
      argumentList = list(
        optimize = TRUE,
        fit = FALSE,
        scaleTime = 1,
        scaleTI = FALSE
      )
    ),
    class = "CoTiMAFit"
  )
  fit <- structure(
    list(
      source_object = setup,
      optimization = list(loglik = -123.5),
      likelihood = list(julia = -130),
      coefficients = list(raw = c(0.1, -0.2, 0.3))
    ),
    class = "ctmaJuliaFit"
  )
  list(setup = setup, fit = fit)
}
