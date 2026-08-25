# Vignette sources and generated guides

The canonical documentation sources are the five `.Rnw` files in
`vignettes/`:

- `CoTiMA_User_Guide_1.0.4.Rnw`
- `CoTiMA_Julia_Usage.Rnw`
- `CoTiMA_Random_Effects.Rnw`
- `CoTiMA_RI_Random_Effects.Rnw`
- `Julia_Install.Rnw`

The nine PNG files in `vignettes/figures/` are source assets used by the
general user guide.

The correspondingly named PDF files in `inst/doc/` are generated outputs.
They are retained deliberately so that the guides are included in source
tarballs built with `--no-build-vignettes`. They must not be edited directly.
Regenerate them from the `.Rnw` sources with `tools/build-CoTiMA2.R` whenever
the published PDF content needs updating.
