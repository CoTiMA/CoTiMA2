# Vignette sources and generated guides

The canonical documentation sources are the six `.Rnw` files in
`vignettes/`:

- `CoTiMA2_User_Guide_0.14.14.Rnw`
- `CoTiMA2_Julia_Usage.Rnw`
- `CoTiMA2_Julia_Fit-Diagnostics.Rnw`
- `CoTiMA2_Random_Effects.Rnw`
- `CoTiMA2_RI_Random_Effects.Rnw`
- `Julia_Install.Rnw`

Package-document names use `CoTiMA2`; the scientific method may still be
called CoTiMA in explanatory prose. Before every new tarball, the version in
the User Guide filename, title page, and vignette metadata must match the
`Version` field in `DESCRIPTION`. `tools/build-CoTiMA2.R` enforces this rule and
stops before building when the documentation version is stale.

The nine PNG files in `vignettes/figures/` are source assets used by the
general user guide.

The correspondingly named PDF files in `inst/doc/` are generated outputs.
They are retained deliberately so that the guides are included in source
tarballs built with `--no-build-vignettes`. They must not be edited directly.
Regenerate them from the `.Rnw` sources with `tools/build-CoTiMA2.R` whenever
the published PDF content needs updating.
