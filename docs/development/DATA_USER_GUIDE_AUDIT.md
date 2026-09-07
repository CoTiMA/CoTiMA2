# Data objects and the general CoTiMA user guide

Audit date: 2026-08-25

## Scope and method

All `.rda` files in `data/` were loaded separately into temporary
environments. Before the resulting cleanup, the directory contained 126
files and exactly one named object per file.

Object names were matched case-sensitively and with identifier boundaries
against `vignettes/CoTiMA2_User_Guide_0.14.14.Rnw`. This prevents, for example,
`CoTiMAFullFit_6` from being counted merely because
`CoTiMAFullFit_6_new` occurs in the guide.

## Result

The initial audit found two objects not named in the general user guide:

- `alphas313`
- `pubList_8`

## Resolution

- The `ctmaCorRel()` example now uses the already documented `empcov128` and
  `alphas128` objects. The redundant `empcov313` and `alphas313` package-data
  objects and their data documentation were removed.
- The Appendix B description of `ctmaGetPub()` now refers readers to
  `summary(pubList_8)`.

After this cleanup, `data/` contains 124 files with one object per file, and
all 124 objects are named in `CoTiMA2_User_Guide_0.14.14.Rnw`.
