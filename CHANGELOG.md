# Changelog

All notable changes to the Eventdisplay_v4 project will be documented in this file.
Changes for upcoming releases can be found in the [docs/changes](docs/changes) directory.
Note that changes before release v492.0 are not documented here, but can be found in the
[GitHub repository](https://github.com/VERITAS-Observatory/VTS-SimPipe/releases).

This changelog is generated using [Towncrier](https://towncrier.readthedocs.io/).

<!-- towncrier release notes start -->

## [v5.18.1](https://github.com/VERITAS-Observatory/EventDisplay_v4/releases/tag/v5.18.1) - 2026-02-11

### Bug Fixes

- Fix prod6 run script to select hyper-array layout applying case-insensitive file name checks.

### New Feature

## [v5.18.0](https://github.com/VERITAS-Observatory/EventDisplay_v4/releases/tag/v5.18.0) - 2026-02-10

### New Feature

- Improved prod6 container building to take different HESSIO flags for North / South grid productions into account. ([#188](https://github.com/VERITAS-Observatory/EventDisplay_v4/issues/188))


## [v5.17.0](https://github.com/VERITAS-Observatory/EventDisplay_v4/releases/tag/v5.17.0) - 2026-01-29

### New Feature

- Add support for reading XGB stereo reconstruction results (as friend trees). Effective area / IRF production can use now XGB stereo direction and energy. ([#186](https://github.com/VERITAS-Observatory/EventDisplay_v4/issues/186))


## [v5.16.0](https://github.com/VERITAS-Observatory/EventDisplay_v4/releases/tag/v5.16.0) - 2025-08-28

### New Feature

- Improvements to small-array analysis:

  - acceptance of loss values > 0.2 and introduction of mitigating steps
  - dispBDTEnergy trains on logE instead of `TMath::Power( 10., fdisp_energy_T[i] * log10( img_size[i] ) );`
  - energy reconstruction (dispBDTEnergy) uses median instead of weighted mean
  - introduce `TRAIN_RECO_QUALITY` and `TRAIN_RECO_METHOD` BDT machine learners
  - improved weight expression for training of dispBDTs
  - (many typo and code layout fixes)

  ([#171](https://github.com/VERITAS-Observatory/EventDisplay_v4/issues/171))
- Improve dispBDT analysis for small arrays.

  Re-calculate stereo parameters for the array configuration under consideration. This was actually a bug,
  as the stereo parameters calculated for the full (hyper) array were used to calculate the emission height
  and the direction coordinates used for the `cross` calculation. Recalculate using the simple stereo
  reconstructor.

  Change scale to logarithmic for several variables used by the dispBDTs (width, tgrad_x, length, cross) to improve
  regression performance. ([#178](https://github.com/VERITAS-Observatory/EventDisplay_v4/issues/178))

### Maintenance

- Enable towncrier change logs. ([#176](https://github.com/VERITAS-Observatory/EventDisplay_v4/issues/176))
- Fix GitHub build workflows after updates of several dependencies. ([#181](https://github.com/VERITAS-Observatory/EventDisplay_v4/issues/181))
