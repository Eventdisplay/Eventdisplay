# Eventdisplay - An Analysis and Reconstruction Package for Ground-based Gamma-ray Astronomy

[![DOI](https://zenodo.org/badge/221222023.svg)](https://zenodo.org/badge/latestdoi/221222023)

Original developers: Gernot Maier and Jamie Holder

## Overview

Eventdisplay is a reconstruction and analysis pipline for data of
Imaging Atmospheric Cherenkov Telescopes (IACT).
It has been primarily developed for VERITAS and CTA analysis.

In case Eventdisplay is used in a research project, please cite the 
following publication:

Maier, G.; Holder, J., Eventdisplay: An Analysis and Reconstruction Package for 
Ground-based Gamma-ray Astronomy,  35th International Cosmic Ray Conference.
10-20 July, 2017. Bexco, Busan, Korea, Proceedings of Science, Vol. 301.
Online at [https://pos.sissa.it/cgi-bin/reader/conf.cgi?confid=301], id.747
[https://arxiv.org/abs/1708.04048]

For guidelines on installation, see INSTALL. For further information, 
see files in README directory

The package consists of several analysis steps and tools:

1. evndisp (calibrate and parametrize images, event reconstruction, stereo analysis)
2. mscw_energy (use lookup tables to produce msw, msl, and energies)
3. anasum (produce maps and calculate analysis results)
4. shared library tools and macros (see EVNDISP/lib/libVAnaSum.so and EVNDISP/macros/)  
   (produce the energy spectrum and integral fluxes, plot maps, etc.) 
5. makeEffectiveArea (calculate effective areas)
6. makeOptimizeBoxCutsTMVA (tools to optimize cuts)
7. ...

## Documentation

- INSTALL: information on installation the analysis package, dependencies, environmental variables
- README.CTA: description of a typical CTA analysis
- README.VERITAS.quick_summary: description of a typical VERITAS analysis
- AUTHORS: author description

Description and command line options for the different software parts:

- README.EVNDISP
- README.EVNDISP.commandline
- README.MSCW_ENERGY
- README.ANASUM
- README.EFFECTIVEAREA
- README.ANALYSISLIBRARY
- README.SCRIPTS
- README.MACROS

## Licence

License: BSD-3 (see LICENCE file)

## Contact information:

Gernot Maier (DESY)

## The Eventdisplay Ecosystem

Reconstruction and analysis can be complex; it requires inputs from different sources and execution of several indedependent stages.
Eventdisplay is in use since roughly 2004 and an ecosystem of libaries and repositories grew around the core code base. 
Below an overview of those repositories. 
Some are internal to VERITAS and not accessible to the general public.
There are differences between the Eventdisplay code applied in VERITAS and that made available in the
public repository (see below).

For almost every use case, Eventdisplay consists of at least three major components: 
- the code (Eventdisplay), 
- a library of scripts,
- a set of auxiliary files.

Care should be taken in using the correct versions (releases, tags, branches) combining these three types of repositories.
A blending of different versions of components will lead to incorrect results.

### Code, tools, library

The core library consist of all code, tools, and libraries required to run the analysis.
This includes the following main reconstruction tools (incomplete list):

1. `evndisp`: calibrate and parametrize images, event reconstruction, stereo analysis
2. `trainTMVAforAngularReconstruction`: train boosted decision trees for direction and energy reconstruction
3. `mscw_energy`: fill and use lookup tables for mean scaled with and lenght calculation, energy reconstruction, stereo reconstruction
4. `trainTMVAforGammaHadronSeparation`: train boosted decision trees for gamma/hadron separation
5. `makeEffectiveArea`: calculation of the instrument response functions (effective areas, angular point-spread function, energy resolution)
6. `makeRadialAcceptance`: calculation of radial camera acceptance from data files
7. `anasum`: analysis to calculate sky maps and spectral energy distribution
8. `libVAnaSum`: shared library tools (to be used with [ROOT](https://root.cern/) to e.g., plot instrument response function, spectral energy distributions, light curves, sky maps

Public code repository: <https://github.com/Eventdisplay/Eventdisplay>

Private (VERITAS) repository (to be used for VERITAS analysis): <https://github.com/VERITAS-Observatory/EventDisplay_v4>

### Analysis scripts

Typical use cases for Eventdisplay require the processing of many files (tens to several 100,000 in the case of CTA).
A library of scripts for the efficient execution is available and recommended to be used as the main access to the tools described in the section above.

The analysis scripts depend on the specific observatory use cases. 
The following repositories of scripts are available:

- for **CTA**: <https://github.com/Eventdisplay/Eventdisplay_AnalysisScripts_CTA>
- for **VERITAS** (private repository): <https://github.com/VERITAS-Observatory/Eventdisplay_AnalysisScripts_VTS> (to be used for version v485 and new; scripts are included in the [code library](https://github.com/VERITAS-Observatory/EventDisplay_v4) for older versions)
- for the **CTA pSCT** (incomple): <https://github.com/Eventdisplay/Eventdisplay_AnalysisScripts_pSCT>

### Auxiliary files for parameters, definitions, calibration values

The reconstruction of analysis requires information on the instrument (e.g., telescope positions), access information to data bases, parameters for the analysis (e.g., image cleaning parameters or instruction for the gamma/hadron separation), or basic calibration values.

This information is accessible through repositories for auxiliary files, again dependent on the observatory of interest:

- for **CTA**: <https://github.com/Eventdisplay/Eventdisplay_AnalysisFiles_CTA>
- for **VERITAS** (private repository): <https://github.com/VERITAS-Observatory/Eventdisplay_AnalysisFiles_VTS>

### Release tests

Reconstruction and analysis are complex and a series of tests are required before the release of a new version of Eventdisplay.

Release tests depend on the observatory of interest:

- for **VERITAS** (private repository): <https://github.com/VERITAS-Observatory/EventDisplay_ReleaseTests>

### Container applications

Docker files and images are provided for some Eventdisplay use cases.
Docker files are collected in <https://github.com/Eventdisplay/Eventdisplay_Docker>

### Converters

Tools to convert event lists into DL2 or DL3 format are collected in the converter repositories:

- for **CTA** (DL2): <https://github.com/Eventdisplay/Converters>
- for **VERITAS** (DL3, private repository): <https://github.com/VERITAS-Observatory/V2DL3>



