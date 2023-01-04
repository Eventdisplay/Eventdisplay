# Eventdisplay - An Analysis and Reconstruction Package for Ground-based Gamma-ray Astronomy


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6814319.svg)](https://doi.org/10.5281/zenodo.6814319)
[![ascl:2212.002](https://img.shields.io/badge/ascl-2212.002-blue.svg?colorB=262255)](https://ascl.net/2212.002)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)
[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/6753/badge)](https://bestpractices.coreinfrastructure.org/projects/6753)

[![CI](https://github.com/Eventdisplay/Eventdisplay/actions/workflows/ci.yml/badge.svg)](https://github.com/Eventdisplay/Eventdisplay/actions/workflows/ci.yml)
[![CTA-Prod5 Docker Image](https://github.com/Eventdisplay/Eventdisplay/actions/workflows/packages-cta-prod5.yml/badge.svg)](https://github.com/Eventdisplay/Eventdisplay/actions/workflows/packages-cta-prod5.yml)
[![CTA-slib Docker Image](https://github.com/Eventdisplay/Eventdisplay/actions/workflows/packages-cta-slib.yml/badge.svg)](https://github.com/Eventdisplay/Eventdisplay/actions/workflows/packages-cta-slib.yml)


* Authors and contributors: [CITATION.cff](CITATION.cff)
* Licence: [LICENSE](LICENSE)


## Overview

Eventdisplay is a reconstruction and analysis pipline for data of
Imaging Atmospheric Cherenkov Telescopes (IACT).
It has been primarily developed for VERITAS and CTA analysis and used in
many VERITAS and CTA publications. 
This repository contains the Eventdisplay version used for CTA and other analyses.
For the VERITAS version, please go to [https://github.com/VERITAS-Observatory/EventDisplay_v4](https://github.com/VERITAS-Observatory/EventDisplay_v4)

Original developers: Gernot Maier and Jamie Holder

In case Eventdisplay is used in a research project, please cite this repository and
the following publication:

Maier, G.; Holder, J., Eventdisplay: An Analysis and Reconstruction Package for 
Ground-based Gamma-ray Astronomy,  35th International Cosmic Ray Conference.
10-20 July, 2017. Bexco, Busan, Korea, Proceedings of Science, Vol. 301.
Online at [https://pos.sissa.it/cgi-bin/reader/conf.cgi?confid=301], id.747
[https://arxiv.org/abs/1708.04048]

For guidelines on installation, see [INSTALL.md](INSTALL.md). For further information, 
see files in [README](./README) directory

The package consists of several analysis steps and tools:

1. `evndisp`: calibrate and parametrize images, event reconstruction, stereo analysis
2. `trainTMVAforAngularReconstruction`: train boosted decision trees for direction and energy reconstruction
3. `mscw_energy`: fill and use lookup tables for mean scaled with and lenght calculation, energy reconstruction, stereo reconstruction
4. `trainTMVAforGammaHadronSeparation`: train boosted decision trees for gamma/hadron separation
5. `makeEffectiveArea`: calculation of the instrument response functions (effective areas, angular point-spread function, energy resolution)
6. `makeRadialAcceptance`: calculation of radial camera acceptance from data files
7. `anasum`: analysis to calculate sky maps and spectral energy distribution
8. `libVAnaSum`: shared library tools (to be used with [ROOT](https://root.cern/) to e.g., plot instrument response function, spectral energy distributions, light curves, sky maps

## Documentation

- [INSTALL.md](INSTALL.md): information on installation the analysis package, dependencies, environmental variables
- [README.CTA](README/README.CTA): description of a typical CTA analysis
- README.VERITAS.quick_summary: description of a typical VERITAS analysis
- AUTHORS: author description

Description and command line options for the different software parts:

- [README.EVNDISP](README/README.EVNDISP)
- [README.MSCW_ENERGY](README/README.MSCW_ENERGY)
- [README.EFFECTIVEAREA](README/README.EFFECTIVEAREA)
- [README.ANALYSISLIBRARY](README/README.ANALYSISLIBRARY)

## The Eventdisplay Ecosystem

Reconstruction and analysis can be complex; it requires inputs from different sources and execution of several indedependent stages.
Eventdisplay is in use since roughly 2004 and an ecosystem of libaries and repositories grew around the core code base. 
Below an overview of those repositories. 
Some are internal to VERITAS and not accessible to the general public.

For almost every use case, Eventdisplay consists of at least three major components: 
- the code (Eventdisplay), 
- a library of scripts,
- a set of auxiliary files.

Care should be taken in using the correct versions (releases, tags, branches) combining these three types of repositories.
A blending of different versions of components will lead to incorrect results.

## Docker images

Docker images are made available for the following use cases.

- CTA prod5 analysis: analysis of prod5 CTA simulations; [Dockerfile](dockerfiles/Dockerfile-cta-prod5); docker image available from [container page](https://github.com/Eventdisplay/Eventdisplay/pkgs/container/eventdisplay) with tag `cta-prod5`
- CTA slib: analysis library used for CTA; [Dockerfile](dockerfiles/Dockerfile-cta-slib); docker image available from [container page](https://github.com/Eventdisplay/Eventdisplay/pkgs/container/eventdisplay) with tag `cta-slib`

See [dockerfiles/README.md](dockerfiles/README.md) on usage.

(for the time being, there is some overlap with the [Eventdisplay container](https://github.com/Eventdisplay/Eventdisplay_Docker) repository.

### Code, tools, library

The core library consist of all code, tools, and libraries required to run the analysis.

### Analysis scripts

Typical use cases for Eventdisplay require the processing of many files (tens to several 100,000 in the case of CTA).
A library of scripts for the efficient execution is available and recommended to be used as the main access to the tools described in the section above.

The analysis scripts depend on the specific observatory use cases. 
The following repositories of scripts are available:

- for **CTA**: [https://github.com/Eventdisplay/Eventdisplay_AnalysisScripts_CTA](https://github.com/Eventdisplay/Eventdisplay_AnalysisScripts_CTA)
- for the **CTA pSCT** (incomplete): [https://github.com/Eventdisplay/Eventdisplay_AnalysisScripts_pSCT](https://github.com/Eventdisplay/Eventdisplay_AnalysisScripts_pSCT)

### Auxiliary files for parameters, definitions, calibration values

The reconstruction of analysis requires information on the instrument (e.g., telescope positions), access information to data bases, parameters for the analysis (e.g., image cleaning parameters or instruction for the gamma/hadron separation), or basic calibration values.

This information is accessible through repositories for auxiliary files, again dependent on the observatory of interest:

- for **CTA**: [https://github.com/Eventdisplay/Eventdisplay_AnalysisFiles_CTA](https://github.com/Eventdisplay/Eventdisplay_AnalysisFiles_CTA)

### Container applications

Docker files and images are provided for some Eventdisplay use cases.
Docker files are collected in [https://github.com/Eventdisplay/Eventdisplay_Docker](https://github.com/Eventdisplay/Eventdisplay_Docker)

### Converters

Tools to convert event lists into DL2 or DL3 format are collected in the converter repositories:

- for **CTA** (DL2): [https://github.com/Eventdisplay/Converters](https://github.com/Eventdisplay/Converters)
- for **VERITAS** (converter to GADF DL3 Fromat): [https://github.com/VERITAS-Observatory/V2DL3](https://github.com/VERITAS-Observatory/V2DL3)

## Contact

For any questions, contact Gernot Maier (gernot.maier@desy.de)
