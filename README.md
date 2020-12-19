# Eventdisplay - An Analysis and Reconstruction Package for Ground-based Gamma-ray Astronomy

Original developers: Gernot Maier and Jamie Holder

[![DOI](https://zenodo.org/badge/221222023.svg)](https://zenodo.org/badge/latestdoi/221222023)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c4203cc5d45f4db2b30affccd6e0c641)](https://www.codacy.com/manual/GernotMaier/Eventdisplay?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Eventdisplay/Eventdisplay&amp;utm_campaign=Badge_Grade)

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

