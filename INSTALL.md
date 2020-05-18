#  INSTALLATION

## Prerequisites

- ROOT must be installed 
  version >= 6.14
  To compile root, use '-Dbuiltin_cfitsio=ON -Dbuiltin_gsl=ON'

- SOFA library [http://www.iausofa.org/current_C.html] must be installed. Use the script in the $EVNDISPSYS directory:
```
./install_sofa.sh
```
Set the following environmental variable:  SOFASYS=$EVNDISPSYS/sofa

### VERITAS analysis

(see VERITAS internal wiki for all details

### CTA analysis

(2.0) HESSIO libraries needed for the analysis of CTA Monte Carlo can be found here:
   [http://www.mpi-hd.mpg.de/hfm/CTA/internal/MC/Software/]

   (note: this side is password protected, usual CTA details)

### Optional

(3.0) for all FITS related output, cfitsio is needed (see [http://heasarc.gsfc.nasa.gov/fitsio/])

(3.1) GSL libraries (needed for FROGS image template method; ROOT included gsl is fine)
      [http://www.gnu.org/software/gsl/]

## Environmental Variables

### Compiling and linking

ROOTSYS :   (required) ROOT installation; add $ROOTSYS/lib to $LD_LIBRARY_PATH and $ROOTSYS/bin to $PATH 
            (root should be compiled with minuit2, mysql, xml)

SOFASYS:    (required) Astronomy library from Sofa (use ./install_sofa.sh) for installation)

HESSIOSYS : (optional) HESSIO libraries (for CTA analysis only); add $HESSIOSYS/lib to $LD_LIBRARY_PATH (not needed for VERITAS analysis)
            note: make sure to use the same compiler flags as used in the simulations, see Makefile)

VBSYS :     (optional) VBF libraries (for VERITAS analysis only); add $VBFSYS/bin to $PATH and $VBFSYS/lib to $LD_LIBRARY_PATH
            (VBF is not needed for CTA analysis)

FITSSYS :   (optional) FITS libraries (optional, not needed in most cases)

GSLSYS :    (optional) GSL libraries (needed for FROGS image template method)

### Analysis

EVNDISPSYS : EVNDISP directory (scripts expect binaries in $EVNDISPSYS/bin and libraries in $EVNDISPSYS/lib) 

For root versions >=6.xx: add $EVNDISPSYS/obj to LD_LIBRARY_PATH

### Data directories

(different auxiliary data files are needed for the analysis (e.g. calibration files, detector geometry, etc))

The following setup are for an environment where several users for example are analysis the same raw data
(would be in $VERITAS_DATA_DIR, or $CTA_DATA_DIR)), but having their analysis results written to their own
directories ($VERITAS_USER_DATA_DIR). 
Note that e.g. $VERITAS_DATA_DIR and $VERITAS_USER_DATA_DIR can point to the same directory.

For CTA analysis, the following directories are needed:

CTA_EVNDISP_AUX_DIR:      directory with all auxiliary data like calibration files, lookup tables, effective areas, etc
CTA_DATA_DIR :            global data directory: containing raw data or input simulation files
CTA_USER_DATA_DIR :       user data directory: containing output files from this analysis package (most of them in root format)

For VERITAS analysis, the following directories are needed:

VERITAS_EVNDISP_AUX_DIR:  directory with all auxiliary data like calibration files, lookup tables, effective areas, etc
VERITAS_DATA_DIR :        directory containing the raw telescope data or input simulation files 
VERITAS_USER_DATA_DIR :   user data directory: containing output files from this analysis package (most of them in root format)
VERITAS_IRFPRODUCTION_DIR : directory used in IRF production (not needed for normal data analysis; ignore it if you do not know it)

To switch settings for the different observatory (in the $EVNDISPSYS directory): 
```
./setObservatory.sh CTA
```

or
```
./setObservatory.sh VERITAS
```

# Compiling

Check your system's configuration:
```
make config
```
Compiling and installing:

in $EVNDISPSYS:

for VERITAS analysis type:
```
make VTS
```

or, for CTA analysis type:
```
make CTA
```

or, for CTA and VERITAS analysis type:
```
make all
```

## Makefile targets

- all	make all executables and libraries
- clean
- install
- CTA	make all CTA relevant binaries/libraries
- VTS	make all VERITAS relevant binaries/libraries
