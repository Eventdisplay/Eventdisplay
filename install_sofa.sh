#!/bin/bash
# 
# install the sofa package into the $EVNDISPSYS directory
#
# see https://www.iausofa.org for a description
#
#

echo "Installation of sofa into $EVNDISPSYS "

CURDIR=`pwd`
cd "$EVNDISPSYS"

echo "Checking for existing sofa installation " 

if [ -d "sofa" ] && [ -d "sofa/lib" ]
then
    echo "Error, sofa directory exists. Please remove."
    cd "$CURDIR"
    exit
fi

mkdir sofa
cd sofa

# get sofa package from the web page and install
SOFAD="20210512"
SOFA="sofa_c-${SOFAD}.tar.gz"
wget --no-check-certificate https://www.iausofa.org/2021_0512_C/${SOFA}
if [ ! -e ${SOFA} ]
then
    echo "error in downloading sofa package"
    exit
fi
tar -xvzf ${SOFA}
rm -f ${SOFA}

##########################
# prepare make file
cd sofa/${SOFAD}/c/src/
sed -i -- "s/\$(HOME)/\$(EVNDISPSYS)\/sofa/" makefile
# use clang on OSX
OS=`uname -s`
echo "$OS"
if [ "$OS" = "Darwin" ]
then
    sed -i -- 's/gcc/clang/' makefile
else
    sed -i -- 's/pedantic/pedantic -fPIC/' makefile
fi

# make and install
make
make install
make clean
cd ../../../../
rm -rf sofa

echo "Installation completed"
echo "Please set the following environmental variable: "
echo "export SOFASYS=$EVNDISPSYS/sofa"

cd "$CURDIR"
