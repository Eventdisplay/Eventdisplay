#!/bin/bash
# 
# install the sofa package into the $EVNDISPSYS directory
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
wget http://www.iausofa.org/2019_0722_C/sofa_c-20190722.tar.gz
if [ ! -e sofa_c-20190722.tar.gz ]
then
    echo "error in downloading sofa package"
    exit
fi
tar -xvzf sofa_c-20190722.tar.gz
rm -f sofa_c-20190722.tar.gz

##########################
# prepare make file
cd sofa/20190722/c/src/
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

echo "Installation completed"
echo "Please set the following environmental variable: "
echo "export SOFASYS=$EVNDISPSYS/sofa"

cd "$CURDIR"
