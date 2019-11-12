#!/bin/bash
# This "script" is just a wrapper that calls a C program
# Calling this script with no arguments will print a help string

$EVNDISPSYS/bin/VTS.getRunListFromDB $*

exit
