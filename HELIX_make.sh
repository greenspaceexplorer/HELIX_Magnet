#!/bin/bash

cd code/python
rm helix_mag*
echo FILES IN PYTHON CODE DIRECTORY
ls

cd ../fortran
echo FILES IN FORTRAN CODE DIRECTORY
ls
make clean -f Makefile
cd src/
f2py -c -m helix_mag CoilDataMod.f90 RotationMod.f90 MagnetLib.f90
cp helix_mag* ../../python
