#!/bin/bash

cd code/fortran/src
f2py -c -m helix_mag CoilDataMod.f90 RotationMod.f90 MagnetLib.f90
