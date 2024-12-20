#!/bin/bash

# Path to FATES src

FC='gfortran'

F_OPTS="-shared -fPIC -O3 -llapack"

#F_OPTS="-shared -fPIC -O0 -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"

MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod

# Build dgesv from lapack

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/dgesvMod.o f90_src/lapack-3.10/dgesvMod.f


# Build the new file with constants

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesConstantsMod.o ../../main/FatesConstantsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/WrapShrMod.o ../leaf_biophys/f90_src/WrapShrMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesUtilsMod.o ../../main/FatesUtilsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/LeafBiophysicsMod.o ../../biogeophys/LeafBiophysicsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/LeafBiophysSuppMod.o  ../leaf_biophys/f90_src/LeafBiophysSuppMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesRadiationMemMod.o ../../radiation/FatesRadiationMemMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/TwoStreamMLPEMod.o  ../../radiation/TwoStreamMLPEMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/RadiationWrapMod.o ../radiation/f90_src/RadiationWrapMod.F90




