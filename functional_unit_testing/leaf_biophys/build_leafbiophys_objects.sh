#!/bin/bash

# Path to FATES src

FC='gfortran'

F_OPTS="-shared -fPIC -g -O0 -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"


MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod

# Build the new file with constants

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesConstantsMod.o ../../main/FatesConstantsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/WrapShrMod.o f90_src/WrapShrMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesUtilsMod.o ../../main/FatesUtilsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/LeafBiophysicsMod.o ../../biogeophys/LeafBiophysicsMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/LeafBiophysSuppMod.o  f90_src/LeafBiophysSuppMod.F90





