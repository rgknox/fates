#!/bin/bash

# Path to FATES src

FC='gfortran'
F_OPTS="-shared -fPIC -O3 -llapack -ffree-line-length-none"
#F_OPTS="-shared -fPIC -O0 -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"

#FC='ifort'
#F_OPTS="-shared -fPIC"

MOD_FLAG="-J./bld"
INC_FLAG="-I./bld"

#MOD_FLAG="-I"

rm -f bld/*.o
rm -f bld/*.mod

# Build dgesv from lapack

${FC} ${F_OPTS} ${INC_FLAG} ${MOD_FLAG} -o bld/dgesvMod.o f90_src/lapack-3.10/dgesvMod.f

# Build the new file with constants

${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesConstantsMod.o ../../main/FatesConstantsMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/WrapShrMod.o f90_src/WrapShrMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/TwoStreamMLPEMod.o  ../../radiation/TwoStreamMLPEMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesGlobals.o ../../main/FatesGlobals.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/JSONParameterUtilsMod.o ../../main/JSONParameterUtilsMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesUtilsMod.o ../../main/FatesUtilsMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesParametersInterface.o ../../main/FatesParametersInterface.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/LeafBiophysicsMod.o ../../biogeophys/LeafBiophysicsMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/LeafBiophysSuppMod.o  ../leaf_biophys/f90_src/LeafBiophysSuppMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesRadiationMemMod.o ../../radiation/FatesRadiationMemMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/RadiationWrapMod.o ../radiation/f90_src/RadiationWrapMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesLitterMod.o  ../../biogeochem/FatesLitterMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/PRTParametersMod.o ../../parteh/PRTParametersMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/PRTGenericMod.o ../../parteh/PRTGenericMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesInterfaceTypesMod.o ../../main/FatesInterfaceTypesMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/EDParamsMod.o ../../main/EDParamsMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/EDPftvarcon.o ../../main/EDPftvarcon.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesFuelClassesMod.o ../../fire/FatesFuelClassesMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/SFParamsMod.o ../../fire/SFParamsMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesParameterDerivedMod.o ../../main/FatesParameterDerivedMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/EDParamsDerivedSupp.o f90_src/EDParamsDerivedSuppMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/AllometrySuppMod.o ../allometry/f90_src/AllometrySuppMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/DamageMainMod.o ../../biogeochem/DamageMainMod.F90
${FC} ${F_OPTS} $INC_FLAG ${MOD_FLAG} -o bld/FatesAllometryMod.o ../../biogeochem/FatesAllometryMod.F90



