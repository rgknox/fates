#!/bin/bash

# Path to FATES src

FC='gfortran'

#F_OPTS="-fPIC -O3 -llapack"
F_OPTS="-g -fPIC"
F_OBJ_OPTS="-shared"

FATES_PATH='/home/rgknox/Models/CTSM/src/fates/'

#F_OPTS="-fPIC -O0 -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"

MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod
rm -f bld/*.a

# Build dgesv from lapack
${FC} ${F_OPTS} -shared -I bld/ -J./bld/ -o bld/FatesConstantsMod.so ${FATES_PATH}/main/FatesConstantsMod.F90
${FC} ${F_OPTS} -shared -I bld/ -J./bld/ -o bld/WrapShrMod.so ${FATES_PATH}/functional_unit_testing/leaf_biophys/f90_src/WrapShrMod.F90
${FC} ${F_OPTS} -shared -I bld/ -J./bld/ -o bld/FatesUtilsMod.so ${FATES_PATH}/main/FatesUtilsMod.F90
${FC} ${F_OPTS} -shared -I bld/ -J./bld/ -o bld/LeafBiophysicsMod.so ${FATES_PATH}/biogeophys/LeafBiophysicsMod.F90
${FC} ${F_OPTS} -shared -I bld/ -J./bld/ -o bld/LeafBiophysSuppMod.so ${FATES_PATH}/functional_unit_testing/leaf_biophys/f90_src/LeafBiophysSuppMod.F90
ar -rcs bld/photo_infer.a bld/FatesConstantsMod.so bld/WrapShrMod.so bld/FatesUtilsMod.so bld/LeafBiophysicsMod.so bld/LeafBiophysSuppMod.so 
${FC} ${F_OPTS} -o photo_infer PhotoNNInfer.F90  -L${FTorch_DIR}/build/ -lftorch -ltorch -ltorch_cpu -I${FTorch_DIR}/build/modules/ -L${Torch_DIR}/lib/ -L./bld -I bld/ -J./bld/ bld/photo_infer.a -Wl,-rpath,bld/
