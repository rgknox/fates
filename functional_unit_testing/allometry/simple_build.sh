#!/bin/bash

FC='gfortran -g -shared -fPIC'

# First copy over the FatesConstants file, but change the types of the fates_r8 and fates_int
# --------------------------------------------------------------------------------------------

old_fates_r8_str=`grep -e integer ../../main/FatesConstantsMod.F90 | grep fates_r8 | sed 's/^[ \t]*//;s/[ \t]*$//'`
new_fates_r8_str='use iso_c_binding, only: fates_r8  => c_double'

old_fates_int_str=`grep -e integer ../../main/FatesConstantsMod.F90 | grep fates_int | sed 's/^[ \t]*//;s/[ \t]*$//'`
new_fates_int_str='use iso_c_binding, only: fates_int => c_int'

# Add the new lines (need position change, don't swap)

sed "/implicit none/i $new_fates_r8_str" ../../main/FatesConstantsMod.F90 > f90src/FatesConstantsMod.F90
sed -i "/implicit none/i $new_fates_int_str" f90src/FatesConstantsMod.F90

# Delete the old lines

sed -i "/$old_fates_r8_str/d" f90src/FatesConstantsMod.F90
sed -i "/$old_fates_int_str/d" f90src/FatesConstantsMod.F90

sed -i "/private/d" f90src/FatesConstantsMod.F90

# Now lest copy over EDTypesMod.F90's dinc_vai and nlevleaf to AllometryMod.F90 so we
# dont need to copy over the whole module

old_edtypes_str=`cat ../../biogeochem/FatesAllometryMod.F90 | grep EDTypes` #  | sed 's/^[ \t]*//;s/[ \t]*$//'`
#nlevleaf_str=`cat ../../main/EDTypesMod.F90 | grep nlevleaf | grep public`
#dinc_str=`cat ../../main/EDTypesMod.F90 | grep dinc_vai | grep public` # | sed 's/^[ \t]*//;s/[ \t]*$//'`

nlevleaf_str='integer, parameter :: nlevleaf = 1'
dinc_str='real(r8) :: dinc_vai(1) = (/1000/)'

sed "/implicit none/a $nlevleaf_str" ../../biogeochem/FatesAllometryMod.F90 > f90src/FatesAllometryMod.F90

sed -i "/implicit none/a $dinc_str" f90src/FatesAllometryMod.F90

sed -i "/$old_edtypes_str/d" f90src/FatesAllometryMod.F90
#sed -i '/use EDTypesMod , only : nlevleaf, dinc_vai/d' f90src/FatesAllometryMod.F90

#grep -v "$old_edtypes_str" f90src/FatesAllometryMod.F90 > temp && mv temp f90src/FatesAllometryMod.F90

# This re-writes the wrapper so that it uses all the correct parameters
# in FatesAllometryMod.F90
python AutoGenVarCon.py


# Procedure for auto-generating AllomUnitWrap
# 1) scan FatesAllometry and create list of EDPftVarcon_inst variables
# 2) scan EDpftVarcon and get the name of the in-file parameter names associated
#    with these variables




rm -f include/*.o
rm -f include/*.mod


# Build the new file with constants

${FC} -I include/ -J include/ -o include/FatesConstantsMod.o  f90src/FatesConstantsMod.F90

#${FC} -I include/ -J include/ -o include/PRTParametersMod.o  ../../parteh/PRTParametersMod.F90

${FC} -I include/ -J include/ -o include/UnitWrapMod.o f90src/UnitWrapMod.F90

${FC} -I include/ -J include/ -o include/FatesAllometryMod.o f90src/FatesAllometryMod.F90


#${FC} -g -o include/FatesConstantsMod.o  ../main/FatesConstantsMod.F90

#gfortran -shared -fPIC -g -o include/EDTypesMod.o ../main/EDTypesMod.F90




#gfortran
