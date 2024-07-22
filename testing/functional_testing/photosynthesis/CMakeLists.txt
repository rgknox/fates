
add_library(test_photo_lib
          FatesTestPhotosynthesisMod.F90)

set(gas_exchange_sources 
      FatesTestGasExchange.F90)

set(biophys_sources 
      FatesTestBiophysicalRates.F90)

set(leaf_photo_sources 
      FatesTestLeafPhoto.F90)

set(NETCDF_C_DIR ${NETCDF_C_PATH})
set(NETCDF_FORTRAN_DIR ${NETCDF_F_PATH})

FIND_PATH(NETCDFC_FOUND libnetcdf.a ${NETCDF_C_DIR}/lib)
FIND_PATH(NETCDFF_FOUND libnetcdff.a ${NETCDF_FORTRAN_DIR}/lib)

include_directories(${NETCDF_C_DIR}/include
                    ${NETCDF_FORTRAN_DIR}/include)

link_directories(${NETCDF_C_DIR}/lib
                ${NETCDF_FORTRAN_DIR}/lib
                ${PFUNIT_TOP_DIR}/lib)
                
add_executable(FATES_gasexch_exe ${gas_exchange_sources})
add_executable(FATES_biophysrates_exe ${biophys_sources})
add_executable(FATES_leaf_photo_exe ${leaf_photo_sources})

target_link_libraries(FATES_gasexch_exe
                  netcdf
                  netcdff
                  fates
                  csm_share
                  funit
                  test_photo_lib)
                  
target_link_libraries(FATES_biophysrates_exe
                  netcdf
                  netcdff
                  fates
                  csm_share
                  funit)
                  
target_link_libraries(FATES_leaf_photo_exe
                  netcdf
                  netcdff
                  fates
                  csm_share
                  funit
                  test_photo_lib)