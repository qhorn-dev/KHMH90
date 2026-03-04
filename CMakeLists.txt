cmake_minimum_required(VERSION 3.18)
project(khmh_s2_fortran LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

option(KHMH_USE_NETCDF "Enable NetCDF I/O" ON)

if(KHMH_USE_NETCDF)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(NETCDFF REQUIRED netcdf-fortran)
endif()

add_executable(khmh_s2
  src/mod_kinds.f90
  src/mod_config.f90
  src/mod_grid.f90
  src/mod_stats_types.f90
  src/mod_sf_terms.f90
  src/mod_linear_transport_terms.f90
  src/mod_interscale_transfer_terms.f90
  src/mod_interspace_transfer_terms.f90
  src/mod_pdelta_terms.f90
  src/mod_px_terms.f90
  src/mod_netcdf_io.f90
  src/mod_output_netcdf.f90
  src/main.f90
)

set_source_files_properties(src/main.f90 PROPERTIES Fortran_PREPROCESS ON)

if(KHMH_USE_NETCDF)
  target_include_directories(khmh_s2 PRIVATE ${NETCDFF_INCLUDE_DIRS})
  target_link_directories(khmh_s2 PRIVATE ${NETCDFF_LIBRARY_DIRS})
  target_link_libraries(khmh_s2 PRIVATE ${NETCDFF_LIBRARIES})
  target_compile_definitions(khmh_s2 PRIVATE KHMH_USE_NETCDF)
endif()

find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(khmh_s2 PRIVATE OpenMP::OpenMP_Fortran)
endif()
