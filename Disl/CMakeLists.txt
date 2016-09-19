cmake_minimum_required(VERSION 2.6)

project(CADD_Disl)
enable_language(Fortran)
enable_language(CXX OPTIONAL)

set(CADD_Disl_MODULES)

# list of source files
set(SRC_FILES
  disl_cg.f90
  disl_dd.f90
  disl_field.f90
  disl_parameters.f90
  disl_routines.f90
  dumpdisl.f90
  fem_alan.f90
  fem_elastic.f90
  fem_movepad.f90
  fem_parameters.f90
  fem_routines.f90
  fem_service.f90
  mod_dd_slip.f90
  )

#FILE(GLOB INC_FILES *.par)

INCLUDE_DIRECTORIES(
  lammps
  ${CMAKE_BINARY_DIR}/lammps
  Disl
  ${CMAKE_BINARY_DIR}/Disl
  )




### Checking to see if compile flags are inherited
# Create Library for linking into main executable
add_library(dislo SHARED ${SRC_FILES} ${INC_FILES})

#add_library(dislo STATIC ${SRC_FILES} ${INC_FILES})




#if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#  message(FATAL_ERROR "Sorry qsdd does not compile with gfortran. Please choose an other compiler by setting the variable CMAKE_Fortran_COMPILER")
#endif()
