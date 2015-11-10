cmake_minimum_required(VERSION 2.6)

project(CADD)
enable_language(Fortran)
enable_language(CXX OPTIONAL)

# list of source files
INCLUDE_DIRECTORIES(  lammps
  ${CMAKE_BINARY_DIR}/lammps
  Disl
  ${CMAKE_BINARY_DIR}/Disl
  Modular
  ${CMAKE_BINARY_DIR}/Modular
 )

#FILE(GLOB SRC_FILES *.f90)
#FILE(GLOB INC_FILES *.par)
set(SRC_FILES
  lammps_init.f90
  mesh.f90
  field.f90
  mod_crack.f90
  move_crack.f90
  )

SUBDIRS(lammps
  Disl 
  Modular)

#===============================================================================
#Macro to add some compilation flags
macro(add_ff_flags flag)

  set(type ${ARGN})
  if(type)
    set(_var CMAKE_Fortran_FLAGS_${type})
  else()
    set(_var CMAKE_Fortran_FLAGS)
  endif()
  if(NOT ${_var} MATCHES "${flag} ")
    set(${_var} "${flag} ${${_var}}"
      CACHE STRING "Flags used by the compiler during all build types." FORCE)
  endif()
endmacro()

macro(remove_ff_flags flag)
  set(type ${ARGN})
  if(type)
    set(_var CMAKE_Fortran_FLAGS_${type})
  else()
    set(_var CMAKE_Fortran_FLAGS)
  endif()
  string(REPLACE "${flag} " "" ${_var} "${${_var}}")
  set(${_var} "${${_var}}"
    CACHE STRING "Flags used by the compiler during all build types." FORCE)
endmacro()

macro(add_sh_flags flag)
  set(type ${ARGN})
  if(type)
    set(_var CMAKE_SHARED_LINKER_FLAGS_${type})
  else()
    set(_var CMAKE_SHARED_LINKER_FLAGS)
  endif()
  if(NOT ${_var} MATCHES "${flag} ")
    set(${_var} "${flag} ${${_var}}"
      CACHE STRING "Flags used by the compiler during all build types." FORCE)
  endif() 
endmacro()

macro(add_exe_flags flag)
  set(type ${ARGN})
  if(type)
    set(_var CMAKE_EXE_LINKER_FLAGS_${type})
  else()
    set(_var CMAKE_EXE_LINKER_FLAGS)
  endif()
  if(NOT ${_var} MATCHES "${flag} ")
    set(${_var} "${flag} ${${_var}}"
      CACHE STRING "Flags used by the compiler during all build types." FORCE)
  endif() 
endmacro()

#===============================================================================

#  message("adding flags")
include_directories(${CMAKE_BINARY_DIR})

# to add ff flags for a given build type
add_ff_flags(-fPIC)
add_ff_flags(-cxxlib)
#add_ff_flags(-lm)
add_sh_flags(-lz)
#add_sh_flags(-lstdc++)
#add_sh_flags(-lm)

#add_exe_flags(-lstdc++)
#add_exe_flags(-lm)
add_exe_flags(-lz)

set(lammpsf90 ${CMAKE_SOURCE_DIR}/lammps/LAMMPS.o)
set(lammpscpp ${CMAKE_SOURCE_DIR}/lammps/LAMMPS-wrapper.o)
set(lammps /home/srinath/lammps.git/src/liblammps_intel.so)
set(stubs /home/srinath/lammps.git/src/STUBS/libmpi_stubs.a)


# Add executable 
add_executable(CADD ${SRC_FILES} ${INC_FILES})
target_link_libraries(CADD dislo mdlro lammpsf90 lammpscpp ${lammps} ${stubs})

#if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#  message(FATAL_ERROR "Sorry CADD does not compile with gfortran. Please choose an other compiler by s#etting the variable CMAKE_Fortran_COMPILER")
#endif()
