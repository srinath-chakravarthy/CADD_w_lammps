cmake_minimum_required(VERSION 2.6)

project(CADD)
enable_language(Fortran)
enable_language(CXX)

# list of source files
set(FORTRAN_SRC_FILES
  lammps.f90)
set(CPP_SRC_FILES
  lammps-wrapper.cpp)

INCLUDE_DIRECTORIES(
  /home/srinath/lammps.git/src
  /home/srinath/lammps.git/src/STUBS
  )


# Create Library for linking into main executable

add_library(lammpsf90 SHARED ${FORTRAN_SRC_FILES} ${INC_FILES})
add_library(lammpscpp SHARED ${CPP_SRC_FILES} ${INC_FILES})

