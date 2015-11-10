cmake_minimum_required(VERSION 2.6)

project(CADD_Modular)
enable_language(Fortran)
enable_language(CXX OPTIONAL)

# list of source files
#FILE(GLOB SRC_FILES *.f90)
set(SRC_FILES
  atomistic.f90
  bandwidth.f90
  cgma05.f90
  common.f90
  contri.f90
  delaunay36.f90
  dislocation.f90
  energygr.f90
  feaplib.f90
  md.f90
  md_lammps.f90
  mod_lammps.f90
  mod_boundary.f90
  mod_cluster.f90
  mod_dynamo.f90
  mod_file.f90
  mod_global.f90
  mod_grain.f90
  mod_interstitial.f90
  mod_lammps.f90
  mod_main.f90
  mod_material.f90
  mod_output.f90
  mod_poten_eam.f90
  pmacr1.f90
  pmacr3.f90
  pmain.f90
  pmesh.f90
  polygrain.f90
  randgen.f90
  restart.f90
  tecplotnl.f90
  )

INCLUDE_DIRECTORIES(../lammps
  ${CMAKE_BINARY_DIR}/lammps
  ../Disl
  ${CMAKE_BINARY_DIR}/Disl
  /home/srinath/lammps.git/src/STUBS)


# Create Library for linking into main executable
set(lammps /home/srinath/lammps.git/src/liblammps_intel.so)
set(stubs /home/srinath/lammps.git/src/STUBS/libmpi_stubs.a)

add_library(mdlro SHARED ${SRC_FILES} ${INC_FILES})
target_link_libraries(mdlro lammpsf90 lammpscpp ${lammps} ${stubs})

# set the linker flags that are not libraries

