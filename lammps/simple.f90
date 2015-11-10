program simple

!!$   use MPI
   use LAMMPS

   ! The following line is unnecessary, as I have included these three entities
   ! with the LAMMPS module, but I leave them in anyway to remind people where
   ! they came from
   use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int

   implicit none

   ! Notes:
   !  * If LAMMPS returns a scalar that is allocated by the library interface
   !     (see library.cpp), then that memory is deallocated automatically and
   !     the argument to lammps_extract_fix must be a SCALAR.
   !  * If LAMMPS returns a pointer to an array, consisting of internal LAMMPS
   !     data, then the argument must be an interoperable Fortran pointer.
   !     Interoperable means it is of type INTEGER (C_INT) or of type
   !     REAL (C_DOUBLE) in this context.
   !  * Pointers should NEVER be deallocated, as that would deallocate internal
   !     LAMMPS data!
   !  * Note that just because you can read the values of, say, a compute at
   !     any time does not mean those values represent the "correct" values.
   !     LAMMPS will abort you if you try to grab a pointer to a non-current
   !     entity, but once it's bound, it's your responsibility to check that
   !     it's current before evaluating.
   !  * IMPORTANT:  Two-dimensional arrays (such as 'x' from extract_atom)
   !     will be transposed from what they might look like in C++.  This is
   !     because of different bookkeeping conventions between Fortran and C
   !     that date back to about 1970 or so (when C was written).
   !  * Arrays start from 1, EXCEPT for mass from extract_atom, which
   !     starts from 0.  This is because the C array actually has a blank
   !     first element (and thus mass[1] corresponds to the mass of type 1)

   type (C_ptr) :: lmp
   real (C_double), pointer :: compute => NULL()
   real (C_double) :: fix, fix2
   real (C_double), dimension(:), pointer :: compute_v => NULL()
   real (C_double), dimension(:,:), pointer :: x => NULL()
   real (C_double), dimension(:), pointer :: mass => NULL()
   integer, dimension(:), allocatable :: types
   double precision, dimension(:), allocatable :: r
   integer :: error, narg, me, nprocs
   character (len=1024) :: command_line

!!$   call MPI_Init (error)
!!$   call MPI_Comm_rank (MPI_COMM_WORLD, me, error)
!!$   call MPI_Comm_size (MPI_COMM_WORLD, nprocs, error)

   ! You are free to pass any string you like to lammps_open or
   ! lammps_open_no_mpi; here is how you pass it the command line
   !call get_command (command_line)
   !call lammps_open (command_line, MPI_COMM_WORLD, lmp)

   ! And here's how to to it with a string constant of your choice
   call lammps_open_no_mpi ('lmp -log log.simple', lmp)
   call test_lammps(lmp)

   call lammps_close (lmp)

!!$   call MPI_Finalize (error)

end program simple
