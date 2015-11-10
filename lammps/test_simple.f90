subroutine test_lamms(lmp)
  use LAMMPS
  implicit none
  type(C_ptr) :: lmp
  character (len=1024) :: command_line
  real (C_double), pointer :: compute => NULL()
  real (C_double) :: fix, fix2
  real (C_double), dimension(:), pointer :: compute_v => NULL()
  real (C_double), dimension(:,:), pointer :: x => NULL()
  real (C_double), dimension(:), pointer :: mass => NULL()
  integer, dimension(:), allocatable :: types
  double precision, dimension(:), allocatable :: r
  integer :: error, narg, me, nprocs


  call lammps_file (lmp, 'in.simple')
  call lammps_command (lmp, 'run 500')

  ! This extracts f_2 as a scalar (the last two arguments can be arbitrary)
  call lammps_extract_fix (fix, lmp, '2', 0, 0, 1, 1)
  print *, 'Fix is ', fix

!!$   ! This extracts f_4[1][1] as a scalar
!!$   call lammps_extract_fix (fix2, lmp, '4', 0, 2, 1, 1)
!!$   print *, 'Fix 2 is ', fix2

  ! This extracts the scalar compute of compute thermo_temp
  call lammps_extract_compute (compute, lmp, 'thermo_temp', 0, 0)
  print *, 'Compute is ', compute

  ! This extracts the vector compute of compute thermo_temp
  call lammps_extract_compute (compute_v, lmp, 'thermo_temp', 0, 1)
  print *, 'Vector is ', compute_v

  ! This extracts the masses
  call lammps_extract_atom (mass, lmp, 'mass')
  print *, 'Mass is ', mass(1:)

  ! Extracts a pointer to the arrays of positions for all atoms
  call lammps_extract_atom (x, lmp, 'x')
  if ( .not. associated (x) ) print *, 'x is not associated'
  print *, 'x is ', x(:,1)  ! Prints x, y, z for atom 1

  ! Extracts pointer to atom types
  call lammps_gather_atoms (lmp, 'type', 1, types)
  print *, 'types is ', types(1:3)

  ! Allocates an array and assigns all positions to it
  call lammps_gather_atoms (lmp, 'x', 3, r)
  print *, 'size(r) = ', size(r)
  print *, 'r is ', r(1:6)

  ! Puts those position data back
  call lammps_scatter_atoms (lmp, 'x', r)

  call lammps_command (lmp, 'run 1')
  print *, 'x is ', x(:,1)  ! Note that the position updates!
  print *, 'Compute is ', compute  ! This did only because "temp" is part of
  ! the thermo output; the vector part did
  ! not, and won't until we give LAMMPS a
  ! thermo output or other command that
  ! requires its value


end subroutine test_lamms
