subroutine write_lammps_data(Id, X, Ix, F, B, Itx, xmin, xmax, ymin, ymax)
  !> Subroutine to write initial atomic config into lammps format for lammps reading
  !>      Identifies atom types, free atoms --> type 1, pad atoms --> 2
  ! TODO (Srinath#1#): Check to ensure data format can be read by lammps ...
  !    1) write the data file
  !    2) Check with independent lammps script.
  !    3) Include type for interface atoms
  !    4) All atom types from CADD are identified by a max number and the indenter etc will
  !        other atom types.

  USE MOD_GRAIN
  USE MOD_GLOBAL
  USE MOD_FILE
  USE MOD_BOUNDARY
  USE MOD_CRACK
  USE MOD_MATERIAL
  USE MOD_DD_SLIP
  USE MOD_DISL_PARAMETERS
  implicit none
  !     Input variables
  INTEGER Id , Ix , Itx(3,*)
  DOUBLE PRECISION X , F , B
  DIMENSION Id(NDF,1) , X(NXDm,1) , Ix(NEN1,1) , F(NDF,1) , B(NDF,1)
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  double precision :: zmin, zmax, yparticle
  double precision :: atom1_xmin, atom1_xmax, atom1_ymin, atom1_ymax
  integer :: nsteps
  integer :: atomType
  logical :: top, bot, left, right, itop, ibot, ileft, iright

  integer :: i, j, k, l, natoms, npad, n

  OPEN (UNIT=200,FILE='md.inp',STATUS='old')
  READ (200,*) stadium_width, exclude_top, exclude_bot, exclude_left, exclude_right
  READ (200,*) damp_coeff , damp_ratio
  READ (200,*) lammps_temperature
  READ (200,*) lammps_timestep
  READ (200,*) fem_update_steps
  READ (200,*) num_md_steps
  READ (200,*) lammps_output_steps
  READ (200,*) num_restart_steps
  READ (200,*) num_initial_equil
  READ (200,*) particle_velocity
  READ (200,*) particle_radius
  READ (200,*) particle_height
  READ (200,*) particle_rotation
  READ (200,*) impact_angle
  CLOSE (200)


  top = .true.
  bot = .true.
  left = .true.
  right = .true. 

  if (exclude_top > 0) then
     top = .false.
  end if

  if (exclude_bot > 0) then
     bot = .false.
  end if

  if (exclude_left > 0) then
     left = .false.
  end if

  if (exclude_right > 0) then
     right = .false.
  end if


!!!! --- Assume units metal in lammps
  lammps_timestep = lammps_timestep/1.0d-12
  tstart = lammps_temperature
  tstop = lammps_temperature

!!! TODO add exclusion zones for stadium thermostat such as top surface, bottom surface, left, right
!!! --- mainly for impact problem

  if (damp_coeff > 0.0) then 
     damp_coeff = 1.0/(damp_coeff)
     ! TODO error handler 
  end if
!!! Once again assuming metal units in lammps
  damp_coeff = damp_coeff / 1.0d-12


  if ( NGRains/=1 ) stop 'hardwired here for 1 grain for zmin/zmax'

  !!-if material is hex, hardwire the z-spacing
  !!-if material is fcc, read in z-spacing from grain data
  if (nmaterials == 1) then 
     if (material(1)%structure == 'hex') then 
        !!> TODO get zmin and zmax from mat file automatically
        zmin = 0.0d0
        zmax = 2.9573845299d0
     else 
        !! TODO get zmin and zmax from mat file automatically
        !! --- For now assume single grain to be modified for multiple grains
        !!zmax = 0.0d0
        zmax = grains(1)%dcell(3)/2.0
        zmin = -grains(1)%dcell(3)/2.0
     end if
  end if

  !!		JM: ***Must be careful to allow enough y space to add particle***
  !!		JM: ***over free surface***        
  yparticle = -ymin

  open(unit=1010, file='cadd_atoms.dat', status='UNKNOWN')
  natoms = 0
  atom1_xmin = 1.0d20
  atom1_xmax = -1.0d20
  atom1_ymin = 1.0d20
  atom1_ymax = -1.0d20

  do i = 1, numnp
     if (isRelaxed(i) /= 0) then
        natoms = natoms + 1
        if (isRelaxed(i) /= -1) then 
           if (X(1,i) < atom1_xmin) then 
              atom1_xmin = X(1,i)
           end if
           if (X(2,i) < atom1_ymin) then 
              atom1_ymin = X(2,i)
           end if
           if (X(1,i) > atom1_xmax) then 
              atom1_xmax = X(1,i)
           end if
           if (X(2,i) > atom1_ymax) then 
              atom1_ymax = X(2,i)
           end if
        end if
     end if
     if (isRelaxed(i) == -1) then
        npad = npad + 1
     end if
  end do
  stadium_xmax = atom1_xmax
  stadium_xmin = atom1_xmin
  stadium_ymin = atom1_ymin
  stadium_ymax = atom1_ymax


  if (particle_radius > abs(atom1_xmax - atom1_xmin)/2.0) then
     call error_handler('Particle radius must be smaller than substrate dimensions')
  end if

  if (particle_radius > abs(atom1_ymax - atom1_ymin)/2.0) then
     call error_handler('Particle radius must be smaller than substrate dimensions')
  end if

  write(1010,*) "CADD input atoms"
  write(1010,*)
  write(1010, fmt='(I7,1X,A10)')  natoms, 'atoms'
  write(1010, fmt='(I3,1X,A15)')  5, 'atom types'
  write(1010, fmt='(2(1X,F15.8),1X,A15)')  xmin, xmax, 'xlo xhi '
  write(1010, fmt='(2(1X,F15.8),1X,A15)')  ymin, yparticle, 'ylo yhi '
  write(1010, fmt='(2(1X,F15.8),1X,A15)')  zmin, zmax, 'zlo zhi '
  write(1010, *)
  write(1010,*) 'Atoms'
  write(1010,*)
  n = 0
  do i = 1, numnp
     if (isRelaxed(i) /=0 ) then
        n = n + 1
        if (isRelaxed(i) == -1) then
           atomType = 2
        elseif (isRelaxed(i) == 2) then
           atomType = 3
        else
           itop = (x(2,i) > atom1_ymax - stadium_width)
           ibot = (x(2,i) < atom1_ymin + stadium_width)
           ileft = (x(1,i) < atom1_xmin + stadium_width)
           iright = (x(1,i) > atom1_xmax - stadium_width)
           if ((itop .and. top) .or. (ibot .and. bot) .or. (ileft .and. left) .or. (iright .and. right)) then 
!!$		 if (X(1,i) > atom1_xmax - stadium_width .or. X(1,i) < atom1_xmin+stadium_width & 
!!$		    .or. X(2,i) > atom1_ymax - stadium_width .or. X(2,i) < atom1_ymin+stadium_width) then 
              atomType = 4
           else 
              atomType = 1
           end if
        end if
        write(1010,fmt='(I7,1X,I3,1X,3(1X,F15.8))') n, atomType, X(1,i), X(2,i), X(3,i)
     end if
  end do

!!$ no stadium on free surface
  stadium_ymax = 2.0d0*stadium_width

  write(1010,*)
  close(1010)
end subroutine write_lammps_data


subroutine initialize_lammps(Id,X,Ix,F,B,Itx,lmp)
  !> Subroutine to initialize lammps
  !>  This performs the following functions
  !>     Initialize lammps pointer
  !>     Initialize dimension to 2
  !>     Initialize potential for use
  !>     Reads the initial mesh from mesh.f90 through a data file
  ! TODO (Srinath#1#): 1) Make data file name avaialble to this routine  ...
  !2) Fixed atoms identified with atom type
  !3) Compare initial configs between this and CADD

  use lammps
  USE MOD_GRAIN
  USE MOD_GLOBAL
  USE MOD_FILE
  USE MOD_BOUNDARY
  USE MOD_CRACK
  USE MOD_MATERIAL
  USE MOD_DD_SLIP
  USE MOD_DISL_PARAMETERS

  implicit none

  type (c_ptr) :: lmp
  !     Input variables
  INTEGER Id , Ix , Itx(3,*)
  DOUBLE PRECISION X , F , B
  DIMENSION Id(NDF,1) , X(NXDm,1) , Ix(NEN1,1) , F(NDF,1) , B(NDF,1)
!!$        double precision, intent(in) :: xmin, xmax, ymin, ymax, padwidth

  double precision :: zmin, zmax, theta, xoffset, Pi, ia_rad, zbound, pbc_tol
  double precision :: xperiod, yperiod, zperiod, a0
  character(1024):: command_line
  integer :: iatom, i,j,k,l

  real (C_double), pointer :: xlo => NULL()
  real (C_double), pointer :: xhi => NULL()
  real (C_double), pointer :: ylo => NULL()
  real (C_double), pointer :: yhi => NULL()


!!$ --- For now assume single grain/material to be modified for multiple grains
  if (nmaterials == 1) then 
     if (material(1)%structure == 'hex') then 
        zperiod = grains(1)%dcell(3)
     else 
        xperiod = grains(1)%dcell(1)
        yperiod = grains(1)%dcell(2)
        zperiod = grains(1)%dcell(3)
     end if
  end if

!!$        call lammps_open_no_mpi('lmg -log log.CADD', lmp)
!!$     If replacing this entire subroutine by reading input file
!!$        delete all lines in this file and replace by a lammps input file

!!$     call lammps_file(lmp, 'filename')

  call lammps_command(lmp, 'units metal')
  call lammps_command(lmp, 'atom_style atomic')
  call lammps_command(lmp, 'dimension 3')
  call lammps_command(lmp, 'boundary ss sm pp')
!!$        call lammps_command(lmp, 'boundary ff ff pp')
  call lammps_command(lmp,'atom_modify sort 0 0.0 map array')

  !!JM: hcp lattice with short z cylinder region length gives a planar hex lattice
  !!JM: with a little hack using the burger's vector as the lattice constant a
  if (nmaterials == 1) then 
     a0 = material(1)%a0

     if (material(1)%structure == 'hex') then 
        !---- This is a 2d Problem so the temperature compute is restricted to partial in the xy plane
        write(command_line,'(A,1F15.5)') 'lattice hex ', a0
        call lammps_command(lmp, command_line)
        !call lammps_command(lmp, 'lattice hex 2.85105')
     else 
        !!call lammps_command(lmp, 'lattice fcc 4.032')
        !!write(command_line,'(A,1F15.5)') 'lattice hex ', a0
        !!call lammps_command(lmp, command_line)
        !!FCC orientation to prevent screw dislocations
        !!call lammps_command(lmp, 'lattice fcc 4.032 orient x 1 1 -2 orient y 1 1 1 orient z 1 -1 0')
        !!call lammps_command(lmp, 'lattice fcc 4.032 orient x 1 -1 0 orient y 1 1 1 orient z -1 -1 2')
        !!write(command_line,'(A,1F15.5,A,3F15.3)') 'lattice fcc ', a0, ' orient x 1 -1 0 orient y 1 1 1 orient z -1 -1 2 spacing ', &
        !!   xperiod, yperiod, zperiod
        print *, "Lattice constant in lammps = ", a0
        write(command_line,'(A,1F15.5,A,3F15.3)') 'lattice fcc ', a0, ' orient x 1 1 -2 orient y 1 1 1 orient z 1 -1 0 spacing ', &
           xperiod, yperiod, zperiod
        !!write(command_line,'(A,1F15.5,A,3F15.3)') 'lattice fcc ', a0, ' orient x 1 0 0 orient y 0 1 0 orient z 0 0 1 spacing ', &
        !!     xperiod, yperiod, zperiod
        call lammps_command(lmp, command_line)

     end if
  end if

  call lammps_command(lmp, 'read_data cadd_atoms.dat')

!!!!!!!!!!!!!FIX THIS!!!!!!!!!!!!!!!!
  !store number of zperiods from CADD into lammps variable 
  !to be used when creating the particle
  print*, 'numperiodz is: ', NUMperiodz
  !write(command_line,'(A,1F15.3)') 'variable zperiods equal ', NUMperiodz
  write(command_line,'(A,1F15.3)') 'variable zperiods equal 1.0'
  call lammps_command(lmp,command_line)

  !for no particle impact angle, differentiate hex vs. fcc
  if (nmaterials == 1) then 
     if (material(1)%structure == 'hex') then 
        write(command_line,'(A,2F15.3,A)') 'region 1 cylinder z 0.0 ', &
             particle_height + particle_radius, particle_radius, ' 0.0 0.1 units box'
        call lammps_command(lmp, command_line)
     else 
!!!!!!FIGURE OUT THIS 2.0 vs 1.9 tolerance issue to line up CADD and LAMMPS lattices!!
        write(command_line,'(A,4F15.3,A)') 'region 1 cylinder z 0.0 ', &
             particle_height + particle_radius, particle_radius, -grains(1)%dcell(3)/2.0-1.d-3, grains(1)%dcell(3)/2.0+1.d-3, ' units box'
        call lammps_command(lmp, command_line)
     end if
  end if


  !! --- creating particle atoms, type 5---
  call lammps_command(lmp, 'create_atoms 5 region 1')   
  call lammps_command(lmp, 'group particle_atoms type 5')

  ! --- Create groups of atoms for fixes and computes ----
  call lammps_command(lmp, "group md_atoms type 1 3 4 5")
  call lammps_command(lmp, "group sub_atoms type 1 3 4")
  call lammps_command(lmp, "group free_atoms type 1")
  call lammps_command(lmp, "group pad_atoms type 2")
  call lammps_command(lmp, "group interface_atoms type 3")
  call lammps_command(lmp, "group langevin_atoms type 3 4")

  !---------------------------------------------------------------------------------------------       
  ! ------ Code to correctly compute the differential displacement ----
  ! --- Currently only z direction (periodic direction is implemented) ----
  ! --- Can be easily extended to other dimensions

  !!!call lammps_command(lmp,'fix dz all property/atom d_dz')
  call lammps_command(lmp,'compute vz_all all property/atom vz')
  call lammps_command(lmp,'compute vz_inter interface_atoms property/atom vz')
  call lammps_command(lmp,'variable dz_all atom "c_vz_all*dt"')
  call lammps_command(lmp,'variable dz_inter atom "c_vz_inter*dt"')
  !call lammps_command(lmp, 'compute dz sub_atoms property/atom d_dz')
  !---------------------------------------------------------------------------------------------       

  
  ! ------- EAM potentials
  call lammps_command(lmp, "pair_style eam/alloy")
!!$        call lammps_command(lmp, "pair_coeff	* * /home/srinath/lammps_potentials/Al-LEA_hex.eam.alloy Al Al Al")
  call lammps_command(lmp, "pair_coeff	* * Al_adams.eam.alloy Al Al Al Al Al")

  call lammps_command(lmp, "neighbor 0.1 bin ")
  call lammps_command(lmp, "neigh_modify delay 0 every 1 check yes")

  ! ---------- Various Fixes ----------------------------------------------
  write(command_line,*) "variable mytemp equal", lammps_temperature
  call lammps_command(lmp, command_line) 

  !!--- initialize system temperature to 2*desired (NVE ensembles << melting temp)
  if (nmaterials == 1) then 
     if (material(1)%structure == 'hex') then 
        call lammps_command(lmp, "velocity sub_atoms create $(2.0*v_mytemp) 829863 dist uniform mom yes rot yes")
        call lammps_command(lmp, "velocity particle_atoms create $(2.0*v_mytemp) 163103 dist uniform mom yes rot yes")
        !!--- If hex material, then set z velocity to 0
        call lammps_command(lmp, "velocity sub_atoms set NULL NULL 0.0 units box")
        call lammps_command(lmp, "velocity particle_atoms set NULL NULL 0.0 units box")
     else 
        call lammps_command(lmp, "velocity sub_atoms create $(2.0*v_mytemp) 426789 dist uniform mom yes rot yes")
        call lammps_command(lmp, "velocity particle_atoms create $(2.0*v_mytemp) 849823 dist uniform mom yes rot yes")
     end if
  end if

  !! --- equilibrate temperature in two NVT ensembles (1 for particle, 1 for substrate)
  call lammps_command(lmp, "fix int_sub sub_atoms nve/limit 0.5")
  call lammps_command(lmp, "fix int_part particle_atoms nve/limit 0.5")

  write(command_line, fmt='(A38,3(1X,F15.6),I7, A10, 7(1X,F15.6))') "fix fix_temp langevin_atoms langevin ", &
       tstart, tstop, damp_coeff, 699483, " stadium ", stadium_xmin, stadium_xmax, stadium_ymin, stadium_ymax, &
       -1000.00, 1000.00, stadium_width
  call lammps_command(lmp, command_line)

  ! ----------------------------------------------------------------------------------

  ! ---- Temperature Computes and temperature variance for testing  -------------------------------------
  if (nmaterials == 1) then 
     if (material(1)%structure == 'hex') then 
        ! ---- This is a 2d Problem so the temperature compute is restricted to partial in the xy plane
        call lammps_command(lmp, "compute md_temp md_atoms temp/partial 1 1 0")
        call lammps_command(lmp, "compute sub_temp sub_atoms temp/partial 1 1 0")
        call lammps_command(lmp, "compute free_temp free_atoms temp/partial 1 1 0")
        call lammps_command(lmp, "compute stadium_temp langevin_atoms temp/partial 1 1 0")
        call lammps_command(lmp, "compute part_temp particle_atoms temp/partial 1 1 0")
     else 
        call lammps_command(lmp, "compute md_temp md_atoms temp/com")
        call lammps_command(lmp, "compute sub_temp sub_atoms temp/com")
        call lammps_command(lmp, "compute free_temp free_atoms temp/com")
        call lammps_command(lmp, "compute stadium_temp langevin_atoms temp/com")
        call lammps_command(lmp, "compute part_temp particle_atoms temp/com")
     end if
  end if

  ! ---- average the variance of the temperature 
!!$        call lammps_command(lmp, "fix free_variance free_atoms ave/time 1 1000 1000 v_delt_free ave running")
!!$        call lammps_command(lmp, "fix stadium_variance free_atoms ave/time 1 1000 1000 v_delt_stadium ave running")

!!$        call lammps_command(lmp, 'fix print_variance all print 1000 "Temperature Variance =  &
!!$	   $(f_free_variance) $(f_stadium_variance) $(f_free_variance/v_canonical_free) $(f_stadium_variance/v_canonical_stadium)"')


  ! ------------ Energies -------------------------
  call lammps_command(lmp, "compute com_pe free_atoms pe/atom")
  call lammps_command(lmp, "compute pe free_atoms reduce sum c_com_pe")
  call lammps_command(lmp, "compute ke free_atoms ke")
  call lammps_command(lmp, "variable tot_energy equal c_pe+c_ke")
  ! ------------------------------------------------------------------------

  !---- Pad atoms always have zero force so this is fixed here to 0 
  call lammps_command(lmp, "fix fix_zeroforce pad_atoms setforce 0.0 0.0 0.0")

!!$ 2D material fixes (hex)    
  if (nmaterials == 1) then 
     if (material(1)%structure == 'hex') then 
        call lammps_command(lmp, "fix fix_2d all setforce NULL NULL 0.0")
        call lammps_command(lmp, "fix fix_2d all enforce2d")
     end if
  end if


  ! ------------- Various computes -------------------------------
  call lammps_command(lmp, "thermo 1")
  call lammps_command(lmp,"thermo_style custom step c_md_temp c_free_temp c_stadium_temp c_part_temp c_sub_temp c_pe c_ke vol")

  ! --------- Compute differential displacement from original position
  call lammps_command(lmp, "set group all image 0 0 0")
  call lammps_command(lmp, "compute dx_free free_atoms displace/atom")

  call lammps_command(lmp, "compute dx_sub sub_atoms displace/atom")

  ! ---- Compute used for average displacement of interface atoms  -----
  call lammps_command(lmp, "compute dx_inter interface_atoms displace/atom")


  ! ----- Now define a fix to actually calculate the average for interface atoms
  write(command_line, '(A38, 2(1X,I3), A15)') "fix dx_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_inter[1]"
  call lammps_command(lmp, command_line)
  write(command_line,  '(A38, 2(1X,I3), A15)') "fix dy_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_inter[2]"
  call lammps_command(lmp, command_line)
  write(command_line,  '(A38, 2(1X,I3), A15)') "fix dz_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " v_dz_inter"
  call lammps_command(lmp, command_line)

  ! ----- Now define a fix to actually calculate the average for all substrate atoms
  write(command_line, '(A38, 2(1X,I3), A15)') "fix dx_sub_ave sub_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_sub[1]"
  call lammps_command(lmp, command_line)
  write(command_line,  '(A38, 2(1X,I3), A15)') "fix dy_sub_ave sub_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_sub[2]"
  call lammps_command(lmp, command_line)
  write(command_line,  '(A38, 2(1X,I3), A15)') "fix dz_sub_ave sub_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " v_dz_all"
  call lammps_command(lmp, command_line)




  ! ---- Compute used for virial stress on atoms
  call lammps_command(lmp, "compute compute_stress free_atoms stress/atom NULL virial")


  ! ----- Now define a fix to actually calculate the stress average for all atoms
  call lammps_command(lmp, "fix stress_ave_xx free_atoms ave/atom 1 25 25 c_compute_stress[1]")
  call lammps_command(lmp, "fix stress_ave_yy free_atoms ave/atom 1 25 25 c_compute_stress[2]")
  call lammps_command(lmp, "fix stress_ave_zz free_atoms ave/atom 1 25 25 c_compute_stress[3]")
  call lammps_command(lmp, "fix stress_ave_xy free_atoms ave/atom 1 25 25 c_compute_stress[4]")
  call lammps_command(lmp, "fix stress_ave_zx free_atoms ave/atom 1 25 25 c_compute_stress[5]")
  call lammps_command(lmp, "fix stress_ave_yz free_atoms ave/atom 1 25 25 c_compute_stress[6]")

  call lammps_extract_global(xlo, lmp, 'boxxlo')
  call lammps_extract_global(xhi, lmp, 'boxxhi')
  call lammps_extract_global(ylo, lmp, 'boxylo')
  call lammps_extract_global(yhi, lmp, 'boxyhi')


  ! ---- Dump data file 
  write(command_line, '(A18,I3,A91)') "dump 1 all custom ", lammps_output_steps, " atom_lmp*.cfg id type x y z c_dx_sub[1] c_dx_sub[2] c_dx_sub[3] fx fy fz vx vy vz v_dz_all"
  call lammps_command(lmp, command_line)
!!$        call lammps_command(lmp, "dump 1 all custom 200 atom_lmp*.cfg id type x y z c_dx_all[1] c_dx_all[2] fx fy fz")       

  call lammps_command(lmp, "run 0")
  call lammps_command(lmp, "write_restart restart.cadd")

end subroutine initialize_lammps
