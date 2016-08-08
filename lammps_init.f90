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
        double precision :: zmin, zmax
        double precision :: atom1_xmin, atom1_xmax, atom1_ymin, atom1_ymax
        integer :: nsteps
        integer :: atomType
        logical :: top, bot, left, right, itop, ibot, ileft, iright

        integer :: i, j, k, l, natoms, npad, n
        
	OPEN (UNIT=200,FILE='md.inp',STATUS='old')
	READ (200,*) stadium_width, exclude_top, exclude_bot, exclude_left, exclude_right
	READ (200,*) damp_coeff , damp_ratio
	READ (200,*) nh_dampcoeff
	READ (200,*) lammps_temperature
	READ (200,*) TIMestep1
	READ (200,*) INDextimeh
	READ (200,*) fem_update_steps
	READ (200,*) nsteps
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
	lammps_timestep = TIMestep1/1.0d-12
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
	lammps_output_steps = nsteps
        
        
        !! TODO get zmin and zmax from mat file automatically
        !!! --- For now assume single grain to be modified for multiple grains
	zmin = 0.0d0
	zmax = grains(1)%dcell(3)
        
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
        
        write(1010,*) "CADD input atoms"
        write(1010,*)
        write(1010, fmt='(I7,1X,A10)')  natoms, 'atoms'
        write(1010, fmt='(I3,1X,A15)')  4, 'atom types'
        write(1010, fmt='(2(1X,F15.8),1X,A15)')  xmin, xmax, 'xlo xhi '
        write(1010, fmt='(2(1X,F15.8),1X,A15)')  ymin, ymax, 'ylo yhi '
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

        double precision :: zmin, zmax
        character(1024):: command_line
        integer :: iatom, i,j,k,l

        real (C_double), pointer :: xlo => NULL()
        real (C_double), pointer :: xhi => NULL()
        real (C_double), pointer :: ylo => NULL()
        real (C_double), pointer :: yhi => NULL()

        zmin = 0.5
        zmax = 0.5



!!$        call lammps_open_no_mpi('lmg -log log.CADD', lmp)
!!$     If replacing this entire subroutine by reading input file
!!$        delete all lines in this file and replace by a lammps input file
        
!!$     call lammps_file(lmp, 'filename')

        call lammps_command(lmp, 'units metal')
        call lammps_command(lmp, 'atom_style atomic')
        call lammps_command(lmp, 'dimension 3')
        call lammps_command(lmp, 'boundary ss ss pp')
        call lammps_command(lmp,'atom_modify sort 0 0.0 map array')

        
        call lammps_command(lmp, 'read_data cadd_atoms.dat')

        ! --- Create groups of atoms for fixes and computes ----
        call lammps_command(lmp, "group md_atoms type 1 3 4")
        call lammps_command(lmp, "group free_atoms type 1")
        call lammps_command(lmp, "group pad_atoms type 2")
        call lammps_command(lmp, "group interface_atoms type 3")
        call lammps_command(lmp, "group langevin_atoms type 3 4")

        ! ------- EAM potentials
        call lammps_command(lmp, "pair_style eam/alloy")
!!$        call lammps_command(lmp, "pair_coeff	* * /home/srinath/lammps_potentials/Al-LEA_hex.eam.alloy Al Al Al")
        call lammps_command(lmp, "pair_coeff	* * Al_adams.eam.alloy Al Al Al Al")

        call lammps_command(lmp, "neighbor 0.1 bin ")
        call lammps_command(lmp, "neigh_modify delay 0 every 1 check yes")

        ! ---------- Various Fixes ----------------------------------------------
        write(command_line,*) "variable mytemp equal", lammps_temperature
        call lammps_command(lmp, command_line)
        !!$call lammps_command(lmp, "variable mytemp equal 500.0")
        call lammps_command(lmp, "velocity md_atoms create $(2.0*v_mytemp) 426789 dist uniform rot yes mom yes")
        if (nmaterials == 1) then 
	  if (material(1)%structure == 'hex') then 
	  ! ---- Set Z-velocity to zero for 2d hex problem
	    call lammps_command(lmp, "velocity md_atoms set NULL NULL 0.0 units box")
	  end if
	end if
	    
	  

	! ---- Temperature fixes -----------------------------------------------------------
	
!!$        call lammps_command(lmp, "fix fix_temp free_atoms nvt temp 1.0 1.0 100.0")


	! --- Todo put in langevin in the input file to say that there is a stadium ....
	call lammps_command(lmp, "fix fix_integ md_atoms nve")
    
!!$	write(command_line, fmt='(A38,3(1X,F15.6),I7, A10, 5(1X,F15.6))') "fix fix_temp langevin_atoms langevin ", &
!!$	  tstart, tstop, damp_coeff, 699483, " stadium ", stadium_xmin, stadium_xmax, stadium_ymin, stadium_ymax, stadium_width
!!$	call lammps_command(lmp, command_line)
    
    write(command_line, fmt='(A38,3(1X,F15.6),I7, A10, 7(1X,F15.6))') "fix fix_temp langevin_atoms langevin ", &
       tstart, tstop, damp_coeff, 699483, " stadium ", stadium_xmin, stadium_xmax, stadium_ymin, stadium_ymax, &
       -100000.00, 1000000.00, stadium_width
    call lammps_command(lmp, command_line)

	! ----------------------------------------------------------------------------------

	! ---- Temperature Computes and temperature variance for testing  -------------------------------------
	if (nmaterials == 1) then 
	  if (material(1)%structure == 'hex') then 
	    ! ---- This is a 2d Problem so the temperature compute is restricted to partial in the xy plane
	    call lammps_command(lmp, "compute free_temp free_atoms temp/partial 1 1 0")
	    call lammps_command(lmp, "compute stadium_temp langevin_atoms temp/partial 1 1 0")
	  else 
	    call lammps_command(lmp, "compute free_temp free_atoms temp")
	    call lammps_command(lmp, "compute stadium_temp langevin_atoms temp")	  
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
        !call lammps_command(lmp, "fix fix_2d all enforce2d")
        if (nmaterials == 1) then 
	  if (material(1)%structure == 'hex') then 
	    call lammps_command(lmp, "fix fix_2d all setforce NULL NULL 0.0")
	  end if 
	end if
        

        ! ------------- Various computes -------------------------------
        call lammps_command(lmp, "thermo 1")
        call lammps_command(lmp,"thermo_style custom step c_free_temp c_stadium_temp c_pe v_tot_energy")
        
        ! --------- Compute differential displacement from original position
        call lammps_command(lmp, "compute dx_free free_atoms displace/atom")

        call lammps_command(lmp, "compute dx_all all displace/atom")

        ! ---- Compute used for average displacement of interface atoms  -----
        call lammps_command(lmp, "compute dx_inter interface_atoms displace/atom")

        
        ! ----- Now define a fix to actually calculate the average for interface atoms
        write(command_line, '(A38, 2(1X,I3), A15)') "fix dx_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_inter[1]"
        call lammps_command(lmp, command_line)
!!$        call lammps_command(lmp, "fix dx_ave interface_atoms ave/atom 1 25 25 c_dx_inter[1]")
        write(command_line,  '(A38, 2(1X,I3), A15)') "fix dy_ave interface_atoms ave/atom 1 ", fem_update_steps, fem_update_steps, " c_dx_inter[2]"
        call lammps_command(lmp, command_line)
!!$        call lammps_command(lmp, "fix dy_ave interface_atoms ave/atom 1 25 25 c_dx_inter[2]")
!!$        call lammps_command(lmp, "fix dz_ave interface_atoms ave/atom 1 25 25 c_dx_inter[3]")


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
        write(command_line, '(A18,I3,A74)') "dump 1 all custom ", lammps_output_steps, " atom_lmp*.cfg id type x y z c_dx_all[1] c_dx_all[2] fx fy fz c_dx_all[4]"
        call lammps_command(lmp, command_line)
!!$        call lammps_command(lmp, "dump 1 all custom 200 atom_lmp*.cfg id type x y z c_dx_all[1] c_dx_all[2] fx fy fz")       
        
        call lammps_command(lmp, "run 0")
       
      end subroutine initialize_lammps