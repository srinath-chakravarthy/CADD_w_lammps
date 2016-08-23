      SUBROUTINE GETFEM_FORCES(Id,Atomcoord,Ix,F,Atomdispl,&
     &                                Avedispl,Atomforce,Atommass,&
     &                                Systemenergy,Moveatoms,Movedisl,&
     &                                Fullfield,Solvefem,Straine0,Ifem)

      USE MOD_GLOBAL
      USE MOD_TIMMING
      USE LAMMPS
      USE MOD_LAMMPS
      use mod_file
      IMPLICIT NONE
!*--GETENERGIESANDFORCES13


      DOUBLE PRECISION Atomdispl(NDF,*) , Atomcoord(NDF,*) , F(NDF,*) , &
     &                 Atomforce(NDF,*) , Systemenergy , Atommass , &
     &                 Avedispl(NDF,*)
      INTEGER Id(NDF,*) , Ix(NEN1,*)
      LOGICAL Moveatoms , Movedisl , Fullfield , Solvefem , changetime
      INTEGER iatom , j , i , Ifem
      DOUBLE PRECISION displ_old(3) , displ_new(3)
      DOUBLE PRECISION Straine0
      integer :: logic
      character*80 :: filename

!!$      filename = 'out/atom_temp_fem0.cfg'
!!$      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
!!$      CALL DUMP_ATOM(Atomcoord,Atomdispl,logic)
!!$      CLOSE (logic)

      CALL CPU_TIME(CT2)
      IF ( Solvefem==.TRUE. ) THEN
!!	   Solve FEM
         CALL VAFUNCMD(Id,Atomcoord,Ix,F,Avedispl,Atomforce,&
     &                 Systemenergy,Moveatoms,Movedisl,Fullfield,&
     &                 Straine0,Ifem,MOVed)
!!$         filename = 'out/atom_temp_fem1.cfg'
!!$         CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
!!$         CALL DUMP_ATOM(Atomcoord,Atomdispl,logic)
!!$         CLOSE (logic)

         
!!		Get Forces and displacements, specifically on PAD atoms
         DO iatom = 1 , NUMnp
            Atomforce(1:NDF,iatom) = -Atomforce(1:NDF,iatom)
            IF ( ISRelaxed(iatom)==INDexcontinuum .OR. ISRelaxed(iatom) ==INDexpad ) THEN
               Atomdispl(1:NDF,iatom) = Avedispl(1:NDF,iatom)
!                if (isrelaxed(iatom) == INDexpad) then
!                   if (abs(atomcoord(2,iatom)) < 5.0) then
!                      if (atomcoord(1,iatom) < 0) then
!                         write(*, '(A25, I7, 4(1X,E15.8))'),'Pad atom displacement = ', iAtom, &
!                              atomcoord(1:2,iatom), atomdispl(1:2, iatom)
!                      end if
!                   end if
!                end if
            END IF              
         ENDDO

!!$         filename = 'out/atom_temp_fem2.cfg'
!!$         CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
!!$         CALL DUMP_ATOM(Atomcoord,Atomdispl,logic)
!!$         CLOSE (logic)

         
         !!       Now set atom displacements
         !---- Keeping in mind that FE does not affect these displacements
         ! --- Alternatively we can just make lammps do the job for these atoms
         ! ---- This is the strategy to follow i think.
         ! --- 1) Create a variable in the lammps module mod_lammps.f90 --- flag to indicate a fixed atom
         ! ----2) Change the x-coordinate of this either by multiplying it here by fem time
         !  --------------   or   -----------------
         !        give a move command in lammps that is independent of the fem time scale etc
         ! ---- if the second choice is made then


!!	   Zero out forces for fixed nodes
         DO i = 1 , NDF
            DO j = 1 , NUMnp
               IF ( IDTemp(i,j) ) Atomforce(i,j) = 0
            ENDDO
         ENDDO
      ENDIF

!!	Zero out forces on all MD atoms in the atomistic region
      !	call InitialiseEnergy(.true.,.true.,id,atomForce,f)
      ! This is not really necesary since lammps will take care of it...
      ! IT is left for testing purposes only
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom) ==INDexinterface ) THEN
            Atomforce(1:3,iatom) = 0.D0
         END IF
      ENDDO


    END SUBROUTINE GETFEM_FORCES






      SUBROUTINE DOSTEPS_lammps(N,Atomdispl,Atomforce,Db,Fl,Epps,Iprint,Dsmax,&
     &                   Rseed,Dfn,Atomid,Atomcoord,Ix,F,Itx,Checkslip,&
     &                   Addedslip,Lostslip,Movedisl,Moveatoms, lmp)



!

      USE MOD_GLOBAL
      USE MOD_POTEN
      USE MOD_DYNAMO
      USE MOD_OUTPUT
      USE MOD_MATERIAL
      USE MOD_PARALLEL
      USE MOD_TIMMING
      USE MOD_COMMON
      USE MOD_FILE
      USE MOD_DISL_PARAMETERS

      USE LAMMPS
      USE MOD_LAMMPS
      IMPLICIT NONE
!*--DOSTEPS1325
      INTEGER MAXATOMS
      PARAMETER (MAXATOMS=100000)
!
      INTEGER N , Ix(*) , Atomid(NDF,*) , Iprint , Itx(*) , Rseed , &
     &        numdis , numh
      DOUBLE PRECISION Atomdispl(NDF,NUMnp) , Db(N) , Atomforce(3,*) , &
     &                 Fl , Epps , Dsmax , Dfn , Atomcoord(NDF,*) , &
     &                 F(*) , mass
      LOGICAL convonfn , Addedslip , Checkslip , Movedisl , Lostslip , &
     &        Moveatoms , DISLCHECK
      LOGICAL dislpass
      LOGICAL printing , DEBug , plot , fullfield , solvefem , &
     &        restartaveraging
      COMMON /DEBUGGER/ DEBug
      COMMON /MD    / OLDacceleration , OLDvelocity , ACCeleration , &
     &                VELocity , AVEdispl , RAMplevel , NUMmdatoms
      ! Lammps variables
      type(c_ptr) :: lmp

!
!--	Local variables
      INTEGER isdamped(MAXATOMS)
      INTEGER nsteps , i , iatom , j , rc , intatom1 , intatom2 , &
     &        intatom3 , intatom4
      INTEGER ii , jj
      INTEGER NUMmdatoms , femsteps , ndis_checked
      INTEGER mdsteps , nnsteps , maxmdsteps
      INTEGER femstepmin , femsteprange , femstepcounter , istep , &
     &        readneighbors
!--	Indexes for H atoms, using the fact that they are adjacent
      INTEGER :: hindex_init = 0 , hindex_final = -1 , indexmin , &
     &           indexmax
      INTEGER :: numinterstitial = 0

      DOUBLE PRECISION atommass , damped_width , langevincoeff , &
     &                 lcnhdampcoeff , requiredtemp , currenttemp , &
     &                 picosecond , systemenergy , timestep , &
     &                 tinterior , stadiumtemp , atomtemp
!
      DOUBLE PRECISION OLDacceleration(3,MAXATOMS) , &
     &                 OLDvelocity(3,MAXATOMS) , &
     &                 ACCeleration(3,MAXATOMS) , VELocity(3,MAXATOMS) ,&
     &                 AVEdispl(3,MAXATOMS) , RAMplevel(MAXATOMS)

      DOUBLE PRECISION xmin , xmax , ymin , ymax , plottime
      CHARACTER*80 atomfilename , energyfilename , tempfilename , &
     &             neighborfilename

      LOGICAL newmd , finishmd

        !!! JS comment For output Avgstress !!!
      DOUBLE PRECISION tempavgstress(3,3)
      INTEGER ifem

      CHARACTER*80 filename
      INTEGER logic , nstepsorig , npass
      DOUBLE PRECISION :: dtol
!
      TYPE (REGION) simulationcell
      TYPE (MD_THERMOSTAT) thermostat
!

!--	Functions
      INTEGER GETSYSTEMSIZE , FINDMDATOMS
      DOUBLE PRECISION GETTEMPERATURE , random , ZBQLU01 , GETRAMP , straine0

      ! ---- Lammps related local variables
      integer :: total_lammps_steps, fem_call_back_steps, lammps_loop, jstep
      character(1024) :: command_line
      logical :: update_pad, update_all

      DATA mdsteps/0/
      DATA nnsteps/0/
      DATA intatom1/0/
      DATA intatom2/0/
      DATA intatom3/0/
      DATA intatom4/0/
      DATA femstepcounter/0/
      DATA ndis_checked/0/

      npass = 0

      IF ( MAXATOMS<NUMnp ) THEN
         PRINT * , 'increase size of maxAtoms in md.f'
         STOP
      ENDIF
!
      plot = .FALSE.
      restartaveraging = .FALSE.
!
      CALL LOSTSLIPINIT(Lostslip)
!
      picosecond = 1.0D-12      ! s
      atommass = AMAss(1)*picosecond**2         ! eV s^2/A^2
      timestep = 1.0D-15                        ! seconds
      femsteps = 1
      fullfield = .TRUE.
      femstepmin = 5
      femsteprange = 0
!
      dtol = 1.D-3

      PRINT * , 'Entering dosteps'
!
!	Read Data
!C--Jun Song NumMDRescale is # of MD steps
!C--doing temperature rescaling
!C--Reading neighborlist update parameter
      OPEN (UNIT=200,FILE='md.inp',STATUS='old')
      READ (200,*) damped_width, exclude_top, exclude_bot, exclude_left, exclude_right
      READ (200,*) langevincoeff , LVScaleratio
      READ (200,*) lcnhdampcoeff
      READ (200,*) requiredtemp
      READ (200,*) TIMestep1
      READ (200,*) INDextimeh
      READ (200,*) femstepmin
      READ (200,*) nsteps
      READ (200,*) thermostat%TYPE
      READ (200,*) thermostat%DAMPING_MODE
      READ (200,*) NUMmdrescale
      READ (200,*) TWIndow
      READ (200,*) NHRescale , MAXhtemp , HNVeflag
      READ (200,*) maxmdsteps
      CLOSE (200)





!! Jun Song comments: output step, energy and temperature per steps
!!(specified when writing to the file)
      OPEN (5800,FILE='MDlog.CADD',STATUS='unknown')
!! Jun Song comment: Initialization

!	Write Data
      WRITE (*,*) 'Base langevinCoeff: ' , langevincoeff
      WRITE (*,*) 'LangevinCoeff scale ratio' , LVScaleratio
!C--Jun Song: Care with Langevin Coefficient!!
      WRITE (*,*) 'Langevin Coefficient' , langevincoeff*LVScaleratio
      WRITE (6,99002) 'requiredTemp: ' , requiredtemp
      WRITE (6,99002) 'timestep1: ' , TIMestep1
      WRITE (6,99002) 'timestepH: ' , TIMestep1/INDextimeh
      WRITE (6,99003) 'FEMStepMin: ' , femstepmin
      WRITE (6,99003) 'Nsteps: ' , nsteps
      WRITE (6,99002) 'atomic mass: ' , atommass
      WRITE (6,99004) 'Thermostat: ' , thermostat%TYPE
      WRITE (6,99004) 'Damping Mode: ' , thermostat%DAMPING_MODE
      WRITE (*,*) '# T Rescale MD steps: ' , NUMmdrescale
      WRITE (*,*) 'The windown for T rescale ' , TWIndow
      WRITE (*,*) '**********For H atom only**********'
      WRITE (*,*) 'H stablizer Steps' , NHRescale , 'for T>' , MAXhtemp
      WRITE (*,*) 'Exclude H from Thermostat? ' , HNVeflag

      ! ---- Lammps is run for fem_call_back_steps for a total of total_lammps_steps
      ! ---- so that the main loop is only for lammps_loop
      fem_call_back_steps = femstepmin
      total_lammps_steps = nsteps
      lammps_loop = total_lammps_steps/fem_call_back_steps
      if (lammps_loop < 1) then
         call lammps_close(lmp)
         call error_handler("lammps loop cannot be less than 1 closing lammps and quitting")
      end if


      finishmd = .FALSE.
      newmd = .FALSE.
      timestep = TIMestep1

      SYStemp = requiredtemp
      PRINT * , 'MDSteps: ' , mdsteps
      IF ( mdsteps==0 ) newmd = .TRUE.
!c	if (MDSteps .eq. MAXMDSteps) finishMD = .true.

      IF ( finishmd ) GOTO 99999
                                ! Finish MD simulation

!C--Here comes inilization if newMD
      IF ( newmd ) THEN
      
!        Equilibrate initialized temperatures...
!        Needed separately since otherwise the mapping
!        gets mess up
!         call equilibrate_lammps(lmp)
      
!        Added here in order give a particle an impact
!        velocity after it's temperature has been equilibrated
!         call add_fix_lammps(lmp)
         
         SIMstep = 0
         ALLOCATE (AVEvirst(3,3,NUMnp))
         ALLOCATE (virst(3,3,NUMnp))
         avedispl = 0.0d0
         ! ---- Initialize FEM
         if (nnsteps == 0) then 
	    AVEdispl(1:NDF,1:NUMnp) = Atomdispl(1:NDF,1:NUMnp)
            AVEvirst(1:3,1:3,1:NUMnp) = 0.D0
            virst(1:3,1:3,1:NUMnp) = 0.D0

         end if
         solveFEM = .true.
         ifem = 1
         CALL GETFEM_FORCES(Atomid,Atomcoord,Ix,F,Atomdispl,&
     &                                AVEdispl,Atomforce,atommass,&
     &                                systemenergy,Moveatoms,Movedisl,&
     &                                fullfield,solvefem,straine0,ifem)

!         update_pad = .true.
         update_all = .false.
         
         call update_lammps_coords(AtomCoord, AtomDispl, update_pad, update_all, lmp)
!!$         filename = 'out/atom.cfg'
!!$         CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
!!$         CALL DUMP_ATOM(Atomcoord,Atomdispl,logic)
!!$         CLOSE (logic)

         ! --- Initial Lammps Run for zero steps to initialize everything
!!$         call lammps_command(lmp, "undump 1")
!!$         write(command_line, fmt='(A18,I4,A55)') "dump 1 all custom ", total_lammps_steps, " atom_lmp*.cfg id type x y z c_dx_all[1], c_dx_all[2]"
!!$         call lammps_command(lmp, command_line)
         
!!$         call lammps_command(lmp, "run 0 pre yes post no")
        !    --Initialize data
!!$         DO iatom = 1 , NUMnp
!!$            IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom)==INDexinterface ) THEN
!!$               AVEdispl(1:NDF,iatom) = 0.D0
!!$            ENDIF
!!$         ENDDO
         newmd = .false.
         solveFEM = .false.
         ifem = 0

!!$      equilibrate for 10k steps
         call lammps_command(lmp,'run 10000')

!!$      recalculate fem forces after equilibration
         CALL GETFEM_FORCES(Atomid,Atomcoord,Ix,F,Atomdispl,&
     &                                AVEdispl,Atomforce,atommass,&
     &                                systemenergy,Moveatoms,Movedisl,&
     &                                fullfield,solvefem,straine0,ifem)

         call update_lammps_coords(AtomCoord, AtomDispl, update_pad, update_all, lmp)
      
	  ENDIF
!     end of if new md loop
!C--New MD initialization ends


!	Assign oldVelocity = velocity
      OLDvelocity(1:NDF,1:NUMnp) = VELocity(1:NDF,1:NUMnp)
!	Main Loop
!--	get the min and max of H indexes
      indexmin = NUMnp
      indexmax = 1

      femsteps = femstepmin
      ifem = 0
      dislpass = .FALSE.
      nstepsorig = nsteps

      DO istep = 1, lammps_loop
         ! ----- Run lammps one step at a time
         ! ---- Here we can choose, but temporarily to maintain structure of current code
         ! ---- Lammps is basically run one step at a time
         do jstep = 1, 1
            if (istep < lammps_loop) then
               write(command_line, *) "run ", fem_call_back_steps, " pre yes post no"

            else
               write(command_line, *) "run ", fem_call_back_steps, " pre yes post yes"
            end if
            if (mod(femstepcounter, femsteps) == 0) then
               ifem = ifem + 1
               solveFem = .TRUE.
!!$               call update_from_lammps(AtomDispl, AtomCoord, AtomForce,  AveDispl, Velocity, virst, avevirst, lmp)

               CALL GETFEM_FORCES(Atomid,Atomcoord,Ix,F,Atomdispl,&
     &                                AVEdispl,Atomforce,atommass,&
     &                                systemenergy,Moveatoms,Movedisl,&
     &                                fullfield,solvefem,straine0,ifem)
!               update_pad = .true.
               update_all = .false. 
               call update_lammps_coords(AtomCoord, AtomDispl, update_pad, update_all, lmp)
               call lammps_command(lmp,'run 0 pre yes post no')


               restartAveraging = .true.
               print *, 'ifem = ', ifem
               solveFem = .false.
            else
               solvefem = .false.
            end if

            ! --- Perform Velocity Verlet using LAMMPS
            call lammps_command(lmp, command_line)

            ! ---- Updates Atom Displacements and forces
            call update_from_lammps(AtomDispl, AtomCoord, AtomForce,  AveDispl, Velocity, virst, avevirst, lmp)

!!$           ! ---- Update Virial Stress from lammps and also obtains average virial stress
!!$            call update_virial_stress_from_lammps(N, Virst, AveVirst, lmp)

            femstepcounter = femstepcounter + fem_call_back_steps
!!$         Now check passing of dislocations
!!$            Realistically depending on the nature of the problem and temperature
!!$              the need is to check for dislocations every md step to ensure that
!!$              any fast moving dislocations are captured
!!$              for testing lammps is run for fem_call_back_steps and then dislocation
!!$              passing is checked.
!!$         Provided facility to do this every step as well
            IF ( dislpass ) THEN
               filename = 'out/atom_pass.cfg'
               CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
               CALL DUMP_ATOM(Atomcoord,Atomdispl,AtomForce,logic)
               CLOSE (logic)

               filename = 'out/atom_pass.vtk'
               CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
               CALL DUMP_MESH(Atomcoord,Atomdispl,Ix,logic)
               CLOSE (logic)

            ENDIF
            dislpass = .FALSE.
!!$
!!$         hard coded off JM
!!$            IF ( DISLCHECK(Checkslip,Lostslip,Addedslip,Movedisl,Ix,&
!!$                 &        Atomcoord,Atomdispl,Itx,ISRelaxed,NUMnp,NDF,NXDm,NUMel,&
!!$                 &        NEN1,NEWmesh,plottime,dislpass,npass) ) THEN
               IF ( dislpass ) PRINT * , 'Dislocation removed from atomistics'
!!$	     if (npass > 1) then
!!$		Nsteps = NstepsOrig*2
!!$		npass = 0
!!$	     end if

               WRITE (*,*) 'Disl checked!'
!!$               mdsteps = mdsteps + 1
               WRITE (*,*) 'Disl checked by rank' , RANk
               ndis_checked = ndis_checked + 1
!!$            END IF
            nnsteps = nnsteps + 1
!!$	  if (Moved) then
!!$       Output the new_atom config after moving atom_displacements
!!$	     filename='out/moved_atoms.cfg'
!!$	     call iofile(filename,'formatted  ',logic,.false.)
!!$	     call dump_atom(atomCoord, atomDispl, logic)
!!$	     close(logic)
!!$	     Moved = .false.
!!$	  end if

!!$         JM, hard coded off
            MOVemesh = .false.
            IF ( MOVemesh ) THEN
!!$               if (Moved) then
!!$                  Moved = .false.
!!$               end if
               IF ( istep==1 ) THEN
                  CALL GET_CRACK_TIP(Atomcoord,Atomdispl)
                  IF ( XTIp(1)>0.0D0 ) THEN
                     IF ( ABS(XTIp(1)-XTIp_init(1))>(X_Move_mesh) ) THEN
                        DO i = 1 , 2
                           XTIp_actual(i) = ABS(XTIp(i)-XTIp_init(i))
                           IF ( XTIp(i)>XTIp_init(i) ) THEN
                              X_Tip_dir(i) = 1.0D0
                           ELSE
                              X_Tip_dir(i) = -1.0D0
                           ENDIF
                        ENDDO
                        PRINT * , 'xtip_actual = ' , XTIp_actual , &
                             &                     X_Tip_dir
                        IF ( X_Tip_dir(1)>0.0D0 ) THEN
                           IF ( INT(XTIp_actual(1)/X_Move_mesh)>0 ) THEN
                              CALL MOVE_ATOMISTIC_CRACK(Atomcoord,Ix,&
                                   &                        Atomdispl)
                              MOVed = .TRUE.
                           ENDIF
                        ENDIF
                        !!$		   return
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
!!$            if (ndisl > 4) then
!!$               if (mod(iStep,10) .eq. 0) then
!!$                  filename='out/atom_pass2.cfg'
!!$                  call iofile(filename,'formatted  ',logic,.false.)
!!$                  call dump_atom(atomCoord, atomDispl, logic)
!!$                  close(logic)
!!$               end if
!!$            end if


			!CHECK THAT THIS IS WORKING
			!Check the update to the B array inside dislocation 
!            update_all = .true.
!            update_pad = .true. 
!            call update_lammps_coords(AtomCoord, AtomDispl,update_all, update_pad, lmp)
         ENDDO
         
         mdsteps = mdsteps + fem_call_back_steps


      end do
99002 FORMAT (a16,2x,1pe11.4)
99003 FORMAT (a16,2x,i5)
99004 FORMAT (a16,2x,a16)
99005 FORMAT (2(a6,2x,f11.4,2x))
 
      
99999      END subroutine dosteps_lammps

