!*==getenergiesandforces.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Find forces on atoms
      SUBROUTINE GETENERGIESANDFORCES(Id,Atomcoord,Ix,F,Atomdispl,&
     &                                Avedispl,Atomforce,Atommass,&
     &                                Systemenergy,Moveatoms,Movedisl,&
     &                                Fullfield,Solvefem,Straine0,Ifem)
 
      USE MOD_GLOBAL
      USE MOD_TIMMING
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
!!$      CALL DUMP_ATOM(Atomcoord,Atomdispl,atomforce,logic)
!!$      CLOSE (logic)

 
      CALL CPU_TIME(CT2)
      IF ( Solvefem==.TRUE. ) THEN
!!	   Solve FEM
         CALL VAFUNCMD(Id,Atomcoord,Ix,F,Avedispl,Atomforce,&
     &                 Systemenergy,Moveatoms,Movedisl,Fullfield,&
     &                 Straine0,Ifem,MOVed)
 
!!$      filename = 'out/atom_temp_fem1.cfg'
!!$      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
!!$      CALL DUMP_ATOM(Atomcoord,Atomdispl,atomforce, logic)
!!$      CLOSE (logic)

      !!		Get Forces and displacements, specifically on PAD atoms
      
         DO iatom = 1 , NUMnp
            Atomforce(1:NDF,iatom) = -Atomforce(1:NDF,iatom)
            IF ( ISRelaxed(iatom)==INDexcontinuum .OR. ISRelaxed(iatom) ==INDexpad ) THEN
               Atomdispl(1:NDF,iatom) = Avedispl(1:NDF,iatom)
            END IF
         ENDDO
         !!       Now set atom displacements
         do iatom = i, numnp
            if (isrelaxed(iatom) == 1) then
               do j = 1, ndf
                  if (id(j, iatom) == 1) then
                     Atomdispl(1:ndf, iatom) = avedispl(1:ndf, iatom)
                  end if
               end do
            end if
         end do
!!$      filename = 'out/atom_temp_fem2.cfg'
!!$      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
!!$      CALL DUMP_ATOM(Atomcoord,Atomdispl,atomforce,logic)
!!$      CLOSE (logic)

         
         
!!	   Zero out forces for fixed nodes
         DO i = 1 , NDF
            DO j = 1 , NUMnp
               IF ( IDTemp(i,j) ) Atomforce(i,j) = 0
            ENDDO
         ENDDO
      ENDIF
      CALL CPU_TIME(CT3)
      CT7 = CT7 + CT3 - CT2
!
!!	Zero out forces on all MD atoms in the atomistic region
!	call InitialiseEnergy(.true.,.true.,id,atomForce,f)
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom) ==INDexinterface ) THEN
            Atomforce(1:3,iatom) = 0.D0
         END IF
      ENDDO
!
!!	Find forces from EAM
      CALL CPU_TIME(CT2)
      CALL PROCESSCLUMP(Id,Atomcoord,Ix,F,Atomdispl,Atomforce,.TRUE.,.TRUE.)
      CALL CPU_TIME(CT3)
      CT4 = CT4 + CT3 - CT2
 
      END SUBROUTINE GETENERGIESANDFORCES
!*==velocityverlet.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--	Move atoms according to a velocity verlet scheme
!C--Jun Song Comment: Add new variables and embed
!C--getEnergiesAndForces subroutine
      SUBROUTINE VELOCITYVERLET(Atomcoord,Atomdispl,Atomforce,Atomid,&
     &                          Isdamped,Oldacceleration,Acceleration,&
     &                          Oldvelocity,Velocity,Timestep,Atommass,&
     &                          Langevincoeff,Requiredtemp,Currenttemp,&
     &                          Simulationcell,Ramplevel,Thermostat,Ix,&
     &                          F,Avedispl,Systemenergy,Moveatoms,&
     &                          Movedisl,Fullfield,Solvefem,Straine0,&
     &                          Ifem)
 
      USE MOD_GLOBAL
      USE MOD_GRAIN
      USE MOD_POTEN
      USE MOD_COMMON
      IMPLICIT NONE
!*--VELOCITYVERLET91
!
!	type(graintype), dimension(:), pointer :: grains
!
      INTEGER n , iprint , maxfn , numdis , Atomid(NDF,*)
      DOUBLE PRECISION Atomdispl(NDF,*) , Atomforce(NDF,*) , &
     &                 Atomcoord(NDF,*) , Ramplevel(*) , &
     &                 Oldacceleration(NDF,*) , Acceleration(NDF,*) , &
     &                 Oldvelocity(NDF,*) , Velocity(NDF,*) , Timestep ,&
     &                 Atommass , maxdispl , displ , Langevincoeff , &
     &                 dampcoeff , interiorcoeff , Currenttemp , &
     &                 Requiredtemp , xmin , xmax , ymin , ymax , &
     &                 damped_width , zetadot , zeta
 
!C--Jun Song: new local variables
      DOUBLE PRECISION F(NDF,*) , Systemenergy , Avedispl(NDF,*) , &
     &                 Straine0 , newdispl(NDF) , currentcoord(NDF) , &
     &                 timestephh
      INTEGER Ix(NEN1,*)
      LOGICAL Moveatoms , Movedisl , Fullfield , Solvefem
!C--Mod Ends
 
      LOGICAL convonfn , addedslip , checkslip , lostslip , dislcheck
      LOGICAL printing , DEBug , plot
      CHARACTER(LEN=16) :: damping_mode
      TYPE (REGION) Simulationcell
      TYPE (MD_THERMOSTAT) Thermostat
      COMMON /DEBUGGER/ DEBug
      DATA zeta/0.D0/
!C--initialize zeta, only called once!!
!
!--	Local variables
      INTEGER iatom , j , natoms , Isdamped(*) , location
!C--constrain motion to be 2D
      INTEGER ndf2d
      INTEGER selectthermostat
      DOUBLE PRECISION tatom , tinterior , rampvalue , meansqdispl , &
     &                 position , size , rmsforce , tdampregion , ramp ,&
     &                 stadiumtemp , neighbortemp , tdamp , randforce , &
     &                 ZBQLU01
!--	Functions
      DOUBLE PRECISION GETKINETICTEMP , GETNEARNEIGHBORTEMP , getramp , &
     &                 BERENDSENDAMPCOEFF , gettemperature , &
     &                 NOSEHOOVERCOEFF
 
!C--Jun Song: half the timestep
!C--dampfactor for thermstat, trsfactor for T rescale
!C--store damp coefficients for each atom to use in final int
!C--hack: Htrsfactor used to stablize H initially (maintain T_H<10000)
!C--Global NHrescale controls # of steps do Htrsfactor
!c--Store the intemediate V for final integration for H atom
!C--Store the random force in initial integration for final integration
      DOUBLE PRECISION dthalf , dampfactor(I_Final) , ntdamp(I_Final)
      DOUBLE PRECISION ndampcoeff(I_Final) , trsfactor , deltatemp
      DOUBLE PRECISION :: htrsfactor = 1.0D0
      DOUBLE PRECISION tempoldvelocity(NDF,I_Final)
      DOUBLE PRECISION randfc(NDF,I_Final)
!C--Mod Ends
!
      LOGICAL usenosehoover , uselangevin , usemarder , usenve
      INTEGER Ifem
!
      dampfactor = 0D0
      deltatemp = 0D0
      natoms = 0
      maxdispl = 0.5D0
!CJSTest: Constrain motion to be 2D
      ndf2d = NDF - 1
 
!--	timestep for H atom
      timestephh = TIMestep1/INDextimeh
 
 
      CALL GETBOXTEMP(Atomcoord,Velocity,Atommass,Simulationcell,&
     &                tinterior,stadiumtemp,Currenttemp)
 
      selectthermostat = 0
     	!! Check for NoseHoover
      usenosehoover = .FALSE.
      location = INDEX(Thermostat%TYPE,'NoseHoover')
      IF ( location>0 ) usenosehoover = .TRUE.
      IF ( location>0 ) selectthermostat = 1
 
      location = 0
   	!! Check for Langevin
      uselangevin = .FALSE.
      location = INDEX(Thermostat%TYPE,'Langevin')
      IF ( location>0 ) uselangevin = .TRUE.
      IF ( location>0 ) selectthermostat = 2
 
      location = 0
   	!! Check for Marder
      usemarder = .FALSE.
      location = INDEX(Thermostat%TYPE,'Marder')
      IF ( location>0 ) usemarder = .TRUE.
      IF ( location>0 ) selectthermostat = 3
 
!C--Jun Song: Adding new thermostat-T rescaling
!C--UseRescale is global parameter, set in dosteps
      IF ( USErescale ) selectthermostat = 4
 
      location = 0
   	!! Check for NVESTAT
      usenve = .FALSE.
      location = INDEX(Thermostat%TYPE,'Nvestat')
      IF ( location>0 ) usenve = .TRUE.
      IF ( location>0 ) selectthermostat = 5
 
!C--Mod Ends
 
!       print*, 'Select Thermostat', selectThermostat
!
      damping_mode = Thermostat%damping_mode
 
!cccccccccccccccccccccccccccccJun Songcccccccccccccccccccccccccccccccc
!C--Replace useNoseHoover with selectThermostat.eq.1
!C--Thus do not update coefficient if doing T rescaling
!C--Use timestepHH because zeta is updated when both H and Ni move
!C--and do MD first for H for SimStep=1 to indextimeH-1
      IF ( selectthermostat==1 ) THEN
         zetadot = NOSEHOOVERCOEFF(Requiredtemp,Currenttemp)
         zeta = zeta + zetadot*timestephh
!		print*, 'Nose Hoover'
      ENDIF
 
!C--Jun Song: opening files to output T and disp for H atoms
      IF ( SIMstep==1 ) THEN
         OPEN (UNIT=9191,FILE='TempHydrogen.dat',STATUS='unknown')
         OPEN (UNIT=9192,FILE='DisplHydrogen.dat',STATUS='unknown')
      ENDIF
!C--Comment End
 
      deltatemp = ABS(stadiumtemp-Requiredtemp)
!-	JS: rescale T if exceeds Twindow
      trsfactor = 1.0D0
      IF ( deltatemp>TWIndow )&
     &     trsfactor = DSQRT(Requiredtemp/Currenttemp)
!
 
!C--Jun Song: Initial integration, update V1 and Disp
      DO iatom = I_Initial , I_Final
         IF ( ATOmspecie(iatom)==2 ) THEN
            Timestep = timestephh
         ELSE
            Timestep = TIMestep1
         ENDIF
!
!-	used to store dampfactor for second integral
         dampfactor(iatom) = 0.0D0
         ndampcoeff(iatom) = 0.0D0
!
         dthalf = 0.5*Timestep
!
!-	Skip near continuum nodes and pad atoms
         Atommass = AMAss(ATOmspecie(iatom))*1.0D-24
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
 
               neighbortemp = GETNEARNEIGHBORTEMP(iatom,Atomcoord,&
     &                        Atomdispl,Isdamped,Oldvelocity,Atommass)
 
!--	output T and Disp for all H atoms
               IF ( MOD(SIMstep,10)==1 ) THEN
                  IF ( ATOmspecie(iatom)==2 ) THEN
!!$                     if (AtomID(2,iatom) == 1) then 
                     tatom = GETKINETICTEMP(Atomcoord,Oldvelocity,iatom,&
     &                       Atommass)
                     WRITE (9191,*) Atommass/1.0365*1.0D+28 , tatom , &
     &                              iatom
                     WRITE (9192,'(3f16.11,i8)') Atomdispl(1,iatom) , &
     &                      Atomdispl(2,iatom) , Atomdispl(3,iatom) , &
     &                      iatom
!!$                     endif
 
                  ENDIF
               ENDIF
 
!-	  Find the damping coeff
               dampcoeff = 0.D0
               IF ( Isdamped(iatom)==0 ) THEN   ! Interior atom
                  dampcoeff = 0.D0
               ELSE
                  tatom = GETKINETICTEMP(Atomcoord,Oldvelocity,iatom,&
     &                    Atommass)
 
!-	Initialize Tdamp to be stadiumTemp at beginning
                  tdamp = stadiumtemp
 
!- 	Update Tdamp and coefficient only when not rescaling
                  IF ( selectthermostat/=4 ) THEN
                     IF ( usenosehoover ) THEN
                        dampcoeff = zeta
                     ELSE
                        rampvalue = Ramplevel(iatom)
                        location = INDEX(damping_mode,'Neighbor')
                        IF ( location>0 ) tdamp = neighbortemp
!
                        location = INDEX(damping_mode,'Stadium')
                        IF ( location>0 ) tdamp = stadiumtemp
!
                        location = INDEX(damping_mode,'Atom')
                        IF ( location>0 ) tdamp = tatom
 
                        location = INDEX(damping_mode,'OFF')
                        IF ( location>0 ) tdamp = Requiredtemp
!-	damp parameter for Marder thermostat
                        dampcoeff = BERENDSENDAMPCOEFF(Requiredtemp,&
     &                              tdamp,Langevincoeff,rampvalue)
                     ENDIF
 
                  ENDIF
               ENDIF
 
 
!C     Initial newdispl vector
               DO j = 1 , NDF
                  newdispl(j) = 0.0D0
!c	     currentCoord(j)=0.0d0
               ENDDO
 
 
!C--Updating for different thermostat
               IF ( selectthermostat==1 ) THEN                         !
                  dampcoeff = zeta
                  ndampcoeff(iatom) = dampcoeff
                  dampfactor(iatom) = EXP(-dthalf*dampcoeff)
!-	  Repeat for each dimension
                  DO j = 1 , ndf2d
!-           Intemediate velocity (same as lammps nvt)
                     Velocity(j,iatom) = Oldvelocity(j,iatom)&
     &                  *dampfactor(iatom) + dthalf*Atomforce(j,iatom)&
     &                  /Atommass
!-	     update displacement
                     newdispl(j) = Velocity(j,iatom)*Timestep
                  ENDDO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c--JS: override to do NVE update for H atom if HNVEFlag==1
                  IF ( (ATOmspecie(iatom)==2) .AND. (HNVeflag==1) ) THEN
                     DO j = 1 , ndf2d
                        Velocity(j,iatom) = Oldvelocity(j,iatom)&
     &                     + dthalf*Atomforce(j,iatom)/Atommass
!c--Store the intemediate V for final integration
                        tempoldvelocity(j,iatom) = Velocity(j,iatom)
                        newdispl(j) = Velocity(j,iatom)*Timestep
 
                     ENDDO
                  ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               ENDIF
 
               IF ( selectthermostat==3 ) THEN
!-	store Tdamp for final integration
                  ntdamp(iatom) = tdamp
 
                  DO j = 1 , ndf2d
!C--Intemediate V (same as lammps nvt)
                     Velocity(j,iatom) = Oldvelocity(j,iatom)&
     &                  + +dthalf*(Atomforce(j,iatom)+DAMpforce(j,iatom)&
     &                  )/Atommass
!c--Updating Disp
                     newdispl(j) = Velocity(j,iatom)*Timestep
                  ENDDO
!
               ENDIF
!
 
               IF ( selectthermostat==2 ) THEN
!C--JS: Real LangevinCoeff!
                  dampcoeff = Langevincoeff*LVScaleratio
                  ndampcoeff(iatom) = dampcoeff
!-	store Tdamp for final integration
                  ntdamp(iatom) = tdamp
 
                  DO j = 1 , ndf2d
!C--Intemediate V
                     Velocity(j,iatom) = Oldvelocity(j,iatom)&
     &                  + +dthalf*(Atomforce(j,iatom)+DAMpforce(j,iatom)&
     &                  )/Atommass
!c--Updating Disp
                     newdispl(j) = Velocity(j,iatom)*Timestep
                  ENDDO
 
               ENDIF
!
!
               IF ( selectthermostat==4 ) THEN                          
                  IF ( stadiumtemp<1E-10 ) THEN
                     WRITE (*,*) 'No rescale with T=0! Exit...'
                     STOP
                  ENDIF
!-	rescale T is exceeds window
                  DO j = 1 , ndf2d
                     IF ( deltatemp>TWIndow ) Velocity(j,iatom)&
     &                    = Oldvelocity(j,iatom)*trsfactor
!
 
                     Velocity(j,iatom) = Velocity(j,iatom)&
     &                  + dthalf*Atomforce(j,iatom)/Atommass
                     newdispl(j) = Velocity(j,iatom)*Timestep
                  ENDDO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c--JS: override to do NVE update for H atom if HNVEFlag==1
                  IF ( (ATOmspecie(iatom)==2) .AND. (HNVeflag==1) ) THEN
                     DO j = 1 , ndf2d
                        Velocity(j,iatom) = Oldvelocity(j,iatom)&
     &                     + dthalf*Atomforce(j,iatom)/Atommass
!c--Store the intemediate V for final integration
                        tempoldvelocity(j,iatom) = Velocity(j,iatom)
                        newdispl(j) = Velocity(j,iatom)*Timestep
                     ENDDO
                  ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               ENDIF
!
               IF ( selectthermostat==5 ) THEN                          
!-	NVE velocity & displacement update
                  DO j = 1 , ndf2d
                     Velocity(j,iatom) = Oldvelocity(j,iatom)&
     &                  + dthalf*Atomforce(j,iatom)/Atommass
                     newdispl(j) = Velocity(j,iatom)*Timestep
                  ENDDO
               ENDIF
 
 
!cccccccccccccccccccccc*H-Stablizer*cccccccccccccccccccccccc
!C--Jun Song: hack for H atom
!c--avoiding Temperature of H being out of control
               IF ( SIMstep<NHRescale ) THEN
                  IF ( ATOmspecie(iatom)==2 ) THEN
                     tatom = GETKINETICTEMP(Atomcoord,Velocity,iatom,&
     &                       Atommass)
                     IF ( tatom>MAXhtemp ) THEN
!C--	The scaling factor of H atom
                        WRITE (*,*) "******H Stablizer******"
                        htrsfactor = DSQRT(Requiredtemp/tatom)
                        DO j = 1 , ndf2d
                           Velocity(j,iatom) = Velocity(j,iatom)&
     &                        *htrsfactor
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
!cccccccccccccccccccccc*H-Stablizer*ccccccccccccccccccccccccc
!
!-	JS: update atomDisp
               DO j = 1 , ndf2d
                  IF ( newdispl(j)>maxdispl ) THEN
!-	JS: warning message
                     WRITE (*,*) 'Warning! DeltaDispl bigger than maxDispl!'
                     newdispl(j) = maxdispl
                  ELSEIF ( newdispl(j)<-maxdispl ) THEN
!-	JS: warning message
                     WRITE (*,*) 'Warning! DeltaDispl less than -maxDispl!'
                     newdispl(j) = -maxdispl
                  ENDIF
!--	    if iAtom is not allowed to move in the j'th direction
                  IF ( Atomid(j,iatom)/=1 ) THEN
 
!--	    Increment atom displacements
                     Atomdispl(j,iatom) = Atomdispl(j,iatom) + newdispl(j)
!
                     rmsforce = rmsforce + Atomforce(j,iatom)**2
                     meansqdispl = meansqdispl + newdispl(j)**2
!!$                  ELSE
!!$                     PRINT *, 'md atom disp = ', iatom, j, atomdispl(j, iatom)
                  ENDIF
!c	    currentCoord(j)=atomCoord(j,iAtom)+atomDispl(j,iAtom)
!
               ENDDO
 
!--Jun Song: put atoms in the simulation box
!c	    if(currentCoord(3) .gt. z_length) then
!c	    	atomDispl(3, iAtom)=atomDispl(3,iAtom)-z_length
!c	    endif
!c	    if(currentCoord(3). lt. 0.0d0) then
!c	    	atomDispl(3, iAtom)=atomDispl(3,iAtom)+z_length
!c	    endif
!
               natoms = natoms + 1
            ENDIF
         ENDIF
      ENDDO
 
 
 
!-	JS: update E and F
      CALL GETENERGIESANDFORCES(Atomid,Atomcoord,Ix,F,Atomdispl,&
     &                          Avedispl,Atomforce,Atommass,&
     &                          Systemenergy,Moveatoms,Movedisl,&
     &                          Fullfield,Solvefem,Straine0,Ifem)
 
 
!-    JS: Final integration, update V2
      DO iatom = I_Initial , I_Final
         IF ( ATOmspecie(iatom)==2 ) THEN
            Timestep = timestephh
         ELSE
            Timestep = TIMestep1
         ENDIF
 
!-    Again, need to use the right timestep for different specie
         dthalf = 0.5*Timestep
!
!-	  Skip near continuum nodes and pad atoms
         Atommass = AMAss(ATOmspecie(iatom))*1.0D-24
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
 
!C--Updating V2 for different thermostat
               IF ( selectthermostat==1 ) THEN                         !
                  DO j = 1 , ndf2d
                     Velocity(j,iatom) = Velocity(j,iatom)&
     &                  *dampfactor(iatom) + dthalf*dampfactor(iatom)&
     &                  *Atomforce(j,iatom)/Atommass
                  ENDDO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c--JS: override to do NVE update for H atom if HNVEFlag==1
                  IF ( (ATOmspecie(iatom)==2) .AND. (HNVeflag==1) ) THEN
                     DO j = 1 , ndf2d
                        Velocity(j,iatom) = tempoldvelocity(j,iatom)&
     &                     + dthalf*Atomforce(j,iatom)/Atommass
                     ENDDO
                  ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               ENDIF
 
 
               IF ( selectthermostat==3 ) THEN                         !
 
                  IF ( (Isdamped(iatom)>0) ) THEN
 
                     dampcoeff = Langevincoeff*rampvalue
                     DO j = 1 , ndf2d
                        randforce = (-1+2*ZBQLU01(0.0D0))&
     &                              *DSQRT(6*dampcoeff*Atommass*&
     &                              BOLtzmannconst*ntdamp(iatom)&
     &                              /Timestep)
!c--JS: override to do NVE update for H atom if HNVEFlag==1
                        IF ( (ATOmspecie(iatom)==2) .AND. (HNVeflag==1)&
     &                       ) THEN
                           randforce = 0.0D0
                           dampcoeff = 0.0D0
                        ENDIF
!c--JS: update Dampforce for damped atoms
                        DAMpforce(j,iatom) = randforce - &
     &                     dampcoeff*Velocity(j,iatom)*Atommass
                     ENDDO
 
                  ENDIF
!c--JS: final integration: update velocities
                  DO j = 1 , ndf2d
                     Velocity(j,iatom) = Velocity(j,iatom)&
     &                  + dthalf*(Atomforce(j,iatom)+DAMpforce(j,iatom))&
     &                  /Atommass
                  ENDDO
 
               ENDIF
 
 
               IF ( selectthermostat==2 ) THEN
                  dampcoeff = ndampcoeff(iatom)
                  DO j = 1 , ndf2d
                     randforce = (-1+2*ZBQLU01(0.0D0))&
     &                           *DSQRT(6*dampcoeff*Atommass*&
     &                           BOLtzmannconst*ntdamp(iatom)/Timestep)
!c--JS: override to do NVE update for H atom if HNVEFlag==1
                     IF ( (ATOmspecie(iatom)==2) .AND. (HNVeflag==1) )&
     &                    THEN
                        randforce = 0.0D0
                        dampcoeff = 0.0D0
                     ENDIF
!c--JS: update Dampforce for damped atoms
                     DAMpforce(j,iatom) = randforce - &
     &                  dampcoeff*Velocity(j,iatom)*Atommass
!c--JS: final update V
                     Velocity(j,iatom) = Velocity(j,iatom)&
     &                  + dthalf*(Atomforce(j,iatom)+DAMpforce(j,iatom))&
     &                  /Atommass
                  ENDDO
!
               ENDIF
 
!
               IF ( selectthermostat==4 ) THEN                          
                  DO j = 1 , ndf2d
                     Velocity(j,iatom) = Velocity(j,iatom)&
     &                  + dthalf*Atomforce(j,iatom)/Atommass
                  ENDDO
               ENDIF
 
               IF ( selectthermostat==5 ) THEN                          
                  DO j = 1 , ndf2d
                     Velocity(j,iatom) = Velocity(j,iatom)&
     &                  + dthalf*Atomforce(j,iatom)/Atommass
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
 
!
      ENDDO
 
 
      natoms = 0
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)>=1 ) THEN
            DO j = 1 , ndf2d
               Oldacceleration(j,iatom) = Acceleration(j,iatom)
               Oldvelocity(j,iatom) = Velocity(j,iatom)
            ENDDO
            natoms = natoms + 1
         ENDIF
      ENDDO
 
      END SUBROUTINE VELOCITYVERLET
!*==getsystemsize.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!
 
 
      INTEGER FUNCTION GETSYSTEMSIZE(Atomcoord,Atomid,Simulationcell)
 
 
      USE MOD_GLOBAL
      USE MOD_COMMON
      IMPLICIT NONE
!*--GETSYSTEMSIZE616
      DOUBLE PRECISION Atomcoord(NDF,*)
      DOUBLE PRECISION xmin , xmax , ymin , ymax , x , y , size
      INTEGER iatom , j , flag , Atomid(NDF,*) , natoms
      TYPE (REGION) Simulationcell
 
      xmax = -1.0D10
      xmin = 1.0D10
      ymax = xmax
      ymin = xmin
 
      natoms = 0
      DO iatom = 1 , NUMnp
!
         IF ( ISRelaxed(iatom)>=1 ) THEN
 
 
            natoms = natoms + 1
 
            x = Atomcoord(1,iatom)
            y = Atomcoord(2,iatom)
 
            IF ( x>xmax ) THEN
               xmax = x
            ELSEIF ( x<xmin ) THEN
               xmin = x
            ENDIF
 
            IF ( y>ymax ) THEN
               ymax = y
            ELSEIF ( y<ymin ) THEN
               ymin = y
            ENDIF
         ENDIF
 
      ENDDO
 
      GETSYSTEMSIZE = 0
 
      Simulationcell%xmin = xmin
      Simulationcell%xmax = xmax
      Simulationcell%ymin = ymin
      Simulationcell%ymax = ymax
!
      PRINT * , 'nAtoms: ' , natoms
      END FUNCTION GETSYSTEMSIZE
!*==distance.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
 
      DOUBLE PRECISION FUNCTION DISTANCE(X0,Y0,X1,Y1)
      IMPLICIT NONE
!*--DISTANCE669
      DOUBLE PRECISION X0 , Y0 , X1 , Y1
 
      DISTANCE = SQRT((X1-X0)**2+(Y1-Y0)**2)
 
      END FUNCTION DISTANCE
!*==findmdatoms.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!
!
!
      INTEGER FUNCTION FINDMDATOMS(Atomcoord,Isdamped,Simulationcell)
      USE MOD_GLOBAL
      USE MOD_COMMON
      IMPLICIT NONE
!*--FINDMDATOMS684
      DOUBLE PRECISION Atomcoord(NDF,*) , damped_width
      DOUBLE PRECISION xmin , xmax , ymin , ymax
      DOUBLE PRECISION x , y , rampvalue , GETRAMP
      INTEGER iatom , Isdamped(*) , numatoms , numdampedatoms
      TYPE (REGION) Simulationcell
 
      numatoms = 0
      numdampedatoms = 0
      DO iatom = 1 , NUMnp
         Isdamped(iatom) = -1
      ENDDO
 
      DO iatom = 1 , NUMnp
 
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
!
               x = Atomcoord(1,iatom)
               y = Atomcoord(2,iatom)
               rampvalue = GETRAMP(Atomcoord,iatom,Simulationcell)
 
               IF ( rampvalue>=1.D-5 ) THEN
                  Isdamped(iatom) = 1
                  numdampedatoms = numdampedatoms + 1
               ELSE
                  Isdamped(iatom) = 0
               ENDIF
!
               numatoms = numatoms + 1
            ENDIF
         ENDIF
 
      ENDDO
 
      PRINT * , 'Number of Damped Atoms: ' , numdampedatoms
      FINDMDATOMS = numatoms
      END FUNCTION FINDMDATOMS
!*==gettemperature.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      DOUBLE PRECISION FUNCTION GETTEMPERATURE(Atomcoord,Isdamped,&
     &   Velocity,Atommass)
 
      USE MOD_GLOBAL
      USE MOD_POTEN
      IMPLICIT NONE
!*--GETTEMPERATURE731
      DOUBLE PRECISION Atomcoord(NDF,*) , Isdamped(*) , Atommass , &
     &                 Velocity(NDF,*)
      DOUBLE PRECISION kineticenergy , v
      DOUBLE PRECISION temperature
      INTEGER iatom , numatoms , j
 
 
      kineticenergy = 0.D0
      numatoms = 0
      DO iatom = 1 , NUMnp
         Atommass = AMAss(ATOmspecie(iatom))*1.0D-24
         IF ( ISRelaxed(iatom)>=1 ) THEN
 
!		if (isDamped(iAtom) .eq. -1) goto 10
!		if (isDamped(iAtom) .eq. 0) goto 10
 
            numatoms = numatoms + 1
            DO j = 1 , NDF
               v = Velocity(j,iatom)
               kineticenergy = kineticenergy + v*v
! 			print*, iAtom, v , isRelaxed(iAtom)
!
            ENDDO
         ENDIF
 
      ENDDO
!
 
      kineticenergy = 0.5*kineticenergy*Atommass
 
      temperature = kineticenergy/(1.5*BOLtzmannconst*numatoms)
 
!	print*, 'getBathTemperature: ', kineticEnergy, numDampedAtoms
      GETTEMPERATURE = temperature
 
      END FUNCTION GETTEMPERATURE
!*==getboxtemp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!
!
!
      SUBROUTINE GETBOXTEMP(Atomcoord,Velocity,Atommass,Simulationcell,&
     &                      Interiortemp,Stadiumtemp,Currenttemp)
 
 
      USE MOD_GLOBAL
      USE MOD_POTEN
      USE MOD_COMMON
      IMPLICIT NONE
!*--GETBOXTEMP781
      DOUBLE PRECISION Atomcoord(NDF,*) , Atommass , Velocity(NDF,*)
      DOUBLE PRECISION keinterior , kestadium , v , keatom , ketotal
      DOUBLE PRECISION temperature , xmin , xmax , ymin , ymax
      DOUBLE PRECISION damped_width , rampvalue , GETRAMP
      DOUBLE PRECISION Interiortemp , Stadiumtemp , Currenttemp
      INTEGER iatom , numintatoms , numstadiumatoms , j , numatoms
      TYPE (REGION) Simulationcell
 
      keinterior = 0.D0
      kestadium = 0.D0
!
      numintatoms = 0
      numstadiumatoms = 0
      numatoms = 0
!
      DO iatom = 1 , NUMnp
 
         IF ( ISRelaxed(iatom)>=1 ) THEN
            Atommass = AMAss(ATOmspecie(iatom))*1.0D-24
            rampvalue = GETRAMP(Atomcoord,iatom,Simulationcell)
 
            keatom = 0.0
            DO j = 1 , NDF
               v = Velocity(j,iatom)
               keatom = keatom + v*v*Atommass
            ENDDO
 
            IF ( rampvalue>1.D-5 ) THEN
               numstadiumatoms = numstadiumatoms + 1
               kestadium = kestadium + keatom
            ELSE
               numintatoms = numintatoms + 1
               keinterior = keinterior + keatom
            ENDIF
!
            ketotal = ketotal + keatom
            numatoms = numatoms + 1
         ENDIF
!
      ENDDO
!
 
      kestadium = 0.5*kestadium
      keinterior = 0.5*keinterior
      ketotal = 0.5*ketotal
!
      Interiortemp = keinterior/(1.5*BOLtzmannconst*numintatoms)
      Stadiumtemp = kestadium/(1.5*BOLtzmannconst*numstadiumatoms)
      Currenttemp = ketotal/(1.5*BOLtzmannconst*numatoms)
 
      END SUBROUTINE GETBOXTEMP
!*==getnearneighbortemp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct
 
 
!
 
!
!
!	get the temperature of the surrounding region
      DOUBLE PRECISION FUNCTION GETNEARNEIGHBORTEMP(Iatom,Atomcoord,&
     &   Atomdispl,Isdamped,Velocity,Atommass)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      USE MOD_POTEN
      IMPLICIT NONE
!*--GETNEARNEIGHBORTEMP847
 
      DOUBLE PRECISION Atomcoord(NDF,*) , Atommass , Velocity(NDF,*) , &
     &                 Atomdispl(NDF,*)
      INTEGER Isdamped(*) , Iatom
!
!--	Local variables
      DOUBLE PRECISION kineticenergy , v , rcut , weight
      DOUBLE PRECISION x1 , y1 , z1 , x2 , y2 , z2 , dist , value
      INTEGER totneighbors , jatom , natoms , i , j , k
!
!--	Functions
      DOUBLE PRECISION CUTOFFRADIUS
 
 
 
      rcut = CUTOFFRADIUS(1)
      x1 = Atomcoord(1,Iatom)
      y1 = Atomcoord(2,Iatom)
      z1 = Atomcoord(3,Iatom)
 
 
!	Find the kinetic energy of the surrounding atoms
      natoms = 0
      kineticenergy = 0.D0
      weight = 0.D0
      DO k = 1 , NUMneighbors(Iatom)
         jatom = NEIghborlist(k,Iatom)
         if (jatom == 0 .or. jatom > numnp) cycle
         IF ( ISRelaxed(jatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(jatom)/=INDexpad ) THEN
!
! 		if (jAtom .eq. iAtom) then
! C			print*, 'Warning; jAtom = iAtom'
! 		endif
 
               x2 = Atomcoord(1,jatom)
               y2 = Atomcoord(2,jatom)
               z2 = Atomcoord(3,jatom)
               Atommass = AMAss(ATOmspecie(jatom))*1.0D-24
               dist = SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
               DO j = 1 , NDF
                  v = Velocity(j,jatom)
                  value = v*v*Atommass
                  kineticenergy = kineticenergy + value
               ENDDO
               natoms = natoms + 1
            ENDIF
         ENDIF
      ENDDO
!
!	Add the KE of the atom itself
      Atommass = AMAss(ATOmspecie(Iatom))*1.0D-24
      DO j = 1 , NDF
         v = Velocity(j,Iatom)
         kineticenergy = kineticenergy + v*v*Atommass
      ENDDO
      natoms = natoms + 1
!
!
!	Get the temperature
      kineticenergy = 0.5*kineticenergy
      GETNEARNEIGHBORTEMP = kineticenergy/(1.5*BOLtzmannconst*natoms)
 
!      	write(6,9) iAtom, nAtoms
! 9	format('iAtom:', 1x, i5, 2x, 'neighboring atoms:', 1x, i5)
      END FUNCTION GETNEARNEIGHBORTEMP
!*==setatomvelocity.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
!
!
 
!
!
      SUBROUTINE SETATOMVELOCITY(Atomcoord,Isdamped,Velocity,Atommass,&
     &                           Requiredtemp,Rseed, Atomid)
 
      USE MOD_GLOBAL
      USE MOD_PARALLEL
      USE MOD_POTEN
      IMPLICIT NONE
!*--SETATOMVELOCITY927
      DOUBLE PRECISION Atomcoord(NDF,*) , Velocity(NDF,*) , Atommass , &
     &                 kineticenergy , Requiredtemp , v , tmpvalue , &
     &                 dfvelocity , avevelocity(NDF) , totmass, &
     &                 Atomid(NDF,*)
      INTEGER iatom , j , Isdamped(*)
      INTEGER numatoms , Rseed
      DOUBLE PRECISION zbqlu01 , zbqluab , random , GETRANDOMNUMBER
!
      numatoms = 0
 
      DO j = 1 , 3
         avevelocity(j) = 0.D0
      ENDDO
 
      totmass = 0.D0
 
!--	Initialize random number generator
!	call ZBQLINI(2478945)
      CALL ZBQLINI(RANk+Rseed)
      PRINT * , RANk , 'random seed = ' , RANk + Rseed
 
      kineticenergy = BOLtzmannconst*Requiredtemp
!c	v = sqrt(kineticEnergy/atomMass)
!
      DO iatom = 1 , NUMnp
         Atommass = AMAss(ATOmspecie(iatom))*1.0D-24
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
               v = SQRT(kineticenergy/Atommass)
               DO j = 1 , NDF
                  random = GETRANDOMNUMBER()
                  dfvelocity = v*random
                  Velocity(j,iatom) = dfvelocity
                  avevelocity(j) = avevelocity(j) + dfvelocity*Atommass
                  totmass = totmass + Atommass
               ENDDO
!
               numatoms = numatoms + 1
            ENDIF
         ENDIF
 
      ENDDO
!
 
      DO j = 1 , NDF
         avevelocity(j) = avevelocity(j)/totmass
         WRITE (6,99001) j , avevelocity(j)
99001    FORMAT ('Average velocity[',i1,']:',2x,1pe11.4)
      ENDDO
!
 
!	Subract Average velocities
      DO iatom = 1 , NUMnp
!
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
!
               DO j = 1 , NDF
                  Velocity(j,iatom) = Velocity(j,iatom) - avevelocity(j)
               ENDDO
            ENDIF
         ENDIF
!
      ENDDO
 
      PRINT * , 'setAtomVelocity: Total MD Atoms: ' , numatoms
      END SUBROUTINE SETATOMVELOCITY
!*==getrandomnumber.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
 
 
 
 
!
      DOUBLE PRECISION FUNCTION GETRANDOMNUMBER()
      IMPLICIT NONE
!*--GETRANDOMNUMBER1003
      DOUBLE PRECISION random , ZBQLU01
!
      random = -1.D0 + 2.D0*ZBQLU01(0.0D0)
!
      IF ( random>0 ) THEN
         random = 1.D0
      ELSEIF ( random<0 ) THEN
         random = -1.D0
      ELSEIF ( DABS(random)<1.0D-6 ) THEN
         random = 0.D0
      ENDIF
!
      GETRANDOMNUMBER = random
      END FUNCTION GETRANDOMNUMBER
!*==getkinetictemp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
 
!
 
!
!
      DOUBLE PRECISION FUNCTION GETKINETICTEMP(Atomcoord,Velocity,Iatom,&
     &   Atommass)
      USE MOD_GLOBAL
      USE MOD_POTEN
      IMPLICIT NONE
!*--GETKINETICTEMP1030
      DOUBLE PRECISION temp , ke , Atommass , Atomcoord(NDF,*) , &
     &                 Velocity(NDF,*) , v , x , y
      INTEGER j , Iatom
!
 
      ke = 0.D0
      DO j = 1 , NDF
         v = Velocity(j,Iatom)
         ke = ke + v*v
      ENDDO
      Atommass = AMAss(ATOmspecie(Iatom))*1.0D-24
      ke = ke*0.5*Atommass
      temp = ke/(1.5*BOLtzmannconst)
 
      GETKINETICTEMP = temp
!
      END FUNCTION GETKINETICTEMP
!*==getramp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!
!
 
 
      DOUBLE PRECISION FUNCTION GETRAMP(Atomcoord,Iatom,Simulationcell)
      USE MOD_GLOBAL
      USE MOD_COMMON
      IMPLICIT NONE
!*--GETRAMP1058
      DOUBLE PRECISION Atomcoord(NDF,*) , xmin , xmax , ymin , ymax
      INTEGER Iatom
      DOUBLE PRECISION v1 , v2 , v3 , v4 , x , y , damp_width , MINVALUE
      TYPE (REGION) Simulationcell
!
!
      x = Atomcoord(1,Iatom)
      y = Atomcoord(2,Iatom)
 
      v1 = ABS(x-Simulationcell%xmin)
      v2 = ABS(x-Simulationcell%xmax)
      v3 = ABS(y-Simulationcell%ymin)
      v4 = ABS(y-Simulationcell%ymax)
 
 
      GETRAMP = MINVALUE(v1,v2,v3,v4)
      IF ( GETRAMP>Simulationcell%DAMPED_WIDTH ) THEN
         GETRAMP = 0.D0
      ELSE
         GETRAMP = 1.D0 - GETRAMP/Simulationcell%DAMPED_WIDTH
      ENDIF
!	getRamp = 1.0
 
!	print*, 'getRamp: ', getRamp
 
      END FUNCTION GETRAMP
!*==berendsendampcoeff.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 
 
 
!
!
 
      DOUBLE PRECISION FUNCTION BERENDSENDAMPCOEFF(Requiredtemp,&
     &   Kinetictemp,Langevincoeff,Rampvalue)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--BERENDSENDAMPCOEFF1095
      DOUBLE PRECISION Langevincoeff , Requiredtemp , Kinetictemp , &
     &                 damped_width , deltat , t0 , t , Rampvalue
!
 
      deltat = 1.D0/20.D0
      t = Kinetictemp
      t0 = Requiredtemp
!
      BERENDSENDAMPCOEFF = Langevincoeff*(t-t0)/SQRT(t**2+(deltat)**2)&
     &                     *Rampvalue
 
 
      END FUNCTION BERENDSENDAMPCOEFF
!*==nosehoovercoeff.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
!
!
      DOUBLE PRECISION FUNCTION NOSEHOOVERCOEFF(Requiredtemp,&
     &   Kinetictemp)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--NOSEHOOVERCOEFF1116
      DOUBLE PRECISION Requiredtemp , Kinetictemp , zetadot , omegae , &
     &                 q , prefactor , unitpico
      INTEGER numdampedatoms
 
!c	OmegaE = 8.0e14		!Einstein Frequency s^-1
!c	prefactor = 0.005
!c	Q = boltzmannConst * requiredTemp/OmegaE**2
!c	zetaDot = 3.d0/Q * boltzmannConst *(kineticTemp - requiredTemp)
!c 	NoseHooverCoeff = zetaDot * prefactor
 
      unitpico = 1.0E-12
      IF ( NHDampcoeff<0.0D0 ) THEN
         WRITE (*,*) "Nose Hoover dampcoeff less than 0!!!"
         STOP
      ENDIF
 
      NOSEHOOVERCOEFF = (Kinetictemp/Requiredtemp-1.0D0)/unitpico**2
      NOSEHOOVERCOEFF = NOSEHOOVERCOEFF/NHDampcoeff**2
!
      END FUNCTION NOSEHOOVERCOEFF
!*==minvalue.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!
!	minimum of four real numbers, modified to calculate the minimum of fir
!	two
      DOUBLE PRECISION FUNCTION MINVALUE(A,B,C,D)
      IMPLICIT NONE
!*--MINVALUE1144
      DOUBLE PRECISION A , B , C , D , value(4)
      DOUBLE PRECISION min_val
      INTEGER i
 
      value(1) = A
      value(2) = B
      value(3) = C
      value(4) = D
!
      DO i = 1 , 4
         IF ( i==1 ) THEN
            min_val = value(1)
         ELSEIF ( value(i)<min_val ) THEN
            min_val = value(i)
         ENDIF
      ENDDO
 
      MINVALUE = min_val
      END FUNCTION MINVALUE
!*==boxmuller.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!
 
 
!
      SUBROUTINE BOXMULLER(Atomcoord,Isdamped,Velocity,Atommass,&
     &                     Requiredtemp)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--BOXMULLER1175
!
      DOUBLE PRECISION Atomcoord(NDF,*) , Velocity(NDF,*) , &
     &                 avevelocity(NDF) , Atommass , Requiredtemp , v , &
     &                 dfvelocity
      INTEGER iatom , j , Isdamped(*) , numatoms
      DOUBLE PRECISION ZBQLU01 , r1 , r2 , sigma , random , pi
!
!
!
      pi = DACOS(-1.D0)
!
      DO j = 1 , 3
         avevelocity(j) = 0.D0
      ENDDO
 
 
!	Initialize random number generator with a chosen random seed
      CALL ZBQLINI(2478945)
!
!	The variance in velocities
      sigma = SQRT(BOLtzmannconst*Requiredtemp/Atommass)
 
      DO iatom = 1 , NUMnp
 
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
!
               DO j = 1 , NDF
!
!	Get two random numbers
                  r1 = ZBQLU01(0.0D0)
                  r2 = ZBQLU01(0.0D0)
!	Derive the random number used
                  random = SQRT(-2.0*LOG(r1))*COS(2.0*pi*r2)
!
                  dfvelocity = sigma*random
                  Velocity(j,iatom) = dfvelocity
                  avevelocity(j) = avevelocity(j) + dfvelocity
               ENDDO
!
               numatoms = numatoms + 1
            ENDIF
         ENDIF
 
      ENDDO
!
 
      DO j = 1 , NDF
         avevelocity(j) = avevelocity(j)/DFLOAT(numatoms)
         WRITE (6,99001) j , avevelocity(j)
99001    FORMAT ('Average velocity[',i1,']:',2x,1pe11.4)
      ENDDO
 
!
!	Subract Average velocities
      DO iatom = 1 , NUMnp
!
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
            IF ( ISRelaxed(iatom)/=INDexpad ) THEN
!
               DO j = 1 , NDF
                  Velocity(j,iatom) = Velocity(j,iatom) - avevelocity(j)
               ENDDO
            ENDIF
         ENDIF
!
      ENDDO
 
      PRINT * , 'BoxMuller: Total MD Atoms: ' , numatoms
      END SUBROUTINE BOXMULLER
!*==nhstress.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! New subroutine for calculate the average stresses !!!!!
!!!!! Based on the subroutine hstress - Jun Song comment!!!!!
 
      SUBROUTINE NHSTRESS(Tempavgstress,Pnump,Pisrelaxed,Pvirst)
      IMPLICIT NONE
!*--NHSTRESS1255
      INTEGER i , j , k , pnatoms
      INTEGER Pnump
      INTEGER Pisrelaxed(Pnump)
      DOUBLE PRECISION Pvirst(3,3,Pnump)
      DOUBLE PRECISION phydro_stress
      DOUBLE PRECISION pavg_stress(3,3)
      DOUBLE PRECISION ptstress
      DOUBLE PRECISION Tempavgstress(3,3)
 
 
 
      phydro_stress = 0.0
      pavg_stress(1:3,1:3) = 0.0
      pnatoms = 0
      DO i = 1 , Pnump
         IF ( Pisrelaxed(i)==1 ) THEN
            pnatoms = pnatoms + 1
            DO j = 1 , 3
               phydro_stress = phydro_stress + Pvirst(j,j,i)
               DO k = 1 , 3
                  pavg_stress(j,k) = pavg_stress(j,k) + Pvirst(j,k,i)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
 
      pavg_stress(1:3,1:3) = pavg_stress(1:3,1:3)/pnatoms
      phydro_stress = phydro_stress*1.0/3.0/pnatoms
 
      ptstress = 0.5*(pavg_stress(1,1)+pavg_stress(3,3))
 
 
      Tempavgstress(1:3,1:3) = pavg_stress(1:3,1:3)
 
      END SUBROUTINE NHSTRESS
!*==dosteps.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Jun Song comment end			!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The subroutine called by ma06 to do MD steps in cgma05.f!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
! ---- use timeplot for plot view
 
      SUBROUTINE DOSTEPS(N,Atomdispl,Atomforce,Db,Fl,Epps,Iprint,Dsmax,&
     &                   Rseed,Dfn,Atomid,Atomcoord,Ix,F,Itx,Checkslip,&
     &                   Addedslip,Lostslip,Movedisl,Moveatoms)
 
 
 
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
      DOUBLE PRECISION GETTEMPERATURE , random , ZBQLU01 , GETRAMP , &
     &                 straine0
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
      READ (200,*) damped_width
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
      WRITE (6,99002) 'damped_width: ' , damped_width
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
 
      finishmd = .FALSE.
      newmd = .FALSE.
      timestep = TIMestep1
 
!C--Jun Song: set SysTemp and CUTFACT to user defined value
!C--	set Nose Hoover damping coefficient as from input
      NHDampcoeff = lcnhdampcoeff
      SYStemp = requiredtemp
!c--	Stop if CUTFACT value if too small or too big
      IF ( CUTfact<1.0D0 .OR. CUTfact>2.0D0 ) THEN
         WRITE (*,*) "*****************************************"
         WRITE (*,*) "Neigh-cutoff not suitable! Set to default"
         WRITE (*,*) "*****************************************"
         STOP
      ENDIF
 
      PRINT * , 'MDSteps: ' , mdsteps
      IF ( mdsteps==0 ) newmd = .TRUE.
!c	if (MDSteps .eq. MAXMDSteps) finishMD = .true.
 
      USErescale = .FALSE.
      IF ( mdsteps<NUMmdrescale ) USErescale = .TRUE.
!
      IF ( finishmd ) GOTO 99999
                                ! Finish MD simulation
 
!C--Here comes inilization if newMD
      IF ( newmd ) THEN
         SIMstep = 0
         ALLOCATE (AVEvirst(3,3,NUMnp))
 
 
!      Find the atom# of two interface atom along crack line
         DO iatom = 1 , NUMnp
            IF ( ISRelaxed(iatom)==INDexinterface .AND. &
     &           Atomcoord(2,iatom)==0 ) THEN
               IF ( Atomcoord(1,iatom)<0 ) intatom1 = iatom
               IF ( Atomcoord(1,iatom)>0 ) intatom2 = iatom
            ENDIF
            IF ( ISRelaxed(iatom)==INDexinterface .AND. &
     &           Atomcoord(1,iatom)==0 ) THEN
               IF ( Atomcoord(2,iatom)<0 ) intatom3 = iatom
               IF ( Atomcoord(2,iatom)>0 ) intatom4 = iatom
            ENDIF
         ENDDO
 
!         Get System Size information
         rc = GETSYSTEMSIZE(Atomcoord,Atomid,simulationcell)
         simulationcell%damped_width = damped_width
 
         NUMmdatoms = FINDMDATOMS(Atomcoord,isdamped,simulationcell)
 
 
         PRINT * , 'Simulation Cell Dimensions:'
         WRITE (6,99005) 'xmin:' , simulationcell%xmin , 'xmax:' , simulationcell%xmax
         WRITE (6,99005) 'ymin:' , simulationcell%ymin , 'ymax:' , simulationcell%ymax
         PRINT * , 'numMDAtoms: ' , NUMmdatoms
 
 
!    --Initialize data
         DO iatom = 1 , NUMnp
            IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom)==INDexinterface ) THEN
               VELocity(1:NDF,iatom) = 0.D0
               OLDvelocity(1:NDF,iatom) = 0.D0
               OLDacceleration(1:NDF,iatom) = 0.D0
               AVEdispl(1:NDF,iatom) = 0.D0
               RAMplevel(iatom) = 0.D0
               DAMpforce(1:NDF,iatom) = 0.D0
            ENDIF
         ENDDO
 
 
!!	Set atom velocities
!C	Only do this if newMD
         CALL SETATOMVELOCITY(Atomcoord,isdamped,VELocity,atommass,&
     &                        requiredtemp*2.D0,Rseed)
!
!
!!	Get Ramp Value
         DO iatom = 1 , NUMnp
            IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom)&
     &           ==INDexinterface ) RAMplevel(iatom)&
     &           = GETRAMP(Atomcoord,iatom,simulationcell)
         ENDDO
 
      ENDIF
!     end of if new md loop
!C--New MD initialization ends
 
 
!	Assign oldVelocity = velocity
      OLDvelocity(1:NDF,1:NUMnp) = VELocity(1:NDF,1:NUMnp)
      CALL GETBOXTEMP(Atomcoord,VELocity,atommass,simulationcell,&
     &                tinterior,stadiumtemp,currenttemp)
 
      PRINT * , 'Total Temperature:' , currenttemp , 'iT ' , tinterior ,&
     &      'Stemp ' , stadiumtemp
 
!	Main Loop
!--	get the min and max of H indexes
      indexmin = NUMnp
      indexmax = 1
 
      femsteps = femstepmin
      DO iatom = 1 , NUMnp
         IF ( ATOmspecie(iatom)==2 ) THEN
            numh = iatom
            WRITE (*,*) "*********numH*********" , numh
            IF ( iatom>indexmax ) indexmax = iatom
            IF ( iatom<indexmin ) indexmin = iatom
         ENDIF
      ENDDO
 
 
!--	when there is H atom. Set HIndex_Init and Final
      IF ( indexmax>=indexmin ) THEN
         hindex_init = indexmin
         hindex_final = indexmax
         numinterstitial = indexmax - indexmin + 1
!	Output the min and max of H indexes
         WRITE (*,*) "******HIndex_Init is " , hindex_init
         WRITE (*,*) "******HIndex_Final is " , hindex_final
      ENDIF
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C--Jun Song: No multiple timestep if no interstitial!
      IF ( (numinterstitial<=0) .AND. (INDextimeh>1) ) THEN
         WRITE (*,*) "Only use mult-timestep with >0 intersitial!"
         STOP
      ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      ifem = 0
      dislpass = .FALSE.
      nstepsorig = nsteps
      DO istep = 0 , nsteps
!C--Increment SimStep for each simulation step
         SIMstep = SIMstep + 1
!$$$	    if (iStep == Nsteps-f) then
!$$$	       MoveDisl = .true.
!$$$	    else
!$$$	       MoveDisl = .false.
!$$$	    endif
         IF ( MOD(SIMstep,INDextimeh)==0 ) THEN
!	    write(*,*)"iStep_: ",iStep
            I_Initial = 1
            I_Final = NUMnp
         ELSE
            I_Initial = hindex_init
            I_Final = hindex_final
         ENDIF
 
 
         IF ( MOD(istep,101)==0 ) WRITE (*,*) "Simlation Step: " , istep
 
 
         IF ( nnsteps==0 ) THEN
!	!! At the very first step, aveDispl = atomDispl for all nodes
            AVEdispl(1:NDF,1:NUMnp) = Atomdispl(1:NDF,1:NUMnp)
            AVEvirst(1:3,1:3,1:NUMnp) = 0.D0
         ENDIF
!
!   !!Find Average positions of free and interface atoms
         IF ( MOD(femstepcounter,femsteps)==0 ) THEN
            ifem = ifem + 1
            solvefem = .TRUE.
            restartaveraging = .TRUE.
            PRINT * , 'Ifem = ' , ifem
         ELSE
!$$$	     if (dislpass) then
!$$$		solveFEM = .true.
!$$$		restartAveraging = .true.
!$$$		iFem=1
!$$$	     else
            solvefem = .FALSE.
!$$$	     end if
         ENDIF
!$$$	  if (iStep == NSteps-FEM) then
!$$$	     solveFEM = .true.
!$$$	  endif
          !! Get Energies and Forces on MD atoms
!C--Jun Song: Get forces and energies for first runs
         IF ( newmd ) THEN
            CALL GETENERGIESANDFORCES(Atomid,Atomcoord,Ix,F,Atomdispl,&
     &                                AVEdispl,Atomforce,atommass,&
     &                                systemenergy,Moveatoms,Movedisl,&
     &                                fullfield,solvefem,straine0,ifem)
            newmd = .FALSE.
         ENDIF
         filename = 'out/atom_temp_anewmd.cfg'
         CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
         CALL DUMP_ATOM(Atomcoord,Atomdispl,atomforce,logic)
         CLOSE (logic)

         
!         !!Get Current temperature of the MD region
         currenttemp = GETTEMPERATURE(Atomcoord,isdamped,VELocity,&
     &                 atommass)
 
!         !! Increment atom positions
         CALL VELOCITYVERLET(Atomcoord,Atomdispl,Atomforce,Atomid,&
     &                       isdamped,OLDacceleration,ACCeleration,&
     &                       OLDvelocity,VELocity,timestep,atommass,&
     &                       langevincoeff,requiredtemp,currenttemp,&
     &                       simulationcell,RAMplevel,thermostat,Ix,F,&
     &                       AVEdispl,systemenergy,Moveatoms,Movedisl,&
     &                       fullfield,solvefem,straine0,ifem)

            write(filename,fmt='(A14,I5,A4)') 'out/atom_temp_', istep,'.cfg'
!!$	    filename = 'out/atom_temp.cfg'
            CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
            CALL DUMP_ATOM(Atomcoord,Atomdispl,atomforce,logic)
            CLOSE (logic)

          !! Restart averaging displacements for MD atoms
         IF ( restartaveraging==.TRUE. ) THEN
            DO iatom = 1 , NUMnp
               IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom) ==INDexinterface ) THEN
                  AVEdispl(1:NDF,iatom) = 0.D0
               END IF
            ENDDO
!
            random = ZBQLU01(0.0D0)
            femsteps = femstepmin + INT(random*femsteprange)
            IF ( femsteps>nsteps ) femsteps = nsteps + 1
!      !print*, 'FEMSteps: ', FEMSteps
            femstepcounter = 0
            restartaveraging = .FALSE.
         ENDIF
 
!  !!Calculate average displacements for interface and free atoms
         DO iatom = 1 , NUMnp
            IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom) ==INDexinterface ) THEN
               AVEdispl(1:NDF,iatom) = AVEdispl(1:NDF,iatom) + Atomdispl(1:NDF,iatom)/femsteps
            END IF
         ENDDO
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C--Jun Song Test
!c	  do iAtom=1, numnp
!c	     if(atomSpecie(iAtom).eq.2) then
!c	       write(*,*) "# of Steps and H neighbors", SimStep,
!c     $	       NumNeighbors(iAtom)
!c	     endif
!c	  enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
          !! calculate average virial stresses
         IF ( istep>=nsteps/2 ) THEN
            DO iatom = 1 , NUMnp
               atommass = AMAss(ATOmspecie(iatom))*1.0D-24
               IF ( ISRelaxed(iatom)==INDexatom .OR. ISRelaxed(iatom)==INDexinterface ) THEN
                 !! add kinetic term
                  DO ii = 1 , 3
                     DO jj = 1 , 3
                        VIRst(ii,jj,iatom) = VIRst(ii,jj,iatom)&
     &                     - VELocity(ii,iatom)*VELocity(jj,iatom)&
     &                     *atommass/((MATerial(1)%A0**3)/4.D0)
                     ENDDO
                  ENDDO
 
                  AVEvirst(1:3,1:3,iatom)&
     &               = (AVEvirst(1:3,1:3,iatom)*(istep-nsteps/2)&
     &               +VIRst(1:3,1:3,iatom))/(istep-nsteps/2+1)
               ENDIF
            ENDDO
         ENDIF
 
         femstepcounter = femstepcounter + 1
!	  print *, 'FEMSTepCounter =', FEMStepCounter
          !!Check for emitted Dislocations
         plottime = nnsteps*TIMestep1*1.D12
 
         CALL CPU_TIME(CT2)
!	  print *, 'Calling dislcheck'
         IF ( dislpass ) THEN
            filename = 'out/atom_pass.cfg'
            CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
            CALL DUMP_ATOM(Atomcoord,Atomdispl,atomforce,logic)
            CLOSE (logic)
 
            filename = 'out/atom_pass.vtk'
            CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
            CALL DUMP_MESH(Atomcoord,Atomdispl,Ix,logic)
            CLOSE (logic)
 
         ENDIF
         dislpass = .FALSE.
         IF ( DISLCHECK(Checkslip,Lostslip,Addedslip,Movedisl,Ix,&
     &        Atomcoord,Atomdispl,Itx,ISRelaxed,NUMnp,NDF,NXDm,NUMel,&
     &        NEN1,NEWmesh,plottime,dislpass,npass) ) THEN
            IF ( dislpass ) PRINT * , &
     &                            'Dislocation removed from atomistics'
!$$$	     if (npass > 1) then
!$$$		Nsteps = NstepsOrig*2
!$$$		npass = 0
!$$$	     end if
 
            WRITE (*,*) 'time = ' , plottime , 'ps'
            WRITE (*,*) 'Disl checked!'
            mdsteps = mdsteps + 1
            WRITE (*,*) 'Disl checked by rank' , RANk
            ndis_checked = ndis_checked + 1
!             if (ndis_checked.eq.1) then
!               call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
!               stop
!                goto 500
!             endif
         ENDIF
         CALL CPU_TIME(CT3)
         CT5 = CT5 + CT3 - CT2
 
         CALL GETBOXTEMP(Atomcoord,VELocity,atommass,simulationcell,&
     &                   tinterior,stadiumtemp,currenttemp)
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C--Jun Song: output step, energy, temperature, stress per steps
         tempavgstress(1:3,1:3) = 0.0
         IF ( MOD(istep,INDextimeh)==0 ) THEN
!! Jun Song:  get the average stress of the system
!C	     call Nhstress(tempAvgStress,numnp,IsRelaxed,avevirst)
!! Can output average stresses if needed
            WRITE (5800,99001) SIMstep , currenttemp , stadiumtemp , &
     &                         tinterior , TOTal_energymd , ENErgyh
 
99001       FORMAT (i8,3F10.3,2F12.4)
         ENDIF
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
         nnsteps = nnsteps + 1
!$$$	  if (Moved) then
!$$$c       Output the new_atom config after moving atom_displacements
!$$$	     filename='out/moved_atoms.cfg'
!$$$	     call iofile(filename,'formatted  ',logic,.false.)
!$$$	     call dump_atom(atomCoord, atomDispl, logic)
!$$$	     close(logic)
!$$$	     Moved = .false.
!$$$	  end if
 
         IF ( MOVemesh ) THEN
!$$$	     if (Moved) then
!$$$		Moved = .false.
!$$$	     end if
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
!$$$		   return
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         write(filename, '(a,i6.6)') 'dump_inter', simstep
         call iofile(filename,'formatted  ',logic,.false.)
         call dump_atom(atomcoord, atomdispl, atomforce, logic)
         close(logic)

!$$$	  if (ndisl > 4) then
!$$$	     if (mod(iStep,10) .eq. 0) then
!$$$		filename='out/atom_pass2.cfg'
!$$$		call iofile(filename,'formatted  ',logic,.false.)
!$$$		call dump_atom(atomCoord, atomDispl, logic)
!$$$		close(logic)
!$$$	     end if
!$$$	  end if
      ENDDO
 
      mdsteps = mdsteps + 1
 
      PRINT * , 'total_energyMD' , TOTal_energymd
      IF ( finishmd==.FALSE. ) RETURN
!
!! Finish MD
      STOP
99002 FORMAT (a16,2x,1pe11.4)
99003 FORMAT (a16,2x,i5)
99004 FORMAT (a16,2x,a16)
99005 FORMAT (2(a6,2x,f11.4,2x))
 
 
 
99999 END SUBROUTINE DOSTEPS
 
 
