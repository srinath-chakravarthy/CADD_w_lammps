!*==mod_output.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      MODULE MOD_OUTPUT
      IMPLICIT NONE
!*--MOD_OUTPUT6
!
      INTEGER :: atomfile = 10 , energyfile = 11 , tempfile = 17
!
      CONTAINS
!
      SUBROUTINE INITOUTPUT(Atomfilename,Energyfilename,Tempfilename)
      IMPLICIT NONE
!*--INITOUTPUT14
      CHARACTER*80 Atomfilename , Energyfilename , Tempfilename
!
      OPEN (UNIT=ATOmfile,FILE=Atomfilename,STATUS='unknown')
      OPEN (UNIT=ENErgyfile,FILE=Energyfilename,STATUS='unknown')
      OPEN (UNIT=TEMpfile,FILE=Tempfilename,STATUS='unknown')
 
!
      END SUBROUTINE INITOUTPUT
 
!
!
      SUBROUTINE WRITEDATA(Steps,Atomcoord,Velocity,Atomdispl,Atommass,&
     &                     Isdamped)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--WRITEDATA30
      DOUBLE PRECISION Atomcoord(NDF,*) , Velocity(NDF,*) , Atommass , &
     &                 Atomdispl(NDF,*)
      DOUBLE PRECISION cellsize , systemenergy , sysenergy75 , temp75 , &
     &                 temp55 , x0 , y0 , x1 , y1 , x , y , e0 , e1 , &
     &                 e2 , size , tempwhole , sysenergy55 , &
     &                 tempintermediate , temp45 , temp35 , temp25 , &
     &                 temp15
 
      INTEGER i0 , i1 , i2 , i3 , i , Steps , Isdamped(*) , &
     &        numatoms200 , numatoms75 , numatoms55 , numatoms , &
     &        numatoms45 , numatoms35 , numatoms25 , numatoms15
 
! 	size = 55
! 	rc = findCornerAtoms(atomCoord, atomID,
!      $		i0, i1, i2, i3, size)
 
!---	Print Data
      cellsize = 200.0
      systemenergy = GETSYSTEMENERGY(Atomcoord,Velocity,Atommass,&
     &               cellsize)
 
      cellsize = 75.0
      sysenergy75 = GETSYSTEMENERGY(Atomcoord,Velocity,Atommass,&
     &              cellsize)
 
      cellsize = 55.0
      sysenergy55 = GETSYSTEMENERGY(Atomcoord,Velocity,Atommass,&
     &              cellsize)
 
      size = 200.0
      tempwhole = GETTEMPREGION(Atomcoord,Isdamped,Velocity,Atommass,&
     &            size,numatoms200)
 
!
      size = 55.0
      temp55 = GETTEMPREGION(Atomcoord,Isdamped,Velocity,Atommass,size,&
     &         numatoms55)
      size = 45.0
      temp45 = GETTEMPREGION(Atomcoord,Isdamped,Velocity,Atommass,size,&
     &         numatoms45)
      size = 35.0
      temp35 = GETTEMPREGION(Atomcoord,Isdamped,Velocity,Atommass,size,&
     &         numatoms35)
      size = 25.0
      temp25 = GETTEMPREGION(Atomcoord,Isdamped,Velocity,Atommass,size,&
     &         numatoms25)
      size = 15.0
      temp15 = GETTEMPREGION(Atomcoord,Isdamped,Velocity,Atommass,size,&
     &         numatoms15)
 
      numatoms = numatoms200 - numatoms55
      tempintermediate = (tempwhole*numatoms200-temp55*numatoms55)&
     &                   /DFLOAT(numatoms)
!	write(*,*)numAtoms55, numAtoms45,numAtoms35,
!     $            numAtoms25,numAtoms15,numAtoms200
!	stop
 
      WRITE (ENErgyfile,99001) Steps , systemenergy , sysenergy75 , &
     &                         sysenergy55
99001 FORMAT ('steps: ',1x,i5,2x,'energy:',3(1x,1pe16.9))
 
      WRITE (TEMpfile,99002) Steps , tempwhole , temp55 , temp45 , &
     &                       temp35 , temp15
99002 FORMAT ('steps: ',1x,i5,1x,'temp: ',1x,5(f12.4,1x))
 
!
!	if (steps .ge. 5000 .and. mod(steps, 1000) .eq. 0) then
!		size = 75.0
!		call writeAtomConfig(atomCoord, atomDispl, size)
!	endif
 
!
      END SUBROUTINE WRITEDATA
!
 
!
      SUBROUTINE WRITEATOMCONFIG(Atomcoord,Atomdispl,Size)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--WRITEATOMCONFIG111
      DOUBLE PRECISION Atomcoord(NDF,*) , Atomdispl(NDF,*)
      DOUBLE PRECISION , POINTER :: initcoord(:,:)
      DOUBLE PRECISION dx , dy , dz , x , y , z , Size
      INTEGER iatom , atomtype
      LOGICAL firsttime
      DATA firsttime/.TRUE./
 
      IF ( firsttime==.TRUE. ) THEN
         ALLOCATE (initcoord(NDF,NUMnp))
 
         DO iatom = 1 , NUMnp
            initcoord(1,iatom) = Atomcoord(1,iatom)
            initcoord(2,iatom) = Atomcoord(2,iatom)
            initcoord(3,iatom) = Atomcoord(3,iatom)
         ENDDO
 
         firsttime = .FALSE.
      ENDIF
 
!
      DO iatom = 1 , NUMnp
         atomtype = ISRelaxed(iatom)
!
         x = Atomcoord(1,iatom)
         y = Atomcoord(2,iatom)
         z = Atomcoord(3,iatom)
 
         dx = Atomdispl(1,iatom)
         dy = Atomdispl(2,iatom)
         dz = Atomdispl(3,iatom)
 
 
         IF ( atomtype/=INDexcontinuum ) THEN
            IF ( DABS(x)<=Size ) THEN
               IF ( DABS(y)<=Size ) THEN
                  WRITE (ATOmfile,99001) x , y , z , dx , dy , dz , &
     &                   ISRelaxed(iatom)
 
99001             FORMAT (3(f16.9,2x),3(1pe16.9,2x),i3)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      WRITE (ATOmfile,*)
 
      END SUBROUTINE WRITEATOMCONFIG
 
 
!
!
!	Finds the energy of an annular rectangular box with dimensions
!	xmin, ymin, xmax, ymax.
      DOUBLE PRECISION FUNCTION GETENERGYREGION(Atomcoord,Velocity,&
     &   Atommass,Xmin,Xmax,Ymin,Ymax)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--GETENERGYREGION168
      DOUBLE PRECISION Atomcoord(NDF,*) , Velocity(NDF,*) , Atommass , &
     &                 ke , pe , v , size , x , y , Xmin , Xmax , Ymin ,&
     &                 Ymax
 
      INTEGER iatom , j , numatoms
 
      ke = 0.0                  !kinetic energy
      pe = 0.0                  !potential energy
      BOLtzmannconst = 8.62906E-5       ! eV/K
      numatoms = 0
 
      DO iatom = 1 , NUMnp
 
         IF ( ISRelaxed(iatom)<1 ) CYCLE
 
         x = Atomcoord(1,iatom)
         y = Atomcoord(2,iatom)
 
         IF ( (DABS(x)>=Xmin) .OR. (DABS(y)>=Ymin) ) THEN
            IF ( (DABS(x)<=Xmax) .AND. (DABS(y)<=Ymax) ) THEN
!
               DO j = 1 , NDF
                  v = Velocity(j,iatom)
                  ke = ke + v*v
               ENDDO
 
!		write(6,9) x, y, xmin, xmax, ymin, ymax
99001          FORMAT (6(2x,f11.4))
               pe = pe + ENErgy(iatom)
               numatoms = numatoms + 1
            ENDIF
         ENDIF
      ENDDO
 
      ke = 0.5*Atommass*ke                      ! gives kinetic energy i
      GETENERGYREGION = pe + ke         ! gives total energy
!	print*, 'numatoms: ', numAtoms, ' size: ', xmin, xmax, ymin, ymax
      END FUNCTION GETENERGYREGION
 
!
!
!
      DOUBLE PRECISION FUNCTION GETTEMPREGION(Atomcoord,Isdamped,&
     &   Velocity,Atommass,Size,Numatoms)
 
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--GETTEMPREGION216
      DOUBLE PRECISION Atomcoord(NDF,*) , Atommass , Velocity(NDF,*)
      DOUBLE PRECISION kineticenergy , v
      DOUBLE PRECISION temperature , Size , x , y
      INTEGER iatom , Numatoms , j , Isdamped(*)
 
 
      kineticenergy = 0.0
      Numatoms = 0
      DO iatom = 1 , NUMnp
 
         IF ( ISRelaxed(iatom)>=1 ) THEN
 
            x = Atomcoord(1,iatom)
            y = Atomcoord(2,iatom)
 
            IF ( ABS(x)<=Size .AND. ABS(y)<=Size ) THEN
 
               Numatoms = Numatoms + 1
               DO j = 1 , NDF
                  v = Velocity(j,iatom)
                  kineticenergy = kineticenergy + v*v
               ENDDO
            ENDIF
         ENDIF
 
      ENDDO
!
 
      kineticenergy = 0.5*kineticenergy*Atommass
      temperature = kineticenergy/(1.5*BOLtzmannconst*Numatoms)
 
      GETTEMPREGION = temperature
      END FUNCTION GETTEMPREGION
 
 
!
      DOUBLE PRECISION FUNCTION GETSYSTEMENERGY(Atomcoord,Velocity,&
     &   Atommass,Size)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--GETSYSTEMENERGY257
      DOUBLE PRECISION Atomcoord(NDF,*) , Velocity(NDF,*) , Atommass , &
     &                 ke , pe , v , Size , x , y
      INTEGER iatom , j , numatoms
 
      ke = 0.0                  !kinetic energy
      pe = 0.0                  !potential energy
      numatoms = 0
 
      DO iatom = 1 , NUMnp
 
         IF ( ISRelaxed(iatom)<1 ) CYCLE
 
         x = Atomcoord(1,iatom)
         y = Atomcoord(2,iatom)
 
         IF ( ABS(x)<=Size .AND. ABS(y)<=Size ) THEN
 
            DO j = 1 , NDF
               v = Velocity(j,iatom)
               ke = ke + v*v
            ENDDO
 
            pe = pe + ENErgy(iatom)
            numatoms = numatoms + 1
         ENDIF
      ENDDO
 
      ke = 0.5*Atommass*ke                      ! gives kinetic energy i
      GETSYSTEMENERGY = pe + ke         ! gives potential energy
      END FUNCTION GETSYSTEMENERGY
!
!
!
      INTEGER FUNCTION FINDCORNERATOMS(Atomcoord,Atomid,I0,I1,I2,I3,&
     &                                 Size)
 
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--FINDCORNERATOMS296
      DOUBLE PRECISION Atomcoord(NDF,*)
      DOUBLE PRECISION xmin , xmax , ymin , ymax , x , y , Size
      INTEGER iatom , j , flag , Atomid(NDF,*) , natoms
      INTEGER I0 , I1 , I2 , I3
      DOUBLE PRECISION d0 , d1 , d2 , d3 , d
      DOUBLE PRECISION DISTANCE
 
 
      d0 = 10.0
      d1 = 10.0
      d2 = 10.0
      d3 = 10.0
 
      natoms = 0
      DO iatom = 1 , NUMnp
!
         IF ( ISRelaxed(iatom)>=1 ) THEN
 
            natoms = natoms + 1
 
            x = Atomcoord(1,iatom)
            y = Atomcoord(2,iatom)
 
            d = DISTANCE(x,y,-Size,-Size)
            IF ( d<d0 ) THEN
               d0 = d
               I0 = iatom
            ENDIF
 
            d = DISTANCE(x,y,Size,-Size)
            IF ( d<d1 ) THEN
               d1 = d
               I1 = iatom
            ENDIF
 
            d = DISTANCE(x,y,-Size,Size)
            IF ( d<d2 ) THEN
               d2 = d
               I2 = iatom
            ENDIF
 
            d = DISTANCE(x,y,Size,Size)
            IF ( d<d3 ) THEN
               d3 = d
               I3 = iatom
            ENDIF
         ENDIF
      ENDDO
 
      FINDCORNERATOMS = 0
! 	print*, 'Corner Atoms:  ', i0, i1, i2, i3
! 	print*, 'distances: ', d0, d1, d2, d3
 
      END FUNCTION FINDCORNERATOMS
!
!
      SUBROUTINE WRITENEIGHBORLIST(Neighborfilename)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--WRITENEIGHBORLIST357
      INTEGER iatom , lowindex , highindex , range , j
      CHARACTER*80 Neighborfilename
!
      OPEN (UNIT=14,FILE=Neighborfilename,STATUS='unknown')
!	Write the list of indices
      WRITE (14,*) NNIndx(0)
      DO iatom = 1 , NUMnp
!		if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
         WRITE (14,*) iatom , NNIndx(iatom)
      ENDDO
 
!	Write the list of neighbors
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
!
!c		lowIndex = nnindx(iAtom - 1)
            lowindex = 0
            highindex = NNIndx(iatom)
            range = highindex - lowindex
!
            DO j = 1 , range
               WRITE (14,*) iatom , NNLst(iatom,j+lowindex)
            ENDDO
            WRITE (14,*)
         ENDIF
      ENDDO
      CLOSE (14)
!
!
      END SUBROUTINE WRITENEIGHBORLIST
!
 
!	Read the neighborlist
      SUBROUTINE READNEIGHBORLIST(Neighborfilename)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--READNEIGHBORLIST395
      INTEGER iatom , lowindex , highindex , range , j , k
      CHARACTER*80 Neighborfilename
!
      OPEN (UNIT=14,FILE=Neighborfilename,STATUS='old')
!
!	Read the list of indices
      READ (14,*) NNIndx(0)
      DO iatom = 1 , NUMnp
!		if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
         READ (14,*) j , NNIndx(iatom)
         IF ( j/=iatom ) THEN
            PRINT * , 'Stopping.. Error in reading index list'
            STOP
         ENDIF
      ENDDO
 
 
 
!	Read the list of neighbors
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)/=INDexcontinuum ) THEN
!
!c		lowIndex = nnindx(iAtom - 1)
            lowindex = 0
            highindex = NNIndx(iatom)
            range = highindex - lowindex
!
 
            DO j = 1 , range
               READ (14,*) k , NNLst(k,j+lowindex)
!
               IF ( k/=iatom ) THEN
                  PRINT * , &
     &                  'Stopping.. Error in reading list of neighbors'
                  PRINT * , k , j , iatom , range
                  STOP
               ENDIF
!
! 		print*, iAtom, k, nnlst(iAtom, j+lowIndex)
            ENDDO
! 		print*
            READ (14,*)
         ENDIF
      ENDDO
      CLOSE (14)
!
      END SUBROUTINE READNEIGHBORLIST
!
!
      SUBROUTINE CLOSEOUTPUT()
      IMPLICIT NONE
!*--CLOSEOUTPUT447
      CLOSE (ATOmfile)
      CLOSE (ENErgyfile)
      CLOSE (TEMpfile)
      END SUBROUTINE CLOSEOUTPUT
!
!
!	subroutine writeRestartFile(atomCoord, velocity, acceleration)
!
      END MODULE MOD_OUTPUT
