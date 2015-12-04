!*==vafuncmd.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE VAFUNCMD(Id,X,Ix,F,B,Dr,Totener2,Moveatoms,Movedisl,&
     &                    Fullfield,Straine0,Ifem)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--VAFUNCMD8
!
      DOUBLE PRECISION B(NDF,*) , F(NDF,*) , Dr(NDF,*) , X(NXDm,*) , &
     &                 Totener2
      INTEGER Id(NDF,*) , Ix(4,*) , idf , Ifem
      LOGICAL Moveatoms , Movedisl , Fullfield
!
!---- Calculate out-of-balance force residual, strain energy, external
!     work,
!---- and total potential energy associated with mesh
!----
!
      INTEGER i , j
      DOUBLE PRECISION ftot , sed , work , etot , edd , strainenergy
      DOUBLE PRECISION TOTener , STRener , Straine0
      COMMON /ENERDAT/ TOTener , STRener
!
!     Calculate energy and force residual
!
!
!     apply boundary conditions
!
      DO i = 1 , NDF
         DO j = 1 , NUMnp
            IF ( Id(i,j)==1 ) then 
               B(i,j) = TIMe*F(i,j)
            endif
         ENDDO
      ENDDO
!
!
!
 
      edd = 0
      CALL INITIALISEENERGY(.TRUE.,.TRUE.,Id,Dr,F)
!      print *, 'Ifem in VafuncMD is', iFem
      CALL FEM_MOVE_PAD(X,B,Ix,F,TIMe,Z_Length,Id,ISRelaxed,edd,Dr,&
     &                  Fullfield,Movedisl,Straine0,NUMel,AVEvirst,&
     &                  SYStemp,Ifem,MOVed)
!
!     initialize
!
!CCCC
      IF ( Fullfield ) THEN
         DO i = 1 , NUMnp
            IF ( ISRelaxed(i)==2 ) Dr(1:NDF,i) = 0.0D0
         ENDDO
      ELSE
         CALL INITIALISEENERGY(.TRUE.,.TRUE.,Id,Dr,F)
      ENDIF
      print *, 'Ending Vafuncmd'
!
!     atomistic energy and forces.
!
      END SUBROUTINE VAFUNCMD
!*==ma06.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!*******************************************************************
      SUBROUTINE MA06(Id,X,Ix,F,B,Dr,Db,Input,Itx, lmp)
      USE MOD_GLOBAL
      USE MOD_DISLOCATION
      USE MOD_TIMMING
      use lammps
      IMPLICIT NONE
!*--MA0669
!
      DOUBLE PRECISION X(NXDm,NUMnp) , B(NDF,NUMnp) , Dr , Db , F
      INTEGER Id(NDF,NUMnp) , Ix(NEN1,NUMel) , Itx
      LOGICAL addedslip , lostslip , movedisl
      CHARACTER*80 Input , filename
!
      INTEGER lower , upper , NEXT , maxfn , idum , iprint , n , i , j ,&
     &        iat , ii(3) , k , maxfn2 , natoms , iprint2 , rseed
      DOUBLE PRECISION dum , dsmax , dsmax2 , dfn , tolm , ener , &
     &                 LENGTHTOL , xc(2) , xx(2) , cutoffr2 , rcutsq , &
     &                 told
      PARAMETER (LENGTHTOL=1.D-4)
      INTEGER MINflag
      COMMON /CGFLAG/ MINflag
      CHARACTER(LEN=16) :: md_status
 
      LOGICAL DEBug , md
      type (C_ptr) :: lmp
      COMMON /DEBUGGER/ DEBug
 
 
 
!     calling format:   ma06,,rseed,,
!     rseed:   integer, seed for random # generator
!
 
      CALL CPU_TIME(CT1)
 
      lower = 4
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,rseed,dum,1)
      lower = upper
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,idum,dsmax,2)
      lower = upper
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,iprint,dum,1)
 
      IDTemp = .FALSE.
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)<=0 ) THEN
            IDTemp(1:NDF,i) = .TRUE.
         ELSE
            DO j = 1 , NDF
               IDTemp(j,i) = (Id(j,i)==1)
            ENDDO
         ENDIF
      ENDDO
      DEBug = .FALSE.
 
!!$      CALL DOSTEPS(n,B,Dr,Db,ener,tolm,iprint,dsmax,rseed,dfn,Id,X,Ix,F,&
!!$     &             Itx,.TRUE.,addedslip,lostslip,.TRUE.,.TRUE.)
!!$ 
      CALL DOSTEPS_lammps(n,B,Dr,Db,ener,tolm,iprint,dsmax,rseed,dfn,Id,X,Ix,F,&
     &             Itx,.TRUE.,addedslip,lostslip,.TRUE.,.TRUE.,lmp)
 
 
 
      CALL CPU_TIME(CT3)
      CT6 = CT6 + CT3 - CT1
      PRINT * , ' '
      PRINT * , '****amount of cpu time devoted to various tasks '
      PRINT * , ' in detection band = ' , CT5
      PRINT * , ' in md = ' , CT4
      PRINT * , ' in fem = ' , CT7
      PRINT * , ' in ma06 = ' , CT6
      PRINT * , ' '
      IF ( MOVed ) THEN
         PRINT * , 'Crack tip moved to ' , XTIp_init(1)
         PRINT * , 'Repeating calculation without increasing load'
!$$$           Moved = .false.
!$$$           goto 10
      ENDIF
 
      END SUBROUTINE MA06
 
