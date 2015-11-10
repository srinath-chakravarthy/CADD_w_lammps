!*==rest.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**   rest : read/write the restart file
!**
      SUBROUTINE REST(Id,X,Ix,Itx,F,B,Logic,Key)
 
      USE MOD_GLOBAL
      USE MOD_BOUNDARY
      USE MOD_GRAIN
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--REST12
!
      DOUBLE PRECISION X(*) , F(*) , B(*)
      INTEGER Id(*) , Ix(*) , Logic , Itx(*)
      CHARACTER*4 Key
      LOGICAL repsuptodate_dummy
!
!---- read/write restart files
!
      CHARACTER*4 xkey1 , xkey2
      DOUBLE PRECISION XCRit
      COMMON /CCRIT / XCRit(2)
      LOGICAL BOUndary
      INTEGER ICRit , ISIgn , NODecrit , IDFcrit
      DOUBLE PRECISION ARCrate
      COMMON /CSTEP / BOUndary , ICRit , ARCrate , ISIgn , NODecrit , &
     &                IDFcrit
      INTEGER neqtot , i , j , k , nsdim , nxsjdim , niddim , nixdim , &
     &        ntdim , njddim , ng , nshpdim , nxdim , nfdim , nby
 
      ! get number of grains
!
      xkey1 = 'read'
      xkey2 = 'writ'
      IF ( Key==xkey1 ) THEN
         REWIND Logic
         READ (Logic) NUMnp , NUMel , NEQ
         IF ( NUMnp>MAXnp ) THEN
            WRITE (6,&
     &'(///''  **error** detected by subroutine rest''           /'' num&
     &ber of nodes exceeds maximum'')')
            WRITE (6,'(''  **required = '',i10)') NUMnp
            WRITE (6,'(''  **maximum  = '',i10)') MAXnp
            STOP
         ENDIF
         IF ( NUMel>MAXel ) THEN
            WRITE (6,&
     &'(///''  **error** detected by subroutine rest''           /'' num&
     &ber of elements exceeds maximum'')')
            WRITE (6,'(''  **required = '',i10)') NUMel
            WRITE (6,'(''  **maximum  = '',i10)') MAXel
            STOP
         ENDIF
         READ (Logic) TIMe , DT , TIMe , TIMeol
         READ (Logic) BOUndary , ICRit , ISIgn , NODecrit , IDFcrit , &
     &                XCRit(1) , XCRit(2)
         neqtot = NEQ + NSTad*NUMel
         READ (Logic) (B(i),i=1,neqtot)
         nsdim = NSDm*NQUad*NUMel
         nshpdim = NSHpdm*NQUad*NUMel
         nxsjdim = NQUad*NUMel
         niddim = NDF*NUMnp
         nxdim = NXDm*NUMnp
         nixdim = NEN1*NUMel
         nfdim = NDF*NUMnp
         ntdim = NUMnp
         njddim = NDF*NUMnp
         READ (Logic) (Id(i),i=1,niddim)
         READ (Logic) (X(i),i=1,nxdim)
         READ (Logic) (Ix(i),i=1,nixdim)
         READ (Logic) (Itx(i),i=1,3*NUMel)
         READ (Logic) (F(i),i=1,nfdim)
         READ (Logic) NCE , NCB
         IF ( NCE>NCEmax ) CALL INCREASEELIST(NCE-NCEmax+100)
         READ (Logic) ((ELIst(i,j),i=1,2),j=1,NCE)
!--Representative atom data
         READ (Logic) ng
         IF ( ng/=NGRains ) THEN
            WRITE (6,&
     &'(///''  **error** detected by subroutine rest''           /'' num&
     &ber of grains in file different than current value'')')
            WRITE (6,'(''  **in file  = '',i10)') ng
            WRITE (6,'(''  **current  = '',i10)') NGRains
            STOP
         ENDIF
         READ (Logic) repsuptodate_dummy
         READ (Logic) (ENErgy(j),j=1,NUMnp)
         READ (Logic) (ISRelaxed(j),j=1,NUMnp)
!--   dislocations
         READ (Logic) NDIsl
         IF ( NDIsl>MAX_DISL ) STOP 'increase max_disl'
         READ (Logic) BURgers(1:3,1:NDIsl) , R_Disl(1:3,1:NDIsl)
         READ (Logic) THEta_e(1:NDIsl) , THEta_s(1:NDIsl)
         READ (Logic) BURg_length(1:NDIsl)
         READ (Logic) ELEm_disl(1:NDIsl)
         DO i = 1 , NDIsl
            IF ( ELEm_disl(i)/=0 )&
     &           CALL SLIPRANGE(BURgers(1:3,i),R_Disl(1:3,i),&
     &           DISl_range(1:2,i),DISl_index(i))
         ENDDO
      ENDIF
 
      IF ( Key==xkey2 ) THEN
         REWIND Logic
         WRITE (Logic) NUMnp , NUMel , NEQ
         WRITE (Logic) TIMe , DT , TIMe , TIMeol
         WRITE (Logic) BOUndary , ICRit , ISIgn , NODecrit , IDFcrit , &
     &                 XCRit(1) , XCRit(2)
         neqtot = NEQ + NSTad*NUMel
         WRITE (Logic) (B(i),i=1,neqtot)
         nsdim = NSDm*NQUad*NUMel
         nshpdim = NSHpdm*NQUad*NUMel
         nxsjdim = NQUad*NUMel
         niddim = NDF*NUMnp
         nxdim = NXDm*NUMnp
         nixdim = NEN1*NUMel
         nfdim = NDF*NUMnp
         ntdim = NUMnp
         njddim = NDF*NUMnp
         WRITE (Logic) (Id(i),i=1,niddim)
         WRITE (Logic) (X(i),i=1,nxdim)
         WRITE (Logic) (Ix(i),i=1,nixdim)
         WRITE (Logic) (Itx(i),i=1,3*NUMel)
         WRITE (Logic) (F(i),i=1,nfdim)
         WRITE (Logic) NCE , NCB
         WRITE (Logic) ((ELIst(i,j),i=1,2),j=1,NCE)
         WRITE (Logic) NGRains
!--Representative atom data
         WRITE (Logic) repsuptodate_dummy
         WRITE (Logic) (ENErgy(j),j=1,NUMnp)
         WRITE (Logic) (ISRelaxed(j),j=1,NUMnp)
!--   dislocations
         WRITE (Logic) NDIsl
         WRITE (Logic) BURgers(1:3,1:NDIsl) , R_Disl(1:3,1:NDIsl)
         WRITE (Logic) THEta_e(1:NDIsl) , THEta_s(1:NDIsl)
         WRITE (Logic) BURg_length(1:NDIsl)
         WRITE (Logic) ELEm_disl(1:NDIsl)
      ENDIF
      END SUBROUTINE REST
 
 
