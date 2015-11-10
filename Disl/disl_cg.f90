!*==cg_setup.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
! $Id: disl_cg.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE CG_SETUP
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_SETUP8
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
      DATA I_Cg/0/
!
      INTEGER i
      CHARACTER*80 error_message
!
      IF ( I_Cg/=0 ) THEN
         error_message = 'cg_setup: improper exit from CG!'
         CALL ERROR_HANDLER(error_message)
      ELSE
         I_Cg = 1
      ENDIF
      IF ( I_Disl/=1 ) THEN
         error_message = 'cg_setup: call disl_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      CALL CG_RESET
!
      END SUBROUTINE CG_SETUP
!*==cg_exit.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE CG_EXIT
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_EXIT39
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_exit: exit from CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ELSE
         I_Cg = 0
      ENDIF
!
      END SUBROUTINE CG_EXIT
!*==cg_reset.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE CG_RESET
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_RESET62
      DOUBLE PRECISION cg_gdotg
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      INTEGER i
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_reset: CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , NDIsl
         PK_b(i) = 0.0D0
      ENDDO
      END SUBROUTINE CG_RESET
!*==cg_gdotg.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      FUNCTION CG_GDOTG()
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_GDOTG88
      DOUBLE PRECISION CG_GDOTG
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      INTEGER i
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_gdotg: CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      CG_GDOTG = 0.0D0
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) CG_GDOTG = CG_GDOTG + PK_f(i)*PK_f(i)
      ENDDO
      END FUNCTION CG_GDOTG
!*==cg_z.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE CG_Z(Z,Cgsmax)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_Z115
      DOUBLE PRECISION Z , Cgsmax
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      INTEGER i
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_z: CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      Cgsmax = 0.0D0
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) THEN
            PK_b(i) = PK_b(i)*BURg_length(i)*Z + PK_f(i)
            IF ( DABS(PK_b(i))>Cgsmax ) Cgsmax = DABS(PK_b(i))
            PK_b(i) = PK_b(i)/BURg_length(i)
         ELSE
            PK_b(i) = 0.0D0
         ENDIF
      ENDDO
      END SUBROUTINE CG_Z
!*==cg_saxpy.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE CG_SAXPY(Alpha)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_SAXPY148
      DOUBLE PRECISION Alpha
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      INTEGER FE_LOCATE , i , elem_old
      DOUBLE PRECISION b
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_saxpy: CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)/=0 ) THEN
            R_Disl(1,i) = R_Disl(1,i) + Alpha*PK_b(i)*BURgers(1,i)
            R_Disl(2,i) = R_Disl(2,i) + Alpha*PK_b(i)*BURgers(2,i)
            elem_old = ELEm_disl(i)
            ELEm_disl(i) = FE_LOCATE(R_Disl(1,i),ELEm_disl(i))
            IF ( ELEm_disl(i)==0 ) THEN
               IF ( elem_old>0 ) THEN
                  ELEm_disl(i) = -elem_old
               ELSE
                  ELEm_disl(i) = elem_old
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      END SUBROUTINE CG_SAXPY
!*==cg_gdots.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      DOUBLE PRECISION FUNCTION CG_GDOTS()
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_GDOTS185
      DOUBLE PRECISION alpha
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      INTEGER i
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_gdots: CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      CG_GDOTS = 0.D0
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) CG_GDOTS = CG_GDOTS - PK_b(i)&
     &        *BURg_length(i)*PK_f(i)
      ENDDO
      END FUNCTION CG_GDOTS
!*==cg_sdots.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      DOUBLE PRECISION FUNCTION CG_SDOTS()
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--CG_SDOTS213
!
      DOUBLE PRECISION PK_b(MAX_DISL)
      INTEGER I_Cg
      COMMON /CG_COMMON/ PK_b , I_Cg
!
      INTEGER i
      CHARACTER*80 error_message
!
      IF ( I_Cg/=1 ) THEN
         error_message = 'cg_sdots: CG without entering!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      CG_SDOTS = 0.D0
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) CG_SDOTS = CG_SDOTS + PK_b(i)&
     &        *BURg_length(i)*PK_b(i)*BURg_length(i)
      ENDDO
      END FUNCTION CG_SDOTS
 
 
!
! $Log: disl_cg.f,v $
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.5  2002/03/0503:00:45  shilkrot
! Moved pk_b out of dislocations.par into the cg common block.
!
! Revision 1.4  2002/02/2822:43:08  shilkrot
! Removed the condition elem_disl != 0 for moving dislocations in cg_sax
!
! Revision 1.3  2002/02/2520:54:05  shilkrot
! Added a routine disl_reset.
! Changed the code to use burg_length.
!
! Revision 1.2  2001/12/1307:31:24  shilkrot
! Implemented breadth first search to find the element number for
! a dislocation. Changed the interface of fe_locate to use the starting
! element for the search. Old fe_locate is in fem_services.
! Changed the interface of fem_setup. Now two arrays used as temp space 
! passed from outside as the last two parameters.
!
! Revision 1.1  2001/11/1303:50:37  shilkrot
! Routines computing scalar products and moving dislocations for the
! conjugate gradient.
!
!
