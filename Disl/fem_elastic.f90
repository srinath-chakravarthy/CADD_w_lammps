!*==fe_elastic.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
! $Id: fem_elastic.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE FE_ELASTIC(Elas_const)
!
!12345678901234567890123456789012345678901234567890123456789012345678901
!         1         2         3         4         5         6         7
!
      IMPLICIT NONE
!*--FE_ELASTIC11
      DOUBLE PRECISION Elas_const(6,6)
      DOUBLE PRECISION CC(6,6)
      DOUBLE PRECISION XE , XNU , XLAmbda , XMU
      INTEGER I_Elas
!
      COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
      DATA I_Elas/0/
!
      DOUBLE PRECISION xinv1 , xinv2
      INTEGER i , j
      CHARACTER*80 error_message
!
      IF ( I_Elas==0 ) THEN
         I_Elas = 1
      ELSE
         error_message = 'Second call to fe_elastic'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO j = 1 , 6
         DO i = 1 , 6
            CC(i,j) = Elas_const(i,j)
         ENDDO
      ENDDO
!
!  Hirth and Lothe
!
!cHex        xinv1 = cc(1,1)+cc(2,2)+cc(3,3)+cc(4,4)+cc(4,4)+cc(5,5)+cc(
!cHex     &        +cc(6,6)+cc(6,6)
!cHex        xinv2 = cc(1,1)+cc(1,2)+cc(1,3)+cc(2,1)+cc(2,2)+cc(2,3)
!cHex     &        +cc(3,1)+cc(3,2)+cc(3,3)
!
!cHex        xlambda = (2.0d0*xinv2-xinv1)/15.0d0
!cHex        xmu = (3.0d0*xinv1-xinv2)/30.0d0
!
!cHex        xe = xmu*(2.0d0+xlambda/(xlambda+xmu))
!cHex        xnu = xlambda/2.0d0/(xlambda+xmu)
 
!c--JS: 2D plane strain elastic coefficient
      XNU = (CC(1,1)+3.0D0*CC(1,2)-2.0D0*CC(6,6))&
     &      /(4.0D0*CC(1,1)+4.0D0*CC(1,2))
      XE = (1.0D0+XNU)*(CC(1,1)-CC(1,2)+2.0D0*CC(6,6))/2.0D0
      XMU = XE/(2.0D0*(1.0D0+XNU))
      XLAmbda = XE*XNU/((1.0D0+XNU)*(1-2.0D0*XNU))
!$$$	e_dd = xe
!$$$        nu_dd = xnu
!$$$        mu_dd = xmu
      WRITE (*,*) "xnu,xe,xlambda,xmu" , XNU , XE , XLAmbda , XMU
!
      END SUBROUTINE FE_ELASTIC
 
!
! $Log: fem_elastic.f,v $
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.2  2001/07/1206:36:49  shilkrot
! Elastic constants are now passed to fem_setup and further to fe_elasti
!
! Revision 1.1  2001/06/2502:38:09  shilkrot
! Sets up the stiffness matrix and the isotropic moduli
!
