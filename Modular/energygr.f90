!*==strainenergy.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**  StrainEnergy : computes total strain energy in the mesh
!**
!**   Non-Obvious Parameters : NONE
!**
!**   Algorithm :-
!**        Loop over the repatoms and add appropriate nodal energy
!**        depending on local or nonlocal status.
!--
      DOUBLE PRECISION FUNCTION STRAINENERGY()
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--STRAINENERGY14
      INTEGER i
      DOUBLE PRECISION sed
      sed = 0.0
      DO i = 1 , NUMnp
!         if(IsRelaxed(i).ge.1) then
         IF ( ISRelaxed(i)/=0 ) sed = sed + ENErgy(i)
      ENDDO
      IF ( NUMnpp1>NUMnp ) sed = sed + ENErgy(NUMnpp1)
      STRAINENERGY = sed
      END FUNCTION STRAINENERGY
!***********************************************************************
