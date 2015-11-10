!*==mod_interstitial.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
!
!
!
!
 
      MODULE MOD_INTERSTITIAL
      IMPLICIT NONE
!*--MOD_INTERSTITIAL9
 
!* Type Defintions
      TYPE INTERSTITIAL_ATOMS
         INTEGER NATOMS , MAXATOMS
         DOUBLE PRECISION , DIMENSION(:,:) , POINTER :: R
         INTEGER , DIMENSION(:) , POINTER :: ELEMENT
      END TYPE INTERSTITIAL_ATOMS
 
      CONTAINS
 
!- Initializes the collection of interstitials
      SUBROUTINE INIT(Interstitial)
      IMPLICIT NONE
!*--INIT23
 
      TYPE (INTERSTITIAL_ATOMS) Interstitial
      INTEGER maxatoms
 
      maxatoms = 100
 
      ALLOCATE (Interstitial%R(3,maxatoms))
      ALLOCATE (Interstitial%ELEMENT(maxatoms))
      Interstitial%maxatoms = maxatoms
 
      END SUBROUTINE INIT
 
 
!     - Reads in interstitial data
      SUBROUTINE READINTERSTITIALDATA(Interstitial)
      IMPLICIT NONE
!*--READINTERSTITIALDATA40
 
      TYPE (INTERSTITIAL_ATOMS) Interstitial
      CHARACTER(LEN=80) :: input_file
      DOUBLE PRECISION x , y , z
      INTEGER ielement , iatoms , n
 
      input_file = 'interstitial.dat'
      OPEN (UNIT=10,FILE=input_file,STATUS='old')
 
      n = 0
      DO iatoms = 1 , Interstitial%MAXATOMS
         READ (10,*,END=100) x , y , z , ielement
         n = n + 1
         Interstitial%R(1,n) = x
         Interstitial%R(2,n) = y
         Interstitial%R(3,n) = z
         Interstitial%ELEMENT(n) = ielement
      ENDDO
!lose(10)
 
 100  Interstitial%NATOMS = n
 
      END SUBROUTINE READINTERSTITIALDATA
 
 
      END MODULE MOD_INTERSTITIAL
