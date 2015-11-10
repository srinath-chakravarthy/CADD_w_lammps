!*==mod_material.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*************************************************************
!**
!**   MODULE mod_material : contains definition of model materials and
!**   related
!**   routines.
!**
!**
!**   Variable Definitions:
!**   ---------------------
!**
!**   Integer Variables:
!**   nmaterials           -  Number of materials
!**   nbasismax            -  Maximum number of basis atoms
!**
!**   Type(bravaismat) Variables:
!**   material(nmaterials) -  All data relevant to the material
!**   structure:
!**   %name -  Materials name (character string)
!**   %bvec(:,i) -  Components of Bravais vector i
!**   %nbasis -  Number of atoms per Bravais site
!**   %basis(:,i) -  Coordinates of basis atom i
!**   %ispec(nbasis) -  Atomic species of basis atoms (use periodic
!**   table numbers)
!**   %bmag(:) -  Magnitude of Bravais vectors
!**   %volume -  Bravais cell volume
!**   %structure -  3 chr string identifying the structure
!**   (fcc, hcp, hex, etc)cubic lattice constant
!**   %a0     -  cubic lattice constant
!**   %cc     -  cubic elastic constants
!**
!**   Contains Routines:
!**   ------------------
!**
!**   ReadMaterials   - Reads in the number of materials and their
!**   properties
!**   OutputMaterials - Prints out information on the materials read in
!**
!**********************************************
 
      MODULE MOD_MATERIAL
      IMPLICIT NONE
!*--MOD_MATERIAL43
 
!     * Type Defintions
      TYPE BRAVAISMAT
         CHARACTER(LEN=20) :: NAME
         DOUBLE PRECISION BVEC(3,3)
         INTEGER NBASIS
         DOUBLE PRECISION , DIMENSION(:,:) , POINTER :: BASIS
         INTEGER , DIMENSION(:) , POINTER :: ISPEC
         DOUBLE PRECISION BMAG(3)
         DOUBLE PRECISION VOLUME , A0 , CC(6,6)
         CHARACTER*3 STRUCTURE
      END TYPE BRAVAISMAT
 
!     * Variable Definintions
      INTEGER nmaterials , nbasismax
      TYPE (BRAVAISMAT)  , DIMENSION(:) , POINTER::material
 
      CONTAINS
 
!------------------------------------------------------------------
! ReadMaterials -- Read in relevant data for all model materials
!                  and store in the structured array material
!
!      Passed Parameters :
!                      none
!
!      Module Parameters :
!          nmaterials (out) : number of materials
!          nbasismax  (out) : maximum number of basis atoms
!         material(:) (out) : material data (see defn above)
!
!      Algorithm :
!            Read in number of materials and then for each material
!            read in all relevant data as defined above.
!
!      Notes :
!            Does not currently allow for non-stoichiometric structures
!
!      Author :
!            E.B.Tadmor (12/31/97)
!
!      Revisions :
!
!--
      SUBROUTINE READMATERIALS(Ndf,Key)
 
      USE MOD_FILE
      IMPLICIT NONE
!*--READMATERIALS92
      INTEGER Ndf
      CHARACTER*4 Key
 
!** Local Variables **!
      INTEGER i , j , iunit , icc , jcc , ncc
      DOUBLE PRECISION DET33 , ctmp
      LOGICAL error
 
      CHARACTER(LEN=80) datafile
 
      READ(input_file_unit, *)  NMAterials
      PRINT * , 'number of materials' , NMAterials
      IF ( NMAterials<1 ) THEN
         PRINT * , '***ERROR: Illegal number of materials specified.'
         STOP
      ENDIF
      ALLOCATE (MATerial(NMAterials))
      NBAsismax = 0
      DO i = 1 , NMAterials
         IF ( Key=='dire' ) THEN
            iunit = 5
         ELSE
            READ (input_file_unit,'(a)') datafile
            IF ( .NOT.FILEEXISTS(datafile,.TRUE.) ) STOP
            CALL IOFILE(datafile,'formatted  ',iunit,.TRUE.)
         ENDIF
         READ (iunit,*) MATerial(i)%NAME
         READ (iunit,*) (MATerial(i)%BVEC(:,j),j=1,3)
         READ (iunit,*) MATerial(i)%NBASIS
         ALLOCATE (MATerial(i)%BASIS(3,MATerial(i)%NBASIS))
         MATerial(i)%BASIS = 0.
         ALLOCATE (MATerial(i)%ISPEC(MATerial(i)%NBASIS))
         MATerial(i)%ISPEC = 0
         READ (iunit,*) (MATerial(i)%BASIS(:,j),MATerial(i)%ISPEC(j),&
     &                  j=1,MATerial(i)%NBASIS)
         READ (iunit,'(a3)') MATerial(i)%STRUCTURE
         READ (iunit,*) MATerial(i)%A0
         READ (iunit,*) ncc
         MATerial(i)%CC = 0.D0
         DO j = 1 , ncc
            READ (iunit,*) icc , jcc , ctmp
            MATerial(i)%CC(icc,jcc) = ctmp
            MATerial(i)%CC(jcc,icc) = ctmp
         ENDDO
         DO j = 1 , 3
            MATerial(i)%BMAG(j) = DSQRT(DOT_PRODUCT(MATerial(i)%BVEC(:,j&
     &                            ),MATerial(i)%BVEC(:,j)))
         ENDDO
         MATerial(i)%VOLUME = DET33(MATerial(i)%BVEC)
         IF ( MATerial(i)%VOLUME<=0.D0 ) THEN
            WRITE (*,*) &
     &'** ERROR: unit cell volume is negative                 or zero fo&
     &r material' , i
            WRITE (*,*) 'Make sure Bravais lattice is right-handed'
            STOP
         ENDIF
         IF ( NBAsismax<MATerial(i)%NBASIS ) NBAsismax = MATerial(i)&
     &        %NBASIS
         IF ( Key/='dire' ) CLOSE (iunit)
      ENDDO
 
      error = .FALSE.
      IF ( Ndf<3*NBAsismax ) THEN
         WRITE (*,*) '***ERROR: ndf is too small to accomodate the'
         WRITE (*,*) '          defined materials.'
         WRITE (*,*) '          Increase ndf to:' , 3*NBAsismax
         error = .TRUE.
      ENDIF
      IF ( error ) STOP
      END SUBROUTINE READMATERIALS
 
   !------------------------------------------------------------------
   ! OutputMaterials -- Prints information on loaded materials in an
   !                    orderly fashion
   !
   !      Passed Parameters :
   !                      none
   !
   !      Module Parameters :
   !          nmaterials (out) : number of materials
   !         material(:) (out) : material data (see defn above)
   !
   !      Algorithm :
   !            self explanatory
   !
   !      Notes :
   !          none
   !
   !      Author :
   !            E.B.Tadmor (12/31/97)
   !
   !      Revisions :
   !              none
   !
   !--
      SUBROUTINE OUTPUTMATERIALS()
 
      IMPLICIT NONE
!*--OUTPUTMATERIALS191
 
      !** Local Variables **!
      INTEGER i , j
 
      PRINT *
      PRINT '(a)' , 'MATERIAL INFORMATION'
      DO i = 1 , NMAterials
         PRINT *
         PRINT '("Mat #",i2," : ",A)' , i , MATerial(i)%NAME
         PRINT '(a)' , 'Bravais lattice vectors:'
         PRINT '("a",i1," = ",3f10.5)' , (j,MATerial(i)%BVEC(:,j),j=1,3)
         PRINT '(a)' , 'Basis atoms coordinates and species:'
         PRINT '("Atom #",i2," : ",3f10.5,5x,i3)' , &
     &         (j,MATerial(i)%BASIS(:,j),MATerial(i)%ISPEC(j),j=1,&
     &         MATerial(i)%NBASIS)
      ENDDO
      PRINT *
 
      END SUBROUTINE OUTPUTMATERIALS
 
      END MODULE MOD_MATERIAL
