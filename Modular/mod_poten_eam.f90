!*==mod_poten.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
!**
!**  MODULE mod_potential : Embedded Atom Method (EAM) Potential Module
!**
!**
!**  Variable Definitions:
!**  ---------------------
!**
!**  Dynamo EAM variables
!**
!**  Contains Required Routines:
!**  ---------------------------
!**
!**  ReadConstitutive - Reads in potential specific data
!**  Nonlocal()       - Computes energy, stress and moduli using nonloca
!**                     limit of formulation
!**
!***********************************************************************
 
      MODULE MOD_POTEN
        IMPLICIT NONE
!*--MOD_POTEN23
 
!--This file is a modified form of "dyn87.inc" which is the use file
!  for dynamo v8.7.  This is to store the interatomic potentials as done
!  in dynamo
 
      INTEGER , PARAMETER :: NGRID = 1000
      INTEGER , PARAMETER :: NGRIDAR = 1000
      INTEGER , PARAMETER :: NELMAX = 3
      INTEGER ntypes , nrhoar , nrar , nrho , nr
      INTEGER ielement(NELMAX) , netype(NELMAX)
      DOUBLE PRECISION amass(NELMAX)
      DOUBLE PRECISION frho(NGRID,NELMAX) , z2r(NGRID,NELMAX,NELMAX) , &
     &                 rhor(NGRID,NELMAX) , drho , drad , rcutsq , sqrtrcutsq
      DOUBLE PRECISION frhoar(NGRIDAR,NELMAX) , frhoar1(NGRIDAR,NELMAX)&
     &                 , frhoar2(NGRIDAR,NELMAX) , &
     &                 frhoar3(NGRIDAR,NELMAX) , frhoar4(NGRIDAR,NELMAX)&
     &                 , frhoar5(NGRIDAR,NELMAX) , &
     &                 frhoar6(NGRIDAR,NELMAX) , frhoar7(NGRIDAR,NELMAX)&
     &                 , rhorar(NGRIDAR,NELMAX) , &
     &                 rhorar1(NGRIDAR,NELMAX) , rhorar2(NGRIDAR,NELMAX)&
     &                 , rhorar3(NGRIDAR,NELMAX) , &
     &                 rhorar4(NGRIDAR,NELMAX) , rhorar5(NGRIDAR,NELMAX)&
     &                 , rhorar6(NGRIDAR,NELMAX) , &
     &                 rhorar7(NGRIDAR,NELMAX) , &
     &                 z2rar(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar1(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar2(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar3(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar4(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar5(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar6(NGRIDAR,NELMAX,NELMAX) , &
     &                 z2rar7(NGRIDAR,NELMAX,NELMAX) , drhoar , drar
      INTEGER mapspecies(103)
 
      DOUBLE PRECISION CON1
      PARAMETER (CON1=100.D0)
 
      CONTAINS
 
                   !!! REQUIRED ROUTINES !!!
 
!**------------------------------------------------------------------
!** ReadConstitutive : Read in potential specific data
!**
!--
      SUBROUTINE READCONSTITUTIVE()
 
      USE MOD_GRAIN
      USE MOD_MATERIAL
      IMPLICIT NONE
!*--READCONSTITUTIVE75
 
      !** Local Variables **!
      INTEGER i , imat , ibasis
 
      ! Make sure Grains have been defined
      IF ( NGRains<1 ) THEN
         WRITE (*,*) '***ERROR: Grain definitions must come before'
         WRITE (*,*) '          constitutive information'
         WRITE (*,*)
         STOP
      ENDIF
 
      ! Initialize species mapping vector
      MAPspecies = 0
      ! Read Potentials
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*) '-------------here comes some dynamo output ------'
      WRITE (*,*)
      CALL READPOTEN()
      WRITE (*,*)
      WRITE (*,*) '-------------back to our output -------------'
      WRITE (*,*)
      WRITE (*,*)
 
      ! Construct species mapping vector
      DO i = 1 , NTYpes
         MAPspecies(IELement(i)) = i
      ENDDO
 
      ! Check that all species in loaded materials are handled by pot
      DO imat = 1 , NMAterials
         DO ibasis = 1 , MATerial(imat)%NBASIS
            IF ( MAPspecies(MATerial(imat)%ISPEC(ibasis))==0 ) THEN
               WRITE (*,*)
               WRITE (*,99001) imat , MATerial(imat)%ISPEC(ibasis)
99001          FORMAT ('***ERROR: Material ',i2,' contains species ',i2)
               WRITE (*,*) &
     &                   '   But only the following species are defined'
               WRITE (*,*) '   in the EAM files:'
               DO i = 1 , 103
                  IF ( MAPspecies(i)/=0 ) WRITE (*,*) i
               ENDDO
               WRITE (*,*)
               WRITE (*,*) '*** Run Stopped due to an ERROR'
               STOP
            ENDIF
         ENDDO
      ENDDO
 
      END SUBROUTINE READCONSTITUTIVE
 
!**---------------------------------------------------------------------
!** Nonlocal : find energy, force and stiffness contribution due to
!**            a nonlocal representative atom
!**
!**    Parameters :
!**          Look at the element routine
!**
!**
!**    Algorithm :
!**          Make Representative crystallite.
!**
!**
!**
!**    Notes :
!**           Specific to embedded atom.
!**
!**
!**
!**    WARNING : This will work only for NDA = 1
!**
!**
!--
      SUBROUTINE NONLOCAL(Id,X,Ix,F,B,Dr,Fls,Flw,Irep)
 
      USE MOD_MATERIAL
      USE MOD_GLOBAL
      USE MOD_GRAIN
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--NONLOCAL157
 
!--Variables transferred
      INTEGER Irep
 
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*) , F(NDF,*) , Dr(NDF,*)
 
      INTEGER Id(NDF,*) , Ix(NEN1,*)
 
      LOGICAL Fls , Flw , wn
 
!--Local Variables
 
      DOUBLE PRECISION x1(2) , x2(2) , x3(2) , s0(3)
      DOUBLE PRECISION , ALLOCATABLE :: p1(:) , p2(:)
 
      INTEGER  :: ineigh , idf , idx , i0 , isp1 , isp0 , i , j , k , node , l , m , ndnn
 
      DOUBLE PRECISION :: rhosum , phisum , r , ef , vol , df , d2f , ph 
      DOUBLE PRECISION :: dp , d2p , em , dem , d2em , druai , drubj , sed
 
      DOUBLE PRECISION force(3)
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
 
!     hack: assume all atoms are the same species for now, assume only
!     one material defined
      IF ( Irep/=NUMnpp1 ) THEN
         isp0 = ATOmspecie(Irep)
      ELSE
         isp0 = -1
      ENDIF
 
!     Ensure there is sufficient storage allocated for out-of-balance
!     force vector and stiffness matrix and initialize them.
      ndnn = NDF*(NNEips+1)
      IF ( Fls ) THEN
         ALLOCATE (p1(ndnn),p2(ndnn))
         p1 = 0.D0
         p2 = 0.D0
      ENDIF
!
!     Loop over all atoms in local crystal summing energy
!     and force contributions
      rhosum = 0.D0
      phisum = 0.D0
      DO ineigh = 1 , NNEips
         IF ( JNEigh(ineigh)/=NUMnpp1 ) THEN
            isp1 = ATOmspecie(JNEigh(ineigh))
         ELSE
            isp1 = -1
         ENDIF
         r = SQRT(RNEigh(ineigh))
!         call edens(r,isp1,ef,df,d2f,.true.,fls,.false.)
!         call pair(r,isp0,isp1,ph,dp,d2p,.true.,fls,.false.)
         DEBugflag = 0
         CALL EDENS(r,isp1,ef,df,d2f,.TRUE.,Fls,.FALSE.)
         CALL PAIR(r,isp0,isp1,ph,dp,d2p,.TRUE.,Fls,.FALSE.)
 
         rhosum = rhosum + ef
         phisum = phisum + 0.5*ph
 
         IF ( Fls ) THEN
            DO idf = 1 , NDF
               i0 = NNEips*NDF + idf
               idx = (ineigh-1)*NDF + idf
               druai = DNEigh(idf,ineigh)/r
               p1(i0) = p1(i0) + df*druai
               p2(i0) = p2(i0) + 0.5D0*dp*druai
               p1(idx) = p1(idx) - df*druai
               p2(idx) = p2(idx) - 0.5D0*dp*druai
            ENDDO
         ENDIF
      ENDDO
 
 
      !Compute embedding energy
      CALL EMBED(rhosum,em,dem,d2em,isp0,Fls,.FALSE.)
 
      !Compute  energy    !CHANGED for nalpha (May 10,96)
                          !Changed for atomistic (Aug 25, 96)
      ENErgy(Irep) = em + phisum
 
!     Compute contribution to out-of-balance force
      IF ( Fls ) THEN
         DO ineigh = 1 , NNEips + 1
            IF ( ineigh>NNEips ) THEN
               node = Irep
            ELSE
               node = JNEigh(ineigh)
            ENDIF
            DO idf = 1 , NDF
               idx = (ineigh-1)*NDF + idf
               force(idf) = (dem*p1(idx)+p2(idx))
               Dr(idf,node) = Dr(idf,node) - force(idf)
            ENDDO
         ENDDO
	 if (irep == 4082) then 
	  write(*, '(A25,I5,1X,2(F15.8))') 'Force on atom in nonloc 1', irep, dr(1:2, irep)
	 end if	
	
!     Compute Virial Stress FCC material Only
 
!CJSTest         vol=(material(1)%a0**3)/4.d0
!CJSTest: for Hex Al only
         vol = 20.5052
 
         DO ineigh = 1 , NNEips
            node = JNEigh(ineigh)
            DO idf = 1 , NDF
               idx = (ineigh-1)*NDF + idf
               force(idf) = (dem*p1(idx)+p2(idx))
               VIRst(idf,1,node) = VIRst(idf,1,node) - force(idf)&
     &                             *DNEigh(1,ineigh)/vol
               VIRst(idf,2,node) = VIRst(idf,2,node) - force(idf)&
     &                             *DNEigh(2,ineigh)/vol
               VIRst(idf,3,node) = VIRst(idf,3,node) - force(idf)&
     &                             *DNEigh(3,ineigh)/vol
            ENDDO
         ENDDO
!
      ENDIF
 	 if (irep == 4082) then 
	  write(*, '(A25,I5,1X,2(F15.8))') 'Force on atom in nonloc 2', irep, dr(1:2, irep)
	 end if	

      IF ( Fls ) DEALLOCATE (p1,p2)
      END SUBROUTINE NONLOCAL
 
!**---------------------------------------------------------------------
!     Calculate Electron Density as a function of Distance between Atoms
      SUBROUTINE EDENS(R,Ispec,F,Df,D2f,Lf,Ldf,Ld2f)
 
      IMPLICIT NONE
!*--EDENS286
 
      !** Transferred Variables **!
      DOUBLE PRECISION , INTENT(IN) :: R
      INTEGER , INTENT(IN) :: Ispec
      DOUBLE PRECISION , INTENT(OUT) :: F , Df , D2f
      LOGICAL , INTENT(IN) :: Lf , Ldf , Ld2f
 
      !local variables
      DOUBLE PRECISION cutoffradius
      F = 0.D0
      Df = 0.D0
      D2f = 0.D0
      CALL RHFL(R,Ispec,F,Df,D2f,Lf,Ldf,Ld2f)
 
      END SUBROUTINE EDENS
 
!**---------------------------------------------------------------------
!     Pair Potential at Distance R
      SUBROUTINE PAIR(R,Isp1,Isp2,F,Df,D2f,Lf,Ldf,Ld2f)
 
      IMPLICIT NONE
!*--PAIR308
 
      !** Transferred Variables **!
      DOUBLE PRECISION , INTENT(IN) :: R
      INTEGER , INTENT(IN) :: Isp1 , Isp2
      DOUBLE PRECISION , INTENT(OUT) :: F , Df , D2f
      LOGICAL , INTENT(IN) :: Lf , Ldf , Ld2f
 
      !local variables
 
      F = 0.D0
      Df = 0.D0
      D2f = 0.D0
      CALL V2FL(R,Isp1,Isp2,F,Df,D2f,Lf,Ldf,Ld2f)
 
      END SUBROUTINE PAIR
 
!**---------------------------------------------------------------------
!     Embedding function and derivatives wrt rho
      SUBROUTINE FEMB(Rho,Ispec,F,Df,D2f,Lf,Ldf,Ld2f)
 
      IMPLICIT NONE
!*--FEMB330
 
      DOUBLE PRECISION , INTENT(IN) :: Rho
      INTEGER , INTENT(IN) :: Ispec
      DOUBLE PRECISION , INTENT(OUT) :: F , Df , D2f
      LOGICAL , INTENT(IN) :: Lf , Ldf , Ld2f
 
      CALL UUFL(Rho,Ispec,F,Df,D2f,Lf,Ldf,Ld2f) !VBS 23 Jan 96
 
      END SUBROUTINE FEMB
 
!**---------------------------------------------------------------------
!     Calculate Energy to Embed and Atom in Energy Density RHO
      SUBROUTINE EMBED(Rho,Ff,Scon,Dcon,Ispec,Fls,Flt)
      IMPLICIT NONE
!*--EMBED345
      DOUBLE PRECISION , INTENT(IN) :: Rho
      DOUBLE PRECISION , INTENT(OUT) :: Ff , Scon , Dcon
      INTEGER , INTENT(IN) :: Ispec
      LOGICAL , INTENT(IN) :: Fls , Flt
      CALL FEMB(Rho,Ispec,Ff,Scon,Dcon,.TRUE.,(Fls .OR. Flt),Flt)
      END SUBROUTINE EMBED
 
!**---------------------------------------------------------------------
!** rhfl - gives electron density
!**
!**     Parameters :
!**              r  (in) : distance r
!**            isp  (in) : species of the atom
!**              f (out) : electron density
!**             df (out) : first derivative of electron density
!**            d2f (out) : second derivative of electron density
!**             lf  (in) : true to compute f
!**            ldf  (in) : true to compute df
!**           ld2f  (in) : true to compute ld2f
!**
      SUBROUTINE RHFL(R,Isp,F,Df,D2f,Lf,Ldf,Ld2f)
      IMPLICIT NONE
!*--RHFL368
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION D2f , Df , F , one , p , R
      INTEGER Isp , k
!*** End of declarations inserted by SPAG
      DATA one/1.0/
      LOGICAL Lf , Ldf , Ld2f
 
      IF ( Isp<0 ) THEN
         IF ( Lf ) F = 0.D0
         IF ( Ldf ) Df = 0.D0
         IF ( Ld2f ) D2f = 0.D0
         RETURN
      ENDIF
 
      p = R/DRAr + 1.0
      k = p
      k = MIN0(k,NRAr-1)
      p = p - k
      p = MIN(p,one)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Lf ) F = ((RHOrar3(k,Isp)*p+RHOrar2(k,Isp))*p+RHOrar1(k,Isp))&
     &              *p + RHOrar(k,Isp)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      IF ( Ldf ) Df = (RHOrar6(k,Isp)*p+RHOrar5(k,Isp))&
     &                *p + RHOrar4(k,Isp)
 
      IF ( Ld2f ) D2f = RHOrar7(k,Isp)&
     &                  + (RHOrar7(k+1,Isp)-RHOrar7(k,Isp))*p
 
      END SUBROUTINE RHFL
 
!**---------------------------------------------------------------------
!** v2fl - gives pair potential
!**
!**     Parameters :
!**              r  (in) : distance r
!**            ity  (in) : species of the atom1
!**            jty  (in) : species of the atom2
!**              f (out) : pair potential
!**             df (out) : first derivative of pair potential
!**            d2f (out) : second derivative of pair potential
!**             lf  (in) : true to compute f
!**            ldf  (in) : true to compute df
!**           ld2f  (in) : true to compute ld2f
!**
      SUBROUTINE V2FL(R,Ity,Jty,F,Df,D2f,Lf,Ldf,Ld2f)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--V2FL418
      DOUBLE PRECISION one , F , Df , D2f , R , p
      INTEGER Ity , Jty , k
 
      DATA one/1.D0/
      LOGICAL Lf , Ldf , Ld2f
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
 
      IF ( Ity<0 .AND. Jty<0 ) THEN
         IF ( Lf ) F = 0.D0
         IF ( Ldf ) Df = 0.D0
         IF ( Ld2f ) D2f = 0.D0
         RETURN
      ELSEIF ( Ity<0 .OR. Jty<0 ) THEN
         IF ( R>INDrad ) STOP 'oops'
         IF ( Lf ) F = CON1*(R-INDrad)**2
         IF ( Ldf ) Df = 2*CON1*(R-INDrad)
         IF ( Ld2f ) D2f = 2*CON1
         RETURN
      ENDIF
 
      p = R/DRAr + 1.0
      k = p
      k = MIN0(k,NRAr-1)
      p = p - k
      p = MIN(p,one)
 
      IF ( R<=1.0D-6 ) THEN
         PRINT * , '**WARNING -- To small r in pair potential' , R
!C--Jun Song Test
         DEBugflag = 1
      ENDIF
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Lf .OR. Ldf .OR. Ld2f ) THEN
         F = ((Z2Rar3(k,Ity,Jty)*p+Z2Rar2(k,Ity,Jty))&
     &       *p+Z2Rar1(k,Ity,Jty))*p + Z2Rar(k,Ity,Jty)
         F = F/R
      ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Ldf .OR. Ld2f ) THEN
         Df = (Z2Rar6(k,Ity,Jty)*p+Z2Rar5(k,Ity,Jty))&
     &        *p + Z2Rar4(k,Ity,Jty)
         Df = (Df-F)/R
      ENDIF
 
      IF ( Ld2f ) THEN
         D2f = Z2Rar7(k,Ity,Jty)&
     &         + (Z2Rar7(k+1,Ity,Jty)-Z2Rar7(k,Ity,Jty))*p
         D2f = (D2f-2.*Df)/R
      ENDIF
 
      END SUBROUTINE V2FL
 
 
!**---------------------------------------------------------------------
!** uufl - gives embedding energy
!**
!**     Parameters :
!**             rh  (in) : electron density
!**            isp  (in) : species of the atom
!**              f (out) : electron density
!**             df (out) : first derivative of electron density
!**            d2f (out) : second derivative of electron density
!**             lf  (in) : true to compute f
!**            ldf  (in) : true to compute df
!**           ld2f  (in) : true to compute ld2f
!**
      SUBROUTINE UUFL(Rh,Isp,F,Df,D2f,Lf,Ldf,Ld2f)
      IMPLICIT NONE
!*--UUFL489
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION D2f , Df , F , one , p , Rh
      INTEGER Isp , k
!*** End of declarations inserted by SPAG
      DATA one/1.0/
      LOGICAL Lf , Ldf , Ld2f
 
      IF ( Isp<0 ) THEN
         IF ( Lf ) F = 0.D0
         IF ( Ldf ) Df = 0.D0
         IF ( Ld2f ) D2f = 0.D0
         RETURN
      ENDIF
 
      p = Rh/DRHoar + 1.0
      k = p
      k = MIN0(k,NRAr-1)
      p = p - k
      p = MIN(p,one)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Lf ) F = ((FRHoar3(k,Isp)*p+FRHoar2(k,Isp))*p+FRHoar1(k,Isp))&
     &              *p + FRHoar(k,Isp)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Ldf ) Df = (FRHoar6(k,Isp)*p+FRHoar5(k,Isp))&
     &                *p + FRHoar4(k,Isp)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( Ld2f ) D2f = FRHoar7(k,Isp)&
     &                  + (FRHoar7(k+1,Isp)-FRHoar7(k,Isp))*p
 
      END SUBROUTINE UUFL
 
!**---------------------------------------------------------------------
!** NumAtomSpec() : get the number of atomic species
!**
      INTEGER FUNCTION NUMATOMSPEC()
      IMPLICIT NONE
!*--NUMATOMSPEC526
 
      NUMATOMSPEC = NTYpes
 
      END FUNCTION NUMATOMSPEC
 
!**--------------------------------------------------------------------
!** ReadPoten - Read the inter atomic potentials from specified file
!**
!**     Parameters :-
!**
!**     Algorithm :-
!**                 Whatever is there in the dynamo code!!
!**
!**
      SUBROUTINE READPOTEN()
      USE MOD_FILE
      IMPLICIT NONE
      
!*--READPOTEN545
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION alat , blat , cof1 , cof2 , cof3 , cof4 , &
     &                 conmas , drhoin , drin , fheader ,  &
     &                 frhoin , p , r , rcut , rcutall , rhomax , &
     &                 rhomaxi , rhorin , rhotmp
      DOUBLE PRECISION rmax , rmaxi , two , zrin , zrtemp
      INTEGER i , i1 , i2 , ipinter , iunit , j , jr , jrho , k , &
     &        nrhoin , nrin , numelmts
!*** End of declarations inserted by SPAG
      DIMENSION frhoin(NGRID,NELMAX) , drhoin(NELMAX) , nrhoin(NELMAX)
      DIMENSION rhorin(NGRID,NELMAX) , zrin(NGRID,NELMAX) , drin(NELMAX)&
     &          , nrin(NELMAX)
      DIMENSION zrtemp(NGRID,NELMAX)
      DIMENSION fheader(10,NELMAX) , rcut(NELMAX) , blat(NELMAX)
      CHARACTER*8 llat(NELMAX) , latty
      CHARACTER*80 sheader(3)
      DATA conmas/1.0365E-4/
      DATA two/2.0/
      LOGICAL setflag
      CHARACTER*80 funcfl(100) , setfl
      DATA funcfl/100*'none'/ , setfl/'none'/
 
 
 
!--Read in the function file name
      ipinter = 0
      READ (input_file_unit,*) setflag
 
      IF ( .NOT.setflag ) THEN
         READ (input_file_unit,*) numelmts
         IF ( numelmts>NELMAX ) THEN
            PRINT * , '**ERROR : Too many func files'
            STOP
         ENDIF
         DO i = 1 , numelmts
            READ (input_file_unit,99019) funcfl(i)
         ENDDO
         i = 0
         DO
            i = i + 1
            IF ( i>NELMAX ) EXIT
            IF ( funcfl(i)(1:4)=='NONE' .OR. funcfl(i)(1:4)=='none' )&
     &           EXIT
            IF ( .NOT.FILEEXISTS(funcfl(i),.TRUE.) ) STOP
            CALL IOFILE(funcfl(i),'formatted  ',iunit,.TRUE.)
            READ (iunit,99001) (fheader(j,i),j=1,10)
99001       FORMAT (10A8)
            READ (iunit,99020) IELement(i) , AMAss(i) , blat(i) , &
     &                         llat(i)
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C--Jun Song: use unformated input for more flexibility
            READ (iunit,*) nrhoin(i) , drhoin(i) , nrin(i) , drin(i) , &
     &                     rcut(i)
!c9901    format(i5,e24.16,i5,2e24.16)
!
!  assume that z(r) and rho(r) grids coincide
!
            READ (iunit,*) (frhoin(j,i),j=1,nrhoin(i))
            READ (iunit,*) (zrin(j,i),j=1,nrin(i))
            READ (iunit,*) (rhorin(j,i),j=1,nrin(i))
!C9902    format(5e24.16)
!
! close the interaction file
!
            CLOSE (iunit)
         ENDDO
         NTYpes = i - 1
         IF ( NTYpes>NELMAX ) THEN
            WRITE (6,*) 'error: number of types greater than nelmax'
            STOP
         ENDIF
!
! determine common grid spacings and number
!
!
! take largest grid spacing
! take smallest maximum
!
         DRAd = drin(1)
         DRHo = drhoin(1)
         rmax = (nrin(1)-1)*drin(1)
         rhomax = (nrhoin(1)-1)*drhoin(1)
         DO i1 = 2 , NTYpes
            DRAd = MAX(DRAd,drin(i1))
            DRHo = MAX(DRHo,drhoin(i1))
            rmaxi = (nrin(i1)-1)*drin(i1)
            rhomaxi = (nrhoin(i1)-1)*drhoin(i1)
            rmax = MAX(rmax,rmaxi)
            rhomax = MAX(rhomax,rhomaxi)
         ENDDO
         NR = NINT(rmax/DRAd)
         NRHo = NINT(rhomax/DRHo)
!
! set up the z(r) and rho(r) grids
!
         DO i1 = 1 , NTYpes
            DO j = 1 , NR
               r = (j-1)*DRAd
!
!  do four-point lagrange interpolation
!
               p = r/drin(i1) + 1.0
               k = p
               k = MIN0(k,nrin(i1)-2)
               k = MAX0(k,2)
               p = p - k
!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and rhor = last value of rhor
               p = MIN(p,two)
               cof1 = -0.166666667*p*(p-1.)*(p-2.)
               cof2 = 0.5*(p**2-1.)*(p-2.)
               cof3 = -0.5*p*(p+1.)*(p-2.)
               cof4 = 0.166666667*p*(p**2-1.)
               RHOr(j,i1) = cof1*rhorin(k-1,i1) + cof2*rhorin(k,i1)&
     &                      + cof3*rhorin(k+1,i1) + cof4*rhorin(k+2,i1)
               zrtemp(j,i1) = cof1*zrin(k-1,i1) + cof2*zrin(k,i1)&
     &                        + cof3*zrin(k+1,i1) + cof4*zrin(k+2,i1)
            ENDDO
         ENDDO
!
! set up the f(rho) grid
!
         DO i1 = 1 , NTYpes
            DO j = 1 , NRHo
               r = (j-1)*DRHo
!
!  do four-point lagrange interpolation
!
               p = r/drhoin(i1) + 1.0
               k = p
               k = MIN0(k,nrhoin(i1)-2)
               k = MAX0(k,2)
               p = p - k
!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and rhor = last value of rhor
               p = MIN(p,two)
               cof1 = -0.166666667*p*(p-1.)*(p-2.)
               cof2 = 0.5*(p**2-1.)*(p-2.)
               cof3 = -0.5*p*(p+1.)*(p-2.)
               cof4 = 0.166666667*p*(p**2-1.)
               FRHo(j,i1) = cof1*frhoin(k-1,i1) + cof2*frhoin(k,i1)&
     &                      + cof3*frhoin(k+1,i1) + cof4*frhoin(k+2,i1)
            ENDDO
         ENDDO
!
! set up the z2 grid (zi*zj)
!
         DO i1 = 1 , NTYpes
            DO i2 = 1 , NTYpes
               DO j = 1 , NR
                  Z2R(j,i1,i2) = 27.2*0.529*zrtemp(j,i1)*zrtemp(j,i2)
               ENDDO
            ENDDO
         ENDDO
         RCUtsq = 0.0
         DO i = 1 , NTYpes
            RCUtsq = MAX(RCUtsq,rcut(i))
         ENDDO
         IF ( RCUtsq==0.0 ) RCUtsq = 5.0
!       print out types
         WRITE (6,99021)
         WRITE (6,99022) NTYpes
         WRITE (6,99023)
         WRITE (6,99024) (i,IELement(i),AMAss(i),blat(i),llat(i),i=1,&
     &                   NTYpes)
         DO i = 1 , NTYpes
!       print out header
            WRITE (6,99002)
99002       FORMAT ('   type  file name  header',/,&
     &              '   ----  ---------  ---------------')
            WRITE (6,99003) i , funcfl(i) , (fheader(j,i),j=1,10)
99003       FORMAT (1x,i4,5x,a80,3x,10A8)
            WRITE (6,99025) rcut(i)
            IF ( ipinter>0 ) THEN
               WRITE (6,99004)
99004          FORMAT (&
     &           '    r          z        rho              rho        f'&
     &           ,/,&
     &        ' --------  --------  ----------       --------  --------'&
     &        )
               DO j = 1 , NR
                  jr = MIN0(j,NR)
                  jrho = MIN0(j,NRHo)
                  r = (jr-1)*DRAd
                  rhotmp = (jrho-1)*DRHo
                  WRITE (6,99005) r , zrtemp(jr,i) , RHOr(jr,i) , &
     &                            rhotmp , FRHo(jrho,i)
99005             FORMAT (1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7)
               ENDDO
            ENDIF
         ENDDO
 
 
      ELSE
         READ (input_file_unit,99019) setfl
         IF ( .NOT.FILEEXISTS(setfl,.TRUE.) ) STOP
         CALL IOFILE(setfl,'formatted  ',iunit,.TRUE.)
         READ (iunit,'(a80)') sheader(1)
         READ (iunit,'(a80)') sheader(2)
         READ (iunit,'(a80)') sheader(3)
         READ (iunit,99006) NTYpes
99006    FORMAT (i5)
         IF ( NTYpes>NELMAX ) THEN
            WRITE (6,*) 'error: number of types greater than nelmax'
            STOP
         ENDIF
         READ (iunit,*) NRHo , DRHo , NR , DRAd , rcutall
         RCUtsq = rcutall
         DO i = 1 , NTYpes
            READ (iunit,99020) IELement(i) , AMAss(i) , blat(i) , &
     &                         llat(i)
            READ (iunit,*) (FRHo(j,i),j=1,NRHo)
            READ (iunit,*) (RHOr(j,i),j=1,NR)
         ENDDO
         DO i1 = 1 , NTYpes
            DO i2 = 1 , i1
               READ (iunit,*) (Z2R(j,i1,i2),j=1,NR)
            ENDDO
         ENDDO
         DO i1 = 1 , NTYpes
            DO i2 = i1 + 1 , NTYpes
               DO j = 1 , NR
                  Z2R(j,i1,i2) = Z2R(j,i2,i1)
               ENDDO
            ENDDO
         ENDDO
         CLOSE (iunit)
         WRITE (6,99021)
         WRITE (6,99022) NTYpes
         WRITE (6,99023)
         WRITE (6,99024) (i,IELement(i),AMAss(i),blat(i),llat(i),i=1,&
     &                   NTYpes)
         WRITE (6,99007) setfl
99007    FORMAT (1x,'    file name  ',a80)
!  print out header
         WRITE (6,99008)
99008    FORMAT ('   header',/,'   ________________')
         WRITE (6,'(4x,a80)') (sheader(i),i=1,3)
99009    FORMAT (4x,10A8)
         WRITE (6,99025) RCUtsq
         IF ( ipinter>0 ) THEN
            DO i = 1 , NTYpes
               WRITE (6,99026) i
               WRITE (6,99010)
99010          FORMAT ('      r          rho              rho        f',&
     &                 /,&
     &                '   --------  ----------       --------  --------'&
     &                )
               DO j = 1 , MAX0(NR,NRHo)
                  jr = MIN0(j,NR)
                  jrho = MIN0(j,NRHo)
                  r = (jr-1)*DRAd
                  rhotmp = (jrho-1)*DRHo
                  WRITE (6,99011) r , RHOr(jr,i) , rhotmp , FRHo(jrho,i)
99011             FORMAT (1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7)
               ENDDO
            ENDDO
            DO i1 = 1 , NTYpes
               DO i2 = 1 , i1
                  WRITE (6,99027) i1 , i2
                  WRITE (6,99012)
99012             FORMAT ('      r            z**2',/,&
     &                    '   _______      _____________')
                  DO j = 1 , NR
                     r = (j-1)*DRAd
                     WRITE (6,99013) r , Z2R(j,i1,i2)
99013                FORMAT (1x,g15.7,1x,g15.7)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
 
 
 
 
!
!       this is the only place where the mass is in amu
!       here we convert to eV-psec**2/angstrom**2
!       this is the unit used throughout the program
!       restart assumes this mass unit
!
      DO i = 1 , NTYpes
         AMAss(i) = conmas*AMAss(i)
      ENDDO
      WRITE (6,99014) RCUtsq
99014 FORMAT ('  use this cut-off distance: ',g15.5)
      SQRtrcutsq = RCUtsq
      RCUtsq = RCUtsq**2
!
!       set the lattice constant to that for type 1 by default
      alat = blat(1)
      latty = llat(1)
!
!
!  now set up the dense grids
!
      NRHoar = NRHo
      DRHoar = DRHo
      NRAr = NR
      DRAr = DRAd
!
!  compute the value and slope at the computational grid points
!  the methodology is as follows.
!  At each point, the slope is estimated from a 5-point Lagrange interpo
!  of the data.
!  Between each pair of points, a cubic polynomial is used which is fit 
!  value and slope at the two points.
!  This yields an interpolation which gives continuous value and first d
!  without introducing the long-range effects of glitches in the data th
!  results from using splines.  (i.e. the procedure is local)
!
      DO i = 1 , NTYpes
         DO j = 1 , NRHoar
            FRHoar(j,i) = FRHo(j,i)
         ENDDO
         DO j = 1 , NRAr
            RHOrar(j,i) = RHOr(j,i)
         ENDDO
         FRHoar1(1,i) = FRHoar(2,i) - FRHoar(1,i)
         FRHoar1(2,i) = 0.5*(FRHoar(3,i)-FRHoar(1,i))
         FRHoar1(NRHoar-1,i) = 0.5*(FRHoar(NRHoar,i)-FRHoar(NRHoar-2,i))
         FRHoar1(NRHoar,i) = FRHoar(NRHoar,i) - FRHoar(NRHoar-1,i)
         FRHoar7(1,i) = (FRHoar(3,i)+FRHoar(1,i)-2.*FRHoar(2,i))&
     &                  /(DRHoar*DRHoar)
         FRHoar7(2,i) = (FRHoar(3,i)+FRHoar(1,i)-2.*FRHoar(2,i))&
     &                  /(DRHoar*DRHoar)
         FRHoar7(NRHoar-1,i) = (FRHoar(NRHoar-2,i)+FRHoar(NRHoar,i)-2.*&
     &                         FRHoar(NRHoar-1,i))/(DRHoar*DRHoar)
         FRHoar7(NRHoar,i) = (FRHoar(NRHoar-2,i)+FRHoar(NRHoar,i)-2.*&
     &                       FRHoar(NRHoar-1,i))/(DRHoar*DRHoar)
 
         DO j = 3 , NRHoar - 2
            FRHoar1(j,i) = ((FRHoar(j-2,i)-FRHoar(j+2,i))+8.*(FRHoar(j+1&
     &                     ,i)-FRHoar(j-1,i)))/12.
            FRHoar7(j,i) = (16.*(FRHoar(j-1,i)+FRHoar(j+1,i))-(FRHoar(j+&
     &                     2,i)+FRHoar(j-2,i))-30.*FRHoar(j,i))&
     &                     /(12.*DRHoar*DRHoar)
         ENDDO
         DO j = 1 , NRHoar - 1
            FRHoar2(j,i) = 3.*(FRHoar(j+1,i)-FRHoar(j,i))&
     &                     - 2.*FRHoar1(j,i) - FRHoar1(j+1,i)
            FRHoar3(j,i) = FRHoar1(j,i) + FRHoar1(j+1,i)&
     &                     - 2.*(FRHoar(j+1,i)-FRHoar(j,i))
         ENDDO
         FRHoar2(NRHoar,i) = 0.
         FRHoar3(NRHoar,i) = 0.
         DO j = 1 , NRHoar
            FRHoar4(j,i) = FRHoar1(j,i)/DRHoar
            FRHoar5(j,i) = 2.*FRHoar2(j,i)/DRHoar
            FRHoar6(j,i) = 3.*FRHoar3(j,i)/DRHoar
         ENDDO
         RHOrar1(1,i) = RHOrar(2,i) - RHOrar(1,i)
         RHOrar1(2,i) = 0.5*(RHOrar(3,i)-RHOrar(1,i))
         RHOrar1(NRAr-1,i) = 0.5*(RHOrar(NRAr,i)-RHOrar(NRAr-2,i))
         RHOrar1(NRAr,i) = 0.
         RHOrar7(1,i) = (RHOrar(3,i)+RHOrar(1,i)-2.*RHOrar(2,i))&
     &                  /(DRAr*DRAr)
         RHOrar7(2,i) = (RHOrar(3,i)+RHOrar(1,i)-2.*RHOrar(2,i))&
     &                  /(DRAr*DRAr)
         RHOrar7(NRAr-1,i) = (RHOrar(NRAr-2,i)+RHOrar(NRAr,i)-2.*RHOrar(&
     &                       NRAr-1,i))/(DRAr*DRAr)
         RHOrar7(NRAr,i) = (RHOrar(NRAr-2,i)+RHOrar(NRAr,i)-2.*RHOrar(&
     &                     NRAr-1,i))/(DRAr*DRAr)
         DO j = 3 , NRAr - 2
            RHOrar1(j,i) = ((RHOrar(j-2,i)-RHOrar(j+2,i))+8.*(RHOrar(j+1&
     &                     ,i)-RHOrar(j-1,i)))/12.
            RHOrar7(j,i) = (16.*(RHOrar(j-1,i)+RHOrar(j+1,i))-(RHOrar(j+&
     &                     2,i)+RHOrar(j-2,i))-30.*RHOrar(j,i))&
     &                     /(12.*DRAr*DRAr)
         ENDDO
 
         DO j = 1 , NRAr - 1
            RHOrar2(j,i) = 3.*(RHOrar(j+1,i)-RHOrar(j,i))&
     &                     - 2.*RHOrar1(j,i) - RHOrar1(j+1,i)
            RHOrar3(j,i) = RHOrar1(j,i) + RHOrar1(j+1,i)&
     &                     - 2.*(RHOrar(j+1,i)-RHOrar(j,i))
         ENDDO
         RHOrar2(NRAr,i) = 0.
         RHOrar3(NRAr,i) = 0.
         DO j = 1 , NRAr
            RHOrar4(j,i) = RHOrar1(j,i)/DRAr
            RHOrar5(j,i) = 2.*RHOrar2(j,i)/DRAr
            RHOrar6(j,i) = 3.*RHOrar3(j,i)/DRAr
         ENDDO
         i1 = i
         DO i2 = 1 , NTYpes
            DO j = 1 , NRAr
               Z2Rar(j,i1,i2) = Z2R(j,i1,i2)
            ENDDO
            Z2Rar1(1,i1,i2) = Z2Rar(2,i1,i2) - Z2Rar(1,i1,i2)
            Z2Rar1(2,i1,i2) = 0.5*(Z2Rar(3,i1,i2)-Z2Rar(1,i1,i2))
            Z2Rar1(NRAr-1,i1,i2) = 0.5*(Z2Rar(NRAr,i1,i2)-Z2Rar(NRAr-2,&
     &                             i1,i2))
            Z2Rar1(NRAr,i1,i2) = 0.
 
            Z2Rar7(1,i1,i2) = (Z2Rar(3,i1,i2)+Z2Rar(1,i1,i2)-2.*Z2Rar(2,&
     &                        i1,i2))/(DRAr*DRAr)
            Z2Rar7(2,i1,i2) = (Z2Rar(3,i1,i2)+Z2Rar(1,i1,i2)-2.*Z2Rar(2,&
     &                        i1,i2))/(DRAr*DRAr)
            Z2Rar7(NRAr-1,i1,i2) = (Z2Rar(NRAr-2,i1,i2)+Z2Rar(NRAr,i1,i2&
     &                             )-2.*Z2Rar(NRAr-1,i1,i2))/(DRAr*DRAr)
            Z2Rar7(NRAr,i1,i2) = (Z2Rar(NRAr-2,i1,i2)+Z2Rar(NRAr,i1,i2)&
     &                           -2.*Z2Rar(NRAr-1,i1,i2))/(DRAr*DRAr)
 
            DO j = 3 , NRAr - 2
               Z2Rar1(j,i1,i2) = ((Z2Rar(j-2,i1,i2)-Z2Rar(j+2,i1,i2))+8.&
     &                           *(Z2Rar(j+1,i1,i2)-Z2Rar(j-1,i1,i2)))&
     &                           /12.
               Z2Rar7(j,i1,i2) = (16.*(Z2Rar(j-1,i1,i2)+Z2Rar(j+1,i1,i2)&
     &                           )-(Z2Rar(j+2,i1,i2)+Z2Rar(j-2,i1,i2))&
     &                           -30.*Z2Rar(j,i1,i2))/(12.*DRAr*DRAr)
            ENDDO
 
            DO j = 1 , NRAr - 1
               Z2Rar2(j,i1,i2) = 3.*(Z2Rar(j+1,i1,i2)-Z2Rar(j,i1,i2))&
     &                           - 2.*Z2Rar1(j,i1,i2)&
     &                           - Z2Rar1(j+1,i1,i2)
               Z2Rar3(j,i1,i2) = Z2Rar1(j,i1,i2) + Z2Rar1(j+1,i1,i2)&
     &                           - 2.*(Z2Rar(j+1,i1,i2)-Z2Rar(j,i1,i2))
            ENDDO
            Z2Rar2(NRAr,i1,i2) = 0.
            Z2Rar3(NRAr,i1,i2) = 0.
            DO j = 1 , NRAr
               Z2Rar4(j,i1,i2) = Z2Rar1(j,i1,i2)/DRAr
               Z2Rar5(j,i1,i2) = 2.*Z2Rar2(j,i1,i2)/DRAr
               Z2Rar6(j,i1,i2) = 3.*Z2Rar3(j,i1,i2)/DRAr
            ENDDO
         ENDDO
      ENDDO
 
!
!   output dense grids
!
      IF ( ipinter>1 ) THEN
         DO i = 1 , NTYpes
            WRITE (6,99026) i
            WRITE (6,99015)
99015       FORMAT ('      r          rho         rhop',&
     &              '        rho        f         fp',/,&
     &              '   --------  ----------    --------',&
     &              '      --------  --------  ----------')
            DO j = 1 , MAX0(NRAr,NRHoar)
               jr = MIN0(j,NRAr)
               jrho = MIN0(j,NRHoar)
               r = (jr-1)*DRAr
               rhotmp = (jrho-1)*DRHoar
               WRITE (6,99016) r , RHOrar(jr,i) , RHOrar1(jr,i) , &
     &                         rhotmp , FRHoar(jrho,i) , FRHoar1(jrho,i)
99016          FORMAT (1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7,1x,&
     &                 g15.7)
            ENDDO
         ENDDO
         DO i1 = 1 , NTYpes
            DO i2 = 1 , i1
               WRITE (6,99027) i1 , i2
               WRITE (6,99017)
99017          FORMAT ('      r            z**2          z**2p',/,&
     &                 '   _______      _____________  ____________')
               DO j = 1 , NRAr
                  r = (j-1)*DRAr
                  WRITE (6,99018) r , Z2Rar(j,i1,i2) , Z2Rar1(j,i1,i2)
99018             FORMAT (1x,g15.7,1x,g15.7,1x,g15.7)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
99019 FORMAT (a)
99020 FORMAT (i5,2G15.5,a8)
99021 FORMAT (/,/' ******   interactions defined  ')
99022 FORMAT (1x,i5,' particle types')
99023 FORMAT ('   type  element       amass              alat        ',&
     &        ' lattype',/,&
     &        '   ----  --------  --------------  ----------------    ',&
     &        '----------')
99024 FORMAT (1x,i4,i9,g15.5,5x,g15.5,10x,a8)
99025 FORMAT ('   cut-off distance =',g15.5)
99026 FORMAT (' type ',i5,/,' _____________')
99027 FORMAT (' types ',i5,i5,/,' ______________________')
      END SUBROUTINE READPOTEN
      END MODULE MOD_POTEN
!*==cutoffradius.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
      DOUBLE PRECISION FUNCTION CUTOFFRADIUS(I)
      USE MOD_POTEN
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--CUTOFFRADIUS1033
      INTEGER I
 
      IF ( I>0 ) THEN
         CUTOFFRADIUS = SQRtrcutsq
      ELSE
         CUTOFFRADIUS = INDrad
      ENDIF
 
      END FUNCTION CUTOFFRADIUS
!*==cutoffr2.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!**---------------------------------------------------------------------
!** CutoffR2() : Get the cut off radius
!**
      DOUBLE PRECISION FUNCTION CUTOFFR2(I)
      USE MOD_POTEN
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--CUTOFFR21052
      INTEGER I
 
      IF ( I>0 ) THEN
         CUTOFFR2 = RCUtsq
      ELSE
         CUTOFFR2 = INDradsq
      ENDIF
      END FUNCTION CUTOFFR2
 
