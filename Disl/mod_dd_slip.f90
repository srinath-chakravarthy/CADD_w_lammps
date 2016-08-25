!!> --- This module contains the entire DD methodology
!! --- Generate slip planes
!! --- Generate sources to a given density
!! --- Generate obstacles to a given spacing lobs
!! --- Calculate the P-K force (direct or multipole) -- currently
!!     only direct
!! --- Calculate the displacement at any location
!! --- Calculate the stress at any given point
!! --- Handling of pinned dislocations

MODULE MOD_DD_SLIP
      IMPLICIT NONE
      INTEGER, PARAMETER :: MXNSLP = 50000 !< @var   mxnslp --> max. number of slip planes
      INTEGER, PARAMETER :: MXNSYS = 12 !<    mxnsys --> max. number of slip systems (default = 3)

      INTEGER  :: nsys , nslp
      !!<    nsys   --> Current number of slip systems (default = 3)
      !!<    nslp   --> Current number of slip planes
      DOUBLE PRECISION :: dtime_dd , TNUC , tincr , TINCR_SAV
      PARAMETER (TNUC=1.D-10)
      PARAMETER (TINCR_SAV=1.D-11)
      !!<    tincr ---> DD timer (default = 1.d-11)
      !!<    tnuc --->  Nucleation timer (default = 1.d-9).
      !!<    dtime_dd ---> running DD time step
      !!     Currently tnuc and tincr needs to be changed in the code

      INTEGER :: locphi(MXNSLP) !< @var slipplane identification

      INTEGER :: nntot , tot_source , n_active_slip
      !! 1/100 of all slip pla
      DOUBLE PRECISION :: process_xmin , process_xmax , process_ymin , process_ymax, process_area
      DOUBLE PRECISION :: source_den , avg_source_str , sd_source_str
      DOUBLE PRECISION :: atom_xmin , atom_xmax , atom_ymin , atom_ymax
      DOUBLE PRECISION :: e_dd , mu_dd , nu_dd
      DOUBLE PRECISION :: phisys , phislp(MXNSYS) , sinphi(MXNSYS) , cosphi(MXNSYS)
      DOUBLE PRECISION , ALLOCATABLE :: xslp(:) , yslp(:) , xendslp(:) , yendslp(:)
      INTEGER , ALLOCATABLE :: nshift_slip(:)
      INTEGER , ALLOCATABLE :: elem_slip(:) !!<     For elements in the DB this contains the slip plane number

      INTEGER ::  inslip(3) !< in burgers vectors
      INTEGER, PARAMETER :: SLIPSPC=100
      DOUBLE PRECISION :: dx_slip , dy_slip
      DOUBLE PRECISION :: sdis_out_tol , disl_range_tol
!!$ ----------------------------------------------------------------
!!$    Dislocation Nculeation related variables
!!$ ----------------------------------------------------------------  
      INTEGER, PARAMETER :: MXNNUC = 100

      INTEGER , ALLOCATABLE :: nnuc(:)
      DOUBLE PRECISION :: xle , d , blen , tobs , dstar , tref , htinc , epsx , epss

      DOUBLE PRECISION , ALLOCATABLE :: t_fr(:,:) , xlnuc(:,:) ,  tnlaps(:,:) , taui(:,:) 
      DOUBLE PRECISION , ALLOCATABLE :: snuc(:,:) , rnuc(:,:,:)
                                !$!      double precision, allocatable :: nuc_pk_stress(:,:), nuc_pk_force
      INTEGER , ALLOCATABLE :: elem_source(:,:)
!!$ ---------------------------------------------------------------- 
!!$      Obstacles common block
!!$ ---------------------------------------------------------------- 
      INTEGER :: tot_obs , MXNOBS
      PARAMETER (MXNOBS=100)
      INTEGER , ALLOCATABLE :: nobs(:)
      DOUBLE PRECISION , ALLOCATABLE :: tau_obs(:,:) , sobs(:,:)
      DOUBLE PRECISION :: lobs , lobs_max , lobs_min
      !!     (obstacle spacing Angstroms)

!!$ ----------------------------------------------------------------  
!!$     Dislocations Common block
!!$ ----------------------------------------------------------------
      INTEGER, PARAMETER :: MXNDIS = 20
      DOUBLE PRECISION ,PARAMETER :: VCUTOFF=1E13  ! A/s
      DOUBLE PRECISION, PARAMETER :: BDRAG=5.D-15  ! ev/A^3 s

      INTEGER :: npile , nd , nloop , tot_disl
      INTEGER , ALLOCATABLE :: ndis(:)
      DOUBLE PRECISION , ALLOCATABLE :: sdis(:,:) , b_dd(:,:) , vprev(:,:)
      DOUBLE PRECISION , ALLOCATABLE :: bout(:,:) , sdis_out(:,:) , sdis_out1(:,:)
      DOUBLE PRECISION , ALLOCATABLE :: rdis(:,:,:) , rold(:,:,:)
      DOUBLE PRECISION :: rpile , dlim
      INTEGER , ALLOCATABLE :: iobpin(:,:) , jcnpin(:,:) , idple(:,:) , elem_dis(:,:)

      INTEGER :: tot_size       !> Total size of dislocations + nucleati

      DOUBLE PRECISION :: range_ !> Range for multipole calculations

      TYPE DD
         !!<     Discrete dislocation object
         !!<     Contains all the variables relating to dislocations
         DOUBLE COMPLEX :: XY   !< (x,y) coordinate of dislocation
         DOUBLE PRECISION :: SLIPANGLE
         !!< Angle of slip plane
         DOUBLE PRECISION , DIMENSION(3) :: BURGERS_DD
         !!< Burgers vector of 
         DOUBLE PRECISION , DIMENSION(3) :: STRESS
         !!< Stress on dislocation
         !!<     stress(1) --> s11, stress(2) --> s22, stress(3) --> s12
         DOUBLE COMPLEX :: FORCE , FORCED , DISP , DISPD , VELOCITY , FORCEDG
         !!<     force --> PK force on dislocation
         !!<     forced --> PK force on dislocation due to direct calculation
         !!<     disp ---> displacement due to dislocation
         !!<     forceg ---> gradient of PK force along the slip plane
         DOUBLE COMPLEX :: ALPHA , ALPHA1 , ALPHA2 , ALPHA3 , ALPHA4 , ALPHA5 , ALPHA6 , ALPHA7
         !! --- Gradient Direct terms
         DOUBLE COMPLEX :: ALPHA1G , ALPHA2G , BIALPHAG , BETA1G , BETA2G , BIBETAG
         DOUBLE COMPLEX :: ALPHAD1 , BETAD1 , ALPHAD2 , BETAD2
         DOUBLE COMPLEX :: BETA , BETA1 , BETA2 , BETA3 , BETA4 , BETA5 , BETA6 , BETA7
         DOUBLE COMPLEX :: BIALPHA , BIBETA , PHI1 , PHI2 , PHI11 , PHID1
         !!>    Multipole terms
         !!> --- Gradient terms
         DOUBLE COMPLEX :: PHI1G , PHI2G , PHI11G
         DOUBLE COMPLEX :: DPHI1 , DPHI2 , DPHI11 , DPHID1
         INTEGER :: STYPE , SLIPPLANE_NO
      END TYPE DD




CONTAINS
!!$      SUBROUTINE GEN_SLIP_PLANES1(INTERFACE_ATOMS)
!!$            IMPLICIT NONE
!!$            DOUBLE PRECISION, INTENT(IN) :: INTERFACE_ATOMS(:)
!!$            INTEGER :: I, J, K
!!$            INTEGER :: ISLP, II, JSLP, LI, LJ
!!$            DOUBLE PRECISION :: XI, XJ, m , c, m1, c1, tant
!!$            integer :: num_interface_atoms
!!$            
!!$            do i = 1, nsys
!!$               cosphi(i) = cos(phislp(i))
!!$               sinphi(i) = sin(phislp(i))
!!$            end do
!!$
!!$            num_interface_atoms = size(interface_atoms)
!!$            nslp = num_interface_atoms -1
!!$
!!$            ALLOCATE (XSLp(NSLp),YSLp(NSLp),XENdslp(NSLp),YENdslp(NSLp))
!!$            ALLOCATE (NSHift_slip(NSLp))
!!$            ALLOCATE (ELEm_slip(NSLp))
!!$            XSLp = 0.0D0
!!$            YSLp = 0.0D0
!!$            XENdslp = 0.0D0
!!$            YENdslp = 0.0D0
!!$
!!$            islp = 0
!!$            do i = 1, num_interface_atoms-1
!!$               xstart = (interface_atoms(1,i) + interface_atoms(1,i+1))/2.0d0
!!$               ystart = (interface_atoms(2,i) + interface_atoms(1,2+1))/2.0d0
!!$               do isys = 1, nsys
!!$                  islp = islp + 1
!!$                  locphi(islp) = isys
!!$                  xslp(islp) = xstart
!!$                  yslp(islp) = ystart
!!$                  tant = tan(phisys(locphi(islp)))
!!$  !--- Intesection to the box formed by xmax, xmin, ymax, ymin
!!$                  !! eqn of slip plane
!!$                  if (abs(sinphi(islp)) < 1.d-6) then
!!$                     xend1 = xstart + (ymax-ystart)/tant)
!!$                     xend2 = xstart + (ymin-ystart)/tant)
!!$                  end if
!!$                  if (abs(xend1) > abs(xend2)) then
!!$                     xend = xend1
!!$                  else
!!$                     xend = xend2
!!$                  end if
!!$                  if (xend > xmax) then
!!$                     xend = xmax
!!$                     yend = ystart + (xend - xstart)*tant
!!$                  else if (xend < xmin) then
!!$                     xend = xmin
!!$                     yend = ystart + (xend - xstart)*tant
!!$                  else
!!$                     yend = ymin
!!$                  end if
!!$                  
!!$               end do
!!$            end do
!!$      END SUBROUTINE GEN_SLIP_PLANES1



      SUBROUTINE GEN_SLIP_PLANES(Xmin,Xmax,Ymin,Ymax,Atomxmin,Atomxmax,&
           &                           Atomymin,Atomymax,Slip_angle,Burg_len,&
           &                           Xslip_start,Yslip_start,Dxslip,Dyslip,&
           &                           Pad_width)
            IMPLICIT NONE
!!$*--GEN_SLIP_PLANES152

            DOUBLE PRECISION :: Xmin , Xmax , Ymin , Ymax , Burg_len
            DOUBLE PRECISION :: Atomxmin , Atomymin , Atomxmax , Atomymax
            DOUBLE PRECISION :: Xslip_start , Yslip_start , Dxslip , Dyslip
            DOUBLE PRECISION :: Slip_angle(3)
!!$     Local Variables
            INTEGER :: i , j , k , l , islp , ii , nshift , jslp , li , lj
            DOUBLE PRECISION :: xstart , ystart , xend , yend , tmp
            DOUBLE PRECISION :: xi , xj
            DOUBLE PRECISION :: Pad_width



!!$      disl_range_tol = pad_width
!!$      sdis_out_tol = pad_width/2.0d0

            DISl_range_tol = 6.0D0
            SDIs_out_tol = DISl_range_tol

            DTIme_dd = 0.0D0
            NSYs = 3
            DX_slip = Dxslip
            DY_slip = Dyslip
            PRINT * , '---------------------------------------'
            PRINT * , 'Generating slip planes'
            PRINT * , Xmin , Xmax , Atomxmin , Atomxmax
            PRINT * , Ymin , Ymax , Atomymin , Atomymax
            PRINT * , Slip_angle(1) , Slip_angle(2) , Slip_angle(3)
            PRINT * , Burg_len
            PRINT * , Dxslip , Dyslip
            PRINT * , Xslip_start , Yslip_start
            PRINT * , '---------------------------------------'
!!$     Generate bottom and top halves separately ?
!!$     Currently it is set up so that both are generated separately
!!$     to account for the crack face

            PROcess_xmin = 0.75*Xmin
            PROcess_xmax = 0.75*Xmax
            PROcess_ymin = 0.75*Ymin
            PROcess_ymax = 0.75*Ymax

            Xmin = 0.75*Xmin
            Xmax = 0.75*Xmax
            Ymin = 0.75*Ymin
            Ymax = 0.75*Ymax

            PROcess_area = (PROcess_xmax-PROcess_xmin)&
                 &               *(PROcess_ymax-PROcess_ymin)
!!$ Angstrom^2

            ATOm_xmin = Atomxmin
            ATOm_xmax = Atomxmax
            ATOm_ymin = Atomymin
            ATOm_ymax = Atomymax

            BLEn = Burg_len
            XLE = 6.D0*BLEn
            NSLp = 0
            DO i = 1 , MXNSYS
               PHIslp(i) = Slip_angle(i)
               COSphi(i) = COS(Slip_angle(i))
               SINphi(i) = SIN(Slip_angle(i))
               PRINT * , 'Slip angle' , PHIslp(i)
               IF ( i/=3 ) THEN
                  tmp = Ymax/ABS(TAN(Slip_angle(i)))
                  INSlip(i) = 2*INT((Xmax-Xmin)/Dxslip)
               ELSE
                  INSlip(i) = 2*INT((Ymax)/(Dyslip))
               ENDIF
               NSLp = NSLp + INSlip(i)
               PRINT * , 'No. of slip systems =' , INSlip(i)
            ENDDO
            PRINT * , 'Total no. of slip systems' , NSLp

            ALLOCATE (XSLp(NSLp),YSLp(NSLp),XENdslp(NSLp),YENdslp(NSLp))
            ALLOCATE (NSHift_slip(NSLp))
            ALLOCATE (ELEm_slip(NSLp))
            XSLp = 0.0D0
            YSLp = 0.0D0
            XENdslp = 0.0D0
            YENdslp = 0.0D0


!!$     Now actually generate slip_planes through the whole sample
!!$     Count the actual no. of slip planes and compare with nslp
            DO i = 1 , 3
               IF ( i/=3 ) THEN
                  IF ( i==1 ) THEN
                     xstart = Xmin + Dxslip
!!$     ystart is hard coded to be 0.0 for crack
                     ystart = 0.01D-3
                     islp = 0
                     DO WHILE ( xstart<Xmax )
                        xend = xstart + (Ymax-ystart)/TAN(Slip_angle(i))
                        IF ( xend>Xmax ) THEN
                           xend = Xmax
                           yend = ystart + (xend-xstart)*TAN(Slip_angle(i))
                        ELSE
                           yend = Ymax
                        ENDIF
                        islp = islp + 1
                        LOCphi(islp) = 1
                        XSLp(islp) = xstart
                        XENdslp(islp) = xend
                        YSLp(islp) = ystart
                        YENdslp(islp) = yend
                        xstart = xstart + Dxslip
!!$                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
!!$$$$                  if (mod(islp,25) .eq. 0) then
!!$$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
!!$$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp
!!$$$$     $                    yendslp(islp)
!!$$$$                  endif
                     ENDDO
                     PRINT * , 'End of first slip =' , islp
                     INSlip(i) = islp
                  ENDIF
                  IF ( i==2 ) THEN
                     xstart = Xmin + Dxslip
!!$     ystart is hard coded to be 0.0 for crack
                     ystart = 0.01D-3
!!$               xend =xstart+(ymax-ystart)/tan(slip_angle(i))
!!$               print *, xstart, xend
!!$islp = 0
                     j = 0
                     DO WHILE ( xstart<Xmax )
                        j = j + 1
                        xend = xstart + (Ymax-ystart)/TAN(Slip_angle(i))
                        IF ( xend<Xmin ) THEN
                           xend = Xmin
                           yend = ystart + (xend-xstart)*TAN(Slip_angle(i))
                        ELSE
                           yend = Ymax
                        ENDIF
                        islp = islp + 1
                        LOCphi(islp) = 2
                        XSLp(islp) = xstart
                        XENdslp(islp) = xend
                        YSLp(islp) = ystart
                        YENdslp(islp) = yend
                        xstart = xstart + Dxslip
!!$                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
!!$$$$                  if (mod(islp,25) .eq. 0) then
!!$$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
!!$$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp
!!$$$$     $                    yendslp(islp)
!!$$$$                  endif

!!$print *, islp, xslp(islp), xendslp(islp)
                     ENDDO
                     PRINT * , 'End of second slip =' , islp
                     INSlip(i) = islp
                  ENDIF
               ELSE
                  xstart = Xmin
                  xend = Xmax
                  ystart = Dyslip
!!$$$$            nqslp = islp
                  j = 0
                  DO WHILE ( ystart<Ymax )
                     islp = islp + 1
                     j = j + 1
                     LOCphi(islp) = 3
                     XSLp(islp) = xstart
                     XENdslp(islp) = xend
                     YSLp(islp) = ystart
                     YENdslp(islp) = ystart
                     ystart = ystart + Dyslip
!!$$$$                  if (mod(islp,25) .eq. 0) then
!!$$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
!!$$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp
!!$$$$     $                    yendslp(islp)
!!$$$$                  endif

                  ENDDO

                  PRINT * , 'End of third slip =' , islp , j
                  INSlip(i) = islp
               ENDIF
            ENDDO
            PRINT * , INSlip(1) , INSlip(2) , INSlip(3)
            nshift = INSlip(1) + INSlip(2) + INSlip(3)
!!$     Now mirror slip planes
!!$----------------------------------------------------------------------

            DO i = 1 , 3
               IF ( i/=3 ) THEN
                  IF ( i==1 ) THEN
                     xstart = Xmin + Dxslip
!!$     ystart is hard coded to be 0.0 for crack
                     ystart = -0.01D-3
                     DO WHILE ( xstart<Xmax )
                        xend = xstart + (Ymin-ystart)/TAN(Slip_angle(i))
                        IF ( xend<Xmin ) THEN
                           xend = Xmin
                           yend = ystart + (xend-xstart)*TAN(Slip_angle(i))
                        ELSE
                           yend = Ymin
                        ENDIF
                        islp = islp + 1
                        LOCphi(islp) = 1
                        XSLp(islp) = xstart
                        XENdslp(islp) = xend
                        YSLp(islp) = ystart
                        YENdslp(islp) = yend
                        xstart = xstart + Dxslip
!!$                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
!!$$$$                  if (mod(islp,25) .eq. 0) then
!!$$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
!!$$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp
!!$$$$     $                    yendslp(islp)
!!$$$$                  endif
                     ENDDO
                     PRINT * , 'End of first slip =' , islp
                     INSlip(i) = islp
                  ENDIF
                  IF ( i==2 ) THEN
                     xstart = Xmin + Dxslip
!!$     ystart is hard coded to be 0.0 for crack
                     ystart = -0.01D-3
!!$               xend =xstart+(ymax-ystart)/tan(slip_angle(i))
!!$               print *, xstart, xend
!!$islp = 0
                     j = 0
                     DO WHILE ( xstart<Xmax )
                        j = j + 1
                        xend = xstart + (Ymin-ystart)/TAN(Slip_angle(i))
                        IF ( xend>Xmax ) THEN
                           xend = Xmax
                           yend = ystart - (xstart-xend)*TAN(Slip_angle(i))
                        ELSE
                           yend = Ymin
                        ENDIF
                        islp = islp + 1
                        LOCphi(islp) = 2
                        XSLp(islp) = xstart
                        XENdslp(islp) = xend
                        YSLp(islp) = ystart
                        YENdslp(islp) = yend
                        xstart = xstart + Dxslip
!!$                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
!!$$$$                  if (mod(islp,25) .eq. 0) then
!!$$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
!!$$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp
!!$$$$     $                    yendslp(islp)
!!$$$$                  endif

!!$print *, islp, xslp(islp), xendslp(islp)
                     ENDDO
                     PRINT * , 'End of second slip =' , islp
                     INSlip(i) = islp
                  ENDIF
               ELSE
                  xstart = Xmin
                  xend = Xmax
                  ystart = Dyslip
!!$$$$            nqslp = islp
                  j = 0
                  DO WHILE ( ystart<Ymax )
                     islp = islp + 1
                     j = j + 1
                     LOCphi(islp) = 3
                     XSLp(islp) = xstart
                     XENdslp(islp) = xend
                     YSLp(islp) = -ystart
                     YENdslp(islp) = -ystart
                     ystart = ystart + Dyslip
!!$$$$                  if (mod(islp,25) .eq. 0) then
!!$$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
!!$$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp
!!$$$$     $                    yendslp(islp)
!!$$$$                  endif


                  ENDDO

                  PRINT * , 'End of third slip =' , islp , j
                  INSlip(i) = islp
               ENDIF
            ENDDO



            NSLp = islp
            N_Active_slip = NSLp/SLIPSPC
            PRINT * , 'No. of active slip systems = ' , N_Active_slip , NSLp
!!$$$$     Calculate shifting parameters for slip planes
!!$$$$      Shifting is used to control dislocations leaving upper and
!!$$$$     $     entering lowere
            NSHift_slip = 0
            DO islp = 1 , NSLp
               xi = XSLp(islp)
               li = LOCphi(islp)
               IF ( li<3 ) THEN
                  DO jslp = 1 , NSLp
                     xj = XSLp(jslp)
                     lj = LOCphi(jslp)
                     IF ( lj<3 ) THEN
                        IF ( islp/=jslp ) THEN
                           IF ( li==lj ) THEN
                              IF ( xi==xj ) NSHift_slip(islp) = jslp
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
!!$$$$         print *, 'nshift', islp, nshift_slip(islp)
            ENDDO
!!$stop
      END SUBROUTINE GEN_SLIP_PLANES


      SUBROUTINE REMOVE_DISL_GLOBAL(Idisl_slip)
            USE MOD_DISL_PARAMETERS
            IMPLICIT NONE
!!$*--REMOVE_DISL_GLOBAL468
            INTEGER :: i , j , k , l , islp , Idisl_slip
            k = 0
            DO islp = 1 , NSLp
               IF ( NDIs(islp)>0 ) THEN
                  DO i = 1 , NDIs(islp)
                     k = k + 1
!!$                  call rmdisl_continuum(i,islp)
!!$                  call put_in_bucket(i,islp)
!!$ Dislocation is removed from the continuum on slip plane level and
!!$ global level..
                     IF ( k==Idisl_slip ) CALL RMDISL_G(i,islp)
                  ENDDO
               ENDIF
            ENDDO
            CALL ASSIGN_DISLOC_GLOBAL
      END SUBROUTINE REMOVE_DISL_GLOBAL

      SUBROUTINE PUT_IN_BUCKET(I,Islp)
            IMPLICIT NONE
!!$*--PUT_IN_BUCKET488
            INTEGER :: I , j , k , l , Islp , jslp
            DOUBLE PRECISION :: hold
            BOUt(1,Islp) = BOUt(1,Islp) + B_Dd(I,Islp)
            SDIs(I,Islp) = SDIs_out(1,Islp)
            RDIs(1,I,Islp) = XSLp(Islp) + SDIs_out(1,Islp)&
                 &                 *COSphi(LOCphi(Islp))
            RDIs(2,I,Islp) = YSLp(Islp) + SDIs_out(1,Islp)&
                 &                 *SINphi(LOCphi(Islp))
            ELEm_dis(I,Islp) = 0
            PRINT * , 'Dislocation' , I , 'on' , Islp , 'put in bucket'
      END SUBROUTINE PUT_IN_BUCKET


      SUBROUTINE CALC_BOUT
            IMPLICIT NONE
!!$*--CALC_BOUT504
            INTEGER :: i , j , k , l , islp , bplus , bminus
            DO islp = 1 , NSLp
               bplus = 0
               bminus = 0
               IF ( NDIs(islp)>0 ) THEN
                  DO i = 1 , NDIs(islp)
                     IF ( B_Dd(i,islp)>0 ) THEN
                        bplus = bplus + 1
                     ELSE
                        bminus = bminus + 1
                     ENDIF
                  ENDDO
                  BOUt(1,islp) = -(bplus-bminus)*BLEn
                  BOUt(2,islp) = 0
                  PRINT * , 'Bucket total =' , BOUt(1,islp)
               ENDIF
            ENDDO
      END SUBROUTINE CALC_BOUT


      SUBROUTINE RMDISL(I,Islp,Ar1,Ar2,Ar3,Ar4)
            IMPLICIT NONE
!!$*--RMDISL527
            INTEGER :: I , j , k , l , Islp , jslp
            DOUBLE PRECISION :: hold
            DOUBLE PRECISION :: Ar1(MXNDIS,NSLp) , Ar2(MXNDIS,NSLp) , &
                 &                    Ar3(MXNDIS,NSLp) , Ar4(MXNDIS,NSLp)
            hold = SDIs(I,Islp)
            DO j = I + 1 , NDIs(Islp)
               SDIs(j-1,Islp) = SDIs(j,Islp)
               B_Dd(j-1,Islp) = B_Dd(j,Islp)
               RDIs(1:3,j-1,Islp) = RDIs(1:3,j,Islp)
               ELEm_dis(j-1,Islp) = ELEm_dis(j,Islp)
               VPRev(j-1,Islp) = VPRev(j,Islp)
               IOBpin(j-1,Islp) = IOBpin(j,Islp)
               Ar1(j-1,Islp) = Ar1(j,Islp)
               Ar2(j-1,Islp) = Ar2(j,Islp)
               Ar3(j-1,Islp) = Ar3(j,Islp)
               Ar4(j-1,Islp) = Ar4(j,Islp)
            ENDDO
            PRINT * , 'Removing dislocation' , I , ' on ' , Islp
            NDIs(Islp) = NDIs(Islp) - 1
            TOT_disl = TOT_disl - 1
            CALL ASSIGN_DISLOC_GLOBAL
      END SUBROUTINE RMDISL

      SUBROUTINE RMDISL_G(I,Islp)
            IMPLICIT NONE
!!$*--RMDISL_G553
            INTEGER :: I , j , k , l , Islp , jslp
            DOUBLE PRECISION :: hold , roldq(3)
            DOUBLE PRECISION :: v(MXNDIS,NSLp) , c(MXNDIS,NSLp) , &
                 &                    tau(MXNDIS,NSLp) , sold(MXNDIS,NSLp) , bb
            PRINT * , 'Removing dislocation' , I , ' on ' , Islp
            roldq(1:3) = RDIs(1:3,I,Islp)
            hold = SDIs(I,Islp)
            DO j = I + 1 , NDIs(Islp)
               SDIs(j-1,Islp) = SDIs(j,Islp)
               B_Dd(j-1,Islp) = B_Dd(j,Islp)
               RDIs(1:3,j-1,Islp) = RDIs(1:3,j,Islp)
               ELEm_dis(j-1,Islp) = ELEm_dis(j,Islp)
               VPRev(j-1,Islp) = VPRev(j,Islp)
               IOBpin(j-1,Islp) = IOBpin(j,Islp)
            ENDDO
            NDIs(Islp) = NDIs(Islp) - 1
            TOT_disl = TOT_disl - 1
            CALL ASSIGN_DISLOC_GLOBAL

      END SUBROUTINE RMDISL_G

      SUBROUTINE ASSIGN_NEW_GLOBAL(Xtip_move)
            USE MOD_DISL_PARAMETERS
            IMPLICIT NONE
!!$*--ASSIGN_NEW_GLOBAL578
            INTEGER :: i , j , islp , ii , k , li , btot , FE_LOCATE , &
                 &           elem_old
            DOUBLE PRECISION :: pi , xend , yend , cphi , sphi
            DOUBLE PRECISION :: xstart , ystart , bb , bsign , ss , ss1 , sd
            DOUBLE PRECISION :: Xtip_move
            INTEGER :: ndisl1
            pi = 2.D0*ASIN(1.0D0)
            TOT_disl = 0

            DO islp = 1 , NSLp
               DO i = 1 , NNUc(islp)
                  RNUc(1,i,islp) = RNUc(1,i,islp) - Xtip_move
                  elem_old = ELEm_source(i,islp)
                  ELEm_source(i,islp) = FE_LOCATE(RNUc(:,i,islp),elem_old)
               ENDDO
            ENDDO

            DO islp = 1 , NSLp
               IF ( NDIs(islp)>0 ) TOT_disl = TOT_disl + NDIs(islp)
            ENDDO

            CALL CALC_BOUT
!!$ --- Set all quantities to zero
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cphi = COSphi(li)
               sphi = SINphi(li)

               DO i = 1 , NDIs(islp)
                  sd = SDIs(i,islp)
                  IF ( li>=3 ) THEN
                     RDIs(1,i,islp) = XSLp(islp) + sd*cphi
                     RDIs(2,i,islp) = YSLp(islp) + sd*sphi
                  ELSEIF ( YENdslp(islp)>0 ) THEN
                     RDIs(1,i,islp) = XSLp(islp) + sd*cphi
                     RDIs(2,i,islp) = YSLp(islp) + sd*sphi
                  ELSE
                     IF ( li==1 ) THEN
                        RDIs(1,i,islp) = XSLp(islp) - sd*cphi
                     ELSE
                        RDIs(1,i,islp) = XSLp(islp) + sd*cphi
                     ENDIF
                     RDIs(2,i,islp) = -(YSLp(islp)+sd*sphi)
                  ENDIF
                  ROLd(1,i,islp) = ROLd(1,i,islp) - Xtip_move
                  ROLd(2:3,i,islp) = ROLd(2:3,i,islp)
                  elem_old = ELEm_dis(i,islp)
                  ELEm_dis(i,islp) = FE_LOCATE(RDIs(:,i,islp),elem_old)
               ENDDO
            ENDDO
            CALL ASSIGN_DISLOC_GLOBAL

      END SUBROUTINE ASSIGN_NEW_GLOBAL

      SUBROUTINE ASSIGN_DISLOC_GLOBAL
            USE MOD_DISL_PARAMETERS
            IMPLICIT NONE
!!$*--ASSIGN_DISLOC_GLOBAL636
            INTEGER :: i , j , islp , ii , k , li , btot , fe_locate
            DOUBLE PRECISION :: pi , xend , yend , cphi , sphi
            DOUBLE PRECISION :: xstart , ystart , bb , bsign , ss , ss1
            INTEGER :: ndisl1
            pi = 2.D0*ASIN(1.0D0)
            TOT_disl = 0
            DO islp = 1 , NSLp
               IF ( NDIs(islp)>0 ) TOT_disl = TOT_disl + NDIs(islp)
            ENDDO
            CALL CALC_BOUT
!!$ --- Set all quantities to zero
            BURgers = 0.0D0
            BURg_length = BLEn
            THEta_e = 0.0D0
            THEta_s = 0.0D0
            R_Disl = 0.0D0
            R_Old = 0.0D0
            PK_force = 0.0D0
            PK_f = 0.0D0
            ELEm_disl = 0
            NDIsl = 0
            NDIsl_dd = 0
            ndisl1 = 0
            k = 0
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               sphi = SINphi(li)
               cphi = COSphi(li)
               xend = XENdslp(islp)
               yend = YENdslp(islp)
               xstart = XSLp(islp)
               ystart = YSLp(islp)
               ss = SDIs_out1(1,islp) + DISl_range_tol
!!$$$$         if (sdis_out(1,islp) .ne. sdis_out(1,islp)) then
!!$$$$            ss1 = sdis_out1(1,islp) -6.d0*disl_range_tol
!!$$$$         else
               ss1 = SDIs_out(1,islp)
!!$$$$         end if
               IF ( ABS(BOUt(1,islp))>0.0D0 ) THEN
                  bb = BLEn
                  btot = NINT(ABS(BOUt(1,islp))/BLEn)
                  IF ( BOUt(1,islp)<0 ) bb = -BLEn
!!$            print *, 'Bucket dislocations'
                  DO i = 1 , btot
                     k = k + 1
                     BURg_length(k) = bb
                     BURgers(1,k) = bb*cphi
                     BURgers(2,k) = bb*sphi
                     BURgers(3,k) = 0.0D0
                     IF ( bb>0 ) THEN
                        THEta_e(k) = -pi
                     ELSE
                        THEta_e(k) = 0.0D0
                     ENDIF
!!$$$$               if (bb > 0) then
!!$$$$                  if (yend > 0.0d0) then
!!$$$$                     theta_e(k) = pi
!!$$$$                  else
!!$$$$                     theta_e(k) = -pi
!!$$$$                  endif
!!$$$$               else
!!$$$$                  if (yend > 0.0d0) then
!!$$$$                     theta_e(k) = 0.0d0
!!$$$$                  else
!!$$$$                     theta_e(k) = 0.0d0
!!$$$$                  endif
!!$$$$               endif
!!$     Try making the position of the bucket dislocation to be at the
!!$     edge of the atomistic region
!!$     This means replacing the next 3 lines here with sdis_out1
                     R_Disl(1,k) = xstart + (ss1)*cphi
                     R_Disl(2,k) = ystart + (ss1)*sphi
                     R_Disl(3,k) = 0.0D0


                     R_Old(1:3,k) = R_Disl(1:3,k)
                     ELEm_disl(k) = 0
                     DISl_index(k) = 2
                     DISl_range(1,k) = ATOm_xmax + DISl_range_tol
                     DISl_range(2,k) = ATOm_ymax + DISl_range_tol
                     PRINT * , 'Buck' , i , islp , R_Disl(1,k) , R_Disl(2,k)
                  ENDDO
               ENDIF

               IF ( NDIs(islp)>0 ) THEN
                  DO i = 1 , NDIs(islp)
                     k = k + 1
                     ndisl1 = ndisl1 + 1
                     NDIsl_dd(k) = ndisl1
                     BURg_length(k) = B_Dd(i,islp)
                     BURgers(1,k) = B_Dd(i,islp)*COSphi(LOCphi(islp))
                     BURgers(2,k) = B_Dd(i,islp)*SINphi(LOCphi(islp))
                     BURgers(3,k) = 0.0D0
                     bb = B_Dd(i,islp)
!!$print *, 'qqq',k, ndisl1, ndisl_dd(k)
!!$     Do not switch cut plane directions in switching from bottom to top
!!$     This is only to ensure that the currect disl_u subroutine is
!!$     satisfied
!!$     Also since there is always going to be a bucket dislocation
!!$     The direction of the cut plane is not relevant (???)
                     IF ( bb>0 ) THEN
                        THEta_e(k) = -pi
                     ELSE
                        THEta_e(k) = 0.0D0
                     ENDIF
!!$$$$               if (bb > 0) then
!!$$$$                  if (yend > 0.0d0) then
!!$$$$                     theta_e(k) = -pi
!!$$$$                  else
!!$$$$                     theta_e(k) = 0.0d0
!!$$$$                  endif
!!$$$$               else
!!$$$$                  if (yend > 0.0d0) then
!!$$$$                     theta_e(k) = 0.0d0
!!$$$$                  else
!!$$$$                     theta_e(k) = -pi
!!$$$$                  endif
!!$$$$               endif

                     R_Disl(1:3,k) = RDIs(1:3,i,islp)
                     R_Old(1:3,k) = ROLd(1:3,i,islp)
                     ELEm_disl(k) = ELEm_dis(i,islp)
                     DISl_index(k) = 2
                     DISl_range(1,k) = xstart + ABS(ss*cphi)
                     DISl_range(2,k) = ABS(ystart) + ABS(ss*sphi)
                  ENDDO
               ENDIF
            ENDDO
            NDIsl = k
!!$      ndisl = k-1
            PRINT * , 'Global Disloc total = ' , NDIsl
            PRINT * , 'Gloabal Dislocations are'
            DO i = 1 , NDIsl
               WRITE (*,FMT='(I7,6(1X,E15.8),1X,2I7,3(1X,E15.8))') i , &
                    &          R_Old(1:2,i) , R_Disl(1:2,i) , DISl_range(1:2,i) , &
                    &          ELEm_disl(i) , NDIsl_dd(i) , BURgers(1:2,i) , THEta_e(i)
            ENDDO
      END SUBROUTINE ASSIGN_DISLOC_GLOBAL


      SUBROUTINE MOVE_DISLOC(Rhs)
            USE MOD_DISL_PARAMETERS
            IMPLICIT NONE
!!$*--MOVE_DISLOC780
            INTEGER :: i , j , k , l , ii , islp , li , elem_old , FE_LOCATE
            DOUBLE PRECISION :: si , bi , sini , cosi , tauij , taug , dsi , &
                 &                    sout(3) , fg , fgg , snew , roldq(3) , phi2 , &
                 &                    Rhs(*)
            DOUBLE PRECISION :: sinit , sfinal , tauii , tau , xi , yi
            DOUBLE PRECISION :: dist , eps , dlim , flim , dstar , s1 , &
                 &                    xstart , ystart
            DOUBLE PRECISION , ALLOCATABLE :: tau1(:,:) , v1(:,:) , c(:,:) , &
                 &                                  sold(:,:)
!!$$$$      double precision :: tau1(mxndis,nslp), v1(mxndis,nslp),
!!$$$$     $     c(mxndis,nslp), sold(mxndis,nslp)
            DOUBLE PRECISION :: cnew , ddt , dt , dx , dv , ddtm , vi , vpp , &
                 &                    vj , xqone
            DOUBLE PRECISION :: TLCE , rpile , ds , temp_s(MXNDIS,NSLp)
            LOGICAL :: move_bottom

!!$ Include obstacle parameters here
            TYPE (DD) ::dd1(TOT_disl)
            dt = TINcr
            ddt = dt
            eps = 1.0D-3*BLEn
            flim = eps
            dlim = flim/5.0D0
            dstar = 2.D0*BLEn
            rpile = dstar

            ALLOCATE (tau1(MXNDIS,NSLp),v1(MXNDIS,NSLp),c(MXNDIS,NSLp))
            ALLOCATE (sold(MXNDIS,NSLp))
            tau1 = 0.0D0
            v1 = 0.0D0
            c = 0.0D0
            sold = 0.0D0
            IF ( TOT_disl>0 ) THEN
               CALL ASSIGN_DISLOC_ONLY(dd1)
               CALL CALCFORCE(dd1,1,TOT_disl)
            ENDIF
            k = 1
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cosi = COSphi(li)
               sini = SINphi(li)
               sinit = SDIs_out1(1,islp)
               sfinal = SDIs_out(2,islp)
               IF ( NDIs(islp)>0 ) THEN
                  DO i = 1 , NDIs(islp)
                     phi2 = 2.D0*PHIslp(LOCphi(islp))
                     bi = B_Dd(i,islp)
                     si = SDIs(i,islp)
                     sold(i,islp) = si
                     elem_old = ELEm_dis(i,islp)
                     roldq(1:3) = RDIs(1:3,i,islp)
                     fg = REAL(dd1(k)%FORCED)
                     fgg = REAL(dd1(k)%FORCEDG)
                     IF ( ELEm_dis(i,islp)>0 )&
                          &              CALL FE_STRESS(ELEm_dis(i,islp),Rhs,sout)
                     phi2 = 2.D0*PHIslp(LOCphi(islp))
                     taug = 0.5D0*(sout(2)-sout(1))*SIN(phi2) + sout(3)&
                          &                *COS(phi2)
                     tauij = fg/bi
                     tauii = tauij + taug
                     vi = tauii*bi/BDRAG

                     vpp = vi/(1-fgg*dt/BDRAG)
                     IF ( vi/vpp>0 ) THEN
                        vi = vpp
                     ELSEIF ( ABS(vi)>ABS(vpp) ) THEN
                        vi = vi + vpp
                     ENDIF
                     IF ( ABS(vi)>VCUTOFF ) vi = SIGN(VCUTOFF,vi)
                     IF ( YENdslp(islp)<0.0D0 ) vi = -vi
                     tau1(i,islp) = tauii
                     dsi = vi*dt
                     snew = si + dsi
                     IF ( IOBpin(i,islp)==0 )&
                          &               WRITE (*,FMT='(A5,1X,2I7,1X,5(1X,E15.8))')&
                          &              'Disl.' , i , islp , vi , si , dsi , snew , dt
!!$ Check against obstacles
                     IF ( IOBpin(i,islp)>0 ) THEN
                        snew = si
                     ELSE
                        DO j = 1 , NOBs(islp)
                           dist = SOBs(j,islp) - si
                           IF ( ABS(dist)>=eps ) THEN
                              IF ( dsi/dist>1.0D0 ) THEN
                                 snew = SOBs(j,islp)
                                 IOBpin(i,islp) = -j
                                 PRINT * , 'Disloc ' , i , ' on ' , islp , &
                                      &                           'about to be pinned at obstacle' , j
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDIF
                     v1(i,islp) = (snew-si)/dt
                     k = k + 1
                  ENDDO
               ENDIF
            ENDDO
!!$     Determine underrelaxation parameters --- not really necessary but
!!$     still there
!!$      call undrex(v1,c)
            PRINT * , 'Underrelaxation factors'
            c = 1.0D0
!!$v1 = vi
            DO islp = 1 , NSLp
               IF ( NDIs(islp)>0 ) THEN
                  DO
                     DO i = 2 , NDIs(islp)
                        PRINT * , 'undrex' , i , islp , c(i,islp) , v1(i,islp)
                        dx = SDIs(i,islp) - SDIs(i-1,islp)
                        IF ( dx<0.0D0 ) THEN
                           WRITE (*,*) 'Wrong order of dislocations' , dx , &
                                &                           i , islp
                           STOP
                        ENDIF
                        vi = c(i,islp)*v1(i,islp)
                        vj = c(i-1,islp)*v1(i-1,islp)
                        dv = vi - vj
                        IF ( dv/=0.0D0 ) THEN
                           ddtm = -dx/dv
                           IF ( ddtm>0.0D0 ) THEN
                              IF ( B_Dd(i,islp)*B_Dd(i-1,islp)>0.0D0 ) THEN
                                 xqone = 1.0D0 - 1.0D-10
                                 IF ( dx/rpile<xqone ) PRINT * , &
                                      &                          'too short in pileup' , dx/rpile , i , &
                                      &                          islp
                                 dx = dx - rpile
                                 IF ( dx<0.0D0 .AND. TLCE(ABS(dx)/eps)<1.0D0 )&
                                      &                          dx = 0.0D0
                                 ddtm = -dx/dv
                              ENDIF
                              ddtm = ABS(ddtm)
                              IF ( TLCE(ddtm/ddt)<=1.0D0 ) THEN
                                 ddt = ddt/2.0D0
                                 cnew = 0.8D0*ddtm/dt
                                 c(i,islp) = c(i,islp)*cnew
                                 c(i-1,islp) = c(i-1,islp)*cnew
                                 GOTO 20
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                     EXIT
20                ENDDO
               ENDIF
            ENDDO

            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cosi = COSphi(li)
               sini = SINphi(li)
               IF ( SDIs_out1(1,islp)/=SDIs_out(1,islp) ) THEN
                  sinit = SDIs_out1(1,islp)
               ELSE
                  sinit = SDIs_out(1,islp)
               ENDIF
               sfinal = SDIs_out(2,islp)
               IF ( NDIs(islp)>0 ) THEN
                  PRINT * , 'Sinit = ' , sinit , ' Sfinal = ' , sfinal
                  DO i = 1 , NDIs(islp)
                     s1 = SDIs(i,islp) + c(i,islp)*v1(i,islp)*ddt
                     bi = B_Dd(i,islp)
                     dsi = c(i,islp)*v1(i,islp)*ddt
                     SDIs(i,islp) = SDIs(i,islp) + dsi
                     PRINT * , 'Disloc' , i , 'Old = ' , SDIs(i,islp) , &
                          &               ' new = ' , s1
                     RDIs(3,i,islp) = 0.0D0
                     snew = SDIs(i,islp)
                     temp_s(i,islp) = snew
!!$     Pin disloc at the  interface
!!$               if (elem_dis(i,islp) > 0) then
                     IF ( snew<sinit .OR. snew>sfinal ) THEN
                        IF ( snew<sinit ) THEN
                           PRINT * , &
                                &'Dislocation close to interface                                   &
                                &or atomisitcs' , snew
                           snew = sinit
                        ELSE
                           snew = sfinal
                        ENDIF
                     ENDIF
!!$               endif

                     SDIs(i,islp) = snew
                     ROLd(1:3,i,islp) = RDIs(1:3,i,islp)
                     IF ( li>=3 ) THEN
                        RDIs(1,i,islp) = XSLp(islp) + snew*cosi
                        RDIs(2,i,islp) = YSLp(islp) + snew*sini
                     ELSEIF ( YENdslp(islp)>0 ) THEN
                        RDIs(1,i,islp) = XSLp(islp) + snew*cosi
                        RDIs(2,i,islp) = YSLp(islp) + snew*sini
                     ELSE
                        IF ( li==1 ) THEN
                           RDIs(1,i,islp) = XSLp(islp) - snew*cosi
                        ELSE
                           RDIs(1,i,islp) = XSLp(islp) + ABS(snew*cosi)
                        ENDIF
                        RDIs(2,i,islp) = -(YSLp(islp)+snew*sini)
                     ENDIF

                     ELEm_dis(i,islp) = FE_LOCATE(RDIs(:,i,islp),elem_old)
                     PRINT * , 'Dislocation' , i , 'on' , islp , 'with' , bi ,&
                          &               'moved from' , ROLd(1,i,islp) , ROLd(2,i,islp) , &
                          &               'to' , RDIs(1,i,islp) , RDIs(2,i,islp) , si , dsi ,&
                          &               ELEm_dis(i,islp)


!!$$$$               print *, 'pk-force =', i, islp, tauii*bi,
!!$$$$     $              real(dd1(k)%forced)
!!$$$$  print *, 'Dislocation ', i, ' on ', islp, 'moved from',
!!$$$$  $                elem_old, ' to ', elem_dis(i,islp)
!!$$$$               rold(1:3,i,islp) = rdis(i:3,i,islp)
                     k = k + 1
                  ENDDO
               ENDIF
            ENDDO


            DO islp = 1 , NSLp
               IF ( NDIs(islp)>1 ) THEN
                  i = 2
                  DO WHILE ( i<=NDIs(islp) )
                     IF ( B_Dd(i,islp)*B_Dd(i-1,islp)<0.0D0 ) THEN
                        IF ( TLCE(SDIs(i-1,islp)/SDIs(i,islp))>1.0D0 ) THEN
                           si = SDIs(i,islp)
                           WRITE (*,*) 'Wrong order' , i , islp , SDIs(i,islp)&
                                &                           , SDIs(i-1,islp)
                        ELSE
                           ds = SDIs(i,islp) - SDIs(i-1,islp)
                           IF ( ds<=XLE ) THEN
                              CALL RMDISL(i,islp,v1,tau1,c,sold)
                              CALL RMDISL(i-1,islp,v1,tau1,c,sold)
                              i = i - 1
                           ENDIF
                        ENDIF
                     ENDIF
                     i = i + 1
                  ENDDO
               ENDIF
            ENDDO

            DO islp = 1 , NSLp
               DO i = 1 , NDIs(islp)
                  DO j = 1 , NOBs(islp)
                     IF ( IOBpin(i,islp)+j==0 ) THEN
                        dist = DABS(SOBs(j,islp)-SDIs(i,islp))
                        IF ( dist<eps ) THEN
                           IOBpin(i,islp) = j
                           WRITE (*,*) 'Dislocation ' , i , ' on ' , islp , &
                                &                           ' pinned at obstacle ' , j
                        ELSE
                           IOBpin(i,islp) = 0
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO


            CALL RLSDIS(tau1)

!!$     --- Pass dislocations from top to bottom or simply put in bucket
!!$     if xstart < atom_xmin + tol
            move_bottom = .TRUE.
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cosi = COSphi(li)
               sini = SINphi(li)
               xstart = XSLp(islp)
               ystart = YSLp(islp)
               IF ( SDIs_out1(1,islp)/=SDIs_out(1,islp) ) THEN
                  move_bottom = .FALSE.
                  sinit = SDIs_out(1,islp)
               ELSE
                  sinit = SDIs_out(1,islp)
                  move_bottom = .TRUE.
               ENDIF
               sfinal = SDIs_out(2,islp)
               i = 1
               k = 1
               DO WHILE ( i<NDIs(islp) )
                  IF ( NDIs(islp)>0 ) THEN
                     IF ( temp_s(i,islp)>=sinit ) THEN
                        i = i + 1
                     ELSEIF ( move_bottom ) THEN
                        l = islp + NSHift_slip(islp)
                        PRINT * , 'Removing from top and putting in bottom' , &
                             &                  temp_s(i,islp) , sinit
                        IF ( XSLp(islp)<ATOm_xmin ) THEN
                           CALL RMDISL_G(i,islp)
                        ELSE
                           CALL CRDISL(temp_s(i,islp),sinit,l,B_Dd(i,islp),1)
                           CALL RMDISL_G(i,islp)
                        ENDIF
                     ELSE
                        i = i + 1
                     ENDIF
                     k = k + 1
                  ENDIF
               ENDDO
            ENDDO


            DO islp = 1 , NSLp
               IF ( NDIs(islp)>0 ) THEN
                  DO i = 1 , NDIs(islp)
                     VPRev(i,islp) = v1(i,islp)
                  ENDDO
               ENDIF
            ENDDO

            CALL ASSIGN_DISLOC_GLOBAL
!!$$$$      deallocate(tau1,v1,c)
!!$$$$      deallocate(sold)
      END SUBROUTINE MOVE_DISLOC




      SUBROUTINE NUCLEATE_DISL(Rhs)
            IMPLICIT NONE
!!$*--NUCLEATE_DISL1101
            INTEGER islp , jslp , i , j , k , l
            DOUBLE PRECISION s , xi , xlh , bsign , sout(3) , Rhs(*) , phi2 , &
                 &                 tau
            DOUBLE PRECISION ev_convert , fact , tau_source , tauij_s , taug_s
            DOUBLE PRECISION qsign
            TYPE (DD) ::dd1(TOT_size)
            LOGICAL create
            bsign = 1.0D0
            create = .FALSE.
            ev_convert = 1.602176462
            fact = 1.D0/ev_convert/1.D-5
            CALL ASSIGN_DISLOC_ONLY(dd1)
            CALL ASSIGN_SOURCE_ONLY(dd1)
            CALL CALCFORCE(dd1,TOT_disl+1,TOT_size)
            k = TOT_disl + 1
            DO islp = 1 , NSLp
               IF ( NNUc(islp)>0 ) THEN
                  DO i = 1 , NNUc(islp)
                     IF ( ELEm_source(i,islp)>0 )&
                          &              CALL FE_STRESS(ELEm_source(i,islp),Rhs,sout)
                     phi2 = 2.D0*PHIslp(LOCphi(islp))
                     taug_s = 0.5D0*(sout(2)-sout(1))*SIN(phi2) + sout(3)&
                          &                  *COS(phi2)
                     tauij_s = REAL(dd1(k)%FORCED)/BLEn
                     tau_source = taug_s + tauij_s
!!$$$$                 taui(i,islp) = taui(i,islp) + real(dd1(k)%forced)/b
!!$$$$     $                + tau

                     sout = sout*fact
                     WRITE (*,FMT='(3I7,1X,3(1X,E15.7))') i , islp , k , &
                          &                tau_source , REAL(dd1(k)%FORCED)/BLEn*fact , &
                          &                TAUi(i,islp)*fact


                     IF ( TNLaps(i,islp)<0.5*TINcr ) TNLaps(i,islp) = 0.0D0
                     qsign = 1.0D0

                     IF ( tau_source*TAUi(i,islp)<0.0D0 ) qsign = -1.0D0

                     IF ( ABS(tau_source)>=T_Fr(i,islp) ) THEN
                        TNLaps(i,islp) = TNLaps(i,islp) + qsign*TINcr
                        PRINT * , 'Source ' , i , ' on ' , islp , 'tnlaps' , &
                             &                  TNLaps(i,islp)
                     ELSEIF ( TNLaps(i,islp)>0 ) THEN
                        TNLaps(i,islp) = TNLaps(i,islp) - TINcr
                     ENDIF
                     IF ( TNLaps(i,islp)<0.5*TINcr ) TNLaps(i,islp) = 0.0D0
                     k = k + 1
                     TAUi(i,islp) = tau_source
                  ENDDO
               ENDIF
            ENDDO
            DO islp = 1 , NSLp
               IF ( NNUc(islp)>0 ) THEN
                  DO i = 1 , NNUc(islp)
                     IF ( TNLaps(i,islp)>TNUC ) THEN
                        xlh = XLNuc(i,islp)/2.D0
                        xi = SNUc(i,islp)
                        IF ( TAUi(i,islp)<0.0D0 ) bsign = -1.0D0
                        IF ( YENdslp(islp)<0.0D0 ) bsign = -bsign
!!$print *, 'Generating dislocations on slip plane',islp
                        CALL CRDISL(xi-xlh,xi,islp,-BLEn*bsign,i)
                        CALL CRDISL(xi+xlh,xi,islp,BLEn*bsign,i)
                        create = .TRUE.
                        TNLaps(i,islp) = 0.0D0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
!!$if (create) stop
            CALL ASSIGN_DISLOC_GLOBAL
      END SUBROUTINE NUCLEATE_DISL

      SUBROUTINE NUCLEATE_ATOMISTIC_DISL(Islp,S_disl,Th_e,Th_s,Bv)
!!$     Nucleate a dislocation passed in from atomistics
!!$     Automatically handles dislocations to be put into the bucket

!!$     islp -- Slip plane to nucleate dislocation
!!$     s_disl --- Position on slip plane to nucleate dislocation
!!$     bv -- Burgers vector of dislocation
!!$     th_e, th_s --- Angles of dislocation to be inserted
            IMPLICIT NONE
!!$*--NUCLEATE_ATOMISTIC_DISL1184
            INTEGER :: Islp , i , j , k , l , li
            DOUBLE PRECISION :: S_disl , Th_e , Th_s , Bv(3)
            DOUBLE PRECISION :: bsign , bmag
            DOUBLE PRECISION , PARAMETER :: DTOL = 1.D-3
!!$     Angle is already checked when finding slip plane
!!$     Now we only re-check with different criteria

!!$     Assign correct sign for dislocation

            li = LOCphi(Islp)
            IF ( ABS(Bv(1)-BLEn*COSphi(li))<DTOL .AND. &
                 &     ABS(Bv(2)-BLEn*SINphi(li))<DTOL ) THEN
               bsign = 1.0D0
            ELSE
               bsign = -1.0D0
            ENDIF


            CALL CRDISL(S_disl,S_disl-6.D0*BLEn,Islp,bsign*BLEn)
            CALL ASSIGN_DISLOC_GLOBAL
      END SUBROUTINE NUCLEATE_ATOMISTIC_DISL

      SUBROUTINE CRDISL(Si,Slim,Islp,Bi,Nuc_num)
            IMPLICIT NONE
!!$*--CRDISL1209
            INTEGER i , j , k , l , m , n , Islp , FE_LOCATE , el_source
            INTEGER , OPTIONAL :: Nuc_num
            DOUBLE PRECISION Si , Slim , Bi , shold , xcl
            DOUBLE PRECISION xint , dstr , dstar , sold , stest
            dstar = 2*BLEn
            n = NDIs(Islp)
            shold = Si
            xcl = Si - Slim
            xint = 3.D0*2.D0*BLEn
            dstr = 2.D0*dstar
            RPIle = 2.D0*BLEn
            i = 0
            IF ( PRESENT(Nuc_num) ) el_source = ELEm_source(Nuc_num,Islp)
!!$ --- Make sure dislocations are in order
!!$TODO !!$ -- Make sure dislocations are atleast 6*b away from each other
            i = 1
            DO WHILE ( shold>SDIs(i,Islp) .AND. i<=NDIs(Islp) )
               i = i + 1
            ENDDO
            k = i
!!$print *, 'Insertion location', k
            NDIs(Islp) = NDIs(Islp) + 1
            TOT_disl = TOT_disl + 1
            TOT_size = TOT_size + 1
            DO i = NDIs(Islp) , k + 1 , -1
               SDIs(i,Islp) = SDIs(i-1,Islp)
               B_Dd(i,Islp) = B_Dd(i-1,Islp)
               RDIs(1:3,i,Islp) = RDIs(1:3,i-1,Islp)
               ELEm_dis(i,Islp) = ELEm_dis(i-1,Islp)
               ROLd(1:3,i,Islp) = ROLd(1:3,i-1,Islp)
               VPRev(i,Islp) = VPRev(i-1,Islp)
               IOBpin(i,Islp) = IOBpin(i-1,Islp)
            ENDDO

            SDIs(k,Islp) = shold
            B_Dd(k,Islp) = Bi
            IOBpin(k,Islp) = 0
            IF ( YENdslp(Islp)>0 ) THEN
               RDIs(1,k,Islp) = XSLp(Islp) + Si*COSphi(LOCphi(Islp))
               RDIs(2,k,Islp) = YSLp(Islp) + Si*SINphi(LOCphi(Islp))
            ELSE
               IF ( LOCphi(Islp)==1 ) THEN
                  RDIs(1,k,Islp) = XSLp(Islp) - Si*COSphi(LOCphi(Islp))
               ELSE
                  RDIs(1,k,Islp) = XSLp(Islp) + ABS(Si*COSphi(LOCphi(Islp)))
               ENDIF
               RDIs(2,k,Islp) = -(YSLp(Islp)+Si*SINphi(LOCphi(Islp)))
            ENDIF
            RDIs(3,k,Islp) = 0.0D0
            ELEm_dis(k,Islp) = FE_LOCATE(RDIs(1:3,k,Islp),1)
            WRITE (*,FMT='(A30,1X,I7,3(1X,E15.7),1X, A15,I7,1X,E15.8)')&
                 &        'Dislocation Generated on' , Islp , SDIs(k,Islp) , &
                 &       RDIs(1,k,Islp) , RDIs(2,k,Islp) , 'in element' , &
                 &       ELEm_dis(k,Islp) , B_Dd(k,Islp)
      END SUBROUTINE CRDISL

      SUBROUTINE GEN_SLIP_ENDS
            IMPLICIT NONE
!!$*--GEN_SLIP_ENDS1268
            INTEGER :: i , j , k , l , islp , li
            DOUBLE PRECISION :: xstart , xend , ystart , yend , xend1 , yend1
            DOUBLE PRECISION :: sphi , cphi , tphi
            DOUBLE PRECISION :: dlmax , dlmax1
            DOUBLE PRECISION :: xmax1 , xmin1 , ymin1 , ymax1

!!$     Initialize the dislocation end variables

            SDIs_out = 0.0D0
            SDIs_out1 = 0.0D0

!!$$$$      xmax1 = atom_xmax + sdis_out_tol
!!$$$$      xmin1 = (atom_xmin - sdis_out_tol)/2.0d0
!!$$$$      ymax1 = atom_ymax + sdis_out_tol
!!$$$$      ymin1 = atom_ymin - sdis_out_tol

            xmax1 = ATOm_xmax
            xmin1 = ATOm_xmin
            ymax1 = ATOm_ymax
            ymin1 = ATOm_ymin

            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cphi = COSphi(li)
               sphi = SINphi(li)
               tphi = sphi/cphi
               xstart = XSLp(islp)
               ystart = YSLp(islp)
               xend = XENdslp(islp)
               yend = YENdslp(islp)
               IF ( li/=3 ) THEN
                  dlmax = ABS(yend-ystart)/sphi
               ELSE
                  dlmax = xend - xstart
               ENDIF
               SDIs_out(1,islp) = 0.0D0 + 1.D-3
               SDIs_out(2,islp) = dlmax - 1.D-3
               IF ( xstart>xmin1 .AND. xstart<xmax1 ) THEN
                  IF ( li/=3 ) THEN
                     IF ( yend>0 ) THEN
                        xend1 = xstart + (ymax1-ystart)/tphi
                        yend1 = ATOm_ymax
                        IF ( xend1>xmax1 ) THEN
                           xend1 = xmax1
                           yend1 = ABS(xend1-xstart)*ABS(tphi)
                        ENDIF
                        IF ( xend1<xmin1 ) THEN
                           xend1 = xmin1
                           yend1 = ABS(xend1-xmin1)*ABS(tphi)
                        ENDIF
                     ENDIF
                     dlmax1 = SQRT((xstart-xend1)**2+(ystart-yend1)**2)&
                          &                  + SDIs_out_tol
                     SDIs_out1(1,islp) = dlmax1 - 2.0D0
                     WRITE (*,FMT='(2I7,5(1X,E15.8))')  , li , islp , &
                          &                XSLp(islp) , xend1 , yend1 , SDIs_out1(1,islp)
                  ENDIF
               ENDIF
            ENDDO
      END SUBROUTINE GEN_SLIP_ENDS

      SUBROUTINE GEN_SOURCES
            IMPLICIT NONE
!!$*--GEN_SOURCES1332
            DOUBLE PRECISION :: XE , XNU , XLAmbda , XMU , pi
            DOUBLE PRECISION :: CC(6,6)
            INTEGER :: I_Elas
            COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
            INTEGER :: i , j , k , l , islp , iseed1 , iseed2 , fe_locate
            LOGICAL :: accept
            DOUBLE PRECISION :: ev_convert , e , mu , nu , lambda , lnuc
            DOUBLE PRECISION , ALLOCATABLE :: source_dist(:)
!!$ --- Obtain Properties of the material from main input
!!$---- Calculate total no. of sources based on density and area
!!$ --- Generate an array of source strengths
!!$ --- Assign each source the strength and compute L_nuc
!!$ --- Place the source according to L_nuc

            DOUBLE PRECISION :: dlmax , dlengthslp , xstart , xend , ystart , &
                 &                    yend
            DOUBLE PRECISION :: qsrc , rand(2) , rand1(2) , sphi , cphi
            INTEGER :: nsrc , ii , li
            DOUBLE PRECISION :: atomistic_exclusion , lnuc_exclusion , &
                 &                    factor , tol



!!$     currently hard coded -- need to use in input file
            pi = 3.1415926535898
            ev_convert = 1.602176462
            SOUrce_den = 66 !!$ per /um^2
            SOUrce_den = 66.0D-8  !!$ per Angstrom^2
            TOT_source = INT(PROcess_area*SOUrce_den) + 10
            qsrc = TOT_source*1.D0/(N_Active_slip*1.D0) !!$ -- no. of sources pe
            PRINT * , 'QSRC =' , qsrc , TOT_source



!!$      tot_source = 1000
            AVG_source_str = 50 !!$ MPa
            SD_source_str = 5 !!$MPa
            AVG_source_str = AVG_source_str*ev_convert*1.D-5
            SD_source_str = SD_source_str*ev_convert*1.D-5
            PRINT * , 'Source strength =' , AVG_source_str , XMU , XNU
            PRINT * , 'Source Strength MPa' , &
                 &      AVG_source_str/ev_convert/1.D-5 , XMU/ev_convert/1.D-2 , &
                 &      XE/ev_convert/1.D-2
            factor = XMU*BLEn/(2.D0*pi*(1.D0-XNU))
            lnuc = factor/AVG_source_str
            tol = lnuc/4.D0
            PRINT * , 'Lnuc =' , lnuc , BLEn , TOT_source
            iseed1 = 1234
            iseed2 = 567

            CALL RMARIN(iseed1,iseed2)
            ALLOCATE (source_dist(1000))
            CALL RGAUSS(source_dist,1000,AVG_source_str,SD_source_str)
!!$$$$      do i = 1, tot_source
!!$$$$         write(*,fmt='(I7,1X,2E15.8)'), i, source_dist(i),
!!$$$$     $        avg_source_str
!!$$$$      end do
!!$ Generate sources
!!$ For now assume that qsrc always less than 1
!!$$$$      if (qsrc > 1.0d0) then
!!$$$$         nsrc = int(qsrc*10.d0)/10
!!$$$$         qsrc = qsrc-nsrc
!!$$$$      endif

!!$      atom_exclusion = (atom_xmax-atom_xmin)/tan(phislp(1))
            ALLOCATE (NNUc(NSLp))
            ALLOCATE (SNUc(MXNNUC,NSLp))
            ALLOCATE (RNUc(3,MXNNUC,NSLp))
            ALLOCATE (T_Fr(MXNNUC,NSLp))
            ALLOCATE (XLNuc(MXNNUC,NSLp))
            ALLOCATE (TNLaps(MXNNUC,NSLp))
            ALLOCATE (TAUi(MXNNUC,NSLp))
            ALLOCATE (ELEm_source(MXNNUC,NSLp))
!!$      allocate(nuc_pk_stress(3,1,nslp))

            ALLOCATE (NDIs(NSLp))
            ALLOCATE (SDIs(MXNDIS,NSLp))
            ALLOCATE (B_Dd(MXNDIS,NSLp))
            ALLOCATE (VPRev(MXNDIS,NSLp))
            ALLOCATE (IOBpin(MXNDIS,NSLp))
            ALLOCATE (RDIs(3,MXNDIS,NSLp))
            ALLOCATE (ROLd(3,MXNDIS,NSLp))
            ALLOCATE (ELEm_dis(MXNDIS,NSLp))
            ALLOCATE (BOUt(2,NSLp))
            ALLOCATE (SDIs_out(2,NSLp))
            ALLOCATE (SDIs_out1(2,NSLp))

            ALLOCATE (NOBs(NSLp))
            ALLOCATE (SOBs(MXNOBS,NSLp))
            ALLOCATE (TAU_obs(MXNOBS,NSLp))



            NDIs = 0
            SDIs = 0.0D0
            B_Dd = 0.0D0
            VPRev = 0.0D0
            IOBpin = 0
            RDIs = 0.0D0
            ROLd = 0.0D0
            ELEm_dis = 0

            NNUc = 0
            RNUc = 0.0D0
            SNUc = 0.0D0
            T_Fr = 0.0D0
            XLNuc = 0.0D0
            TNLaps = 0.0D0
            TAUi = 0.0D0
            BOUt = 0.0D0
            SDIs_out = 0.0D0
            SDIs_out1 = 0.0D0

            NOBs = 0
            SOBs = 0.0D0
            TAU_obs = 0.0D0
            ii = 0
            PRINT * , 'QSRC = ' , qsrc

            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cphi = COSphi(LOCphi(islp))
               sphi = SINphi(LOCphi(islp))
               xstart = XSLp(islp)
               ystart = YSLp(islp)
               xend = XENdslp(islp)
               yend = YENdslp(islp)
               IF ( li<3 ) THEN
                  dlmax = ABS(yend-ystart)/ABS(sphi)
               ELSE
                  dlmax = xend - xstart
               ENDIF
!!$$$$         write(*,fmt='(A20,1X,2I7,5(1X,E15.8))'), 'Slip plane',
!!$$$$     $        islp, li, xstart, xend, ystart, yend, dlmax
            ENDDO
            PRINT * , 'Source generation'

            islp = 2
            PRINT * , 'Total slip planes = ' , NSLp
!!$$$$      do while (islp <= nslp)
!!$!!$$      DO islp = 1 , NSLp
!!$!!$$         li = LOCphi(islp)
!!$!!$$         cphi = COSphi(li)
!!$!!$$         sphi = SINphi(li)
!!$!!$$         xstart = XSLp(islp)
!!$!!$$         ystart = YSLp(islp)
!!$!!$$         xend = XENdslp(islp)
!!$!!$$         yend = YENdslp(islp)
!!$!!$$         IF ( li<3 ) THEN
!!$!!$$            dlmax = ABS(yend-ystart)/ABS(sphi)
!!$!!$$         ELSE
!!$!!$$            dlmax = xend - xstart
!!$!!$$         ENDIF
!!$!!$$         CALL RANMAR(rand,2)
!!$!!$$ 
!!$!!$$  !!$        Hard code one source on +60 slip plane a little distacnce   a
!!$!!$$         if (locphi(islp) .ne. 3) then
!!$!!$$            if (xslp(islp) > 5.d0*blen-1.d0 .and.  xslp(islp) < 5.d0*blen+1.d0) then
!!$!!$$                  if (yend > 0.0d0) then
!!$!!$$                     if (yendslp(islp) > 0.0d0) then
!!$!!$$                        if (locphi(islp) == 1) then
!!$!!$$                           dlengthslp = lnuc + 100.0d0 + rand(1)*lnuc
!!$!!$$                        else
!!$!!$$                           dlengthslp = 2.d0*lnuc + 100.0d0 + rand(1)*lnuc
!!$!!$$                        end if
!!$!!$$                     else
!!$!!$$                        dlengthslp = 4.d0*lnuc + 100.0d0 + rand(1)*lnuc
!!$!!$$                     end if
!!$!!$$                     if (locphi(islp) == 1) then
!!$!!$$                        if (yendslp(islp) < 0.0d0) then
!!$!!$$                           cycle
!!$!!$$                        end if
!!$!!$$                     end if
!!$!!$$                     k = 1
!!$!!$$                     nnuc(islp) = k
!!$!!$$                     t_FR(k,islp) = avg_source_str
!!$!!$$                     xLnuc(k,islp) = lnuc
!!$!!$$                     ii = ii + 1
!!$!!$$                     snuc(k,islp) = dlengthslp
!!$!!$$                     if (yend > 0) then
!!$!!$$                        rnuc(1,k,islp) = xstart + dlengthslp*cphi
!!$!!$$                        rnuc(2,k,islp) = ystart + dlengthslp*sphi
!!$!!$$                        rnuc(3,k,islp) = 0.0
!!$!!$$                     else
!!$!!$$                        if (locphi(islp) ==1) then
!!$!!$$                           rnuc(1,k,islp) = xstart - dlengthslp*cphi
!!$!!$$                        else
!!$!!$$                           rnuc(1,k,islp) = xstart+abs(dlengthslp*cphi)
!!$!!$$                        endif
!!$!!$$                        rnuc(2,k,islp) = -(ystart + dlengthslp*sphi)
!!$!!$$                        rnuc(3,k,islp) = 0.0
!!$!!$$                     endif
!!$!!$$                     elem_source(k,islp) = fe_locate(rnuc(:,k,islp),1)                     
!!$!!$$                     write( *,fmt='(A10,1X,3(1x,E15.7))') 'Coords are', snuc(k,islp),rnuc(1,k,islp), rnuc(2,k,islp)
!!$!!$$                     print *, 'Located in Element ', elem_source(k,islp)
!!$!!$$                     print *, 'Islp = ', islp, locphi(islp), cphi, sphi
!!$exit
!!$!!$$                  endif
!!$!!$$            else
!!$!!$$               cycle
!!$!!$$            endif
!!$!!$$         else
!!$!!$$            cycle
!!$!!$$         endif

!!$$$$!!$     Generate source randomly from gaussian distribution
!!$$$$
!!$$$$         if (rand(1) <= qsrc) then
!!$$$$            accept = .false.
!!$$$$            k = 1
!!$$$$            ii = ii + 1
!!$$$$            if (yend > 0.0d0) then
!!$$$$               t_FR(k,islp) = source_dist(ii)
!!$$$$            else
!!$$$$               t_FR(k,islp) = source_dist(ii)*10
!!$$$$            endif
!!$$$$            xLnuc(k,islp) = factor/t_FR(k,islp)
!!$$$$            if (dlmax < 2.d0*xLnuc(k,islp)) then
!!$$$$c$$$  print *, 'Rejecting slip plane', ii, islp, li,
!!$$$$c$$$  $              dlmax,(2.d0*xLnuc(k,islp)), t_FR(k,islp)
!!$$$$               t_FR(k,islp) = 0.0d0
!!$$$$               xLnuc(k,islp) = 0.0d0
!!$$$$               ii = ii -1
!!$$$$               islp = islp + 1
!!$$$$               goto 10
!!$$$$            endif
!!$$$$            nnuc(islp) = k
!!$$$$
!!$$$$            do while (accept == .false.)
!!$$$$               call ranmar(rand1,2)
!!$$$$               dlengthslp = rand1(1)*(dlmax)
!!$$$$               if (locphi(islp) == 3) then
!!$$$$                  if (abs(ystart) > atom_ymax + xLnuc(k,islp)*2.0d0)
!!$$$$     $                 then
!!$$$$                     if (xstart + dlengthslp < atom_xmin .or.
!!$$$$     $                    xstart + dlengthslp > atom_xmax) then
!!$$$$                        accept = .true.
!!$$$$                     endif
!!$$$$                  endif
!!$$$$               endif
!!$$$$               if (dlengthslp > xLnuc(k,islp)/2.d0 + tol/2.d0) then
!!$$$$                  if (dlengthslp < dlmax - xLnuc(k,islp)/2.d0-tol/2.
!!$$$$     $                 then
!!$$$$                     accept = .true.
!!$$$$                  endif
!!$$$$               endif
!!$$$$            end do
!!$$$$            snuc(k,islp) = dlengthslp
!!$$$$            print *, 'Source no. generated = ', ii, ' on ', islp,
!!$$$$     $           li, t_FR(k,islp)/ev_convert/1.d-5, snuc(k,islp)
!!$$$$
!!$$$$
!!$$$$            if ( locphi(islp) .ne. 3) then
!!$$$$               if (yend > 0) then
!!$$$$                  rnuc(1,k,islp) = xstart + dlengthslp*cphi
!!$$$$                  rnuc(2,k,islp) = ystart + dlengthslp*sphi
!!$$$$                  rnuc(3,k,islp) = 0.0
!!$$$$               else
!!$$$$                  if (locphi(islp) == 1) then
!!$$$$                     rnuc(1,k,islp) = xstart - dlengthslp*cphi
!!$$$$                  else
!!$$$$                     rnuc(1,k,islp) = xstart + abs(dlengthslp*cphi)
!!$$$$                  endif
!!$$$$                  rnuc(2,k,islp) = -(ystart + dlengthslp*sphi)
!!$$$$                  rnuc(3,k,islp) = 0.0
!!$$$$               endif
!!$$$$            else
!!$$$$               rnuc(1,k,islp) = xstart + dlengthslp*cphi
!!$$$$               rnuc(2,k,islp) = ystart
!!$$$$            endif
!!$$$$
!!$$$$            elem_source(k,islp) = fe_locate(rnuc(:,k,islp),1)
!!$$$$            write( *,fmt='(A10,1X,3(1x,E15.7))') 'Coords are',
!!$$$$     $           snuc(k,islp),rnuc(1,k,islp), rnuc(2,k,islp)
!!$$$$            print *, 'Located in Element ', elem_source(k,islp)
!!$$$$            islp = islp + 100
!!$$$$         else
!!$$$$            islp = islp + 100
!!$$$$         endif
!!$!!$$      ENDDO

            DO islp = 1 , NSLp
               IF ( NNUc(islp)>0 ) THEN
                  DO i = 1 , NNUc(islp)
                     PRINT * , 'Source = ' , i , islp , ELEm_source(i,islp)
                  ENDDO
               ENDIF
            ENDDO
            PRINT * , 'Total no. of sources = ' , ii

            NNTot = ii
            TOT_source = ii
            TOT_size = TOT_source

            CALL GEN_SLIP_ENDS

!!$$$$      stop
      END SUBROUTINE GEN_SOURCES

      SUBROUTINE PLOT_SLIP
            IMPLICIT NONE
!!$*--PLOT_SLIP1638
            INTEGER :: i , j , k , l , li , islp , ifile
            INTEGER :: kk(3)
            DOUBLE PRECISION :: x , y , xstart , xend , ystart , yend , sn , &
                 &                    x1 , y1
            DOUBLE PRECISION :: cphi , sphi , tphi
            CHARACTER*80 filename
            INTEGER :: tot_obs , tot_source , tot_pts

!!$$$$      filename='slip_lines.plt'
!!$$$$      ifile = 999
!!$$$$      open(unit=ifile,file=filename,status='unknown')
!!$$$$      write(999,*) 'Variables = X Y'
!!$$$$      do islp = 1,nslp
!!$$$$         if (nnuc(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            cphi = cosphi(li)
!!$$$$            sphi = sinphi(li)
!!$$$$            tphi = sphi/cphi
!!$$$$            xstart = xslp(islp)
!!$$$$            xend = xendslp(islp)
!!$$$$            ystart = yslp(islp)
!!$$$$            yend = yendslp(islp)
!!$$$$            write(999,*) "Zone T = line"
!!$$$$            write(999,fmt='(2(1X,E18.11))') xstart, ystart
!!$$$$            write(999,fmt='(2(1X,E18.11))') xend, yend
!!$$$$         endif
!!$$$$      end do
!!$$$$      close(ifile)
            filename = ''
            filename = 'source_obs.vtk'
            tot_obs = 0
            tot_source = 0
            tot_pts = 0
            DO islp = 1 , NSLp
               tot_source = tot_source + NNUc(islp)
               tot_obs = tot_obs + NOBs(islp)
               tot_pts = tot_pts + NNUc(islp) + NOBs(islp)
            ENDDO

            OPEN (UNIT=ifile,FILE=filename,STATUS='unknown')
            WRITE (ifile,FMT='(A)') '# vtk DataFile Version 2.0'
            WRITE (ifile,FMT='(A)') 'Strains from CADD'
            WRITE (ifile,FMT='(A)') 'ASCII'
            WRITE (ifile,FMT='(A)') 'DATASET UNSTRUCTURED_GRID'
            WRITE (ifile,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , tot_pts , 'float'
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               cphi = COSphi(li)
               sphi = SINphi(li)
               yend = YENdslp(islp)
!!$     Sources are written first
               DO i = 1 , NNUc(islp)
                  sn = SNUc(i,islp)
                  IF ( yend>0.0D0 ) THEN
                     x1 = XSLp(islp) + sn*cphi
                     y1 = YSLp(islp) + sn*sphi
                  ELSE
                     IF ( li==1 ) THEN
                        x1 = XSLp(islp) - sn*cphi
                     ELSE
                        x1 = XSLp(islp) + ABS(sn*cphi)
                     ENDIF
                     y1 = -(YSLp(islp)+sn*sphi)
                  ENDIF
                  WRITE (ifile,'(3(1X,E14.6))') x1 , y1 , 0.0D0
               ENDDO
!!$     Obstacles are written next
               DO i = 1 , NOBs(islp)
                  sn = SOBs(i,islp)
                  IF ( yend>0.0D0 ) THEN
                     x1 = XSLp(islp) + sn*cphi
                     y1 = YSLp(islp) + sn*sphi
                  ELSE
                     IF ( li==1 ) THEN
                        x1 = XSLp(islp) - sn*cphi
                     ELSE
                        x1 = XSLp(islp) + ABS(sn*cphi)
                     ENDIF
                     y1 = -(YSLp(islp)+sn*sphi)
                  ENDIF
                  WRITE (ifile,'(3(1X,E14.6))') x1 , y1 , 0.0D0
               ENDDO
            ENDDO
            WRITE (ifile,*)
            WRITE (ifile,'(A5,1X,I7,1X,I7)') 'CELLS' , tot_pts , 2*tot_pts
            DO i = 1 , tot_pts
               WRITE (ifile,'(2(1X,I7))') 1 , i - 1
            ENDDO
            WRITE (ifile,*)
            WRITE (ifile,'(A10,1X,I7)') 'CELL_TYPES' , tot_pts
            DO i = 1 , tot_pts
               WRITE (ifile,FMT='(5(1x,I7))') 1
            ENDDO
            WRITE (ifile,*) 'POINT_DATA' , tot_pts
            WRITE (ifile,*) 'SCALARS SO integer 1'
            WRITE (ifile,*) 'LOOKUP_TABLE default'
            DO islp = 1 , NSLp
               DO i = 1 , NNUc(islp)
                  WRITE (ifile,'(I7)') 1
               ENDDO
               DO i = 1 , NOBs(islp)
                  WRITE (ifile,'(I7)') 2
               ENDDO
            ENDDO


!!$$$$      write(999,*) 'Variables = X Y'
!!$$$$
!!$$$$      kk = 0
!!$$$$      write(ifile, *) 'Zone T = nuc1'
!!$$$$      do islp = 1,nslp
!!$$$$         if (nnuc(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            if (li .eq. 1) then
!!$$$$               cphi = cosphi(li)
!!$$$$               sphi = sinphi(li)
!!$$$$               tphi = sphi/cphi
!!$$$$               xstart = xslp(islp)
!!$$$$               xend = xendslp(islp)
!!$$$$               ystart = yslp(islp)
!!$$$$               yend = yendslp(islp)
!!$$$$               do i = 1,nnuc(islp)
!!$$$$                  si = snuc(i,islp)
!!$$$$                  if ( yend > 0) then
!!$$$$                     x = xstart + cphi*si
!!$$$$                     y = ystart + sphi*si
!!$$$$                  else
!!$$$$                     x = xstart - cphi*si
!!$$$$                     y = ystart - sphi*si
!!$$$$                  endif
!!$$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
!!$$$$               end do
!!$$$$            endif
!!$$$$         endif
!!$$$$      end do
!!$$$$      write(ifile, *) 'Zone T = nuc2'
!!$$$$
!!$$$$      do islp = 1,nslp
!!$$$$         if (nnuc(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            if (li .eq. 2) then
!!$$$$               cphi = cosphi(li)
!!$$$$               sphi = sinphi(li)
!!$$$$               tphi = sphi/cphi
!!$$$$               xstart = xslp(islp)
!!$$$$               xend = xendslp(islp)
!!$$$$               ystart = yslp(islp)
!!$$$$               yend = yendslp(islp)
!!$$$$               do i = 1,nnuc(islp)
!!$$$$                  si = snuc(i,islp)
!!$$$$                  if ( yend > 0) then
!!$$$$                     x = xstart + cphi*si
!!$$$$                     y = ystart + sphi*si
!!$$$$                  else
!!$$$$                     x = xstart - cphi*si
!!$$$$                     y = ystart - sphi*si
!!$$$$                  endif
!!$$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
!!$$$$               end do
!!$$$$            endif
!!$$$$         endif
!!$$$$      end do
!!$$$$      write(ifile, *) 'Zone T = nuc3'
!!$$$$
!!$$$$      do islp = 1,nslp
!!$$$$         if (nnuc(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            if (li .eq. 3) then
!!$$$$               cphi = cosphi(li)
!!$$$$               sphi = sinphi(li)
!!$$$$               tphi = sphi/cphi
!!$$$$               xstart = xslp(islp)
!!$$$$               xend = xendslp(islp)
!!$$$$               ystart = yslp(islp)
!!$$$$               yend = yendslp(islp)
!!$$$$               do i = 1,nnuc(islp)
!!$$$$                  si = snuc(i,islp)
!!$$$$                     x = xstart + cphi*si
!!$$$$                     y = ystart + sphi*si
!!$$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
!!$$$$               end do
!!$$$$            endif
!!$$$$         endif
!!$$$$      end do
!!$$$$
!!$$$$
!!$$$$      write(ifile, *) 'Zone T = obs1'
!!$$$$      do islp = 1,nslp
!!$$$$         if (nobs(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            if (li .eq. 1) then
!!$$$$               cphi = cosphi(li)
!!$$$$               sphi = sinphi(li)
!!$$$$               tphi = sphi/cphi
!!$$$$               xstart = xslp(islp)
!!$$$$               xend = xendslp(islp)
!!$$$$               ystart = yslp(islp)
!!$$$$               yend = yendslp(islp)
!!$$$$               do i = 1,nobs(islp)
!!$$$$                  si = sobs(i,islp)
!!$$$$                  if ( yend > 0) then
!!$$$$                     x = xstart + cphi*si
!!$$$$                     y = ystart + sphi*si
!!$$$$                  else
!!$$$$                     x = xstart - cphi*si
!!$$$$                     y = ystart - sphi*si
!!$$$$                  endif
!!$$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
!!$$$$               end do
!!$$$$            endif
!!$$$$         endif
!!$$$$      end do
!!$$$$      write(ifile, *) 'Zone T = obs2'
!!$$$$
!!$$$$      do islp = 1,nslp
!!$$$$         if (nobs(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            if (li .eq. 2) then
!!$$$$               cphi = cosphi(li)
!!$$$$               sphi = sinphi(li)
!!$$$$               tphi = sphi/cphi
!!$$$$               xstart = xslp(islp)
!!$$$$               xend = xendslp(islp)
!!$$$$               ystart = yslp(islp)
!!$$$$               yend = yendslp(islp)
!!$$$$               do i = 1,nobs(islp)
!!$$$$                  si = sobs(i,islp)
!!$$$$                  if ( yend > 0) then
!!$$$$                     x = xstart + cphi*si
!!$$$$                     y = ystart + sphi*si
!!$$$$                  else
!!$$$$                     x = xstart - cphi*si
!!$$$$                     y = -(ystart + sphi*si)
!!$$$$                  endif
!!$$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
!!$$$$               end do
!!$$$$            endif
!!$$$$         endif
!!$$$$      end do
!!$$$$      write(ifile, *) 'Zone T = obs3'
!!$$$$
!!$$$$      do islp = 1,nslp
!!$$$$         if (nobs(islp) > 0) then
!!$$$$            li =locphi(islp)
!!$$$$            if (li .eq. 3) then
!!$$$$               cphi = cosphi(li)
!!$$$$               sphi = sinphi(li)
!!$$$$               tphi = sphi/cphi
!!$$$$               xstart = xslp(islp)
!!$$$$               xend = xendslp(islp)
!!$$$$               ystart = yslp(islp)
!!$$$$               yend = yendslp(islp)
!!$$$$               do i = 1,nobs(islp)
!!$$$$                  si = sobs(i,islp)
!!$$$$                     x = xstart + cphi*si
!!$$$$                     y = ystart + sphi*si
!!$$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
!!$$$$               end do
!!$$$$            endif
!!$$$$         endif
!!$$$$      end do








            CLOSE (ifile)
      END SUBROUTINE PLOT_SLIP

      SUBROUTINE GEN_OBSTACLES
            IMPLICIT NONE
!!$*--GEN_OBSTACLES1913
            DOUBLE PRECISION :: XE , XNU , XLAmbda , XMU , pi
            DOUBLE PRECISION :: CC(6,6)
            INTEGER :: I_Elas , inuc
            COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
            INTEGER :: i , j , k , l , li , ii , islp , iseed1 , iseed2 , &
                 &           fe_locate
            LOGICAL :: accept
            DOUBLE PRECISION :: ev_convert , e , mu , nu , lambda , lnuc
            DOUBLE PRECISION :: tobs , dlmax , dlengthslp
            DOUBLE PRECISION :: x , x1 , rand(2)
            DOUBLE PRECISION :: sphi , cphi , xstart , ystart , xend , yend
            DOUBLE PRECISION :: x_lobs_nuc , xhold , xtemp
            INTEGER :: nobs_slip

            pi = 3.1415926535898D0
            ev_convert = 1.602176462D0
            tobs = 1500.0D0 !!$ MPa
            tobs = tobs*ev_convert*1.D-5
            LOBs = 2000.0D0 !!$ Angstrom
!!$     take care the lobs > 2*Lnuc
            LOBs_max = 1.5D0*LOBs
            LOBs_min = 0.5D0*LOBs
            ii = 0
            iseed1 = 111
            iseed2 = 222
            CALL RMARIN(iseed1,iseed2)

            DO islp = 1 , NSLp
               i = 0
               li = LOCphi(islp)
               cphi = COSphi(li)
               sphi = SINphi(li)
               xstart = XSLp(islp)
               ystart = YSLp(islp)
               xend = XENdslp(islp)
               yend = YENdslp(islp)
               IF ( LOCphi(islp)<3 ) THEN
                  dlmax = ABS(yend-ystart)/ABS(sphi)
               ELSE
                  dlmax = xend - xstart
               ENDIF
!!$ ------ Put special conditions for slip-plane 3
               IF ( li/=3 ) THEN
                  x = (ATOm_ymax+15.0D0)/sphi
               ELSE
                  x = 0.0D0
               ENDIF
               IF ( NNUc(islp)>0 ) THEN
                  IF ( dlmax>2.D0*LOBs ) THEN
                     DO WHILE ( x<dlmax )
!!$call ranmar(rand,2)
                        accept = .FALSE.
                        k = 0
                        DO WHILE ( accept==.FALSE. )
                           CALL RANMAR(rand,2)
                           IF ( i>0 ) THEN
                              x1 = LOBs + LOBs*(rand(1)-0.5)
                           ELSE
                              x1 = LOBs_min*rand(1)
                           ENDIF
                           x = x + x1
                           accept = .TRUE.
                           CALL XMIN1(x,SNUc(1:NNUc(islp),islp),NNUc(islp),&
                                &                          xhold,inuc)
                           x_lobs_nuc = MAX(LOBs/3.0D0,XLNuc(inuc,islp))
!!$$$$                     write( *,fmt='(A20,2I7,5(1X,E15.9))')
!!$$$$     $                    'obstacle trial ', islp,k, x,
!!$$$$     $                    xhold,snuc(inuc,islp), x_lobs_nuc
                           IF ( xhold<x_lobs_nuc ) THEN
                              accept = .FALSE.
                              x = x - x1
                              k = k + 1
                              IF ( k>20 ) x = x + LOBs_min
                           ENDIF
                           IF ( li==3 ) THEN
                              IF ( ABS(ystart)<ATOm_ymax+15.0D0 ) THEN
                                 xtemp = xstart + x
                                 IF ( xtemp>=ATOm_xmin .AND. &
                                      &                          xtemp<=ATOm_xmax ) accept = .FALSE.
                              ENDIF
                           ENDIF
                        ENDDO
                        IF ( x>dlmax ) EXIT
                        i = i + 1
                        NOBs(islp) = NOBs(islp) + 1
                        SOBs(i,islp) = x
                        TAU_obs(i,islp) = tobs
                        PRINT * , islp , i , x , TAU_obs(i,islp)
                     ENDDO
                  ENDIF
                  PRINT * , 'Total obstatcles on ' , islp , ' = ' , NOBs(islp)
               ENDIF
            ENDDO

      END SUBROUTINE GEN_OBSTACLES

      SUBROUTINE XMIN1(X,Xlist,N,Xhold,Inuc)
            IMPLICIT NONE
!!$*--XMIN12012
            DOUBLE PRECISION :: Xlist(N) , X , Xhold , d
            INTEGER :: Inuc , N , i
            Xhold = ABS(X-Xlist(1))
            Inuc = 1
            DO i = 2 , N
               d = ABS(X-Xlist(i))
               IF ( d<Xhold ) THEN
                  Xhold = d
                  Inuc = i
               ENDIF
            ENDDO

      END SUBROUTINE XMIN1


      SUBROUTINE RGAUSS(X,Len,Xmean,Sd)
!!$     !!$$  --------------------------------------------------------------
!!$     -------------
!!$     !!$$  Generate a random array of variables x with a Gaussian
!!$     distribution
!!$     !!$$  with a mean value xmean and a standard deviation sd. The
!!$     routine uses a trick
!!$     !!$$  found in
!!$     !!$$  Roes, P.B.M Van Oorschot, H.J.L : Kansrekening and statistek.
!!$     DUM 1978
!!$     !!$$  The routine makes use of ranmar to generate random numbers and
!!$     assumes
!!$     !!$$  ramrin has been called before to initialize ranmar
!!$     !!$$  --------------------------------------------------------------
!!$     -------------
            IMPLICIT NONE
!!$*--RGAUSS2044
!!$*** Start of declarations inserted by SPAG
            INTEGER i , ix , n
            DOUBLE PRECISION scale , v1 , v2
!!$*** End of declarations inserted by SPAG
            DOUBLE PRECISION :: pi
            INTEGER :: Len
            DOUBLE PRECISION :: rand(2) , Xmean , Sd , X(Len)
            DOUBLE PRECISION :: pi2
            pi = 3.1415926535898
            pi2 = 2.D0*pi
            n = Len/2
            IF ( (2*n)/=Len ) n = n + 1
            ix = 0
            DO i = 1 , n
               CALL RANMAR(rand,2)
               v1 = rand(1)
               v2 = rand(2)
               scale = Sd*SQRT(-2.0D0*LOG(v1))
               ix = ix + 1
               X(ix) = Xmean + scale*COS(pi2*v2)
               ix = ix + 1
               IF ( ix>Len ) RETURN
               X(ix) = Xmean + scale*SIN(pi2*v2)
            ENDDO
      END SUBROUTINE RGAUSS

!!$$$$      function xmin1(x,xlist,n)
!!$$$$      implicit double precision (a-h, o-z)
!!$$$$      double precision :: xlist(n), x
!!$$$$!!$     !!$$    print *, 'xmin1'
!!$$$$!!$     !!$$    do i = 1,n
!!$$$$!!$     !!$$       print *, i, xlist(i)
!!$$$$!!$     !!$$    end do
!!$$$$      xmin1 = abs(x-xlist(1))
!!$$$$      do i = 2, n
!!$$$$         d = abs(x-xlist(i))
!!$$$$         if (d < dmin) then
!!$$$$            xmin1 = d
!!$$$$         end if
!!$$$$      end do
!!$$$$      return
!!$$$$      end function xmin1



!!$ This random number generator originally appeared in "Toward a Universa
!!$ Random Number Generator" by George Marsaglia and Arif Zaman.
!!$ Florida State University Report: FSU-SCRI-87-50 (1987)
!!$
!!$ It was later modified by F. James and published in "A Review of Pseudo
!!$ random Number Generators"
!!$
!!$ THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
!!$       (However, a newly discovered technique can yield
!!$         a period of 10^600. But that is still in the development stage
!!$
!!$ It passes ALL of the tests for random number generators and has a peri
!!$   of 2^144, is completely portable (gives bit identical results on all
!!$   machines with at least 24-bit mantissas in the floating point
!!$   representation).
!!$
!!$ The algorithm is a combination of a Fibonacci sequence (with lags of 9
!!$   and 33, and operation "subtraction plus one, modulo one") and an
!!$   "arithmetic sequence" (using subtraction).
!!$
!!$ On a Vax 11/780, this random number generator can produce a number in
!!$    13 microseconds.
!!$=======================================================================

      SUBROUTINE RMARIN(Ij,Kl)
            IMPLICIT NONE
!!$*--RMARIN2116
!!$*** Start of declarations inserted by SPAG
            INTEGER i , ii , Ij , j , jj , k , Kl , l , m
            REAL s , t
!!$*** End of declarations inserted by SPAG
!!$ This is the initialization routine for the random number generator RAN
!!$ NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!!$                                                      0 <= KL <= 30081
!!$ The random number sequences created by these two seeds are of sufficie
!!$ length to complete an entire calculation with. For example, if sveral
!!$ different groups are working on different parts of the same calculatio
!!$ each group could be assigned its own IJ seed. This would leave each gr
!!$ with 30000 choices for the second seed. That is to say, this random
!!$ number generator can create 900 million different subsequences -- with
!!$ each subsequence having a length of approximately 10^30.


            DOUBLE PRECISION :: U(97) , C , CD , CM
            INTEGER I97 , J97
            LOGICAL TESt
            COMMON /RASET1/ U , C , CD , CM , I97 , J97 , TESt

            TESt = .FALSE.

            IF ( Ij<0 .OR. Ij>31328 .OR. Kl<0 .OR. Kl>30081 ) THEN
               PRINT '(A)' , &
                    &' The first random number seed must have a value between 0 and 313&
                    &28'
               PRINT '(A)' , &
                    &         ' The second seed must have a value between 0 and 30081'
               STOP
            ENDIF

            i = MOD(Ij/177,177) + 2
            j = MOD(Ij,177) + 2
            k = MOD(Kl/169,178) + 1
            l = MOD(Kl,169)

            DO ii = 1 , 97
               s = 0.0
               t = 0.5
               DO jj = 1 , 24
                  m = MOD(MOD(i*j,179)*k,179)
                  i = j
                  j = k
                  k = m
                  l = MOD(53*l+1,169)
                  IF ( MOD(l*m,64)>=32 ) s = s + t
                  t = 0.5*t
               ENDDO
               U(ii) = s
            ENDDO

            C = 362436.0/16777216.0
            CD = 7654321.0/16777216.0
            CM = 16777213.0/16777216.0

            I97 = 97
            J97 = 33

            TESt = .TRUE.
      END SUBROUTINE RMARIN


      SUBROUTINE RANMAR(Rvec,Len)
            IMPLICIT NONE
!!$*--RANMAR2182
!!$*** Start of declarations inserted by SPAG
            INTEGER Len
            REAL uni
!!$*** End of declarations inserted by SPAG
!!$ This is the random number generator proposed by George Marsaglia in
!!$ Florida State University Report: FSU-SCRI-87-50
!!$ It was slightly modified by F. James to produce an array of pseudorand
!!$ numbers.
            DOUBLE PRECISION :: Rvec(*)
            DOUBLE PRECISION :: U(97) , C , CD , CM
            INTEGER I97 , J97
            LOGICAL TESt
            COMMON /RASET1/ U , C , CD , CM , I97 , J97 , TESt

            INTEGER ivec

            IF ( .NOT.TESt ) THEN
               PRINT '(A)' , &
                    &         ' Call the init routine (RMARIN) before calling RANMAR'
               STOP
            ENDIF

            DO ivec = 1 , Len
               uni = U(I97) - U(J97)
               IF ( uni<0.0 ) uni = uni + 1.0
               U(I97) = uni
               I97 = I97 - 1
               IF ( I97==0 ) I97 = 97
               J97 = J97 - 1
               IF ( J97==0 ) J97 = 97
               C = C - CD
               IF ( C<0.0 ) C = C + CM
               uni = uni - C
               IF ( uni<0.0 ) uni = uni + 1.0
               Rvec(ivec) = uni
            ENDDO
      END SUBROUTINE RANMAR

      SUBROUTINE ASSIGN_DISLOC_ONLY(Dd1)
            IMPLICIT NONE
!!$*--ASSIGN_DISLOC_ONLY2223
            TYPE (DD)  , DIMENSION(TOT_disl)::Dd1

!!$-----------------------------------------------------
!!$ Local variables
!!$------------------------
            INTEGER :: islp , i , j , k , l , ntot , ntotmax , iib
            DOUBLE PRECISION :: x , y , s , si
            DOUBLE PRECISION :: maxx , maxy , minx , miny , pii
!!$-------------------------------------------------------
            maxx = 0.D0
            RANge_ = 2.D0*PROcess_xmax + 1000.0D0
            pii = 2.D0*ASIN(1.D0)
            ntotmax = MXNSLP*(MXNDIS+MXNNUC)
            PRINT * , 'Size of dd1 entering assign is ' , SIZE(Dd1) , TOT_disl
            k = 1
            PRINT * , 'Dislocations are '
!!$ ----- Read in the dislocations-------------------------------------
            DO islp = 1 , NSLp
               IF ( NDIs(islp)>0 ) THEN
                  DO i = 1 , NDIs(islp)
!!$               if (elem_dis(i,islp) .ne. 0) then
                     iib = LOCphi(islp)
                     si = SDIs(i,islp)
!!$$$$               x = xslp (islp) + si * cosphi (iib)
!!$$$$               y = yslp (islp) + si * sinphi (iib)
                     x = RDIs(1,i,islp)
                     y = RDIs(2,i,islp)
                     Dd1(k)%XY = DCMPLX(x,y)
                     Dd1(k)%SLIPANGLE = PHIslp(iib)
                     Dd1(k)%BURGERS_DD(1) = B_Dd(i,islp)
                     Dd1(k)%BURGERS_DD(2) = B_Dd(i,islp)*COSphi(iib)
                     Dd1(k)%BURGERS_DD(3) = B_Dd(i,islp)*SINphi(iib)
                     Dd1(k)%SLIPPLANE_NO = islp
                     Dd1(k)%STYPE = 1
                     Dd1(k)%VELOCITY = VPRev(i,islp)&
                          &                           *CMPLX(COSphi(iib),SINphi(iib))
                     PRINT * , i , islp , SDIs(i,islp) , IOBpin(i,islp)
                     k = k + 1
!!$               if (k >= tot_disl ) print *,k
!!$               end do
                  ENDDO
               ENDIF
            ENDDO
            PRINT * , 'End Dislocations'
            PRINT * , &
                 &      '-------------------------------------------------------'
!!$tot_disloc = k-1
      END SUBROUTINE ASSIGN_DISLOC_ONLY

      SUBROUTINE ASSIGN_SOURCE_ONLY(Dd1)
            IMPLICIT NONE
!!$*--ASSIGN_SOURCE_ONLY2275
            TYPE (DD)  , DIMENSION(TOT_size)::Dd1
!!$----------------------------------------------------
!!$ Local variables
!!$-----------------------------------------------------
            INTEGER :: islp , i , j , k , l , ntot , ntotmax , iib
            DOUBLE PRECISION :: x , y , s , si
!!$----------------------------------------------------
            k = TOT_disl + 1
            DO islp = 1 , NSLp
               IF ( NNUc(islp)>0 ) THEN
                  DO i = 1 , NNUc(islp)
                     iib = LOCphi(islp)
                     si = SNUc(i,islp)
!!$$$$               x = xslp (islp) + si * cosphi (iib)
!!$$$$               y = yslp (islp) + si * sinphi (iib)
                     x = RNUc(1,i,islp)
                     y = RNUc(2,i,islp)
!!$$$$               if (yendslp(islp) > 0.0d0) then
!!$$$$                  x = xslp(islp) + si*cosphi(iib)
!!$$$$                  y = yslp(islp) + si*sinphi(iib)
!!$$$$               else
!!$$$$                  if (iib == 1) then
!!$$$$                     x = xslp(islp) - si*cosphi(iib)
!!$$$$                  else
!!$$$$                     x = xslp(islp) + abs(si*cosphi(iib))
!!$$$$                  end if
!!$$$$                  y = -(yslp(islp) + si*sinphi(iib))
!!$$$$               end if
                     IF ( y==0 ) y = 1.E-9
                     IF ( x==0 ) x = 1.E-9
                     Dd1(k)%XY = CMPLX(x,y)
                     Dd1(k)%SLIPANGLE = PHIslp(iib)
                     Dd1(k)%BURGERS_DD(1) = BLEn
                     Dd1(k)%BURGERS_DD(2) = BLEn*COSphi(iib)
                     Dd1(k)%BURGERS_DD(3) = BLEn*SINphi(iib)
                     Dd1(k)%SLIPPLANE_NO = islp
                     Dd1(k)%STYPE = 2
                     Dd1(k)%VELOCITY = 0.D0
                     k = k + 1
                  ENDDO
               ENDIF
            ENDDO
            PRINT * , 'Total size =' , TOT_size , k - 1
!!$tot_size = k -1;
      END SUBROUTINE ASSIGN_SOURCE_ONLY

      SUBROUTINE CALCFORCE(Dd1,Start,Finish)
            IMPLICIT NONE
!!$*--CALCFORCE2324
            DOUBLE PRECISION :: XE , XNU , XLAmbda , XMU
            DOUBLE PRECISION :: CC(6,6)
            INTEGER :: I_Elas
            COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
            TYPE (DD)  , DIMENSION(:) , INTENT(INOUT)::Dd1
            DOUBLE COMPLEX :: b , bc , t , ct , t1 , ct1 , tmp1 , im , tc , &
                 &                  zph , zmh , czmh , czph , zstar , z2 , z3 , zx
            DOUBLE COMPLEX :: phi1 , phi2 , phi11 , z , z0 , zc , z0c , &
                 &                  phi2c , zi , z0i , zic , z0ic , th , tmp2 , dzds
            DOUBLE PRECISION :: sig11 , sig12 , sig22 , bi , sin2i , cos2i , &
                 &                    fg , fc , fgg , s11 , s12 , s22
            DOUBLE COMPLEX :: phi1z , phi2z , phi11z , tmp1z , tmp2z , gz
            DOUBLE PRECISION :: sig11z , sig12z , sig22z
            INTEGER :: i , j , k , Start , start1 , Finish , finish1
            DOUBLE PRECISION :: e , mu , nu , factor , h , rcore , x , y , &
                 &                    x0 , pii
            DOUBLE PRECISION , DIMENSION(2) :: zstar1(2) , zstar2(2) , z1(2)

            rcore = 2.D0*BLEn
            pii = 2.D0*ASIN(1.D0)
!!$$$$      start = 1
!!$$$$      finish = size(dd1)
            start1 = 1
            finish1 = TOT_disl
            im = CMPLX(0.,1.)
            mu = XMU
            nu = XNU
            factor = mu/(4.*pii*(1-nu))
!!$     print *, 'Factor = ', factor
            DO i = Start , Finish
               Dd1(i)%FORCED = 0.D0
               Dd1(i)%FORCEDG = 0.D0
               Dd1(i)%phi1 = 0.D0
               Dd1(i)%DPHI1 = 0.D0
               phi1 = 0.D0
               phi11 = 0.0d0;
               phi1z = 0.D0
               sig12z = 0.D0

               z = Dd1(i)%XY
               zc = CONJG(z)
!!$     write (*,fmt='(I2,1x,2E12.5,2x,2e12.5)')i, z,
!!$     cmplx(dd1(i)%burgers_dd(2),dd1(i)%burgers_dd(3))
               sig11 = 0.
               DO j = start1 , finish1
                  tmp1z = 0.D0
                  z0 = Dd1(j)%XY
                  z0c = CONJG(z0)
                  b = CMPLX(Dd1(j)%BURGERS_DD(2),Dd1(j)%BURGERS_DD(3))
                  bc = CONJG(b)

                  IF ( ABS(z-z0)>2*BLEn ) THEN
                     phi1 = -im*factor*b/(z-z0)
                     phi11 = im*factor*b/((z-z0)*(z-z0))
                     phi2 = im*factor*bc/(z-z0)
                     tmp1 = 2.D0*(phi1+CONJG(phi1))
                     tmp2 = -2.D0*((z-z0)*CONJG(phi11)+CONJG(phi2))


                     phi1z = im*factor*b/((z-z0)*(z-z0))
                     phi11z = -2.D0*im*factor*b/((z-z0)*(z-z0)*(z-z0))
                     phi2z = -im*factor*bc/((z-z0)*(z-z0))
                     tmp1z = 2.D0*(phi1z+CONJG(phi1z))
                     tmp2z = -2.D0*((z-z0)*CONJG(phi11z)+CONJG(phi2z))

                     sig11z = sig11z + REAL(0.5D0*(tmp1z+tmp2z))
                     sig22z = sig22z + REAL(0.5D0*(tmp1z-tmp2z))
                     sig12z = sig12z + AIMAG(0.5D0*(tmp2z))
                     sig11 = sig11 + REAL(0.5D0*(tmp1+tmp2))
                     sig22 = sig22 + REAL(0.5D0*(tmp1-tmp2))
                     sig12 = sig12 + AIMAG(0.5D0*(tmp2))

                  ELSE
                     t = z - z0
                     IF ( i/=j ) THEN
                        IF ( REAL(t)==0 ) THEN
                           th = DCMPLX(0,ASIN(SIN(Dd1(j)%SLIPANGLE)))
                        ELSE
                           th = DCMPLX(0,AIMAG(LOG(t)))
                        ENDIF
                        zstar = rcore*EXP(th)
                        phi1 = -im*factor*b/zstar
                        phi11 = im*factor*b/(zstar*zstar)
                        phi2 = im*factor*bc/(zstar)
                        tmp1 = 2.D0*(phi1+CONJG(phi1))
                        tmp2 = -2.D0*((z-z0)*CONJG(phi11)+CONJG(phi2))
                        sig11 = sig11 + REAL(0.5D0*(tmp1+tmp2))
                        sig22 = sig22 + REAL(0.5D0*(tmp1-tmp2))
                        sig12 = sig12 + AIMAG(0.5D0*(tmp2))
                     ENDIF
                  ENDIF
               ENDDO
               sin2i = SIN(2*Dd1(i)%SLIPANGLE)
               cos2i = COS(2*Dd1(i)%SLIPANGLE)
               bi = Dd1(i)%BURGERS_DD(1)
               fg = Dd1(i)%BURGERS_DD(1)*((sig22-sig11)*0.5*sin2i+sig12*cos2i)
               Dd1(i)%FORCED = Dd1(i)%FORCED + fg
               fgg = Dd1(i)%BURGERS_DD(1)&
                    &         *((sig22z-sig11z)*0.5*sin2i+sig12z*cos2i)
               Dd1(i)%FORCEDG = Dd1(i)%FORCEDG + fgg
!!$print *, 'xxxxx', i, dd1(i)%burgers_dd,dd1(i)%forced
            ENDDO
      END SUBROUTINE CALCFORCE


!!$
      SUBROUTINE DISLP(R,U)
!!$-----------------------------------------------------------------------
!!$     Determine displacement vector u at node n caused by dislocations
!!$
!!$     bout(1,islp) = total length of Burgers vectors of dislocations hav
!!$     moved out of the bottom side (1) of the cell
!!$-----------------------------------------------------------------------
!!$     Edition log:
!!$     06/10/94 The field of escaped dislocations is corrected to b/4
!!$-----------------------------------------------------------------------
            IMPLICIT NONE
!!$*--DISLP2440
            DOUBLE PRECISION U(3) , R(3)
            DOUBLE PRECISION dxx , dyy , r2 , rcore , dxi , dyi , cosi , &
                 &                 sini , r1
            DOUBLE PRECISION dx2 , dy2 , a , pi , t , u1 , u2 , dx0 , dy0 , &
                 &                 x , y , dxy2
            DOUBLE PRECISION XNU , XMU , CC(6,6) , XLAmbda , XE
            INTEGER I_Elas
            COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
            INTEGER i , j , islp , li
            DOUBLE PRECISION sign1(2)
            x = R(1)
            y = R(2)
            pi = 2.D0*ASIN(1.0D0)
            PRINT * , 'dislp' , NSLp , XNU

!!$
!!$
!!$.....For each node, loop over all active dislocations inside the cell
            U(1) = 0.0D0
            U(2) = 0.0D0
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               sini = SINphi(li)
               cosi = COSphi(li)
               dx0 = x - XSLp(islp)
               dy0 = y - YSLp(islp)
               dxi = dx0*cosi + dy0*sini
               dyi = -dx0*sini + dy0*cosi
               u1 = 0.0D0
               u2 = 0.0D0
               DO i = 1 , NDIs(islp)
                  rcore = 2.0D0*DABS(B_Dd(i,islp))
                  dxx = dxi - SDIs(i,islp)
                  dyy = dyi
                  dx2 = dxx*dxx
                  dy2 = dyy*dyy
                  r2 = dx2 + dy2
!!$     Cut off when closer than core radius
                  IF ( r2<rcore**2 ) THEN
                     r1 = SQRT(r2)
                     dxx = dxx/r1*rcore
                     dyy = dyy/r1*rcore
                     dx2 = dxx*dxx
                     dy2 = dyy*dyy
                     r2 = rcore**2
                  ENDIF
                  dxy2 = dxx**2 + dyy**2
                  a = B_Dd(i,islp)/2.0D0/pi/(1.0D0-XNU)
                  IF ( dyy==0.0D0 ) THEN
                     t = DSIGN(pi/2.0D0,dxx)
                  ELSE
                     t = DATAN(dxx/dyy)
                  ENDIF
!!$.....first the local displacements
                  u1 = u1 + a*(dxx*dyy/2.0D0/dxy2-(1.0D0-XNU)*t)
                  u2 = u2 + a*(dyy**2/2.0D0/dxy2-(1.0D0-2.0D0*XNU)&
                       &           /4.0D0*LOG(dxy2/B_Dd(i,islp)**2))
               ENDDO
!!$.....transform to global displacements
               U(1) = U(1) + cosi*u1 - sini*u2
               U(2) = U(2) + sini*u1 + cosi*u2
               U(3) = 0.0D0
            ENDDO

!!$.....Dislocations that have escaped from the cell previously:
!!$.....The displacement field generated by a + dislocation moved out at t
!!$.....right-hand side (2) is
!!$.....          b/4 for y>yd
!!$.....   u(1)=
!!$.....         -b/4 for y<yd
!!$.....Dislocations having moved out at the bottom side (1) generate
!!$.....displacements with the opposite sign.
            sign1(1) = -1.0D0
            sign1(2) = 1.0D0
            DO islp = 1 , NSLp
               li = LOCphi(islp)
               sini = SINphi(li)
               cosi = COSphi(li)
               dx0 = x - XSLp(islp)
               dy0 = y - YSLp(islp)
               dxi = dx0*cosi + dy0*sini
               dyi = -dx0*sini + dy0*cosi
               u1 = 0.0D0
               u2 = 0.0D0
!!$     sout(1,islp) = -1.d10
!!$     sout(2,islp) = 1.d10
!!$      soutq = 1.d10
               DO i = 1 , 2
                  IF ( BOUt(i,islp)/=0 ) THEN
!!$soutq = sign1(i)*soutq
                     rcore = 2.D0*BLEn
                     dxx = dxi - SDIs_out(i,islp)
                     dyy = dyi
                     dx2 = dxx*dxx
                     dy2 = dyy*dyy
                     r2 = dx2 + dy2
!!$       Cut off when closer than core radius
                     IF ( r2<rcore**2 ) THEN
                        r1 = DSQRT(r2)
                        dxx = dxx/r1*rcore
                        dyy = dyy/r1*rcore
                        dx2 = dxx*dxx
                        dy2 = dyy*dyy
                        r2 = rcore**2
                     ENDIF
                     dxy2 = dxx**2 + dyy**2
                     a = BOUt(i,islp)/2.0D0/pi/(1.0D0-XNU)
                     IF ( dyy==0.0D0 ) THEN
                        t = DSIGN(pi/2.0D0,dxx)
                     ELSE
                        t = DATAN(dxx/dyy)
                     ENDIF
!!$.....first the local displacements
                     IF ( dyy>0.0D0 ) THEN
                        u1 = u1 + sign1(i)*BOUt(i,islp)/4.0D0
                     ELSE
                        u1 = u1 - sign1(i)*BOUt(i,islp)/4.0D0
                     ENDIF
                     u2 = u2 + a*(dyy**2/2.0D0/dxy2-(1.0D0-2.0D0*XNU)&
                          &              /4.0D0*DLOG(dxy2/BOUt(i,islp)**2))
                  ENDIF
               ENDDO
!!$.....Now transform to global displacements
               U(1) = U(1) + cosi*u1 - sini*u2
               U(2) = U(2) + sini*u1 + cosi*u2

            ENDDO

      END SUBROUTINE DISLP


      SUBROUTINE RLSDIS(Tau)
            IMPLICIT NONE
!!$*--RLSDIS2574
!!$-----------------------------------------------------------------------
!!$ Release pinned dislocations when the resolved shear stress `tau'
!!$ exceeds a critical value
!!$ 08/04/01 dislocations pinned at junctions skipped at present
!!$-----------------------------------------------------------------------
            DOUBLE PRECISION :: Tau(MXNDIS,NSLp) , tauii
            INTEGER :: islp , i , j

            DO islp = 1 , NSLp
               DO i = 1 , NDIs(islp)
                  IF ( IOBpin(i,islp)>0 ) THEN
                     tauii = DABS(Tau(i,islp))
                     IF ( tauii>TAU_obs(IOBpin(i,islp),islp) ) THEN
                        WRITE (6,99001) i , islp , IOBpin(i,islp)
99001                   FORMAT (' Disloc. ',i3,' on slip plane ',i4,&
                             &                    ' released from obstacle ',i2)
                        IOBpin(i,islp) = 0
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
      END SUBROUTINE RLSDIS

      SUBROUTINE UNDREX(V,C)
!!$-----------------------------------------------------------------------
!!$ Determine underrelaxation factor for updating dislocation positions.
!!$ Currently, underrelaxation with a factor 0.75 is used when the
!!$ velocity of two like-signed dislocations (in pile-ups) changes sign
!!$ as compared to the previous increment.
!!$ In addition, underrelaxation with a factor 0.5 is used for ANY
!!$ sign change in velocity
!!$-----------------------------------------------------------------------
!!$ Edition log
!!$ 02/09/96 Underrelaxation implemented also for any change of sig
!!$-----------------------------------------------------------------------
            IMPLICIT NONE
!!$*--UNDREX2611
            DOUBLE PRECISION :: V(MXNDIS,NSLp) , C(MXNDIS,NSLp)
            INTEGER :: i , j , k , l , islp
            DOUBLE PRECISION :: bi
            DO islp = 1 , NSLp
               DO i = 1 , NDIs(islp)
                  C(i,islp) = 1.0D0
                  IF ( (V(i,islp)*VPRev(i,islp))<0.0D0 ) THEN
                     C(i,islp) = 0.5D0*C(i,islp)
                     bi = B_Dd(i,islp)
                     j = i - 1
                     IF ( j>=1 ) THEN
                        IF ( (bi*B_Dd(j,islp))>0.0D0 ) C(i,islp)&
                             &                 = 0.75D0*C(i,islp)
                     ENDIF
                     j = i + 1
                     IF ( j<=NDIs(islp) ) THEN
                        IF ( (bi*B_Dd(j,islp))>0.0D0 ) C(i,islp)&
                             &                 = 0.75D0*C(i,islp)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
      END SUBROUTINE UNDREX


END MODULE MOD_DD_SLIP
!!$*==tlce.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
FUNCTION TLCE(Xx)
      IMPLICIT NONE
!!$*--TLCE2641
!!$*** Start of declarations inserted by SPAG
      DOUBLE PRECISION TLCE , Xx
!!$*** End of declarations inserted by SPAG
      TLCE = DINT(Xx*1.0D6)/1.0D6
END FUNCTION TLCE
