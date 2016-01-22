!*==ma01.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
!     Qu modified on 02/01/2005---the displacement boundary conditions
!     are treated
!     according to  anisotropic elasticity theory. The major
!     modification is
!     employed in subroutine kfiled_displ()
      SUBROUTINE MA01(Id,X,Ix,F,B,Dr,Db,Input)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--MA0110
!
      CHARACTER Input*80
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*) , F(NDF,*) , Dr(*) , Db(*)
      INTEGER Id(NDF,*) , Ix(NEN1,*)
!
!     local variables
!
      DOUBLE PRECISION btmp(NDF) , lambda , mu
      INTEGER i , j
      LOGICAL fixed
 
      CALL FINDAVERAGEELASTICCONST(lambda,mu)
 
!
!     applied displacement b.c.
!
      PRINT * , '!!!!!!!!Calling k-field!!!!!!'
 
!!$      DO i = 1 , NUMnp
!!$!$$$         do j = 1, ndf
!!$!$$$            if (id(j,i) .eq. 1) then
!!$!$$$               if (isRelaxed(i) .eq. 1) then
!!$!$$$                  btmp(1) = 0.0
!!$!$$$                  btmp(2) = 10.0d0
!!$!$$$               end if
!!$!$$$!cc   --Only used when doing FEM, when b and f are obtained by scali
!!$!$$$!cc   --the current ones by time (that's why 1.0d0 is used below)
!!$!$$$               !call kfield_displ(1.0d0, x(1, i), btmp, lambda, mu, 
!!$!$$$
!!$!$$$!cc   --End of Mod
!!$!$$$            endif
!!$!$$$         end do
!!$         DO j = 1 , NDF
!!$            IF ( Id(j,i)==1 ) THEN
!!$               IF ( ISRelaxed(i)==1 ) THEN
!!$			     print *, 'atomdisp on ', i,j
!!$                  F(j,i) = 10.0D0
!!$               ELSE
!!$                  F(j,i) = 0.0
!!$               ENDIF
!!$            ENDIF
!!$         ENDDO
!!$      ENDDO
      do i = 1, numnp
         do j = 1, ndf
            if (id(j,i) == 1) then
               call kfield_displ(1.0d0, x(1,i), btmp, lambda, mu, 0.0d0)
            end if
         end do
         do j = 1, ndf
            if (id(j,i) == 1) then
               f(j,i) = btmp(j)
            end if
         end do
      end do
      END SUBROUTINE MA01
!*==precrack.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
 
 
      SUBROUTINE PRECRACK(Id,X,B,F,Input)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--PRECRACK60
!     c
      INTEGER Id(NDF,*)
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*) , F(NDF,*)
      DOUBLE PRECISION lambda , mu
!CC   --Hack parameter
      DOUBLE PRECISION yhshift
      INTEGER i
      CHARACTER*80 Input
!     c
      INTEGER lower , upper , NEXT , idum
      DOUBLE PRECISION dflag , kinit
      DOUBLE PRECISION btmp(NDF) , temp , btmp1(NDF)
      
!!$   Pulse parameters
      double precision :: a0, A, sigma, rc, delr, uc, r
      lower = 4
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,idum,kinit,2)
      lower = upper
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,idum,dflag,2)
      PRINT * , 'Time in precrack is' , TIMe , kinit , dflag

!!$!     c
!!$!     c
!!$!     Hard coded pulse for stadium langevin

      a0 = 4.032
      A = 0.86
      sigma = 5.0*a0/sqrt(2.0)
      rc = 10.0*a0;

      do i = 1, numnp
         if (isrelaxed(i) /=0) then
            if (isrelaxed(i) /= -1) then
               r = sqrt(X(1,i)**2 + X(2,i)**2)
               if (r < rc) then
                  uc = A*exp(-(rc/sigma)**2)
                  delr = A*(A*exp(-(r/sigma)**2)-uc)/(A-uc)
                  b(1,i) = b(1,i) + delr
                  b(2,i) = b(2,i) + delr
               end if 
            end if
         end if
      end do
      
!     c
!!$
!!$
!!$      CALL FINDAVERAGEELASTICCONST(lambda,mu)
!!$      DO i = 1 , NUMnp
!!$!CC   --JS Hack for H atoms
!!$         yhshift = 0.0D0
!!$         IF ( (ATOmspecie(i)==2) .AND. (X(1,i)<0.1) ) THEN
!!$            IF ( (X(2,i)>-0.2) .AND. (X(2,i)<0.2) ) yhshift = 0.1D0
!!$            IF ( (X(2,i)>1.0) .AND. (X(2,i)<2.0) ) yhshift = -0.1D0
!!$         ENDIF
!!$!CC   --END
!!$         IF ( dflag==0.0 ) THEN
!!$            TIMe = kinit
!!$            CALL KFIELD_DISPL(kinit,X(1,i),btmp,lambda,mu,yhshift)
!!$            B(1,i) = btmp(1)
!!$            B(2,i) = btmp(2)
!!$            IF ( Id(1,i)==1 ) F(1,i) = B(1,i)/TIMe
!!$            IF ( Id(2,i)==1 ) F(2,i) = B(2,i)/TIMe
!!$            IF ( Id(1,i)==1 .OR. Id(2,i)==1 ) THEN
!!$!     print *, 'K_Field in prec', i, f(1,i), f(2,i), x(1,i), x(2,i)
!!$            ENDIF
!!$         ELSE
!!$            CALL KFIELD_DISPL(kinit,X(1,i),btmp1,lambda,mu,yhshift)
!!$            temp = kinit + dflag
!!$            CALL KFIELD_DISPL(temp,X(1,i),btmp,lambda,mu,yhshift)
!!$            B(1,i) = B(1,i) - btmp1(1) + btmp(1)
!!$            B(2,i) = B(2,i) - btmp1(2) + btmp(2)
!!$            IF ( Id(1,i)==1 ) F(1,i) = B(1,i)/(TIMe)
!!$            IF ( Id(2,i)==1 ) F(2,i) = B(2,i)/(TIMe)
!!$         ENDIF
!!$ 
!!$      ENDDO
!!$ 

      !     c
      END SUBROUTINE PRECRACK
!*==pdelcalc.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
 
!**********************************************************************
!**********************************************************************
!**********************************************************************
 
!     c**---------------------------------------------------------------
!     c**   pdelcalc : compute the pdelta curve
!     c**
      SUBROUTINE PDELCALC(F,Dr,Id,Force,X)
      USE MOD_BOUNDARY
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--PDELCALC130
!     c
      DOUBLE PRECISION F(NDF,*) , Dr(NDF,*) , Force , X(NXDm,*)
      INTEGER Id(NDF,*)
!     c
      INTEGER i , j
!     c
      WRITE (6,*) 'p-delta is not implemented'
      END SUBROUTINE PDELCALC
!*==kfield_displ.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
 
 
 
      SUBROUTINE KFIELD_DISPL(K_u,Xin,Bout,Lambda,Mu,Yshift)
      USE MOD_GLOBAL
      USE MOD_CRACK
      IMPLICIT NONE
!*--KFIELD_DISPL147
      DOUBLE PRECISION M_PI
      PARAMETER (M_PI=3.141592653589793D0)
!     C
!     C  This is valid for the hexagonal lattice only. Elastic constants
!     should be
!     C  exported from the code!
!     C
      DOUBLE PRECISION Lambda
!     C      parameter (lambda = 0.5747d0) ! for 2D Hexagonal Al.
!     C	parameter (lambda = 0.3742) ! for 3D FCC Al.
!     C      parameter (lambda = 0.3784) ! for 3D FCC Al, Baskes and Daw
!     potential.
!     C      parameter (lambda = 0.7961) ! for 3D FCC Ni, Baskes and Daw
!     Potential
!     C
      DOUBLE PRECISION K_u , Xin(NXDm) , Bout(NDF)
!     c
      DOUBLE PRECISION Mu , nu , e , k , Yshift
      DOUBLE PRECISION x , y , r , theta
      DOUBLE PRECISION k_i , k_ii , ux1 , ux2 , uy1 , uy2
      DOUBLE PRECISION s1x , s1y , s2x , s2y , p1x , p1y , p2x , p2y
      DOUBLE PRECISION q1x , q1y , q2x , q2y
      DOUBLE PRECISION xnu , xmu
      COMPLEX*16 s1 , s2 , p1 , p2 , q1 , q2 , b1 , b2
 
      Mu = (0.9581D0-Lambda)/2.0D0 ! for 2D Hexagonal Al.
!      print *, 'mu = ', mu, ' lambda = ', lambda
!     C      mu = (0.7371 - lambda)/2.0d0 ! for 3D FCC Al.
!     C      mu = (0.7126 - lambda)/2.0d0 ! for 3D FCC Al, Baskes and
!     Daw potential.
!     C      mu = 0.3728 ! for 3D FCC Ni, Baskes and Daw Potential
 
      nu = Lambda/2.0D0/(Lambda+Mu)
      e = Mu*2.D0*(1.D0+nu)
      k = 3.D0 - 4.D0*nu
 
!      print *, 'mu =  ', mu,' nu = ', nu, ' K_u = ',K_u
!     C      stop
 
 
      x = Xin(1) - X0Crack
      y = Xin(2) - Y0Crack - Yshift
 
      r = DSQRT(x*x+y*y)
      theta = DATAN2(y,x)
!     c Qu modification begins
      k_i = K_u
      k_ii = 0.D0
      OPEN (UNIT=12,FILE='anisoEig.inp',STATUS='old')
      READ (12,*) s1x , s1y
      READ (12,*) s2x , s2y
      READ (12,*) p1x , p1y
      READ (12,*) p2x , p2y
      READ (12,*) q1x , q1y
      READ (12,*) q2x , q2y
      CLOSE (12)
      s1 = DCMPLX(s1x,s1y)
      s2 = DCMPLX(s2x,s2y)
      p1 = DCMPLX(p1x,p1y)
      p2 = DCMPLX(p2x,p2y)
      q1 = DCMPLX(q1x,q1y)
      q2 = DCMPLX(q2x,q2y)
      b1 = CDSQRT(DCOS(theta)+s1*DSIN(theta))
      b2 = CDSQRT(DCOS(theta)+s2*DSIN(theta))
      CALL ANISODISPL(s1,s2,p1,p2,q1,q2,b1,b2,ux1,ux2,uy1,uy2)
      Bout(1) = DSQRT(2.D0*r/M_PI)*(k_i*ux1+k_ii*ux2)
      Bout(2) = DSQRT(2.D0*r/M_PI)*(k_i*uy1+k_ii*uy2)
 
!     bout(1)=0.0*1e6/(2.0*.1925*100e9)*xin(2)
!     bout(2)=0.01*xin(2) + 0.0*1e6/(2.0*.1925*100e9)*xin(1)
 
!     xmu=mu*100e9/1e30/1.602e-19
!     xnu=nu
!     bout(1)=K_u/xmu*dsqrt(r/(2.0*M_PI))*
!     &     dcos(theta/2.0)*(1.0-2.0*xnu+dsin(theta/2.0)*dsin(theta
!     /2.0))
!     bout(2)=K_u/xmu*dsqrt(r/(2.0*M_PI))*
!     &     dsin(theta/2.0)*(2.0-2.0*xnu-dcos(theta/2.0)*dcos(theta
!     /2.0))
 
 
!     c
!     bout(1)=K_u/2.d0/e*dsqrt(r/2.d0/M_PI)*(1.d0+nu)
!     &     *((2.d0*k-1.0d0)*dcos(theta/2.0d0)-dcos(1.5d0*theta))
!     bout(2)=K_u/2.d0/e*dsqrt(r/2.d0/M_PI)*(1.d0+nu)
!     &     *((2.d0*k+1.0d0)*dsin(theta/2.0d0)-dsin(1.5d0*theta))
!     c
!     c Qu modification ends
      END SUBROUTINE KFIELD_DISPL
!*==anisodispl.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
 
 
      SUBROUTINE ANISODISPL(S1,S2,P1,P2,Q1,Q2,B1,B2,Ux1,Ux2,Uy1,Uy2)
      IMPLICIT NONE
!*--ANISODISPL242
!     c normalized displacement fields around a crack tip in anisotropic
!     material
!     c    ux=ux/(sqrt(2*r/PI)/K_u)
!     c    ux1--x-component under mod_I load
!     c    ux2--x-component under mod_II load
!     c    uy1--y-component under mod_I load
!     c    uy2--y-component under mod_II load
!     c
      COMPLEX*16 S1 , S2 , P1 , P2 , Q1 , Q2 , B1 , B2
      DOUBLE PRECISION Ux1 , Ux2 , Uy1 , Uy2
 
      Ux1 = DREAL((S1*P2*B2-S2*P1*B1)/(S1-S2))
      Ux2 = DREAL((P2*B2-P1*B1)/(S1-S2))
      Uy1 = DREAL((S1*Q2*B2-S2*Q1*B1)/(S1-S2))
      Uy2 = DREAL((Q2*B2-Q1*B1)/(S1-S2))
 
      END SUBROUTINE ANISODISPL
!*==findaverageelasticconst.spg  processed by SPAG 6.70Rc at 12:37 on 29
 
      SUBROUTINE FINDAVERAGEELASTICCONST(Lambda,Mu)
 
      USE MOD_MATERIAL
      IMPLICIT NONE
!*--FINDAVERAGEELASTICCONST266
      DOUBLE PRECISION Lambda , nu , Mu , cel(6,6)
!     C	type(bravaismat), dimension(:), pointer :: material
 
      cel = MATerial(1)%CC
 
      Lambda = 1.0/15.0*(cel(1,1)+cel(3,3)+5.0*cel(1,2)+8.0*cel(1,3)&
     &         -4.0*cel(4,4))
 
      Mu = 1.0/30.0*(7.0*cel(1,1)-5.0*cel(1,2)+2.0*cel(3,3)&
     &     +12.0*cel(4,4)-4.0*cel(1,3))
 
      END SUBROUTINE FINDAVERAGEELASTICCONST
 
 
