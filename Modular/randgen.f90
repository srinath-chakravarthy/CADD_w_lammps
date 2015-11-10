!*==zbqlbd01.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!******************************************************************
!*******	FILE: randgen.f				***********
!*******	AUTHORS: Richard Chandler		***********
!*******		 (richard@stats.ucl.ac.uk)	***********
!*******		 Paul Northrop 			***********
!*******		 (northrop@stats.ox.ac.uk)	***********
!*******	LAST MODIFIED: 1/7/03			***********
!*******	See file randgen.txt for details	***********
!******************************************************************
 
      BLOCKDATA ZBQLBD01
      IMPLICIT NONE
!*--ZBQLBD0114
!
!       Initializes seed array etc. for random number generator.
!       The values below have themselves been generated using the
!       NAG generator.
!
      COMMON /ZBQL0001/ ZBQlix , B , C
      INTEGER ZBQlix(43) , i , B , C
      DATA (ZBQlix(i),i=1,43)/8001441 , 55321801 , 169570999 , &
     &      288589940 , 291581871 , 103842493 , 79952507 , 381202335 , &
     &      311575334 , 402878631 , 249757109 , 115192595 , 210629619 , &
     &      399952890 , 412280521 , 133873288 , 71345525 , 223467704 , &
     &      282934796 , 99756750 , 168564303 , 286817366 , 114310713 , &
     &      347045253 , 93762426 , 109670477 , 320029657 , 326369301 , &
     &      9441177 , 353244738 , 244771580 , 159804337 , 207319904 , &
     &      337342907 , 375423178 , 70893571 , 426059785 , 395854390 , &
     &      20081010 , 59250059 , 162176640 , 320429173 , 263576576/
      DATA B/424967291/
      DATA C/0/
      END BLOCKDATA ZBQLBD01
!*==zbqlini.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
!*****************************************************************
!*****************************************************************
      SUBROUTINE ZBQLINI(Seed)
      IMPLICIT NONE
!*--ZBQLINI40
!*****************************************************************
!       To initialize the random number generator - either
!       repeatably or nonrepeatably. Need double precision
!       variables because integer storage can't handle the
!       numbers involved
!*****************************************************************
!	ARGUMENTS
!	=========
!	SEED	(integer, input). User-input number which generates
!		elements of the array ZBQLIX, which is subsequently used
!		in the random number generation algorithm. If SEED=0,
!		the array is seeded using the system clock if the
!		FORTRAN implementation allows it.
!*****************************************************************
!	PARAMETERS
!	==========
!	LFLNO	(integer). Number of lowest file handle to try when
!		opening a temporary file to copy the system clock into.
!		Default is 80 to keep out of the way of any existing
!		open files (although the program keeps searching till
!		it finds an available handle). If this causes problems,
!               (which will only happen if handles 80 through 99 are
!               already in use), decrease the default value.
!*****************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
!*****************************************************************
!	VARIABLES
!	=========
!	SEED	See above
!	ZBQLIX	Seed array for the random number generator. Defined
!		in ZBQLBD01
!	B,C	Used in congruential initialisation of ZBQLIX
!	SS,MM,}	System clock secs, mins, hours and days
!	HH,DD }
!	FILNO	File handle used for temporary file
!	INIT	Indicates whether generator has already been initialised
!
      INTEGER Seed , ZBQlix(43) , B , C , ss , mm , hh , dd , filno , i
      INTEGER init
      DOUBLE PRECISION tmpvar1 , dss , dmm , dhh , ddd , db
 
      COMMON /ZBQL0001/ ZBQlix , B , C
      SAVE init
 
!
!	Ensure we don't call this more than once in a program
!
      IF ( init>=1 ) THEN
         IF ( init==1 ) THEN
            WRITE (*,99001)
 
99001       FORMAT (//5X,&
     &              '****WARNING**** You have called routine ZBQLINI ',&
     &              'more than',/5X,&
     &              'once. I''m ignoring any subsequent calls.',//)
            init = 2
         ENDIF
         RETURN
      ELSE
         init = 1
      ENDIF
!
!       If SEED = 0, cat the contents of the clock into a file
!       and transform to obtain ZQBLIX(1), then use a congr.
!       algorithm to set remaining elements. Otherwise take
!       specified value of SEED.
!
      db = DBLE(B)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>	NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
!>>>>>>>	(NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
!>>>>>>>	THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
!>>>>>>>	OF THE FOLLOWING IF BLOCK SHOULD BE	 >>>>>>>
!>>>>>>>	COMMENTED OUT.				 >>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
!>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       CALL SYSTEM(' date +%S%M%H%j > zbql1234.tmp')
!*
!*       Try all file numbers for LFLNO to 999
!*
!       FILNO = LFLNO
! 10    OPEN(FILNO,FILE='zbql1234.tmp',ERR=11)
!       GOTO 12
! 11    FILNO = FILNO + 1
!       IF (FILNO.GT.999) THEN
!        WRITE(*,2)
!        RETURN
!       ENDIF
!       GOTO 10
! 12    READ(FILNO,'(3(I2),I3)') SS,MM,HH,DD
!       CLOSE(FILNO)
!       CALL SYSTEM('rm zbql1234.tmp')
!       DSS = DINT((DBLE(SS)/6.0D1) * DBLE(B))
!       DMM = DINT((DBLE(MM)/6.0D1) * DBLE(B))
!       DHH = DINT((DBLE(HH)/2.4D1) * DBLE(B))
!       DDD = DINT((DBLE(DD)/3.65D2) * DBLE(B))
!       TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,DB)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<
!<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF ( Seed/=0 ) tmpvar1 = DMOD(DBLE(Seed),db)
      ZBQlix(1) = INT(tmpvar1)
      DO i = 2 , 43
         tmpvar1 = DBLE(ZBQlix(i-1))*3.0269D4
         tmpvar1 = DMOD(tmpvar1,db)
         ZBQlix(i) = INT(tmpvar1)
      ENDDO
99002 FORMAT (//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',&
     &        ' find an',/5X,&
     &    'available file number. To rectify the problem, decrease the '&
     &    ,'value of',/5X,&
     &    'the parameter LFLNO at the start of this routine (in file ',&
     &    'randgen.f)',/5X,&
     &    'and recompile. Any number less than 100 should work.')
      END SUBROUTINE ZBQLINI
!*==zbqlu01.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLU01(Dummy)
      IMPLICIT NONE
!*--ZBQLU01165
!
!       Returns a uniform random number between 0 & 1, using
!       a Marsaglia-Zaman type subtract-with-borrow generator
!
      DOUBLE PRECISION ZBQLU01 , Dummy
      INTEGER C , B , ZBQlix(43) , x , curpos , id22 , id43
 
      COMMON /ZBQL0001/ ZBQlix , B , C
      SAVE /ZBQL0001/ 
      SAVE curpos , id22 , id43
      DATA curpos , id22 , id43/1 , 22 , 43/
 
      x = ZBQlix(id22) - ZBQlix(id43) - C
      IF ( x<0 ) THEN
         x = x + B
         C = 1
      ELSE
         C = 0
      ENDIF
      ZBQLU01 = DBLE(x)/DBLE(B)
      ZBQlix(id43) = x
!
!     Update array pointers. Do explicit check for bounds of each to
!     avoid expense of modular arithmetic. If one of them is 0 the other
!     won't be
!
      curpos = curpos - 1
      id22 = id22 - 1
      id43 = id43 - 1
      IF ( curpos==0 ) THEN
         curpos = 43
      ELSEIF ( id22==0 ) THEN
         id22 = 43
      ELSEIF ( id43==0 ) THEN
         id43 = 43
      ENDIF
 
      END FUNCTION ZBQLU01
!*==zbqluab.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLUAB(A,B)
      IMPLICIT NONE
!*--ZBQLUAB208
!
!       Returns a random number uniformly distributed on (A,B)
!
      DOUBLE PRECISION A , B , ZBQLU01 , ZBQLUAB
 
!
!       Even if A > B, this will work as B-A will then be -ve
!
      IF ( A/=B ) THEN
         ZBQLUAB = A + ((B-A)*ZBQLU01(0.0D0))
      ELSE
         ZBQLUAB = A
         WRITE (*,99001)
99001    FORMAT (/5X,'****WARNING**** (function ZBQLUAB) Upper and ',&
     &           'lower limits on uniform',/5X,&
     &           'distribution are identical',/)
      ENDIF
      END FUNCTION ZBQLUAB
!*==zbqlexp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLEXP(Mu)
      IMPLICIT NONE
!*--ZBQLEXP231
!
!       Returns a random number exponentially distributed with
!       mean MU
!
      DOUBLE PRECISION Mu , ZBQLEXP , ZBQLU01
 
      ZBQLEXP = 0.0D0
 
      IF ( Mu<0.0D0 ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLEXP',/)
         RETURN
      ENDIF
 
      ZBQLEXP = -DLOG(ZBQLU01(0.0D0))*Mu
 
      END FUNCTION ZBQLEXP
!*==zbqlnor.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLNOR(Mu,Sigma)
      IMPLICIT NONE
!*--ZBQLNOR255
!
!       Returns a random number Normally distributed with mean
!       MU and standard deviation |SIGMA|, using the Box-Muller
!       algorithm
!
      DOUBLE PRECISION theta , r , ZBQLNOR , ZBQLU01 , pi , Mu , Sigma
      DOUBLE PRECISION spare
      INTEGER status
      SAVE status , spare , pi
      DATA status/ - 1/
 
      IF ( status==-1 ) pi = 4.0D0*DATAN(1.0D0)
 
      IF ( status<=0 ) THEN
         theta = 2.0D0*pi*ZBQLU01(0.0D0)
         r = DSQRT(-2.0D0*DLOG(ZBQLU01(0.0D0)))
         ZBQLNOR = (r*DCOS(theta))
         spare = (r*DSIN(theta))
         status = 1
      ELSE
         ZBQLNOR = spare
         status = 0
      ENDIF
 
      ZBQLNOR = Mu + (Sigma*ZBQLNOR)
 
      END FUNCTION ZBQLNOR
!*==zbqlbin.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLBIN(N,P)
      IMPLICIT NONE
!*--ZBQLBIN287
!
!       Returns a random number binomially distributed (N,P)
!
      DOUBLE PRECISION P , ZBQLBET1
      DOUBLE PRECISION pp , ppp , g , y , tiny
      INTEGER N , ZBQLBIN , ZBQLGEO , iz , nn
 
      tiny = 1.0D-8
      ZBQLBIN = 0
 
      IF ( (P<0.0D0) .OR. (P>1.0D0) ) THEN
         WRITE (*,99001)
         RETURN
      ELSEIF ( N<=0 ) THEN
         WRITE (*,99001)
         RETURN
      ENDIF
!
!	First step: if NP > 10, say, things will be expensive, and
!	we can get into the right ballpark by guessing a value for
!	ZBQLBIN (IZ, say), and simulating Y from a Beta distribution
!	with parameters IZ and NN-IZ+1 (NN starts off equal to N).
!	If Y is less than PP (which starts off as P) then the IZth order
!	statistic from NN U(0,1) variates is less than P, and we know
!	that there are at least IZ successes. In this case we focus on
!	the remaining (NN-IZ) order statistics and count how many are
!	less than PP, which is binomial (NN-IZ,(PP-Y)/(1-Y)).
!	Otherwise, if Y is greater than PP there must be less
!	than IZ successes, so we can count the number of order statistics
!	under PP, which is binomial (IZ-1,P/Y). When we've got NN*PP
!	small enough, we go to the next stage of the algorithm and
!	generate the final bits directly.
!
      nn = N
      pp = P
      DO
         iz = INT(DBLE(nn)*pp) + 1
         IF ( (iz>10) .AND. (iz<nn-10) ) THEN
            y = ZBQLBET1(DBLE(iz),DBLE(nn-iz+1))
            IF ( y<pp ) THEN
               ZBQLBIN = ZBQLBIN + iz
               nn = nn - iz
               pp = (pp-y)/(1.0D0-y)
            ELSE
               nn = iz - 1
               pp = pp/y
            ENDIF
            CYCLE
         ENDIF
!
!	PP is the probability of the binomial we're currently
!	simulating from. For the final part, we simulate either number
!	of failures or number of success, depending which is cheaper.
!
         IF ( pp>0.5 ) THEN
            ppp = 1.0D0 - pp
         ELSE
            ppp = pp
         ENDIF
 
         g = 0
         iz = 0
!
!     ZBQLGEO falls over for miniscule values of PPP, so ignore these
!     (tiny probability of any successes in this case, anyway)
!
         IF ( ppp>tiny ) THEN
            DO
               g = g + ZBQLGEO(ppp)
               IF ( g<=nn ) THEN
                  iz = iz + 1
                  CYCLE
               ENDIF
               EXIT
            ENDDO
         ENDIF
 
         IF ( pp>0.5 ) iz = nn - iz
         ZBQLBIN = ZBQLBIN + iz
         EXIT
      ENDDO
 
99001 FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &        ' ZBQLBIN',/)
      END FUNCTION ZBQLBIN
!*==zbqlgeo.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLGEO(P)
      IMPLICIT NONE
!*--ZBQLGEO377
!
!       Returns a random number geometrically distributed with
!       parameter P ie. mean 1/P
!
 
      DOUBLE PRECISION P , ZBQLU01 , u , tiny
      INTEGER ZBQLGEO
 
      tiny = 1.0D-12
      ZBQLGEO = 0
 
      IF ( (P<0.0D0) .OR. (P>1.0D0) ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLGEO',/)
         RETURN
      ENDIF
 
      IF ( P>0.9D0 ) THEN
         DO
            ZBQLGEO = ZBQLGEO + 1
            u = ZBQLU01(0.0D0)
            IF ( u<=P ) EXIT
         ENDDO
      ELSE
         u = ZBQLU01(0.0D0)
!
!	For tiny P, 1-p will be stored inaccurately and log(1-p) may
!	be zero. In this case approximate log(1-p) by -p
!
         IF ( P>tiny ) THEN
            ZBQLGEO = 1 + INT(DLOG(u)/DLOG(1.0D0-P))
         ELSE
            ZBQLGEO = 1 + INT(-DLOG(u)/P)
         ENDIF
      ENDIF
      END FUNCTION ZBQLGEO
!*==zbqlpoi.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLPOI(Mu)
      IMPLICIT NONE
!*--ZBQLPOI420
!
!       Returns a random number Poisson distributed with mean MU
!
 
      DOUBLE PRECISION ZBQLU01 , x , y , Mu , pi
      DOUBLE PRECISION ZBQLLG , ZBQLGAM , mu1 , tmp1 , tmp2 , t
      INTEGER ZBQLPOI , ZBQLBIN , k , init
      SAVE init , pi
      DATA init/0/
 
      IF ( init==0 ) THEN
         pi = 4.0D0*DATAN(1.0D0)
         init = 1
      ENDIF
 
      ZBQLPOI = 0
 
      IF ( Mu<0.0D0 ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLPOI',/)
         RETURN
      ENDIF
!
!      For small MU, generate exponentials till their sum exceeds 1
!      (equivalently, uniforms till their product falls below e^-MU)
!
      IF ( Mu<=1.0D3 ) THEN
         mu1 = Mu
         DO
!
!     For values of MU less than 1000, use order statistics - the Kth
!     event in a Poisson process of rate MU has a Gamma distribution
!     with parameters (MU,K); if it's greater than 1 we know that there
!     are less than K events in (0,1) (and the exact number is binomial)
!     and otherwise the remaining number is another Poisson. Choose K so
!     that we'll get pretty close to 1 in the first go but are unlikely
!     to overshoot it.
!
            IF ( mu1>1.0D1 ) THEN
               k = INT(mu1-DSQRT(mu1))
               y = ZBQLGAM(DBLE(k),mu1)
               IF ( y>1.0D0 ) THEN
                  ZBQLPOI = ZBQLPOI + ZBQLBIN(k-1,(1.0D0/y))
                  RETURN
               ENDIF
               ZBQLPOI = ZBQLPOI + k
               mu1 = mu1*(1.0D0-y)
               CYCLE
            ENDIF
            y = DEXP(-mu1)
            x = 1.0D0
            EXIT
         ENDDO
         DO
            x = x*ZBQLU01(0.0D0)
            IF ( x>y ) THEN
               ZBQLPOI = ZBQLPOI + 1
               CYCLE
            ENDIF
            EXIT
         ENDDO
!
!     For really huge values of MU, use rejection sampling as in
!     Press et al (1992) - large numbers mean some accuracy may be
!     lost, but it caps the execution time.
!
      ELSE
         tmp1 = DSQRT(2.0D0*Mu)
         tmp2 = ZBQLLG(Mu+1.0D0) - (Mu*DLOG(Mu))
         DO
            y = DTAN(pi*ZBQLU01(0.0D0))
            ZBQLPOI = INT(Mu+(tmp1*y))
            IF ( ZBQLPOI>=0 ) THEN
               x = DBLE(ZBQLPOI)
               t = (x*DLOG(Mu)-ZBQLLG(x+1.0D0)) + tmp2
               IF ( DABS(t)<1.0D2 ) THEN
                  t = 0.9D0*(1.0D0+(y*y))*DEXP(t)
                  IF ( ZBQLU01(0.0D0)<=t ) EXIT
               ELSE
                  t = DLOG(0.9D0*(1.0D0+(y*y))) + t
                  IF ( DLOG(ZBQLU01(0.0D0))<=t ) EXIT
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      END FUNCTION ZBQLPOI
!*==zbqlgam.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*****************************************************************
      FUNCTION ZBQLGAM(G,H)
      IMPLICIT NONE
!*--ZBQLGAM513
!
!       Returns a random number with a gamma distribution with mean
!       G/H and variance G/(H^2). (ie. shape parameter G & scale
!       parameter H)
!
      DOUBLE PRECISION c , d , r , ZBQLGAM , ZBQLU01 , G , H , a , z1 , &
     &                 z2 , b1 , b2 , m
      DOUBLE PRECISION u1 , u2 , u , v , test , x
      DOUBLE PRECISION c1 , c2 , c3 , c4 , c5 , w
 
      ZBQLGAM = 0.0D0
 
      IF ( (G<=0.0D0) .OR. (H<0.0D0) ) THEN
         WRITE (*,99001)
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLGAM',/5X,'(both parameters must be positive)',/)
         RETURN
      ENDIF
 
      IF ( G<1.0D0 ) THEN
         DO
            u = ZBQLU01(0.0D0)
            v = ZBQLU01(0.0D0)
            IF ( u>EXP(1.0D0)/(G+EXP(1.0D0)) ) THEN
               ZBQLGAM = -LOG((G+EXP(1.0D0))*(1.0D0-u)/(G*EXP(1.0D0)))
               IF ( v<=ZBQLGAM**(G-1.0) ) EXIT
            ELSE
               ZBQLGAM = ((G+EXP(1.0D0))*u/EXP(1.0D0))**(1.0D0/G)
               IF ( v<=EXP(-ZBQLGAM) ) EXIT
            ENDIF
         ENDDO
         ZBQLGAM = ZBQLGAM/H
         RETURN
      ELSEIF ( G<2.0D0 ) THEN
         m = 0.0D0
      ELSEIF ( G>10.0D0 ) THEN
         c1 = G - 1.0D0
         c2 = (G-1.0D0/(6.0D0*G))/c1
         c3 = 2.0D0/c1
         c4 = c3 + 2.0D0
         c5 = 1.0D0/SQRT(G)
         DO
            u = ZBQLU01(0.0D0)
            v = ZBQLU01(0.0D0)
            IF ( G>2.50D0 ) u = v + c5*(1.0D0-1.860D0*u)
            IF ( u>0.0D0 .AND. u<1.0D0 ) THEN
               w = c2*v/u
               IF ( c3*u+w+1.0D0/w>c4 ) THEN
                  IF ( c3*LOG(u)-LOG(w)+w>=1.0D0 ) CYCLE
               ENDIF
               ZBQLGAM = c1*w/H
               RETURN
            ENDIF
         ENDDO
      ELSE
         m = -(G-2.0D0)
      ENDIF
      r = 0.50D0
      a = ((G-1.0D0)/EXP(1.0D0))**((G-1.0D0)/(r+1.0D0))
      c = (r*(m+G)+1.0D0)/(2.0D0*r)
      d = m*(r+1.0D0)/r
      z1 = c - DSQRT(c*c-d)
      z2 = c + DSQRT(c*c-d)
      b1 = (z1*(z1-m)**(r*(G-1.0D0)/(r+1.0D0)))&
     &     *DEXP(-r*(z1-m)/(r+1.0D0))
      b2 = (z2*(z2-m)**(r*(G-1.0D0)/(r+1.0D0)))&
     &     *DEXP(-r*(z2-m)/(r+1.0D0))
      DO
         u1 = ZBQLU01(0.0D0)
         u2 = ZBQLU01(0.0D0)
         u = a*u1
         v = b1 + (b2-b1)*u2
         x = v/(u**r)
         IF ( x>m ) THEN
            test = ((x-m)**((G-1)/(r+1)))*EXP(-(x-m)/(r+1.0D0))
            IF ( u>test ) CYCLE
            ZBQLGAM = (x-m)/H
            EXIT
         ENDIF
      ENDDO
 
      END FUNCTION ZBQLGAM
!*==zbqlbet1.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**************************************************************
      FUNCTION ZBQLBET1(Nu1,Nu2)
      IMPLICIT NONE
!*--ZBQLBET1600
!
!       Returns a random number, beta distributed with degrees
!       of freedom NU1 and NU2.
!
      DOUBLE PRECISION Nu1 , Nu2 , ZBQLGAM , ZBQLBET1 , ZBQLU01 , x1 , &
     &                 x2
 
      ZBQLBET1 = 0.0D0
 
      IF ( (Nu1<=0.0) .OR. (Nu2<=0.0) ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLBET1',/5X,&
     &           '(both degrees of freedom must be positive)',/)
         RETURN
      ENDIF
!
!       If parameters are too small, gamma subroutine tends to return ze
!       as all the probability goes to the origin and we get rounding
!       errors, even with double precision. In this case, we use Johnk's
!       method, suitably scaled to avoid rounding errors as much as poss
!
 
      IF ( (Nu1<0.9D0) .AND. (Nu2<0.9D0) ) THEN
         DO
            x1 = ZBQLU01(0.0D0)
            x2 = ZBQLU01(0.0D0)
            IF ( (x1**(1.0D0/Nu1))+(x2**(1.0D0/Nu2))<=1.0D0 ) THEN
               x1 = (DLOG(x2)/Nu2) - (DLOG(x1)/Nu1)
               ZBQLBET1 = (1.0D0+DEXP(x1))**(-1)
               IF ( ZBQLBET1<=1.0D0 ) EXIT
            ENDIF
         ENDDO
      ELSE
         x1 = ZBQLGAM(Nu1,1.0D0)
         x2 = ZBQLGAM(Nu2,1.0D0)
         ZBQLBET1 = x1/(x1+x2)
      ENDIF
 
      RETURN
 
      END FUNCTION ZBQLBET1
!*==zbqlwei.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**************************************************************
      FUNCTION ZBQLWEI(A,B)
      IMPLICIT NONE
!*--ZBQLWEI648
!
!       Returns a random number, Weibull distributed with shape paramete
!       A and location parameter B, i.e. density is
!	f(x) = ( A/(B**A) ) * x**(A-1) * EXP( -(x/B)**A )
!
      DOUBLE PRECISION A , B , ZBQLU01 , ZBQLWEI , u
 
      ZBQLWEI = 0.0D0
 
      IF ( (A<=0.0) .OR. (B<=0.0) ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLWEI',/5X,'(both parameters must be positive)',/)
         RETURN
      ENDIF
 
      u = ZBQLU01(0.0D0)
      ZBQLWEI = B*((-DLOG(u))**(1.0D0/A))
      END FUNCTION ZBQLWEI
!*==zbqlnb.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**************************************************************
      FUNCTION ZBQLNB(R,P)
      IMPLICIT NONE
!*--ZBQLNB673
!
!       Returns a pseudo-random number according to a Negative
!	Binomial distribution with parameters (R,P). NB these are
!	both DOUBLE - it copes with non-integer R as well. The
!       form of the distribution is *not* the no. of trials to
!       the Rth success - see documentation for full spec.
!
      DOUBLE PRECISION R , P , ZBQLGAM , y
      INTEGER ZBQLNB , ZBQLPOI
 
      ZBQLNB = 0
 
      IF ( (R<=0.0D0) .OR. (P<=0.0D0) .OR. (P>=1.0D0) ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLNB')
         RETURN
      ENDIF
 
      y = ZBQLGAM(R,1.0D0)
      y = y*P/(1.0D0-P)
      ZBQLNB = ZBQLPOI(y)
      END FUNCTION ZBQLNB
!*==zbqlpar.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**************************************************************
      FUNCTION ZBQLPAR(A,B)
      IMPLICIT NONE
!*--ZBQLPAR702
!
!     Returns a random number, Pareto distributed with parameters
!     A and B. The density is A*(B**A) / (B+X)**(A+1) for X > 0.
!     (this is slightly nonstandard - see documentation in
!     randgen.txt). The algorithm is straightforward - it uses the
!     inverse CDF method.
!
      DOUBLE PRECISION A , B , ZBQLPAR , ZBQLU01 , u
 
      ZBQLPAR = 0.0D0
 
      IF ( (A<=0.0D0) .OR. (B<=0.0D0) ) THEN
         WRITE (*,99001)
 
99001    FORMAT (/5X,'****ERROR**** Illegal parameter value in ',&
     &           ' ZBQLPAR',/5X,'(both parameters must be positive)',/)
         RETURN
      ENDIF
 
      u = ZBQLU01(0.0D0)
      ZBQLPAR = B*(u**(-1.0D0/A)-1.0D0)
      END FUNCTION ZBQLPAR
!*==zbqllg.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**************************************************************
      FUNCTION ZBQLLG(X)
      IMPLICIT NONE
!*--ZBQLLG729
!
!     Returns log(G(X)) where G is the Gamma function. The algorithm is
!     that given in Press et al (1992), Section 6.1, although this
!     version also allows for arguments less than 1.
!
      DOUBLE PRECISION X , z , z2 , ZBQLLG , pi , rln2p , c(0:6) , tmp ,&
     &                 sum
      INTEGER init , i
      SAVE init , c , rln2p , pi
      DATA init/0/
      DATA (c(i),i=0,6)/1.000000000190015D0 , 76.18009172947146D0 , &
     &      -86.50532032941677D0 , 24.01409824083091D0 , &
     &      -1.231739572450155D0 , 0.1208650973866179D-2 , &
     &      -0.5395239384953D-5/
 
      IF ( init==0 ) THEN
         pi = 4.0D0*DATAN(1.0D0)
         rln2p = 0.5D0*DLOG(2.0D0*pi)
         init = 1
      ENDIF
!
!     Compute for x > 1, then use transformation if necessary. Z is
!     our working argument.
!
      IF ( X>=1.0D0 ) THEN
         z = X
      ELSE
         z = 2.0D0 - X
         z2 = 1.0D0 - X
      ENDIF
 
      IF ( DABS(z-1.0D0)<1.0D-12 ) THEN
         ZBQLLG = 0.0D0
         RETURN
      ENDIF
 
      tmp = z + 4.5D0
      tmp = ((z-0.5D0)*DLOG(tmp)) - tmp + rln2p
 
      sum = c(0)
      DO i = 1 , 6
         sum = sum + (c(i)/(z+DBLE(i-1)))
      ENDDO
      ZBQLLG = tmp + DLOG(sum)
!
!     Transformation required if X<1
!
      IF ( X<1.0D0 ) THEN
         tmp = pi*z2
         ZBQLLG = DLOG(tmp/DSIN(tmp)) - ZBQLLG
      ENDIF
 
      END FUNCTION ZBQLLG
