!*==pmesh.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**    pmesh  :  control routine for the macro "feap".
!**
!**   Non-Obvious Parameters : NONE
!**
!--
      SUBROUTINE PMESH(Id,X,Ix,F,B,Dr,Itx, lmp)
      USE MOD_GLOBAL
      USE MOD_POTEN
      USE MOD_GRAIN
      USE MOD_MATERIAL
      USE MOD_BOUNDARY
      USE MOD_FILE
      use lammps
      IMPLICIT NONE
!*--PMESH15
!
      !---- data input routine for mesh description
      type (c_ptr) :: lmp
!
      INTEGER Id(NDF,*) , Ix(NEN1,*) , Itx(3,*)
      DOUBLE PRECISION B(NDF,*) , F(NDF,*) , X(NXDm,*) , Dr(NDF,*)
!
      LOGICAL PCOMP , init , die
      CHARACTER*4 va(2) , xm(3) , xv(3) , wd(25) , cd(3) , te(3) ,&
           & fd(3) , u0(3) , cc(2)
      CHARACTER*4 bl
      CHARACTER*80 geomfile , input
 
      INTEGER list , n , i , ii , lower , upper , NEXT , ldif
 
 
      DATA wd/'xxxx' , 'xxxx' , 'mate' , 'xxxx' , 'xxxx' , 'xxxx' ,&
           & 'end ' , 'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' ,&
           & 'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , 'mp01' , 'xxxx' ,&
           & 'xxxx' , 'xxxx' , 'cons' , 'grai' , 'xxxx' , 'xxxx' ,&
           & 'xxxx'/
      DATA bl/'blan'/ , list/25/
      DATA va/' val' , 'ue  '/
      DATA xm/' mas' , 'ses ' , '    '/ , xv/' vel' , 'ocit' , 'ies '&
           &/ , cd/' coo' , 'rdin' , 'ates'/ , te/' tem' , 'pera' ,&
           & 'ture'/ , fd/' for' , 'ce/d' , 'ispl'/ , u0/' dis' ,&
           & 'pl. ' , '    '/
!
!
!---- initialize arrays
      init = .TRUE.
      DO n = 1 , NUMnp
         DO i = 1 , NXDm
            IF ( X(i,n)/=0. ) init = .FALSE.
         ENDDO
      ENDDO
      IF ( init ) THEN
         DO n = 1 , NUMnp
            X(1,n) = 0.0
         ENDDO
      ENDIF
 100  DO
!---- read macro cards
         READ (input_file_unit,'(a80)') input
         upper = 0
         DO ii = 1 , 2
            lower = upper
            upper = NEXT(lower,input)
            ldif = upper - lower - 1
            IF ( ldif==0 ) THEN
               cc(ii) = '    '
            ELSE
               CALL KEYEQUAL(cc(ii),input,lower+1,upper-1)
            ENDIF
         ENDDO
         DO i = 1 , list
            IF ( PCOMP(cc(1),wd(i)) ) GOTO 200
         ENDDO
      ENDDO
!---- process macros
 200  SELECT CASE (i)
      CASE (1)
         GOTO 100
      CASE (2)
         GOTO 100
      CASE (3)
!
!---- macro 'mate'
!---- material data input
!--Read Material Property table (Ellad's routine)
         WRITE (*,*) '** Reading Material Data'
         WRITE (*,*)
         CALL READMATERIALS(NDF,cc(2))
         CALL OUTPUTMATERIALS()
!
!     Verify user input is correct for this element
!
         IF ( NDM/=2 ) THEN
            PRINT * , '***ERROR: Incorrect spatial dimension (ndm=2)'
            STOP
         ENDIF
         IF ( NEN/=3 ) THEN
            PRINT * , '***ERROR: Incorrect #nodes per element (nen=3)'
            STOP
         ENDIF
         IF ( NSDm/=6 ) THEN
            PRINT * , '***ERROR: Incorrect stress dimension (nsdm=6)'
            STOP
         ENDIF
         IF ( NQUad/=1 ) THEN
            PRINT * , '***ERROR: Incorrect #quadrature points (nquad=1)'
            STOP
         ENDIF
         IF ( NAD/=0 ) THEN
            PRINT * , '***ERROR: Incorrect #internal nodes (nad=0)'
            STOP
         ENDIF
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!---- define the model boundaries
!
!---- macro 'xxxx'
         GOTO 100
      CASE (4)
         GOTO 100
      CASE (5)
         GOTO 100
      CASE (6)
         GOTO 100
      CASE (7)
!
!---- macro 'end '
!---- terminate mesh input
!
! check that some basic stuff is properly defined
!
         die = .FALSE.
         IF ( NMAterials==0 ) THEN
            WRITE (*,*) '**ERROR: no materials are defined'
            die = .TRUE.
         ENDIF
         IF ( NGRains==0 ) THEN
            WRITE (*,*) '**ERROR: no grains are defined'
            die = .TRUE.
         ENDIF
         DO i = 1 , 103
            IF ( MAPspecies(i)/=0 ) GOTO 250
         ENDDO
         WRITE (*,*) '**ERROR: no constitutive information found'
         die = .TRUE.
 250     IF ( die ) STOP
         WRITE (6,'(/2x,''**end of mesh definition**''/)')
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
         RETURN
      CASE (8)
         GOTO 100
      CASE (9)
         GOTO 100
      CASE (10)
         GOTO 100
      CASE (11)
         GOTO 100
      CASE (12)
         GOTO 100
      CASE (13)
         GOTO 100
      CASE (14)
         GOTO 100
      CASE (15)
         GOTO 100
      CASE (16)
         GOTO 100
      CASE (17)
!
!---- macro 'mp01'
!---- user supplied macro
         CALL MP01(Id,X,Ix,F,B,Itx, lmp)
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
         GOTO 100
      CASE (18)
         GOTO 100
      CASE (19)
         GOTO 100
      CASE (20)
         GOTO 100
      CASE (21)
!
!---- macro 'cons'
!---- read in potential specific constitutive input
         CALL READCONSTITUTIVE()
         GOTO 100
      CASE (22)
!
!---- macro 'grai'
!---- user supplied macro
!--Read in filename
         IF ( cc(2)=='dire' ) THEN
            geomfile = ' '
         ELSE
            READ (input_file_unit,99001) geomfile
99001       FORMAT (a)
         ENDIF
 
!--Read data form  file
         CALL READGRAINDATA(geomfile,cc(2))
!
         CALL PROCESSGRAINS
!--Print out grain structure
         CALL OUTPUTGRAINDATA()
         IF ( NGRains>1 ) STOP 'hardwired for 1 grain'
         CALL ROTATEBURGERS(GRAins(1)%XLATVECT(1:3,1:3),GRAins(1)&
     &                      %ROTMAT(1:2),MATerial(GRAins(1)%MATGRAIN)&
     &                      %A0,MATerial(GRAins(1)%MATGRAIN)%STRUCTURE)
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
         GOTO 100
!
99002    FORMAT (/5x,'material set',i3,' for element type',i2)
99003    FORMAT (/' ',a80//5x,'material properties')
99004    FORMAT (i10,9E13.3)
      CASE (23)
         GOTO 100
      CASE (24)
         GOTO 100
      CASE DEFAULT
         GOTO 100
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
      END SELECT
      END SUBROUTINE PMESH
 
 
