!*==contri.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
      SUBROUTINE CONTRI(Npts,N,Nce,Ncb,Elist,Xx,Nxdm,List,W,V,T,Ntri)
      IMPLICIT NONE
!*--CONTRI4
!*** Start of declarations inserted by SPAG
      INTEGER Nxdm
!*** End of declarations inserted by SPAG
!**********************************************************************
!
!     PURPOSE:
!     --------
!
!     Assemble constrained Delaunay triangulation for collection of
!     points in the plane. This code has a total memory requirement
!     equal to 4*NPTS + 13*N + 2*NCE + 18
!
!     INPUT:
!     ------
!
!     NPTS   - Total number of points in data set (NPTS ge N)
!     N      - Total number of points to be triangulated (N ge 3)
!     NCE    - Total number of constrained edges which must be
!              present, including those which define boundaries
!            - NCE=0 indicates triangulation is unconstrained so that
!              ELIST is empty and the code will produce a
!              triangulation of the convex hull
!     NCB    - Total number of constrained edges which define one
!              external boundary and any internal boundaries (holes)
!            - NCB=0 indicates there are no boundary edge constraints
!              and the code will produce a triangulation of the convex
!              hull
!            - NCB must not be greater than NCE
!            - If NCB gt 0, then at least one of the boundaries
!              specified must be external
!     ELIST  - List of edges which must be present in triangulation
!            - These may be internal edges or edges which define
!              boundaries
!            - Edge I defined by vertices in ELIST(1,I) and ELIST(2,I)
!            - Edges which define boundaries must come first in ELIST
!              and thus occupy the first NCB columns
!            - Edges which define an external boundary must be listed
!              anticlockwise but may be presented in any order
!            - Edges which define an internal boundary (hole) must be
!              listed clockwise but may be presented in any order
!            - An internal boundary (hole) cannot be specified unless
!              an external boundary is also specified
!            - All boundaries must form closed loops
!            - An edge may not appear more than once in ELIST
!            - An external or internal boundary may not cross itself
!              and may not share a common edge with any other boundary
!            - Internal edges, which are not meant to define boundaries
!              but must be present in the final triangulation, occupy
!              columns NCB+1,... ,NCE of ELIST
!     X      - X-coords of all points in data set
!            - X-coord of point I given by X(I)
!            - Last three locations are used to store x-coords of
!              supertriangle vertices in subroutine delaun
!     Y      - Y-coords of all points in data set
!            - Y-coord of point I given by Y(I)
!            - Last three locations are used to store y-coords of
!              supertriangle vertices in subroutine delaun
!     LIST   - List of points to be triangulated
!            - If N eq NPTS, set LIST(I)=I for I=1,2,...,NPTS
!              prior to calling this routine
!            - No point in LIST may lie outside any external boundary
!              or inside any internal boundary
!     W      - Not defined, workspace of length ge 2*(NPTS+3)
!     V      - Not defined
!     T      - Not defined
!     NTRI   - Not defined
!
!     OUTPUT:
!     -------
!
!     NPTS   - Unchanged
!     N      - Unchanged
!     NCE    - Unchanged
!     NCB    - Unchanged
!     ELIST  - Unchanged
!     X      - Unchanged, except that last three locations contain
!              normalised x-coords of supertriangle vertices
!     Y      - Unchanged, except that last three locations contain
!              normalised y-coords of supertriangle vertices
!     LIST   - List of points to be triangulated
!            - Points ordered such that consecutive points are close
!              to one another in the x-y plane
!     W      - Not defined
!     V      - Vertex array for triangulation
!            - Vertices listed in anticlockwise sequence
!            - Vertices for triangle J are found in V(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - First vertex is at point of contact of first and third
!              adjacent triangles
!     T      - Adjacency array for triangulation
!            - Triangles adjacent to J are found in T(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - Adjacent triangles listed in anticlockwise sequence
!            - Zero denotes no adjacent triangle
!     NTRI   - Total number of triangles in final triangulation
!
!     SUBROUTINES  CALLED:  BSORT, DELAUN, EDGE, TCHECK, BCLIP
!     --------------------
!
!     PROGRAMMER:             Scott Sloan
!     -----------
!
!     LAST MODIFIED:          3 march 1991        Scott Sloan
!     --------------
!
!     COPYRIGHT 1990:         Scott Sloan
!     ---------------         Department of Civil Engineering
!                             University of Newcastle
!                             NSW 2308
!
!**********************************************************************
      INTEGER i , j , N , p
      INTEGER jw , nb , vi , vj
      INTEGER Ncb , Nce , nef
      INTEGER Npts , Ntri
      INTEGER T(3,2*N+1) , V(3,2*N+1) , W(Npts+3,2)
      INTEGER List(N)
      INTEGER Elist(2,*)
!
      DOUBLE PRECISION dmax , xmin , xmax , ymin , ymax
      DOUBLE PRECISION C00001
      DOUBLE PRECISION fact
      DOUBLE PRECISION Xx(Nxdm,*)
      DOUBLE PRECISION , POINTER :: x(:) , y(:)
!
      PARAMETER (C00001=1.0)
!
! Ron's addition: just putting xx(nxdm,1) into X and Y
!
      ALLOCATE (x(Npts+3),y(Npts+3))
      DO i = 1 , Npts
         x(i) = Xx(1,i)
         y(i) = Xx(2,i)
      ENDDO
!---------------------------------------------------------------------
!     Check input for obvious errors
!
      IF ( Npts<3 ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''LESS THAN 3 POINTS IN DATASET'',                           &
     &   /,''CHECK VALUE OF NPTS'')')
         STOP
      ENDIF
      IF ( N<3 ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''LESS THAN 3 POINTS TO BE TRIANGULATED'',                   &
     &   /,''CHECK VALUE OF N'')')
         STOP
      ELSEIF ( N>Npts ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''TOO MANY POINTS TO BE TRIANGULATED'',                      &
     &   /,''N MUST BE LESS THAN OR EQUAL TO NPTS'',                    &
     &   /,''CHECK VALUES OF N AND NPTS'')')
         STOP
      ENDIF
      IF ( Ncb>Nce ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''NCB GREATER THAN NCE'',                                    &
     &   /,''CHECK BOTH VALUES'')')
         STOP
      ENDIF
!---------------------------------------------------------------------
!     Check for invalid entries in LIST
!
      DO i = 1 , N
         p = List(i)
         IF ( p<1 .OR. p>Npts ) THEN
            WRITE (6,99001) i , p
!---------------------------------------------------------------------
99001       FORMAT (//,'***ERROR IN CONTRI***',/,&
     &              'ILLEGAL VALUE IN LIST',/,'LIST(',I5,' )=',I5)
            STOP
         ENDIF
         W(p,1) = 0
      ENDDO
      DO i = 1 , N
         p = List(i)
         W(p,1) = W(p,1) + 1
      ENDDO
      DO i = 1 , N
         p = List(i)
         IF ( W(p,1)>1 ) THEN
            WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''POINT'',I5,'' OCCURS'',I5,'' TIMES IN LIST'',              &
     &   /,''POINT SHOULD APPEAR ONLY ONCE'')') p , W(p,1)
            STOP
         ENDIF
      ENDDO
!---------------------------------------------------------------------
!     Check for invalid entries in ELIST
!
      IF ( Nce>0 ) THEN
         DO i = 1 , Npts
            W(i,1) = 0
         ENDDO
         DO i = 1 , N
            W(List(i),1) = 1
         ENDDO
         DO i = 1 , Nce
            vi = Elist(1,i)
            IF ( vi<1 .OR. vi>Npts ) THEN
               WRITE (6,99002) 1 , i , vi
               STOP
            ELSEIF ( W(vi,1)/=1 ) THEN
               WRITE (6,99003) 1 , i , vi
               STOP
            ENDIF
            vj = Elist(2,i)
            IF ( vj<1 .OR. vj>Npts ) THEN
               WRITE (6,99002) 2 , i , vj
               STOP
            ELSEIF ( W(vj,1)/=1 ) THEN
               WRITE (6,99003) 2 , i , vj
               STOP
            ENDIF
         ENDDO
      ENDIF
!---------------------------------------------------------------------
!     Check that all boundaries (if there are any) form closed loops
!     Count appearances in ELIST(1,.) and ELIST(2,.) of each node
!     These must match if each boundary forms a closed loop
!
      IF ( Ncb>0 ) THEN
         DO i = 1 , Npts
            W(i,1) = 0
            W(i,2) = 0
         ENDDO
         DO i = 1 , Ncb
            vi = Elist(1,i)
            vj = Elist(2,i)
            W(vi,1) = W(vi,1) + 1
            W(vj,2) = W(vj,2) + 1
         ENDDO
         DO i = 1 , Ncb
            vi = Elist(1,i)
            IF ( W(vi,1)/=W(vi,2) ) THEN
               WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''BOUNDARY SEGMENTS DO NOT FORM A '',                        &
     &     ''CLOSED LOOP'',                                             &
     &   /,''CHECK ENTRIES IN COLS 1...NCB OF ELIST '',                 &
     &     ''FOR NODE'',I5)') vi
               STOP
            ENDIF
         ENDDO
      ENDIF
!---------------------------------------------------------------------
!     Compute min and max coords for x and y
!     Compute max overall dimension
!
      p = List(1)
      xmin = x(p)
      xmax = xmin
      ymin = y(p)
      ymax = ymin
      DO i = 2 , N
         p = List(i)
         xmin = MIN(xmin,x(p))
         xmax = MAX(xmax,x(p))
         ymin = MIN(ymin,y(p))
         ymax = MAX(ymax,y(p))
      ENDDO
      IF ( xmin==xmax ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''ALL POINTS HAVE SAME X-COORD'',                            &
     &   /,''ALL POINTS ARE COLLINEAR'')')
         STOP
      ENDIF
      IF ( ymin==ymax ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN CONTRI***'',                                   &
     &   /,''ALL POINTS HAVE SAME Y-COORD'',                            &
     &   /,''ALL POINTS ARE COLLINEAR'')')
         STOP
      ENDIF
      dmax = MAX(xmax-xmin,ymax-ymin)
!---------------------------------------------------------------------
!     Normalise x-y coords of points
!
      fact = C00001/dmax
      DO i = 1 , N
         p = List(i)
         x(p) = (x(p)-xmin)*fact
         y(p) = (y(p)-ymin)*fact
      ENDDO
!---------------------------------------------------------------------
!     Sort points into bins using a linear bin sort
!     This call is optional
!
      CALL BSORT(N,Npts,x,y,xmin,xmax,ymin,ymax,dmax,W,List,W(1,2))
!---------------------------------------------------------------------
!     Compute Delaunay triangulation
!
      CALL DELAUN(Npts,N,x,y,List,W,V,T,Ntri)
!---------------------------------------------------------------------
!     For each node, store any triangle in which it is a vertex
!     Include supertriangle vertices
!
      DO j = 1 , Ntri
         DO i = 1 , 3
            vi = V(i,j)
            W(vi,1) = j
         ENDDO
      ENDDO
!---------------------------------------------------------------------
!     Constrain triangulation by forcing edges to be present
!
      nef = 0
      jw = (Npts+3)/2
      DO i = 1 , Nce
         vi = Elist(1,i)
         vj = Elist(2,i)
         CALL EDGE(vi,vj,Npts,Ntri,nef,jw,x,y,V,T,W,W(1,2))
      ENDDO
!---------------------------------------------------------------------
!     Clip triangulation and check it
!
      CALL BCLIP(Npts,Ncb,Elist,W,V,T,Ntri,nb)
      CALL TCHECK(Npts,N,x,y,List,V,T,Ntri,nef,nb,Nce,Ncb,Elist,W)
!---------------------------------------------------------------------
!     Reset x-y coords to original values
!
      DO i = 1 , N
         p = List(i)
         x(p) = x(p)*dmax + xmin
         y(p) = y(p)*dmax + ymin
      ENDDO
      DEALLOCATE (x,y)
99002 FORMAT (//,'***ERROR IN CONTRI***',/,'ILLEGAL VALUE IN ELIST',/,&
     &        'ELIST(',I5,',',I5,' )=',I5)
99003 FORMAT (//,'***ERROR IN CONTRI***',/,'ILLEGAL VALUE IN ELIST',/,&
     &        'ELIST(',I5,',',I5,' )=',I5,/,'THIS POINT IS NOT IN LIST')
      END SUBROUTINE CONTRI
!*==bclip.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE BCLIP(Npts,Ncb,Blist,Tn,V,T,Ntri,Nb)
      IMPLICIT NONE
!*--BCLIP348
!**********************************************************************
!
!     PURPOSE:
!     --------
!
!     Clip constrained Delaunay triangulation to boundaries in BLIST
!     If BLIST is empty, then the triangulation is clipped to a convex
!     hull by removing all triangles that are formed with the
!     supertriangle vertices
!     The triangulation MUST initially use the supertriangle
!     vertices
!
!     INPUT:
!     ------
!
!     NPTS   - Total number of points in data set (NPTS ge N)
!     NCB    - Total number of constrained edges which define one
!              external boundary and any internal boundaries (holes)
!            - NCB=0 indicates there are no boundary edge constraints
!              and the code will produce a triangulation of the convex
!              hull
!     BLIST  - List of edges which must be present in triangulation
!              and define boundaries
!            - Edge I defined by vertices in BLIST(1,I) and BLIST(2,I)
!              where I ranges from 1,...,NCB
!            - Edges which define an external boundary must be listed
!              anticlockwise but may be presented in any order
!            - Edges which define an internal boundary (hole) must be
!              listed clockwise but may be presented in any order
!     TN     - List of triangle numbers such that vertex I can be
!              found in triangle TN(I)
!     V      - Vertex array for triangulation
!            - Vertices listed in anticlockwise sequence
!            - Vertices for triangle J are found in V(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - First vertex is at point of contact of first and third
!              adjacent triangles
!     T      - Adjacency array for triangulation
!            - Triangles adjacent to J are found in T(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - Adjacent triangles listed in anticlockwise sequence
!            - Zero denotes no adjacent triangle
!     NTRI   - Number of triangles, including those formed with
!              vertices of supertriangle
!     NB     - Not defined
!
!     OUTPUT:
!     -------
!
!     NPTS   - Unchanged
!     NCB    - Unchanged
!     BLIST  - Unchanged
!     TN     - Not defined
!     V      - Updated vertex array for triangulation
!            - Vertices listed in anticlockwise sequence
!            - Vertices for triangle J are found in V(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - First vertex is at point of contact of first and third
!              adjacent triangles
!     T      - Updated adjacency array for triangulation
!            - Triangles adjacent to J are found in T(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - Adjacent triangles listed in anticlockwise sequence
!            - Zero denotes no adjacent triangle
!     NTRI   - Updated number of triangles in final triangulation
!     NB     - Number of boundaries defining the triangulation
!            - NB=1 for a simply connected domain with no holes
!            - NB=H+1 for a domain with H holes
!
!     PROGRAMMER:           Scott Sloan
!     -----------
!
!     LAST MODIFIED:        3 march 1991        Scott Sloan
!     --------------
!
!     COPYRIGHT 1990:       Scott Sloan
!     ---------------       Department of Civil Engineering
!                           University of Newcastle
!                           NSW 2308
!
!**********************************************************************
      INTEGER a , c , e , i , j , l , r , s
      INTEGER Nb , ts , v1 , v2
      INTEGER Ncb , nev
      INTEGER Npts , Ntri
      INTEGER nntri , npts1 , ntri1
      INTEGER T(3,Ntri) , V(3,Ntri)
      INTEGER Tn(Npts+3)
      INTEGER Blist(2,*)
!---------------------------------------------------------------------
!     Skip to triangle deletion step if triangulation has no
!     boundary constraints
!
      IF ( Ncb==0 ) THEN
         Nb = 1
         GOTO 100
      ENDIF
!---------------------------------------------------------------------
!     Mark all edges which define the boundaries
!     S=triangle in which search starts
!     R=triangle to right of edge V1-V2
!
      DO i = 1 , Ncb
         v1 = Blist(1,i)
         v2 = Blist(2,i)
         s = Tn(v1)
         r = s
         DO
!
!       Circle anticlockwise round V1 until edge V1-V2 is found
!
            IF ( V(1,r)==v1 ) THEN
               IF ( V(3,r)==v2 ) THEN
                  T(3,r) = -T(3,r)
                  EXIT
               ELSE
                  r = ABS(T(3,r))
               ENDIF
            ELSEIF ( V(2,r)==v1 ) THEN
               IF ( V(1,r)==v2 ) THEN
                  T(1,r) = -T(1,r)
                  EXIT
               ELSE
                  r = ABS(T(1,r))
               ENDIF
            ELSEIF ( V(2,r)==v2 ) THEN
               T(2,r) = -T(2,r)
               EXIT
            ELSE
               r = ABS(T(2,r))
            ENDIF
            IF ( r==s ) THEN
               WRITE (6,&
     &'(//,''***ERROR IN BCLIP***'',                                    &
     &   /,''CONSTRAINED BOUNDARY EDGE NOT FOUND'',                     &
     &   /,''V1='',I5,'' V2='',I5,                                      &
     &   /,''CHECK THAT THIS EDGE IS NOT CROSSED'',                     &
     &   /,''BY ANOTHER BOUNDARY EDGE'')') v1 , v2
               STOP
            ENDIF
         ENDDO
      ENDDO
!--------------------------------------------------------------------
!     Mark all triangles which are to right of boundaries by
!     circling anticlockwise around the outside of each boundary
!     Loop while some boundary edges have not been visited
!     S = triangle from which search starts
!     NEV = number of edges visited
!
      s = 0
      Nb = 0
      nev = 0
      ntri1 = Ntri + 1
      npts1 = Npts + 1
      DO WHILE ( nev<Ncb )
         DO
!
!       Find an edge on a boundary
!
            s = s + 1
            IF ( T(1,s)<0 ) THEN
               e = 1
            ELSEIF ( T(2,s)<0 ) THEN
               e = 2
            ELSEIF ( T(3,s)<0 ) THEN
               e = 3
            ELSE
               CYCLE
            ENDIF
!
!       Store and mark starting edge
!       Mark starting triangle for deletion
!       Increment count of edges visited
!       C = current triangle
!
            ts = T(e,s)
            T(e,s) = ntri1
            V(1,s) = npts1
            nev = nev + 1
            c = s
!
!       Increment to next edge
!
            e = MOD(e+1,3) + 1
            EXIT
         ENDDO
         DO
!
!       Loop until trace of current boundary is complete
!
            IF ( T(e,c)/=ntri1 ) THEN
               IF ( T(e,c)<0 ) THEN
!
!           Have found next boundary edge
!           Increment count of boundary edges visited, unmark the edge
!           and move to next edge
!
                  nev = nev + 1
                  T(e,c) = -T(e,c)
                  e = MOD(e+1,3) + 1
               ELSE
!
!           Non-boundary edge, circle anticlockwise around boundary
!           vertex to move to next triangle, and mark next
!           triangle for deletion
!
                  l = c
                  c = T(e,l)
                  IF ( T(1,c)==l ) THEN
                     e = 3
                  ELSEIF ( T(2,c)==l ) THEN
                     e = 1
                  ELSE
                     e = 2
                  ENDIF
                  V(1,c) = npts1
               ENDIF
               CYCLE
            ENDIF
!
!       Trace of current boundary is complete
!       Reset marked edge to original value and check for any more
!       boundaries
!
            T(e,c) = -ts
            Nb = Nb + 1
            EXIT
         ENDDO
      ENDDO
!---------------------------------------------------------------------
!     Remove all triangles which have been marked for deletion
!     Any triangle with a vertex number greater than NPTS is deleted
!
 100  nntri = Ntri
      Ntri = 0
      DO j = 1 , nntri
         IF ( MAX(V(1,j),V(2,j),V(3,j))>Npts ) THEN
!
!         Triangle J is to be deleted
!         Update adjacency lists for triangles adjacent to J
!
            DO i = 1 , 3
               a = T(i,j)
               IF ( a/=0 ) THEN
                  IF ( T(1,a)==j ) THEN
                     T(1,a) = 0
                  ELSEIF ( T(2,a)==j ) THEN
                     T(2,a) = 0
                  ELSE
                     T(3,a) = 0
                  ENDIF
               ENDIF
            ENDDO
         ELSE
!
!         Triangle J is not to be deleted
!         Update count of triangles
!
            Ntri = Ntri + 1
            IF ( Ntri<j ) THEN
!
!           At least one triangle has already been deleted
!           Relabel triangle J as triangle NTRI
!
               DO i = 1 , 3
                  a = T(i,j)
                  T(i,Ntri) = a
                  V(i,Ntri) = V(i,j)
                  IF ( a/=0 ) THEN
                     IF ( T(1,a)==j ) THEN
                        T(1,a) = Ntri
                     ELSEIF ( T(2,a)==j ) THEN
                        T(2,a) = Ntri
                     ELSE
                        T(3,a) = Ntri
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      END SUBROUTINE BCLIP
!*==bsort.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE BSORT(N,Npts,X,Y,Xmin,Xmax,Ymin,Ymax,Dmax,Bin,List,W)
      IMPLICIT NONE
!*--BSORT635
!***********************************************************************
!
!     PURPOSE:
!     --------
!
!     Sort points such that consecutive points are close to one another
!     in the x-y plane using a bin sort
!
!     INPUT:
!     ------
!
!     N      - Total number of points to be triangulated (N le NPTS)
!     NPTS   - Total number of points in data set
!     X      - X-coords of all points in data set
!            - If point is in list,the coordinate must be normalised
!              according to X=(X-XMIN)/DMAX
!            - X-coord of point I given by X(I)
!            - Last three locations are used to store x-coords of
!              supertriangle vertices in subroutine delaun
!     Y      - Y-coords of all points in data set
!            - If point is in list,the coordinate must be normalised
!              according to Y=(Y-YMIN)/DMAX
!            - Y-coord of point I given by Y(I)
!            - Last three locations are used to store y-coords of
!              supertriangle vertices in subroutine delaun
!     XMIN   - Min x-coord of points in LIST
!     XMAX   - Max x-coord of points in LIST
!     YMIN   - Min y-coord of points in LIST
!     YMAX   - Max y-coord of points in LIST
!     DMAX   - DMAX=MAX(XMAX-XMIN,YMAX-YMIN)
!     BIN    - Not defined
!     LIST   - List of points to be triangulated
!     W      - Undefined workspace
!
!     OUTPUT:
!     -------
!
!     N      - Unchanged
!     NPTS   - Unchanged
!     X      - Unchanged
!     Y      - Unchanged
!     XMIN   - Unchanged
!     XMAX   - Unchanged
!     YMIN   - Unchanged
!     YMAX   - Unchanged
!     DMAX   - Unchanged
!     BIN    - Not used
!     LIST   - List of points to be triangulated
!            - Points ordered such that consecutive points are close
!              to one another in the x-y plane
!     W      - Not used
!
!     PROGRAMMER:             Scott Sloan
!     -----------
!
!     LAST MODIFIED:          3 march 1991        Scott Sloan
!     --------------
!
!     COPYRIGHT 1990:         Scott Sloan
!     ---------------         Department of Civil Engineering
!                             University of Newcastle
!                             NSW 2308
!
!***********************************************************************
      INTEGER b , i , j , k , N , p
      INTEGER li , lt , nb
      INTEGER ndiv , Npts
      INTEGER W(N)
      INTEGER Bin(Npts)
      INTEGER List(N)
!
      DOUBLE PRECISION Dmax , Xmax , Xmin , Ymax , Ymin
      DOUBLE PRECISION factx , facty
      DOUBLE PRECISION X(Npts+3) , Y(Npts+3)
!---------------------------------------------------------------------
!     Compute number of bins in x-y directions
!     Compute inverse of bin size in x-y directions
!
      ndiv = NINT(REAL(N)**0.25)
      factx = REAL(ndiv)/((Xmax-Xmin)*1.01/Dmax)
      facty = REAL(ndiv)/((Ymax-Ymin)*1.01/Dmax)
!---------------------------------------------------------------------
!     Zero count of points in each bin
!
      nb = ndiv*ndiv
      DO i = 1 , nb
         W(i) = 0
      ENDDO
!---------------------------------------------------------------------
!     Assign bin numbers to each point
!     Count entries in each bin and store in W
!
      DO k = 1 , N
         p = List(k)
         i = INT(Y(p)*facty)
         j = INT(X(p)*factx)
         IF ( MOD(i,2)==0 ) THEN
            b = i*ndiv + j + 1
         ELSE
            b = (i+1)*ndiv - j
         ENDIF
         Bin(p) = b
         W(b) = W(b) + 1
      ENDDO
!---------------------------------------------------------------------
!     Form pointers to end of each bin in final sorted list
!     Note that some bins may contain no entries
!
      DO i = 2 , nb
         W(i) = W(i-1) + W(i)
      ENDDO
!---------------------------------------------------------------------
!     Now perform linear sort
!
      DO i = 1 , N
         IF ( List(i)>0 ) THEN
            li = List(i)
            b = Bin(li)
            p = W(b)
            DO
               IF ( p/=i ) THEN
                  W(b) = W(b) - 1
                  lt = List(p)
                  List(p) = li
                  li = lt
                  b = Bin(li)
                  List(p) = -List(p)
                  p = W(b)
                  CYCLE
               ENDIF
               W(b) = W(b) - 1
               List(i) = li
               EXIT
            ENDDO
         ELSE
            List(i) = -List(i)
         ENDIF
      ENDDO
      END SUBROUTINE BSORT
!*==delaun.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE DELAUN(Npts,N,X,Y,List,Stack,V,T,Ntri)
      IMPLICIT NONE
!*--DELAUN779
!***********************************************************************
!
!     PURPOSE:
!     --------
!
!     Assemble Delaunay triangulation
!
!     INPUT:
!     ------
!
!     NPTS   - Total number of points in data set
!     N      - Total number of points to be triangulated (N le NPTS)
!     X      - X-coords of all points in data set
!            - X-coord of point I given by X(I)
!            - If point is in LIST, coordinate must be normalised
!              such that X=(X-XMIN)/DMAX
!            - Last three locations are used to store x-coords of
!              supertriangle vertices in subroutine delaun
!     Y      - Y-coords of all points in data set
!            - Y-coord of point I given by Y(I)
!            - If point is in LIST, coordinate must be normalised
!              such that Y=(Y-YMIN)/DMAX
!            - Last three locations are used to store y-coords of
!              supertriangle vertices in subroutine delaun
!     LIST   - List of points to be triangulated
!            - Coincident points are flagged by an error message
!            - Points are ordered such that consecutive points are
!              close to one another in the x-y plane
!     STACK  - Not defined
!     V      - Not defined
!     T      - Not defined
!     NTRI   - Not defined
!
!     OUTPUT:
!     -------
!
!     NPTS   - Unchanged
!     N      - Unchanged
!     X      - Unchanged
!     Y      - Unchanged
!     LIST   - Unchanged
!     STACK  - Not defined
!     V      - Vertex array for triangulation
!            - Vertices listed in anticlockwise sequence
!            - Vertices for triangle J are found in V(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - First vertex is at point of contact of first and third
!              adjacent triangles
!     T      - Adjacency array for triangulation
!            - Triangles adjacent to J are found in T(I,J) for I=1,2,3
!              J=1,2,...,NTRI
!            - Adjacent triangles listed in anticlockwise sequence
!            - Zero denotes no adjacent triangle
!     NTRI   - Number of triangles in triangulation, including those
!              formed with the supertriangle vertices
!            - NTRI = 2*N+1
!
!     NOTES:
!     ------
!
!     - This is a faster version of the code appearing in AES 1987 vol 9
!       (small subroutines and functions have been coded in-line)
!     - Also some changes in code to improve efficiency
!     - Saving in cpu-time about 60 percent for larger problems
!     - A test has been implemented to detect coincident points
!     - Triangles formed with supertriangle vertices have not been
!       deleted
!
!     PROGRAMMER:             Scott Sloan
!     -----------
!
!     LAST MODIFIED:          3 march 1991          Scott Sloan
!     --------------
!
!     COPYRIGHT 1990:         Scott Sloan
!     ---------------         Department of Civil Engineering
!                             University of Newcastle
!                             NSW 2308
!
!***********************************************************************
      INTEGER a , b , c , i , j , l , N , p , r
      INTEGER v1 , v2 , v3
      INTEGER Npts , Ntri
      INTEGER topstk
      INTEGER T(3,2*N+1) , V(3,2*N+1)
      INTEGER List(N)
      INTEGER Stack(Npts)
!
      DOUBLE PRECISION d
      DOUBLE PRECISION xp , yp
      DOUBLE PRECISION TOL , x13 , x23 , x1p , x2p , y13 , y23 , y1p , &
     &                 y2p
      DOUBLE PRECISION cosa , cosb
      DOUBLE PRECISION C00000 , C00100
      DOUBLE PRECISION X(Npts+3) , Y(Npts+3)
!
!     TOL is the tolerance used to detect coincident points
!     The square of the distance between any two points must be
!     greater then TOL to avoid an error message
!     This value of TOL is suitable for single precision on most
!     32-bit machines (which typically have a precision of 0.000001)
!
      PARAMETER (TOL=1.E-10)
      PARAMETER (C00000=0.0,C00100=100.0)
!---------------------------------------------------------------------
!     Define vertices for supertriangle
!
      v1 = Npts + 1
      v2 = Npts + 2
      v3 = Npts + 3
!---------------------------------------------------------------------
!     Set coords of supertriangle
!
      X(v1) = -C00100
      X(v2) = C00100
      X(v3) = C00000
      Y(v1) = -C00100
      Y(v2) = -C00100
      Y(v3) = C00100
!---------------------------------------------------------------------
!     Introduce first point
!     Define vertex and adjacency lists for first 3 triangles
!
      p = List(1)
      V(1,1) = p
      V(2,1) = v1
      V(3,1) = v2
      T(1,1) = 3
      T(2,1) = 0
      T(3,1) = 2
      V(1,2) = p
      V(2,2) = v2
      V(3,2) = v3
      T(1,2) = 1
      T(2,2) = 0
      T(3,2) = 3
      V(1,3) = p
      V(2,3) = v3
      V(3,3) = v1
      T(1,3) = 2
      T(2,3) = 0
      T(3,3) = 1
!---------------------------------------------------------------------
!     Loop over each point (except first) and construct triangles
!
      Ntri = 3
      topstk = 0
      DO i = 2 , N
         p = List(i)
         xp = X(p)
         yp = Y(p)
!
!       Locate triangle J in which point lies
!
         j = Ntri
         DO
            v1 = V(1,j)
            v2 = V(2,j)
            IF ( (Y(v1)-yp)*(X(v2)-xp)>(X(v1)-xp)*(Y(v2)-yp) ) THEN
               j = T(1,j)
               CYCLE
            ENDIF
            v3 = V(3,j)
            IF ( (Y(v2)-yp)*(X(v3)-xp)>(X(v2)-xp)*(Y(v3)-yp) ) THEN
               j = T(2,j)
               CYCLE
            ENDIF
            IF ( (Y(v3)-yp)*(X(v1)-xp)>(X(v3)-xp)*(Y(v1)-yp) ) THEN
               j = T(3,j)
               CYCLE
            ENDIF
!
!       Check that new point is not coincident with any previous point
!
            d = (X(v1)-xp)**2
            IF ( d<TOL ) THEN
               d = d + (Y(v1)-yp)**2
               IF ( d<TOL ) WRITE (6,99002) v1 , p
            ENDIF
            d = (X(v2)-xp)**2
            IF ( d<TOL ) THEN
               d = d + (Y(v2)-yp)**2
               IF ( d<TOL ) WRITE (6,99002) v2 , p
            ENDIF
            d = (X(v3)-xp)**2
            IF ( d<TOL ) THEN
               d = d + (Y(v3)-yp)**2
               IF ( d<TOL ) WRITE (6,99002) v3 , p
            ENDIF
!
!       Create new vertex and adjacency lists for triangle J
!
            a = T(1,j)
            b = T(2,j)
            c = T(3,j)
            V(1,j) = p
            V(2,j) = v1
            V(3,j) = v2
            T(1,j) = Ntri + 2
            T(2,j) = a
            T(3,j) = Ntri + 1
!
!       Create new triangles
!
            Ntri = Ntri + 1
            V(1,Ntri) = p
            V(2,Ntri) = v2
            V(3,Ntri) = v3
            T(1,Ntri) = j
            T(2,Ntri) = b
            T(3,Ntri) = Ntri + 1
            Ntri = Ntri + 1
            V(1,Ntri) = p
            V(2,Ntri) = v3
            V(3,Ntri) = v1
            T(1,Ntri) = Ntri - 1
            T(2,Ntri) = c
            T(3,Ntri) = j
!
!       Put each edge of triangle J on STACK
!       Store triangles on left side of each edge
!       Update adjacency lists for triangles B and C
!
            topstk = topstk + 2
            Stack(topstk-1) = j
            IF ( T(1,c)==j ) THEN
               T(1,c) = Ntri
            ELSE
               T(2,c) = Ntri
            ENDIF
            Stack(topstk) = Ntri
            IF ( b/=0 ) THEN
               IF ( T(1,b)==j ) THEN
                  T(1,b) = Ntri - 1
               ELSEIF ( T(2,b)==j ) THEN
                  T(2,b) = Ntri - 1
               ELSE
                  T(3,b) = Ntri - 1
               ENDIF
               topstk = topstk + 1
               Stack(topstk) = Ntri - 1
            ENDIF
            EXIT
         ENDDO
!
!       Loop while STACK is not empty
!
         DO WHILE ( topstk>0 )
!
!         Find triangles L and R which are either side of stacked edge
!         triangle L is defined by V3-V1-V2 and is left of V1-V2
!         triangle R is defined by V4-V2-V1 and is right of V1-V2
 
            r = Stack(topstk)
            topstk = topstk - 1
            l = T(2,r)
!
!         Check if new point P is in circumcircle for triangle L
!
            IF ( T(1,l)==r ) THEN
               v1 = V(1,l)
               v2 = V(2,l)
               v3 = V(3,l)
               a = T(2,l)
               b = T(3,l)
            ELSEIF ( T(2,l)==r ) THEN
               v1 = V(2,l)
               v2 = V(3,l)
               v3 = V(1,l)
               a = T(3,l)
               b = T(1,l)
            ELSE
               v1 = V(3,l)
               v2 = V(1,l)
               v3 = V(2,l)
               a = T(1,l)
               b = T(2,l)
            ENDIF
            x13 = X(v1) - X(v3)
            y13 = Y(v1) - Y(v3)
            x23 = X(v2) - X(v3)
            y23 = Y(v2) - Y(v3)
            x1p = X(v1) - xp
            y1p = Y(v1) - yp
            x2p = X(v2) - xp
            y2p = Y(v2) - yp
            cosa = x13*x23 + y13*y23
            cosb = x2p*x1p + y1p*y2p
            IF ( cosa<C00000 .OR. cosb<C00000 ) THEN
               IF ( (cosa<C00000 .AND. cosb<C00000) .OR. &
     &              ((x13*y23-x23*y13)*cosb<(x1p*y2p-x2p*y1p)*cosa) )&
     &              THEN
!
!             New point is inside circumcircle for triangle L
!             Swap diagonal for convex quad formed by P-V2-V3-V1
!
                  c = T(3,r)
!
!             Update vertex and adjacency list for triangle R
!
                  V(3,r) = v3
                  T(2,r) = a
                  T(3,r) = l
!
!             Update vertex and adjacency list for triangle L
!
                  V(1,l) = p
                  V(2,l) = v3
                  V(3,l) = v1
                  T(1,l) = r
                  T(2,l) = b
                  T(3,l) = c
!
!             Put edges R-A and L-B on STACK
!             Update adjacency lists for triangles A and C
!
                  IF ( a/=0 ) THEN
                     IF ( T(1,a)==l ) THEN
                        T(1,a) = r
                     ELSEIF ( T(2,a)==l ) THEN
                        T(2,a) = r
                     ELSE
                        T(3,a) = r
                     ENDIF
                     topstk = topstk + 1
                     IF ( topstk>Npts ) THEN
                        WRITE (6,99001)
                        STOP
                     ENDIF
                     Stack(topstk) = r
                  ENDIF
                  IF ( b/=0 ) THEN
                     topstk = topstk + 1
                     IF ( topstk>Npts ) THEN
                        WRITE (6,99001)
                        STOP
                     ENDIF
                     Stack(topstk) = l
                  ENDIF
                  T(1,c) = l
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!---------------------------------------------------------------------
!     Check consistency of triangulation
!
      IF ( Ntri/=2*N+1 ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN DELAUN***'',                                   &
     &   /,''INCORRECT NUMBER OF TRIANGLES FORMED'')')
         STOP
      ENDIF
!---------------------------------------------------------------------
99001 FORMAT (//,'***ERROR IN SUBROUTINE DELAUN***',/,'STACK OVERFLOW')
99002 FORMAT (//,'***WARNING IN DELAUN***',/,'POINTS',I5,' AND',I5,&
     &        ' ARE COINCIDENT',/,&
     &        'DELETE EITHER POINT FROM LIST VECTOR')
      END SUBROUTINE DELAUN
!*==edge.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE EDGE(Vi,Vj,Npts,Ntri,Nef,Jw,X,Y,V,T,Tn,W)
      IMPLICIT NONE
!*--EDGE1143
!***********************************************************************
!
!     PURPOSE:
!     --------
!
!     Force edge VI-VJ to be present in Delaunay triangulation
!
!     INPUT:
!     ------
!
!     VI,VJ  - Vertices defining edge to be present in triangulation
!     NPTS   - Total number of points in data set
!     NTRI   - Number of triangles in triangulation
!     NEF    - Running total of edges that have been forced
!            - Set to zero before first call to EDGE
!     JW     - Number of cols in workspace vector W
!            - JW must not be less than the number of edges in the
!              existing triangulation that intersect VI-VJ
!     X      - X-coords of all points in data set
!            - X-coord of point I given by X(I)
!            - Last three locations are used to store x-coords of
!              supertriangle vertices in subroutine delaun
!     Y      - Y-coords of all points in data set
!            - Y-coord of point I given by Y(I)
!            - Last three locations are used to store y-coords of
!              supertriangle vertices in subroutine delaun
!     V      - Vertex array for triangulation
!            - Vertices listed in anticlockwise sequence
!            - Vertices for triangle J are found in V(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - First vertex is at point of contact of first and third
!              adjacent triangles
!     T      - Adjacency array for triangulation
!            - Triangles adjacent to J are found in T(I,J) for I=1,2,3
!              J=1,2,...,NTRI
!            - Adjacent triangles listed in anticlockwise sequence
!            - Zero denotes no adjacent triangle
!     TN     - List of triangle numbers such that vertex I can be
!              found in triangle TN(I)
!     W      - Not defined, used as workspace
!
!     OUTPUT:
!     -------
!
!     VI,VJ  - Unchanged
!     NPTS   - Unchanged
!     NTRI   - Unchanged
!     NEF    - If VI-VJ needs to be forced, NEF is incremented by unity
!            - Else NEF is unchanged
!     JW     - Unchanged
!     X      - Unchanged
!     Y      - Unchanged
!     V      - Vertex array for triangulation updated so that edge
!              V1-V2 is present
!     T      - Adjacency array for triangulation updated so that edge
!              V1-V2 is present
!     TN     - List of triangle numbers updated so that edge
!              V1-V2 is present
!     W      - List of new edges that replace old edges in
!              triangulation
!            - Vertices in W(1,I) and W(2,I) define each new edge I
!
!     NOTES:
!     ------
!
!     - This routine assumes that the edge defined by VI-VJ does not
!       lie on an outer boundary of the triangulation and, thus, the
!       triangulation must use the triangles that are formed with
!       the supertriangle vertices
!
!     PROGRAMMER:             Scott Sloan
!     -----------
!
!     LAST MODIFIED:          3 march 1991          Scott Sloan
!     --------------
!
!     COPYRIGHT 1990:         Scott Sloan
!     ---------------         Department of Civil Engineering
!                             University of Newcastle
!                             NSW 2308
!
!***********************************************************************
      INTEGER a , c , e , i , l , r , s
      INTEGER Jw , nc , v1 , v2 , v3 , v4 , Vi , Vj
      INTEGER elr , erl , Nef
      INTEGER last , Ntri , Npts
      INTEGER first
      INTEGER T(3,Ntri) , V(3,Ntri) , W(2,Jw)
      INTEGER Tn(Npts+3)
!
      DOUBLE PRECISION x1 , x2 , x3 , x4 , xi , xj , y1 , y2 , y3 , y4 ,&
     &                 yi , yj
      DOUBLE PRECISION TOL , x13 , x14 , x23 , x24 , y13 , y14 , y24 , &
     &                 y23
      DOUBLE PRECISION cosa , cosb
      DOUBLE PRECISION C00000 , detij3
      DOUBLE PRECISION X(Npts+3) , Y(Npts+3)
!
      LOGICAL swap
!
      PARAMETER (TOL=1.E-6)
      PARAMETER (C00000=0.0)
!---------------------------------------------------------------------
!     Check data
!
      IF ( Vi<=0 .OR. Vi>Npts ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''ILLEGAL VALUE FOR VI'',                                    &
     &   /,''VI='',I5)') Vi
         STOP
      ENDIF
      IF ( Vj<=0 .OR. Vj>Npts ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''ILLEGAL VALUE FOR VJ'',                                    &
     &   /,''VJ='',I5)') Vj
         STOP
      ENDIF
!----------------------------------------------------------------------
!     Find any triangle which has VI as a vertex
!
      s = Tn(Vi)
      IF ( s<=0 ) THEN
         WRITE (6,99001) Vi
         STOP
      ENDIF
      IF ( Tn(Vj)<=0 ) THEN
         WRITE (6,99001) Vj
         STOP
      ENDIF
      xi = X(Vi)
      yi = Y(Vi)
      xj = X(Vj)
      yj = Y(Vj)
!----------------------------------------------------------------------
!     Find an arc that crosses VI-VJ
!     C=current triangle
!     S=triangle in which search is started
!
      c = s
      DO
!
!     Vertices V1 and V2 are such that V1-V2 is opposite VI
!     Circle anticlockwise round VI until V1-V2 crosses VI-VJ
!
         IF ( V(1,c)==Vi ) THEN
            v2 = V(3,c)
            e = 2
         ELSEIF ( V(2,c)==Vi ) THEN
            v2 = V(1,c)
            e = 3
         ELSE
            v2 = V(2,c)
            e = 1
         ENDIF
!
!     Test if arc VI-VJ already exists
!
         IF ( v2==Vj ) RETURN
!
!     Test if V1-V2 crosses VI-VJ
!
         x2 = X(v2)
         y2 = Y(v2)
         IF ( (xi-x2)*(yj-y2)>(xj-x2)*(yi-y2) ) THEN
!
!       V2 is left of VI-VJ
!       Check if V1 is right of VI-VJ
!
            v1 = V(e,c)
            x1 = X(v1)
            y1 = Y(v1)
            IF ( (xi-x1)*(yj-y1)<(xj-x1)*(yi-y1) ) THEN
!
!         V1-V2 crosses VI-VJ , so edge needs to be forced
!
               Nef = Nef + 1
!-------------------------------------------------------------------
!     Loop to store all arcs which cross VI-VJ
!     Vertices V1/V2 are right/left of   VI-VJ
!
               nc = 0
               DO
                  nc = nc + 1
                  IF ( nc>Jw ) THEN
                     WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''NOT ENOUGH WORKSPACE'',                                    &
     &   /,''INCREASE JW'')')
                     STOP
                  ENDIF
                  W(1,nc) = v1
                  W(2,nc) = v2
                  c = T(e,c)
                  IF ( V(1,c)==v2 ) THEN
                     v3 = V(3,c)
                     e = 2
                  ELSEIF ( V(2,c)==v2 ) THEN
                     v3 = V(1,c)
                     e = 3
                  ELSE
                     v3 = V(2,c)
                     e = 1
                  ENDIF
!
!     Termination test, all arcs crossing VI-VJ have been stored
!
                  IF ( v3==Vj ) THEN
!-------------------------------------------------------------------
!     Swap each arc that crosses VI-VJ if it is a diagonal of a
!     convex quadrilateral
!     Execute all possible swaps, even if newly formed arc also
!     crosses VI-VJ, and iterate until no arcs cross
!
                     last = nc
                     DO
                        IF ( last>0 ) THEN
                           first = 1
                           DO
                              IF ( first<=last ) THEN
!
!         Find triangle L which is left of V1-V2
!         Find triangle R which is right of V1-V2
!
                                 v1 = W(1,first)
                                 v2 = W(2,first)
!
!         Exchange V1 and V2 if V1 is a supertriangle vertex
!
                                 IF ( v1>Npts ) THEN
                                    IF ( v2>Npts ) THEN
                                       WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''ARC BETWEEN VERTICES'',I5,'' AND'',I5,                     &
     &   /,''CROSSES SUPERTRIANGLE BOUNDARY DEFINED '',                 &
     &     ''BY VERTICES'',I5,''AND'',I5)') Vi , Vj , v1 , v2
                                       STOP
                                    ENDIF
                                    W(1,first) = v2
                                    W(2,first) = v1
                                    v2 = v1
                                    v1 = W(1,first)
                                 ENDIF
                                 l = Tn(v1)
                                 DO
                                    IF ( V(1,l)==v1 ) THEN
                                       IF ( V(2,l)==v2 ) THEN
                                         v3 = V(3,l)
                                         elr = 1
                                         r = T(1,l)
                                       ELSE
                                         l = T(3,l)
                                         CYCLE
                                       ENDIF
                                    ELSEIF ( V(2,l)==v1 ) THEN
                                       IF ( V(3,l)==v2 ) THEN
                                         v3 = V(1,l)
                                         elr = 2
                                         r = T(2,l)
                                       ELSE
                                         l = T(1,l)
                                         CYCLE
                                       ENDIF
                                    ELSEIF ( V(1,l)==v2 ) THEN
                                       v3 = V(2,l)
                                       elr = 3
                                       r = T(3,l)
                                    ELSE
                                       l = T(2,l)
                                       CYCLE
                                    ENDIF
!
!         Find vertices V3 and V4 where:
!         triangle L is defined by V3-V1-V2
!         triangle R is defined by V4-V2-V1
!
                                    IF ( T(1,r)==l ) THEN
                                       v4 = V(3,r)
                                       a = T(2,r)
                                       erl = 1
                                    ELSEIF ( T(2,r)==l ) THEN
                                       v4 = V(1,r)
                                       a = T(3,r)
                                       erl = 2
                                    ELSE
                                       v4 = V(2,r)
                                       a = T(1,r)
                                       erl = 3
                                    ENDIF
!
!         Test if quad formed by V3-V1-V4-V2 is convex
!
                                    x3 = X(v3)
                                    y3 = Y(v3)
                                    x4 = X(v4)
                                    y4 = Y(v4)
                                    x1 = X(v1)
                                    y1 = Y(v1)
                                    x2 = X(v2)
                                    y2 = Y(v2)
                                    IF ( (x3-x1)*(y4-y1)<(x4-x1)*(y3-y1)&
     &                                 ) THEN
                                       IF ( (x3-x2)*(y4-y2)>(x4-x2)&
     &                                    *(y3-y2) ) THEN
!
!             Quad is convex so swap diagonal arcs
!             Update vertex and adjacency lists for triangle L
!
                                         IF ( elr==1 ) THEN
                                         V(2,l) = v4
                                         c = T(2,l)
                                         T(1,l) = a
                                         T(2,l) = r
                                         ELSEIF ( elr==2 ) THEN
                                         V(3,l) = v4
                                         c = T(3,l)
                                         T(2,l) = a
                                         T(3,l) = r
                                         ELSE
                                         V(1,l) = v4
                                         c = T(1,l)
                                         T(3,l) = a
                                         T(1,l) = r
                                         ENDIF
!
!             Update vertex and adjacency lists for triangle R
!
                                         IF ( erl==1 ) THEN
                                         V(2,r) = v3
                                         T(1,r) = c
                                         T(2,r) = l
                                         ELSEIF ( erl==2 ) THEN
                                         V(3,r) = v3
                                         T(2,r) = c
                                         T(3,r) = l
                                         ELSE
                                         V(1,r) = v3
                                         T(3,r) = c
                                         T(1,r) = l
                                         ENDIF
!
!             Update adjacency lists for triangles A and C
!
                                         IF ( T(1,c)==l ) THEN
                                         T(1,c) = r
                                         ELSEIF ( T(2,c)==l ) THEN
                                         T(2,c) = r
                                         ELSE
                                         T(3,c) = r
                                         ENDIF
                                         IF ( T(1,a)==r ) THEN
                                         T(1,a) = l
                                         ELSEIF ( T(2,a)==r ) THEN
                                         T(2,a) = l
                                         ELSE
                                         T(3,a) = l
                                         ENDIF
!
!             Update vertex-triangle list
!
                                         Tn(v1) = l
                                         Tn(v2) = r
!
!             Test if new diagonal arc crosses VI-VJ and store it if it
!             does
!
                                         IF ( ((xi-x3)*(yj-y3)-(xj-x3)*(&
     &                                      yi-y3))&
     &                                      *((xi-x4)*(yj-y4)-(xj-x4)&
     &                                      *(yi-y4))<C00000 ) THEN
                                         W(1,first) = v4
                                         W(2,first) = v3
                                         first = first + 1
                                         ELSE
                                         W(1,first) = W(1,last)
                                         W(1,last) = v3
                                         W(2,first) = W(2,last)
                                         W(2,last) = v4
                                         last = last - 1
                                         ENDIF
                                         GOTO 2
                                       ENDIF
                                    ENDIF
!
!         Arc cannot be swapped, so move to next intersecting arc
!
                                    first = first + 1
                                    GOTO 2
                                 ENDDO
                              ENDIF
                              GOTO 4
 2                         ENDDO
                        ENDIF
!----------------------------------------------------------------------
!     Optimise all new arcs (except VI-VJ)
!
                        swap = .TRUE.
                        EXIT
 4                   ENDDO
                     DO
                        IF ( swap ) THEN
                           swap = .FALSE.
                           DO i = 2 , nc
!
!         Find triangle L which is left of V1-V2
!         Find triangle R which is right of V1-V2
!
                              v1 = W(1,i)
                              v2 = W(2,i)
!
!         Exchange V1 and V2 if V1 is a supertriangle vertex
!
                              IF ( v1>Npts ) THEN
                                 IF ( v2>Npts ) THEN
                                    WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''ARC BETWEEN VERTICES'',I5,'' AND'',I5,                     &
     &   /,''CANNOT BE OPTIMISED SINCE IT IS A '',                      &
     &     ''SUPERTRIANGLE BOUNDARY'')') v1 , v2
                                    STOP
                                 ENDIF
                                 W(1,i) = v2
                                 W(2,i) = v1
                                 v2 = v1
                                 v1 = W(1,i)
                              ENDIF
                              l = Tn(v1)
                              DO
                                 IF ( V(1,l)==v1 ) THEN
                                    IF ( V(2,l)==v2 ) THEN
                                       v3 = V(3,l)
                                       elr = 1
                                       r = T(1,l)
                                    ELSE
                                       l = T(3,l)
                                       CYCLE
                                    ENDIF
                                 ELSEIF ( V(2,l)==v1 ) THEN
                                    IF ( V(3,l)==v2 ) THEN
                                       v3 = V(1,l)
                                       elr = 2
                                       r = T(2,l)
                                    ELSE
                                       l = T(1,l)
                                       CYCLE
                                    ENDIF
                                 ELSEIF ( V(1,l)==v2 ) THEN
                                    v3 = V(2,l)
                                    elr = 3
                                    r = T(3,l)
                                 ELSE
                                    l = T(2,l)
                                    CYCLE
                                 ENDIF
!
!         Find vertices V3 and V4 where:
!         triangle L is defined by V3-V1-V2
!         triangle R is defined by V4-V2-V1
!
                                 IF ( T(1,r)==l ) THEN
                                    v4 = V(3,r)
                                    a = T(2,r)
                                    erl = 1
                                 ELSEIF ( T(2,r)==l ) THEN
                                    v4 = V(1,r)
                                    a = T(3,r)
                                    erl = 2
                                 ELSE
                                    v4 = V(2,r)
                                    a = T(1,r)
                                    erl = 3
                                 ENDIF
                                 x13 = X(v1) - X(v3)
                                 y13 = Y(v1) - Y(v3)
                                 x14 = X(v1) - X(v4)
                                 y14 = Y(v1) - Y(v4)
                                 x23 = X(v2) - X(v3)
                                 y23 = Y(v2) - Y(v3)
                                 x24 = X(v2) - X(v4)
                                 y24 = Y(v2) - Y(v4)
                                 cosa = x13*x23 + y23*y13
                                 cosb = x24*x14 + y24*y14
                                 IF ( cosa<C00000 .OR. cosb<C00000 )&
     &                                THEN
                                    IF ( (cosa<C00000 .AND. cosb<C00000)&
     &                                 .OR. &
     &                                 ((x13*y23-x23*y13)*cosb-(x14*y24-&
     &                                 x24*y14)&
     &                                 *cosa<-TOL*SQRT((x13*x13+y13*y13)&
     &                                 *(x23*x23+y23*y23)&
     &                                 *(x24*x24+y24*y24)&
     &                                 *(x14*x14+y14*y14))) ) THEN
!
!             V4 is inside circumcircle for triangle L
!             Swap diagonal for convex quad formed by V3-V1-V4-V2
!             Update vertex and adjacency lists for triangle L
!
                                       swap = .TRUE.
                                       IF ( elr==1 ) THEN
                                         V(2,l) = v4
                                         c = T(2,l)
                                         T(1,l) = a
                                         T(2,l) = r
                                       ELSEIF ( elr==2 ) THEN
                                         V(3,l) = v4
                                         c = T(3,l)
                                         T(2,l) = a
                                         T(3,l) = r
                                       ELSE
                                         V(1,l) = v4
                                         c = T(1,l)
                                         T(3,l) = a
                                         T(1,l) = r
                                       ENDIF
!
!             Update vertex and adjacency lists for triangle R
!
                                       IF ( erl==1 ) THEN
                                         V(2,r) = v3
                                         T(1,r) = c
                                         T(2,r) = l
                                       ELSEIF ( erl==2 ) THEN
                                         V(3,r) = v3
                                         T(2,r) = c
                                         T(3,r) = l
                                       ELSE
                                         V(1,r) = v3
                                         T(3,r) = c
                                         T(1,r) = l
                                       ENDIF
!
!             Update adjacency lists for triangles A and C
!
                                       IF ( T(1,c)==l ) THEN
                                         T(1,c) = r
                                       ELSEIF ( T(2,c)==l ) THEN
                                         T(2,c) = r
                                       ELSE
                                         T(3,c) = r
                                       ENDIF
                                       IF ( T(1,a)==r ) THEN
                                         T(1,a) = l
                                       ELSEIF ( T(2,a)==r ) THEN
                                         T(2,a) = l
                                       ELSE
                                         T(3,a) = l
                                       ENDIF
!
!             Update vertex-triangle list and arc list
!
                                       Tn(v1) = l
                                       Tn(v2) = r
                                       W(1,i) = v3
                                       W(2,i) = v4
                                    ENDIF
                                 ENDIF
                                 EXIT
                              ENDDO
                           ENDDO
                           CYCLE
                        ENDIF
                        GOTO 99999
                     ENDDO
                  ELSE
                     x3 = X(v3)
                     y3 = Y(v3)
                     detij3 = (xi-x3)*(yj-y3) - (xj-x3)*(yi-y3)
                     IF ( detij3<C00000 ) THEN
                        e = MOD(e,3) + 1
                        v1 = v3
                     ELSEIF ( detij3>C00000 ) THEN
                        v2 = v3
                     ELSE
                        WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''VERTEX'',I5,'' IS ON ARC'',                                &
     &   /,''BETWEEN VERTICES'',I5,'' AND'',I5)') v3 , Vi , Vj
                        STOP
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
!
!     No crossing, move anticlockwise around VI to the next triangle
!
         c = T(MOD(e,3)+1,c)
         IF ( c==s ) THEN
            WRITE (6,&
     &'(//,''***ERROR IN EDGE***'',                                     &
     &   /,''VERTEX ADJACENT TO'',I5,'' IS ON ARC'',                    &
     &   /,''BETWEEN VERTICES'',I5,'' AND'',I5)') Vi , Vi , Vj
            STOP
         ENDIF
      ENDDO
!---------------------------------------------------------------------
99001 FORMAT (//,'***ERROR IN EDGE***',/,'   VERTEX',I5,&
     &        ' NOT IN ANY TRIANGLE')
99999 END SUBROUTINE EDGE
!*==tcheck.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE TCHECK(Npts,N,X,Y,List,V,T,Ntri,Nef,Nb,Nce,Ncb,Elist,W)
      IMPLICIT NONE
!*--TCHECK1748
!***********************************************************************
!
!     PURPOSE:
!     --------
!
!     Check Delaunay triangulation which may be constrained
!
!     INPUT:
!     ------
!
!     NPTS   - Total number of points in data set
!     N      - Total number of points to be triangulated (N le NPTS)
!     X      - X-coords of all points in data set
!            - X-coord of point I given by X(I)
!            - Last three locations are used to store x-coords of
!              supertriangle vertices in subroutine delaun
!     Y      - Y-coords of all points in data set
!            - Y-coord of point I given by Y(I)
!            - Last three locations are used to store y-coords of
!              supertriangle vertices in subroutine delaun
!     LIST   - List of points to be triangulated
!     V      - Vertex array for triangulation
!            - Vertices listed in anticlockwise sequence
!            - Vertices for triangle J are found in V(I,J) for I=1,2,3
!              and J=1,2,...,NTRI
!            - First vertex is at point of contact of first and third
!              adjacent triangles
!     T      - Adjacency array for triangulation
!            - Triangles adjacent to J are found in T(I,J) for I=1,2,3
!              J=1,2,...,NTRI
!            - Adjacent triangles listed in anticlockwise sequence
!            - Zero denotes no adjacent triangle
!     NTRI   - Number of triangles in triangulation
!     NEF    - Number of forced edges in triangulation
!     NB     - Number of boundaries defining the triangulation
!            - NB=1 for a simply connected domain with no holes
!            - NB=H+1 for a domain with H holes
!     NCE    - Total number of constrained edges which must be
!              present, including those which define boundaries
!            - NCE=0 indicates triangulation is unconstrained so that
!              ELIST is empty and the code will produce a
!              triangulation of the convex hull
!     NCB    - Total number of constrained edges which define one
!              external boundary and any internal boundaries (holes)
!            - NCB=0 indicates there are no boundary edge constraints
!              and the code will produce a triangulation of the convex
!              hull
!     ELIST  - List of edges which must be present in triangulation
!            - These may be internal edges or edges which define
!              boundaries
!     W      - Undefined, vector used as workspace
!
!     OUTPUT:
!     -------
!
!     NPTS   - Unchanged
!     N      - Unchanged
!     X      - Unchanged
!     Y      - Unchanged
!     LIST   - Unchanged
!     V      - Unchanged
!     T      - Unchanged
!     NTRI   - Unchanged
!     NEF    - Unchanged
!     NB     - Unchanged
!     NCE    - Unchanged
!     NCB    - Unchanged
!     ELIST  - Unchanged
!     W      - Not used
!
!     NOTES:
!     ------
!
!     - This routine performs a number of elementary checks to test
!       the integrity of the triangulation
!     - NTRI=2*(N+NB)-NBOV-4 for a valid triangulation, where NBOV is
!       the number of boundary vertices
!     - NEDG=N+NTRI+NB-2 for a valid triangulation, where NEDG is the
!       number of edges in the triangulation
!     - NOPT le NEF for a valid triangulation, where NOPT is the number
!       of non-optimal edges in the triangulation
!     - The triangulation is tested to ensure that each non-boundary
!       constrained edge (if there are any) is present
!
!     PROGRAMMER:             Scott Sloan
!     -----------
!
!     LAST MODIFIED:          3 march 1991          Scott Sloan
!     --------------
!
!     COPYRIGHT 1990:         Scott Sloan
!     ---------------         Department of Civil Engineering
!                             University of Newcastle
!                             NSW 2308
!
!***********************************************************************
      INTEGER i , j , l , N , r , s
      INTEGER Nb , v1 , v2 , v3 , v4
      INTEGER nbov , nedg , nopt , Npts , Ntri
      INTEGER Ncb , Nce , Nef
      INTEGER T(3,Ntri) , V(3,Ntri) , W(Npts)
      INTEGER List(N)
      INTEGER Elist(2,*)
!
      DOUBLE PRECISION x1 , x2 , x3 , x4 , y1 , y2 , y3 , y4
      DOUBLE PRECISION det , dis , r21 , r31 , rad , TOL , x21 , x31 , &
     &                 xcc , y21 , y31 , ycc
      DOUBLE PRECISION C00000 , CP5000
      DOUBLE PRECISION X(Npts+3) , Y(Npts+3)
!
      PARAMETER (TOL=1.E-5)
      PARAMETER (C00000=0.0,CP5000=0.5)
!-----------------------------------------------------------------------
!     Loop over each triangle and count the following
!     NOPT=number of edges which are not optimal
!     NEDG=number of edges in triangulation
!     NBOV=number of boundary vertices
!
      nopt = 0
      nedg = 0
      nbov = 0
      DO l = 1 , Ntri
!
!       Loop over each side of triangle
!
         DO i = 1 , 3
            r = T(i,l)
            IF ( r==0 ) THEN
               nedg = nedg + 1
               nbov = nbov + 1
            ELSEIF ( r<l ) THEN
               nedg = nedg + 1
               IF ( i==1 ) THEN
                  v1 = V(1,l)
                  v2 = V(2,l)
                  v3 = V(3,l)
               ELSEIF ( i==2 ) THEN
                  v1 = V(2,l)
                  v2 = V(3,l)
                  v3 = V(1,l)
               ELSE
                  v1 = V(3,l)
                  v2 = V(1,l)
                  v3 = V(2,l)
               ENDIF
!
!           Triangle L is left of V1-V2 and is defined by V3-V1-V2
!           Triangle R is right of V1-V2 and is defined by V4-V2-V1
!
               IF ( T(1,r)==l ) THEN
                  v4 = V(3,r)
               ELSEIF ( T(2,r)==l ) THEN
                  v4 = V(1,r)
               ELSE
                  v4 = V(2,r)
               ENDIF
!
!           Find circumcentre of triangle L and its circumcircle radius
!
               x1 = X(v1)
               y1 = Y(v1)
               x21 = X(v2) - x1
               y21 = Y(v2) - y1
               x31 = X(v3) - x1
               y31 = Y(v3) - y1
               det = x21*y31 - x31*y21
               IF ( det<=C00000 ) THEN
                  WRITE (6,&
     &'(//,''***WARNING IN TCHECK***'',                                 &
     &   /,''ZERO OR -VE AREA FOR TRIANGLE'',I5)') l
               ELSE
                  det = CP5000/det
                  r21 = x21*x21 + y21*y21
                  r31 = x31*x31 + y31*y31
                  xcc = det*(r21*y31-r31*y21)
                  ycc = det*(x21*r31-x31*r21)
                  rad = SQRT(xcc*xcc+ycc*ycc)
                  xcc = xcc + x1
                  ycc = ycc + y1
!
!             Check if V4 is inside circumcircle for triangle L
!
                  dis = SQRT((xcc-X(v4))**2+(ycc-Y(v4))**2)
                  IF ( rad-dis>TOL*rad ) nopt = nopt + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     Check triangulation is valid
!
      IF ( nbov<3 ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN TCHECK***'',                                   &
     &   /,''NUMBER BOUNDARY VERTICES LT 3'',                           &
     &   /,''NBOV='',I5)') nbov
         STOP
      ENDIF
      IF ( Ntri/=2*(N+Nb)-nbov-4 ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN TCHECK***'',                                   &
     &   /,''INVALID TRIANGULATION'',                                   &
     &   /,''NTRI IS NOT EQUAL TO 2*(N+NB)-NBOV-4'')')
         WRITE (6,99001)
         WRITE (6,99002)
         STOP
      ENDIF
      IF ( nedg/=N+Ntri+Nb-2 ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN TCHECK***'',                                   &
     &   /,''INVALID TRIANGULATION'',                                   &
     &   /,''NEDG IS NOT EQUAL TO N+NTRI+NB-2'')')
         WRITE (6,99001)
         WRITE (6,99002)
         STOP
      ENDIF
      IF ( nopt>Nef ) THEN
         WRITE (6,&
     &'(//,''***ERROR IN TCHECK***'',                                   &
     &   /,''INVALID TRIANGULATION'',                                   &
     &   /,''TOO MANY NON-OPTIMAL EDGES'')')
         WRITE (6,99001)
         WRITE (6,99002)
         STOP
      ENDIF
      IF ( Ncb>0 ) THEN
         IF ( Ncb/=nbov ) THEN
            WRITE (6,&
     &'(//,''***ERROR IN TCHECK***'',                                   &
     &   /,''INVALID TRIANGULATION'',                                   &
     &   /,''NO. BOUNDARY VERTICES NOT EQUAL TO NBC'')')
            WRITE (6,99001)
            WRITE (6,99002)
            STOP
         ENDIF
      ENDIF
!----------------------------------------------------------------------
!     Check that each node appears in at least one triangle
!
      DO i = 1 , Npts
         W(i) = 0
      ENDDO
      DO j = 1 , Ntri
         DO i = 1 , 3
            W(V(i,j)) = j
         ENDDO
      ENDDO
      DO i = 1 , N
         v1 = List(i)
         IF ( W(v1)==0 ) THEN
            WRITE (6,&
     &'(//,''***ERROR IN TCHECK***'',                                   &
     &   /,''INVALID TRIANGULATION'',                                   &
     &   /,''VERTEX'',I5,'' IS NOT IN ANY TRIANGLE'')') v1
            STOP
         ENDIF
      ENDDO
      IF ( Nce==Ncb ) RETURN
!----------------------------------------------------------------------
!     Check that all non-boundary constrained edges are present
!
      DO i = Ncb + 1 , Nce
         v1 = Elist(1,i)
         v2 = Elist(2,i)
         s = W(v1)
         r = s
         DO
!
!       Circle anticlockwise round V1 until
!       - an edge V1-V2 is found, or,
!       - we are back at the starting triangle, this indicates an error
!         since V1 is interior and all neighbours have been checked, or,
!       - a boundary edge is found, this indicates that V1 is a boundary
!         vertex and the search must be repeated by circling clockwise
!           round V1
!
            IF ( V(1,r)==v1 ) THEN
               IF ( V(3,r)==v2 ) EXIT
               r = T(3,r)
            ELSEIF ( V(2,r)==v1 ) THEN
               IF ( V(1,r)==v2 ) EXIT
               r = T(1,r)
            ELSEIF ( V(2,r)==v2 ) THEN
               EXIT
            ELSE
               r = T(2,r)
            ENDIF
            IF ( r==s ) THEN
               WRITE (6,&
     &'(//,''***WARNING IN TCHECK***'',                                 &
     &   /,''CONSTRAINED EDGE WITH VERTICES'',I5,                       &
     &     '' AND'',I5,'' IS NOT IN TRIANGULATION'')') v1 , v2
!
!         Edge V1-V2 is not in triangulation
!         Check if it crosses any other constrained edges
!
               x1 = X(v1)
               y1 = Y(v1)
               x2 = X(v2)
               y2 = Y(v2)
               DO j = Ncb + 1 , Nce
                  IF ( j/=i ) THEN
                     v3 = Elist(1,j)
                     v4 = Elist(2,j)
                     x3 = X(v3)
                     y3 = Y(v3)
                     x4 = X(v4)
                     y4 = Y(v4)
                     IF ( ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))&
     &                    *((x1-x4)*(y2-y4)-(x2-x4)*(y1-y4))<C00000 )&
     &                    THEN
                        IF ( ((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))&
     &                       *((x3-x2)*(y4-y2)-(x4-x2)*(y3-y2))<C00000 )&
     &                       THEN
!
!               Edge V1-V2 crosses edge V3-V4
!
                           WRITE (6,&
     &'(''IT INTERSECTS ANOTHER CONSTRAINED EDGE '',                    &
     &  ''WITH VERTICES'',I5,'' AND'',I5)') v3 , v4
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               EXIT
            ENDIF
            IF ( r<=0 ) THEN
!
!       V1 must be a boundary vertex
!       Circle clockwise round V1 until
!       - an edge V1-V2 is found, or,
!       - a boundary edge is found which indicates an error
!
               l = s
               DO
                  IF ( V(1,l)==v1 ) THEN
                     IF ( V(2,l)==v2 ) GOTO 100
                     l = T(1,l)
                  ELSEIF ( V(2,l)==v1 ) THEN
                     IF ( V(3,l)==v2 ) GOTO 100
                     l = T(2,l)
                  ELSEIF ( V(1,l)==v2 ) THEN
                     GOTO 100
                  ELSE
                     l = T(3,l)
                  ENDIF
                  IF ( l<=0 ) THEN
                     WRITE (6,&
     &'(//,''***WARNING IN TCHECK***'',                                 &
     &   /,''CONSTRAINED EDGE WITH VERTICES'',I5,                       &
     &     '' AND'',I5,'' IS NOT IN TRIANGULATION'')') v1 , v2
!
!       Edge V1-V2 is not in triangulation
!       Check if it crosses any other constrained edges
!
                     x1 = X(v1)
                     y1 = Y(v1)
                     x2 = X(v2)
                     y2 = Y(v2)
                     DO j = Ncb + 1 , Nce
                        IF ( j/=i ) THEN
                           v3 = Elist(1,j)
                           v4 = Elist(2,j)
                           x3 = X(v3)
                           y3 = Y(v3)
                           x4 = X(v4)
                           y4 = Y(v4)
                           IF ( ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))&
     &                          *((x1-x4)*(y2-y4)-(x2-x4)*(y1-y4))&
     &                          <C00000 ) THEN
                              IF ( ((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))&
     &                             *((x3-x2)*(y4-y2)-(x4-x2)*(y3-y2))&
     &                             <C00000 ) THEN
!
!             Edge V1-V2 crosses edge V3-V4
!
                                 WRITE (6,&
     &'(''IT INTERSECTS ANOTHER CONSTRAINED EDGE '',                    &
     &  ''WITH VERTICES'',I5,'' AND'',I5)') v3 , v4
                                 EXIT
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                     GOTO 100
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
 100  ENDDO
99001 FORMAT (//,'CHECK THE FOLLOWING CONDITIONS:',/,&
     &        '-EDGES WHICH DEFINE BOUNDARIES MUST COME FIRST IN',/,&
     &        ' ELIST AND THUS OCCUPY THE FIRST NCB COLUMNS',/,&
     &        '-EDGES WHICH DEFINE AN EXTERNAL BOUNDARY MUST BE',/,&
     &        ' LISTED ANTICLOCKWISE BUT MAY BE PRESENTED IN ANY',/,&
     &        ' ORDER',/,&
     &        '-EDGES WHICH DEFINE AN INTERNAL BOUNDARY (HOLE) MUST',/,&
     &        ' BE LISTED CLOCKWISE BUT MAY BE PRESENTED IN ANY',/,&
     &        ' ORDER',/,&
     &        '-AN INTERNAL BOUNDARY (HOLE) CANNOT BE SPECIFIED',/,&
     &        ' UNLESS AN EXTERNAL BOUNDARY IS ALSO SPECIFIED')
99002 FORMAT ('-ALL BOUNDARIES MUST FORM CLOSED LOOPS',/,&
     &        '-AN EDGE MAY NOT APPEAR MORE THAN ONCE IN ELIST',/,&
     &        '-AN EXTERNAL OR INTERNAL BOUNDARY MAY NOT CROSS',/,&
     &        ' ITSELF AND MAY NOT SHARE A COMMON EDGE WITH ANY',/,&
     &        ' OTHER BOUNDARY',/,&
     &        '-INTERNAL EDGES, WHICH ARE NOT MEANT TO DEFINE',/,&
     &        ' BOUNDARIES BUT MUST BE PRESENT IN THE FINAL',/,&
     &        ' TRIANGULATION, OCCUPY COLUMNS NCB+1,... ,NCE OF',/,&
     &        ' ELIST',/,&
     &        '-NO POINT IN THE LIST VECTOR MAY LIE OUTSIDE ANY',/,&
     &        ' EXTERNAL BOUNDARY OR INSIDE ANY INTERNAL BOUNDARY')
      END SUBROUTINE TCHECK
