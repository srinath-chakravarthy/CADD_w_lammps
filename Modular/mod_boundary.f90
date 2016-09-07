!*==mod_boundary.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
      MODULE MOD_BOUNDARY
      IMPLICIT NONE
!*--MOD_BOUNDARY4
! elist is the model boundary, see explanation in CONTRI.F
! dbpoly contains the polygons where you want to put one or more
! detection bands.  Just allocate the necessary storage, define the
! polygon (CCW) and the code will find the elements necessary to make
! the band.  If the band should go right up to the free surface (like
! when the atomistic region is on the surface) define the polygon to
! extend out into free space and then come back into the mesh at the
! right place where the detection band should start again.
!
!  ndbpoly: number of detection band polygons.  Usually 1 for each
!           atomistic region
!  ndbvtx(ndbpoly): number of vertices for each det band polygons
!  dbpoly(2,ndbpoly):  coordinates of the vertices defining the det
!                      band polygons
!  eldb(:,ndbpoly):  element global # of the detection band element
!     on idb
!  nelidb(ndbpoly):  # of elements on idb ring
!  dbbound(ndb
      INTEGER nce , ncb , ncemax , ndbpoly
      INTEGER , ALLOCATABLE :: ndbvtx(:)
      INTEGER , ALLOCATABLE :: eldb(:,:)
      INTEGER , ALLOCATABLE :: nelidb(:)
      INTEGER , ALLOCATABLE :: elidb(:)
      INTEGER , POINTER :: elist(:,:)
      DOUBLE PRECISION , ALLOCATABLE :: dbpoly(:,:,:)
      LOGICAL , ALLOCATABLE :: dbbound(:,:) , dbboundnear(:)
      INTEGER , PARAMETER :: DBNMAX = 1800
 
      CONTAINS
!**
!**---------------------------------------------------------------
!**  IncreaseElist : Increase storage for constrained edges
!**
      SUBROUTINE INCREASEELIST(Nadd)
      IMPLICIT NONE
!*--INCREASEELIST40
      INTEGER Nadd
      INTEGER , POINTER :: e(:,:)
      ALLOCATE (e(2,NCEmax))
      e(1:2,1:NCEmax) = ELIst(1:2,1:NCEmax)
      DEALLOCATE (ELIst)
      ALLOCATE (ELIst(2,NCEmax+Nadd))
      ELIst(1:2,1:NCEmax) = e(1:2,1:NCEmax)
      NCEmax = NCEmax + Nadd
      DEALLOCATE (e)
      END SUBROUTINE INCREASEELIST
 
      SUBROUTINE ALLOCATE_DB
      IMPLICIT NONE
!*--ALLOCATE_DB54
      ALLOCATE (NDBvtx(NDBpoly))
      ALLOCATE (DBPoly(2,4,NDBpoly))
      ALLOCATE (DBBoundnear(NDBpoly))
      ALLOCATE (NELidb(NDBpoly))
      ALLOCATE (ELDb(DBNMAX,NDBpoly))
      ALLOCATE (DBBound(DBNMAX,NDBpoly))
 
      DBBoundnear(1:NDBpoly) = .FALSE.
      NDBvtx(1:NDBpoly) = 4
      DBBound = .FALSE.
      END SUBROUTINE ALLOCATE_DB
 
      END MODULE MOD_BOUNDARY
