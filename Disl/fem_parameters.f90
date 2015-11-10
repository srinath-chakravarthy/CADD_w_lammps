!*==mod_fem_parameters.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 
!
! $Id: mod_fem_paramters,v 1.22004/04/2114:29:39 shastry Exp $
!
      MODULE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--MOD_FEM_PARAMETERS7
      INTEGER KNODE , NDOF , MAXGEO , MAXLMN , MAXEQS , MAXBND , MAXFIXED
      INTEGER maxface , MAXSEGM , MAXPAD
      PARAMETER (KNODE=3)
      PARAMETER (NDOF=2)
      PARAMETER (MAXGEO=100000)
      PARAMETER (MAXLMN=100000)
      PARAMETER (MAXEQS=NDOF*MAXGEO)
      PARAMETER (MAXBND=420)
      PARAMETER (MAXFIXED=50000)
      PARAMETER (MAXSEGM=40000)
      PARAMETER (MAXPAD=40000)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                       
! Common blocks                                                         
!  /setup_fem/ flag assuring that the arrays are initialized once       
!  /arrays_fem/:                                                        
!     imap(maxgeo) mapping from fem node_number to global node_number   
!     iconn(knode,maxlmn) connectivity matrix                           
!     iadj(knode, maxlmn) element ajacency matrix                       
!     ifixed(maxfixed) constrained degrees of freedom                   
!     ifixed_hold(2,maxfixed) map to nodes with constrained d.f.        
!     isegm(2,maxsegm) boundary segments                                
!                                                                       
!  /contrl/:                                                            
!      nnodes number of fem nodes                                       
!      nfixed the number of constrained degrees of freedom              
!      nelm number of elements                                          
!      nsegm number boundary segments                                   
!                                                                       
!  /data_fem/:                                                          
!       x0(ndof+1,maxgeo) coordinates of nodal points                   
!       a_stiff(maxbnd, maxeqs) stiffness matrix                        
!       ad_stiff(maxbnd, maxeqs) decomposed stiffness matrix            
!       mbandw matrix bandwidth                                         
!       nequ number of equations                                        
!                                                                       
!  /pad_cntrl/:                                                         
!       npad number of pad atoms                                        
!       padmap(maxpad) mapping from pad atom number to the gloal number 
!       padelmnt(maxpad) element the pad atom lies in                   
!                                                                       
!  /pad_data:                                                           
!       padtricoord(knode,maxpad) triangular coordinates of pad atoms   
!                                                                       
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!23456789012345678901234567890123456789012345678901234567890123456789012
!         1         2         3         4         5         6         7
!
      INTEGER imap(MAXGEO) , iconn(KNODE,MAXLMN) , iadj(KNODE,MAXLMN)
      INTEGER ifixed(MAXFIXED) , ifix_hold(2,MAXFIXED) , isegm(2,MAXSEGM)
      DOUBLE PRECISION x0(NDOF+1,MAXGEO) , absegm(3,MAXSEGM)
      DOUBLE PRECISION a_stiff(MAXBND,MAXEQS) , ad_stiff(MAXBND,MAXEQS)
      INTEGER nnodes , nfixed , nelm , nsegm
      INTEGER mbandw , nequ
      INTEGER i_flag
!
      INTEGER npad , padmap(MAXPAD) , padelmnt(MAXPAD)
      DOUBLE PRECISION padtricoord(KNODE,MAXPAD)
      END MODULE MOD_FEM_PARAMETERS
!
! $Log: mod_fem_paramters,v $
! Revision 1.2  2004/04/2114:29:39  shastry
! vijay-    mod_fem_paramters fem_alan.f: increased storage.
!
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.5  2002/06/0420:31:44  shilkrot
! 1. Rewrote fem solve (changed commons in fem_parameters and energy
! and numerical force computation.
! 2. Introduced negative element # andpenalty.
! 3. Add flag MoveDisl to fem_solve.
!
! Revision 1.4  2001/12/1307:31:24  shilkrot
! Implemented breadth first search to find the element number for
! a dislocation. Changed the interface of fe_locate to use the starting
! element for the search. Old fe_locate is in fem_services.
! Changed the interface of fem_setup. Now two arrays used as temp space 
! passed from outside as the last two parameters.
!
! Revision 1.3  2001/08/2203:18:35  shilkrot
! Fixed the expression for the energy and polished fem_alan a little bit
! This wersion works with dislocation passing.
!
! Revision 1.2  2001/07/1205:16:58  shilkrot
! Change to have implicit none in each subroutine
!
! Revision 1.1  2001/06/1800:22:23  shilkrot
! FEM parameters
!
!
