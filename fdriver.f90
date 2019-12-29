!
!  Driver program for pdnet: 
!              Fortran subroutines for the
!              solution of linear minimum cost network flow problems.
!
!              Inputs a network with pdrdmc(), calls pdnet() to compute
!              minimum cost network flow, deals with error conditions 
!              with pdierr(), and prints output report with pdwrsl().

program pdnet_driver

  implicit none

  !     Authors: Luis Portugal (Univ. Coimbra - Portugal)
  !              Mauricio G. C. Resende (AT&T Labs Research - USA)
  !              Geraldo Veiga (AT&T Labs - USA)
  !              Joaquim Judice (Univ. Coimbra - Portugal)
  !              Joao Patricio (EST Tomar - Portugal)

  !     Reference: L.F. Portugal, M.G.C. Resende, G. Veiga, and 
  !     J.J. Judice, "A truncated primal-infeasible dual-feasible 
  !     interior point network flow method," Networks, vol. 35, 
  !     pp. 91-108, 2000.

  

  !--------------------------------------------------------------
  !   Integer variables:
  ! 
  !   na: number of arcs.
  !   nn: number of nodes.
  !--------------------------------------------------------------
 ! integer ::  na, nn, iprcnd, maxit, mxcgit, optgap,  optmf, optst, prttyp
  
  integer ::  na, nn
  integer, parameter :: maxn = 60000
  integer, parameter :: maxa = 800000

  !--------------------------------------------------------------
  !Array declaration:
  !
  !Problem data: 
  !               b     (1:maxn): array of capacities at each node.
  !               c     (1:maxa): array of costs at each arc.
  !               endn  (1:maxa): end vertex of edge.
  !               l     (1:maxa): lower bound of flow at each arc.
  !               optflo(1:maxa): array with optimal flow at each arc.
  !               strn  (1:maxa): start vertex of edge.
  !               u     (1:maxa): upper bound of flow at each arc.
  !--------------------------------------------------------------
  integer, dimension(1:maxn):: b
  integer, dimension(1:maxa):: c, endn, l, optflo, strn, u

  !--------------------------------------------------------------
  !               info(1:30): some integer statistics.
  !--------------------------------------------------------------
  integer, dimension (1:30) :: info

  !--------------------------------------------------------------
  !               dinfo(1:35): some double precision statistics.
  !--------------------------------------------------------------
  double precision , dimension (1:35) :: dinfo

  !--------------------------------------------------------------
  !     Start of executable portion of main program.
  !
  !     Read the problem data. This procedure reads the 
  !     problem in standard DIMACS format. The network has nn nodes and
  !     na arcs.  Nodes of arc(i) are placed in strn(i) and endn(i).
  !     The cost, upper bound, and lower bound of arc(i) are placed in 
  !     c(i), u(i), and l(i), respectively.  The supply/demand of node(i) 
  !     is placed in b(i).
  !     ------------------------------------------------------------------
  
  call pdnet_read (nn, na, maxn, maxa, strn, endn, b, c, l, u)
  
  !     ------------------------------------------------------------------      
  !     Startup definitions.  These paramters are placed in array info.
  !     
  !     Required info parameters for input are:
  !
  !          info(23) - maximum number of arcs
  !          info(24) - maximum number of nodes
  !          info( 3) = fortran output unit
  !          info( 2) = selected preconditioner 
  !          info(19) = output level
  !          info( 3) = unit of output file (6=standard output)
  !     ------------------------------------------------------------------
  !
  !     ------------------------------------------------------------------
  !     Set preconditioner to be used.
  !
  !     info( 2) = -1 : diagonal preconditioner.
  !     info( 2) = -2 : IQRD preconditioner.
  !     info( 2) = -3 : Spanning tree preconditioner.
  !     info( 2) = -4 : diagonally compensated spanning tree 
  !                     preconditioner.
  !     ------------------------------------------------------------------
  !
  !     ------------------------------------------------------------------
  !     Set output level.
  !
  !     info(19) =  0 : completely silent.
  !     info(19) =  1 : only prints iteration progression.
  !     info(19) =  2 : also prints some statistics at each iteration.
  !     info(19) =  3 : same as previous, but also prints solution.
  !     ------------------------------------------------------------------

  call pdnet_setintparm (info)

  !     ------------------------------------------------------------------
  !     Solve the problem.  If optimal solution is found, optimal flow 
  !     on arc(i) is returned in optflo(i).  The error condition is in
  !     info(50).  On return, array info contains statistics in addition
  !     to parameters:
  !
  !           dinfo( 1): dual objective function value
  !           dinfo( 2): dual objective function value on tree
  !           dinfo( 5): p-d gap
  !           dinfo(10): miu
  !            info(10): p-d iterations
  !           dinfo(15): primal objective function
  !           dinfo(16): primal objective function value on tree
  !            info(15): error condition
  !     ------------------------------------------------------------------ 
 
  call pdnet (b(1:nn), c, dinfo, endn, info, l, na, nn, optflo(1:na), strn, u)
 

  !     ------------------------------------------------------------------
  !     Deal with error conditions.  
  !     ------------------------------------------------------------------

  if (info(15) /= 0) then
     call pdnet_error(info(15))
     stop
  end if

  !     ------------------------------------------------------------------
  !     Write output report.
  !     ------------------------------------------------------------------

  call pdnet_solreport(c, dinfo, endn, info, l, na, strn, u, optflo)

  !     ------------------------------------------------------------------
  !     End of main program.
  !     ------------------------------------------------------------------

end program pdnet_driver
