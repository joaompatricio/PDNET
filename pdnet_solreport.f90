!     ------------------------------------------------------------------
!     pdnet_solreport: Prints summary report.
!     ------------------------------------------------------------------

subroutine pdnet_solreport (c, dinfo, endn, info, l,  na, strn, u, x)

  implicit none

  !     ------------------------------------------------------------------  
  !     Integer input parameters:
  !
  !             na : number of arcs.
  !     ------------------------------------------------------------------ 
 
  integer :: na

  !     ------------------------------------------------------------------ 
  !     Integer input arrays:
  !             c     (1:na): array of costs at each arc.
  !             endn  (1:na): pointer of network data structure.
  !             l     (1:na): lower bound of flow at each arc.
  !             strn  (1:na): adjacency list of network data structure.
  !             u     (1:na): upper bound of flow at each arc.
  !             x     (1:na): array with optimal flow at each arc.
  !             info  (1:30): array with many statistics.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: c, endn, l, strn, u, x
  integer, dimension(1:30) :: info

  !     ------------------------------------------------------------------
  !     Double precision input arrays:
  !
  !             dinfo  (1:    35): array with many statistics.
  !     ------------------------------------------------------------------

  double precision, dimension(1:35) :: dinfo

  !     ------------------------------------------------------------------
  !     Integer work variables:
  !
  !             cgstop: flags the way PCG stopped.
  !             i     : loop counter.
  !             iprcnd: preconditioner used.
  !             itr   : interior-point method iteration count.
  !             iout  : output file.
  !             mfstat: maxflow status. 
  !             mxfe  : edges on optimal flow network. 
  !             mxfv  : vertices on optimal flow network.
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             opt   : flags optimality of solution. 
  !             pcgitr: PCG iteration count.
  !             pf    : iteration level.
  !             prttyp: printing level.
  !    ------------------------------------------------------------------
  
  integer :: cgstop, i, iprcnd, itr, iout, mfstat, mxfe, mxfv, opt, &
       pcgitr, pf, prttyp

  !     ------------------------------------------------------------------
  !     Double precision work variables.
  !    
  !             cgtol : tolerance for PCG stopping criterion.
  !             dobj  : dual objective function value on tree.
  !             dof   : dual objective function value on tree.
  !             gap   : value of the gap.
  !             mfdopt: maxflow value.
  !             miu   : value of interior-point mu parameter.
  !             pcgcos: value of cosine on PCG iteration. 
  !             pcgres: value of residual on PCG iteration.
  !             pcgtol: value of error on PCG iteration.
  !             pobj  : objective function value.
  !             pof   : primal objective function value on tree.
  !             tol1  : lower bound of the tolerance for maxflow.
  !             tol2  : upper bound of the tolerance for maxflow.  
  !     ------------------------------------------------------------------ 

  double precision :: cgtol , dobj  , dof   , gap   , mfdopt, miu   , &
       pcgcos, pcgres, pcgtol, pobj  ,  pof   , tol1  , tol2

  !     ------------------------------------------------------------------ 
  !     Character variables:
  !         
  !                ccgstp*3 : string corresponding to stopping of CG. 
  !                cmfstt*10: string corresponding to maxflow termination.
  !                precnd*4 : string corresponding to preconditioner.
  !     ------------------------------------------------------------------
  
  character(len= 3) :: ccgstp
  character(len=10) :: cmfstt
  character(len= 4) :: precnd

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------

  if (info(19) == 0) return
  
  itr    =  info(10)
  pcgitr =  info(13)
  gap    = dinfo( 5)
  pcgtol = dinfo(27)
  cgtol  = dinfo( 7)
  mfdopt = dinfo( 9)                 
  iprcnd =  info( 2)
  mxfv   =  info( 6)
  mxfe   =  info( 7)
  cgstop =  info(16)
  mfstat =  info(17)
  miu    = dinfo(10)             
  pcgcos = dinfo(29)
  pcgres = dinfo(14)
  dof    = dinfo( 2)
  pof    = dinfo(16)
  pf     =  info(14)
  opt    =  info(12)
  tol1   = dinfo(21)
  tol2   = dinfo(22)         
  pobj   = dinfo(15)
  dobj   = dinfo( 1)
  prttyp =  info(19)
  iout   =  info( 3)
 
  select case (cgstop)
  case(0)
     ccgstp = 'MAX'
  case(1)
     ccgstp = 'NOR'
  case(2)
     ccgstp = 'COS'
  case default
     ccgstp ='???'
  end select

  select case (mfstat)
  case(0)
     cmfstt = 'NONE'
  case(1)
     cmfstt = 'NONOPTIMAL'
  case(2)
     cmfstt = 'OPTIMAL'
  case default
     cmfstt = '?'
  end select

  select case (iprcnd)
  case(1)
     precnd='DIAG'
  case(2)
     precnd='IQRD'
  case(3)
     precnd='SPAN'
  case(4)
     precnd='DMST'
  case default
     precnd='?'
  end select

  write (iout,900)
900 format('----------------------------------------',&
         '----------------------------------------')
  if (prttyp >= 2) then
     write (iout,902) itr,miu
902  format('  CONV>     itr: ',i14,'      miu: ',e14.8,/)
     write (iout,901) pobj,dobj,gap
901  format('  COST>    pobj: ',e14.8,'     dobj: ',e14.8, '  pd-gap: ',e14.8,/)
     write (iout,903) precnd,ccgstp,pcgitr,pcgtol,pcgres,cgtol, pcgcos
903  format('    CG>  precnd: ',a4,7x,'      cgstop: ',a3,3x,&
          '          pcgitr: ',i14,/,&
          '    CG>  pcgtol: ',e14.8,'   pcgres: ',e14.8,/,&
          '    CG>   cgtol: ',e14.8,'   pcgcos: ',e14.8,/)
     if (pf == 1 .and. opt == 1) then
        write (iout,904)
904     format('  TREE> pb-feas: TRUE             pb-opt: TRUE')
     endif
     if (pf == 1 .and. opt == 0) then
        write (iout,905)
905     format('  TREE> pb-feas: TRUE             pb-opt: FALSE')
     endif
     if (pf == 0) then
        write (iout,906)
906     format('  TREE> pb-feas: FALSE')
     endif
  else
     write (iout,902) itr,miu
     if (pf == 1 .and. opt == 1) then
        write (iout,904)
     endif
     if (pf == 1 .and. opt == 0) then
        write (iout,905)
     endif
     if (pf == 0) then
        write (iout,906)
     endif
  endif
      
  if (pf == 1) then
     write (iout,907) pof,dof
907  format('  TREE>     pof: ',e14.8,'      dof: ',e14.8)
  endif
  write (iout,911)
  if (miu < 1.0d0) then
     write (iout,908) mxfv,mxfe
908  format(' MFLOW>    mxfv: ',i14,'     mxfe: ',i14)
  endif
  if (prttyp >= 2) then
     if (miu < 1.0d0) then
        write (iout,909) tol1,tol2,mfdopt
909     format(' MFLOW> pdlotol: ',e14.8,'  pdhitol: ',e14.8,&
             '    dobj: ',e14.8)
     endif
     write (iout,910) cmfstt
910  format(' MFLOW>  mfstat: ',a10)
  endif
  if (miu < 1.0d0) then
     write (iout,911)
911  format(' ')
  endif
  write (iout,900)
  if (prttyp <= 2) return
  write (iout,915)
915 format (' PDNET 1.0:  Primal-feasible,',&
         ' dual-infeasible interior point network flow method',/,&
         ' OPTIMAL primal flow',//,&
         '                from        to       arc     lower',&
         '     upper       arc',/,&
         '       arc      node      node      cost     bound',&
         '     bound      flow',/)
  do i = 1, na
     if (x(i) > 0) then
        write (iout,950) i,strn(i),endn(i),c(i),l(i),u(i)+l(i), x(i)+l(i)
950     format(7i10)
     endif
  end do
  write (iout,900)
  
  !     ------------------------------------------------------------------
  !     End of executable portion of subroutine pdnet_solreport.
  !     ------------------------------------------------------------------

end subroutine pdnet_solreport
