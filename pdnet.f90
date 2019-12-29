subroutine pdnet( b, c, dinfo, endn, info, l, na, nn, optflo, strn, u)
  !     ------------------------------------------------------------------
  !     pdnet: This is the main subroutine for the PDNET package.
  !
  !            The subroutine starts by assigning working vectors to
  !            pointers, with pdpntr(). The control parameters are then
  !            read from the vector info, with pdginf(). Correctness of
  !            the input is checked with pdchkp(), and data structure
  !            created with pdmked(). Internal parameters for methods
  !            are set with pddflt(), and additional structures are
  !            created with pddstr(). Subroutines pdtrfm() and
  !            pdprtb() are called in order to shift the lower bounds to
  !            zero and to transform the data into double precision,
  !            respectively. Subroutines pdpdat() and pdckfs() check the
  !            problem characteristics and inquire if the network has
  !            enough capacity to transport the proposed amount of
  !            commodity. The primal-dual main loop is then started and
  !            the current iteration direction system right-hand-side is
  !            computed by pcrhs(). The maximum spanning tree is computed
  !            by pdhprm(), and its optimality is tested by pdoptm().
  !            Under certain conditions (miu < 1.0d0), the maxflow
  !            stopping criterion is invoked by subroutine pdchmf().
  !            The preconditioner is computed next, and if iqrd=-2,
  !            subroutine pdiqrd() is called to compute the IRQD
  !            preconditioner. The preconditioner switch comes next,
  !            and pdpcgd() is called afterwards, so that the linear
  !            system is solved. Preconditioner update is then made, and
  !            a summary of the iteration is printed by pdsout(). Primal
  !            and dual updates are made by pdupsl() and stopping criteria
  !            check takes place, before returning to the start of the
  !            iteration loop.
  !
  !     Authors: Luis Portugal (Univ. Coimbra - Portugal)
  !              Mauricio G. C. Resende (AT&T Labs Research - USA)
  !              Geraldo Veiga (AT&T Labs - USA)
  !              Joaquim Judice (Univ. Coimbra - Portugal)
  !              Joao Patricio (EST Tomar - Portugal)
  !
  !     Reference: L.F. Portugal, M.G.C. Resende, G. Veiga, and J.J. Judice,
  !     "A truncated primal-infeasible dual-feasible interior point network
  !     flow method," Networks, vol. 35, pp. 91-108, 2000.
  !
  !     ------------------------------------------------------------------
  
  implicit none
  
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !     ------------------------------------------------------------------
  integer :: na, nn
  integer :: i
  
  
  !     Integer input/output arrays:
  !
  !              b     (1:  nn): array of capacities at each node.
  !              c     (1:  na): array of costs at each arc.
  !              endn  (1:  na): pointer of network data structure.
  !              info  (1:  30): array with integer statistics and
  !                              parameters.
  !              l     (1:  na): lower bound of flow at each arc.
  !              optflo(1:  na): array with optimal flow at each arc.
  !              strn  (1:  na): adjacency list of network data structure.
  !              u     (1:  na): upper bound of flow at each arc.
  !     ------------------------------------------------------------------
  
  integer, dimension(1: nn) :: b
  integer, dimension(1: na) :: c, endn, l, optflo, strn, u
  integer, dimension(1: 30) :: info

  !   ------------------------------------------------------------------
  !   Integer allocatable arrays.
  !   ------------------------------------------------------------------

  integer, dimension(:), allocatable ::  bucket, bposit, bstats, father,&
       first, ifirst, innext, link , list  , lost  , lsets , mark  , oufrst,&
       ounext, pred  , rank  , rlink , rsets , rtrth , stedg , stptr , stvrtx,&
       wk2e1 , wrke1 , wrke2 , wrke3 , wrke4 , treath, trearc, unext , balance
  
  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             dinfo (1:    35): array with double precision statistics 
  !                               and parameters.
  !     ------------------------------------------------------------------
  
  double precision, dimension(1:   35) :: dinfo
  
  !   ------------------------------------------------------------------
  !   Real allocatable arrays.
  !   ------------------------------------------------------------------
  
  double precision, dimension(:), allocatable :: av1, av3, av4, av5, av6, av7, av8,&
       db, dc, diag , du  , key, ndiag, p   , q  , rhs, r , sd, s  , theta,&
       wd  , w  , wrkdv1, xd   , yd  , x  , y , zd , z , zz 

  !     ------------------------------------------------------------------
  !     Integer working variables:
  !
  !              bfirst: used for building maximum spanning tree.
  !              bound : used for computing aditional data structures.
  !              cgstop: flags the way PCG stopped.
  !              ierror: flags errors in PDNET.
  !              ioflag: controls the level of output of summary.
  !              iout  : output file unit.
  !              iprcnd: preconditioner used.
  !              larcs : auxiliary variable.
  !              lfirst: auxiliary variable.
  !              maxit : maximum number of primal-dual iterations.
  !              mfstat: maxflow status.
  !              mr    : estimate of maximum number of iterations.
  !              mxcgit: maximum number of PCG iterations.
  !              mxfe  : edges on optimal flow network.
  !              mxfv  : vertices on optimal flow network.
  !              nbound: used for computing aditional data structures.
  !              nbuk  : number of buckets.
  !              npfrst: First bounded arc in upnext list
  !              opt   : flags optimal flow.
  !              optgap: Duality gap optimality indicator flag.
  !              optmf : Maximum flow optimality indicator flag.
  !              optst : Spanning tree optimality indicator flag.
  !              pcgit : number of PCG iterations performed.
  !              pduit : number of primal-dual iterations performed.
  !              pf    : iteration level.
  !              prttyp: printing level.
  !              rc    : flags optimality status.
  !              root  : sets beginning of structure.
  !              uarcs : number of arcs at the upper bound.
  !              ufirst: auxiliary variable.
  !              upfrst: First unbounded arc in upnext list
  !     ------------------------------------------------------------------
  
  integer :: bfirst, bound , cgstop, ierror, ioflag, iout, iprcnd, larcs, &
       lfirst, maxit, mfstat, mr, mxcgit, mxfe  , mxfv  , nbound, nbuk,&
       npfrst, opt, optgap, optmf , optst ,  pcgit , pduit , pf, prttyp, &
       rc, root,  uarcs, ufirst, upfrst
  
  !     ----------------------------------------------------------------
  !     Double precision variables:
  !
  !             dobj  : dual objective function value.
  !             dof   : dual objective function value on tree.
  !             ds    : updating parameter for vector s.
  !             factor: factor for Newton step.
  !             fadd  : auxiliary variable for cost computation.
  !             fcttol: factorization zero tolerance
  !             gap   : value of the gap.
  !             grho0 : factor for rho0 update.
  !             gtolmf: update factor for tol1 and tol2.
  !             huge  : maxdouble precision number.
  !             maxb  : largest absolute value of b.
  !             maxc  : largest absolute value of c.
  !             maxu  : largest absolute value of u.
  !             mfdopt: maxflow value.
  !             minb  : smallest absolute value of b.
  !             minc  : smallest absolute value of c.
  !             minu  : smallest absolute value of u.
  !             miu   : value of interior-point mu parameter.
  !             mntlcg: lower bound for CG tolerance.
  !             mrho0 : lower bound for rho0.
  !             norb  : 1-norm of vector b.
  !             norc  : 1-norm of vector c.
  !             noru  : 1-norm of vector u.
  !             oldcos: value of cosine on previous PCG iteration.
  !             oldtol: value of residual on previous PCG iteration.
  !             pcos  : value of cosine on PCG iteration.
  !             pcgres: value of residual on PCG iteration.
  !             pobj  : objective function value.
  !             pof   : primal objective function value on tree.
  !             ps    : maximum primal stepsize.
  !             rgap  : relative primal-dual gap.
  !             rho0  : parameter from IQRD preconditioner.
  !             s1fctr: factor for miu.
  !             stptol: zero for dual slacks.
  !             stpval: largest slack on dual face.
  !             sw1   : iteration limit on CG for preconditioner switch.
  !             sw2   : another limit for switch.
  !             tol1  : lower bound of the tolerance for maxflow.
  !             tol2  : upper bound of the tolerance for maxflow.
  !             tolcg : tolerance for PCG stopping criterion.
  !             tolcg0: guess for tolerance for PCG stopping criterion.
  !             tolpcg: initial value for tolcg.
  !             weight: tree weight.
  !             zero  : value for zero.
  !     ------------------------------------------------------------------
  double precision :: dobj  , dof   , ds    , factor, fadd  ,fcttol, gap   , grho0 ,&
       gtolcg, gtolmf, huge  ,maxb  , maxc  , maxu  , mfdopt, minb  , minc  , minu,&
       miu   , mntlcg, mrho0 , norb  , norc  ,  noru  , oldcos, oldtol, pcgcos, &
       pcgres, pobj  ,  pof   , ps    , rgap  , rho0  , s1fctr, stptol, stpval, sw1, &
       sw2   , tol1  , tol2  , tolcg , tolcg0, tolpcg, tolslk, tolsw1, tolsw2, &
       weight,  zero

  !     ------------------------------------------------------------------
  !     Start of executable section of the subroutine.
  !     ------------------------------------------------------------------
  ierror = 0

  !     ------------------------------------------------------------------
  !     Set default values for internal parameters
  !     ------------------------------------------------------------------     

  call pdnet_default(dinfo, info, na, nn)

  call pdnet_getinfo( bound,  cgstop, dinfo,  dobj,   dof,    factor, fcttol, &
       gap, grho0,  gtolcg, gtolmf, huge, ierror, info,   ioflag, iout, &
       iprcnd, maxit,  mfdopt, mfstat, miu, mntlcg, mr,  mrho0,  mxcgit, &
       mxfe,   mxfv,   nbound, nbuk, oldcos, oldtol, opt, optgap, optmf, &
       optst, pcgcos, pcgit,  pcgres, pduit,  pf, pobj,  pof, prttyp, rho0, &
       root,   s1fctr, stptol,   stpval, tol1, tol2, tolcg,  tolcg0, tolpcg, &
       tolslk, tolsw1, tolsw2, sw1, sw2, zero )

  !     ------------------------------------------------------------------
  !     Print output report header
  !     ------------------------------------------------------------------
  
  call pdnet_header(dinfo, info, na, nn)

  !     ------------------------------------------------------------------
  !     Check correctness of input data and exit if error is detected.
  !     ------------------------------------------------------------------

  call pdnet_check (endn, ierror, na, nn, prttyp)

  !     ------------------------------------------------------------------
  !     Exit if error is detected.
  !     ------------------------------------------------------------------
  if (ierror /= 0) then
      call pdnet_putinfo( bound,  cgstop, dinfo,  dobj,   dof,    factor, fcttol, &
       gap, grho0,  gtolcg, gtolmf, huge, ierror, info,   ioflag, iout, &
       iprcnd, maxit,  mfdopt, mfstat, miu, mntlcg, mr,  mrho0,  mxcgit, &
       mxfe,   mxfv,   nbound, nbuk, oldcos, oldtol, opt, optgap, optmf, &
       optst, pcgcos, pcgit,  pcgres, pduit,  pf, pobj,  pof, prttyp, rho0, &
       root,   s1fctr, stptol,   stpval, tol1, tol2, tolcg,  tolcg0, tolpcg, &
       tolslk, tolsw1, tolsw2, sw1, sw2, zero )
      call pdnet_error(ierror)
      return
   end if

   !     ------------------------------------------------------------------
   !     Initialize maxflow status to inactive.
   !     ------------------------------------------------------------------
   mfstat = 0

   !   ------------------------------------------------------------------
   !     Variable ioflag controls the level of output of summary.
   !              ioflag=1: only number of arcs and nodes.
   !              ioflag=0: full output.
   !   ------------------------------------------------------------------
   ioflag = 0

   !     ------------------------------------------------------------------
   !     Construct additional data structures.
   !     ------------------------------------------------------------------
   !
   !     ------------------------------------------------------------------
   !     Allocate auxiliary integer vectors.
   !     ------------------------------------------------------------------
   allocate (bucket(0:mr))
   allocate (bposit(1:na))
   allocate (bstats(1:na))
   allocate (father(1:nn))
   allocate (first (1:nn))
   allocate (ifirst(1:nn))
   allocate (innext(1:na))
   allocate (link  (0:nn))
   allocate (list  (1:nn))
   allocate (lost  (1:nn))
   allocate (lsets (1:na))
   allocate (mark  (1:nn))
   allocate (oufrst(1:nn))
   allocate (ounext(1:na))
   allocate (pred  (1:nn))
   allocate (rank  (1:nn))
   allocate (rlink (0:nn))
   allocate (rsets (1:na))
   allocate (rtrth (1:na))
   allocate (stedg (1:na))
   allocate (stptr (1:na))
   allocate (stvrtx(1:na))
   allocate (wk2e1 (1:2*na))
   allocate (wrke1 (1:na))
   allocate (wrke2 (1:na))
   allocate (wrke3 (1:na))
   allocate (wrke4 (1:na))
   allocate (treath(1:nn))
   allocate (trearc(1:nn))
   allocate (unext (1:na))
   allocate (balance(1:na))

   !     ------------------------------------------------------------------
   !     Allocate auxiliary double precision vectors.
   !     ------------------------------------------------------------------
   allocate (av1   (1:na))
   allocate (av3   (1:na))
   allocate (av4   (1:na))
   allocate (av5   (1:na))
   allocate (av6   (1:na))
   allocate (av7   (1:na))
   allocate (av8   (1:na))
   allocate (db    (1:na))
   allocate (dc    (1:na))
   allocate (diag  (1:nn))
   allocate (du    (1:na))
   allocate (key   (1:nn + 1))
   allocate (ndiag (1:nn))
   allocate (p     (1:nn))
   allocate (q     (1:nn))
   allocate (rhs   (1:nn))
   allocate (r     (1:nn))
   allocate (sd    (1:na))
   allocate (s     (1:na))
   allocate (theta (1:na))
   allocate (wd    (1:na))
   allocate (w     (1:na))
   allocate (wrkdv1(1:nn))
   allocate (xd    (1:na))
   allocate (yd    (1:nn))
   allocate (x     (1:na))
   allocate (y     (1:nn))
   allocate (zd    (1:na))
   allocate (z     (1:na))
   allocate (zz    (1:nn))

   !-------------------------------------------------------------
   !Build data structures.
   !-------------------------------------------------------------

   call pdnet_datastruct (bound, bstats, endn, list, ifirst, innext, &
        na, nbound, nn, npfrst, lost, oufrst, ounext, strn, u, &
        upfrst, unext)
   
   do i = 1, na
      db(i) = 0.0d0
   end do

   !     ------------------------------------------------------------------
   !     Transform the problem by shifting the lower bounds to zero.
   !     ------------------------------------------------------------------

   call pdnet_transform (b, bound, c, endn, fadd, l, na, nbound, nn, npfrst,&
        strn, u, upfrst, unext)

   !     ------------------------------------------------------------------
   !     Introduce perturbations and transform the data to double 
   !     precision.
   !     ------------------------------------------------------------------   

   call pdnet_perturb (b, bound, c, db, dc, du, endn, na, nn, s, strn, u,&
        upfrst, unext, w, x, z)

   !     ------------------------------------------------------------------
   !     Find the problem data characteristics.
   !     ------------------------------------------------------------------

   call pdnet_probdata (db, bound, dc, huge, maxb, maxc, maxu, minb, minc, &
        minu, na, nn, norb, norc, noru, du, upfrst, unext)

   !     ------------------------------------------------------------------
   !     Check problem feasibility
   !     ------------------------------------------------------------------

   call pdnet_checkfeas(db, du, endn, ierror, na, nn, strn, wrke1, balance)

  !     ------------------------------------------------------------------
  !     Exit if error is detected.
  !     ------------------------------------------------------------------
  if (ierror /= 0) then
      call pdnet_putinfo( bound,  cgstop, dinfo,  dobj,   dof,    factor, fcttol, &
       gap, grho0,  gtolcg, gtolmf, huge, ierror, info,   ioflag, iout, &
       iprcnd, maxit,  mfdopt, mfstat, miu, mntlcg, mr,  mrho0,  mxcgit, &
       mxfe,   mxfv,   nbound, nbuk, oldcos, oldtol, opt, optgap, optmf, &
       optst, pcgcos, pcgit,  pcgres, pduit,  pf, pobj,  pof, prttyp, rho0, &
       root,   s1fctr, stptol,   stpval, tol1, tol2, tolcg,  tolcg0, tolpcg, &
       tolslk, tolsw1, tolsw2, sw1, sw2, zero )
      call pdnet_error(ierror)
      return
   end if

   !     ------------------------------------------------------------------
   !     Define the starting values.
   !     ------------------------------------------------------------------

   call pdnet_startvalues (av1, key, bound, db, dc, du, endn, gap, maxb, &
        maxc, minb, miu, na, nbound, nn, norb, npfrst, root, s, s1fctr, strn, &
        tolpcg, upfrst, unext, w, x, y, z)

   !     ------------------------------------------------------------------
   !     Initialize interior point iteration counter and define label for 
   !     main loop.
   !     ------------------------------------------------------------------
   
   pduit = 0
   info(19) = 3
   prttyp = 3
   do i = 1, nn
      yd(i) = 0.0d0
   end do

   do i = 1, na
      optflo(i) = 0.0d0
   end do
   
   do 
      
      !     ------------------------------------------------------------------
      !     Compute the right-hand-side for preconditioned conjugate gradient.
      !     ------------------------------------------------------------------

      call pdnet_comprhs(av1, av3, db, bound, bposit, dc, endn, larcs, lfirst, lsets, &
           miu, na, nbound, nn, npfrst, rhs, root, rsets, s, strn, theta, du, uarcs, &
           ufirst, upfrst, unext, w, wd, x, xd, y, z, zd)

      !     ------------------------------------------------------------------
      !     Identify maximum spanning tree and build data structures.
      !     ------------------------------------------------------------------

40    continue
      
      if (iprcnd == 1) then
         go to 50
      end if
      
      call pdnet_heap (bfirst, bposit, bucket,  endn, father, first, huge, ifirst,&
           innext, key, larcs, lfirst, list, link, lost, lsets, mark, mr, na, nn,&
           oufrst, ounext, pred, rank, rlink, root, rsets, rtrth, strn, theta, trearc,&
           treath, uarcs, ufirst, weight)

      !     ------------------------------------------------------------------
      !     Test optimality of the maximum spanning tree.
      !     ------------------------------------------------------------------
      
      call pdnet_optcheck (p, av4, b, bound, bstats, c, dof, endn, father, mark, na, &
           nbound, nn, list, npfrst, opt, optflo, pf, pof, root, rsets, rtrth, strn, &
           trearc, treath, u, uarcs, ufirst, upfrst, unext, lost, y, zz, zero)

      !     ------------------------------------------------------------------
      !     Test max flow stopping criterion if miu < 1.0d0.
      !     ------------------------------------------------------------------
      
      if (miu < 1.0d0) then
         call pdnet_checkmaxflow (av4, du, dc, endn, mfstat, mxfe, mxfv, nbuk, &
              na, nn, mfdopt, rc, s, stedg, stptol, stptr, stpval, strn, stvrtx,&
              db, tol2, tol1, wk2e1, wrkdv1, wrke1, wrke2, wrke3, wrke4, balance,&
              w, optflo, x, y, z)
         
         !        ---------------------------------------------------------------
         !        Update the tolerances used to guess the optimal partition.
         !        ---------------------------------------------------------------
         
         tol1=tol1*gtolmf
         tol2=tol2/gtolmf
         
         !        ---------------------------------------------------------------
         !        Exit on optimal maxflow.
         !        ---------------------------------------------------------------
         
         if (rc == 0) goto 70

      end if

50    continue

      if (pduit >= 1) then
         
         !        --------------------------------------------------------------
         !        Compute objective function.
         !        --------------------------------------------------------------

         call pdnet_compobjectivefun(db, bound, dc, dobj, fadd, na, nn, pobj, &
              root, du, upfrst, unext, w, x, y)

         !        --------------------------------------------------------------
         !        Print iteration summary.
         !        --------------------------------------------------------------

         call pdnet_printout(cgstop, oldcos, dobj,   dof,    gap,    iout, iprcnd,&
              pduit,  mfdopt, mfstat, miu,    mxfe, mxfv,   na,        opt,&
              pcgcos, pcgit, pcgres, oldtol, pf,     pobj,   pof,    prttyp, tol1,&
              tol2,   optflo)

      end if

      !     ------------------------------------------------------------------
      !     Preconditioner is determined by the value of iprcnd.
      !
      !          iprcnd=1 : Diagonal preconditioner.
      !          iprcnd=2 : IQRD preconditioner.
      !          iprcnd=3 : Spanning tree preconditioner.
      !          iprcnd=4 : Diagonally compensated spanning tree
      !                      preconditioner.
      !
      !     Compute preconditioner.
      !     ------------------------------------------------------------------

      if (iprcnd /= 3) then
         call pdnet_iqrd(av4, diag, endn, father, fcttol, iprcnd, na, ndiag,&
              nn, rho0, root, rtrth, strn, theta, trearc)

      end if

60    continue

      !     ------------------------------------------------------------------
      !     Compute the dual direction yd and the gap.
      !     ------------------------------------------------------------------

      if (pduit > 1) then
         rgap = dabs(pobj-dobj)/dabs(pobj)
      else
         rgap = 1
      end if

      !     ------------------------------------------------------------------
      !     Compute the tolerances for PCG stopping criterion.
      !     ------------------------------------------------------------------

      if (rgap > 1.0d1) then
         tolcg=tolcg0
      else
         tolcg=tolcg0*rgap
      end if

      !     ------------------------------------------------------------------
      !     Decide which preconditioner to use.
      !     ------------------------------------------------------------------

      if (pduit > 1 .and.  iprcnd == 4 .and. miu < 1.0d0) iprcnd = 3
      if (tolcg < mntlcg)  tolcg=mntlcg
      rho0=rho0*grho0
      if (rho0 < mrho0) rho0=mrho0

      !     ------------------------------------------------------------------
      !     Solve the system using preconditioned conjugate gradient method.
      !     ------------------------------------------------------------------

      call pdnet_precconjgrd(rhs, cgstop, diag, endn, pcgcos, father, iprcnd, &
           mxcgit, na, ndiag, nn, p, pcgit, q, r, pcgres, root, rtrth, strn,&
           sw1, theta, tolcg, trearc,treath, yd, zz)

      !     ------------------------------------------------------------------
      !     Check if change from diagonal preconditioner is needed.
      !     ------------------------------------------------------------------

      if ((iprcnd == 1) .and. (dble(pcgit) > sw1)) then
         
         !        ---------------------------------------------------------------
         !        Compute objective function.
         !        ---------------------------------------------------------------

         call pdnet_compobjectivefun(db, bound, dc, dobj, fadd, na, nn, pobj, root,&
              du, upfrst, unext, w, x, y)
        
         !        ---------------------------------------------------------------
         !        Print summary of iteration.
         !        ---------------------------------------------------------------

         call pdnet_printout(cgstop, oldcos, dobj, dof, gap, iout, iprcnd, pduit,&
              mfdopt, mfstat, miu,    mxfe, mxfv, na,  opt, pcgcos, pcgit, &
              pcgres, oldtol, pf, pobj, pof, prttyp, tol1, tol2, optflo)

         !        ---------------------------------------------------------------
         !        Switch to maximum spanning tree preconditioner.
         !        ---------------------------------------------------------------

         iprcnd = 3

         !        ---------------------------------------------------------------
         !        Ignore current iteration using diagonal preconditioner.
         !        Redo it using maximum spanning tree preconditioner.
         !        ---------------------------------------------------------------

         go to 40
      
      end if

      !     ------------------------------------------------------------------
      !     Current iteration accepted. Update primal and dual solutions.
      !     ------------------------------------------------------------------
      
      oldtol=tolpcg

      call pdnet_updatesol(av1, key, av3, db, bound,&
           ds, endn, factor, gap, huge, miu, na, nbound, nn, npfrst, ps, &
           root, s, sd, strn, theta, tolpcg, du, upfrst, unext, w, wd, x, xd,&
           y, yd, z, zd)

      !     ------------------------------------------------------------------
      !     Increment iteration counter and save PCG cosine tolerance.
      !     ------------------------------------------------------------------
      
      pduit=pduit+1
      oldcos=tolcg

      !     ------------------------------------------------------------------
      !     Check iteration counter against maximum number of iterations.
      !     ------------------------------------------------------------------ 
      
      if (pduit > maxit) then

         !        ---------------------------------------------------------------
         !        Maximum number of iterations reached.
         !        Set error flag and go to synchronization before returning.
         !        ---------------------------------------------------------------

         ierror = 6
         goto 80
      end if

      !     ------------------------------------------------------------------
      !     Maximum number of iterations not reached yet, Choose the
      !     preconditioner.
      !     ------------------------------------------------------------------
      
      if (iprcnd == 1 .and. pcgit > sw1 .and. pcgit < sw2) then
         iprcnd = 3
      elseif (iprcnd == 0 .and. pcgit < sw2) then
         iprcnd = 3
      elseif (pcgit > sw2) then
         iprcnd = 0
      endif

      !     ------------------------------------------------------------------
      !     End of primal-dual iterative loop.
      !     ------------------------------------------------------------------

   end do

70 continue

   !     ------------------------------------------------------------------
   !     Compute objective value at current iteration.
   !     ------------------------------------------------------------------

   call pdnet_compobjectivefun(db, bound, dc, dobj, fadd, na, nn, pobj, root,&
        du, upfrst, unext, w, x, y)

   !     -----------------------------------------------------------------
   !     Synchronization: set all the info values in vector info.
   !     -----------------------------------------------------------------

80 call pdnet_putinfo(bound,  cgstop, dinfo, dobj,  dof, factor, fcttol, &
        gap,    grho0,  gtolcg, gtolmf, huge, ierror, info, ioflag, iout,&
        iprcnd, maxit,  mfdopt, mfstat, miu, mntlcg, mr,  mrho0,  mxcgit,& 
        mxfe,   mxfv, nbound, nbuk, oldcos, oldtol,  opt, optgap, optmf, &
        optst, pcgcos, pcgit, pcgres, pduit, pf,   pobj,  pof,  prttyp,  &
        rho0, root, s1fctr, stptol, stpval, tol1, tol2,   tolcg,  tolcg0,&
        tolpcg, tolslk, tolsw1, tolsw2, sw1, sw2, zero)

   call pdnet_error(ierror)

   !     ------------------------------------------------------------------
   !     End of subroutine pdnet.
   !     ------------------------------------------------------------------

 end subroutine pdnet
!     ------------------------------------------------------------------
!     pdnet_approxcosine: 
!             This subroutine computes approximate cosine of directions.
!     ------------------------------------------------------------------

subroutine pdnet_approxcosine(b, dcosn, nn, nrmb, q, r, root)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              nn    : order of the vectors.
  !              root  : beginning of data structure.
  !     ------------------------------------------------------------------
  
  integer ::  nn    , root

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !              b(1:nn): first vector of the angle.
  !              q(1:nn): subtraction of the two vectors of the angle.
  !              r(1:nn): second vector of the angle.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: b, q, r

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              i     : loop counter.
  !     ------------------------------------------------------------------
 
  integer ::  i

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !              dcosn : cosine between b and r.
  !              nrmb  : 2-norm of b.
  !     ------------------------------------------------------------------
  
  double precision ::  dcosn , nrmb

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !              dnrm2 : 2-norm (BLAS function).
  !              eps   : inner product of b and q.
  !              nrmq  : 2-norm of q.
  !     ------------------------------------------------------------------
  
  double precision ::  dnrm2 , eps   , nrmq

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  do i=1,nn
     q(i)=b(i)-r(i)
  end do

  nrmq= dnrm2 ( nn, q, 1)

  call pdnet_ddot(eps,    nn,     root,   b,      q)
  
  dcosn = dabs(1.0d0-eps/(nrmb*nrmq))

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_approxcosine.
  !     ------------------------------------------------------------------

end subroutine pdnet_approxcosine
!     ------------------------------------------------------------------
!     pdnet_attx: This subroutine computes the product y=a*(1/teta)*a'x.
!     ------------------------------------------------------------------

subroutine pdnet_attx(endn, na, nn, root, strn, teta, x, y)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------

  integer :: na    , nn    , root

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endn(1:na): pointer of network data structure.
  !             strn(1:na): adjacency list of network data structure.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: endn, strn

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !            teta(1:na): normal equations diagonal matrix.
  !            x   (1:nn): primal iterate x.
  !            y   (1:nn): primal iterate y.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: x, y
  double precision, dimension(1:na) :: teta

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             arc   : current arc.
  !             enod  : current element on endn.
  !             snod  : current element on strn.
  !             i     : counter variable
  !     ------------------------------------------------------------------

  integer :: arc   , enod  , snod, i

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             aux1  : auxiliary variable.
  !     ------------------------------------------------------------------

  double precision :: aux1

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  do i = 1, nn
     y(i) = 0.0d0
  end do

  do arc=1,na
     snod = strn(arc)
     enod = endn(arc)
     aux1 = (x(enod)-x(snod))/teta(arc)
     y(snod) = y(snod)-aux1
     y(enod) = y(enod)+aux1
  end do

  y(root) = 0.0d0

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_attx.
  !     ------------------------------------------------------------------

end subroutine pdnet_attx
!     ------------------------------------------------------------------
!     pdnet_backsub: 
!             Given a spanning tree defined by a list of vertices in 
!             dfs order and parent edges, it computes the solution the
!             system s'x = b (back substitution).  s is a block diagonal
!             incidence matrix with one linear dependent row per block.
!             The rhs is mapped to each vertex, as is s were a square
!             matrix with a "zero" edge added to the beginning of each 
!             block.  Solution is returned on rhs.
!     ------------------------------------------------------------------

subroutine pdnet_backsub(ntree,  nedge,  nvrtx,  strv, endv, stptr, stvrtx,  &
     stedg,  rhs)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             ntree : number of trees in spanning forest.
  !             nedge : number of edges.
  !             nvrtx : number of vertices.
  !     ------------------------------------------------------------------
  
  integer :: ntree , nedge , nvrtx

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             strv  (1:  nedge): start vertex in network data structure.
  !             stptr (1:    ntree + 1): pointer to first vertex in a tree
  !             stvrtx(1:    nvrtx    ): dfs ordering of vertices.
  !             stedg (1:    nedge    ): parent edge for each vertex.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: strv, endv, stedg
  integer, dimension(1:nvrtx) :: stvrtx
  integer, dimension(1:ntree+1) ::  stptr 

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             rhs   (1:    nvrtx    ): rhs and solution of system.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nvrtx) :: rhs

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             edg   : edge index.
  !             tree  : tree counter.
  !             ptr   : temporary pointer.
  !             v1    : vertex index.
  !             v2    : adjacent vertex.
  !             vrtx  : vertex index.
  !     ------------------------------------------------------------------
  
  integer :: edg   , tree  , ptr   , v1    , v2    , vrtx

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Process each tree.
  !     ------------------------------------------------------------------

  do tree = 1, ntree
     
     !        ---------------------------------------------------------------
     !        Set variable corresponding to root to zero.
     !        ---------------------------------------------------------------
     
     rhs(stvrtx(stptr(tree)))=0.d0

     !        ---------------------------------------------------------------
     !        Process remainder of tree.
     !        ---------------------------------------------------------------

     do ptr=stptr(tree)+1,stptr(tree+1)-1
        
        !           ------------------------------------------------------------
        !           Identify vertices and edge.
        !           ------------------------------------------------------------
        
        vrtx=stvrtx(ptr)
        edg=stedg(vrtx)
        v1=strv(edg)
        v2=endv(edg)

        !           ------------------------------------------------------------
        !           Solve for current edge.
        !           ------------------------------------------------------------
        
        if(vrtx == v1) then
           rhs(v1)=rhs(v1)+rhs(v2)
        else
           rhs(v2)=-rhs(v2)+rhs(v1)
        end if
     end do
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_backsub.
  !     ------------------------------------------------------------------

end subroutine pdnet_backsub
!     ------------------------------------------------------------------
!     pdnet_buildstruct: 
!             This subroutine builds a data structure for accessing each
!             individual tree in a spanning forest from subroutine
!             pdmxst. Vertices are listed in ascending order of indices
!             in each tree list.
!     ------------------------------------------------------------------

subroutine pdnet_buildstruct(ntree,  nvrtx,  parent, stptr,  stvrtx)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             ntree : number of trees in spanning forest.
  !             nvrtx : number of vertices.
  !     ------------------------------------------------------------------
  
  integer :: ntree , nvrtx

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             parent(1:nvrtx  ): parent array from pdmxst.
  !             stptr (1:ntree+1): pointer to first vertex in a tree.
  !             stvrtx(1:nvrtx  ): list of vertices in spanning forest.
  !     ------------------------------------------------------------------

  integer, dimension(1:nvrtx)   ::  parent, stvrtx
  integer, dimension(1:ntree+1) ::  stptr

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             tree  : tree index.
  !             vrtx  : vertex index.
  !     ------------------------------------------------------------------
 
  integer :: tree  , vrtx

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Initialize pointers.
  !     ------------------------------------------------------------------

  do tree = 1, ntree
     stptr(tree)=0
  end do

  !     ------------------------------------------------------------------
  !     Count number of vertices in each tree.
  !     ------------------------------------------------------------------
  
  do vrtx = 1, nvrtx
     tree=parent(vrtx)
     stptr(tree)=stptr(tree)+1
  end do

  !     ------------------------------------------------------------------
  !     Compute pointers.
  !     ------------------------------------------------------------------ 

  stptr(ntree+1)=nvrtx-stptr(ntree)+1
  
  do tree = ntree, 2, -1
     stptr(tree)=stptr(tree+1)-stptr(tree-1)
  end do
  stptr(1)=1

  !     ------------------------------------------------------------------
  !     Build data structure holding spanning forest.
  !     ------------------------------------------------------------------
  
  do vrtx = 1, nvrtx
     tree=parent(vrtx)+1
     stvrtx(stptr(tree))=vrtx
     stptr(tree)=stptr(tree)+1
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_buildstruct.
  !     ------------------------------------------------------------------

end subroutine pdnet_buildstruct
!     ------------------------------------------------------------------
!     pdnet_check: This subroutine checks validity of parameters.
!     ------------------------------------------------------------------


subroutine pdnet_check (endn, ierror, na, nn, prttyp)

  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             ierror: flags errors in input.
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             prttyp: printing level.
  !     ------------------------------------------------------------------
  
  integer :: ierror, na, nn, prttyp
  integer :: ipos

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endn  (1:  na): pointer of network data structure.
  !     ------------------------------------------------------------------
  integer, dimension (1:na) :: endn

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Set error condition for no error.
  !     ------------------------------------------------------------------
  
  ierror = 0

  !     ------------------------------------------------------------------
  !     Check if the number of each node does not exceed the number of
  !     nodes.
  !     ------------------------------------------------------------------
  ipos = 1
  do while (ipos <= na)
     if (endn(ipos) > nn) then
        ierror = 3
        return
     end if
     ipos = ipos + 1
  end do

  !     ------------------------------------------------------------------
  !     Check if output level was correctly set.
  !     ------------------------------------------------------------------

  if (prttyp < 0 .or. prttyp > 3) then
     ierror = 5
  end if

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_check.
  !     ------------------------------------------------------------------

end subroutine pdnet_check
!     ---------------------------------------------------------------------------
!     pdnet_checkfeas : This subroutine checks if the network has enough capacity 
!                       for the proposed flow.
!     ---------------------------------------------------------------------------

subroutine pdnet_checkfeas( b, cap, endv,  ierror, nedge,  nvrtx,strv, wrke1, &
     balance)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             ierror: error code returned (0=no error).
  !             nedge : number of arcs.
  !             nvrtx : number of nodes.
  !     ------------------------------------------------------------------

  integer :: ierror, nedge, nvrtx

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             strv  (1:  nedge): start vertex in network data structure.
  !             wrke1 (1:  nedge): temporary workspace.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: endv, strv, wrke1, balance

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             b     (1:  nvrtx): array of supply at each node.
  !             cap   (1:  nedge): upper bound of flow at each edge
  !     ------------------------------------------------------------------

  double precision, dimension(1:nvrtx) :: b
  double precision, dimension(1:nedge) :: cap

  !     ------------------------------------------------------------------
  !     Integer variables
  !
  !             i     : loop counter.
  !             mflow : maximum flow supplied.
  !             mfreq : maximum flow requested.
  !             mxfe  : edges on optimal flow network.
  !             mxfv  : vertices on optimal flow network.
  !             nactmf: number of edges in active list.
  !             ncapmf: number of edges at capacity.
  !     ------------------------------------------------------------------
  
  integer   ::   i, mflow , mfreq , mxfe  , mxfv  , nactmf, ncapmf
  
  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             adem  : overall supply.
  !             asup  : overall demand.
  !             bi    : capacity at node i.
  !     ------------------------------------------------------------------

  double precision :: adem  , asup  , bi
  
  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------
  
  ierror = 0

  !     ------------------------------------------------------------------
  !     Check if the overall supply is equal to the overall demand.
  !     ------------------------------------------------------------------
  asup = 0.0d0
  adem = 0.0d0  
  do i = 1, nvrtx
     bi = b(i)
     if (bi < 0.0d0) then
        adem = adem - bi
     end if
     if (bi > 0.0d0) then
        asup = asup + bi
     end if
  end do

  if (asup /= adem) then
     ierror = 6
     return
  end if

  !     ------------------------------------------------------------------
  !     Place all edges in the active list.
  !     ------------------------------------------------------------------

  nactmf = 0
  ncapmf = 0 

  do i = 1, nedge
     nactmf = nactmf + 1
     wrke1(nactmf) = nedge - i + 1     
  end do

  !     ------------------------------------------------------------------
  !     Solve maxflow problem.
  !     ------------------------------------------------------------------  

  call pdnet_feasible (wrke1, cap, endv, mflow, mfreq, mxfe, mxfv, nactmf,&
       ncapmf, nedge, nvrtx, strv, b, balance)
       

  if (mflow /= int(asup)) then
     ierror = 7
  end if

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_checkfeas.
  !     ------------------------------------------------------------------

end subroutine pdnet_checkfeas

  
!     ------------------------------------------------------------------
!     pdnet_checkmaxflow: 
!             This subroutine checks optimality of flow by using max
!             flow criterion:
!
!             1. Identify active, upper bound and lower bound edges.
!
!             2. Compute spanning tree (maximal forest) for active edges.
!
!             3. Build tentative dual optimal solution.
!
!             4. Identify dual face.
!
!             5. Solve max flow for edges defining optimal face.
!     ------------------------------------------------------------------ 

subroutine pdnet_checkmaxflow (av4, cap, cost,  endv,  mfstat, mxfe, mxfv, &
     nbuk,   nedge,  nvrtx,  optdl,  rc,ssol,   stedg,  stptol, stptr,  stpval, &
     strv, stvrtx, supply, tolhi,  tollo,  wrk2e1, wrkdv1, wrke1,  wrke2,  wrke3,&
     wrke4,  balance, wsol, xopt,   xsol,   ysol,   zsol)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             mxfe  : edges on optimal flow network.
  !             mxfv  : vertices on optimal flow network.
  !             nbuk  : number of buckets.
  !             nedge : number of arcs.
  !             nvrtx : number of nodes.
  !             rc    : flags optimality status.
  !     ------------------------------------------------------------------
  
  integer :: mxfe  , mxfv  , nbuk  , nedge , nvrtx , rc

  !     ------------------------------------------------------------------
  !     Integer arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             stedg (1:  nedge): edges on spanning tree.
  !             stptr (1:  nedge): pointer for edges on spanning tree.
  !             stvrtx(1:  nedge): vertices on the spanning forest.
  !             strv  (1:  nedge): start vertex in network data structure.
  !             wrk2e1(1:2*nedge): working array of order 2*nedge.
  !             wrke1 (1:  nedge): working array of order nedge.
  !             wrke2 (1:  nedge): working array of order nedge.
  !             wrke3 (1:  nedge): working array of order nedge.
  !             wrke4 (1:  nedge): working array of order nedge.
  !             xopt  (1:  nedge): optimal solution of maxflow problem.
  !     ------------------------------------------------------------------

  integer, dimension(1:2*nedge) :: wrk2e1
  integer, dimension(1:  nedge) ::  endv, stedg, stptr, stvrtx, strv, wrke1, &
       wrke2, wrke3, wrke4, balance, xopt

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !              av4   (1:nedge): auxiliary vector.
  !              cap   (1:nedge): capacity at each arc.
  !              cost (1:nedge): cost at each node.
  !              ssol  (1:nedge): solution for s.
  !              supply(1:nedge): supply at each node.
  !              wrkdv1(1:nvrtx): work array. 
  !              wsol  (1:nedge): solution for w.
  !              xsol  (1:nedge): solution for x.
  !              ysol  (1:nvrtx): solution for y.
  !              zsol  (1:nedge): solution for z.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nvrtx) ::  wrkdv1, ysol
  double precision, dimension(1:nedge) ::  supply, cost, av4, cap, ssol, wsol, xsol, zsol

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             mflow : maximum flow supplied.
  !             mfreq : maximum flow requested.
  !             mfstat: optimality status.
  !             nact  : number of active edges.
  !             nactmf: number of edges in active list.
  !             ncap  : number of edges at the upper bound.
  !             ncapmf: number of edges at capacity.
  !             ntree : number of trees in spanning forest.
  !     ------------------------------------------------------------------
  
  integer :: i , mflow , mfreq , mfstat, nact  , nactmf, ncap  , ncapmf, ntree

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             optdl : dual objective value.
  !             stptol: zero for dual slacks.
  !             stpval: largest slack on dual face.
  !             tolhi : upper bound of the tolerance for maxflow.
  !             tollo : lower bound of the tolerance for maxflow.
  !     ------------------------------------------------------------------

  double precision :: optdl , stptol, stpval, tolhi , tollo
  
  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------
  !
  !
  !     ------------------------------------------------------------------
  !     Select active edges.
  !     ------------------------------------------------------------------

  call pdnet_tapia (wrke1,  wrke2,  nact,   ncap,   nedge,  ssol, tolhi,&
       tollo,  wsol,   xsol,   zsol )

  !     ------------------------------------------------------------------
  !     Sort active edges for spanning tree.
  !     ------------------------------------------------------------------

  do i = 1, nedge
     av4(i)=zsol(i)/xsol(i) + wsol(i)/ssol(i)
  end do
  call pdnet_sortedges(wrke1,  wrke4,  nedge, nact,   nbuk,   wrke2,  wrke3, av4)

  !     ------------------------------------------------------------------
  !     Identify spanning tree for active edges.
  !     ------------------------------------------------------------------

  call pdnet_kruskal(endv, nact, ntree, nvrtx, nedge, wrke2,  wrke3, &
       wrke4, stedg, strv)
  
  !     ------------------------------------------------------------------
  !     Build data structure for spanning forest.
  !     ------------------------------------------------------------------

  call pdnet_buildstruct (ntree,  nvrtx,  wrke3,  stptr,  stvrtx)
  
  !     ------------------------------------------------------------------
  !     Order spanning forest.
  !     ------------------------------------------------------------------
  
  call pdnet_orderstruct(endv,  wrk2e1, wrke3,  nedge,  ntree,  nvrtx, &
       wrke4,  stedg,  stptr,  strv, stvrtx)

  !     ------------------------------------------------------------------
  !     Build tentative dual optimal.
  !     ------------------------------------------------------------------

  call pdnet_dualoptimal (cost,  endv,  nedge,  ntree,  nvrtx,  stedg, &
       stptr,  strv, stvrtx, wrkdv1, ysol)

  !     ------------------------------------------------------------------
  !     Identify optimal dual face.
  !     ------------------------------------------------------------------

  call pdnet_optimaldualface (wrke1,  cost,   endv,  stpval, nactmf, &
       ncapmf, nedge,  nvrtx,   strv,  stptol, wrkdv1)

  !     ------------------------------------------------------------------
  !     Compute dual objective value.
  !     ------------------------------------------------------------------

  call pdnet_dualobjective (cap,  cost,  endv,  nedge,  nvrtx,  supply, &
       strv, wrkdv1, optdl)
  
  !     ------------------------------------------------------------------
  !     Solve maxflo problem.
  !     ------------------------------------------------------------------
  
  call pdnet_maxflow ( wrke1,  cap,      endv,  mflow,  mfreq,  mxfe,  mxfv,&
       nactmf, ncapmf,  nedge,  nvrtx,  strv, supply, balance, xopt)

  !     ------------------------------------------------------------------
  !     Return optimality status.
  !     ------------------------------------------------------------------
    
  if(mfreq == mflow) then
     rc = 0
     mfstat = 2
  else
     rc = 1
     mfstat = 1
  end if

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_checkmaxflow.
  !     ------------------------------------------------------------------

end subroutine pdnet_checkmaxflow
!     ------------------------------------------------------------------
!     pdnet_compobjectivefun: 
!             This subroutine computes the current value of the
!             objective function.
!     ------------------------------------------------------------------

subroutine pdnet_compobjectivefun(b,      bound,  c,      dobj,   fadd,&
     na, nn,     pobj,   root,   u, upfrst, upnext, w,      x,      y)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              bound : number of bounded arcs.
  !              na    : number of arcs.
  !              nn    : number of nodes.
  !              root  : sets beginning of structure.
  !              upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------
  
  integer  :: bound , na    , nn    , root  , upfrst
  
  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: upnext

  !     ------------------------------------------------------------------
  !     Double precision arrays:
  !             b(1:nn): array of capacities at each node.
  !             c(1:nn): array of costs at each arc.
  !             u(1:na): upper bound of flow at each arc.
  !             w(1:na): dual iterate w.
  !             x(1:na): primal iterate x.
  !             y(1:nn): dual iterate y.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: b, y
  double precision, dimension(1:na) :: c, u, w, x

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !              arc   : current arc.
  !              i     : loop counter.
  !              node  : another loop counter.
  !     ------------------------------------------------------------------
    
  integer  :: arc   , i     , node

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !              dobj  : dual objective function value.
  !              fadd  : auxiliary variable for cost computation.
  !              pobj  : objective function value.
  !     ------------------------------------------------------------------
  
  double precision ::  dobj  , fadd  ,  pobj

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------
  
  dobj = 0d0
  pobj = 0d0
  
  do node=1,root-1
     dobj = dobj+b(node)*y(node)
  end do

  do node=root+1,nn
     dobj = dobj+b(node)*y(node)
  end do

  arc = upfrst

  do i=1,bound
     dobj = dobj-w(arc)*u(arc)
     arc = upnext(arc)
  end do

  do arc=1,na
     pobj = pobj+c(arc)*x(arc)
  end do

  pobj = pobj+fadd
  dobj = dobj+fadd

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_compobjectivefunction.
  !     ------------------------------------------------------------------

end subroutine pdnet_compobjectivefun

!     ----------------------------------------------------------------------
!     pdnet_comprhs : This subroutine computes the right-hand side of the 
!                     normal equations system to be solved in the network 
!                     primal-dual interior-point algorithm.  This subroutine 
!                     predefines the sets l and u.
!     ----------------------------------------------------------------------

subroutine pdnet_comprhs (av1,    av3,    b,      bound,  bposit, c, endn, &
     larcs,  lfirst, lsets,  miu,    na, nbound, nn,     npfrst, rhs, root,&
     rsets,   s,      strn,   teta,   u, uarcs,  ufirst, upfrst, upnext, w,&
     wd,     x,      xd,  y,      z,      zd )

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             bound : number of bounded arcs.
  !             larcs : number of arcs at the lower bound.
  !             lfirst: first arc at the lower bound.
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             nbound: number of unbounded arcs.
  !             npfrst: next bounded arc to be considered.
  !             root  : sets beginning of structure.
  !             uarcs : number of arcs at the upper bound.
  !             ufirst: first arc at the upper bound.
  !             upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------

  integer :: bound, larcs, lfirst, na, nn, nbound, npfrst, root, uarcs,&
       ufirst, upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !          bposit(1:na): flags arcs according to bounds.
  !          endn  (1:na): pointer of network data structure.
  !          lsets (1:na): list of arcs at lower bound.
  !          rsets (1:na): list of arcs at upper bound.
  !          strn  (1:na): adjacency list of network data structure.
  !          upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: bposit, endn, lsets, rsets, strn, upnext

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !          av1 (1:na): auxilliary vector 1.
  !          av3 (1:na): auxilliary vector 3.
  !          b   (1:nn): array of capacities at each node.
  !          c   (1:na): array of costs at each arc.
  !          rhs (1:nn): right-hand side of the linear system.
  !          s   (1:na): dual iterate s.
  !          teta(1:na): normal equations diagonal matrix.
  !          u   (1:na): upper bound of flow at each arc.
  !          w   (1:na): dual iterate w.
  !          wd  (1:na): direction dw.
  !          x   (1:na): primal iterate x.
  !          xd  (1:na): direction dx.
  !          y   (1:nn): primal iterate y.
  !          z   (1:na): primal iterate z.
  !          zd  (1:na): direction dz.
  !     -----------------------------------------------------------------

  double precision, dimension(1:nn) :: b, rhs, y
  double precision, dimension(1:na) :: av1, av3, c, s, teta, u, w, wd, x, xd, &
       z, zd

  !     ------------------------------------------------------------------
  !     Integer variables.
  !
  !             arc   : current arc.
  !             enod  : current element on endn.
  !             i     : loop counter.
  !             node  : another loop counter.
  !             luxarc: current arc on lower bound.
  !             snod  : current element on strn.
  !             uuxarc: current arc on upper bound.
  !     ------------------------------------------------------------------

  integer :: arc   , enod  , i     , node  , luxarc, snod  , uuxarc

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             miu: parameter mu.
  !     ------------------------------------------------------------------

  double precision :: miu

  !     ------------------------------------------------------------------
  !     Double precision variables.
  !
  !            aux1-7: auxilliary variables.
  !     ------------------------------------------------------------------

  double precision :: aux1  , aux2  , aux3  , aux4  ,  aux5  , aux6  , aux7

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     
  !     Setting lfirst and ufirst equal to zero.
  !     ------------------------------------------------------------------

  lfirst = 0
  ufirst = 0
  larcs  = 0
  uarcs  = 0

  !     ------------------------------------------------------------------
  !     Initialize the right-hand side vector rhs.
  !     ------------------------------------------------------------------

  do node = 1, root - 1
     rhs(node) = b(node)
  end do

  do node = root + 1, nn
     rhs(node) = b(node)
  end do

  rhs(root) = 0.0d0

  !     ------------------------------------------------------------------
  !     Compute the right-hand side for pcg.
  !     Define sets b, l and u.
  !     
  !     
  !     Arcs without upper bound.
  !     ------------------------------------------------------------------

  arc = npfrst
  
  do i = 1, nbound
     
     !        ---------------------------------------------------------------
     !        Compute working parameters.
     !        ---------------------------------------------------------------  

     snod = strn(arc)
     enod = endn(arc)
     aux1 = z(arc)/x(arc)
     aux3 = miu/x(arc)
     aux6 = (aux3-y(snod)+y(enod)-c(arc))/aux1
     aux7 = aux6+x(arc)

     !        ---------------------------------------------------------------
     !        Compute the right-hand side for pcg.
     !        ---------------------------------------------------------------  

     rhs(snod) = rhs(snod)+aux7
     rhs(enod) = rhs(enod)-aux7     

     !        ---------------------------------------------------------------
     !        Predefine set l.
     !        ---------------------------------------------------------------

     if (lfirst == 0) then
        lfirst         = arc
        lsets(lfirst)  = 0
        rsets(lfirst)  = 0
        luxarc        = lfirst
     else
        lsets(arc)     = luxarc
        rsets(arc)     = 0
        rsets(luxarc) = arc
        luxarc        = arc
     end if

     bposit(arc) = -1
     larcs = larcs+1    

     !        ---------------------------------------------------------------
     !        Compute auxilliary vectors.
     !        ---------------------------------------------------------------

     xd(arc)   = aux6
     zd(arc)   = aux3
     av1(arc)  = aux1
     teta(arc) = aux1

     !        ---------------------------------------------------------------
     !        Compute the next arc.
     !        ---------------------------------------------------------------
     
     arc = upnext(arc)     
  end do

  !     ------------------------------------------------------------------
  !     Arcs with upper bound
  !     ------------------------------------------------------------------
    
  arc = upfrst

  do i = 1, bound

     !        ---------------------------------------------------------------
     !        Compute working parameters.
     !        ---------------------------------------------------------------     

     snod = strn(arc)
     enod = endn(arc)
     aux1 = z(arc)/x(arc)
     aux2 = w(arc)/s(arc)
     aux3 = miu/x(arc)
     aux4 = miu/s(arc)
     aux5 = aux1+aux2
     aux6 = (aux3-aux4-y(snod)+y(enod)-c(arc)+ aux2*(u(arc)-x(arc))-w(arc))/aux5
     aux7 = aux6+x(arc)

     !        ---------------------------------------------------------------
     !        Compute the right-hand side for pcg.
     !        ---------------------------------------------------------------   

     rhs(snod) = rhs(snod)+aux7
     rhs(enod) = rhs(enod)-aux7    

     !        ---------------------------------------------------------------
     !        Predefine sets l and u.
     !        ---------------------------------------------------------------

     if (aux1 > aux2) then
        if (lfirst == 0) then
           lfirst         = arc
           lsets(lfirst)  = 0
           rsets(lfirst)  = 0
           luxarc         = lfirst    
        else
           lsets(arc)     = luxarc
           rsets(arc)     = 0
           rsets(luxarc) = arc
           luxarc        = arc  
        end if
        larcs             = larcs+1
        bposit(arc)       = -1
     else
        if (ufirst == 0) then
           ufirst         = arc
           lsets(ufirst)  = 0
           rsets(ufirst)  = 0
           uuxarc        = ufirst          
        else
           lsets(arc)     = uuxarc
           rsets(arc)     = 0
           rsets(uuxarc) = arc
           uuxarc        = arc     
        end if
        uarcs             = uarcs+1
        bposit(arc)       = 1
     end if

     !        ---------------------------------------------------------------
     !        Compute auxilliary vectors.
     !        ---------------------------------------------------------------     

     xd(arc)   = aux6
     zd(arc)   = aux3
     wd(arc)   = aux4
     av1(arc)  = aux1
     av3(arc)  = aux2
     teta(arc) = aux5   

     !        ---------------------------------------------------------------
     !        Compute the next arc.
     !        ---------------------------------------------------------------
     
     arc = upnext(arc)

  end do

  !     ------------------------------------------------------------------
  !     Setting rhs(root) and av2(root) equal to zero.
  !     ------------------------------------------------------------------  
 
  rhs(root) = 0.0d0

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_comprhs.
  !     ------------------------------------------------------------------

end subroutine pdnet_comprhs
!     -------------------------------------------------------------------------
!     pdnet_datastruct: This subroutine constructs an additional data structure
!                       used by the primal-dual interior-point algorithm.
!     -------------------------------------------------------------------------


subroutine pdnet_datastruct( bound, bstats, endn, inaux, infrst, innext, na,&
     nbound, nn, npfrst, ouaux,  oufrst, ounext, strn,   u, upfrst, upnext)

  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables.
  !
  !             bound : number of bounded arcs.
  !             na    : number of arcs.
  !             nbound: number of unbounded arcs.
  !             nn    : number of nodes.
  !             npfrst: First bounded arc in upnext list
  !             upfrst: First unbounded arc in upnext list
  !     ------------------------------------------------------------------

  integer ::  bound , na    , nbound,  nn    , npfrst, upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             bstats(1:na): flags if each arc has an upper bound or not.
  !             endn  (1:na): pointer of network data structure.
  !             inaux (1:nn): pointer for a in-adjc network linked list.
  !             infrst(1:nn): points to the first node in in-adcncy list.
  !             innext(1:na): points to the last node in in-adcncy list.
  !             ouaux (1:nn): pointer for a out-adjc network linked list.
  !             oufrst(1:nn): points to the first node in out-adcncy list.
  !             ounext(1:na): points to the last node in out-adcncy list.
  !             strn  (1:na): adjacency list of network data structure.
  !             u     (1:na): upper bound of flow at each arc.
  !             upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: bstats, endn, innext, ounext, strn, u, upnext
  integer, dimension(1:nn) :: inaux, infrst, ouaux, oufrst
  
  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             arc   : current arc.
  !             auxnup: current unbounded position on data stucture.
  !             auxup : current bounded position on data stucture.
  !             enod  : current element on endn.
  !             i     : loop counter.
  !             snod  : current element on strn.
  !     ------------------------------------------------------------------
  
  integer :: arc   , auxnup, auxup , enod  , i     , snod

  !     ------------------------------------------------------------------
  !     Start of executable section of the subroutine.
  !     ------------------------------------------------------------------
  bound  = 0
  upfrst = 0
  nbound = 0
  npfrst = 0
  do i = 1, nn
     infrst(i) = 0
     oufrst(i) = 0
  end do
  do i = 1, na
     innext(i) = 0
     ounext(i) = 0
     upnext(i) = 0
  end do

  do arc = 1, na
     snod = strn(arc)
     enod = endn(arc)   
     if (infrst(enod) == 0) then
        infrst(enod)       = arc
        inaux(enod)        = arc
     else
        innext(inaux(enod)) = arc
        inaux(enod)         = arc
     end if
     if (oufrst(snod) == 0) then
        oufrst(snod)       = arc
        ouaux(snod)        = arc
     else
        ounext(ouaux(snod)) = arc
        ouaux(snod)         = arc
     end if
     if (u(arc) > 0) then
        if (upfrst == 0) then
           upfrst          = arc
           auxup           = arc
        else
           upnext(auxup)   = arc
           auxup    = arc
        end if
        bound = bound+1
        bstats(arc)        = 1
     else
        if (npfrst == 0) then
           npfrst           = arc
           auxnup           = arc
        else
           upnext(auxnup)   = arc
           auxnup           = arc
        endif
        nbound    = nbound+1
        bstats(arc)        = 0
     end if
  end do
  
  !---------------------------------------
  ! End of subroutine pdnet_datastruct.
  !---------------------------------------
end subroutine pdnet_datastruct
  
!     ------------------------------------------------------------------
!     pdnet_ddot: 
!             This subroutine computes the scalar product alfa=x'y.
!             The vectors x and y have both dimension dim.  This 
!             subroutine is based on the dblas ddot routine.
!     ------------------------------------------------------------------

subroutine pdnet_ddot (alfa,   dim,    root,   x,      y)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             dim   : dimension of vectors.
  !             root  : beginning of data structure.
  !     ------------------------------------------------------------------
  
  integer :: dim   , root

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             x(1:dim): vector x.
  !             y(1:dim): vector y.
  !     ------------------------------------------------------------------

  double precision, dimension(1:dim) :: x, y

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             m     : size of unrolled loop.
  !             mp1   : m+1.
  !     ------------------------------------------------------------------
  
  integer :: i     , m     , mp1

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             alfa  : inner product.
  !             dtemp : temporary variable.
  !     ------------------------------------------------------------------

  double precision :: alfa  , dtemp

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  dtemp = 0d0
  alfa  = 0d0
  x(root) = 0d0
  y(root) = 0d0
  m = mod(dim,5)

  if (m /= 0) then
     do i = 1,m
        dtemp = dtemp + x(i)*y(i)
     end do
     if (dim < 5) then
        alfa = dtemp
        return
     end if
  end if

  mp1 = m + 1

  do i = mp1,dim,5
     dtemp = dtemp + x(i)*y(i) + x(i + 1)*y(i + 1) + &
          x(i + 2)*y(i + 2) + x(i + 3)*y(i + 3)&
          + x(i + 4)*y(i + 4)
  end do

  alfa = dtemp

!     ------------------------------------------------------------------
!     End of subroutine pdnet_ddot.
!     ------------------------------------------------------------------

end subroutine pdnet_ddot
!     ------------------------------------------------------------------
!     pdnet_dualobjective: 
!           This subroutine computes dual objective value for a given
!             solution of an LP with upper bounds and 0 lower bounds.
!     ------------------------------------------------------------------

subroutine pdnet_dualobjective (cap,    cost,   endv,  nedge,  nvrtx, &
     supply, strv,   ysol, pddlob)

  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nedge : number of edges.
  !             nvrtx : number of vertices.
  !     ------------------------------------------------------------------
     
  integer :: nedge , nvrtx

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             edg   : edge index.
  !     ------------------------------------------------------------------
  
  integer :: edg

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             strv  (1:  nedge): start vertex in network data structure.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: endv, strv

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             cost  (1:nedge): cost coefficients.
  !             supply(1:nvrtx): flow supply.
  !             cap   (1:nedge): flow capacity.
  !             ysol  (1:nvrtx): dual solution.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nedge) :: cost, cap
  double precision, dimension(1:nvrtx) :: supply,  ysol
  
  !     ------------------------------------------------------------------
  !     double precision variables
  !
  !             tmp   : temporary variable.
  !             tobj  : temporary dual objective.
  !             ddot  : inner product function (DBLAS).
  !             pddlob: objective function value.
  !     ------------------------------------------------------------------
  
  double precision ::  tmp   , tobj  , ddot, pddlob

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Initialize dual objective value.
  !     ------------------------------------------------------------------
  
  tobj=-ddot(nvrtx,supply,1,ysol,1)

  do edg=1,nedge
     
     !        ---------------------------------------------------------------
     !        Compute reduced cost.
     !        ---------------------------------------------------------------
     
     tmp=cost(edg)-ysol(strv(edg))+ysol(endv(edg))
     if(tmp < 0.d0) tobj=tobj+tmp*cap(edg)

  end do

  pddlob = tobj
  
  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_dualobjective.
  !     ------------------------------------------------------------------

end subroutine pdnet_dualobjective
!     ------------------------------------------------------------------
!     pdnet_dualoptimal: 
!             (Orthogonally) project current interior dual solution 
!             onto hyperplane defined by spanning forest.
!     ------------------------------------------------------------------

subroutine pdnet_dualoptimal(cost,   endv,  nedge,  ntree,  nvrtx,  stedg, &
     stptr,  strv, stvrtx, yopt,   ysol)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nedge : number of edges.
  !             ntree : number of trees in spanning forest.
  !             nvrtx : number of vertices.
  !     ------------------------------------------------------------------
     
  integer :: nedge , ntree , nvrtx
  
  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             stedg (1:    nedge    ): parent edge for each vertex.
  !             stptr (1:    ntree + 1): pointer to first vertex in a tree
  !             strv  (1:  nedge): start vertex in network data structure.
  !             stvrtx(1:    nvrtx    ): triangular ordering of vertices.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: stedg, endv, strv
  integer, dimension(1:nvrtx) :: stvrtx
  integer, dimension(1:ntree + 1) :: stptr

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             cost(1:nedge): cost coefficients.
  !             yopt(1:nvrtx): boundary dual solution.
  !             ysol(1:nvrtx): interior dual solution.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nedge) :: cost
  double precision, dimension(1:nvrtx) :: yopt, ysol

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             edg   : edge index
  !             ptr   : temporary pointer
  !             tree  : tree counter
  !             vrtx  : vertex index
  !     ------------------------------------------------------------------
  
  integer :: edg   , ptr   , tree  , vrtx

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             tmp   : temporary variable.
  !     ------------------------------------------------------------------

  double precision :: tmp

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Build rhs for non-homogeneous system.
  !     ------------------------------------------------------------------
  
  do tree = 1, ntree
     yopt(stvrtx(stptr(tree)))=0.d0
     do ptr=stptr(tree)+1,stptr(tree+1)-1
        vrtx=stvrtx(ptr)
        edg=stedg(vrtx)
        yopt(vrtx)=cost(edg)
     end do
  end do

  !     ------------------------------------------------------------------
  !     Compute non-homogeneous solution.
  !     ------------------------------------------------------------------

  call pdnet_backsub(ntree,  nedge,  nvrtx,  strv,  endv, stptr,  stvrtx, &
       stedg,  yopt)

  !     ------------------------------------------------------------------
  !     Compute orthogonal projection for each tree.
  !     ------------------------------------------------------------------
  
  do tree=1,ntree
     
     !     ------------------------------------------------------------------
     !     Compute factor for homogeneous solution.
     !     ------------------------------------------------------------------
     
     tmp=0.d0

     do ptr=stptr(tree),stptr(tree+1)-1
        vrtx=stvrtx(ptr)
        tmp=tmp-ysol(vrtx)-yopt(vrtx)
     end do

     tmp=tmp/dble(stptr(tree+1)-stptr(tree))

     !        ---------------------------------------------------------------
     !        Compute projection for current tree.
     !        ---------------------------------------------------------------

     do ptr=stptr(tree),stptr(tree+1)-1
        vrtx=stvrtx(ptr)
        yopt(vrtx)=yopt(vrtx)+tmp
     end do

  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_dualoptimal.
  !     ------------------------------------------------------------------

end subroutine pdnet_dualoptimal
!     ------------------------------------------------------------------
!     pdierr: Prints error message corresponding to code ierror.
!     ------------------------------------------------------------------

subroutine pdnet_error(ierror)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output parameter:
  !
  !                ierror: error code.
  !     ------------------------------------------------------------------
  
  integer :: ierror

  !     ------------------------------------------------------------------
  !     Local data:
  !
  !       errorlist(): 	list of error messages error code.
  !       nerr:		number of of items in errlst()		
  !     ------------------------------------------------------------------
  
  integer, parameter :: nerr = 7

  character(len=60), dimension (1:nerr) :: errorlist

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------
  !
  !     ------------------------------------------------------------------
  !     Assign error messages.
  !     ------------------------------------------------------------------
  errorlist(1) = 'Node exceeds allowed label.'
  errorlist(2) = 'Maximum interior-point iterations reached.'
  errorlist(3) = 'Supply does not equal demand.'
  errorlist(4) = 'Problem instance is infeasible.'
  errorlist(7) = 'Network unable to provide the requested flow.'

  !     ------------------------------------------------------------------
  !     Print error message.
  !     ------------------------------------------------------------------  
  if (ierror > 0 .and. ierror <= nerr) then
     print *, 'PDNET runtime error : ', ierror
     print *, errorlist(ierror)
  end if

  !     ------------------------------------------------------------------ 
  !     End of subroutine pdnet_error.
  !     ------------------------------------------------------------------
end subroutine pdnet_error
  
!     ------------------------------------------------------------------
!     pdnet_father: 
!             This subroutine returns the canonical element of the set
!             containing  x. Each set is represented by a rooted tree.
!             The nodes of the tree are the elements of the set; the
!             canonical element is the root of the tree.  Each node x
!             has a pointer p(x) to its parent in the tree; the root
!             points to itself.
!
!     Reference:
!     R.E. Tarjan, "Data structures and network algorithms",SIAM,
!     Philadelphia, PA, 1983.
!     ------------------------------------------------------------------

subroutine pdnet_father(nvrtx, outpt,  parent, x)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             outpt : canonical element of the set containing x.
  !             nvrtx : number of vertices.
  !             x     : set element.
  !     ------------------------------------------------------------------
  
  integer :: nvrtx, outpt , x

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             parent(1:nvrtx): parent array for tree representation.
  !     ------------------------------------------------------------------

  integer, dimension(1:nvrtx) :: parent

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             c     : canonical element.
  !             l     : pointer to last element.
  !             t     : temporary pointer in search.
  !     ------------------------------------------------------------------
  
  integer :: c     , l     , t

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Find canonical element.
  !     ------------------------------------------------------------------

  t = x
  do while (t /= parent(t))
     t = parent(t)
  end do
  
  !     ------------------------------------------------------------------
  !     Canonical element found.
  !     ------------------------------------------------------------------
  
  c=t
  outpt=t
  
  !     ------------------------------------------------------------------
  !     Path compression.
  !     ------------------------------------------------------------------
  
  t=x

  do while (t /= parent(t))
     l=t
     t=parent(t)
     parent(l)=c
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_father.
  !     ------------------------------------------------------------------

end subroutine pdnet_father
!     ------------------------------------------------------------------
!     pdginf: Gets values of all parameters and statistics from arrays
!             info() and dinfo().
!     ------------------------------------------------------------------
subroutine pdnet_getinfo (bound , cgstop, dinfo,  dobj,   dof, factor, fcttol, &
     gap, grho0,  gtolcg, gtolmf, huge,   ierror, info,   ioflag,  iout, &
     iprcnd, maxit, mfdopt, mfstat, miu, mntlcg, mr, mrho0, mxcgit, mxfe, &
     mxfv, nbound, nbuk,  oldcos, oldtol, opt, optgap, optmf, optst , pcgcos, &
     pcgit,  pcgres, pduit,  pf, pobj, pof, prttyp, rho0,  root, s1fctr,stptol, &
     stpval, tol1, tol2, tolcg, tolcg0, tolpcg, tolslk, tolsw1, tolsw2, sw1, sw2, zero)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             info  (1:30): array with integer statistics and parameters.
  !     ------------------------------------------------------------------

  integer, dimension (1:30) :: info

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             bound : used for computing aditional data structures.
  !             cgstop: flags the way PCG stopped.
  !             ierror: flags errors in input.
  !             ioflag: controls the level of output of summary.
  !             iout  : output file unit.
  !             iprcnd: preconditioner used.
  !             maxit : maximum number of primal-dual iterations.
  !             mfstat: maxflow status.
  !             mr    : estimate of maximum number of iterations.
  !             mxcgit: maximum number of PCG iterations.
  !             mxfe  : edges on optimal flow network.
  !             mxfv  : vertices on optimal flow network.
  !             nbound: used for computing aditional data structures.
  !             nbuk  : number of buckets.
  !             opt   : flags optimal flow.
  !             optgap: Duality gap optimality indicator flag.
  !             optmf : Maximum flow optimality indicator flag.
  !             optst : Spanning tree optimality indicator flag.
  !             pcgit : number of PCG iterations performed.
  !             pduit : number of primal-dual iterations performed.
  !             pf    : iteration level.
  !             prttyp: printing level.
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------

  integer ::  bound , cgstop, ierror, ioflag, iout, iprcnd, &
       maxit , mfstat, mr    , mxcgit, mxfe  , mxfv  , &
       nbound, nbuk  , opt   , optgap, optmf , optst , pcgit, &
       pduit ,  pf    , prttyp, root 

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             dinfo  (1:35): array with double precision statistics and 
  !                            parameters.
  !     ------------------------------------------------------------------

  double precision, dimension(1:35) :: dinfo
  
  !     ----------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             dobj  : dual objective function value.
  !             dof   : dual objective function value on tree.
  !             factor: factor for Newton step.
  !             fcttol: factorization zero tolerance.
  !             gap   : value of the gap.
  !             grho0 : factor for rho0 update.
  !             gtolmf: update factor for tol1 and tol2.
  !             huge  : maxdouble precision number.
  !             mfdopt: maxflow value.
  !             miu   : value of interior-point mu parameter.
  !             mntlcg: lower bound for CG tolerance.
  !             mrho0 : lower bound for rho0.
  !             oldcos: value of cosine on previous PCG iteration.
  !             oldtol: value of residual on previous PCG iteration.
  !             pcgcos: value of cosine on PCG iteration.
  !             pcgres: value of residual on PCG iteration.
  !             pobj  : objective function value.
  !             pof   : primal objective function value on tree.
  !             rho0  : parameter from IQRD preconditioner.
  !             s1fctr: factor for miu.
  !             stptol: zero for dual slacks.
  !             stpval: largest slack on dual face.
  !             sw1   : iteration limit on CG for preconditioner switch.
  !             sw2   : another limit for switch.
  !             tol1  : lower bound of the tolerance for maxflow.
  !             tol2  : upper bound of the tolerance for maxflow.
  !             tolcg : tolerance for PCG stopping criterion.
  !             tolcg0: guess for tolerance for PCG stopping criterion.
  !             tolpcg: initial value for tolcg.
  !             zero  : value for zero.
  !     ------------------------------------------------------------------

  double precision :: dobj  , dof   , factor,  fcttol, gap   , grho0,&
       gtolcg, gtolmf, huge  ,  mfdopt,&
       miu   , mntlcg, mrho0 ,   oldcos, oldtol, &
       pcgcos, pcgres, pobj  ,  pof   , rho0  , s1fctr, stptol,&
       stpval, sw1   , sw2   , tol1  ,  tol2  , tolcg , tolcg0, tolpcg, tolslk,&
       tolsw1, tolsw2,  zero
  
  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------
  !
  !     ------------------------------------------------------------------
  !     Integer statistics assignment.
  !     ------------------------------------------------------------------
  bound  =  info( 1)
  iprcnd =  info( 2)
  iout   =  info( 3)
  maxit  =  info( 4)
  mxcgit =  info( 5)
  mxfv   =  info( 6)
  mxfe   =  info( 7)
  nbound =  info( 8)
  nbuk   =  info( 9)
  pduit  =  info(10)
  root   =  info(11)
  opt    =  info(12)
  pcgit  =  info(13)
  pf     =  info(14)
  ierror =  info(15)
  cgstop =  info(16)
  mfstat =  info(17)
  ioflag =  info(18)
  prttyp =  info(19)
  optgap =  info(20)
  optmf  =  info(21)
  optst  =  info(22)
  mr     =  info(23)
  
  !     ------------------------------------------------------------------
  !     Double precision statistics assignment.
  !     ------------------------------------------------------------------
  dobj   = dinfo( 1)
  dof    = dinfo( 2)
  factor = dinfo( 3)
  fcttol = dinfo( 4)
  gap    = dinfo( 5)
  grho0  = dinfo( 6)
  gtolcg = dinfo( 7)
  gtolmf = dinfo( 8)
  mfdopt = dinfo( 9)
  miu    = dinfo(10)
  mntlcg = dinfo(11)
  mrho0  = dinfo(12)
  oldtol = dinfo(13)
  pcgres = dinfo(14)
  pobj   = dinfo(15)
  pof    = dinfo(16)
  rho0   = dinfo(17)
  s1fctr = dinfo(18)
  stptol = dinfo(19)
  stpval = dinfo(20)
  tol1   = dinfo(21)
  tol2   = dinfo(22)
  tolcg  = dinfo(23)
  tolcg0 = dinfo(24)
  tolslk = dinfo(25)
  huge   = dinfo(26)
  tolpcg = dinfo(27)
  oldcos = dinfo(28)
  pcgcos = dinfo(29)
  zero   = dinfo(30)
  tolsw1 = dinfo(31)
  tolsw2 = dinfo(32)
  sw1    = dinfo(33)
  sw2    = dinfo(34)
  
  
  !    ------------------------------------------------------------------
  !     End of subroutine pdnet_getinfo.
  !     ------------------------------------------------------------------
end subroutine pdnet_getinfo
!     -------------------------------------------------------------------------
!     pdnet_header: Print report header with program identification, problem
!                   instance information, run-time and compile-time parameters.
!     -------------------------------------------------------------------------

subroutine pdnet_header (dinfo, info, na, nn)

  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nn    : number of nodes.  
  !             na    : number of arcs.
  !     ------------------------------------------------------------------

  integer :: nn, na, iout
  
  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             info  (1:30): array with many integer statistics.
  !     ------------------------------------------------------------------

  integer, dimension(1:30) :: info

  !     ------------------------------------------------------------------
  !     Double precision input arrays:
  !
  !             dinfo  (1:35): array with many double precision statistics.
  !     ------------------------------------------------------------------
  double precision, dimension (1:35) :: dinfo

  !     ------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------
  !
  !     --------------------
  !     Write problem specs.
  !     --------------------

  iout = info(3)
  write (iout,'(A)')  '-----------------------------------------------'
  write (iout, '(A)') 'PDNET 2.0     Release date: 2005-06-01'
  write (iout, '(A)') 'Primal-feasible dual-infeasible interior point '
  write (iout, '(A)') 'network flow method'
  write (iout, '(A)') '-----------------------------------------------'
  write (iout, '(A)') 'Problem Size : '
  write (iout, '(2x,a,i14)') 'Number of Nodes : ', nn
  write (iout, '(2x,a,i14)') 'Number of Arcs  : ', na

  !     ------------------
  !     Runtime parameters
  !     ------------------
  write (iout, '(A)') '-----------------------------------------------'
  write (iout, '(a)') 'Runtime parameters : '
  write (iout, '(2(2x,a,i14))') '      CG Preconditioner:', info( 2), &
                                '           Output level:', info(19)
  write (iout, '(2(2x,a,i14))') '     Maximum iterations:', info( 4), &
                                '  Maximum CG iterations:', info( 5)
  write (iout, '(2(2x,a,i14))') ' Duality Gap Optimality:', info(20), &
                                '   Max. Flow Optimality:', info(21)
  write (iout, '(  2x,a,i14 )') '  Span. Tree Optimality:', info(22)

  !     -----------------------
  !     Compile-time parameters
  !     -----------------------
  write (iout, '(a)') '-----------------------------------------------'
  write (iout, '(a)') 'Compile-time parameters : '
  write (iout, '(  2x,a,i14 )')  '      Output file unit:', iout
  write (iout, '(2(2x,a,i14))')  '    Spanning tree root:', info(11), &
                                 'Number of sort buckets:', info( 9)
  write (iout, '(2x,a,e14.8,2x,a,i14)') &
                                 'Factor for Newton step:',dinfo( 3), &
                                 'Number of sort buckets:', info( 9)
  !     ----------------------------------------------------------------
  !    End of subroutine pdnet_header.
  !   ------------------------------------------------------------------
end subroutine pdnet_header
!     -------------------------------------------------------------------
!     pdnet_heap: This subroutine picks-up a minimum spanning tree from a
!                 directed network. It uses a Fibonacci heap data structure
!                 to maintain the set of nodes that have not yet been 
!                 included in the spanning tree.
!     ------------------------------------------------------------------

subroutine pdnet_heap (bfirst, bposit, bucket, endn,  father, first, huge,&
     infrst, innext, key, larcs, lfirst, list, llink, lost, lsets, mark, mr,&
     na,  nn, oufrst, ounext, pred, rank, rlink, root, rsets, rtreth, strn, &
     teta, trearc, treath, uarcs, ufirst, weight)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              bfirst: used for building maximum spanning tree.
  !              larcs : auxiliary variable.
  !              lfirst: auxiliary variable.
  !              mr    : estimate of maximum number of iterations.
  !              na    : number of arcs.
  !              nn    : number of nodes.
  !              root  : sets beginning of tree structure.
  !     ------------------------------------------------------------------
  
  integer :: bfirst, larcs, mr, na, nn, root

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             bposit(1:na): flags arcs according to bounds.
  !             bucket(0:mr): stores buckets.
  !             endn  (1:na): pointer of network data structure.
  !             father(1:nn): marks the father of each node.
  !             first (1:nn): first node.
  !             infrst(1:nn): points to the first node in in-adcncy list.
  !             innext(1:na): points to the last node in in-adcncy list.
  !             list  (1:nn): temporary list of nodes.
  !             llink (0:nn): part of the node left linked list.
  !             lost  (1:nn): number of lost nodes per node.
  !             lsets (1:na): list of arcs at lower bound.
  !             mark  (1:nn): selects nodes to belong to the tree.
  !             pred  (1:nn): list of predecessors.
  !             rank  (1:nn): rank of nodes.
  !             rlink (0:nn): part of the node right linked list.
  !             rsets (1:na): list of arcs at upper bound.
  !             rtreth(1:nn): adjacent to a given node.
  !             strn  (1:na): adjacency list of network data structure.
  !             trearc(1:nn): arc of tree corresponding to node.
  !             treath(1:nn): node that is adjacent to given.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: father, first, infrst, list, lost, mark,&
       oufrst, pred, rank, rtreth, trearc, treath
  integer, dimension(1:na) :: bposit, endn, innext, lsets, ounext, rsets,&
       strn
  integer, dimension(0:nn) :: llink, rlink
  integer, dimension(0:mr) :: bucket

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !              key (0:nn): array of key of nodes.
  !              teta(1:na): array of teta values of arcs.
  !     ------------------------------------------------------------------

  double precision, dimension(0:nn) :: key
  double precision, dimension(1:na) :: teta

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !              arc   : current arc.
  !              larc  : auxiliary variable.
  !              lauxrc: variable used for setting up rsets.
  !              lfirst: first arc at the lower bound.
  !              minkey: minimum key.
  !              node  : loop counter.
  !              oldmk : variable used for setting up mark.
  !              rarc  : auxiliary variable.
  !              snod  : stores a node.
  !              uarcs : arcs at the upper bound.
  !              ufirst: first arc at the upper bound.
  !     ------------------------------------------------------------------
  
  integer :: arc   , enod, larc  , lauxrc, lfirst, minkey, node, oldmk , &
       rarc  , snod  , uarcs , ufirst

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !              huge  : very huge number.
  !              value : auxiliary variable.
  !              weight: weight of the tree.
  !     ------------------------------------------------------------------ 

  double precision ::  huge  , value , weight

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !
  !     Initializations.
  !     ------------------------------------------------------------------

  do node = 1, nn
     key(node)     = huge
     pred(node)    = 0
     first(node)   = 0
     llink(node)   = 0
     rlink(node)   = 0
     rank(node)    = 0
     lost(node)    = 0
     list(node)    = 0
     father(node)  = 0
     mark(node)    = 0
     treath(node)  = 0
     rtreth(node)  = 0
     trearc(node)  = 0
  end do

  do node = 0, mr
     bucket(node) = 0
  end do

  key(0)        = huge
  llink(0)      = 0
  rlink(0)      = 0
  weight        = 0d0
  bfirst  = 0  

  !     ------------------------------------------------------------------
  !     Define minkey.
  !     ------------------------------------------------------------------
  
  minkey = root
  
  !     ------------------------------------------------------------------
  !     Insert nodes in the heap.
  !     ------------------------------------------------------------------  

  do node = 1, nn
     call pdnet_insert (nn,     mr,     node,   pred,   first,  llink, &
          rlink,  rank,   bucket, key,    minkey )
  end do

  !     ------------------------------------------------------------------
  !     Prim's algorithm.
  !     ------------------------------------------------------------------  

  node = 1
  do while (node /= nn)
     node = node + 1
     
     !     ------------------------------------------------------------------
     !     Select node minkey to belong to the tree.
     !     ------------------------------------------------------------------
     
     oldmk = minkey
     mark(oldmk) = node-1
     
     !     ------------------------------------------------------------------
     !     Delete node minkey from the heap.
     !     ------------------------------------------------------------------
     
     call pdnet_min (bucket, first,  huge,   key,    llink,  lost, minkey, &
          mr,     nn,     pred,   rank,   rlink)
     
     !     ------------------------------------------------------------------
     !     Decreasing the key of the nodes adjacent to the node oldmk.
     !     
     !     Incoming arcs:
     !     ------------------------------------------------------------------
     
     arc = infrst(oldmk)
     
     do while (arc /= 0)
        snod = strn(arc)
        if (key(snod) > teta(arc) .and. mark(snod) == 0) then
           trearc(snod) = arc
           value = teta(arc)
           father(snod) = oldmk
           call pdnet_key (bucket, first, snod, key, list, &
                llink,lost, minkey, mr, nn, pred, rank, rlink, value)
        end if
        arc = innext(arc)
     end do

     !     ------------------------------------------------------------------
     !     Outgoing arcs:
     !     ------------------------------------------------------------------

     arc = oufrst(oldmk)
     do while (arc /= 0)
        enod = endn(arc)
        if (key(enod) > teta(arc) .and. mark(enod) == 0) then
           trearc(enod) = arc
           value = teta(arc)
           father(enod) = oldmk
           call pdnet_key (bucket, first, enod, key, list, llink, lost, minkey,&
                mr, nn, pred, rank, rlink, value)
        end if
        arc = ounext(arc)
     end do


     
     !     ------------------------------------------------------------------
     !     Update vectors treath, rtreth, lsets and rsets.
     !     ------------------------------------------------------------------
     
     treath(oldmk) = minkey
     rtreth(minkey) = oldmk
     arc  = trearc(minkey)
     if (arc /= 0) then
        larc = lsets(arc)
        rarc = rsets(arc)
        if (rarc /= 0) lsets(rarc) = larc
        if (larc /= 0) rsets(larc) = rarc
        if (arc == lfirst) lfirst=rarc
        if (arc == ufirst) ufirst=rarc
        lsets(arc) = 0
        rsets(arc) = 0
        if (bfirst == 0) then
           bfirst         = arc
           lauxrc        = bfirst
        else
           rsets(lauxrc) = arc
           lsets(arc)     = lauxrc
           lauxrc        = arc     
        end if
        if (bposit(arc) == -1) larcs=larcs-1
        if (bposit(arc) ==  1) uarcs=uarcs-1
        bposit(arc) = 0
     else
        stop 'disconnected network'
     end if

     !     ------------------------------------------------------------------
     !     Update tree weight.
     !     ------------------------------------------------------------------
     
     weight = weight + key(minkey)

  end do

  !     ------------------------------------------------------------------
  !     Update vectors treath and rtreth.
  !     ------------------------------------------------------------------

  treath(minkey) = root
  rtreth(root)  = minkey
  mark(minkey)   = node

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_heap.
  !     ------------------------------------------------------------------

end subroutine pdnet_heap
!     ----------------------------------------------------------------------
!     pdnet_insert: This subroutine inserts a new tree in the Fibonacci heap 
!                   and performs multilinking to restore the heap invariants.
!     ----------------------------------------------------------------------

subroutine pdnet_insert (nn,     mr,     node,   pred,   first,  llink, &
     rlink,  rank,   bucket, key,    minkey)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             minkey: minimum key.
  !             mr    : estimate of maximum number of iterations.
  !             nn    : number of nodes.
  !             node  : node to insert.
  !     ------------------------------------------------------------------

  integer :: minkey, mr    , nn    , node
  
  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             bucket(0:mr): stores buckets.
  !             first (1:nn): first node.
  !             llink (0:nn): part of the node left linked list.
  !             pred  (1:nn): list of predecessors.
  !             rank  (1:nn): rank of nodes.
  !             rlink (0:nn):part of the node right linked list.
  !     ------------------------------------------------------------------  

  integer, dimension(1:nn) :: first, pred, rank
  integer, dimension(0:nn) :: llink, rlink
  integer, dimension(0:mr) :: bucket

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             key (0:nn): array of key of nodes.
  !     ------------------------------------------------------------------

  double precision, dimension(0:nn) :: key
  
  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             cnode : copy of node.
  !             croot : copy of root.
  !             root  : bucket of a given node.
  !     ------------------------------------------------------------------
 
  integer :: cnode , croot , root

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------

  cnode = node

  do
     root = bucket(rank(cnode))
     if (root == 0) then
        bucket(rank(cnode)) = cnode
        return
     else
        bucket(rank(cnode)) = 0
        if (key(root) > key(cnode) .or. cnode == minkey) then
           croot = root
           root  = cnode
           cnode = croot
        end if
        call pdnet_link(nn, root, cnode, pred, first, llink, &
             rlink,  rank)
        cnode = root
     end if
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_insert.
  !     ------------------------------------------------------------------

end subroutine pdnet_insert
!     ------------------------------------------------------------------
!     pdnet_iqrd: 
!             This subroutine constructs the diagonal matrix diag and
!             the vector ndiag wich stores the nondiagonal elements of
!             the matrix f.
!     ------------------------------------------------------------------

subroutine pdnet_iqrd(av4, diag, endn, father, fcttol, iprcnd, na, ndiag,&
     nn, rho0, root, rtreth, strn, teta, trearc)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             iprcnd: preconditioner used.
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------
  
  integer :: iprcnd, na    , nn    , root

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endn  (1:na): pointer of network data structure.
  !             father(1:nn): marks the father of each node.
  !             rtreth(1:nn): adjacent to a given node.
  !             strn  (1:na): adjacency list of network data structure.
  !             trearc(1:nn): arc of tree corresponding to node.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: father, rtreth, trearc
  integer, dimension(1:na) :: endn, strn 

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             av4   (1:na): auxiliary vector.
  !             diag  (1:nn): diagonal of matrix.
  !             ndiag (1:nn): diagonal of preconditioner.
  !             teta  (1:na): diagonal theta.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: diag,  ndiag
  double precision, dimension(1:na) :: av4, teta
  
  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             arc   : current arc.
  !             enod  : current element on endn.
  !             epred : father of enod.
  !             i     : loop counter.
  !             node  : another loop counter.
  !             pred  : father of node.
  !             snod  : current element on strn.
  !             spred : father of snode.
  !     ------------------------------------------------------------------
  
  integer :: arc, enod, epred, i, node, pred, snod, spred
  
  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             fcttol: factorization zero tolerance.
  !             rho0  : parameter from IQRD preconditioner.
  !     ------------------------------------------------------------------
  
  double precision :: fcttol, rho0

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             aux1  : auxiliary variable 1.
  !             aux2  : auxiliary variable 2.
  !             aux3  : auxiliary variable 3.
  !             aux4  : auxiliary variable 4.
  !             maxav4: maximum element of vector av4.
  !             minav4: minimum element of vector av4.
  !             rho   : parameter from IQRD preconditioner.
  !     ------------------------------------------------------------------

  double precision :: aux1  , aux2  , aux3  , aux4  , maxav4, minav4, rho

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  do node=1,nn
     diag(node)  = 0d0
     ndiag(node) = 0d0
  end do

  minav4=1.0d64
  maxav4=0.0d0

  do arc=1,na
     av4(arc)    = 1d0/teta(arc)
     snod  = strn(arc)
     enod  = endn(arc)
     if (arc == trearc(snod) .or. arc == trearc(enod)) then
        maxav4=dmax1(maxav4,av4(arc))
        minav4=dmin1(minav4,av4(arc))
     end if
  end do

  if (maxav4 /= 0.0d0) then
     rho=rho0*dsqrt(minav4/maxav4)
  else
     rho=rho0
  end if

  !     ------------------------------------------------------------------
  !     Compute diag.
  !     ------------------------------------------------------------------

  if (iprcnd == 1) then
     do arc=1,na
        snod  = strn(arc)
        enod  = endn(arc)
        diag(enod)  = diag(enod)+av4(arc)
        diag(snod)  = diag(snod)+av4(arc)
     end do
     diag(root) = 0.0d0
  end if

  if (iprcnd == 0) then
     do arc=1,na
            snod  = strn(arc)
            enod  = endn(arc)
            spred = father(snod)
            epred = father(enod)
            if (snod == epred) then
               diag(enod)  =  diag(enod)+av4(arc)
               if (snod == root) goto 40
               ndiag(enod) = ndiag(enod)+av4(arc)
            else
               if (enod == spred) then
                  diag(snod)  =  diag(snod)+av4(arc)
                  if (enod == root) goto 40
                  ndiag(snod) = ndiag(snod)+av4(arc)
               else
                  diag(enod)  = diag(enod)+av4(arc)
                  diag(snod)  = diag(snod)+av4(arc)
               end if
            end if
40       end do

         diag(root) = 0d0

         do  node = 1,root-1
            ndiag(node) = ndiag(node)/diag(node)
         end do

         do node = root+1,nn
            ndiag(node) = ndiag(node)/diag(node)
         end do

         ndiag(root) = 0d0
      end if
      
      if (iprcnd == 4) then
         do arc = 1,na
            snod = strn(arc)
            enod = endn(arc)
            spred = father(snod)
            epred = father(enod)
            if (arc == trearc(snod) .or. arc == trearc(enod)) then
               diag(enod) = diag(enod)+av4(arc)
               diag(snod) = diag(snod)+av4(arc)
            else
               diag(enod) = diag(enod)+rho*av4(arc)
               diag(snod) = diag(snod)+rho*av4(arc)
            endif
            if (snod == epred) then
               if (snod == root) goto 70
               if (arc /= trearc(enod)) goto 70
               ndiag(enod) = ndiag(enod)+av4(arc)
            else
               if (enod == spred) then
                  if (enod == root) goto 70
                  if (arc /= trearc(snod)) goto 70
                  ndiag(snod) = ndiag(snod)+av4(arc)
               end if
            end if
70       end do
         
         diag(root) = 0d0
         ndiag(root) = 0d0
         node = rtreth(root)
         
         do i = 2, nn
            pred = father(node)
            aux1=diag(node)
            if (aux1 < fcttol) then
               aux1 = 1.d50
            endif
            aux2 = ndiag(node)
            aux3 = diag(pred)
            aux4 = aux2/aux1
            ndiag(node) = aux4
            diag(pred) = aux3 - aux2*aux4
            node = rtreth(node)
         end do
      end if

!     ------------------------------------------------------------------
!     End of subroutine pdnet_irqd.
!     ------------------------------------------------------------------

    end subroutine pdnet_iqrd
!     ---------------------------------------------------------------------
!     pdnet_iqrdsolv: This subroutine solves a system with the IQRD matrix.
!     ---------------------------------------------------------------------

subroutine pdnet_iqrdsolv(b, diag, father, iprcnd, na, ndiag,  nn, root, &
     rtreth, teta, trearc, treath, x)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             iprcnd: preconditioner used.
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------
  
  integer :: iprcnd, na    , nn    , root

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             father(1:nn): marks the father of each node.
  !             rtreth(1:nn): adjacent to a given node.
  !             trearc(1:nn): arc of tree corresponding to node.
  !             treath(1:nn): node that is adjacent to given.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: father, rtreth, trearc,treath

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             b    (1:nn): array of capacities at each node.
  !             diag (1:nn): diagonal of matrix.
  !             ndiag(1:nn): diagonal of preconditioner.
  !             teta (1:na): diagonal theta.
  !             x    (1:nn): primal iteration x.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: b, diag, ndiag, x
  double precision, dimension(1:na) :: teta

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             node  : another loop counter.
  !             pred  : father of node.
  !     ------------------------------------------------------------------
  
  integer :: i     , node  , pred

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !
  !     Initializing the vector x.
  !     ------------------------------------------------------------------

  do node=1,root-1
     x(node) = b(node)
  end do

  do node=root+1,nn
     x(node) = b(node)
  end do

  x(root) = 0.0d0

  !     ------------------------------------------------------------------
  !     Reverse treath traversal of the spanning tree.
  !     ------------------------------------------------------------------

  if (iprcnd == 3) then
     node = rtreth(root)
     do i = 2, nn
        pred = father(node)
        x(pred) = x(pred)+x(node)
        node = rtreth(node)
     end do
  end if

  if (iprcnd == 0 .or. iprcnd == 4) then
     node = rtreth(root)
     do i=2,nn
        pred = father(node)
        x(pred) = x(pred)+x(node)*ndiag(node)
        node = rtreth(node)
     end do
  end if

  x(root) = 0.0d0
  
  !     ------------------------------------------------------------------
  !     Dividing x by diag.
  !     ------------------------------------------------------------------

  if (iprcnd == 0 .or. iprcnd == 1 .or. iprcnd == 4) then
     do node=1,root-1
        x(node) = x(node)/diag(node)
     end do
     do node=root+1,nn
        x(node) = x(node)/diag(node)
     end do
  end if
  if (iprcnd == 3) then
     do node=1,root-1
        x(node) = x(node)*teta(trearc(node))
     end do
     do node=root+1,nn
        x(node) = x(node)*teta(trearc(node))
     end do
  end if

  !     ------------------------------------------------------------------
  !     Treath traversal of the spanning tree.
  !     ------------------------------------------------------------------
      
  if (iprcnd == 3) then
     node = treath(root)
     do i=2,nn
        pred = father(node)
        x(node) = x(node)+x(pred)
        node = treath(node)
     end do
  end if

  if (iprcnd == 0 .or. iprcnd == 4) then
     node = treath(root)
     do i=2,nn
        pred = father(node)
        x(node) = x(node)+x(pred)*ndiag(node)
        node = treath(node)
     end do
  end if

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_iqrdsolv.
  !     ------------------------------------------------------------------

end subroutine pdnet_iqrdsolv
!     ------------------------------------------------------------------
!     pdnet_key: This subroutine decreases the key of a node in the
!             Fibonacci heap and performs multispliting and multilinking
!             to restore the invariants of the heap.
!     ------------------------------------------------------------------

subroutine pdnet_key (bucket, first,  k, key, list, llink, lost, minkey, &
     mr, nn, pred, rank, rlink, value)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              k     : node whose key is to be decreased.
  !              minkey: minimum key.
  !              mr    : estimate of maximum number of iterations.
  !              nn    : number of nodes.
  !     ------------------------------------------------------------------
  
  integer :: k     , minkey, mr    , nn

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !              bucket(0:mr): stores buckets.
  !              first (1:nn): first node.
  !              list  (1:nn): temporary list of nodes.
  !              llink (0:nn): part of the node left linked list.
  !              lost  (1:nn): number of lost nodes per node.
  !              pred  (1:nn): list of predecessors.
  !              rank  (1:nn): rank of nodes.
  !              rlink (0:nn): part of the node right linked list.
  !     ------------------------------------------------------------------  

  integer, dimension(1:nn) :: first, list, lost, pred, rank
  integer, dimension(0:nn) :: llink, rlink
  integer, dimension(0:mr) :: bucket

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !              key   (0:nn): array of key of nodes.
  !     ------------------------------------------------------------------

  double precision, dimension(0:nn) :: key

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !              father: father of node k.
  !              i     : loop counter.
  !              itotal: total of nodes processed.
  !              kcopy : copy of k.
  !     ------------------------------------------------------------------
  
  integer :: father, i     , itotal, kcopy

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !              value : value the key is decreased by.
  !     ------------------------------------------------------------------  

  double precision :: value

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------

  kcopy = k
  key(kcopy) = value
  father = pred(kcopy)  
  
  if (key(minkey) > value) then
     minkey = kcopy
  end if

  if (father == 0 .or. key(father) <= value) then
     return
  end if

  !     ------------------------------------------------------------------
  !     Multisplitting.
  !     ------------------------------------------------------------------
 
  itotal = 0
  do
     itotal = itotal + 1
     list(itotal) = kcopy
     call pdnet_split (father, first, llink, lost, nn, pred, rank, rlink, kcopy)
     if (lost(father) < 2) then
        exit
     end if
     kcopy = father
     father = pred(kcopy)
  end do

  if (pred(father) == 0) then
     itotal = itotal+1
     list(itotal) = father
     bucket(rank(father)+1) = 0
  end if

  !     ------------------------------------------------------------------
  !     Multilinking.
  !     ------------------------------------------------------------------ 

  do i = 1, itotal
     kcopy = list(i)
     call pdnet_insert(nn, mr, kcopy, pred, first, llink, rlink, rank, &
          bucket, key, minkey)
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pddkey.
  !     ------------------------------------------------------------------

end subroutine pdnet_key
!     ------------------------------------------------------------------
!     pdnet_kruskal: Given a graph g and a ordering of the edges, this
!             subroutine applies Kruskal's algorithm and builds a
!             maximum spanning forest.
!
!             An ordered list with "nord" indices of allowable edges is
!             given in the input array "ordedg".  Edges are identified
!             by their start and end vertices in arrays strv and  endv.
!
!             Output is expressed by 3 items:
!
!                   1. ntree:  number of trees in spanning forest.
!                   2. stedg:  list of (nvrtx - ntree) edges in spanning
!                              forest.
!                   3. parent: tree index for each vertex.
!
!             Arrays "stedg" and "ordedg" can point to the same memory 
!             location.
!     ------------------------------------------------------------------

subroutine pdnet_kruskal(endv,  nord,  ntree,  nvrtx,  nedge, ordedg, &
     parent, rank, stedg,   strv)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nord  : number of edges in ordered list.
  !             ntree : number of trees in spanning forest.
  !             nvrtx : number of vertices.
  !             nedge : number of edges.
  !     ------------------------------------------------------------------
  
  integer :: nord  , ntree , nvrtx, nedge

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             ordedg(1:  nord ): edge ordering.
  !             parent(1:  nvrtx): parent array for tree representation.
  !             rank  (1:  nvrtx): upper bound on height of element.
  !             stedg (1:  nord) : edges in spanning forest.
  !             strv  (1:  nedge): start vertex in network data structure.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: endv, strv
  integer, dimension(1:nord ) :: ordedg, stedg 
  integer, dimension(1:nvrtx) :: parent, rank

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             e     : edge index.
  !             ntedg : count of edges in spanning forest.
  !             ord   : ordered edge list index.
  !             partv : canonical element for set containing v.
  !             partw : canonical element for set containing w.
  !             v     : vertex index.
  !             vrtx  : vertex index.
  !             w     : vertex index.
  !     ------------------------------------------------------------------
  
  integer :: e     , ntedg , ord   , partv , partw , v     , vrtx  , w

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Make sets -- initialize parent and rank.
  !     ------------------------------------------------------------------

  do vrtx = 1, nvrtx
     parent(vrtx)=vrtx
     rank(vrtx)=0
  end do

  !     ------------------------------------------------------------------
  !     Initialize number of edges in spanning tree (forest).
  !     ------------------------------------------------------------------
  
  ntedg=0
  do ord = 1, nord

     !        ---------------------------------------------------------------
     !        Identify next edge.
     !        ---------------------------------------------------------------
     
     e=ordedg(ord)
     v=strv(e)
     w=endv(e)
     
     !        ---------------------------------------------------------------
     !        If v and w are from different partitions join both partitions.
     !        ---------------------------------------------------------------

     call pdnet_father( nvrtx, partv,  parent, v )
     call pdnet_father( nvrtx, partw,  parent, w )
     
     if (partv /= partw) then
        call pdnet_linkstruct(nvrtx, partv,  partw,  parent, rank)

        !           ------------------------------------------------------------
        !           Add e=(v,w) to max weighted spanning tree t.
        !           ------------------------------------------------------------
  
        ntedg=ntedg+1
        stedg(ntedg)=e
     end if

     !           ------------------------------------------------------------
     !           Break if tree is complete.
     !           ------------------------------------------------------------
     
     if ((partv /= partw) .and. (ntedg == nvrtx - 1)) exit
  end do

  !     ------------------------------------------------------------------
  !     Initialize number of trees in spanning forest.
  !     ------------------------------------------------------------------

  ntree = 0

  !     ------------------------------------------------------------------
  !     Associate vertices with roots, count vertices for each root.
  !     ------------------------------------------------------------------

  do vrtx = 1, nvrtx
     call pdnet_father(nvrtx,  v,      parent, vrtx )
     if (v == vrtx) then
        ntree=ntree+1
        rank(vrtx)=ntree
     end if
  end do

  !     ------------------------------------------------------------------
  !     Translate entries in parent array and count.
  !     ------------------------------------------------------------------
  
  do vrtx = 1, nvrtx
     parent(vrtx)=rank(parent(vrtx))
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_kruskal.
  !     ------------------------------------------------------------------

end subroutine pdnet_kruskal
!     ------------------------------------------------------------------
!     pdnet_l1norm: 
!             This subroutine computes the l1 norm of z. The vector z
!             has dimension dim.
!     ------------------------------------------------------------------

subroutine pdnet_l1norm(alfa,   dim,    z)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             dim: dimension of vector.
  !     ------------------------------------------------------------------
  
  integer :: dim

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             z(1:dim): vector whose norm is to be computed.
  !     ------------------------------------------------------------------

  double precision, dimension(1:dim) :: z

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             m     : size of unrolled loop.
  !             mp1   : m+1.
  !     ------------------------------------------------------------------
  
  integer :: i     , m     , mp1

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             alfa  : l-1 norm of vector.
  !     ------------------------------------------------------------------

  double precision :: alfa

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------
  
  alfa = 0d0
  m = mod(dim,4)
  
  if (m /= 0) then
     do i = 1,m
        alfa = alfa+dabs(z(i))
     end do
     if (dim < 4) return
  end if

  mp1 = m + 1
  do i = mp1,dim,4
     alfa=alfa+dabs(z(i))+dabs(z(i+1))+dabs(z(i+2))+dabs(z(i+3))
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_l1norm.
  !     ------------------------------------------------------------------

end subroutine pdnet_l1norm
!     --------------------------------------------------------------------
!     pdnet_link: This subroutine links two subtrees in a Fiboncci heap by
!                 updating the arrays pred, rank, llink, rlink  and first.
!     --------------------------------------------------------------------

subroutine pdnet_link (nn, root, nsoon, pred, first, llink, rlink, rank)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nn    : number of nodes.
  !             nsoon : son of node.
  !             root  : root node.
  !     ------------------------------------------------------------------
  
  integer :: nn    , nsoon , root

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             first (1:nn): first node.
  !             llink (0:nn): part of the node left linked list.
  !             pred  (1:nn): list of predecessors.
  !             rank  (1:nn): rank of nodes.
  !             rlink (0:nn): part of the node right linked list.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: first, pred, rank
  integer, dimension(0:nn) :: llink, rlink

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             osoon : auxiliary variable.
  !     ------------------------------------------------------------------
  
  integer osoon

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  osoon = first(root)
  rlink(nsoon) = osoon
  llink(osoon) = nsoon
  pred(nsoon) = root
  rank(root) = rank(root) + 1
  first(root) = nsoon

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_link.
  !     ------------------------------------------------------------------  

end subroutine pdnet_link
  
!     ------------------------------------------------------------------
!     pdlkst: This subroutine forms a new set that is the union of the 
!             two sets whose canonical elements are  x and y, destroying 
!             the two old sets.  select and return the canonical element 
!             of the new set. This operation assumes x != y.  This 
!             routine uses a union by rank heiristic.
!
!     ARGUMENTS:
!
!     n         dimension of vectors parent and rank.
!     x         canonical element in first set
!     y         canonical element in second set
!     parent()  parent array for tree representation
!     rank()    upper bound on height of element.
!
!     local variables:
!
!     c   canonical element
!
!     Reference:
!     R.E. Tarjan, "Data structures and network algorithms",SIAM,
!     Philadelphia, PA, 1983
!     ------------------------------------------------------------

subroutine pdnet_linkstruct(n, x,      y,      parent, rank)

  implicit none

  !     ------------------------------------------------------------------
  !     Storage definition.
  !     ------------------------------------------------------------------

  integer :: n, x, y
  integer, dimension(1:n) :: parent, rank

  !-------------------------------------------------------------------

  if (rank(x) < rank(y)) then
     parent(x) = y
     return
  else if (rank(x) > rank(y)) then
     parent(y) = x
     return
  else
     parent(y) = x
     rank(x) = rank(x) + 1
     return
  end if

  !-------------------------------------------------------------------

end subroutine pdnet_linkstruct
!     -------------------------------------------------------------------
!     pdnet_min: This subroutine deletes the node with minimum key in the
!                 Fibonnacci heap.
!     -------------------------------------------------------------------

subroutine pdnet_min (bucket, first,  huge,   key,    llink,  lost, &
     minkey, mr,     nn,     pred,   rank,   rlink)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              minkey: minimum key.
  !              mr    : estimate of maximum number of iterations.
  !              nn    : number of nodes.
  !     ------------------------------------------------------------------

  integer ::  minkey, mr    , nn

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             bucket(0:mr): stores buckets.
  !             first (1:nn): first node.
  !             llink (0:nn): part of the node left linked list.
  !             lost  (1:nn): number of lost nodes per node.
  !             pred  (1:nn): list of predecessors.
  !             rank  (1:nn): rank of nodes.
  !             rlink (0:nn): part of the node right linked list.
  !     -----------------------------------------------------------------

  integer, dimension(1:nn) :: first, lost, pred, rank
  integer, dimension(0:nn) :: llink, rlink
  integer, dimension(0:mr) :: bucket

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !              key   (0:nn) : array of key of nodes.
  !     ------------------------------------------------------------------

  double precision, dimension(0:nn) :: key

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !              node  : element from array first.
  !              nsoons: loop counter.
  !              rnode : element from array rlink.
  !     ------------------------------------------------------------------
  
  integer :: node  , nsoons, rnode

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !              auxkey : auxiliary variable.
  !              huge   : very huge number.
  !     ------------------------------------------------------------------  

  double precision :: auxkey, huge

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------
      
  bucket(rank(minkey)) = 0

  !     ------------------------------------------------------------------
  !     Restore invariants of the fibonacci heap.
  !     ------------------------------------------------------------------
      
  node = first(minkey)

  do while (node /= 0)
     rnode = rlink(node)
     call pdnet_split (minkey, first, llink, lost, nn, pred, &
          rank, rlink, node)
     call pdnet_insert (nn, mr, node, pred, first, llink, &
          rlink, rank, bucket, key, minkey)
     node = rnode
  end do

  !     ------------------------------------------------------------------
  !     Recompute minkey.
  !     ------------------------------------------------------------------

  auxkey = huge
  do nsoons = 0, mr
     node = bucket(nsoons)
     if (node /= 0 .and. key(node) <= auxkey) then
        minkey = node
        auxkey = key(minkey)
     end if
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pddmin.
  !     ------------------------------------------------------------------

end subroutine pdnet_min
!     ------------------------------------------------------------------
!     pdoptm: This subroutine checks the optimality of the sets b, l,
!             and u indicated by the primal-dual interior point method.
!             If the sets are optimal, then the optimal flows are
!             computed.
!     ------------------------------------------------------------------

subroutine pdnet_optcheck (alfa, av4, b, bound, bstats, c, dof, endn, father, &
     mark, na, nbound, nn, nnodes, npfrst, opt, optflo, pf, pof, root, rsets, &
     rtreth, strn, trearc, treath, u, uarcs, ufirst, upfrst, upnext, x, y, y0, zero)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             bound : used for computing aditional data structures.
  !             na    : number of arcs.
  !             nbound: used for computing aditional data structures.
  !             nn    : number of nodes.
  !             npfrst: auxiliary variable.
  !             opt   : flags optimal flow.
  !             pf    : iteration level.
  !             root  : sets beginning of structure.
  !             uarcs : number of arcs at the upper bound.
  !             ufirst: auxiliary variable.
  !             upfrst: arc with upper bound.
  !     ------------------------------------------------------------------

  integer :: bound, na, nbound, nn, npfrst, opt, pf, root, uarcs, ufirst, upfrst 

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             b     (1:nn): array of capacities at each node.
  !             bstats(1:na): flags if each arc has an upper bound or not.
  !             c     (1:na): array of costs at each arc.
  !             endn  (1:na): pointer of network data structure.
  !             father(1:nn): marks the father of each node.
  !             mark  (1:nn): selects nodes to belong to the tree.
  !             nnodes(1:nn): number of nodes adjacent to node.
  !             optflo(1:na): array with optimal flow at each arc.
  !             rsets (1:na): list of arcs at upper bound.
  !             rtreth(1:nn): adjacent to a given node.
  !             strn  (1:na): adjacency list of network data structure.
  !             trearc(1:nn): arc of tree corresponding to node.
  !             treath(1:nn): node that is adjacent to given.
  !             u     (1:na): upper bound of flow at each arc.
  !             upnext(1:na): next arc with positive upper bound.
  !             x     (1:nn): current solution.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: b, father, mark, nnodes, trearc, treath, x
  integer, dimension(1:na) :: bstats, c, endn, optflo, rsets, rtreth, strn, &
       u, upnext

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             alfa(1:nn): auxiliary vector.
  !             av4 (1:na): another auxiliary vector.
  !             y   (1:nn): variable y.
  !             y0  (1:nn): copy of y.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: alfa, y, y0
  double precision, dimension(1:na) :: av4

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             arc   : arc being processed.
  !             enod  : current element on endn.
  !             i     : loop counter.
  !             mnod  : mark of node.
  !             node  : another loop counter.
  !             pred  : predecessor (father) of node.
  !             snod  : current element on strn.
  !     ------------------------------------------------------------------
  
  
  integer :: arc   , enod  , i     , mnod  , node  , pred  , snod

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             dof   : dual objective function value on tree.
  !             pof   : primal objective function value on tree.
  !     ------------------------------------------------------------------
  
  double precision :: dof   , pof

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             zero  : value for zero.
  !     ------------------------------------------------------------------
  
  double precision :: zero

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Initializations.
  !     ------------------------------------------------------------------

  pf  = 0
  opt = 0
  pof = 0d0
  dof = 0d0
  
  !     ------------------------------------------------------------------
  !     Initialize flow variables x.
  !     ------------------------------------------------------------------

  do node = 1, root - 1
     x(node) =  b(node)
  end do

  do node = root + 1, nn
     x(node) =  b(node)
  end do
   
  x(root) = 0
  arc=ufirst
  do i = 1, uarcs
     snod = strn(arc)
     enod = endn(arc)
     x(snod) = x(snod)+u(arc)
     x(enod) = x(enod)-u(arc)
     arc = rsets(arc)
  end do
  x(root) = 0

  !     ------------------------------------------------------------------
  !     Compute flow variables x belonging to the set b by
  !     performing a reverse treath traversal of the spanning
  !     tree. This procedure computes also the total cost of
  !     the basic variables.
  !     ------------------------------------------------------------------
  
  pf = 1
  node = rtreth(root)
  
  do i = 2, nn
     arc  = trearc(node)
     snod = strn(arc)
     pred = father(node)
     x(pred) = x(pred)+x(node)
     if (node == snod) x(node)=-x(node)
     pof = pof + dble(x(node))*dble(c(arc))
     if (bstats(arc)  == 1) then
        if (x(node)< 0 .or. x(node) > u(arc)) then
           pf = 0
           return
        endif
     else
        if (x(node) < 0) then
           pf = 0
           return
        endif
     endif
     node = rtreth(node)     
  end do
  x(root) = 0

  !     ------------------------------------------------------------------
  !     Add the total cost of the flow variables at the upper bound.
  !     ------------------------------------------------------------------
  
  arc=ufirst
  do i = 1, uarcs
     snod = strn(arc)
     enod = endn(arc)
     pof = pof + dble(c(arc))*dble(u(arc))
     arc  = rsets(arc)
  end do

  !     ------------------------------------------------------------------
  !     Project y onto the dual face complementary to x.
  !     ------------------------------------------------------------------
  
  y0(root)   = 0d0
  mark(root) = root
  node = treath(root)

  do i = 2, nn
     arc = trearc(node)
     if (x(node) == 0 .or. &
          (bstats(arc) == 1 .and. x(node) == u(arc))) then
        y0(node)     = 0d0
        mark(node)   = node
        nnodes(node) = 1
        alfa(node)   = 0d0
     else
        pred         = father(node)
        arc          = trearc(node)
        y0(node)     = dble(c(arc))
        mark(node)   = mark(pred)
        mnod         = mark(node)
        nnodes(mnod) = nnodes(mnod)+1
     endif
     node = treath(node)
  end do

  node = treath(root)

  do i = 2, nn
     arc  = trearc(node)
     snod = strn(arc)
     enod = endn(arc)
     mnod = mark(node)
     if ( (bstats(arc) == 0 .and. x(node) > 0) .or. &
          (bstats(arc) == 1 .and. x(node) > 0 .and. &
          x(node) < u(arc))) then
        if (node == enod) then
           y0(enod) = y0(enod)+y0(snod)
        else
           y0(snod) = y0(enod)-y0(snod)
        endif
     endif
     alfa(mnod) = alfa(mnod)+y(node)-y0(node)
     node = treath(node) 
  end do

  do node = 1, root - 1
     mnod = mark(node)
     y0(node) = y0(node)+alfa(mnod)/dble(nnodes(mnod))
  end do

  do node = root + 1, nn
     mnod = mark(node)
     y0(node) = y0(node)+alfa(mnod)/dble(nnodes(mnod))
  end do

  !     ------------------------------------------------------------------
  !     Estimate dual objective function value.
  !     ------------------------------------------------------------------

  do arc = 1, na
     snod = strn(arc)
     enod = endn(arc)
     av4(arc) = y0(enod)-y0(snod)-dble(c(arc))     
  end do

  arc = npfrst

  do i = 1, nbound
     if (av4(arc) > zero) return
     arc = upnext(arc)
  end do

  do node = 1, root - 1
     dof = dof+dble(b(node))*y0(node)
  end do

  do node = root + 1, nn
     dof = dof+dble(b(node))*y0(node)
  end do
  
  arc = upfrst
  
  do i = 1, bound
     snod = strn(arc)
     enod = endn(arc)
     dof  = dof-dble(u(arc))*dmax1(av4(arc),0.0d0)
     arc  = upnext(arc)    
  end do

  if (dabs(pof-dof) < 1.0d0) then
     
     !        ----------------------------------------------------------------
     !        Dump the optimal flow to optflo and exit.
     !        ----------------------------------------------------------------     
     
     opt = 1
     do i = 1, na
        optflo(i) = 0
     end do
     do i = 1, nn
        if (trearc(i) > 0) optflo(trearc(i)) = x(i)
     end do
     return
  end if

  !     ------------------------------------------------------------------
  !     End of subroutine pdoptm.
  !     ------------------------------------------------------------------

end subroutine pdnet_optcheck
!     ------------------------------------------------------------------
!     pdnet_optimaldualface: 
!             Based on a tentative optimal dual solution, determines 
!             the corresponding optimal face. assumes the given dual 
!             solution is interior to the optimal face.
!     ------------------------------------------------------------------

subroutine pdnet_optimaldualface(actlst, cost,     endv,  maxslk, nact,&
     ncap, nedge,  nvrtx,    strv,  tolslk, ysol)

  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nact  : number of active edges.
  !             ncap  : number of dropped edges at capacity.
  !             nedge : number of edges.
  !             nvrtx : number of vertices.
  !     ------------------------------------------------------------------
  
  integer :: nact  , ncap  , nedge , nvrtx

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             actlst(1:    nedge): list of active edges.
  !             endv  (1:  nedge): end vertex in network data structure
  !             strv  (1:  nedge): start vertex in network data structure.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: actlst, endv, strv 

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             cost  (1:nedge): cost coefficients.
  !             ysol  (1:nvrtx): interior dual solution.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nedge) :: cost
  double precision, dimension(1:nvrtx) :: ysol

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             edg   : edge index.
  !             pact  : pointer to actlst.
  !             pcap  : pointer to actlst.
  !             plwb  : pointer to actlst.
  !     ------------------------------------------------------------------
  
  integer :: edg   , pact  , pcap  , plwb

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             maxslk: largest slack on dual.
  !             tolslk: zero for dual slacks.
  !     ------------------------------------------------------------------

  double precision :: maxslk, tolslk

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             nrm   : normalizing factor.
  !             tmp   : temporary variable.
  !     ------------------------------------------------------------------

  double precision :: nrm   , tmp

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Count number of edges in each class.
  !     ------------------------------------------------------------------
  
  maxslk=0.d0
  nact=0
  ncap=0

  do edg=1,nedge
     
     !        ---------------------------------------------------------------
     !        Compute normalizing factor.
     !        ---------------------------------------------------------------
     
     nrm=dabs(cost(edg))
     
     if(nrm < 1.d0) nrm=1.d0

     !        ---------------------------------------------------------------
     !        Compute reduced cost.
     !        ---------------------------------------------------------------
     
     tmp=cost(edg)-ysol(strv(edg))+ysol(endv(edg))
     
     !        ---------------------------------------------------------------
     !        Active edge.
     !        ---------------------------------------------------------------

     if(dabs(tmp) < nrm*tolslk) then
        nrm=dabs(tmp)/nrm
        if(maxslk < nrm) maxslk=nrm
        nact=nact+1
        
        !        ---------------------------------------------------------------
        !        Capacity ( relative cost < 0 ).
        !        ---------------------------------------------------------------
        
     else if (tmp < 0.0d0) then
        ncap = ncap + 1
     end if
  end do
  
  !     ------------------------------------------------------------------
  !     Assign edges to class.
  !     ------------------------------------------------------------------
  
  pact=1
  pcap=nact+1
  plwb=nact+ncap+1
  
  do edg=1,nedge
     
     !        ---------------------------------------------------------------
     !        Compute normalizing factor.
     !        ---------------------------------------------------------------
     
     nrm=dabs(cost(edg))
     if(nrm < 1.d0) nrm=1.d0
     
     !        ---------------------------------------------------------------
     !        Compute reduced cost.
     !        ---------------------------------------------------------------
     
     tmp=cost(edg)-ysol(strv(edg))+ysol(endv(edg))
     
     if(dabs(tmp) < nrm*tolslk) then
        actlst(pact)=edg
        pact=pact+1
        
        !           ------------------------------------------------------------
        !           Capacity ( aux < 0 ).
        !           ------------------------------------------------------------
        
     else if(tmp < 0d0) then
        actlst(pcap)=edg
        pcap=pcap+1
        
        !           ------------------------------------------------------------
        !           Lower bound ( aux > 0 ).
        !           ------------------------------------------------------------
        
     else
        actlst(plwb)=edg
        plwb=plwb+1
     end if
  end do

!     ------------------------------------------------------------------
!     End of subroutine pdnet_optimaldualface.
!     ------------------------------------------------------------------

end subroutine pdnet_optimaldualface
!     ------------------------------------------------------------------
!     pdnet_orderstruct: 
!             Given a spanning tree or forest, this subroutine finds
!             permutations for vertices and tree edges that yield an
!             incidence matrix with a block triangular structure.
!     ------------------------------------------------------------------

subroutine pdnet_orderstruct (endv, edglst, edgptr, nedge,  ntree,  &
     nvrtx, stack,  stedg,  stptr,  strv, stvrtx)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nedge : number of edges.
  !             ntree : number of trees.
  !             nvrtx : number of vertices.
  !     ------------------------------------------------------------------
  
  integer :: nedge , ntree , nvrtx

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endv  (1:  nedge): end vertex in network data structure
  !             edglst(1:2 * (nvrtx - ntree): list for edges in tree adjacent
  !                                           to each vertex.
  !             edgptr(1:    nvrtx + 1): pointer to beginning edge list for 
  !                                      each vertex.
  !             stack (1:    nvrtx    ): stack of eligible vertices in dfs
  !             stedg (1:    nedge    ): edge linking vertex to parent in
  !                                      a rooted tree representation.
  !             stptr (1:    ntree + 1): pointer to first vertex in
  !                                      a tree.
  !             strv  (1:  nedge): start vertex in network data structure.
  !             stvrtx(1:    nvrtx    ): list of vertices in spanning
  !                                      forest. Returns a triagular
  !                                      ordering of vertices.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) ::  endv, stedg, strv
  integer, dimension(1:2 * nvrtx -1) :: edglst
  integer, dimension(1:nvrtx + 1) :: edgptr
  integer, dimension(1:ntree + 1) :: stptr
  integer, dimension(1:nvrtx    ) :: stack, stvrtx

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             edg   : edge index.
  !             eptr  : temporary pointer to the st-edge data structure.
  !             ntedge: number of edges in spanning forest.
  !             nxtvtx: adjacent vertex.
  !             tedg  : index for tree edges.
  !             top   : poiter to top of stack.
  !             tree  : tree counter.
  !             vptr  : temporary pointer to the vertex data structure.
  !             vrtx  : vertex index.
  !     ------------------------------------------------------------------
  
  integer :: edg, eptr, ntedge, nxtvtx, tedg, top, tree, vptr, vrtx


  !     ------------------------------------------------------------------
  !     Initialize vertex degrees in edgptr array.  This array will be
  !     used as pointers to the edglst array, indicating the beginning of 
  !     incoming and outgoing edge list  for each vertex.  
  !     ------------------------------------------------------------------

  do vrtx = 1, nvrtx
     edgptr(vrtx)=0
  end do

  !     ------------------------------------------------------------------
  !     Compute degree of each edge.
  !     ------------------------------------------------------------------

  ntedge=nvrtx-ntree
  do tedg = 1, ntedge
     edg=stedg(tedg)
     vrtx=strv(edg)
     edgptr(vrtx)=edgptr(vrtx)+1
     vrtx=endv(edg)
     edgptr(vrtx)=edgptr(vrtx)+1
  end do

  !     ------------------------------------------------------------------
  !     Set edge list pointers after the end of list (+1) for each vertex.
  !     ------------------------------------------------------------------

  edgptr(1)=edgptr(1)+1
  do vrtx = 2, nvrtx
     edgptr(vrtx)=edgptr(vrtx)+edgptr(vrtx-1)
  end do

  !     ------------------------------------------------------------------
  !     Build lists of edges using stored pointers.  After this operation,
  !     the pointers indicate the beginning of the edge list for each 
  !     vertex.
  !     ------------------------------------------------------------------

  do tedg = ntedge, 1, -1
     edg=stedg(tedg)
     vrtx=strv(edg)
     edgptr(vrtx)=edgptr(vrtx)-1
     edglst(edgptr(vrtx))=edg
     vrtx=endv(edg)
     edgptr(vrtx)=edgptr(vrtx)-1
     edglst(edgptr(vrtx))=edg
  end do

  edgptr(nvrtx+1)=2*ntedge+1

  !     ------------------------------------------------------------------
  !     Initialize parent edges.
  !     ------------------------------------------------------------------

  do vrtx = 1, nvrtx
     stedg(vrtx)=0
  end do

  !     ------------------------------------------------------------------
  !     Traverse each tree in depth first order saving dfs rank and
  !     parent edge.
  !     ------------------------------------------------------------------

  do tree = 1, ntree
     
     !        ---------------------------------------------------------------
     !        Initialize current vertex with tree root (parent edge 0).
     !        ---------------------------------------------------------------
       
     vptr=stptr(tree)
     vrtx=stvrtx(vptr)

     !        ---------------------------------------------------------------
     !        Initialize stack of eligible vertices.
     !        ---------------------------------------------------------------

     top=1
     stack(top)=vrtx

     !        ---------------------------------------------------------------
     !        Position pointer to first edge in list for vrtx.
     !        ---------------------------------------------------------------
     
     eptr=edgptr(vrtx)

     !        ---------------------------------------------------------------
     !        Start Do until with empty.stack
     !        ---------------------------------------------------------------

     do
        
        !        ---------------------------------------------------------------
        !        Test for end of edge list.
        !        ---------------------------------------------------------------
        
        if(eptr >= edgptr(vrtx+1)) then
           top=top-1
           
           !           ------------------------------------------------------------
           !           Skip to next tree if stack is empty.
           !           ------------------------------------------------------------
           
           if (top <= 0) exit
           
           !           ------------------------------------------------------------
           !           Identify vertex and saved pointer to edge list.
           !           ------------------------------------------------------------
           
           vrtx=stack(top)
           eptr=edglst(edgptr(vrtx)) 
           cycle
        end if
        
        !           ------------------------------------------------------------
        !           Identify next successor edge.
        !           ------------------------------------------------------------
        
        edg=edglst(eptr)

        !           ------------------------------------------------------------
        !           Advance and save edge pointer for current vertex.
        !           ------------------------------------------------------------
        
        eptr=eptr+1
        edglst(edgptr(vrtx))=eptr     
        
        !           ------------------------------------------------------------
        !           Skip if edge is parent.
        !           ------------------------------------------------------------
        
        if(edg == stedg(vrtx)) cycle

        !           ------------------------------------------------------------
        !           Identify next vertex (edg=(vrtx,nxtvtx) in rooted tree).
        !           ------------------------------------------------------------
        
        nxtvtx=strv(edg)
        if(nxtvtx == vrtx) nxtvtx=endv(edg)

        !           ------------------------------------------------------------
        !           Save dfs order and parent edge.
        !           ------------------------------------------------------------
        
        vptr=vptr+1
        stvrtx(vptr)=nxtvtx
        stedg(nxtvtx)=edg

        !           ------------------------------------------------------------
        !           Add visited vertex to stack of eligible vertices.
        !           ------------------------------------------------------------
        
        top=top+1
        vrtx=nxtvtx
        stack(top)=vrtx

        !           ------------------------------------------------------------
        !           Position pointer to first edge in list for vrtx.
        !           ------------------------------------------------------------
        
        eptr=edgptr(vrtx)

     end do

     !        ---------------------------------------------------------------
     !        End of "do until".
     !        ---------------------------------------------------------------
     
     
  end do
  
  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_orderstruct.
  !     ------------------------------------------------------------------

end subroutine pdnet_orderstruct


!     ------------------------------------------------------------------
!     pdnet_perturb: This subroutine perturbs and transforms the data to
!                    double precision.
!     ------------------------------------------------------------------


subroutine pdnet_perturb (b,      bound,  c,      db,     dc,     du,&
     endn,   na,     nn,     s,      strn,   u,  upfrst, upnext, w, &
     x,      z )

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables.
  !
  !                   bound : number of bounded arcs.
  !                   na    : number of arcs.
  !                   nn    : number of nodes.
  !                   upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------

  integer :: bound , na    , nn    , upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             b     (1:nn): array of capacities at each node.
  !             c     (1:na): array of costs at each arc.
  !             endn  (1:na): pointer of network data structure.
  !             strn  (1:na): adjacency list of network data structure.
  !             u     (1:na): upper bound of flow at each arc.
  !             upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: b
  integer, dimension(1:na) :: c, endn, strn, u, upnext

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             db    (1:nn): double precision copy of int vector b.
  !             dc    (1:na): double precision copy of int vector c.
  !             du    (1:na): double precision copy of int vector u.
  !             s     (1:na): dual iterate s.
  !             w     (1:na): dual iterate w.
  !             x     (1:na): primal iterate x.
  !             z     (1:na): primal iterate z.
  !     ------------------------------------------------------------------

  double precision, dimension (1:nn) :: db
  double precision, dimension (1:na) :: dc, du, s, w, x, z

  !     ------------------------------------------------------------------
  !     Integer variables.
  !
  !                   arc   : current arc.
  !                   enod  : current element on endn.
  !                   i     : loop counter.
  !                   node  : another loop counter.
  !                   snod  : current element on strn.
  !     ------------------------------------------------------------------

  integer :: arc   , enod  , i     , node  , snod

  !     ------------------------------------------------------------------
  !     Start of executable section of the subroutine.
  !     ------------------------------------------------------------------

  do arc = 1, na
     x(arc) = 0.0d0
     s(arc) = 0.0d0
     z(arc) = 0.0d0
     w(arc) = 0.0d0
  end do

  do node = 1, nn
     db(node) = dble(b(node))
  end do

  do arc = 1, na
     snod = strn(arc)
     enod = endn(arc)
     db(snod) = db(snod)+x(arc)
     db(enod) = db(enod)-x(arc)
     dc(arc) = dble(c(arc))+z(arc)-w(arc)    
  end do

  arc = upfrst

  do i = 1, bound
     du(arc)  = dble(u(arc))+x(arc)+s(arc)
     arc = upnext(arc)     
  end do

!     ------------------------------------------------------------------
!     End of subroutine pdnet_perturb.
!     ------------------------------------------------------------------

end subroutine pdnet_perturb
!     ------------------------------------------------------------------
!     pdnet_precconjgrd: 
!             This sobroutine solves the system a*teta*at x = b by means
!             of a preconditioned conjugate gradient (pcg) algorithm.
!     ------------------------------------------------------------------

subroutine pdnet_precconjgrd(b,      cgstop, diag,   endn,   err,    &
     father,iprcnd, mxcgit, na,     ndiag,  nn,     p, pcgit,  q,    &
     r,      res,    root,   rtreth, strn,   sw1,    teta,   tolcg,  &
     trearc, treath,  x,      z)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             cgstop: flags the way PCG stopped.
  !             iprcnd: preconditioner used.
  !             mxcgit: maximum number of PCG iterations.
  !             na    : number of arcs.
  !             nn    : number of nodes.
  !             pcgit : number of PCG iterations performed.
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------
  
  integer :: cgstop, iprcnd, mxcgit, na    , nn    , pcgit  , root

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endn  (1:na): pointer of network data structure.
  !             father(1:nn): marks the father of each node.
  !             rtreth(1:nn): adjacent to a given node.
  !             strn  (1:na): adjacency list of network data structure.
  !             trearc(1:nn): arc of tree corresponding to node.
  !             treath(1:nn): node that is adjacent to given.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: father, rtreth, trearc, treath
  integer, dimension(1:na) :: endn, strn

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             b     (1:nn): array of capacities at each node.
  !             diag  (1:nn): diagonal of matrix.
  !             ndiag (1:nn): diagonal of preconditioner.
  !             p     (1:nn): vector p.
  !             q     (1:nn): vector q.
  !             r     (1:nn): vector r.
  !             teta  (1:na): diagonal theta.
  !             x     (1:nn): primal solution x.
  !             z     (1:nn): dual solution z.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: b, diag, ndiag, p, q, r, x, z
  double precision, dimension(1:na) :: teta

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             err   : value of cosine on PCG iteration.
  !             res   : value of residual on PCG iteration.
  !             tolcg : tolerance for PCG stopping criterion.
  !     ------------------------------------------------------------------

  double precision :: err   , res   , tolcg
  
  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             alfa  : parameter from PCG algorithm.
  !             beta  : parameter from PCG algorithm.
  !             dnrm2 : 2-norm of a vector (BLAS function).
  !             eps   : temporary storage for inner product.
  !             nrmb  : 2-norm of vector b.
  !             oldeps: value of eps in previous iteration.
  !             sw1   : copy of maximum iteration count.
  !     ------------------------------------------------------------------

  double precision :: alfa  , beta  , dnrm2 , eps   , nrmb  , oldeps, sw1

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !
  !     Initializations.
  !     ------------------------------------------------------------------
  
  pcgit = 0


  call pdnet_attx(endn, na, nn, root, strn, teta, x, r)

  call dcopy (nn,     r,      1,      z,      1)
  call pdnet_xminusy( nn,     root,   b,      z,      r)
  
  call pdnet_iqrdsolv(r, diag, father, iprcnd, na, ndiag, nn, root, rtreth,&
       teta, trearc, treath, z)

  call pdnet_xequaly(nn,     root,   p,      z)

  call pdnet_ddot ( eps,    nn,     root,   r,      z)

  
  call pdnet_l1norm (res,    nn,     r)

  !     ------------------------------------------------------------------
  !     PCG algorithm.
  !     ------------------------------------------------------------------

  nrmb = dnrm2(nn, b, 1)

  do
     pcgit = pcgit + 1
     if (pcgit >= mxcgit) then
        cgstop=0
        return
     end if

     if ((iprcnd == 1).and.(dble(pcgit) > sw1)) then
        cgstop=1
        return
     end if

     call pdnet_attx(endn, na, nn, root, strn, teta, p, q)
     call pdnet_ddot( alfa,   nn,     root,   p,      q )
     alfa = eps / alfa
     call pdnet_update1( alfa,   nn,     root,   x,      p)


     alfa = -alfa
     call pdnet_update1( alfa,   nn,     root,   r,      q )

     call pdnet_iqrdsolv( r, diag, father, iprcnd, na, ndiag, nn, root, &
          rtreth, teta,   trearc, treath, z)
     oldeps = eps
     call pdnet_ddot(eps,    nn,     root,   r,      z )
     beta = eps / oldeps
     call pdnet_update2( beta,   nn,     root,   p,      z)
     call pdnet_l1norm(res,    nn,     r)

     !     ------------------------------------------------------------------
     !     Cosine stopping criterion.
     !     ------------------------------------------------------------------
     

     call pdnet_approxcosine(b, err, nn, nrmb, q, r, root)

     if (err < tolcg) then
        cgstop = 2
        return
     end if
  


  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_precconjgrd.
  !     ------------------------------------------------------------------

end subroutine pdnet_precconjgrd
!     ------------------------------------------------------------------
!     pdnet_printout: 
!             This subroutine prints summary report at the end of each
!             iteration.
!     ------------------------------------------------------------------

subroutine pdnet_printout(cgstop, cgtol,  dobj,   dof,    gap,    iout,&
     iprcnd, itr,    mfdopt, mfstat, miu,    mxfe, mxfv,   na,   &
     opt,    pcgcos, pcgitr, pcgres, pcgtol, pf,     pobj,   pof,  &
     prttyp, tol1,   tol2,   x)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             cgstop: flags the way PCG stopped.
  !             iout  : output file unit.
  !             iprcnd: preconditioner used.
  !             itr   : number of primal-dual iterations performed.
  !             mfstat: maxflow status.
  !             mxfe  : edges on optimal flow network.
  !             mxfv  : vertices on optimal flow network.
  !             na    : number of arcs.
  !             opt   : flags optimal flow.
  !             pcgitr: number of PCG iterations performed.
  !             pf    : iteration level.
  !             prttyp: printing level.
  !     ------------------------------------------------------------------
     
  integer :: cgstop, iout, iprcnd, itr, mfstat, mxfe, mxfv, na, &
       opt, pcgitr, pf, prttyp

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             x(1:na): array with optimal flow at each arc.
  !     ------------------------------------------------------------------
  
  integer, dimension(1:na) :: x

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !     ------------------------------------------------------------------
  
  integer :: i

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             cgtol : value of cosine on previous PCG iteration.
  !             dobj  : dual objective function value.
  !             dof   : dual objective function value on tree.
  !             gap   : value of the gap.
  !             mfdopt: maxflow value.
  !             miu   : value of interior-point mu parameter.
  !             pcgcos: value of cosine on PCG iteration.
  !             pcgres: value of residual on PCG iteration
  !             pcgtol: value of residual on previous PCG iteration.
  !             pobj  : objective function value.
  !             pof   : primal objective function value on tree.
  !             tol1  : lower bound of the tolerance for maxflow.
  !             tol2  : upper bound of the tolerance for maxflow.
  !     ------------------------------------------------------------------
  
  double precision :: cgtol, dobj, dof, gap, mfdopt, miu, pcgcos, pcgres, &
       pcgtol, pobj, pof, tol1, tol2

  !     ------------------------------------------------------------------
  !     String type variables:
  !
  !             ccgstp: type of stopping on CG.
  !             cmfstt: type of optimality state.
  !             precnd: type of preconditioner.
  !     ------------------------------------------------------------------

  character(len= 8) :: ccgstp
  character(len=10) :: cmfstt
  character(len=19) :: precnd

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------
  
  if (prttyp == 0) return
  if (cgstop == 0) then
     ccgstp = 'MAX'
  else
     if (cgstop == 1) then
        ccgstp = 'RESIDUAL'
     else
        if (cgstop == 2) then
           ccgstp = 'COSINE'
        else
           ccgstp ='??????'
        endif
     endif
  endif
  if (mfstat == 0) then
     cmfstt = 'INACTIVE'
  else
     if (mfstat == 1) then
        cmfstt = 'SUBOPTIMAL'
     else
        if (mfstat == 2) then
           cmfstt = 'OPTIMAL'
        else
               cmfstt = '?'
            endif
         endif
      endif
      if (iprcnd == 1) then
         precnd='DIAGONAL'
      else
         if (iprcnd == 3) then
            precnd='SPANNING TREE'
         else
            if (iprcnd == 2) then
               precnd='IQRD'
            else
               if (iprcnd == 4) then
                  precnd='INCOMPLETE CHOLESKY'
               else
                  precnd='?'
               endif
            endif
         endif
      endif

      write (iout,899)
899   format('----------------------------------------',&
           '----------------------------------------')
      
      if (prttyp >= 2) then
         write (iout,902) itr,miu
902      format('  CONV>     itr: ',i14,'       mu: ',e14.8,/)
         write (iout,901) pobj,dobj,gap
901      format('  COST>    pobj: ',e14.8,'     dobj: ',e14.8,&
              '  pd-gap: ',e14.8,/)
         write (iout,903) precnd,ccgstp,pcgitr,pcgtol,pcgres,cgtol,pcgcos
903      format('    CG>  precnd: ',a4,7x,'      cgstop: ',a4,3x,&
              '         pcgitr: ',i14,/,&
              '    CG>  pcgtol: ',e14.8,'   pcgres: ',e14.8,/,&
              '    CG>   cgtol: ',e14.8,'   pcgcos: ',e14.8,/)
         if (pf == 1 .and. opt == 1) then
            write (iout,904)
904         format('  TREE> pb-feas: TRUE             pb-opt: TRUE')
         endif
         if (pf == 1 .and. opt == 0) then
            write (iout,905)
905         format('  TREE> pb-feas: TRUE             pb-opt: FALSE')
         endif
         if (pf == 0) then
            write (iout,906)
906         format('  TREE> pb-feas: FALSE')
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
907      format('  TREE>    pobj: ',e14.8,'     dobj: ',e14.8)
      endif
      write (iout,911)
      
      if (miu < 1.0d0) then
         if (prttyp >= 2) then
            write (iout,908) mxfv,mxfe
908         format(' MFLOW>    mxfv: ',i14,'     mxfe: ',i14)
         endif
      endif
      if (prttyp >= 2) then
         if (miu < 1.0d0) then
            write (iout,909) tol1,tol2,mfdopt
909         format(' MFLOW> pdlotol: ',e14.8,'  pdhitol: ',e14.8,&
                 '    dobj: ',e14.8)
         endif
      endif
      write (iout,910) cmfstt
910   format(' MFLOW>  mfstat: ',a10)
      if (miu < 1.0d0) then
         write (iout,911)
911      format(' ')
      endif
      if (prttyp == 3) then
         open (14, file = 'solution.dat', status = 'unknown',&
              access = 'sequential', form = 'formatted')
         do  i = 1, na
            write (14,*) x(i)
         end do
         close(14)
      endif

      !     ------------------------------------------------------------------
      !     End of subroutine pdnet_printout.
      !     ------------------------------------------------------------------

    end subroutine pdnet_printout


!     ---------------------------------------------------------------------
!     pdnet_probdata: This subroutine finds the data characteristics of the 
!                     problem.
!     ---------------------------------------------------------------------

subroutine pdnet_probdata (b, bound,  c, huge, maxb, maxc,maxu,   minb,   &
     minc,   minu,   na,     nn,norb,   norc,   noru,   u,      upfrst, upnext)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !               bound : number of bounded arcs.
  !               na    : number of arcs.
  !               nn    : number of nodes.
  !               upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------
  
  integer :: bound , na    , nn    , upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !               upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------  

  integer, dimension(1:na) :: upnext

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !               b(1:nn): array of capacities at each node.
  !               c(1:na): array of costs at each arc.
  !               u(1:na): upper bound of flow at each arc.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: b
  double precision, dimension(1:na) :: c, u

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !               arc   : current arc.
  !               i     : loop counter.
  !               node  : another loop counter.
  !     ------------------------------------------------------------------

  integer :: arc   , i     ,  node

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !               huge  : maxdouble precision number.
  !               maxb  : largest absolute value of b.
  !               maxc  : largest absolute value of c.
  !               maxu  : largest absolute value of u.
  !               minb  : smallest absolute value of b.
  !               minc  : smallest absolute value of c.
  !               minu  : smallest absolute value of u.
  !               norb  : 1-norm of vector b.
  !               norc  : 1-norm of vector c.
  !               noru  : 1-norm of vector u.
  !     ------------------------------------------------------------------  

  double precision ::  huge  , maxb  , maxc  , maxu  , minb  , minc  , minu  , &
       norb  , norc  , noru

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------  

  maxc = -huge
  maxu = -huge
  maxb = -huge
  minc =  huge
  minu =  huge
  minb =  huge
  norc = 0d0
  noru = 0d0
  norb = 0d0

  do arc = 1, na
     maxc = dmax1(maxc,c(arc))
     minc = dmin1(minc,c(arc))
     norc = norc+dabs(c(arc))
  end do

  arc = upfrst

  do i = 1, bound
     maxu = dmax1(maxu,u(arc))
     minu = dmin1(minu,u(arc))
     noru = noru+dabs(u(arc))
     arc = upnext(arc)     
  end do

  do node = 1, nn
     maxb = dmax1(maxb,b(node))
     minb = dmin1(minb,b(node))
     norb = norb+dabs(b(node))     
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_probdata.
  !     ------------------------------------------------------------------

end subroutine pdnet_probdata
!     ------------------------------------------------------------------
!     pdginf: Saves values of all parameters and statistics to arrays
!             info() and dinfo().
!     ------------------------------------------------------------------
subroutine pdnet_putinfo (bound , cgstop, dinfo,  dobj,   dof, factor, fcttol, &
     gap, grho0,  gtolcg, gtolmf, huge,   ierror, info,   ioflag,  iout, &
     iprcnd, maxit, mfdopt, mfstat, miu, mntlcg, mr, mrho0, mxcgit, mxfe, &
     mxfv, nbound, nbuk,  oldcos, oldtol, opt, optgap, optmf, optst , pcgcos, &
     pcgit,  pcgres, pduit,  pf, pobj, pof, prttyp, rho0,  root, s1fctr,stptol, &
     stpval, tol1, tol2, tolcg, tolcg0, tolpcg, tolslk, tolsw1, tolsw2, sw1, sw2, zero)

  implicit none
  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             info  (1:30): array with integer statistics and parameters.
  !     ------------------------------------------------------------------

  integer, dimension (1:30) :: info

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             bound : used for computing aditional data structures.
  !             cgstop: flags the way PCG stopped.
  !             ierror: flags errors in input.
  !             ioflag: controls the level of output of summary.
  !             iout  : output file unit.
  !             iprcnd: preconditioner used.
  !             maxit : maximum number of primal-dual iterations.
  !             mfstat: maxflow status.
  !             mr    : estimate of maximum number of iterations.
  !             mxcgit: maximum number of PCG iterations.
  !             mxfe  : edges on optimal flow network.
  !             mxfv  : vertices on optimal flow network.
  !             nbound: used for computing aditional data structures.
  !             nbuk  : number of buckets.
  !             opt   : flags optimal flow.
  !             optgap: Duality gap optimality indicator flag.
  !             optmf : Maximum flow optimality indicator flag.
  !             optst : Spanning tree optimality indicator flag.
  !             pcgit : number of PCG iterations performed.
  !             pduit : number of primal-dual iterations performed.
  !             pf    : iteration level.
  !             prttyp: printing level.
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------

  integer ::  bound , cgstop, ierror, ioflag, iout, iprcnd, &
       maxit , mfstat, mr    , mxcgit, mxfe  , mxfv  , &
       nbound, nbuk  , opt   , optgap, optmf , optst , pcgit, &
       pduit ,  pf    , prttyp, root 

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             dinfo  (1:35): array with double precision statistics and 
  !                            parameters.
  !     ------------------------------------------------------------------

  double precision, dimension(1:35) :: dinfo
  
  !     ----------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             dobj  : dual objective function value.
  !             dof   : dual objective function value on tree.
  !             factor: factor for Newton step.
  !             fcttol: factorization zero tolerance.
  !             gap   : value of the gap.
  !             grho0 : factor for rho0 update.
  !             gtolmf: update factor for tol1 and tol2.
  !             huge  : maxdouble precision number.
  !             mfdopt: maxflow value.
  !             miu   : value of interior-point mu parameter.
  !             mntlcg: lower bound for CG tolerance.
  !             mrho0 : lower bound for rho0.
  !             oldcos: value of cosine on previous PCG iteration.
  !             oldtol: value of residual on previous PCG iteration.
  !             pcgcos: value of cosine on PCG iteration.
  !             pcgres: value of residual on PCG iteration.
  !             pobj  : objective function value.
  !             pof   : primal objective function value on tree.
  !             rho0  : parameter from IQRD preconditioner.
  !             s1fctr: factor for miu.
  !             stptol: zero for dual slacks.
  !             stpval: largest slack on dual face.
  !             sw1   : iteration limit on CG for preconditioner switch.
  !             sw2   : another limit for switch.
  !             tol1  : lower bound of the tolerance for maxflow.
  !             tol2  : upper bound of the tolerance for maxflow.
  !             tolcg : tolerance for PCG stopping criterion.
  !             tolcg0: guess for tolerance for PCG stopping criterion.
  !             tolpcg: initial value for tolcg.
  !             zero  : value for zero.
  !     ------------------------------------------------------------------

  double precision :: dobj  , dof   , factor,  fcttol, gap   , grho0,&
       gtolcg, gtolmf, huge  ,  mfdopt,&
       miu   , mntlcg, mrho0 ,   oldcos, oldtol, &
       pcgcos, pcgres, pobj  ,  pof   , rho0  , s1fctr, stptol,&
       stpval, sw1   , sw2   , tol1  ,  tol2  , tolcg , tolcg0, tolpcg, tolslk,&
       tolsw1, tolsw2,  zero
  

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------
  !
  !     ------------------------------------------------------------------
  !     Integer statistics assignment.
  !     ------------------------------------------------------------------
  info( 1) = bound  
  info( 2) = iprcnd 
  info( 3) = iout   
  info( 4) = maxit  
  info( 5) = mxcgit 
  info( 6) = mxfv   
  info( 7) = mxfe   
  info( 8) = nbound 
  info( 9) = nbuk   
  info(10) = pduit  
  info(11) = root   
  info(12) = opt    
  info(13) = pcgit  
  info(14) = pf     
  info(15) = ierror 
  info(16) = cgstop 
  info(17) = mfstat 
  info(18) = ioflag 
  info(19) = prttyp 
  info(20) = optgap 
  info(21) = optmf  
  info(22) = optst  
  info(23) = mr     
  
  !     ------------------------------------------------------------------
  !     Double precision statistics assignment.
  !     ------------------------------------------------------------------
  dinfo( 1) = dobj   
  dinfo( 2) = dof    
  dinfo( 3) = factor 
  dinfo( 4) = fcttol 
  dinfo( 5) = gap    
  dinfo( 6) = grho0  
  dinfo( 7) = gtolcg 
  dinfo( 8) = gtolmf 
  dinfo( 9) = mfdopt 
  dinfo(10) = miu    
  dinfo(11) = mntlcg 
  dinfo(12) = mrho0  
  dinfo(13) = oldtol 
  dinfo(14) = pcgres 
  dinfo(15) = pobj   
  dinfo(16) = pof    
  dinfo(17) = rho0   
  dinfo(18) = s1fctr 
  dinfo(19) = stptol 
  dinfo(20) = stpval 
  dinfo(21) = tol1   
  dinfo(22) = tol2   
  dinfo(23) = tolcg  
  dinfo(24) = tolcg0 
  dinfo(25) = tolslk 
  dinfo(26) = huge   
  dinfo(27) = tolpcg 
  dinfo(28) = oldcos 
  dinfo(29) = pcgcos 
  dinfo(30) = zero   
  dinfo(31) = tolsw1 
  dinfo(32) = tolsw2 
  dinfo(33) = sw1    
  dinfo(34) = sw2    
  
  
  !    ------------------------------------------------------------------
  !     End of subroutine pdnet_putinfo.
  !     ------------------------------------------------------------------
end subroutine pdnet_putinfo
!-----------------------------------------------------------------------
! pdnet_setintparm: Sets values of external parameters stored in info().
!-----------------------------------------------------------------------

subroutine pdnet_setintparm(info)

  implicit none

  integer :: i
  integer, dimension(1:30) :: info

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  do i = 1, 30
     info(i) = 0
  enddo

  info( 2) = 3
  info( 4) = 1000
  info( 5) = 1000
  info(19) = 0
  info(20) = 1
  info(21) = 1
  info(22) = 1
  info( 3) = 6
  
  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_setintparm.
  !     ------------------------------------------------------------------
end subroutine pdnet_setintparm

  
!     ------------------------------------------------------------------
!     pdnet_sortedges: This subroutine sorts edges according to weights.
!     ------------------------------------------------------------------

subroutine pdnet_sortedges (actlst, buklst, nedge, nact, nbuk, ordedg, ptr, weight)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nact  : number of active edges.
  !             nbuk  : number of buckets for bucket sorting.
  !             nedge : number of edges.
  !     ------------------------------------------------------------------
  
  integer :: nact  , nbuk, nedge
  
  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             actlst(1:nact): list of active edges.
  !             ordedg(1:nact): order of edges for spanning tree.
  !             ptr   (1:nbuk): temporary array for pointers.
  !             buklst(1:nact): bucket for each edge.
  !     ------------------------------------------------------------------

  integer, dimension(1:nact) :: actlst, ordedg, buklst
  integer, dimension(1:nbuk) :: ptr

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !            weight(1:nact): weight vector for active edges.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nedge) :: weight

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             act   : active edge counter.
  !             buk   : bucket index.
  !             edg   : edge index.
  !     ------------------------------------------------------------------
  
  integer :: act   , buk   , edg

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             const : constant term.
  !             nrm   : normalizing factor.
  !             minw  : minimum weight.
  !             maxw  : maximum weight.
  !     ------------------------------------------------------------------

  double precision ::  const ,nrm   ,minw  , maxw, auxwe

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !
  !     Compute min and max weights.
  !     ------------------------------------------------------------------

  minw=1.d20
  maxw=-1.d20
  do act = 1, nact
     edg=actlst(act)
     auxwe=weight(edg)
     if(auxwe > maxw) then
        maxw=auxwe
     else if(auxwe < minw) then
        minw=auxwe
     end if
  end do
  
  minw=dlog(minw)
  maxw=dlog(maxw)

  !     ------------------------------------------------------------------
  !     Initialize head and tail for each bucket.
  !     ------------------------------------------------------------------

  do buk = 1, nbuk
     ptr(buk) = 0
  end do

  !   ------------------------------------------------------------------
  !   Count number of edges in each bucket.
  !   ------------------------------------------------------------------ 

  const=dble(nbuk)*(1.d0+minw/(maxw-minw))
  nrm=dble(nbuk)/(maxw-minw)

  do act = 1, nact
     edg=actlst(act)
     buk=int(const-dlog(weight(edg))*nrm)
     if(buk < 1) then
        buk=1
     else if(buk > nbuk) then
        buk=nbuk
     endif
     ptr(buk)=ptr(buk)+1
     buklst(act)=buk
  end do

  !     ------------------------------------------------------------------
  !     Compute pointers to ordered array.
  !     ------------------------------------------------------------------
  
  ptr(nbuk)=nact-ptr(nbuk)+1

  do buk = nbuk-1, 2, -1
     ptr(buk)=ptr(buk+1)-ptr(buk)
  end do
  ptr(1) = 1

  !     ------------------------------------------------------------------
  !     Build array with ordered edges.
  !     ------------------------------------------------------------------

  do act = 1, nact
     edg=actlst(act)
     buk=buklst(act)
     ordedg(ptr(buk))=edg
     ptr(buk)=ptr(buk)+1
  end do

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_sortedges.
  !     ------------------------------------------------------------------
  
end subroutine pdnet_sortedges
!     ------------------------------------------------------------------
!     pdnet_split: This subroutine splits a tree into two subtrees in a
!             Fibonacci heap by updating the arrays pred, rank, llink,
!             rlink, first and lost.
!     ------------------------------------------------------------------

subroutine pdnet_split (father, first, llink, lost, nn, pred, rank, rlink, soon)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             father: father of node.
  !             nn    : number of nodes.
  !             soon  : son of node.
  !     ------------------------------------------------------------------
  
  integer :: father, nn    , soon

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             first (1:nn): first node.
  !             llink (0:nn): part of the node left linked list.
  !             lost  (1:nn): number of lost nodes per node.
  !             pred  (1:nn): list of predecessors.
  !             rank  (1:nn): rank of nodes.
  !             rlink (0:nn): part of the node right linked list.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: first, lost, pred, rank
  integer, dimension(0:nn) :: llink, rlink

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             lsoon : left son.
  !             rsoon : right son.
  !     ------------------------------------------------------------------
  
  integer :: lsoon , rsoon

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------
      
  lsoon = llink(soon)
  rsoon = rlink(soon)
  rlink(soon) = 0
  llink(soon) = 0
  rlink(lsoon) = rsoon
  llink(rsoon) = lsoon
  pred(soon) = 0
  rank(father) = rank(father)-1

  if (lsoon == 0) then
     first(father)=rsoon
  end if

  lost(soon) = 0

  if (pred(father) /= 0) then
     lost(father)=lost(father)+1
  end if

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_split.
  !     ------------------------------------------------------------------
  
end subroutine pdnet_split
  
!     -----------------------------------------------------------------------
!     pdnet_startvalues: This subroutine computes the starting values for the
!                        primal-dual interior point algorithm.
!     -----------------------------------------------------------------------

subroutine pdnet_startvalues ( av1,    av2,    bound , db,     dc,     du, &
     endn,   gap,    maxb,   maxc,   minb,   miu,    na,     nbound, nn,   &
     norb,   npfrst, root,  s,      s1fctr, strn,   tolpcg, upfrst, upnext, &
     w,      x,      y,      z )

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             bound : number of bounded arcs.
  !             na    : number of arcs.
  !             nbound: number of unbounded arcs.
  !             nn    : number of nodes.
  !             npfrst: next bounded arc to be considered.
  !             root  : sets beginning of structure.
  !             upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------
  
  integer :: bound , na    , nbound, nn    , npfrst, root  ,  upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !          endn  (1:na): pointer of network data structure.
  !          strn  (1:na): adjacency list of network data structure.
  !          upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: endn, strn, upnext

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !           av1(1:na): temporary array 1.
  !           av2(1:nn): temporary array 2.
  !           db (1:nn): array of capacities at each node.
  !           dc (1:na): array of costs at each arc.
  !           du (1:na): upper bound of flow at each arc.
  !           s  (1:na): dual iterate s.
  !           w  (1:na): dual iterate w.
  !           x  (1:na): primal iterate x.
  !           y  (1:nn): dual iterate y.
  !           z  (1:na): primal iterate z.
  !     ------------------------------------------------------------------  

  double precision, dimension(1:nn) :: av2, db, y

  double precision, dimension(1:na) :: av1, dc, du, s, w, x, z

!     ------------------------------------------------------------------
!     Integer variables:
!
!             arc   : current arc.
!             enod  : current element on endn.
!             i     : loop counter.
!             node  : another loop counter.
!             snod  : current element on strn.
!     ------------------------------------------------------------------
 
  integer :: arc   , enod  , i     , node  , snod

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             gap   : value of the gap.
  !             maxb  : largest absolute value of b.
  !             maxc  : largest absolute value of c.
  !             minb  : smallest absolute value of b.
  !             miu   : value of interior-point mu parameter.
  !             norb  : 1-norm of vector b.
  !             s1fctr: factor for miu.
  !             tolpcg: initial value for tolcg.
  !     ------------------------------------------------------------------

  double precision :: gap, maxb, maxc, minb, miu, norb, s1fctr, tolpcg

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !              adema : initial demand.
  !              aux1  : auxiliary variable 1.
  !              aux2  : auxiliary variable 2.
  !              demand: demand.
  !              nav1  : small quantity.
  !     ------------------------------------------------------------------

  double precision :: adema , aux1  , aux2  , demand, nav1

  !     ------------------------------------------------------------------
  !     External functions:
  !
  !             dnrm2  : euclidean norm of a vector (BLAS).
  !     ------------------------------------------------------------------
  
  double precision :: dnrm2

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------ 

  adema  = dmax1(dabs(minb),dabs(maxb))
  demand = norb/2d0  

  !     ------------------------------------------------------------------
  !     Computing auxilliary vector av2.
  !     ------------------------------------------------------------------

  do node = 1, nn
     av2(node) = db(node)
  end do

  !     ------------------------------------------------------------------
  !     Computing initial vector y.
  !     ------------------------------------------------------------------ 

  aux1 = maxc/adema
  do node = 1, root - 1
     y(node) = aux1 * db(node)
  end do

  do node = root + 1, nn
     y(node) = aux1 * db(node)
  end do

  y(root) = 0.0d0

  !     ------------------------------------------------------------------
  !     Computing auxilliary vector av1.
  !     ------------------------------------------------------------------ 

  do arc = 1, na
     snod = strn(arc)
     enod = endn(arc)
     av1(arc) = y(snod)-y(enod)+dc(arc)
  end do

  !     ------------------------------------------------------------------
  !     Computing initial miu.
  !     ------------------------------------------------------------------

  miu      = 0d0
  arc      = npfrst  
  
  do i = 1, nbound
     aux1  = dabs(av1(arc))
     miu   = dmax1(miu,aux1)
     arc   = upnext(arc)    
  end do

  miu      = demand*miu
  arc      = upfrst

  do i = 1, bound
     aux1  = dabs(av1(arc)*du(arc))
     miu   = dmax1(miu,aux1)
     arc   = upnext(arc)     
  end do

  miu = s1fctr*miu

  !     ------------------------------------------------------------------
  !     Computing initial values for variables with upper bound.
  !     ------------------------------------------------------------------

  nav1=1.0d-8 * dnrm2 ( na, av1, 1)
  arc = upfrst  

  do i = 1, bound
     aux1 = du(arc)/2d0
     if (av1(arc) > nav1) then
        aux2   = miu/av1(arc)
        x(arc) = aux1+aux2-dsqrt(aux1*aux1+aux2*aux2)
     else if (av1(arc) < -nav1) then
        aux2   = miu/av1(arc)
        x(arc) = aux1+aux2+dsqrt(aux1*aux1+aux2*aux2)
     else
        x(arc) = aux1
     end if
     s(arc)    = du(arc)-x(arc)
     z(arc)    = miu/x(arc)
     w(arc)    = miu/s(arc)
     arc = upnext(arc)
  end do
  
  !     ------------------------------------------------------------------
  !     Computing initial values for variables without upper bound.
  !     ------------------------------------------------------------------  

  arc    = npfrst

  do i = 1, nbound
     aux1 = demand/2d0
     if (av1(arc) > 0.0d0) then
        aux2   = miu/av1(arc)
        x(arc) = aux1+aux2-dsqrt(aux1*aux1+aux2*aux2)
     else if (av1(arc) < 0.0d0) then
        aux2   = miu/av1(arc)
        x(arc) = aux1+aux2+dsqrt(aux1*aux1+aux2*aux2)      
     else
        x(arc) = aux1
     end if
     z(arc)  = miu/x(arc)
     arc = upnext(arc)
  end do

  !     ------------------------------------------------------------------
  !     Computing auxialliary vector av2.
  !     ------------------------------------------------------------------

  do arc = 1, na
     snod = strn(arc)
     enod = endn(arc)
     av2(snod) = av2(snod)+x(arc)
     av2(enod) = av2(enod)-x(arc)   
  end do

  av2(root) = 0.0d0

  !    ------------------------------------------------------------------
  !    Computing initial values to miu and tolpcg.
  !    ------------------------------------------------------------------

  tolpcg    = 0d0

  do node = 1, nn
     tolpcg = tolpcg+dabs(av2(node))
  end do
  
  gap = miu*dble(na)
  miu =  gap/dble(10*na)
  tolpcg = tolpcg/1d1

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_startvalues.
  !     ------------------------------------------------------------------

end subroutine pdnet_startvalues
  
!     ------------------------------------------------------------------
!     pdnet_tapia: 
!             This subroutine uses Tapia indicators to determine edge
!             classes = (low,active,cap).
!
!             p: min c^tx            d: max b^ty - u^tw
!                      Ax     = b           A^ty -    w + z = c
!                       x + s = u                     w,  z >=0
!                       x,  s >=0
!     ------------------------------------------------------------------

subroutine pdnet_tapia (actlst, edgcls, nact,   ncap,   nedge,  ssol,&
     tolhi,  tollo,  wsol,   xsol,   zsol )

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             nact  : number of active edges.
  !             ncap  : number of edges at capacity.
  !             nedge : number of edges.
  !     ------------------------------------------------------------------
  
  integer :: nact  , ncap  , nedge

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !             actlst(1:nedge): active edges list.
  !             edgcls(1:nedge): edges list.
  !     ------------------------------------------------------------------

  integer, dimension(1:nedge) :: actlst, edgcls

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             ssol(1:nedge): solution for s.
  !             wsol(1:nedge): solution for w.
  !             xsol(1:nedge): solution for x.
  !             zsol(1:nedge): solution for z.
  !     ------------------------------------------------------------------
  
  double precision, dimension(1:nedge) ::   ssol, wsol, xsol, zsol

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !     ------------------------------------------------------------------
  
  integer :: i

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             tolhi : upper bound of the tolerance for maxflow.
  !             tollo : lower bound of the tolerance for maxflow.
  !     ------------------------------------------------------------------
  
  double precision :: tolhi , tollo

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             swrat : ratio between s and w.
  !             xzrat : ratio between x and z.
  !     ------------------------------------------------------------------

  double precision  :: swrat , xzrat

  !     ------------------------------------------------------------------
  !     Start of executable portion of subroutine.
  !     ------------------------------------------------------------------

  nact=0
  ncap=0

  do i = 1, nedge
     xzrat=xsol(i)/zsol(i)
     swrat=ssol(i)/wsol(i)
     if ((xzrat > tolhi).and.(swrat < tollo)) then
        !           ------------------------------------------------------------
        !           Edge is set to upper bound.
        !           ------------------------------------------------------------
        edgcls(i)=2
        ncap=ncap+1 
     else
        if ((xzrat < tollo).and.(swrat > tolhi)) then
           edgcls(i)=0
        else
           !              ---------------------------------------------------------
           !              Edge is set active.
           !              ---------------------------------------------------------
           edgcls(i)=1
           nact=nact+1
           actlst(nact)=i
        end if
     end if
  end do

!     ------------------------------------------------------------------
!     End of subroutine pdnet_tapia.
!     ------------------------------------------------------------------

end subroutine pdnet_tapia
!     ------------------------------------------------------------------
!     pdnet_transform: This subroutine removes the lower bounds from the
!                      formulation.
!     ------------------------------------------------------------------


subroutine pdnet_transform ( b,      bound,  c,      endn,   fadd, &
     l,      na,     nbound, nn,     npfrst, strn,   u, upfrst, upnext )

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables
  !
  !             bound : number of bounded arcs.
  !             na    : number of arcs.
  !             nbound: number of unbounded arcs.
  !             nn    : number of nodes
  !             npfrst: next bounded arc to be considered.
  !             upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------

  integer :: bound , na    , nbound, nn    , npfrst, upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             b     (1:nn): array of capacities at each node.
  !             c     (1:na): array of costs at each arc.
  !             endn  (1:na): pointer of network data structure.
  !             l     (1:na): lower bound of flow at each arc.
  !             strn  (1:na): adjacency list of network data structure.
  !             u     (1:na): upper bound of flow at each arc.
  !             upnext(1:na): pointer to next arc in adjcncy structure.
  !     ------------------------------------------------------------------

  integer, dimension(1:nn) :: b
  integer, dimension(1:na) :: c, endn, l, strn, u, upnext

  !     ------------------------------------------------------------------
  !     Integer variables
  !
  !             arc   : current arc.
  !             enod  : current element on endn.
  !             i     : loop counter.
  !             snod  : current element on strn.
  !     ------------------------------------------------------------------

  integer ::  arc   , enod  , i     , snod

  !     ------------------------------------------------------------------
  !     Double precision variables
  !
  !             fadd: temp value of objective function.
  !     ------------------------------------------------------------------

  double precision fadd

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  fadd = 0.0d0
  arc = upfrst
  do i = 1, bound
     snod    = strn(arc)
     enod    = endn(arc)
     fadd    = fadd+dble(c(arc))*dble(l(arc))
     b(snod) = b(snod)+l(arc)
     b(enod) = b(enod)-l(arc)
     u(arc)  = u(arc)-l(arc)
     arc     = upnext(arc)
  end do

  arc = npfrst

  do i = 1, nbound
     snod    = strn(arc)
     enod    = endn(arc)
     fadd    = fadd+dble(c(arc))*dble(l(arc))
     b(snod) = b(snod)+l(arc)
     b(enod) = b(enod)-l(arc)
     arc     = upnext(arc) 
  end do
  
  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_transform.
  !     ------------------------------------------------------------------

end subroutine pdnet_transform
!     ------------------------------------------------------------------
!     pdnet_update1: 
!             This subroutine updates the vector x by using the formula
!             x = x + alfa*y, where y is a scalar. The vectors x and y
!             have both dimension dim. This subroutine is based on the
!             dblas daxpy routine.
!     ------------------------------------------------------------------

subroutine pdnet_update1 (alfa,   dim,    root,   x,      y)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             dim   : dimension of vectors x and y
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------
  
  integer :: dim   , root

  
  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             x(1:dim): vector x.
  !             y(1:dim): vector y, which will accomodate the result.
  !     ------------------------------------------------------------------

  double precision, dimension(1:dim) :: x, y

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             m     : size of unrolled loop.
  !             mp1   : m+1.
  !     ------------------------------------------------------------------

  integer :: i     , m     , mp1

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             alfa  : factor of y.
  !     ------------------------------------------------------------------

  double precision :: alfa

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  m = mod(dim,4)

  if (m /= 0) then
     do i = 1, m
         x(i) = x(i) + alfa*y(i)
      end do
      if (dim < 4) return
   end if

   mp1 = m + 1
   
   do i = mp1,dim,4
      x(i) = x(i) + alfa*y(i)
      x(i + 1) = x(i + 1) + alfa*y(i + 1)
      x(i + 2) = x(i + 2) + alfa*y(i + 2)
      x(i + 3) = x(i + 3) + alfa*y(i + 3) 
   end do

   x(root) = 0d0

   !     ------------------------------------------------------------------
   !     End of subroutine pdnet_update1.
   !     ------------------------------------------------------------------

 end subroutine pdnet_update1
!     ------------------------------------------------------------------
!     pdnet_update2: 
!             This subroutine updates the vector x by using the formula
!             x = y + alfa*x, wherw alfa is a scalar. The vectors x and
!             y have both dimension dim.
!     ------------------------------------------------------------------


subroutine pdnet_update2 (alfa,   dim,    root,   x,      y)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             dim   : dimension of vectors x and y
  !             root  : sets beginning of structure.
  !     ------------------------------------------------------------------
  
  integer :: dim   , root

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             x(1:dim): vector x, which will accomodate the result.
  !             y(1:dim): vector y.
  !     ------------------------------------------------------------------

  double precision, dimension(1:dim) :: x, y

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             m     : size of unrolled loop.
  !             mp1   : m+1.
  !     ------------------------------------------------------------------
 
  integer :: i     , m     , mp1

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             alfa  : factor of x.
  !     ------------------------------------------------------------------
  
  double precision :: alfa

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  m = mod(dim,4)

  if (m /= 0) then
     do  i = 1,m
        x(i) = y(i) + alfa*x(i)
     end do
     if (dim < 4) return
  end if

  mp1 = m + 1

  do i = mp1,dim,4
     x(i) = y(i) + alfa*x(i)
     x(i + 1) = y(i + 1) + alfa*x(i + 1)
     x(i + 2) = y(i + 2) + alfa*x(i + 2)
     x(i + 3) = y(i + 3) + alfa*x(i + 3)
  end do

  x(root) = 0d0

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_update2.
  !     ------------------------------------------------------------------

end subroutine pdnet_update2
!     ------------------------------------------------------------------
!     pdnet_updatesol: 
!             This subroutine updates the primal and dual solutions in
!             the network primal-dual interior-point algorithm.
!             Furthermore, this subroutine updates the parameter miu.
!     ------------------------------------------------------------------

subroutine pdnet_updatesol (av1,    av2,    av3,       b,      bound,&
     ds,     endn, factor, gap,&
     huge,   miu,    na,     nbound, nn,     npfrst, ps,     root,   s,&
     sd, strn,   teta,   tolpcg, u,      upfrst, upnext, w,      wd,&
     x,      xd,     y,      yd, z,      zd)

  implicit none


  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             bound : used for computing aditional data structures.
  !             na    : number of arcs.
  !             nbound: number of bounds.    
  !             nn    : number of nodes.
  !             npfrst: next bounded arc to be considered.
  !             root  : sets beginning of structure.
  !             upfrst: next unbounded arc to be considered.
  !     ------------------------------------------------------------------
  
  integer :: bound , na    , nbound, nn    , npfrst, root  , upfrst

  !     ------------------------------------------------------------------
  !     Integer input/output arrays:
  !
  !             endn  (1:na): pointer of network data structure.
  !             strn  (1:na): adjacency list of network data structure.
  !             upnext(1:na): next arc with positive upper bound.
  !     ------------------------------------------------------------------

  integer, dimension(1:na) :: endn, strn, upnext
  
  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             av1  (1:na): auxiliary vector 1.
  !             av2  (1:nn): auxiliary vector 2.
  !             av3  (1:na): auxiliary vector 3.
  !             b    (1:nn): array of capacities at each node.
  !             s    (1:na): dual variable s.
  !             sd   (1:na): direction for s.
  !             teta (1:na): diagonal matrix theta.
  !             u    (1:na): upper bound of flow at each arc.
  !             w    (1:na): primal variable w.
  !             wd   (1:na): direction for w.
  !             x    (1:na): primal variable x.
  !             xd   (1:na): direction for x.
  !             y    (1:nn): dual variable y.
  !             yd   (1:nn): direction for y.
  !             z    (1:na): primal variable z.
  !             zd   (1:na): direction for z.
  !     ------------------------------------------------------------------

  double precision, dimension(1:nn) :: av2, y, yd, b
  double precision, dimension(1:na) :: av1, av3,&
       s, sd, teta, u, w, wd, x, xd, z, zd

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             arc   : current arc.
  !             enod  : current element on endn.
  !             i     : loop counter.
  !             node  : another loop counter.
  !             snod  :  current element on strn.
  !     ------------------------------------------------------------------
  
  integer :: arc   , enod  , i     , node  , snod

  !     ------------------------------------------------------------------
  !     Double precision input/output variables:
  !
  !             factor: factor for lambda_max update.
  !             gap   : value of the gap.
  !             huge  : maxdouble precision number.
  !             miu   : value of interior-point mu parameter.
  !             tolpcg: initial value for tolcg.
  !     ------------------------------------------------------------------

  double precision :: factor, gap   , huge  , miu   , tolpcg

  !     ------------------------------------------------------------------
  !     Double precision variables:
  !
  !             aux1  : auxiliary variable.
  !             ds    : dual step.
  !             ps    : primal step.
  !     ------------------------------------------------------------------

  double precision :: aux1  , ds    , ps

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !
  !     Initializations.
  !     ------------------------------------------------------------------

  do node=1,nn
     av2(node) = b(node)
  end do

  gap = 0d0
  ps = huge
  ds = huge

  !     ------------------------------------------------------------------
  !     Compute Newton's direction xd, sd, zd and wd and maximum primal 
  !     and dual stepsizes ps and ds.
  !
  !     Arcs without upper bound.
  !     ------------------------------------------------------------------

  arc = npfrst

  do i=1,nbound
     
     !        ---------------------------------------------------------------
     !        Compute working parameters.
     !        ---------------------------------------------------------------
     
     snod = strn(arc)
     enod = endn(arc)
     aux1 = -yd(snod)+yd(enod)

     !        ---------------------------------------------------------------
     !        Compute the Newton's direction.
     !        ---------------------------------------------------------------

     xd(arc)  = aux1/teta(arc)+xd(arc)
     zd(arc)  = zd(arc)-av1(arc)*xd(arc)-z(arc)

     !        ---------------------------------------------------------------
     !        Compute the primal and dual steps.
     !        ---------------------------------------------------------------

     if (xd(arc) < 0.) then
        aux1 = -x(arc)/xd(arc)
        ps = dmin1(ps,aux1)
     endif
     if (zd(arc) < 0.) then
        aux1 = -z(arc)/zd(arc)
        ds = dmin1(ds,aux1)
     endif

     !        ---------------------------------------------------------------
     !        Compute the next arc.
     !        ---------------------------------------------------------------
     
     arc = upnext(arc)
  end do

  !     ------------------------------------------------------------------
  !     Arcs with upper bound.
  !     ------------------------------------------------------------------

  arc = upfrst
  
  do i=1,bound
     
     !        ---------------------------------------------------------------
     !        Compute working parameters.
     !        ---------------------------------------------------------------
     
     snod = strn(arc)
     enod = endn(arc)
     aux1 = -yd(snod)+yd(enod)

     !        ---------------------------------------------------------------
     !        Compute Newton's direction.
     !        ---------------------------------------------------------------
     
     xd(arc)  = aux1/teta(arc)+xd(arc)
     sd(arc)  = u(arc)-x(arc)-s(arc)-xd(arc)
     zd(arc)  = zd(arc)-av1(arc)*xd(arc)-z(arc)
     wd(arc)  = wd(arc)-av3(arc)*sd(arc)-w(arc)

     !        ---------------------------------------------------------------
     !        Compute primal and dual steps.
     !        ---------------------------------------------------------------
     
     if (xd(arc) < 0.) then
        aux1 = -x(arc)/xd(arc)
        ps = dmin1(ps,aux1)
     endif
     if (sd(arc) < 0.) then
        aux1 = -s(arc)/sd(arc)
        ps = dmin1(ps,aux1)
     endif
     if (zd(arc) < 0.) then
        aux1 = -z(arc)/zd(arc)
        ds = dmin1(ds,aux1)
     endif
     if (wd(arc) < 0.) then
        aux1 = -w(arc)/wd(arc)
        ds = dmin1(ds,aux1)
     endif
     
     !        ---------------------------------------------------------------
     !        Compute the next arc.
     !        ---------------------------------------------------------------
     
     arc = upnext(arc)

  end do
  
  !     ------------------------------------------------------------------
  !     Multiply primal and dual steps by factor.
  !     ------------------------------------------------------------------
  
  ps = ps*factor
  ds = ds*factor

  do node=1,nn
     y(node) = y(node)+ds*yd(node)
  end do

  !     ------------------------------------------------------------------
  !     Set y(root) equal to zero.
  !     ------------------------------------------------------------------
  
  y(root) = 0d0
  
  !     ------------------------------------------------------------------
  !     Update vectors x, w, z, and s.  Compute the duality gap.
  !
  !     Arcs with upper bound
  !     ------------------------------------------------------------------
  
  arc = upfrst

  do i=1,bound
     
     !        ---------------------------------------------------------------
     !        Compute working parameters.
     !        ---------------------------------------------------------------
     
     snod = strn(arc)
     enod = endn(arc)
     
     !        ---------------------------------------------------------------
     !        Update x, w, z, and s.
     !        ---------------------------------------------------------------
     
     x(arc) = x(arc)+ps*xd(arc)
     s(arc) = s(arc)+ps*sd(arc)
     z(arc) = z(arc)+ds*zd(arc)
     w(arc) = w(arc)+ds*wd(arc)
     
     !        ---------------------------------------------------------------
     !        Compute the duality gap.
     !        ---------------------------------------------------------------
     
     gap  = gap+x(arc)*z(arc)+s(arc)*w(arc)
     
     !        ---------------------------------------------------------------
     !        Compute auxilliary vector av2.
     !        ---------------------------------------------------------------
     
     av2(snod) = av2(snod)+x(arc)
     av2(enod) = av2(enod)-x(arc)
     arc = upnext(arc)

  end do

  !     ------------------------------------------------------------------
  !     Arcs without upper bound.
  !     ------------------------------------------------------------------
  
  arc = npfrst

  do i=1,nbound
     
     !        ---------------------------------------------------------------
     !        Compute working parameters.
     !        ---------------------------------------------------------------
     
     snod = strn(arc)
     enod = endn(arc)
     
     !        ---------------------------------------------------------------
     !        Update x, w, z, and s.
     !        ---------------------------------------------------------------
     
     x(arc) = x(arc)+ps*xd(arc)
     z(arc) = z(arc)+ds*zd(arc)
     
     !        ---------------------------------------------------------------
     !        Compute the duality gap.
     !        ---------------------------------------------------------------
     
     gap  = gap+x(arc)*z(arc)
     
     !        ---------------------------------------------------------------
     !        Compute auxilliary vector av2.
     !        ---------------------------------------------------------------
     
     av2(snod) = av2(snod)+x(arc)
     av2(enod) = av2(enod)-x(arc)
     arc = upnext(arc)
  end do
  
  !     ------------------------------------------------------------------
  !     Set av2(root) equal to zero.
  !     ------------------------------------------------------------------
  
  av2(root) = 0d0
  
  !     ------------------------------------------------------------------
  !     Compute the parameter miu and tolerance for conjugate gradients.
  !     ------------------------------------------------------------------
  
  tolpcg    = 0d0

  do node=1,nn
     tolpcg = tolpcg+dabs(av2(node))
  end do

  miu    = gap/dble(5*na)
  tolpcg = tolpcg/1d1

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_updatesol.
  !     ------------------------------------------------------------------

end subroutine pdnet_updatesol
!     ------------------------------------------------------------------
!     pdnet_xequaly: 
!             This subroutine performs the operation x = y.  The vectors 
!             x and y have both dimension dim.  This subroutine is based 
!             on dcopy dblas routine.
!     ------------------------------------------------------------------

subroutine pdnet_xequaly(dim,    root,   x,      y)

  implicit none

  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !             dim   : dimension of vectors.
  !             root  : beginning of data structure.
  !     ------------------------------------------------------------------
  

  integer :: dim   , root

  
  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             x(1:dim): vector x.
  !             y(1:dim): vector y.
  !     ------------------------------------------------------------------

  double precision, dimension(1:dim) :: x, y

  !     ------------------------------------------------------------------
  !     Integer variables:
  !
  !             i     : loop counter.
  !             m     : size of unrolled loop.
  !             mp1   : m + 1.
  !     ------------------------------------------------------------------
  
  integer :: i     , m     , mp1

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  m = mod(dim,7)

  if (m /= 0) then
     do i = 1, m
        x(i) = y(i)
     end do
     if (dim < 7) return
  end if

  mp1 = m + 1
  do i = mp1,dim,7
     x(i) = y(i)
     x(i + 1) = y(i + 1)
     x(i + 2) = y(i + 2)
     x(i + 3) = y(i + 3)
     x(i + 4) = y(i + 4)
     x(i + 5) = y(i + 5)
     x(i + 6) = y(i + 6)
  end do

   x(root) = 0d0

   !     ------------------------------------------------------------------
   !     End of subroutine pdnet_xequaly.
   !     ------------------------------------------------------------------

 end subroutine pdnet_xequaly

!     ------------------------------------------------------------------
!     pdnet_xminusy: 
!             This subroutine performs the operation z=x-y. The vectors
!             x, z and y have dimension dim. This subroutine is based on
!             the dblas daxpy routine.
!     ------------------------------------------------------------------

subroutine pdnet_xminusy(dim,    root,   x,      y,      z)

  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output variables:
  !
  !              dim   : dimension of vectors.
  !              root  : sets beginning of structure.
  !     ------------------------------------------------------------------
  
  integer ::  dim   , root

!     ------------------------------------------------------------------
!     Double precision arrays:
!
!              x(1:dim): first vector in subtraction.
!              y(1:dim): second vector in subtraction.
!              z(1:dim): x-y.
!     ------------------------------------------------------------------

  double precision, dimension(1:dim) :: x, y, z

  !     ------------------------------------------------------------------
  !     Integer  variables:
  !
  !             i     : loop counter.
  !             m     : size of unrolled loop.
  !             mp1   : m+1.
  !     ------------------------------------------------------------------
  
  integer ::  i     , m     , mp1

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  m = mod(dim,4)

  if (m /= 0) then
     do i = 1, m
        z(i) = x(i) - y(i)
     end do
     if (dim < 4) return
  end if
  mp1 = m + 1
  do i = mp1, dim, 4
     z(i) = x(i) - y(i)
     z(i + 1) = x(i + 1) - y(i + 1)
     z(i + 2) = x(i + 2) - y(i + 2)
     z(i + 3) = x(i + 3) - y(i + 3)
  end do
  z(root) = 0.0d0

  !     ------------------------------------------------------------------
  !     End of subroutine pdnet_xminusy.
  !     ------------------------------------------------------------------

end subroutine pdnet_xminusy
