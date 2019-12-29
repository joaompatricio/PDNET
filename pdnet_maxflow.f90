!     -----------------------------------------------------------------
!     pdnet_maxflow: 
!             This subroutine computes a maximum flow across a network
!             comprised of the active edges in the optimal solution.
!
!             1. Correct the flow balances in or out of each vertex, as
!                implied by edges in upper and lower bound.
!
!             2. Build a list of virtual edges from a source vertex into
!                vertices with positive flow balance and from vertices
!                with negative flow balance into a sink vertex.
!
!             3. Solve maximum flow problem with Dinic's algorithm.
!
!             4. Build a primal maxflow solution.
!             
!
!     Authors: Luis Portugal (Univ. Coimbra - Portugal)
!              Mauricio G. C. Resende (AT&T Research Labs - USA)
!              Geraldo Veiga (AT&T Research Labs - USA)
!              Joaquim Judice (Univ. Coimbra - Portugal)
!              Joao Patricio (EST Tomar - Portugal)
!     ------------------------------------------------------------------

subroutine pdnet_maxflow (actlst, cap, endv, mflow, mfreq, mxfe, mxfv, &
     nact, ncap, nedge, nvrtx, strv, supply, balance, xopt)

  implicit none
  
  !----------------------------------------------------------------------
  
  integer :: mflow, mfreq, mxfe, mxfv, nact, ncap, nedge, nvrtx, edg, i,&
       mfedge, mfsnk, mfsrc, mfvrtx, nacut, nncut, numaug, numstg, rc, &
       tmp, v1, v2

  integer, dimension(0:nedge - 1) :: actlst, endv, strv, balance, xopt


  integer, dimension (:), allocatable :: dnfadj
  integer, dimension (:), allocatable :: dnfrom
  integer, dimension (:), allocatable :: dnmap
  integer, dimension (:), allocatable :: dnto
  integer, dimension (:), allocatable :: dnbadj
  integer, dimension (:), allocatable :: dngcap
  integer, dimension (:), allocatable :: dncap
  integer, dimension (:), allocatable :: dnflow
  integer, dimension (:), allocatable :: dnbtof
  integer, dimension (:), allocatable :: dnlist
  integer, dimension (:), allocatable :: dnfapt
  integer, dimension (:), allocatable :: dndist
  integer, dimension (:), allocatable :: dnptrf
  integer, dimension (:), allocatable :: dnbapt
  integer, dimension (:), allocatable :: dnflab
  integer, dimension (:), allocatable :: dnptrb
  
  double precision, dimension(0:nedge - 1) :: cap, supply
  
  !     ----------------------------------------------
  !     Beginning of executable section of subroutine.
  !     ----------------------------------------------
  !
  !     ---------------------------
  !     Compute number of vertices.
  !     ---------------------------

  mfvrtx = nvrtx + 2
  mfsrc  = nvrtx + 1
  mfsnk  = nvrtx + 2
  
  !     -------------------------------------
  !     Initialize capacity of virtual edges.
  !     -------------------------------------

  do i = 0, nvrtx
     balance(i) = -int(supply(i))
  end do

  !     ------------------------------------------
  !     Adjust flow balance for edges at capacity.
  !     ------------------------------------------

  do i = nact, nact + ncap - 1
     edg = actlst(i) - 1
     v1  = strv(edg) - 1
     v2  = endv(edg) - 1
     balance(v1) = balance(v1) - int(cap(edg))
     balance(v2) = balance(v2) + int(cap(edg))
  end do

  !     -----------------------------------------------
  !     Compute number of vertices and edges in maxflow
  !     -----------------------------------------------

  mfedge = nact
  mfreq = 0
  do i = 0, nvrtx - 1
     if (balance(i) /= 0) then
        mfedge = mfedge + 1
        if (balance(i) > 0) mfreq = mfreq + balance(i)
     end if
  end do

  !     -------------------------------------------
  !     Dynamically allocate some auxiliary arrays.
  !     -------------------------------------------

  allocate (dnmap(0:mfedge + 1))
  allocate (dnto(0:mfedge + 1))
  allocate (dnfadj(0:mfedge + 1))
  allocate (dnbadj(0:mfedge + 1))
  allocate (dncap(0:mfedge + 1))
  allocate (dnbtof(0:mfedge + 1))
  allocate (dnflow(0:mfedge + 1))
  allocate (dnfrom(0:mfedge + 1))
  allocate (dngcap(0:mfedge + 1))
  allocate (dnlist(0:mfvrtx + 1))
  allocate (dnfapt(0:mfvrtx + 1))
  allocate (dndist(0:mfvrtx + 1))
  allocate (dnptrf(0:mfvrtx + 1))
  allocate (dnbapt(0:mfvrtx + 1))
  allocate (dnflab(0:mfvrtx + 1))
  allocate (dnptrb(0:mfvrtx + 1))

  !     ----------------------
  !     Generate active edges.
  !     ----------------------

  do i = 0, nact - 1
     edg = actlst(i) - 1
     tmp = int(cap(edg))
     dnfrom(i) = strv(edg)
     dnto(i) = endv(edg)
     dngcap(i) = min(tmp, mfreq)
  end do

  !     -----------------------
  !     Generate virtual edges.
  !     -----------------------
  i = 0
  edg = nact 
  do i = 0, nvrtx - 1
     if (balance(i) < 0) then
        dnfrom(edg) = i + 1
        dnto(edg) = mfsnk
        dngcap(edg) = -balance(i)
        edg = edg + 1
     elseif (balance(i) > 0) then
        dnfrom(edg) = nvrtx + 1
        dnto(edg) = i + 1
        dngcap(edg) = balance(i)
        edg = edg + 1
     end if
  end do

  !     ----------------------------------------------------------
  !     Create input data in fwd adjacencies as required by dnsub.
  !     ----------------------------------------------------------
  
  rc = 0
  call dnfwd(mfvrtx, mfedge, rc, dnlist, dnfapt, dnfadj, &
       dnfrom, dnto, dncap, dngcap, dnmap)
  if (rc /= 0) then
     print *, '-- PDMXFL >> Error in DNFWD.'
     stop
  end if
  
  !     ----------------------
  !     Solve maxflow problem.
  !     ----------------------
  
  mxfv = mfvrtx
  mxfe = mfedge
  call dnsub(mfvrtx, mfedge, mfsrc, mfsnk, mflow, numaug, numstg, &
       nncut, nacut, rc, dnlist, dndist, dnfapt, dnptrf, dnbapt, &
       dnflab, dnptrb, dnfadj,   dnbadj, dncap, &
       dnflow, dnbtof)
  if (rc /= 0) then
     print *, '-- PDMXFL >> Error in DNSUB.'
     stop
  end if

  !     -------------------------
  !     Build a maxflow solution.
  !     -------------------------

  do i = 0, nact - 1
     edg = actlst(i) - 1
     xopt(edg) = dnflow(dnmap(i) - 1)
  end do

  do i = nact, nact + ncap - 1
     edg = actlst(i) - 1
     xopt(edg) = int(cap(edg))
  end do

  do i = nact + ncap, nedge - 1
     edg = actlst(i) - 1
     xopt(edg) = 0
  end do


  deallocate (dnmap)
  deallocate (dnto)
  deallocate (dnfadj)
  deallocate (dnbadj)
  deallocate (dncap)
  deallocate (dnbtof)
  deallocate (dnflow)
  deallocate (dnfrom)
  deallocate (dngcap)
  deallocate (dnlist)
  deallocate (dnfapt)
  deallocate (dndist)
  deallocate (dnptrf)
  deallocate (dnbapt)
  deallocate (dnflab)
  deallocate (dnptrb)

  !     --------------------------------
  !     End of subroutine pdnet_maxflow.
  !     --------------------------------

end subroutine pdnet_maxflow
