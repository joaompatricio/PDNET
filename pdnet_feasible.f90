!     -----------------------------------------------------------------
!     pdfeas: this subroutine is a slight modification of subroutine
!             pdmxfl, for computing a maximum flow across a network
!             comprised of the active edges without computing the
!             maximum flow array. This subroutine is invoked to test
!             feasibility of the active edge network.
!     -----------------------------------------------------------------

subroutine pdnet_feasible  (actlst, cap, endv, mflow, mfreq, mxfe, mxfv,&
     nact, ncap, nedge, nvrtx, strv, supply, balance)

  implicit none

  integer :: mflow, mfreq, mxfe, mxfv, nact, ncap, nedge, nvrtx, edg, i, mfedge,&
       mfsnk, mfsrc, mfvrtx, nacut, nncut, numaug, numstg, rc, tmp, v1, v2

  integer, dimension(0:nedge - 1) :: actlst, endv,  strv, balance

  integer, dimension (:), allocatable :: dnfadj, dnfrom, dnmap, dnto, dnbadj,&
       dngcap, dncap, dnflow, dnbtof, dnlist, dnfapt, dndist, dnptrf, dnbapt,&
       dnflab, dnptrb

  double precision, dimension(0:nedge-1) :: cap
  double precision, dimension(0:nvrtx-1) :: supply

  !     ----------------------------------------------
  !     Start of executable section of subroutine.
  !     ----------------------------------------------

  
  !     ----------------------------------------------
  !     Compute number of vertices.
  !     ----------------------------------------------

  mfvrtx = nvrtx + 2
  mfsrc  = nvrtx + 1
  mfsnk  = nvrtx + 2
  
  !     -------------------------------------
  !     Initialize capacity of virtual edges.
  !     -------------------------------------

  do i = 0, nvrtx-1
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

  mfreq = 0
  mfedge = nact


  do i = 0, nvrtx - 1
     if (balance(i) /= 0) then
        mfedge = mfedge + 1
        if (balance(i) > 0) then
           mfreq = mfreq + balance(i)
        end if
     endif
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
  allocate (dnlist(1:mfvrtx + 2))
  allocate (dnfapt(1:mfvrtx + 2))
  allocate (dndist(1:mfvrtx + 2))
  allocate (dnptrf(1:mfvrtx + 2))
  allocate (dnbapt(1:mfvrtx + 2))
  allocate (dnflab(1:mfvrtx + 2))
  allocate (dnptrb(1:mfvrtx + 2))
  
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
     else if (balance(i) > 0) then
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

  call dnfwd(mfvrtx, mfedge, rc, dnlist, dnfapt, dnfadj, dnfrom, dnto, &
       dncap, dngcap, dnmap)  

  if (rc /= 0) then
     print *,  '-- pdnet_feasible >> Error in DNFWD.'
     stop
  end if

  !     ----------------------
  !     Solve maxflow problem.
  !     ----------------------

  mxfv = mfvrtx
  mxfe = mfedge
  call dnsub(mfvrtx, mfedge, mfsrc, mfsnk, mflow, numaug, numstg,nncut, nacut, &
       rc, dnlist, dndist, dnfapt, dnptrf, dnbapt, dnflab, dnptrb, dnfadj, &
       dnbadj, dncap, dnflow, dnbtof)  

  if (rc /= 0) then
     print *, '-- pdnet_feasible >> Error in DNSUB.'
     stop
  end if

  !     ---------------------------------
  !     End of subroutine pdnet_feasible.
  !     ---------------------------------
end subroutine pdnet_feasible
