!
!     ------------------------------------------------------------------
!     pdnet_default: Sets in info array default values for PDNET internal 
!                    parameters.  Some parameters are set as function of
!                    instance size.
!     ------------------------------------------------------------------ 
subroutine pdnet_default (dinfo, info, na, nn)
  
  implicit none
  
  !     ------------------------------------------------------------------
  !     Integer input/output:
  !
  !             na        : number of arcs.
  !             nn        : number of nodes.
  !             info(1:23): array for integer stats and params.
  !     ------------------------------------------------------------------
  
  integer :: na, nn, i
  integer, dimension (1:30) :: info

  !     ------------------------------------------------------------------
  !     Double precision input/output arrays:
  !
  !             dinfo(1:35): array for double stats and params.
  !     ------------------------------------------------------------------

  double precision, dimension (1:35) :: dinfo
  double precision :: aux1, aux2

  !     ------------------------------------------------------------------
  !     Start of executable section of subroutine.
  !     ------------------------------------------------------------------

  do i = 1, 35
     dinfo(i) = 0.0d0
  end do

  dinfo(26) = 1.0d30
  info ( 3) = 6
  dinfo( 3) = 9.0d-1
  dinfo(30) = 1.0d-20
  dinfo(18) = 2.0d-1
  dinfo(23) = 1.0d-3
  dinfo(24) = 1.0d-3
  dinfo(11) = 1.0d-6
  dinfo( 7) = 9.75d-1
  info ( 9) = na
  dinfo(21) = 1.0d-2
  dinfo(22) = 1.0d2
  dinfo( 8) = 9.5d-1
  dinfo(25) = 1.0d-8
  dinfo(19) = 1.0d-8
  dinfo(20) = 1.0d-8
  dinfo(17) = 1.0d-3
  dinfo( 6) = 1.0d-1
  dinfo(12) = 1.0d-4
  dinfo( 4) = 1.0d-8
  info (11) = 1
  aux1 = dlog(dble(nn))
  aux2 = alog(2.0)
  aux1 = aux1 / aux2
  info (23) = 2*int(aux1) + 2
  dinfo(31) = 4.0d0
  dinfo(32) = 5.0d-2
  dinfo(33) = dsqrt(dble(nn))/dinfo(31)
  dinfo(34) = dsqrt(dble(nn))/dinfo(32)

  !    ------------------------------------------------------------------
  !     End of subroutine pdnet_default.
  !     ------------------------------------------------------------------
end subroutine pdnet_default
