! pdnet_read: this subroutine reads a DIMACS-format network data
!             file, creating the data structure returned via
!             subroutine arguments.

subroutine pdnet_read (nodes, arcs, maxn, maxa, strv, endv, supply, cost, lwbnd, cap)
  
  implicit none

  integer :: arcs, cst, dst, edges, itmp, lines, low, nodes, src, vrtx, i
  integer :: maxn, maxa
  integer, dimension (1:maxn) :: supply
  integer, dimension (1:maxa) :: cap, cost, endv, lwbnd, strv

  character (len=1) desc
  character (len=3) problm

  logical readp

  !     ----------------------------------------------
  !     Beginning of executable section of subroutine.
  !     ----------------------------------------------
  lines = 0
  edges = 0
  readp = .false.
  do 
     read (*, fmt = 1000, end = 500, advance = 'no') desc
     lines = lines + 1
     if (desc == 'p' .or. desc == 'P') then
        if (.not. readp) then
           read *, problm, nodes, arcs
           readp = .true.
           do i = 1, nodes
              supply(i) = 0
           end do
        else
           print *, 'Error:'
           print *, 'P section appears twice in '
           print *, 'DIMACS file.'
           print *, 'Aborting...'
           stop
        endif
     elseif (desc == 'a' .or. desc == 'A') then
        read *, src, dst, low, itmp, cst
        edges = edges + 1
        strv(edges) = src
        endv(edges) = dst
        lwbnd(edges) = low
        cost(edges) = cst
        cap(edges) = itmp
     elseif (desc == 'n' .or. desc == 'N') then
        read *, vrtx, itmp
        supply(vrtx) = -itmp
     elseif (desc == 'c' .or. desc == 'C') then
        read *, problm
     endif
  enddo
  
  !     ------------
  !     End of file.
  !     ------------
500 return
  
  !     --------
  !     Formats.
  !     --------
1000 format(a1)

  !     -------------------------
  !     End of subroutine pdnet_read.
  !     -------------------------
end subroutine pdnet_read
