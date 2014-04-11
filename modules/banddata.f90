MODULE banddata
  !
  use constants
  !
  implicit none
  !
  real(dp), allocatable :: xk(:)
  real(dp), allocatable :: eig(:, :)
  real(dp), allocatable :: proj(:, :, :)
  !
  integer nbnd
  !
  integer nkpt
  !
CONTAINS

SUBROUTINE init_band
  !
  use input,     only : nkvec, kvec, lproj, bvec
  !
  implicit none
  !
  integer ik
  real(dp) distance
  !
  nkpt=nkvec
  !
  allocate(xk(1:nkpt))
  allocate(eig(1:nbnd, 1:nkpt))
  !
  xk(1) = 0.d0
  !
  do ik=2, nkpt
    xk(ik)=xk(ik-1)+distance(kvec(:, ik-1), kvec(:, ik), bvec)
  enddo
  !
  eig(:, :)=0.d0
  !
  if (lproj) then
    allocate(proj(1:nbnd, 1:nbnd, 1:nkpt))
    proj(:, :, :)=0.d0
  endif
  !
END SUBROUTINE

SUBROUTINE output_band
  !
  use constants, only : stdout
  use input,     only : lproj
  !
  implicit none
  !
  integer ik, ib, ip
  !
  do ib=1, nbnd
    do ik=1, nkpt
      if(lproj) then
        write(stdout, '(2F12.4,50F12.4)') xk(ik), eig(ib, ik), proj(:, ib, ik)
      else
        write(stdout, '(2F12.4)') xk(ik), eig(ib, ik)
      endif
    enddo
    write(stdout, " ")
  enddo
  !
END SUBROUTINE

SUBROUTINE finalize_band
  !
  implicit none
  !
  if(allocated(xk)) deallocate(xk)
  if(allocated(eig)) deallocate(eig)
  if(allocated(proj)) deallocate(proj)
  !
END SUBROUTINE

END MODULE

