include 'lapack.f90'

SUBROUTINE interpolate_bands(ik)
  !
  use lapack95,  only : heev
  use input,     only : kvec, lproj
  use constants, only : dp, twopi, cmplx_0, cmplx_i
  use wanndata,  only : rvec, ham, weight, nrpt, norb
  use banddata,  only : eig, proj, nbnd
  !
  implicit none
  !
  integer ik
  real(dp) rdotk
  complex(dp) fact
  complex(dp), allocatable :: work(:, :)
  !
  integer ir, info
  !
  allocate(work(1:nbnd, 1:nbnd))
  work(:,:)=cmplx_0
  !
  do ir=1, nrpt
    rdotk=SUM(kvec(:,ik)*rvec(:, ir))
    fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
    work(:,:)=work(:,:)+ham(:,:,ir)*fact
  enddo ! ir
  !
  call heev(work, eig(:, ik), 'V', 'U', info)
  !
  ! work(ii, io): ii\th element eigenvector of eig(io)
  !
  if(lproj) then
    do ir=1, nbnd
      do info=1, nbnd
        proj(info, ir, ik)=real(work(info, ir)*conjg(work(info, ir)))
      enddo
    enddo
  endif
  !
END SUBROUTINE
