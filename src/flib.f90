SUBROUTINE cross_product(z, x, y)
  !
  use constants,   only : dp
  !
  implicit none
  !
  real(dp), dimension(1:3) :: z, x, y
  !
  z(1)=x(2)*y(3)-x(3)*y(2)
  z(2)=x(3)*y(1)-x(1)*y(3)
  z(3)=x(1)*y(2)-x(2)*y(1)
  !
END SUBROUTINE

SUBROUTINE real2reciprocal(bv, av)
  !
  use constants,   only : dp
  !
  implicit none
  !
  real(dp),dimension(1:3,1:3) :: av, bv
  !
  real(dp) omega
  !
  omega=av(1,1)*av(2,2)*av(3,3)+av(1,2)*av(2,3)*av(3,1)+av(1,3)*av(2,1)*av(3,2)- &
       (av(1,3)*av(2,2)*av(3,1)+av(1,2)*av(2,1)*av(3,3)+av(1,1)*av(2,3)*av(3,2))
  !
  CALL cross_product(bv(:, 1), av(:, 2), av(:, 3))
  CALL cross_product(bv(:, 2), av(:, 3), av(:, 1))
  CALL cross_product(bv(:, 3), av(:, 1), av(:, 2))
  !
  bv(:, :)=bv(:, :)/omega
  !
END SUBROUTINE

FUNCTION distance(k1, k2, bv)
  !
  use constants,   only : dp
  !
  implicit none
  !
  real(dp) distance
  real(dp), dimension(1:3) :: k1, k2
  real(dp), dimension(1:3, 1:3) :: bv
  !
  integer ii
  !
  real(dp), dimension(1:3) :: dk
  !
  do ii=1, 3
    dk(ii)=sum(k1(:)*bv(ii,:))-sum(k2(:)*bv(ii,:))
  enddo
  !
  distance=sqrt(dk(1)**2+dk(2)**2+dk(3)**2)
  !
END FUNCTION
