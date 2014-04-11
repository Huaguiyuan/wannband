MODULE input
  !
  use constants
  !
  implicit none
  !
  character(len=80) seed
  real(dp) bvec(1:3,1:3)
  logical lproj
  integer nkvec
  real(dp), allocatable :: kvec(:, :)
  !
CONTAINS
  !
 SUBROUTINE read_input
  !
  !************ INPUT FILE *************
  !** file name: wannband.inp
  !line 1: seed name
  !line 2: A1
  !line 3: A2
  !line 4: A3
  !line 5: mode
  !line 6: nkseg, nk_per_seg
  !line 7.. : high_symmetry_kpts
  !*************************************
  !
  use constants, only : fin
  use para
  !
  implicit none
  !
  real(dp) avec(1:3, 1:3)
  real(dp), allocatable :: kbnd_vec(:, :)
  integer mode, nksec, nk_per_sec, iksec, ik
  !
  if (inode.eq.0) then
    open(unit=fin, file="wannband.inp")
    !
    read(fin, *) seed
    !
    read(fin, *) avec(:, 1)
    read(fin, *) avec(:, 2)
    read(fin, *) avec(:, 3)
    read(fin, *) mode
    read(fin, *) nksec, nk_per_sec
    !
  endif
  !
  CALL para_sync(avec, 3, 3)
  !
  CALL real2reciprocal(bvec, avec)
  !
  CALL para_sync(mode)
  CALL para_sync(nksec)
  CALL para_sync(nk_per_sec)
  nkvec=(nksec-1)*nk_per_sec+1
  allocate(kvec(1:3, 1:nkvec))
  if (inode.eq.0) then
    allocate(kbnd_vec(1:3, 1:nksec))
    do iksec=1, nksec
      read(fin, *) kbnd_vec(:, iksec)
    enddo
    do iksec=1, nksec-1
      do ik=1, nk_per_sec
        kvec(:, (iksec-1)*nk_per_sec+ik)=kbnd_vec(:, iksec)+(kbnd_vec(:, iksec+1)-kbnd_vec(:, iksec))*(ik-1)/nk_per_sec
      enddo
    enddo
    kvec(:, nkvec)=kbnd_vec(:, nksec)
    deallocate(kbnd_vec)
  endif
  CALL para_sync(kvec, 3, nkvec)
  !
  if (mode.eq.0) then 
    lproj = .false.
  else
    lproj = .true.
  endif
  !
 END SUBROUTINE
 !
 SUBROUTINE finalize_input
   !
   implicit none
   !
   if(allocated(kvec)) deallocate(kvec)
   !
 END SUBROUTINE
  !
END MODULE
