PROGRAM wannband
  !
  use para,     only : init_para, inode, distribute_k, finalize_para, first_k, last_k, para_merge
  use wanndata, only : read_ham, norb, finalize_wann
  use banddata
  use input,    only : read_input, seed, lproj, finalize_input
  !
  implicit none
  !
  integer ik
  !
  CALL init_para
  CALL read_input
  CALL read_ham(seed)
  !
  nbnd=norb
  !
  CALL init_band
  CALL distribute_k(nkpt)
  !
  do ik=first_k, last_k
    !
    CALL interpolate_bands(ik)
    !
  enddo
  !
  CALL para_merge(eig, nbnd, nkpt)
  !
  if (lproj) then
    CALL para_merge(proj, nbnd, nbnd, nkpt)
  endif
  !
  if (inode.eq.0) CALL output_band
  !
  CALL finalize_input
  !
  CALL finalize_wann
  !
  CALL finalize_band
  !
  CALL finalize_para
  !
END PROGRAM
