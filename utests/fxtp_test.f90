program fxtp_test
  use iso_c_binding
  use env_module, only: rk
  !use fxt_fif
  use fxt_fwrap
  use ply_fxt_module, only:  ply_init_fxt

  implicit none
 
  type(fxtf_flptld_type), target :: flpt
 
  real(kind=rk) :: nodal_orig(10)
  real(c_double) :: modal(10)
  real(c_double) :: nodal(10)
  integer(c_int):: modalLen, nodalLen 

  integer :: nPoints        ! number of points
  integer :: maxDegree      ! maximum degree
  real(kind=rk) :: prec

  nPoints = 10
  maxDegree = 9
  prec = 8*epsilon(prec)
  modalLen = nPoints 
  nodalLen = nPoints 

  call ply_init_fxt(    flpt    = flpt,      &
    &                   degree  = maxDegree, &
    &                   nPoints = nPoints,   &
    &                   prec    = prec       )
  
  call random_number(nodal_orig)
  write(*,*) 'orig :', nodal_orig
  nodal = nodal_orig

  write(*,*) "passing nodal values = ", nodal
  write(*,*) "modalLen, nodalLen = ", modalLen, nodalLen  
  ! ther e....
  ! transform from physical to wave space
  call fxtf_flptld_exp( modal, modalLen, flpt%handle, &
    &                   nodal, nodalLen, flpt%work    )

  call fxt_error_print()

  write(*,*) "modal val", modal

  ! ...and back again
  ! transform from wave to physical space
  call fxtf_flptld_evl( nodal, nodalLen, flpt%handle, &
    &                   modal, modalLen, flpt%work    )

  write(*,*) 'trafo:', nodal
  write(*,*) 'Should be the same as orig.'

  if (all(abs(nodal - nodal_orig) < 2.*epsilon(1.0_c_double))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Data does not match after conversion:'
    write(*,*) 'FAILED'
  end if

end program fxtp_test
