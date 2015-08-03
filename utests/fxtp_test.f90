program fxtp_test
  use iso_c_binding
  use env_module, only: rk
  use fxt_fif
  use fxt_fwrap
  use ply_fxt_module, only:  ply_init_fxt

  implicit none
 
  type(fxtf_flptld_type) :: flpt
 
  real(kind=rk) :: v_orig(10)
  real(kind=rk), target :: u(10)
  real(kind=rk), target :: v(10)

  integer :: nPoints        ! number of points
  integer :: maxDegree      ! maximum degree
  real(kind=rk) :: prec

  nPoints = 10
  maxDegree = 9
  prec = 8*epsilon(prec)

  call ply_init_fxt(    flpt    = flpt,      &
    &                   degree  = maxDegree, &
    &                   nPoints = nPoints,   &
    &                   prec    = prec       )
  
  call random_number(v_orig)
  write(*,*) 'orig :', v_orig
  v = v_orig

  ! there ....
  ! transform from physical to wave space
  call fxtf_flptld_exp( c_loc(u), 10_c_int, flpt%handle, &
    &                   c_loc(v), 10_c_int, flpt%work    )

  ! ...and back again
  ! transform from wave to physical space
  call fxtf_flptld_evl( c_loc(v), 10_c_int, flpt%handle, &
    &                   c_loc(u), 10_c_int, flpt%work    )
  write(*,*) 'trafo:', v
  write(*,*) 'Should be the same as orig.'

  if (all(v - v_orig < 2.*epsilon(1.0_c_double))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Data does not match after conversion:'
    write(*,*) 'FAILED'
  end if

end program fxtp_test
