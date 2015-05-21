program test_fxtd
  use iso_c_binding
  use fxt_fif
  use fxt_fwrap

  implicit none
  integer, parameter :: rk = selected_real_kind(15)

  character(c_char) :: fname
  integer(c_long) :: wsize
  type(c_ptr) :: flpt
  type(c_ptr) :: w

  real(kind=c_double) :: v_orig(10)
  real(kind=c_double), target :: u(10)
  real(kind=c_double), target :: v(10)

  integer(8) :: p           ! number of points
  integer(8) :: n           ! maximum degree
  real(kind=rk) :: prec
  
!  load fast spherical harmonic transform
!  flpt = fxt_flptld_init(10_c_long, 9_c_long, epsilon(1.0_c_double)) 
!  write(*,*) 'after flpt init'
!  flush(6)
! 
!  ! size of working array
!  wsize = fxt_flptld_wsize(flpt)
!  write(*,*) 'wsize: ', wsize
!  flush(6)
!
!  ! create a new vector (working vector)
!  w = fxt_vecld_new(wsize)
!  write(*,*) 'allocated working vector'
!  flush(6)

  p = 10
  n = 9
  prec = 1.0

  call fxtf_flptld_init(p, n, prec, flpt, w)
  
  call random_number(v_orig)
  write(*,*) 'orig :', v_orig
  v = v_orig

  ! there ....
  ! transform from physical to wave space
  call fxtf_flptld_exp(c_loc(u), 10_c_int, flpt, c_loc(v), 10_c_int, w)

  ! ...and back again
  ! transform from wave to physical space
  call fxtf_flptld_evl(c_loc(v), 10_c_int, flpt, c_loc(u), 10_c_int, w)  
  write(*,*) 'trafo:', v
  write(*,*) 'Should be the same as orig.'

  if (all(v - v_orig < 2.*epsilon(1.0_c_double))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Data does not match after conversion:'
    write(*,*) 'FAILED'
  end if

  flush(6)

end program test_fxtd
