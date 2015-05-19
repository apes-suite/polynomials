program test_fxtd_falt
  use iso_c_binding
  use fxt_fif
  use fxt_fwrap

  implicit none

  character(c_char) :: fname
  integer(c_long) :: wsize
  integer(c_long) :: m
  type(c_ptr) :: falt
  type(c_ptr) :: w

  real(kind=c_double) :: v_orig(10)
  real(kind=c_double), target :: u(10)
  real(kind=c_double), target :: v(10)

  ! load fast spherical harmonic transform
  falt = fxt_faltld_init(10_c_long, 9_c_long, epsilon(1.0_c_double)) 
  write(*,*) 'after falt init'
  flush(6)
 
  ! size of working array
  wsize = fxt_faltld_wsize(falt)
  write(*,*) 'wsize: ', wsize
  flush(6)

  ! create a new vector (working vector)
  w = fxt_vecld_new(wsize)
  write(*,*) 'allocated working vector'
  flush(6)

  call random_number(v_orig)
  write(*,*) 'orig :', v_orig
  v = v_orig

  ! there ....
  ! transform from physical to wave space
  call fxtf_faltld_exp(c_loc(u), 10_c_int, falt, m, c_loc(v), 10_c_int, w)

  ! ...and back again
  ! transform from wave to physical space
  call fxtf_faltld_evl(c_loc(v), 10_c_int, m, falt, c_loc(u), 10_c_int, w)  
  write(*,*) 'trafo:', v
  write(*,*) 'Should be the same as orig.'

  if (all(v - v_orig < 2.*epsilon(1.0_c_double))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Data does not match after conversion:'
    write(*,*) 'FAILED'
  end if

  flush(6)

end program test_fxtd_falt
