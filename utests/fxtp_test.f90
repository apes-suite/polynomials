program test_fxtd
  use iso_c_binding
  use fxt_fif
  use fxt_fwrap

  implicit none

  character(c_char) :: fname
  integer(c_long) :: wsize
  type(c_ptr) :: flpt
  type(c_ptr) :: w

  real(kind=c_double), target :: u(10)
  real(kind=c_double), target :: v(10)
  
  ! load fast spherical harmonic transform
  flpt = fxt_flptld_init(10_c_long, 9_c_long, epsilon(1.0_c_double)) 
  write(*,*) 'after flpt init'
  flush(6)
  
  ! size of working array
  wsize = fxt_flptld_wsize(flpt)
  write(*,*) 'wsize: ', wsize
  flush(6)

  ! create a new vector (working vector)
  w = fxt_vecld_new(wsize)
  write(*,*) 'allocated working vector'
  flush(6)
   
  ! evaluate fast spherical harmonic transform
  call fxtf_flptld_evl(c_loc(v), 10_c_int, flpt, c_loc(u), 10_c_int, w)  

  ! expand fast spherical harmonic transform
  call fxtf_flptld_exp(c_loc(u), 10_c_int, flpt, c_loc(v), 10_c_int, w)

  write(*,*) 'PASSED'
  flush(6)

end program test_fxtd
