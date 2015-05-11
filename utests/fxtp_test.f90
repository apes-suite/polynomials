program test_fxtd
  use iso_c_binding
  use fxt_fif

  implicit none

  character(c_char) :: fname
  integer(c_long) :: wsize
  type(c_ptr) :: flpt
  type(c_ptr) :: u, v, w
  
  ! load fast spherical harmonic transform
  flpt = fxt_flptld_init(10_c_long, 9_c_long, epsilon(1.0_c_double)) 
  
  ! size of working array
  wsize = fxt_flptld_wsize(flpt)

  ! create a new vector (working vector)
  w = fxt_vecld_new(wsize)
   
  ! evaluate fast spherical harmonic transform
  call fxt_flptld_evl(v, flpt, u, w)  

  ! expand fast spherical harmonic transform
  call fxt_flptld_exp(u, flpt, v, w)

end program test_fxtd
