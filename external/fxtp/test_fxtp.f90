program test_fxtd
  use iso_c_binding
  use fxt_flptld_binding
  use fxt_faltld_binding
  implicit none

  character(c_char) :: fname
  integer(c_long) :: m
  type(c_ptr) :: u, v, w
  
  ! load fast spherical harmonic transform
  falt = fxt_faltld_load(fname) 
  
  ! size of working array
  wsize = fxt_faltld_wsizemax(falt)

  ! create a new vector (working vector)
  w = fxt_vecld_new(wsize)
   
  ! evaluate fast spherical harmonic transform
  call fxt_faltld_evl(v, falt, m, u, w)  

  ! expand fast spherical harmonic transform
  call fxt_faltld_exp(u, falt, m, v, w)

end program test_fxtd
