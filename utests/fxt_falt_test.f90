program test_fxtd_falt
  use iso_c_binding
  use fxt_fif
  use fxt_fwrap

  implicit none

  write(*,*) 'This utest is broken and needs to be fixed or deleted!'
  write(*,*) 'FAILED'
!HK!  character(c_char) :: fname
!HK!  integer(c_long) :: wsize
!HK!  integer(c_long) :: m
!HK!  type(c_ptr) :: falt
!HK!  type(c_ptr) :: w
!HK!
!HK!  real(kind=c_double) :: v_orig(10)
!HK!  real(kind=c_double), target :: u(10)
!HK!  real(kind=c_double), target :: v(10)
!HK!
!HK!  ! load fast spherical harmonic transform
!HK!  falt = fxt_faltld_init(10_c_long, 9_c_long, epsilon(1.0_c_double)) 
!HK!  write(*,*) 'after falt init'
!HK!  flush(6)
!HK! 
!HK!  ! size of working array
!HK!  wsize = fxt_faltld_wsize(falt)
!HK!  write(*,*) 'wsize: ', wsize
!HK!  flush(6)
!HK!
!HK!  ! create a new vector (working vector)
!HK!  w = fxt_vecld_new(wsize)
!HK!  write(*,*) 'allocated working vector'
!HK!  flush(6)
!HK!
!HK!  call random_number(v_orig)
!HK!  write(*,*) 'orig :', v_orig
!HK!  v = v_orig
!HK!
!HK!  ! there ....
!HK!  ! transform from physical to wave space
!HK!  call fxtf_faltld_exp(c_loc(u), 10_c_int, falt, m, c_loc(v), 10_c_int, w)
!HK!
!HK!  ! ...and back again
!HK!  ! transform from wave to physical space
!HK!  call fxtf_faltld_evl(c_loc(v), 10_c_int, m, falt, c_loc(u), 10_c_int, w)  
!HK!  write(*,*) 'trafo:', v
!HK!  write(*,*) 'Should be the same as orig.'
!HK!
!HK!  if (all(v - v_orig < 2.*epsilon(1.0_c_double))) then
!HK!    write(*,*) 'PASSED'
!HK!  else
!HK!    write(*,*) 'Data does not match after conversion:'
!HK!    write(*,*) 'FAILED'
!HK!  end if

  flush(6)

end program test_fxtd_falt
