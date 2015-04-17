module ftx_faltld_binding
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE 

! type, bind(c) :: vecl
!   integer(c_long), value :: n
!   integer(c_long) :: v
! end type vecl
! 
! type, bind(c) :: fxt_faltld
!   integer(c_long), value :: p
!   integer(c_long), value :: n
!   real(c_double), value :: prec
! end type fxt_flptld

interface
  subroutine fxt_faltld_preproc(p, n, mv, prec, fname) bind (c)
    use, intrinsic :: iso_c_binding
    integer(c_long), value :: p
    integer(c_long), value :: n
    real(c_double), value :: prec
    character(c_char) :: fname
    type(c_ptr), value :: mv 
  end subroutine fxt_faltld_preproc

  ! fxt faltld*
 
  ! deallocate fast spherical harmonic transform
  ! void fxt_faltld_del(fxt_faltld *falt);
  subroutine fxt_faltld_del(falt) bind (c)
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: falt
  end subroutine fxt_faltld_del

  ! size of working array
  ! long fxt_faltld_wsize(fxt_faltld *falt, long m); 
  integer(c_long) function fxt_faltld_wsize(falt, m) bind (c)
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: falt
    integer(c_long), value :: m
  end function fxt_faltld_wsize

  !  /*** maximum size of working array ***/
  ! long fxt_faltld_wsizemax(fxt_faltld *falt);
  integer(c_long) function fxt_faltld_wsizemax(falt) bind (c)
     use, intrinsic :: iso_c_binding
     type(c_ptr), value :: falt
  end function fxt_faltld_wsizemax

  ! /*** evaluate fast spherical harmonic transform ***/
  ! void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
  !                     fxt_vecld *u, fxt_vecld *w);
  subroutine fxt_faltld_evl(v, falt, m, u, w) bind(c)
     use, intrinsic :: iso_c_binding
     type(c_ptr), value :: v, falt, u, w
     integer(c_long), value :: m
  end subroutine fxt_faltld_evl

  !  /*** expand fast spherical harmonic transform ***/
  ! void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
  !                     fxt_vecld *v, fxt_vecld *w);
  subroutine fxt_faltld_exp(u, falt, m, v, w) bind(c)
     use, intrinsic :: iso_c_binding
     type(c_ptr), value :: u, falt, v, w
     integer(c_long), value :: m
  end subroutine fxt_faltld_exp


end interface

end module fxt_faltld_binding
