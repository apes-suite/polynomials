module fxt_faltld_binding
  use, intrinsic :: iso_c_binding
  use fxt_types
  IMPLICIT NONE 

interface
  subroutine fxt_faltld_preproc(p, n, mv, prec, fname) bind (c)
    use, intrinsic :: iso_c_binding
    use fxt_types
    integer(c_long), value :: p
    integer(c_long), value :: n
    real(c_double), value :: prec
    character(c_char) :: fname
    type(fxt_vecl), value :: mv 
  end subroutine fxt_faltld_preproc

  ! fxt faltld*
 
  ! deallocate fast spherical harmonic transform
  ! void fxt_faltld_del(fxt_faltld *falt);
  subroutine fxt_faltld_del(falt) bind (c)
    use, intrinsic :: iso_c_binding
    use fxt_types
   type(fxt_faltld), value :: falt
  end subroutine fxt_faltld_del

  ! size of working array
  ! long fxt_faltld_wsize(fxt_faltld *falt, long m); 
  integer(c_long) function fxt_faltld_wsize(falt, m) bind (c)
    use, intrinsic :: iso_c_binding
    use fxt_types
    type(fxt_faltld) :: falt
    integer(c_long), value :: m
  end function fxt_faltld_wsize

  !  /*** maximum size of working array ***/
  ! long fxt_faltld_wsizemax(fxt_faltld *falt);
  integer(c_long) function fxt_faltld_wsizemax(falt) bind (c)
     use, intrinsic :: iso_c_binding
     use fxt_types
     type(fxt_faltld) :: falt
  end function fxt_faltld_wsizemax

  ! /*** evaluate fast spherical harmonic transform ***/
  ! void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
  !                     fxt_vecld *u, fxt_vecld *w);
  subroutine fxt_faltld_evl(v, falt, m, u, w) bind(c)
     use, intrinsic :: iso_c_binding
     use fxt_types
     ! type(c_ptr), value :: v, falt, u, w
     type(fxt_vecld) :: u, v, w
     type(fxt_faltld) :: falt
     integer(c_long), value :: m
  end subroutine fxt_faltld_evl

  !  /*** expand fast spherical harmonic transform ***/
  ! void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
  !                     fxt_vecld *v, fxt_vecld *w);
  subroutine fxt_faltld_exp(u, falt, m, v, w) bind(c)
     use, intrinsic :: iso_c_binding
     use fxt_types
     type(fxt_vecld) :: u, v, w
     type(fxt_faltld) :: falt
     integer(c_long), value :: m
  end subroutine fxt_faltld_exp


end interface

end module fxt_faltld_binding
