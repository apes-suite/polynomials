module fxt_faltld_binding
  use, intrinsic :: iso_c_binding
  use fxt_types
  IMPLICIT NONE 

interface
 
!  /*** load fast spherical harmonic transform ***/
!  fxt_faltld* fxt_faltld_load(char *fname)
   function fxt_faltld_load(fname) result(falt) bind(c)
     use, intrinsic :: iso_c_binding
     character(c_char) :: fname
     type(c_ptr) :: falt
   end function fxt_faltld_load
 
  ! deallocate fast spherical harmonic transform
  ! void fxt_faltld_del(fxt_faltld *falt);
  subroutine fxt_faltld_del(falt) bind (c)
    use, intrinsic :: iso_c_binding
    use fxt_types
   type(c_ptr), value :: falt           !fxt_faltld
  end subroutine fxt_faltld_del

  ! size of working array
  ! long fxt_faltld_wsize(fxt_faltld *falt, long m); 
  integer(c_long) function fxt_faltld_wsize(falt, m) bind (c)
    use, intrinsic :: iso_c_binding
    use fxt_types
    type(c_ptr) :: falt            !fxt_faltld
    integer(c_long), value :: m
  end function fxt_faltld_wsize

  !  /*** maximum size of working array ***/
  ! long fxt_faltld_wsizemax(fxt_faltld *falt);
  integer(c_long) function fxt_faltld_wsizemax(falt) bind (c)
     use, intrinsic :: iso_c_binding
     use fxt_types
     type(c_ptr) :: falt          !fxt_faltld
  end function fxt_faltld_wsizemax

  ! /*** evaluate fast spherical harmonic transform ***/
  ! void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
  !                     fxt_vecld *u, fxt_vecld *w);
  subroutine fxt_faltld_evl(v, falt, m, u, w) bind(c)
     use, intrinsic :: iso_c_binding
     use fxt_types
     ! type(c_ptr), value :: v, falt, u, w
     type(c_ptr) :: u, v, w                 !fxt_vecld
     type(c_ptr) :: falt                   !fxt_faltld 
     integer(c_long), value :: m
  end subroutine fxt_faltld_evl

  !  /*** expand fast spherical harmonic transform ***/
  ! void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
  !                     fxt_vecld *v, fxt_vecld *w);
  subroutine fxt_faltld_exp(u, falt, m, v, w) bind(c)
     use, intrinsic :: iso_c_binding
     use fxt_types
     type(c_ptr) :: u, v, w            !fxt_vecld 
     type(c_ptr) :: falt              !fxt_faltld
     integer(c_long), value :: m
  end subroutine fxt_faltld_exp

end interface

end module fxt_faltld_binding
