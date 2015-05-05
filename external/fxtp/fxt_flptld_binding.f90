module fxt_flptld_binding
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE 

interface
    subroutine fxt_flptld_preproc(p, n, prec, fname) bind (c)
      use, intrinsic :: iso_c_binding
      use fxt_types
      integer(c_long), value :: p
      integer(c_long), value :: n
      real(c_double), value :: prec
      character(c_char) :: fname
    end subroutine fxt_flptld_preproc
  
    ! deallocate fast Legendre polynomial transform
    ! void fxt_flptld_del(fxt_flptld *flpt);
    subroutine fxt_flptld_del(flpt) bind (c)
      use, intrinsic :: iso_c_binding
      use fxt_types
      type(c_ptr) :: flpt             !fxt_flptld
    end subroutine fxt_flptld_del
  
    ! size of working array
    ! long fxt_flptld_wsize(fxt_flptld *flpt); 
    integer(c_long) function fxt_flptld_wsize(flpt) bind (c)
      use, intrinsic :: iso_c_binding
      use fxt_types
      type(c_ptr) :: flpt              !fxt_flptld
    end function fxt_flptld_wsize
   
    ! evaluate fast Legendre Polynomial transform 
    ! void fxt_flptld_evl(fxt_vecld *v, fxt_flptld *flpt,
    !                     fxt_vecld *u, fxt_vecld *w);
    subroutine fxt_flptld_evl(v, flpt, u, w) bind (c)
      use, intrinsic :: iso_c_binding
      use fxt_types
      type(c_ptr) :: v, u, w              ! fxt_vecld
      type(c_ptr) :: flpt               ! fxt_flptld
    end subroutine fxt_flptld_evl
  
    ! expand fast Legendre Polynomial transform 
    ! void fxt_flptld_exp(fxt_vecld *u, fxt_flptld *flpt,
    !                     fxt_vecld *v, fxt_vecld *w);  
    subroutine fxt_flptld_exp(u, flpt, v, w) bind(c)
      use, intrinsic :: iso_c_binding
      use fxt_types
      type(c_ptr) ::  u, v, w        ! fxt_vecld
      type(c_ptr) :: flpt            ! fxt_flptld
    end subroutine fxt_flptld_exp
  end interface
  
  end module fxt_flptld_binding
