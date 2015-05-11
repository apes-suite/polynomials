!> This module provides the ISO_C_Binding interfaces to the fxtpack routines.
module fxt_fif
  implicit none

  interface

    ! FALTLD routines:
    ! /*** load fast spherical harmonic transform ***/
    ! fxt_faltld* fxt_faltld_load(char *fname)
    function fxt_faltld_load(fname) result(falt) bind(c)
      use, intrinsic :: iso_c_binding
      character(c_char) :: fname
      type(c_ptr) :: falt
    end function fxt_faltld_load

    ! deallocate fast spherical harmonic transform
    ! void fxt_faltld_del(fxt_faltld *falt);
    subroutine fxt_faltld_del(falt) bind (c)
      use, intrinsic :: iso_c_binding
     type(c_ptr), value :: falt           !fxt_faltld
    end subroutine fxt_faltld_del

    ! size of working array
    ! long fxt_faltld_wsize(fxt_faltld *falt, long m);
    integer(c_long) function fxt_faltld_wsize(falt, m) bind (c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: falt            !fxt_faltld
      integer(c_long), value :: m
    end function fxt_faltld_wsize

    !  /*** maximum size of working array ***/
    ! long fxt_faltld_wsizemax(fxt_faltld *falt);
    integer(c_long) function fxt_faltld_wsizemax(falt) bind (c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: falt          !fxt_faltld
    end function fxt_faltld_wsizemax

    ! /*** evaluate fast spherical harmonic transform ***/
    ! void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
    !                     fxt_vecld *u, fxt_vecld *w);
    subroutine fxt_faltld_evl(v, falt, m, u, w) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: u, v, w                 !fxt_vecld
       type(c_ptr) :: falt                   !fxt_faltld
       integer(c_long), value :: m
    end subroutine fxt_faltld_evl

    !  /*** expand fast spherical harmonic transform ***/
    ! void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
    !                     fxt_vecld *v, fxt_vecld *w);
    subroutine fxt_faltld_exp(u, falt, m, v, w) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: u, v, w            !fxt_vecld
       type(c_ptr) :: falt              !fxt_faltld
       integer(c_long), value :: m
    end subroutine fxt_faltld_exp
    ! ------------------------------------------------------------------------ !


    ! ........................................................................ !
    ! FLPTLD routines:
    subroutine fxt_flptld_preproc(p, n, prec, fname) bind (c)
      use, intrinsic :: iso_c_binding
      integer(c_long), value :: p
      integer(c_long), value :: n
      real(c_double), value :: prec
      character(c_char) :: fname
    end subroutine fxt_flptld_preproc

    ! deallocate fast Legendre polynomial transform
    ! void fxt_flptld_del(fxt_flptld *flpt);
    subroutine fxt_flptld_del(flpt) bind (c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: flpt             !fxt_flptld
    end subroutine fxt_flptld_del

    ! size of working array
    ! long fxt_flptld_wsize(fxt_flptld *flpt);
    integer(c_long) function fxt_flptld_wsize(flpt) bind (c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: flpt              !fxt_flptld
    end function fxt_flptld_wsize

    ! evaluate fast Legendre Polynomial transform
    ! void fxt_flptld_evl(fxt_vecld *v, fxt_flptld *flpt,
    !                     fxt_vecld *u, fxt_vecld *w);
    subroutine fxt_flptld_evl(v, flpt, u, w) bind (c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: v, u, w              ! fxt_vecld
      type(c_ptr) :: flpt               ! fxt_flptld
    end subroutine fxt_flptld_evl

    ! expand fast Legendre Polynomial transform
    ! void fxt_flptld_exp(fxt_vecld *u, fxt_flptld *flpt,
    !                     fxt_vecld *v, fxt_vecld *w);
    subroutine fxt_flptld_exp(u, flpt, v, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr) ::  u, v, w        ! fxt_vecld
      type(c_ptr) :: flpt            ! fxt_flptld
    end subroutine fxt_flptld_exp
    ! ------------------------------------------------------------------------ !


    ! ........................................................................ !
    ! VECL / VECLD
    type(c_ptr) function fxt_vecl_new(size) bind (c)
      use, intrinsic :: iso_c_binding
      integer(c_long) :: size
    end function fxt_vecl_new

    type(c_ptr) function fxt_vecld_new(size) bind (c)
      use, intrinsic :: iso_c_binding
      integer(c_long) :: size
    end function fxt_vecld_new
    ! ------------------------------------------------------------------------ !

  end interface

end module fxt_fif
