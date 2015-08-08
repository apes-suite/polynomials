module fxt_fwrap
  use, intrinsic :: iso_c_binding
  use env_module, only: rk
  use fxt_fif

  implicit none


  !> This datatype provides a handle to the information that FXTPACK needs
  !! to have about the transformation.
  type fxtf_flptld_type
    !> Handle for the fast Legendre polynomial transformation data in FXTPACK.
    type(c_ptr) :: handle

    !> Pointer to the working array, that is required by the transformations.
    type(c_ptr) :: work
  end type fxtf_flptld_type

!  public :: fxtf_flptld_type
!  public :: fxtf_flptld_evl 
!  public :: fxtf_flptld_exp


  !HK: This is totally unnecessary!
  !HK: The declaration of the interfaces to the C-Routines is in fxt_fif,
  !HK: which is used above!
  !HK:
  !HK: It is absolutely confusing to have two contradicting interface
  !HK: declarations in the same code!
  !HK: Also, I am not sure why the declarations in fxt_fif should be
  !HK: problematic!
  interface
    subroutine fxtf_flptld_evl(v, vn, flpt, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      real(c_double), dimension(*) :: v
      integer(c_int), value :: vn
      type(c_ptr), value :: flpt
      real(c_double), dimension(*) :: u
      integer(c_int), value :: un
      type(c_ptr), value :: w
    end subroutine fxtf_flptld_evl

    subroutine fxtf_flptld_exp(u, un, flpt, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      real(c_double), dimension(*)  ::  u
      integer(c_int), value  :: un
      type(c_ptr), value  :: flpt
      real(c_double), dimension(*)  :: v
      integer(c_int), value :: vn
      type(c_ptr), value  :: w
    end subroutine fxtf_flptld_exp

    subroutine fxtf_faltld_evl(v, vn, falt, m, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: falt
      integer(kind=c_long) :: m
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: w
    end subroutine fxtf_faltld_evl

    subroutine fxtf_faltld_exp(u, un, falt, m, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: falt
      integer(kind=c_long) :: m
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: w
    end subroutine fxtf_faltld_exp

!HK: Not used anywhere...
!HK: I would like to keep the modifications to the upstream code
!HK: to a minimum. So, it would be nice, if we could avoid using it.
!HK!    subroutine fxt_error_print() bind(c)
!HK!
!HK!    end subroutine fxt_error_print

  end interface

end module fxt_fwrap
