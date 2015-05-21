module fxt_fwrap
  use, intrinsic :: iso_c_binding
  use :: fxt_fif

  implicit none

  !> This datatype provides a handle to the information that FXTPACK needs
  !! to have about the transformation.
  type fxtf_flptld
    !> Handle for the fast Legendre polynomial transformation data in FXTPACK.
    type(c_ptr) :: handle

    !> Pointer to the working array, that is required by the transformations.
    type(c_ptr) :: work
  end type fxtf_flptld  


  interface
    subroutine fxtf_flptld_evl(v, vn, fplt, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: fplt
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: w
    end subroutine fxtf_flptld_evl

    subroutine fxtf_flptld_exp(u, un, fplt, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: fplt
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: w
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

  end interface
 

contains
  

  !> Initialize the flpt data structure for fast legendre polynomial
  !! transformation via the fxtpack.
  subroutine fxtf_flptld_init(flpt, degree, nPoints, prec)
    !> Handle to the resulting fast polynomial table.
    type(fxtf_flptld), intent(out) :: flpt

    !> Polynomial degree.
    integer, intent(in) :: degree

    !> Number of points.
    !!
    !! Optional, defaults to degree+1.
    integer, intent(in), optional :: nPoints

    !> Required precision for the transformation.
    !!
    !! Optional, defaults to 8 times the precision of c_double.
    real(kind=c_double), intent(in), optional :: prec

    integer(c_long) :: wsize
    integer(c_long) :: p
    integer(c_long) :: n

    real(kind=c_double) :: lprec

    n = degree
    if (present(nPoints)) then
      p = nPoints
    else
      p = degree + 1
    end if

    if (present(prec)) then
      lprec = prec
    else
      lprec = 8*epsilon(lprec)
    end if

    flpt%handle = fxt_flptld_init(p, n, lprec)
    wsize = fxt_flptld_wsize(flpt%handle)
    flpt%work = fxt_vecld_new(wsize)

  end subroutine fxtf_flptld_init

end module fxt_fwrap
