module ply_fxt_module
  use env_module,                only: rk
  use fxt_fif
  use fxt_fwrap

  implicit none

  private
 
!  !> This datatype provides a handle to the information that FXTPACK needs
!  !! to have about the transformation.
!  type fxtf_flptld_type
!    !> Handle for the fast Legendre polynomial transformation data in FXTPACK.
!    type(c_ptr) :: handle
!
!    !> Pointer to the working array, that is required by the transformations.
!    type(c_ptr) :: work
!  end type fxtf_flptld_type

  type ply_fxt_type
    type(fxtf_flptld_type) :: flpt
    real(kind=rk) :: prec
  end type ply_fxt_type 


  public :: ply_fxt_type  
  public :: ply_fxt_init, ply_fxt_m2n, ply_fxt_n2m

contains

  !****************************************************************************!
   !> Initialize the flpt data structure for fast legendre polynomial
   !! transformation via the fxtpack.
   subroutine ply_fxt_init(flpt, degree, nPoints, prec)
    !--------------------------------------------------------------------------!
     !> Handle to the resulting fast polynomial table.
     type(fxtf_flptld_type), intent(out) :: flpt

     !> Polynomial degree.
     integer, intent(in) :: degree

     !> Number of points.
     !! Optional, defaults to degree+1.
     integer, intent(in), optional :: nPoints

     !> Required precision for the transformation.
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

   end subroutine ply_fxt_init
  !****************************************************************************!


  !****************************************************************************!
  !> Convert modal data to nodal data using flpt.
   !!
   !! This encapsualtes the pure C-Interface, with extraction of the array
   !! sizes and dealing with the flpt data.
   !!
   !! Note: The modal and nodal data array sizes need to match the flpt
   !! definitions, provided in the fxtf_flptld_init call.
   subroutine ply_fxt_m2n(fxt, modal_data, nodal_data, nModes, nNodes)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(:,:)
     !> Modal data
     real(kind=c_double), target :: modal_data(:,:)

     call fxtf_flptld_evl( c_loc(nodal_data), nNodes, fxt%flpt%handle, &
       &                   c_loc(modal_data), nModes, fxt%flpt%work    )

   end subroutine ply_fxt_m2n
  !****************************************************************************!

  !****************************************************************************!
   !> Convert nodal data to modal data using flpt.
   !!
   !! This encapsualtes the pure C-Interface, with extraction of the array
   !! sizes and dealing with the flpt data.
   !!
   !! Note: The modal and nodal data array sizes need to match the flpt
   !! definitions, provided in the fxtf_flptld_init call.
   subroutine ply_fxt_n2m(fxt, nodal_data, modal_data, nNodes, nModes)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(:,:)
     !> Modal data
     real(kind=c_double), target :: modal_data(:,:)

     integer(kind=c_int) :: un
     integer(kind=c_int) :: vn

     vn = size(nodal_data)
     un = size(modal_data)

     call fxtf_flptld_exp( c_loc(modal_data), un, fxt%flpt%handle, &
       &                   c_loc(nodal_data), vn, fxt%flpt%work    )

   end subroutine ply_fxt_n2m 
  !****************************************************************************!

end module ply_fxt_module
