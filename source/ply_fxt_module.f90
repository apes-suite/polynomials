module ply_fxt_module
  use env_module,                only: rk
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
    integer :: ndims
  end type ply_fxt_type 


  public :: ply_fxt_type  
  public :: ply_init_fxt
  public :: ply_fxt_m2n_1D, ply_fxt_m2n_2D,ply_fxt_m2n_3D
  public :: ply_fxt_n2m_1D, ply_fxt_n2m_2D,ply_fxt_n2m_3D

contains

  !****************************************************************************!
   !> Initialize the flpt data structure for fast legendre polynomial
   !! transformation via the fxtpack.
   subroutine ply_init_fxt(flpt, degree, nPoints, prec)
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

   end subroutine ply_init_fxt
  !****************************************************************************!


  !****************************************************************************!
  !> Convert modal data to nodal data using flpt.
   !!
   !! This encapsualtes the pure C-Interface, with extraction of the array
   !! sizes and dealing with the flpt data.
   !!
   !! Note: The modal and nodal data array sizes need to match the flpt
   !! definitions, provided in the fxtf_flptld_init call.
  !> Convert modal data to nodal data in 1D using flpt.
   subroutine ply_fxt_m2n_1D( fxt, modal_data, nodal_data, nModes, nNodes, &
     &                        oversamp_degree                          )
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=c_double), target :: modal_data(nModes)
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine

 
     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = ub + oversamp_degree
       call fxtf_flptld_n2m( flpt        = fxt%flpt,         & 
         &                   nodal_data = nodal_data(lb:ub), &
         &                   modal_data = modal_data(lb:ub), &
         &                   nModes     = nModes,            &
         &                   nNodes     = nNodes             )
     end do

   end subroutine ply_fxt_m2n_1D
  !****************************************************************************!

  !****************************************************************************!
  !> Convert modal data to nodal data in 2D using flpt.
   subroutine ply_fxt_m2n_2D(fxt, modal_data, nodal_data, nModes, nNodes, &
                          &  oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=c_double), target :: modal_data(nModes)
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim 

     nModesPerDim = (oversamp_degree+1)**2

     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = ub + oversamp_degree
       call fxtf_flptld_n2m( flpt        = fxt%flpt,         & 
         &                   nodal_data = nodal_data(lb:ub), &
         &                   modal_data = modal_data(lb:ub), &
         &                   nModes     = nModes,            &
         &                   nNodes     = nNodes             )
     end do

     do iColumn = 1, oversamp_degree+1
       lb = iColumn
       ub = (oversamp_degree*oversamp_degree-1) + iColumn
       call fxtf_flptld_n2m( flpt       = fxt%flpt,                    & 
         &                   nodal_data =                              &
         &                    nodal_data(lb:ub:oversamp_degree+1), &
         &            modal_data = modal_data(nModesPerDim+lb          &
         &                                   :nModesPerDim+ub          &
         &                                   :oversamp_degree+1),  &
         &                   nModes     = nModes,                      &
         &                   nNodes     = nNodes                       )
     end do
   end subroutine ply_fxt_m2n_2D
  !****************************************************************************!

  !****************************************************************************!
  !> Convert modal data to nodal data in 3D using flpt.
   subroutine ply_fxt_m2n_3D(fxt, modal_data, nodal_data, nModes, nNodes, &
                          &  oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=c_double), target :: modal_data(nModes)
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine

     call fxtf_flptld_evl( c_loc(nodal_data), nNodes, fxt%flpt%handle, &
       &                   c_loc(modal_data), nModes, fxt%flpt%work    )

     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = ub + oversamp_degree
       call fxtf_flptld_n2m( flpt        = fxt%flpt,         & 
         &                   nodal_data = nodal_data(lb:ub), &
         &                   modal_data = modal_data(lb:ub), &
         &                   nModes     = nModes,            &
         &                   nNodes     = nNodes             )
     end do
   end subroutine ply_fxt_m2n_3D
  !****************************************************************************!
  !****************************************************************************!

  !****************************************************************************!
   !> Convert nodal data to modal data using flpt.
   !!
   !! This encapsualtes the pure C-Interface, with extraction of the array
   !! sizes and dealing with the flpt data.
   !!
   !! Note: The modal and nodal data array sizes need to match the flpt
   !! definitions, provided in the fxtf_flptld_init call.
   subroutine ply_fxt_n2m_1D(fxt, nodal_data, modal_data, nNodes, nModes)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=c_double), target :: modal_data(nModes)

     call fxtf_flptld_exp( c_loc(modal_data), nModes, fxt%flpt%handle, &
       &                   c_loc(nodal_data), nNodes, fxt%flpt%work    )

   end subroutine ply_fxt_n2m_1D 
  !****************************************************************************!

   !> todo :NA: This routine needs to be adapted to work with the specified dimension
   subroutine ply_fxt_n2m_2D(fxt, nodal_data, modal_data, nNodes, nModes)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=c_double), target :: modal_data(nModes)

     call fxtf_flptld_exp( c_loc(modal_data), nModes, fxt%flpt%handle, &
       &                   c_loc(nodal_data), nNodes, fxt%flpt%work    )
   end subroutine ply_fxt_n2m_2D 
  !****************************************************************************!


   !> todo :NA: This routine needs to be adapted to work with the specified dimension
   subroutine ply_fxt_n2m_3D(fxt, nodal_data, modal_data, nNodes, nModes)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     integer(kind=c_int), intent(in) :: nModes, nNodes
     !> Nodal data
     real(kind=c_double), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=c_double), target :: modal_data(nModes)

     call fxtf_flptld_exp( c_loc(modal_data), nModes, fxt%flpt%handle, &
       &                   c_loc(nodal_data), nNodes, fxt%flpt%work    )

   end subroutine ply_fxt_n2m_3D 
  !****************************************************************************!


end module ply_fxt_module
