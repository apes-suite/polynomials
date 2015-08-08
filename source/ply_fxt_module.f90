!> Fast polynomial transformation using the FXTPACK implementation of a
!! fast multipole method.
module ply_fxt_module
  use env_module, only: rk
  use fxt_fwrap, only: fxtf_flptld_type, fxtf_flptld_evl, fxtf_flptld_exp

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
   !!
   !!\todo Actually implement filling of ply_fxt_type here.
   !!      Basically needs to set precision and dimension, and call
   !!      fxtf_flptld_init
   subroutine ply_init_fxt(flpt, nDims, degree, nPoints, prec)
    !--------------------------------------------------------------------------!
     !> Handle to the resulting fast polynomial table.
     type(ply_fxt_type), intent(out) :: flpt

     !> Number of dimensions to use for this transformation.
     integer, intent(in) :: nDims

     !> Polynomial degree.
     integer, intent(in) :: degree

     !> Number of points.
     !! Optional, defaults to degree+1.
     integer, intent(in), optional :: nPoints

     !> Required precision for the transformation.
     !! Optional, defaults to 8 times the precision of c_double.
     real(kind=rk), intent(in), optional :: prec


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
     !> Nodal data
     real(kind=rk), intent(inout) :: nodal_data(nNodes)
     !> Modal data
     real(kind=rk), intent(inout) :: modal_data(nModes)
     integer, intent(in) :: nModes, nNodes
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine
     !-----------------------------------------------------------!

 
     call fxtf_flptld_evl( nodal_data, nNodes,fxt%flpt%handle, &
       &                   modal_data, nModes,fxt%flpt%work    )


   end subroutine ply_fxt_m2n_1D
  !****************************************************************************!

  !****************************************************************************!
  !> Convert modal data to nodal data in 2D using flpt.
   subroutine ply_fxt_m2n_2D(fxt, modal_data, nodal_data, nModes, nNodes, &
                          &  oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), intent(inout) :: nodal_data(nNodes)
     !> Modal data
     real(kind=rk), intent(inout) :: modal_data(nModes)
     integer, intent(in) :: nModes, nNodes
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq 
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim

     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       call fxtf_flptld_evl( nodal_data(lb:ub), nModesPerDim,   &
         &                   fxt%flpt%handle,modal_data(lb:ub), &
         &                   nModesPerDim,fxt%flpt%work         )
     end do

     do iColumn = 1, oversamp_degree+1
       lb = iColumn
       ub = oversamp_degree +1
       call fxtf_flptld_evl( nodal_data(lb:msq:ub), nModesPerDim ,        &
         &                   fxt%flpt%handle, modal_data(lb : msq : ub),  &
         &                   nModesPerDim,fxt%flpt%work                   )
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
     !> Nodal data
     real(kind=rk), intent(inout) :: nodal_data(nNodes)
     !> Modal data
     real(kind=rk), intent(inout) :: modal_data(nModes)
     integer, intent(in) :: nModes, nNodes
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq, ntotalDofs
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim
     nTotalDofs =  (oversamp_degree+1)**3

     ! The loop for msq stripes for independent x Dir evaluations
     do iLine = 1, msq
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       call fxtf_flptld_evl( nodal_data(lb:ub),                 &
         &                   nModesPerDim,                      &
         &                   fxt%flpt%handle,                   &
         &                   modal_data(lb:ub),                 &
         &                   nModesPerDim,                      &
         &                   fxt%flpt%work                      )
     end do

     ! The loop for msq stripes for independent y Dir evaluations
     do iColumn = 1, msq
       lb = int((iColumn - 1)/ (oversamp_degree + 1))*msq + 1 
       ub = int((iColumn-1)/(oversamp_degree + 1))*msq + msq
       call fxtf_flptld_evl( nodal_data( lb : ub : oversamp_degree + 1),  &
         &                   nModesPerDim ,                               &
         &                   fxt%flpt%handle,                             &
         &                   modal_data(lb : ub : oversamp_degree + 1),   &
         &                   nModesPerDim,                                &
         &                   fxt%flpt%work                                )
     end do

     ! The loop for msq stripes for independent z Dir evaluations
     do iColumn = 1, msq
       lb = iColumn 
       ub = nTotalDofs
       call fxtf_flptld_evl( nodal_data( lb : ub : msq),                  &
         &                   nModesPerDim ,                               &
         &                   fxt%flpt%handle,                             &
         &                   modal_data(lb : ub : msq),                   &
         &                   nModesPerDim,                                &
         &                   fxt%flpt%work                                )
     end do
     


   end subroutine ply_fxt_m2n_3D
  !****************************************************************************!

  !****************************************************************************!
   !> Convert nodal data to modal data using flpt.
   !!
   !! This encapsualtes the pure C-Interface, with extraction of the array
   !! sizes and dealing with the flpt data.
   !!
   !! Note: The modal and nodal data array sizes need to match the flpt
   !! definitions, provided in the fxtf_flptld_init call.
   subroutine ply_fxt_n2m_1D(fxt, nodal_data, modal_data, nNodes, nModes, &
    &                        oversamp_degree                              )
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), intent(inout) :: nodal_data(nNodes)
     !> Modal data
     real(kind=rk), intent(inout) :: modal_data(nModes)
     integer, intent(in) :: nModes, nNodes
     integer, intent(in) :: oversamp_degree 
    !--------------------------------------------------------------------------!

     call fxtf_flptld_exp( modal_data, nModes, fxt%flpt%handle,   &
       &                   nodal_data, nNodes, fxt%flpt%work      )


   end subroutine ply_fxt_n2m_1D 
  !****************************************************************************!

   subroutine ply_fxt_n2m_2D(fxt, nodal_data, modal_data, nNodes, nModes,    &
    &                                                         oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target :: nodal_data(nNodes)
     !> Modal data
     real(kind=rk), target :: modal_data(nModes)
     integer, intent(in) :: nModes, nNodes
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq 
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim

     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       call fxtf_flptld_exp( modal_data(lb:ub), nModesPerDim ,          &
         &                   fxt%flpt%handle, nodal_data(lb:ub),        &
         &                   nModesPerDim ,fxt%flpt%work                )
     end do

     do iColumn = 1, oversamp_degree+1
       lb = iColumn
       ub = oversamp_degree + 1
       call fxtf_flptld_exp( modal_data(lb:msq:ub), nModesPerDim ,      &
         &                   fxt%flpt%handle, nodal_data(lb : msq :ub), &
         &                   nModesPerDim ,fxt%flpt%work                )
     end do
   end subroutine ply_fxt_n2m_2D 
  !****************************************************************************!


   !> todo :NA: This routine needs to be adapted to work with the specified dimension
   subroutine ply_fxt_n2m_3D(fxt, nodal_data, modal_data, nNodes, nModes,  &
    &                        oversamp_degree                               )
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), intent(inout) :: nodal_data(nNodes)
     !> Modal data
     real(kind=rk), intent(inout) :: modal_data(nModes)
     integer, intent(in) :: nModes, nNodes
     integer, intent(in) :: oversamp_degree 
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq, ntotalDofs
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim
     nTotalDofs =  (oversamp_degree+1)**3

     ! The loop for msq stripes for independent x Dir evaluations
     do iLine = 1, msq
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       call fxtf_flptld_exp( nodal_data(lb:ub),                 &
         &                   nModesPerDim,                      &
         &                   fxt%flpt%handle,                   &
         &                   modal_data(lb:ub),                 &
         &                   nModesPerDim,                      &
         &                   fxt%flpt%work                      )
     end do

     ! The loop for msq stripes for independent y Dir evaluations
     do iColumn = 1, msq
       lb = int((iColumn - 1)/ (oversamp_degree + 1))*msq + 1 
       ub = int((iColumn-1)/(oversamp_degree + 1))*msq + msq
       call fxtf_flptld_exp( nodal_data( lb : ub : oversamp_degree + 1),  &
         &                   nModesPerDim ,                               &
         &                   fxt%flpt%handle,                             &
         &                   modal_data(lb : ub : oversamp_degree + 1),   &
         &                   nModesPerDim,                                &
         &                   fxt%flpt%work                                )
     end do

     ! The loop for msq stripes for independent z Dir evaluations
     do iColumn = 1, msq
       lb = iColumn 
       ub = nTotalDofs
       call fxtf_flptld_exp( nodal_data( lb : ub : msq),                  &
         &                   nModesPerDim ,                               &
         &                   fxt%flpt%handle,                             &
         &                   modal_data(lb : ub : msq),                   &
         &                   nModesPerDim,                                &
         &                   fxt%flpt%work                                )
     end do

   end subroutine ply_fxt_n2m_3D 
  !****************************************************************************!


end module ply_fxt_module
