!> Fast polynomial transformation using the FXTPACK implementation of a
!! fast multipole method.
module ply_fxt_module
  use env_module, only: rk
  use fxt_fwrap,  only: fxtf_flptld_type, fxtf_flptld_n2m, fxtf_flptld_m2n,&
    &                  fxtf_flptld_init
  use ply_fxt_header_module,     only: ply_fxt_header_type
  use ply_nodes_module,          only: init_gauss_nodes, ply_faceNodes_type, &
    &                                  init_gauss_nodes_2d, init_gauss_nodes_1d
  use ply_space_integration_module, only: ply_create_surface_gauss_points_cube, &
    &                                     ply_create_surface_gauss_points_cube_2d, &
    &                                     ply_create_surface_gauss_points_cube_1d, &
    &                                     ply_gaussLegPoints
  use ply_modg_basis_module,     only: legendre_1D

  implicit none

  private

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
   subroutine ply_init_fxt(fxt, header, degree, nDims, nodes, faces)
     !-------------------------------------------------------------------------!
     !> Handle to the resulting fast polynomial table.
     type(ply_fxt_type), intent(out) :: fxt
     type(ply_fxt_header_type), intent(in)     :: header
      !> Polynomial degree.
     integer, intent(in)                       :: degree
      !> Number of dimensions to use for this transformation.
     integer, intent(in)                       :: nDims
     real(kind=rk), intent(out), allocatable   :: nodes(:,:)
     type(ply_faceNodes_type), intent(out), allocatable :: faces(:,:)
     !-------------------------------------------------------------------------!
     real(kind=rk), allocatable :: tmp_weights(:)
     real(kind=rk), allocatable :: gaussp1D(:)
     real(kind=rk), allocatable :: weights1D(:)
     real(kind=rk), allocatable :: leg1D_at_gauss(:,:)
     integer :: iDir, iAlign, nPoints, nDofs, iPoint, lb, ub
     !-------------------------------------------------------------------------!
     nPoints = degree + 1
     nDofs = nPoints
     allocate(gaussp1D(npoints))
     allocate(weights1D(nPoints))
     allocate(leg1D_at_gauss(max(2, nDofs), nPoints))

     ! create the quadrature points for volume and on the face
     ! for the oversampled projection
     call ply_gaussLegPoints( x1 = -1.0_rk, x2 = 1.0_rk, nIntP = nPoints, &
       &          w = weights1D, x = gaussp1D                    )

     leg1D_at_gauss = legendre_1D(gaussp1D, degree)

     select case(nDims)
     case(3) 
       call fxtf_flptld_init( flpt    = fxt%flpt,           &
          &                   degree  = degree,             &
          &                   nPoints = degree+1,           &
          &                   prec    = header%prec         )

       allocate( nodes(nPoints**3, 3) )
       do iPoint=1,nPoints
         lb = (iPoint-1)*nPoints + 1
         ub = iPoint*nPoints
         nodes(lb:ub,1) = gaussP1D
         nodes(lb:ub,2) = gaussP1D(iPoint)
       end do
       nodes(:nPoints**2,3) = gaussP1D(1)

       do iPoint=2,nPoints
         lb = (iPoint-1)*nPoints**2 + 1
         ub = iPoint*nPoints**2
         nodes(lb:ub,1) = nodes(:nPoints**2,1)
         nodes(lb:ub,2) = nodes(:nPoints**2,2)
         nodes(lb:ub,3) = gaussP1D(iPoint)
       end do

       allocate( faces(3,2) )
       do iDir = 1,3
         do iAlign = 1,2
           faces(iDir,iAlign)%nquadpoints = nPoints**2
           call ply_create_surface_gauss_points_cube(                 &
             &    num_intp_per_direction = nPoints,                   &
             &    points                 = faces(iDir,iAlign)%points, &
             &    weights                = tmp_weights,               &
             &    refElemMin             = -1.0_rk,                   &
             &    refElemMax             =  1.0_rk,                   &
             &    dir                    = idir,                      &
             &    align                  = iAlign                     )
           deallocate(tmp_weights)
         end do
       end do

     case(2) 
       call fxtf_flptld_init(flpt    = fxt%flpt,            &
          &                   degree  = degree,             &
          &                   nPoints = degree+1,           &
          &                   prec    = header%prec         )  
 
       ! Fill up the nodes and the face with gauss legendre points   
       allocate( nodes(nPoints**2, 3) )
       do iPoint=1,nPoints
         lb = (iPoint-1)*nPoints + 1
         ub = iPoint*nPoints
         nodes(lb:ub,1) = gaussP1D
         nodes(lb:ub,2) = gaussP1D(iPoint)
       end do

       allocate( faces(2,2) )
       do iDir = 1,2
         do iAlign = 1,2
           faces(iDir,iAlign)%nquadpoints = nPoints
           call ply_create_surface_gauss_points_cube_2d(              &
             &    num_intp_per_direction = nPoints,                   &
             &    points                 = faces(iDir,iAlign)%points, &
             &    weights                = tmp_weights,               &
             &    refElemMin             = -1.0_rk,                   &
             &    refElemMax             =  1.0_rk,                   &
             &    dir                    = idir,                      &
             &    align                  = iAlign                     )
           deallocate(tmp_weights)
         end do
       end do

     case(1) 
       call fxtf_flptld_init(flpt    = fxt%flpt,           &
         &                   degree  = degree,             &
         &                   nPoints = degree+1,           &
         &                   prec    = header%prec         )
      allocate( nodes(degree + 1, 3) )
      nodes(:,1) = gaussP1D
      nodes(:,2:) = 0.0_rk

      allocate( faces(1,2) )
      iDir = 1
      do iAlign = 1,2
        faces(iDir,iAlign)%nquadpoints = 1
        call ply_create_surface_gauss_points_cube_1d(   &
          &    points  = faces(iDir,iAlign)%points,     &
          &    weights = tmp_weights,                   &
          &    dir     = idir,                          &
          &    align   = iAlign                         )
        deallocate(tmp_weights)
      end do
     end select

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
   subroutine ply_fxt_m2n_1D( fxt, modal_data, nodal_data, oversamp_degree )
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target, intent(inout) :: nodal_data(:)
     !> Modal data
     real(kind=rk), target, intent(inout) :: modal_data(:)
     integer, intent(in) :: oversamp_degree
     !-----------------------------------------------------------!

     call fxtf_flptld_m2n(flpt       = fxt%flpt,   &
       &                  modal_data = modal_data, &
       &                  nodal_data = nodal_data  )

   end subroutine ply_fxt_m2n_1D
  !****************************************************************************!

  !****************************************************************************!
  !> Convert modal data to nodal data in 2D using flpt.
   subroutine ply_fxt_m2n_2D(fxt, modal_data, nodal_data, oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target, intent(inout) :: nodal_data(:)
     !> Modal data
     real(kind=rk), target, intent(inout) :: modal_data(:)
     integer, intent(in) :: oversamp_degree
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim

     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       call fxtf_flptld_m2n(flpt       = fxt%flpt,          &
         &                  modal_data = modal_data(lb:ub), &
         &                  nodal_data = nodal_data(lb:ub)  )
     end do

     do iColumn = 1, oversamp_degree+1
       lb = iColumn
       call fxtf_flptld_m2n(flpt       = fxt%flpt,                             &
         &                  modal_data = nodal_data(lb:msq:oversamp_degree+1), &
         &                  nodal_data = modal_data(lb:msq:oversamp_degree+1)  )
     end do
     nodal_data = modal_data

   end subroutine ply_fxt_m2n_2D
  !****************************************************************************!

  !****************************************************************************!
  !> Convert modal data to nodal data in 3D using flpt.
   subroutine ply_fxt_m2n_3D(fxt, modal_data, nodal_data, oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target, intent(inout) :: nodal_data(:)
     !> Modal data
     real(kind=rk), target, intent(inout) :: modal_data(:)
     integer, intent(in) :: oversamp_degree
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq, ntotalDofs
     real(kind=rk), pointer :: tmp_in(:), tmp_out(:)
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim
     nTotalDofs =  (oversamp_degree+1)**3
     allocate(tmp_in(nModesPerDim))
     allocate(tmp_out(nModesPerDim))
     tmp_in = -42
     tmp_out = -42

     ! The loop for msq stripes for independent x Dir evaluations
     do iLine = 1, msq
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       tmp_in = modal_data(lb:ub)
       call fxtf_flptld_m2n(flpt       = fxt%flpt, &
         &                  modal_data = tmp_in,   &
         &                  nodal_data = tmp_out   )
       nodal_data(lb:ub) = tmp_out
     end do

     ! The loop for msq stripes for independent y Dir evaluations
     do iColumn = 1, msq
       lb = int((iColumn-1)/nModesPerDim)*msq + mod(iColumn-1,nModesPerDim)+1
       ub = lb + msq - 1
       call fxtf_flptld_m2n(flpt       = fxt%flpt,                       &
         &                  modal_data = nodal_data(lb:ub:nModesPerDim), &
         &                  nodal_data = modal_data(lb:ub:nModesPerDim)  )
     end do

     ! The loop for msq stripes for independent z Dir evaluations
     ub = nTotalDofs
     do iColumn = 1, msq
       lb = iColumn
       call fxtf_flptld_m2n(flpt       = fxt%flpt,              &
         &                  modal_data = modal_data(lb:ub:msq), &
         &                  nodal_data = nodal_data(lb:ub:msq)  )
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
   subroutine ply_fxt_n2m_1D(fxt, nodal_data, modal_data, oversamp_degree )
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target, intent(inout) :: nodal_data(:)
     !> Modal data
     real(kind=rk), target, intent(inout) :: modal_data(:)
     integer, intent(in) :: oversamp_degree
    !--------------------------------------------------------------------------!

     call fxtf_flptld_n2m(flpt       = fxt%flpt,   &
       &                  nodal_data = nodal_data, &
       &                  modal_data = modal_data  )

   end subroutine ply_fxt_n2m_1D
  !****************************************************************************!

   subroutine ply_fxt_n2m_2D(fxt, nodal_data, modal_data, oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target :: nodal_data(:)
     !> Modal data
     real(kind=rk), target :: modal_data(:)
     integer, intent(in) :: oversamp_degree
     !-----------------------------------------------------------!
     integer :: ub, lb, iLine, iColumn, nModesPerDim, msq
     !-----------------------------------------------------------!

     nModesPerDim = (oversamp_degree+1)
     msq = nModesPerDim*nModesPerDim

     do iLine = 1, oversamp_degree+1
       lb = (iLine-1)*(oversamp_degree+1)+1
       ub = lb + oversamp_degree
       call fxtf_flptld_n2m(flpt       = fxt%flpt,          &
         &                  nodal_data = nodal_data(lb:ub), &
         &                  modal_data = modal_data(lb:ub)  )
     end do

     do iColumn = 1, oversamp_degree+1
       lb = iColumn
       call fxtf_flptld_n2m(flpt       = fxt%flpt,          &
         &                  nodal_data = modal_data(lb:msq:oversamp_degree+1), &
         &                  modal_data = nodal_data(lb:msq:oversamp_degree+1)  )
     end do
     modal_data = nodal_data
   end subroutine ply_fxt_n2m_2D
  !****************************************************************************!


   subroutine ply_fxt_n2m_3D(fxt, nodal_data, modal_data, oversamp_degree)
    !--------------------------------------------------------------------------!
     !> Description of the Fast Legendre Polynomial Transform
     type(ply_fxt_type) :: fxt
     !> Nodal data
     real(kind=rk), target, intent(inout) :: nodal_data(:)
     !> Modal data
     real(kind=rk), target, intent(inout) :: modal_data(:)
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
       call fxtf_flptld_n2m(flpt       = fxt%flpt,          &
         &                  nodal_data = nodal_data(lb:ub), &
         &                  modal_data = modal_data(lb:ub)  )
     end do

     ! The loop for msq stripes for independent y Dir evaluations
     do iColumn = 1, msq
       lb = int((iColumn-1)/nModesPerDim)*msq + mod(iColumn-1,nModesPerDim)+1
       ub = lb + msq - 1
       call fxtf_flptld_n2m(flpt       = fxt%flpt,                       &
         &                  nodal_data = modal_data(lb:ub:nModesPerDim), &
         &                  modal_data = nodal_data(lb:ub:nModesPerDim)  )
     end do

     ! The loop for msq stripes for independent z Dir evaluations
     ub = nTotalDofs
     do iColumn = 1, msq
       lb = iColumn
       call fxtf_flptld_n2m(flpt       = fxt%flpt,              &
         &                  nodal_data = nodal_data(lb:ub:msq), &
         &                  modal_data = modal_data(lb:ub:msq)  )
     end do

   end subroutine ply_fxt_n2m_3D
  !****************************************************************************!


end module ply_fxt_module
