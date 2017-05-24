module ply_l2p_module
  use env_module,                only: rk, labelLen

  use aotus_module,              only: flu_State, aot_get_val
  use aot_table_module,          only: aot_table_open, aot_table_close

  use tem_compileconf_module,    only: vlen

  use ply_modg_basis_module,     only: scalProdLeg
  use ply_space_integration_module, only:ply_create_surface_gauss_points_cube, &
    &                                    ply_create_surface_gauss_points_cube_2d, &
    &                                    ply_create_surface_gauss_points_cube_1d, &
    &                                    ply_gaussLegPoints
  use ply_modg_basis_module,     only: legendre_1D
  use ply_l2p_header_module,     only: ply_l2p_header_type
  use ply_nodes_module,          only: init_gauss_nodes, ply_faceNodes_type, &
    &                                  init_gauss_nodes_2d, init_gauss_nodes_1d

  implicit none

  private

  !> Storage of the transformation matrices for the L2 projection method to
  !! convert between modal and nodal values.
  type ply_l2p_type
    real(kind=rk), allocatable :: leg2node(:,:)
    real(kind=rk), allocatable :: node2leg(:,:)
  end type ply_l2p_type
 

  interface assignment(=)
    module procedure Copy_ply_l2p
  end interface

  public :: ply_l2p_type
  public :: ply_init_l2p
  public :: assignment(=)
  public :: ply_l2p_trafo_1D, ply_l2p_trafo_2D, ply_l2p_trafo_3D

contains

  !****************************************************************************!
  subroutine Copy_ply_l2p(left,right)
    !--------------------------------------------------------------------------!
    !> fpt to copy to
    !type(ply_legFpt_2D_type), intent(out) :: left
    type(ply_l2p_type), intent(out) :: left
    !> fpt to copy from
    !type(ply_legFpt_2D_type), intent(in) :: right
    type(ply_l2p_type), intent(in) :: right
    !--------------------------------------------------------------------------!

    left%leg2node = right%leg2node
    left%node2leg = right%node2leg

  end subroutine Copy_ply_l2p
  !****************************************************************************!


  !****************************************************************************!
  !> Initialize the transformations via L2 projections.
  subroutine ply_init_l2p(l2p, degree, nDims, nodes, faces)
    !--------------------------------------------------------------------------!
    type(ply_l2p_type), intent(out)       :: l2p
    integer, intent(in)                       :: degree
    integer, intent(in)                       :: nDims
    real(kind=rk), intent(out), allocatable   :: nodes(:,:)
    type(ply_faceNodes_type), intent(out), allocatable :: faces(:,:)
    !--------------------------------------------------------------------------!
    integer :: iPoint, iDof
    integer :: iDir, iAlign
    integer :: nDofs
    integer :: nPoints
    integer :: lb,ub
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: leg1D_at_gauss(:,:)
    real(kind=rk), allocatable :: weights1D(:)
    real(kind=rk), allocatable :: tmp_weights(:)
    real(kind=rk) :: scalprod_q
    !-------------------------------------------------------------------------!

    nDofs = degree+1
    nPoints = nDofs

    allocate(leg1D_at_gauss(max(2,nDofs), nPoints))
    allocate(gaussp1D(nPoints))
    allocate(weights1D(nPoints))

    ! create the quadrature points for volume and on the face
    ! for the oversampled projection
    call ply_gaussLegPoints( x1 = -1.0_rk, x2 = 1.0_rk, nIntP = nPoints, &
      &          w = weights1D, x = gaussp1D                 )
    leg1D_at_gauss = legendre_1D(gaussp1D, degree)

    allocate(l2p%leg2node(nDofs, nPoints))
    allocate(l2p%node2leg(nPoints, nDofs))

    ! Coefficients to transform legendre to nodal values.
    l2p%leg2node = leg1D_at_gauss(:nDofs, :)

    ! Coefficients to transform nodal to legendre values.
    do iDoF = 1,nDofs
      scalProd_q = 1.0_rk / scalProdLeg(iDoF)
      do iPoint = 1,nPoints
        l2p%node2leg(iPoint, iDoF) = l2p%leg2node(iDoF, iPoint) &
          &                          * weights1D(iPoint)        &
          &                          * scalProd_q
      end do
    end do

    ! Build the Gauss-Legendre nodes on the reference faces
    ! HK: This should probably moved to a separate routine.
    select case(nDims)
    case (1)
      allocate( nodes(nPoints, 3) )
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

    case (2)
      allocate( nodes(nPoints**2, 3) )
      do iPoint=1,nPoints
        lb = (iPoint-1)*nPoints + 1
        ub = iPoint*nPoints
        nodes(lb:ub,1) = gaussP1D
        nodes(lb:ub,2) = gaussP1D(iPoint)
      end do
      nodes(:,3) = 0.0_rk

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

    case (3)
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

    end select

  end subroutine ply_init_l2p
  !****************************************************************************!


  !****************************************************************************!
  !> Actual implementation of the matrix operation to change between nodal
  !! and modal representations.
  !!
  !! The operation is applied simultaneously to nIndeps 1D sections of the
  !! provided original data.
  !! These 1D sections have to run fastest in the original array and will be
  !! transposed (running slowest in the projected array).
  !! The actual direction of the operation depends on the passed matrix.
  !! matrix = l2p%leg2node will do the modal to nodal transformation
  !! matrix = l2p%node2leg will do the nodal to modal transformation
  subroutine ply_l2_projection(nIndeps, projected, original, matrix)
    !--------------------------------------------------------------------------!
    !> Number of values that can be computed independently.
    integer, intent(in) :: nIndeps

    !> Projected data.
    !!
    !! Size has to be nIndeps*size(matrix,1), and the layout is changed here
    !! when compared to the original array, as the projected direction moves
    !! to the end.
    real(kind=rk), intent(out) :: projected(:)

    !> Original data.
    !!
    !! Size has to be size(matrix,1) and the direction to be projected has to
    !! be the fastest running one.
    real(kind=rk), intent(in) :: original(:)

    !> Matrix to apply in this operation.
    !!
    !! The matrix defines wether this is a modal to nodal transformation or the
    !! other way around.
    real(kind=rk), intent(in) :: matrix(:,:)
    !--------------------------------------------------------------------------!
    integer :: i
    integer :: iStrip
    integer :: indep
    integer :: iDof
    integer :: iOrig
    integer :: strip_ub
    integer :: nDofs
    integer :: orig_off
    real(kind=rk) :: res(vlen)
    !--------------------------------------------------------------------------!

    nDofs = size(matrix,1)

    if (nDofs > 1) then

      !$OMP DO
      do iStrip=0,nIndeps-1,vlen

        ! Calculate the upper bound of the current strip
        strip_ub = min(iStrip + vlen, nIndeps) - iStrip

        ! Naive dense matrix vector multiply here.
        ! Maybe use DGEMV from BLAS here instead (though not straight forward
        ! with the memory layout used here).
        do iDof=1,nDofs
          do i=1,strip_ub
            indep = iStrip + i
            orig_off = nDofs*(indep-1)
            res(i) = matrix(1,iDof) * original(orig_off+1)
            do iOrig=2,nDofs-1
              res(i) = res(i) &
                &      + matrix(iOrig,iDof) * original(orig_off+iOrig)
            end do
            projected(indep + nIndeps*(iDof-1)) = res(i) &
              &                                 + matrix(nDofs,iDof) &
              &                                   * original(orig_Off+nDofs)
          end do
        end do

      end do
      !$OMP END DO

    else

      !$OMP WORKSHARE
      projected = matrix(nDofs,1) * original
      !$OMP END WORKSHARE

    end if

  end subroutine ply_l2_projection
  !****************************************************************************!
 

  !****************************************************************************!
  !> Transformation between modal and nodal values in 1D via L2 projection.
  subroutine ply_l2p_trafo_1D(trafo, projected, original)
    !--------------------------------------------------------------------------!
    !> L2 Projection matrix, this determines the direction of the trafo at hand
    !!
    !! l2p%leg2node = modal to nodal
    !! l2p%node2leg = nodal to modal
    real(kind=rk), intent(in) :: trafo(:,:)

    !> Original coefficients to project.
    real(kind=rk), intent(inout) :: original(:)

    !> Projected coefficients.
    real(kind=rk), intent(inout) :: projected(:)
    !--------------------------------------------------------------------------!
    !--------------------------------------------------------------------------!

    call ply_l2_projection( nIndeps   = 1,         &
      &                     projected = projected, &
      &                     original  = original,  &
      &                     matrix    = trafo      )

  end subroutine ply_l2p_trafo_1D
  !****************************************************************************!


  !****************************************************************************!
  !> Transformation between modal and nodal values in 2D via L2 projection.
  subroutine ply_l2p_trafo_2D(trafo, projected, original)
    !--------------------------------------------------------------------------!
    !> L2 Projection matrix, this determines the direction of the trafo at hand
    !!
    !! l2p%leg2node = modal to nodal
    !! l2p%node2leg = nodal to modal
    real(kind=rk), intent(in) :: trafo(:,:)

    !> Original coefficients to project.
    real(kind=rk), intent(inout) :: original(:)

    !> Projected coefficients.
    real(kind=rk), intent(inout) :: projected(:)
    !--------------------------------------------------------------------------!
    integer :: nDofs
    !--------------------------------------------------------------------------!

    nDofs = size(trafo,1)

    ! Transformation in X direction
    call ply_l2_projection( nIndeps   = nDofs,     &
      &                     projected = projected, &
      &                     original  = original,  &
      &                     matrix    = trafo      )

    ! Transformation in Y direction
    call ply_l2_projection( nIndeps   = nDofs,     &
      &                     projected = original,  &
      &                     original  = projected, &
      &                     matrix    = trafo      )

    ! As we reuse the original array in Y-direction to store the projected
    ! values, thus we need to copy those back into the projected array.

    !$OMP WORKSHARE
    projected = original
    !$OMP END WORKSHARE

  end subroutine ply_l2p_trafo_2D
  !****************************************************************************!


  !****************************************************************************!
  !> Transformation between modal and nodal values in 3D via L2 projection.
  subroutine ply_l2p_trafo_3D(trafo, projected, original)
    !--------------------------------------------------------------------------!
    !> L2 Projection matrix, this determines the direction of the trafo at hand
    !!
    !! l2p%leg2node = modal to nodal
    !! l2p%node2leg = nodal to modal
    real(kind=rk), intent(in) :: trafo(:,:)

    !> Original coefficients to project.
    real(kind=rk), intent(inout) :: original(:)

    !> Projected coefficients.
    real(kind=rk), intent(inout) :: projected(:)
    !--------------------------------------------------------------------------!
    integer :: nDofs
    integer :: nDofs_square
    !--------------------------------------------------------------------------!

    nDofs = size(trafo,1)
    nDofs_square = nDofs**2

    ! Transformation in X direction
    call ply_l2_projection( nIndeps   = nDofs_square, &
      &                     projected = projected,    &
      &                     original  = original,     &
      &                     matrix    = trafo         )

    ! Transformation in Y direction
    call ply_l2_projection( nIndeps   = nDofs_square, &
      &                     projected = original,     &
      &                     original  = projected,    &
      &                     matrix    = trafo         )

    ! Transformation in Z direction
    call ply_l2_projection( nIndeps   = nDofs_square, &
      &                     projected = projected,    &
      &                     original  = original,     &
      &                     matrix    = trafo         )

  end subroutine ply_l2p_trafo_3D
  !****************************************************************************!

end module ply_l2p_module
