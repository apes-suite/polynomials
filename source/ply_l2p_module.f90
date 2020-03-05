! Copyright (c) 2012-2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2015,2018,2020 Harald Klimach <harald.klimach@uni-siegen.de.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2012 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014,2016 Verena Krupp
! Copyright (c) 2013-2014,2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Jan Hueckelheim, Melven
! Zoellner and Harald Klimach for German Research School for Simulation
! Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Tobias Girresser, Jiaxing Qi and Nikhil Anand for University of Siegen.
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !

module ply_l2p_module
  use env_module, only: rk

  use tem_compileconf_module, only: vlen

  use ply_modg_basis_module, only: scalProdLeg
  use ply_space_integration_module, only: ply_gaussLegPoints
  use ply_modg_basis_module, only: legendre_1D

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


  ! ************************************************************************ !
  subroutine Copy_ply_l2p( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    !type(ply_legFpt_2D_type), intent(out) :: left
    type(ply_l2p_type), intent(out) :: left
    !> fpt to copy from
    !type(ply_legFpt_2D_type), intent(in) :: right
    type(ply_l2p_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%leg2node = right%leg2node
    left%node2leg = right%node2leg

  end subroutine Copy_ply_l2p
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Initialize the transformations via L2 projections.
  subroutine ply_init_l2p( l2p, degree )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_type), intent(out) :: l2p
    integer, intent(in) :: degree
    ! -------------------------------------------------------------------- !
    integer :: iPoint, iDof
    integer :: nDofs
    integer :: nPoints
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: leg1D_at_gauss(:,:)
    real(kind=rk), allocatable :: weights1D(:)
    real(kind=rk) :: scalprod_q
    ! -------------------------------------------------------------------- !

    nDofs = degree+1
    nPoints = nDofs

    allocate(leg1D_at_gauss(max(2,nDofs), nPoints))
    allocate(gaussp1D(nPoints))
    allocate(weights1D(nPoints))

    ! create the quadrature points for volume and on the face
    ! for the oversampled projection
    call ply_gaussLegPoints( x1    = -1.0_rk,   &
      &                      x2    = 1.0_rk,    &
      &                      nIntP = nPoints,   &
      &                      w     = weights1D, &
      &                      x     = gaussp1D   )
    leg1D_at_gauss = legendre_1D( points = gaussp1D, &
      &                           degree = degree    )

    allocate(l2p%leg2node(nDofs, nPoints))
    allocate(l2p%node2leg(nPoints, nDofs))

    ! Coefficients to transform legendre to nodal values.
    l2p%leg2node = leg1D_at_gauss(:nDofs, :)

    ! Coefficients to transform nodal to legendre values.
    do iDoF = 1,nDofs
      scalProd_q = 1.0_rk / scalProdLeg(iDoF)
      do iPoint = 1,nPoints
        l2p%node2leg(iPoint, iDoF) = l2p%leg2node(iDoF, iPoint) &
          &                            * weights1D(iPoint)        &
          &                            * scalProd_q
      end do
    end do

  end subroutine ply_init_l2p
  ! ************************************************************************ !


  ! ************************************************************************ !
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
  subroutine ply_l2_projection( nDofs, nIndeps, projected, original, matrix )
    !ICE! Directive for Cray compiler to prevent inlining of this routine,
    !ICE! what causes the compiler to fail.
    !dir$ inlinenever ply_l2_projection
    ! -------------------------------------------------------------------- !
    !> Number of degree of freedoms
    integer, intent(in) :: nDofs

    !> Number of values that can be computed independently.
    integer, intent(in) :: nIndeps

    !> Projected data.
    !!
    !! Size has to be nIndeps*size(matrix,1), and the layout is changed here
    !! when compared to the original array, as the projected direction moves
    !! to the end.
    real(kind=rk), intent(out) :: projected(nIndeps, nDofs)

    !> Original data.
    !!
    !! Size has to be size(matrix,1) and the direction to be projected has to
    !! be the fastest running one.
    real(kind=rk), intent(in) :: original(nDofs, nIndeps)

    !> Matrix to apply in this operation.
    !!
    !! The matrix defines wether this is a modal to nodal transformation or the
    !! other way around.
    real(kind=rk), intent(in) :: matrix(nDofs,nDofs)
    ! -------------------------------------------------------------------- !
    integer :: iRow, iCol, iCell, iStrip, strip_ub
    real(kind=rk) :: mval
    ! JQ: on SX-ACE, vlen=nIndeps gives the best performance
    !     on    x86, vlen=256     gives the best performance
    ! integer, parameter :: vlen = nIndeps
    ! -------------------------------------------------------------------- !

    !$OMP PARALLEL DEFAULT(SHARED), &
    !$OMP PRIVATE(iStrip, iRow, iCell, iCol, strip_ub, mval)
    if (nDofs > 1) then

      !$OMP DO
      do iStrip=0,nIndeps-1,vlen

        ! Calculate the upper bound of the current strip
        strip_ub = min(iStrip + vlen, nIndeps) - iStrip

        do iRow = 1, nDofs

          do iCell = iStrip+1, iStrip+strip_ub
            projected(iCell, iRow) = 0.0_rk
          end do
          do iCol = 1, nDofs
            mval =  matrix(iCol,iRow)
            do iCell = iStrip+1, iStrip+strip_ub
              ! on SX-ACE, this can be identified as matrix multiplication
              ! which results in VERY HIGH performance
              projected(iCell, iRow) = projected(iCell, iRow) &
                &                   + mval * original(iCol, iCell)
            end do ! iCell
          end do ! iCol = 1, nCols

        end do ! iRow = 1, nRows
      end do ! iStrip
      !$OMP END DO

    else

      !$OMP SINGLE
      projected = matrix(nDofs,1) * original
      !$OMP END SINGLE

    end if
    !$OMP END PARALLEL

  end subroutine ply_l2_projection
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transformation between modal and nodal values in 1D via L2 projection.
  subroutine ply_l2p_trafo_1D( trafo, projected, original )
    ! -------------------------------------------------------------------- !
    !> L2 Projection matrix, this determines the direction of the trafo at hand
    !!
    !! l2p%leg2node = modal to nodal
    !! l2p%node2leg = nodal to modal
    real(kind=rk), intent(in) :: trafo(:,:)

    !> Original coefficients to project.
    real(kind=rk), intent(inout) :: original(:)

    !> Projected coefficients.
    real(kind=rk), intent(inout) :: projected(:)
    ! -------------------------------------------------------------------- !

    call ply_l2_projection( nIndeps   = 1,             &
      &                     nDofs     = size(trafo,1), &
      &                     projected = projected,     &
      &                     original  = original,      &
      &                     matrix    = trafo          )

  end subroutine ply_l2p_trafo_1D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transformation between modal and nodal values in 2D via L2 projection.
  subroutine ply_l2p_trafo_2D( trafo, projected, original )
    ! -------------------------------------------------------------------- !
    !> L2 Projection matrix, this determines the direction of the trafo at hand
    !!
    !! l2p%leg2node = modal to nodal
    !! l2p%node2leg = nodal to modal
    real(kind=rk), intent(in) :: trafo(:,:)

    !> Original coefficients to project.
    real(kind=rk), intent(inout) :: original(:)

    !> Projected coefficients.
    real(kind=rk), intent(inout) :: projected(:)
    ! -------------------------------------------------------------------- !
    integer :: nDofs
    ! -------------------------------------------------------------------- !

    nDofs = size(trafo,1)

    ! Transformation in X direction
    call ply_l2_projection( nIndeps   = nDofs,     &
      &                     nDofs     = nDofs,     &
      &                     projected = projected, &
      &                     original  = original,  &
      &                     matrix    = trafo      )

    ! Transformation in Y direction
    call ply_l2_projection( nIndeps   = nDofs,     &
      &                     nDofs     = nDofs,     &
      &                     projected = original,  &
      &                     original  = projected, &
      &                     matrix    = trafo      )

    ! As we reuse the original array in Y-direction to store the projected
    ! values, thus we need to copy those back into the projected array.

    projected = original

  end subroutine ply_l2p_trafo_2D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transformation between modal and nodal values in 3D via L2 projection.
  subroutine ply_l2p_trafo_3D( trafo, projected, original )
    ! -------------------------------------------------------------------- !
    !> L2 Projection matrix, this determines the direction of the trafo at hand
    !!
    !! l2p%leg2node = modal to nodal
    !! l2p%node2leg = nodal to modal
    real(kind=rk), intent(in) :: trafo(:,:)

    !> Original coefficients to project.
    real(kind=rk), intent(inout) :: original(:)

    !> Projected coefficients.
    real(kind=rk), intent(inout) :: projected(:)
    ! -------------------------------------------------------------------- !
    integer :: nDofs
    integer :: nDofs_square
    ! -------------------------------------------------------------------- !

    nDofs = size(trafo,1)
    nDofs_square = nDofs**2

    ! Transformation in X direction
    call ply_l2_projection( nIndeps   = nDofs_square, &
      &                     nDofs     = nDofs,        &
      &                     projected = projected,    &
      &                     original  = original,     &
      &                     matrix    = trafo         )

    ! Transformation in Y direction
    call ply_l2_projection( nIndeps   = nDofs_square, &
      &                     nDofs     = nDofs,        &
      &                     projected = original,     &
      &                     original  = projected,    &
      &                     matrix    = trafo         )

    ! Transformation in Z direction
    call ply_l2_projection( nIndeps   = nDofs_square, &
      &                     nDofs     = nDofs,        &
      &                     projected = projected,    &
      &                     original  = original,     &
      &                     matrix    = trafo         )

  end subroutine ply_l2p_trafo_3D
  ! ************************************************************************ !

end module ply_l2p_module
