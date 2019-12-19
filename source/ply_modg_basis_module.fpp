! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2016, 2018 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014 Verena Krupp <v.krupp@grs-sim.de>
! Copyright (c) 2014, 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Jan Hueckelheim, Melven
! Zoellner and Harald Klimach for German Research School for Simulation
! Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Nikhil Anand, Daniel Petró and Tobias Girresser for University of Siegen.
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
?? include "ply_dof_module.inc"
!> Routines and datatypes related to the modal basis functions of the
!! modal discontinuous Galerkin scheme.
!! \author{Jens Zudrop}
module ply_modg_basis_module
  use env_module,                   only: rk

  use ply_dof_module,               only: Q_space, P_space
  use ply_space_integration_module, only: ply_gaussLegPoints

  implicit none

  private

  !> \brief Coefficients for the projections of the elemental basis functions
  !! from coarser to finer elements and vice versa.
  !!
  !! The MODG scheme is defined with Legendre polynomials (ansatz functions)
  !! and modified Legendre polynomials (test functions). Because of our
  !! dimension-by-dimension approach we can consider 1D elements without loss
  !! of generalty. \n
  !! For MODG scheme the reference element is always \f$ [-1,+1] \f$ . In case
  !! of non-conforming element refinement we have the follwing faces overlying
  !! in the 3D case: \n
  !! \n
  !!      faces of refined                       face of non-refined
  !!          cube                                      cube
  !! ------------------------                 ------------------------
  !! |          |           |                 |                      |
  !! |    3     |     4     |                 |                      |
  !! |          |           |                 |                      |
  !! ------------------------   <---------->  |          5           |
  !! |          |           |                 |                      |
  !! |    1     |     2     |                 |                      |
  !! |          |           |                 |                      |
  !! ------------------------                 ------------------------
  !! \n
  !! To enable flux calculations and projections between the two element sizes
  !! we have to tansfer polynomial functions from one element size to another
  !! one. \n
  !! Therefore we have two tasks: \n
  !! 1. Restrict polynomial functions on 5 to each of the fine element 1 to 4.
  !!    The restriction has to deliver a polynomial approximation on 1 to 4
  !!    in terms of fine element's ansatz functions. \n
  !! 2. Approximate polynomial functions on 1 to 4 by a L2-projection on 5.
  !!    Again the approximation on 5 has to be delivered in terms of ansatz
  !!    functions defined on 5. \n
  !! Without loss of generaltiy we can restrict ourself to the following 1D
  !! situation: \n
  !!                                                             \n
  !!                                                             \n
  !!        Coarse face's ref. element                           \n
  !!  f(x)                                                       \n
  !!   |   -------                                               \n
  !!   | /         \         /                                   \n
  !!   |/           \       /                                    \n
  !!   |             \     /                                     \n
  !!   |              -----                                      \n
  !!   |----------------------|------------->                    \n
  !!  -1                      +1             x                   \n
  !!                                                             \n
  !!             /|\          |                                  \n
  !!              |           |                                  \n
  !!     L2 proj. |           | L2 proj.                         \n
  !!     (approx) |           | (exact)                          \n
  !!              |          \|/                                 \n
  !!                                                             \n
  !!         Fine face's ref. element                            \n
  !!  f(x)                  f(x)                                 \n
  !!   |   -------            |                                  \n
  !!   | /                    |\         /                       \n
  !!   |/                     | \       /                        \n
  !!   |                      |  \     /                         \n
  !!   |                      |   -----                          \n
  !!   |----------|-->        |-----------|-->                   \n
  !!  -1          0  x        0           +1  x                  \n
  !! \n
  !! This datatype stores all the coefficients to calculate the necessary
  !! L2 projections to transfer polynomial functions between coarser and
  !! finer elements (, faces and volumes).
  type ply_modg_refine_type

    !! 1st dim: standard anzatz function [-1,1]
    !! 2nd dim: shifted anzatz function
    !! 3rd dim: shifting for coarse basis function
    !!          1: 1/2x|y|z - 1/2
    !!          2: 1/2x|y|z + 1/2
    real(kind=rk), allocatable :: anz_anzShift(:,:,:)

  end type ply_modg_refine_type

  !> Projection coefficients for covolume filtering.
  type ply_modg_covolume_type

    !! 1st dim: standard anzatz function [-1,1]
    !! 2nd dim: shifted anzatz function
    !! 3rd dim: shifting for coarse basis function
    !!          1: 1/2x|y|z - 1/2
    !!          2: 1/2x|y|z + 1/2
    real(kind=rk), allocatable :: anz_anzShift(:,:,:)

  end type ply_modg_covolume_type


  !> Datatype to represent the polynomial basis functions of the modg scheme.
  type ply_modg_basis_type

    !> Projections of ansatz functions of a finer element to a coarser element
    !! and vice versa. These coefficients are required for non-conforming element
    !! refinement in the MODG scheme.
    type(ply_modg_refine_type) :: refineBaseCoeff

    !> Projections of ansatz functions to covolume grid and vice versa.
    !! These coefficients are required for covolume stabilizations.
    type(ply_modg_covolume_type) :: covolumeBaseCoeff

  end type ply_modg_basis_type


  public :: init_modg_multilevelCoeffs, evalLegendreTensPoly, scalProdLeg, &
    &       scalProdDualLeg, scalProdDualLegDiff,scalProdDualLeg_vec,      &
    &       ply_modg_refine_type, ply_modg_covolume_type,                                        &
    &       faceValLeftBndAns, faceValLeftBndTest,                         &
    &       faceValRightBndTest, ply_modg_basis_type, legendre_1D,         &
    &       faceValLeftBndTestGrad, faceValRightBndTestGrad,               &
    &       faceValLeftBndgradTest, faceValRightBndgradTest,               &
    &       faceValLeftBndTestGrad_vec, faceValRightBndTestGrad_vec,       &
    &       faceValLeftBndgradTest_vec, faceValRightBndgradTest_vec,       &
    &       faceValLeftBndAns_vec, faceValLeftBndTest_vec,                 &
    &       faceValRightBndTest_vec,                                       &
    &       faceValLeftBndDiffAns, faceValRightBndDiffAns,                 &
    &       init_modg_covolumeCoeffs, integrateLeg

contains

  ! ************************************************************************ !
  !> Integral of combination of all anzatz functions for
  !! projection onto finer element
  subroutine init_modg_covolumeCoeffs( nPoints, nFunc, integral )
    ! -------------------------------------------------------------------- !
    ! The number of quadrature points to be used
    integer, intent(in) :: nPoints
    ! number of anzatz and test functions
    integer, intent(in) :: nFunc
    ! integration results
    type(ply_modg_covolume_type), intent(out) :: integral
    ! -------------------------------------------------------------------- !
    ! Gaussian points array
    real(kind=rk), allocatable  :: GaussPoints(:)
    real(kind=rk), allocatable  :: GaussPoints_left(:)
    real(kind=rk), allocatable  :: GaussPoints_right(:)
    !! points and weights for gauss-legendre quadrature
    real(kind=rk) :: tempLeft(nPoints), tempRight(nPoints)
    real(kind=rk) :: sumLeft, sumRight
    !! Gaussian weights
    real(kind=rk), allocatable  :: w(:)
    !! legendre polynomila values left on [-1;0]
    real(kind=rk) :: legendre_left(nFunc, nPoints)
    real(kind=rk) :: legendre_left_shifted(nFunc, nPoints)
    !! legendre polynomila values right on [0;+1]
    real(kind=rk) :: legendre_right(nFunc, nPoints)
    real(kind=rk) :: legendre_right_shifted(nFunc, nPoints)
    integer :: iFunc, jFunc
    ! -------------------------------------------------------------------- !

    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(jFunc, iFunc)

    allocate(GaussPoints(nPoints))
    allocate(GaussPoints_left(nPoints))
    allocate(GaussPoints_right(nPoints))
    allocate(w(nPoints))

    ! get GL points and weights on reference element [-1;+1]
    call ply_gaussLegPoints(x1    = -1.0_rk,     &
      &         x2    = 1.0_rk,                  &
      &         x     = GaussPoints,             &
      &         w     = w,                       &
      &         nIntP = nPoints                  )

    ! shift the gauss points to
    ! ... the left integral domain, i.e. [-1;0]
    GaussPoints_left = ( (0.0_rk + 1.0_rk ) / 2.0_rk) * GaussPoints &
      & + ( ( 0.0_rk - 1.0_rk ) / ( 2.0_rk ) )
    ! ... the right integral domain, i.e. [0;+1]
    GaussPoints_right = ( ( 1.0_rk - 0.0_rk ) / ( 2.0_rk ) ) * GaussPoints &
      & + ( ( 1.0_rk + 0.0_rk ) / ( 2.0_rk ) )

    ! Scale the weights for integration over domains of length 1.0
    w(:) = w(:) * 0.5_rk

    deallocate( GaussPoints )

    ! Calculate values of Legendre polynomials
    ! ... on [-1;0]
    legendre_left = legendre_1D(GaussPoints_left, nFunc-1)
    ! ... on [0;+1]
    legendre_right = legendre_1D(GaussPoints_right, nFunc-1)

    ! Calculate values of shifted Legendre polynomials
    ! ... for [-1;0]
    legendre_left_shifted = legendre_1D(GaussPoints_left+1.0_rk, nFunc-1)
    ! ... for [0;+1]
    legendre_right_shifted = legendre_1D(GaussPoints_right-1.0_rk, nFunc-1)

    allocate( integral%anz_anzShift(1:nFunc, 1:nFunc, 2))

    !$OMP DO
    !loop over anzatz functions
    do jFunc = 1, nFunc
      do iFunc = 1, nFunc

        ! ansatz-ansatz with left shift integral (on [-1;0])
        tempLeft = legendre_left(iFunc, :)             &
          &          * legendre_left_shifted(jFunc, :) &
          &          * w(:)
        sumLeft = sum(tempLeft)
        integral%anz_anzShift(iFunc, jFunc, 1) = sumLeft / scalProdLeg(iFunc)


        ! ansatz-ansatz with right shift integral (on [0;+1])
        tempRight = legendre_right(iFunc, :)             &
          &           * legendre_right_shifted(jFunc, :) &
          &           * w(:)
        sumRight = sum(tempRight)
        integral%anz_anzShift(iFunc, jFunc, 2) = sumRight / scalProdLeg(iFunc)

      end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine init_modg_covolumeCoeffs
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Integral of combination of all anzatz functions for
  !! projection onto finer element
  subroutine init_modg_multilevelCoeffs( nPoints, nFunc, integral )
    ! -------------------------------------------------------------------- !
    ! integration results
    type(ply_modg_refine_type), intent(out) :: integral
    ! number of anzatz and test functions
    integer, intent(in) :: nFunc
    integer, intent(in) :: nPoints
    ! -------------------------------------------------------------------- !
    ! Gaussian points array
    real(kind=rk), allocatable  :: GaussPoints(:)
    !! points and weights for gauss-legendre quadrature
    real(kind=rk) :: tempLeft(nPoints), tempRight(nPoints)
    real(kind=rk) :: sumLeft, sumRight
    !! Gaussian weights
    real(kind=rk), allocatable :: w(:)
    !! legendre polynomial values [-1,1]
    real(kind=rk) :: legendre_standard(nFunc, nPoints)
    !! legendre polynomila values left shift
    real(kind=rk) :: legendre_left(nFunc, nPoints)
    !! legendre polynomila values right shift
    real(kind=rk) :: legendre_right(nFunc, nPoints)
    integer :: iFunc, jFunc
    ! -------------------------------------------------------------------- !

    allocate(GaussPoints(nPoints))
    allocate(w(nPoints))

    ! get GL points and weights
    call ply_gaussLegPoints( x1    = -1.0_rk,     &
      &                      x2    = 1.0_rk,      &
      &                      x     = GaussPoints, &
      &                      w     = w,           &
      &                      nIntP = nPoints      )

    ! Calculate values of legendre polynomials
    legendre_standard = legendre_1D(GaussPoints, nFunc-1)
    legendre_left = legendre_1D(GaussPoints/2.0_rk-1.0_rk/2.0_rk, nFunc-1)
    legendre_right = legendre_1D(GaussPoints/2.0_rk+1.0_rk/2.0_rk, nFunc-1)
    allocate( integral%anz_anzShift(1:nFunc, 1:nFunc, 2))

    !loop over anzatz functions
    do jFunc = 1, nFunc
      do iFunc = 1, nFunc
        !left shift

        ! anzatz-anzatz with left shift integral
        tempLeft = legendre_standard(iFunc, :) * legendre_left(jFunc, :) * w(:)
        sumLeft = sum(tempLeft)
        integral%anz_anzShift(iFunc, jFunc, 1) = sumLeft

        ! anzatz-anzatz with right shift integral
        tempRight = legendre_standard(iFunc, :)  &
          &           * legendre_right(jFunc, :) &
          &           * w(:)
        sumRight = sum(tempRight)
        integral%anz_anzShift(iFunc, jFunc, 2) = sumRight

      end do
    end do


  end subroutine init_modg_multilevelCoeffs
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Integrate the integrand function in Legendre basis, and represent the
  !! integral again in the Legendre basis up to the maximal degree.
  !!
  !! The maximal degree needs to be one higher than the maximal degree of the
  !! integrand in order to fully represent the integral, otherwise the result
  !! is truncated and only an approximation up to maxdegree is obtained.
  !!
  !! The integral needs to have a length of maxdegree+1.
  !! maxdegree needs to be non-negative. Thus, integral needs to be an array
  !! with a length of at least 1!
  !!
  !! Implemented property of Legendre Polynomials:
  !! L_i(x) = 1/(2*i + 1) * d/dx [ L_{i+1}(x) - L_{i-1}(x) ]
  pure function integrateLeg( integrand, maxdegree ) result(integral)
    ! -------------------------------------------------------------------- !
    !> Coefficients of the function to integrate in Legendre basis.
    real(kind=rk), intent(in) :: integrand(:)

    !> Maximal polynomial degree for the integral, should be larger than the
    !! degree of the integrand
    integer, intent(in) :: maxdegree

    !> Legendre coefficients of the resulting integral.
    real(kind=rk) :: integral(maxdegree+1)
    ! -------------------------------------------------------------------- !
    integer :: nOrigModes
    integer :: minModes
    integer :: iMode
    ! -------------------------------------------------------------------- !

    nOrigModes = size(integrand)
    minModes = min(nOrigModes-1, maxdegree+1)

    if (nOrigModes >= 2) then
      integral(1) = -1.0_rk/3.0_rk * integrand(2)
    else
      integral(1) = 0.0_rk
    end if

    do iMode=2,minModes
      integral(iMode) = integrand(iMode-1)/real(2*iMode-3, kind=rk) &
        &               - integrand(iMode+1)/real(2*iMode+1, kind=rk)
    end do

    do iMode=max(nOrigModes,2),min(maxDegree+1,nOrigModes+1)
      integral(iMode) = integrand(iMode-1)/real(2*iMode-3, kind=rk)
    end do

    do iMode=nOrigModes+2,maxDegree+1
      integral(iMode) = 0.0_rk
    end do

  end function integrateLeg
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Evaluate all 1D Legendre polynomials at a given set
  !! of points up to the given degree.
  pure function legendre_1D( points, degree ) result(one_dim_eval)
    ! -------------------------------------------------------------------- !
    !> 1D points to evaluate.
    real(kind=rk), intent(in) :: points(:)
    !> Degree up to which to evaluate the polynomials
    integer,intent(in) ::degree
    !> Resulting vector of all mode values at all points
    real(kind=rk) :: one_dim_eval(degree+1, size(points))
    ! -------------------------------------------------------------------- !
    integer :: iDegree
    real(kind=rk) :: n_q
    ! -------------------------------------------------------------------- !

    !> init the first two Legendre polynomials.
    !! ... the first Legendre polynomial is 1
    one_dim_eval(1, :) = 1

    if (degree > 0) then
      !! ... the second Legendre polynomial is x
      one_dim_eval(2, :) = points(:)

      do iDegree = 2, degree
        n_q = 1.0_rk / real(iDegree , kind=rk)
        !> Recursive polynomial evaluation:
        !! \f$ n L_{n}(x)= (2n - 1) x L_{n-1}(x) - (n-1)L_{n-2}(x) \f$
        one_dim_eval(iDegree + 1,:)                  &
          &  = n_q * ( ( 2 * iDegree - 1 )           &
          &              * points(:)                 &
          &              * one_dim_eval(iDegree,:)   &
          &            - ( iDegree - 1 )             &
          &              * one_dim_eval(iDegree-1,:) &
          &          )
      end do
    end if
  end function legendre_1D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Evaluate three-dimensional tensor product Legendre polynomials
  !! (not-normalized) at a given set of coordinates.
  subroutine evalLegendreTensPoly( coords , nCoords, maxPolyDegree, basisType, &
    &                              polyVal )
    ! -------------------------------------------------------------------- !
    !> Array of coordinates (on the reference element) to evaluate the tensor
    !! product polynomials at. First dimension is nCoord, second is 3 for x,y,z
    !! component.
    real(kind=rk), intent(in) :: coords(:,:)

    !> The number of coordinates to evaluate the polynomials at.
    integer, intent(in) :: nCoords

    !> The maximum polynomail degree of the MODG scheme.
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: basisType

    !> The polynomial values. First dimension is the number of tensor product
    !! polynomials and the second dimension is the number of points, i.e.
    !! nCoords.
    real(kind=rk), allocatable, intent(out) :: polyVal(:,:)
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: polyValX(:,:), polyValY(:,:), polyValZ(:,:)
    integer :: iAnsX, iAnsY, iAnsZ, iAns, ansPos, ansPosMax
    real(kind=rk) :: n_q
    ! -------------------------------------------------------------------- !

    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iAns, iAnsX, iAnsY, iAnsZ, n_q, ansPos)

    ! allocate the output array
    select case(basisType)
      case(Q_space)
        allocate( polyVal( (maxPolyDegree+1)**3 ,nCoords) )
      case(P_space)
        allocate( polyVal((maxPolydegree+1) * (maxPolydegree+2) &
          &         * ( maxPolydegree+3) / 6, &
          &       nCoords ) )
    end select

    allocate( polyValX( (maxPolyDegree+1) ,nCoords) )
    allocate( polyValY( (maxPolyDegree+1) ,nCoords) )
    allocate( polyValZ( (maxPolyDegree+1) ,nCoords) )

    ! Evaluate the Legendre polynomials per direction:
    ! ... first Legendere polynmoial is constant
    polyValX(1,:) = 1.0_rk
    polyValY(1,:) = 1.0_rk
    polyValZ(1,:) = 1.0_rk

    if(maxPolyDegree > 0) then
      ! ... second Legendere polynmoial is identity
      polyValX(2,:) = coords(:,1)
      polyValY(2,:) = coords(:,2)
      polyValZ(2,:) = coords(:,3)
      !$OMP DO
      ! ... higher order polynomials are build recursively
      do iAns = 3, maxPolyDegree+1
        n_q = 1.0_rk / real(iAns-1,kind=rk)
        ! x recursion
        polyValX(iAns,:) = ( ( 2 * ( iAns - 1 ) - 1 ) &
          &                    * coords(:,1)          &
          &                    * polyValX(iAns-1,:)   &
          &                  - ( ( iAns - 1 ) - 1 )   &
          &                    * polyValX(iAns-2,:) ) &
          &                  * n_q
        ! y recursion
        polyValY(iAns,:) = ( ( 2 * ( iAns - 1 ) - 1 ) &
          &                    * coords(:,2)          &
          &                    * polyValY(iAns-1,:)   &
          &                  - ( ( iAns - 1 ) - 1 )   &
          &                    * polyValY(iAns-2,:) ) &
          &                  * n_q
        ! z recursion
        polyValZ(iAns,:) = ( ( 2 * ( iAns - 1 ) - 1 ) &
          &                    * coords(:,3)          &
          &                    * polyValZ(iAns-1,:)   &
          &                 - ( ( iAns - 1 ) - 1 )    &
          &                   * polyValZ(iAns-2,:) )  &
          &                 *n_q
      end do
      !$OMP END DO
    end if

    ! Now, build the complete point value.
    select case(basisType)
      case(Q_space)
        !$OMP DO
        do iAnsX = 1, maxPolyDegree+1
          do iAnsY = 1, maxPolyDegree+1
            do iAnsZ = 1, maxPolyDegree+1
              ! get the position of this ansatz function combination.
?? copy :: posOfModgCoeffQTens(iAnsX, iAnsY, iAnsZ, maxPolyDegree, ansPos)
              polyVal(ansPos, :) = polyValX(iAnsX,:)     &
                &                    * polyValY(iAnsY,:) &
                &                    * polyValZ(iAnsZ,:)
            end do
          end do
        end do
        !$OMP END DO
      case(P_space)
        iAnsX = 1
        iAnsY = 1
        iAnsZ = 1
?? copy :: getDofsPTens(maxPolyDegree, ansPosMax)
        !$OMP DO
        do ansPos = 1, ansPosMax
          polyVal(ansPos, :) = polyValX(iAnsX,:)     &
            &                    * polyValY(iAnsY,:) &
            &                    * polyValZ(iAnsZ,:)
?? copy :: nextModgCoeffPTens(iAnsX, iAnsY, iAnsZ)
        end do
        !$OMP END DO
    end select

    !$OMP END PARALLEL

  end subroutine evalLegendreTensPoly
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Returns the value of the non-normalized differentiated Legendre polynomial
  !! at the right boundary of the reference element, i.e. at +1.
  pure function faceValRightBndDiffAns( ansFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first ansatz function has index 1.
    integer, intent(in) :: ansFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    val = ansFunc * ( ansFunc + 1 ) * 0.5_rk

  end function faceValRightBndDiffAns
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Returns the value of the non-normalized Legendre polynomial at the left
  !! boundary of the reference element, i.e. at -1.
  pure function faceValLeftBndAns( ansFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first ansatz function has index 1.
    integer, intent(in) :: ansFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    val = ( -1.0_rk )**( ansFunc - 1 )

  end function faceValLeftBndAns
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the non-normalized Legendre polynomial at the left
  !! boundary of the reference element, i.e. at -1.
  function faceValLeftBndAns_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first ansatz function has index 1.
    integer, intent(in) :: mPD
    integer :: ansFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))
    
    do ansFunc=1, mPD
      val(ansFunc) = ( -1.0_rk )**( ansFunc - 1 )
    end do
  end function faceValLeftBndAns_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the non-normalized differentiated Legendre polynomial
  !! at the leftboundary of the reference element, i.e. at -1.
  pure function faceValLeftBndDiffAns( ansFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first ansatz function has index 1.
    integer, intent(in) :: ansFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if (ansFunc ==1) then
      val = 0.0_rk
    else
      val = ansFunc * ( ansFunc + 1 ) * 0.5_rk * ( ( -1.0_rk )**( ansFunc ) )
    endif

  end function faceValLeftBndDiffAns
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Returns the value of the dual Legendre polynomial at the right
  !! boundary of the reference element, i.e. at +1.
  pure function faceValRightBndTest( testFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: testFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if(testFunc < 3) then
      val = 1.0_rk
    else
      val = 0.0_rk
    end if

  end function faceValRightBndTest
  ! ************************************************************************ !

! ************************************************************************ !
  !> Returns the value of the dual Legendre polynomial at the right
  !! boundary of the reference element, i.e. at +1. Vectorized Version.
  function faceValRightBndTest_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: mPD
    integer :: testFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))
    
    do testFunc = 1, mPD
      if(testFunc < 3) then
        val(testFunc) = 1.0_rk
      else
        val(testFunc) = 0.0_rk
      end if
    end do
  end function faceValRightBndTest_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the gradient of dual Legendre polynomial at the right
  !! boundary of the reference element, i.e. at +1.
  pure function faceValRightBndgradTest( testFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: testFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if(testFunc == 1) then
      val = 0.0_rk
    else
      val = ( testFunc - 2 ) * 2 + 1.0_rk
    end if

  end function faceValRightBndgradTest
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the gradient of dual Legendre polynomial at the right
  !! boundary of the reference element, i.e. at +1. Vectorized version.
  function faceValRightBndgradTest_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: mPD
    integer :: testFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))
    
    do testFunc = 1, mPD
      if(testFunc == 1) then
        val(testFunc) = 0.0_rk
      else
        val(testFunc) = ( testFunc - 2 ) * 2 + 1.0_rk
      end if
    end do
  end function faceValRightBndgradTest_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the dual Legendre polynomial at the left
  !! boundary of the reference element, i.e. at -1.
  pure function faceValLeftBndTest( testFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: testFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if(testFunc == 1) then
      val = 1.0_rk
    elseif(testFunc == 2) then
      val = -1.0_rk
    else
      val = 0.0_rk
    end if

  end function faceValLeftBndTest
  ! ************************************************************************ !

! ************************************************************************ !
  !> Returns the value of the dual Legendre polynomial at the left
  !! boundary of the reference element, i.e. at -1.Vectorized version.
  function faceValLeftBndTest_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: mPD
    integer :: testFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))

    do testFunc = 1, mPD
      if(testFunc == 1) then
        val(testFunc) = 1.0_rk
      elseif(testFunc == 2) then
        val(testFunc) = -1.0_rk
      else
        val(testFunc) = 0.0_rk
      end if
    end do
  end function faceValLeftBndTest_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the gradient of the dual Legendre polynomial at the
  !! left boundary of the reference element, i.e. at -1.
  pure function faceValLeftBndgradTest( testFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: testFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if(testFunc == 1) then
      val = 0.0_rk
    else
      val = ( 2.0_rk * ( testFunc - 2 ) + 1.0_rk ) * ( -1 )**testFunc
    end if

  end function faceValLeftBndgradTest
  ! ************************************************************************ !

 ! ************************************************************************ !
  !> Returns the value of the gradient of the dual Legendre polynomial at the
  !! left boundary of the reference element, i.e. at -1.
  function faceValLeftBndgradTest_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: mPD
    integer :: testFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))

    do testFunc = 1, mPD
      if(testFunc == 1) then
        val(testFunc) = 0.0_rk
      else
        val(testFunc) = ( 2.0_rk * ( testFunc - 2 ) + 1.0_rk ) * ( -1 )**testFunc
      end if
    end do
    
  end function faceValLeftBndgradTest_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the derivaitve of the dual Legendre polynomial at the
  !! left boundary of the reference element, i.e. at -1.
  pure function faceValLeftBndTestGrad( testFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: testFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if(testFunc==1) then
      val = 0.0_rk
    else
      val = (-1.0_rk)**(testFunc)
    end if

  end function faceValLeftBndTestGrad
  ! ************************************************************************ !

! ************************************************************************ !
  !> Returns the value of the derivaitve of the dual Legendre polynomial at the
  !! left boundary of the reference element, i.e. at -1.Vectorized version.
  function faceValLeftBndTestGrad_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: mPD
    integer :: testFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))

    do testFunc = 1, mPD
      if(testFunc==1) then
        val(testFunc) = 0.0_rk
      else
        val(testFunc) = (-1.0_rk)**(testFunc)
      end if
    end do 

  end function faceValLeftBndTestGrad_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Returns the value of the derivaitve of the dual Legendre polynomial at the right
  !! boundary of the reference element, i.e. at +1.
  pure function faceValRightBndTestGrad( testFunc ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: testFunc
    !> The function value.
    real(kind=rk) :: val
    ! -------------------------------------------------------------------- !

    if(testFunc==1) then
      val = 0.0_rk
    else
      val = 1.0_rk
    end if

  end function faceValRightBndTestGrad
  ! ************************************************************************ !

! ************************************************************************ !
  !> Returns the value of the derivaitve of the dual Legendre polynomial at the right
  !! boundary of the reference element, i.e. at +1.vectoized version.
  function faceValRightBndTestGrad_vec( mPD ) result(val)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, first test function has index 1.
    integer, intent(in) :: mPD
    integer :: testFunc
    !> The function value.
    real(kind=rk),allocatable :: val(:)
    ! -------------------------------------------------------------------- !
    allocate(val(mPD))

    do testFunc = 1, mPD
      if(testFunc==1) then
        val(testFunc) = 0.0_rk
      else
        val(testFunc) = 1.0_rk
      end if
    end do
    
  end function faceValRightBndTestGrad_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Function to calculate the L2 scalar product of a Legendre polynomial
  !! with itself on the reference element [-1,+1].
  pure function scalProdLeg( ansFunc ) result(scalProd)
    ! -------------------------------------------------------------------- !
    !> The Legendre polynomial to calculate the scalar product for.
    !! The first Legendre polynomial has index 1.
    integer, intent(in) :: ansFunc
    !> The scalar product on the refenece element [-1,+1].
    real(kind=rk) :: scalProd
    ! -------------------------------------------------------------------- !

    scalProd = 2.0 / ( 2.0_rk * ansFunc - 1.0_rk )

  end function scalProdLeg
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function to calculate the scalar product between a Legendre polynomial
  !! (ansatz function) and a dual Legendre polynomial (test function) on the
  !! reference element [-1;+1].
  pure function scalProdDualLeg( ansFunc, testFunc ) result(scalProd)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, there first ansatz function has index 1.
    integer, intent(in) :: ansFunc
    !> The test function index, there first test function has index 1.
    integer, intent(in) :: testFunc
    !> The scalar product of the two functions.
    real(kind=rk) :: scalProd
    ! -------------------------------------------------------------------- !

    if( ansFunc == testFunc ) then
      scalProd = 2.0_rk / ( 2.0_rk * testFunc - 1.0_rk )
    elseif( ansFunc == testFunc - 2 ) then
      scalProd = ( -2.0_rk ) / ( 2.0_rk * ansFunc - 1.0_rk )
    else
      scalProd = 0.0_rk
    end if

  end function scalProdDualLeg
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Function to calculate the scalar product between a Legendre polynomial
  !! (ansatz function) and a differentiated dual Legendre polynomial (test
  !! function) on the reference element [-1;+1].
  pure function scalProdDualLegDiff( ansFunc, testFunc ) result(scalProd)
    ! -------------------------------------------------------------------- !
    !> The ansatz function index, there first ansatz function has index 1.
    integer, intent(in) :: ansFunc
    !> The test function index, there first test function has index 1.
    integer, intent(in) :: testFunc
    !> The scalar product of the two functions.
    real(kind=rk) :: scalProd
    ! -------------------------------------------------------------------- !

    if(ansFunc == testFunc-1) then
      scalProd = 2.0_rk
    else
      scalProd = 0.0_rk
    end if

  end function scalProdDualLegDiff
  ! ************************************************************************ !
  
  ! ************************************************************************ !
  !> Vectorized Function to calculate the scalar product between a Legendre polynomial
  !! (ansatz function) and a dual Legendre polynomial (test function) on the
  !! reference element [-1;+1] and
  !! to calculate the scalar product between a Legendre polynomial
  !! (ansatz function) and a differentiated dual Legendre polynomial (test
  !! function) on the reference element [-1;+1].
  function scalProdDualLeg_vec( ansFunc, testFunc, mPd) result(scalProd)
  ! -------------------------------------------------------------------------!
    !> maxPolyDegree
    integer, intent(in) :: mPd
    !> The ansatz function index, there first ansatz function has index 1.
    integer, intent(in) :: ansFunc
    !> The test function index, there first test function has index 1.
    integer, intent(in) :: testFunc
    !> The scalar product of the two functions.
    real(kind=rk) :: scalProd(mPd+1)
    integer :: iTest

    do iTest = 1, mPd+1
     if( ansFunc == testFunc ) then
       scalProd(iTest) = 2.0_rk / ( 2.0_rk * iTest - 1.0_rk )
     elseif( ansFunc == testFunc - 2 ) then
       scalProd(iTest) = ( -2.0_rk ) / ( 2.0_rk * (iTest-2) - 1.0_rk )
     elseif(ansFunc == testFunc-1) then
       scalProd(iTest) = 2.0_rk
     else
       scalProd(iTest) = 0.0_rk
     end if
   end do
   
  end function scalProdDualLeg_vec
  
end module ply_modg_basis_module
