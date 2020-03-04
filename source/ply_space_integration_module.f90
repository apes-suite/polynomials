! Copyright (c) 2011-2012 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014,2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012-2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2012 Vyacheslav Korchagin <v.korchagin@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp <v.krupp@grs-sim.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Vyacheslav Korchagin, Melven
! Zoellner and Harald Klimach for German Research School for Simulation
! Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Daniel Petró and Nikhil Anand for University of Siegen.
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

!> Spatial integration with the Gauss-Legendre numerical integration.
module ply_space_integration_module
  use env_module,         only: rk
  use tem_param_module,   only: PI
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  implicit none

  private

  !> Type to specify the space integration of our variational method.
  type ply_space_quadrature_type
    !> Unique characters specifying the space integration type.
    !! If you add an additional space integration method, then you should
    !! add a unique abbreviation for that here. The following relations
    !! hold:
    !! int_type == 'GALEINT' means Gauss-Legendre interpolation
    character(len=32) :: int_type

    !> The number of space integration points we need to evaluate a volume
    !! integral correctly. For example for
    !! Gauss-Legendre quadrature it is (\frac{2 * pm + 1}{2})^3.
    integer :: numVolumePoints = 0

    !> The number of space integration points we need to evaluate a surface
    !! integral correctly. For example for
    !! Gauss-Legendre quadrature it is (\frac{2 * pm + 1}{2})^2.
    integer :: numSurfacePoints = 0

    !> The integration points to evaluate a volume integral over the
    !! reference cubic element.
    !! The first dimension of this array is numVolumePoints.
    !! The second dimension is 3.
    real(kind=rk), allocatable :: integration_points_volume(:, :)

    !> Weights used for quadrature formula of volume integration.
    !! The dimensions is numVolumePoints.
    real(kind=rk), allocatable :: integration_weights_volume(:)

    !> The integration points to evaluate a surface integral over the
    !! reference element.
    !! The first dimension of this array is numSufacePoints.
    !! The second dimens is 2
    real(kind=rk), allocatable :: integration_points_surface(:,:)

    !> Weights used for quadrature formula of surface integration.
    !! The dimensions is numSurfacePoints.
    real(kind=rk), allocatable :: integration_weights_surface(:)

    !> This vector holds all gauss point coordinates for a one dimensional
    !! quadrature. Additionally the start and end point of the unit interval
    !! (i.e. 0.0 and 1.0) are appended.
    real(kind=rk), allocatable :: integration_points_line(:)

  end type ply_space_quadrature_type

  !> interface to create gauss points for pure gauss quadrature on the
  !! reference element. Depending on the given parameters this will create
  !! quadrature points on the reference volume element or the reference
  !! surface element.
  interface ply_create_volume_gauss_points_cube
    module procedure ply_create_volume_gauss_points_cube
  end interface

  public :: ply_space_quadrature_type,               &
    &       ply_create_volume_gauss_points_cube,     &
    &       ply_create_volume_gauss_points_cube_2d,  &
    &       ply_create_volume_gauss_points_cube_1d,  &
    &       ply_create_surface_gauss_points_cube,    &
    &       ply_create_surface_gauss_points_cube_2d, &
    &       ply_create_surface_gauss_points_cube_1d, &
    &       ply_gaussLegPoints


contains


  ! ------------------------------------------------------------------------ !
  !> Create the quadrature points for a Gauss-Legendre quadrature
  !! on the cubic reference volume. The reference cubic element has the space
  !! directions (xi,eta,zeta) (correspond to (x,y,z) in physical space).
  subroutine ply_create_volume_gauss_points_cube(                              &
    &          num_intp_per_direction, points, weights, refElemMin, refElemMax )
    ! -------------------------------------------------------------------- !
    !> Number auf integration points in each direction.
    integer, intent(in) :: num_intp_per_direction

    !> Resulting list of points. First index runs over all points, second
    !! indicates the coordinate dimension (x=1,y=2,z=3).
    real(kind=rk), allocatable, intent(out) :: points(:,:)

    !> Integration weight for each point.
    real(kind=rk), allocatable, intent(out) :: weights(:)

    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin

    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    ! -------------------------------------------------------------------- !
    integer :: i, j, k ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: numQuadPoints
    ! -------------------------------------------------------------------- !

    numQuadPoints = num_intp_per_direction**3
    allocate(points(numQuadPoints,3))
    allocate(weights(numQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call ply_gaussLegPoints( x1    = refElemMin,            &
      &                      x2    = refElemMax,            &
      &                      x     = gaussp1D,              &
      &                      w     = weights1D,             &
      &                      nIntP = num_intp_per_direction )

    pointNumber = 1
    do k = 1, num_intp_per_direction
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          !> here we build all possible combinations of the one-dimensional
          !! quadrature points to get the three dimensional values.
          points(PointNumber, 1 ) = gaussp1D(i)
          points(PointNumber, 2 ) = gaussp1D(j)
          points(PointNumber, 3 ) = gaussp1D(k)

          weights(PointNumber) = weights1D(i) * weights1D(j) * weights1D(k)

          pointNumber = pointNumber + 1
        end do
      end do
    end do

  end subroutine ply_create_volume_gauss_points_cube
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create the quadrature points for a 2D Gauss-Legendre quadrature
  !! on the cubic reference volume. The reference cubic element has the space
  !! directions (xi,eta,zeta) (correspond to (x,y,z) in physical space).
  !! The z-coordinates will be set to 0.
  subroutine ply_create_volume_gauss_points_cube_2d(                           &
    &          num_intp_per_direction, points, weights, refElemMin, refElemMax )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(out) :: points(:,:)
    real(kind=rk),allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    ! -------------------------------------------------------------------- !
    integer :: i, j ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: numQuadPoints
    ! -------------------------------------------------------------------- !

    numQuadPoints = num_intp_per_direction**2
    allocate(points(numQuadPoints,3))
    allocate(weights(numQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call ply_gaussLegPoints( x1    = refElemMin,            &
      &                      x2    = refElemMax,            &
      &                      x     = gaussp1D,              &
      &                      w     = weights1D,             &
      &                      nIntP = num_intp_per_direction )

    pointNumber = 1
    do j = 1, num_intp_per_direction
      do i = 1, num_intp_per_direction
        !> here we build all possible combinations of the one-dimensional
        !! quadrature points to get the three dimensional values.
        points(PointNumber, 1 ) = gaussp1D(i)
        points(PointNumber, 2 ) = gaussp1D(j)
        points(PointNumber, 3 ) = 0.0_rk

        weights(PointNumber) = weights1D(i) * weights1D(j)
        pointNumber = pointNumber + 1
      end do
    end do

  end subroutine ply_create_volume_gauss_points_cube_2d
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create the quadrature points for a 1D Gauss-Legendre quadrature
  !! on the cubic reference volume. The reference cubic element has the space
  !! directions (xi,eta,zeta) (correspond to (x,y,z) in physical space).
  !! The z and y-coordinates will be set to 0.
  subroutine ply_create_volume_gauss_points_cube_1d( &
    &          num_intp_per_direction, points, weights, refElemMin, refElemMax )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(out) :: points(:,:)
    real(kind=rk),allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    ! -------------------------------------------------------------------- !
    integer :: i ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    ! -------------------------------------------------------------------- !

    allocate(points(num_intp_per_direction,3))
    allocate(weights(num_intp_per_direction))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call ply_gaussLegPoints( x1    = refElemMin,            &
      &                      x2    = refElemMax,            &
      &                      x     = gaussp1D,              &
      &                      w     = weights1D,             &
      &                      nIntP = num_intp_per_direction )

    pointNumber = 1
    do i = 1, num_intp_per_direction
      !> here we build all possible combinations of the one-dimensional
      !! quadrature points to get the three dimensional values.
      points(PointNumber, 1 ) = gaussp1D(i)
      points(PointNumber, 2 ) = 0.0_rk
      points(PointNumber, 3 ) = 0.0_rk

      weights(PointNumber) = weights1D(i)
      pointNumber = pointNumber + 1
    end do
  end subroutine ply_create_volume_gauss_points_cube_1d
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create the integration points on the surface of (cubical) elements.
  subroutine ply_create_surface_gauss_points_cube(      &
    &          num_intp_per_direction, points, weights, &
    &          refElemMin, refElemMax,                  &
    &          dir, align                               )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(out) :: points(:,:)
    real(kind=rk),allocatable,  intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: i, j ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: nquadPoints
    ! -------------------------------------------------------------------- !
    nQuadPoints = num_intp_per_direction**2

    allocate(points(nQuadPoints,3))
    allocate(weights(nQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call ply_gaussLegPoints( x1    = refElemMin,            &
      &                      x2    = refElemMax,            &
      &                      x     = gaussp1D,              &
      &                      w     = weights1D,             &
      &                      nIntP = num_intp_per_direction )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          !> here we build all possible combinations of the one-dimensional
          !! quadrature points to get the three dimensional values.
          points(PointNumber, 1) = (-1.0_rk)**align
          points(PointNumber, 2) = gaussp1D(i)
          points(PointNumber, 3) = gaussp1D(j)
          weights(PointNumber) = weights1D(i) * weights1D(j)
          pointNumber = pointNumber + 1
        end do
      end do

    case(2) ! face in y direction, y coord is fixed
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          !> here we build all possible combinations of the one-dimensional
          !! quadrature points to get the three dimensional values.
          points(PointNumber, 1 ) = gaussp1D(i)
          points(PointNumber, 2 ) = (-1.0_rk)**align
          points(PointNumber, 3 ) = gaussp1D(j)
          weights(PointNumber) = weights1D(i) * weights1D(j)
          pointNumber = pointNumber + 1
        end do
      end do

    case(3) ! face in z direction, z coord is fixed
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          !> here we build all possible combinations of the one-dimensional
          !! quadrature points to get the three dimensional values.
          points(PointNumber, 1 ) = gaussp1D(i)
          points(PointNumber, 2 ) = gaussp1D(j)
          points(PointNumber, 3 ) = (-1.0_rk)**align
          weights(PointNumber) = weights1D(i) * weights1D(j)
          pointNumber = pointNumber + 1
        end do
      end do

    case default
      call tem_abort( 'ERROR in create_surface_gauss_points_cube:' &
        & // ' unknown face direction'                             )

    end select

  end subroutine ply_create_surface_gauss_points_cube
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create the integration points on the surface of (cubical) elements.
  subroutine ply_create_surface_gauss_points_cube_2d(   &
    &          num_intp_per_direction, points, weights, &
    &          refElemMin, refElemMax,                  &
    &          dir, align                               )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    real(kind=rk), allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: i ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: nQuadPoints
    ! -------------------------------------------------------------------- !
    ! The number of quadrature points on the boundary of a 2d volume is the
    ! number of quad points in one direction
    nQuadPoints = num_intp_per_direction

    allocate(points(nQuadPoints,3))
    allocate(weights(nQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call ply_gaussLegPoints( x1    = refElemMin,            &
      &                      x2    = refElemMax,            &
      &                      x     = gaussp1D,              &
      &                      w     = weights1D,             &
      &                      nIntP = num_intp_per_direction )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      do i = 1, num_intp_per_direction
        !> here we build all possible combinations of the one-dimensional
        !! quadrature points for 2d case to get the three dimensional values.
        points(PointNumber, 1 ) = (-1.0_rk)**align
        points(PointNumber, 2 ) = gaussp1D(i)
        points(PointNumber, 3 ) = 0.0_rk
        weights(PointNumber) = weights1D(i)
        pointNumber = pointNumber + 1
      end do

    case(2) ! face in y direction, y coord is fixed
      do i = 1, num_intp_per_direction
        !> here we build all possible combinations of the one-dimensional
        !! quadrature points in 2d case to get the three dimensional values.
        points(PointNumber, 1 ) = gaussp1D(i)
        points(PointNumber, 2 ) = (-1.0_rk)**align
        points(PointNumber, 3 ) = 0.0_rk
        weights(PointNumber) = weights1D(i)
        pointNumber = pointNumber + 1
      end do

    case default
      call tem_abort( 'ERROR in create_surface_gauss_points_cube_2d:' &
        & // ' unknown face direction'                                )

    end select

  end subroutine ply_create_surface_gauss_points_cube_2d
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create the integration points on the surface of (cubical) elements.
  subroutine ply_create_surface_gauss_points_cube_1d( points, weights, dir, &
    &                                                 align                 )
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    real(kind=rk), allocatable, intent(out) :: weights(:)
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: pointNumber
    integer :: nQuadPoints
    ! -------------------------------------------------------------------- !
    ! for 1d case, the number of points on the  boundary is 1
    nQuadPoints = 1
    allocate(points(nQuadPoints,3))
    allocate(weights(nQuadPoints))

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
            ! for 1d case, there is no combination of 1d quadrature points
      points(pointNumber, 1 ) = (-1.0_rk)**align
      points(pointNumber, 2 ) = 0.0_rk
      points(pointNumber, 3 ) = 0.0_rk
      ! for 1d case, the weights of the point on the boundary is
      ! automatically 1
      ! weights(pointNumber) = 1.0_rk

    case default
      call tem_abort( 'ERROR in create_surface_cheb_points_cube_1d:' &
        & // ' unknown face direction'                               )

    end select

  end subroutine ply_create_surface_gauss_points_cube_1d
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Create Gauss-Legendre integration points and weights for one-dimensional
  !! integration on the interval [x1,x2].
  subroutine ply_gaussLegPoints( x1, x2, x, w, nIntP )
    ! -------------------------------------------------------------------- !
    !> lower limit of integration interval
    real(kind=rk), intent(in) :: x1
    !> upper limit of integration interval
    real(kind=rk), intent(in) :: x2
    !> The coordinates of the gauss points on the interval [x1,x2].
    !! The array has the length nIntP.
    real(kind=rk), intent(out) :: x(:)
    !> The quadrature weights. The array has the length nIntP.
    real(kind=rk), intent(out) :: w(:)
    !> The number of integration points.
    integer, intent(in) :: nIntP
    ! -------------------------------------------------------------------- !
    !> some working variables
    real(kind=rk) :: z1,z,xm,xl,pp,p3,p2,p1;
    !> the relative precision of the points
    real(kind=rk) :: EPS
    integer :: m, i, j
    ! -------------------------------------------------------------------- !

    EPS = 1.0 / (10.0**(PRECISION(1.0_rk)-2) )
    m = (nIntP+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)

    do i = 1, m

      z = cos(PI*((i-1)+0.75_rk)/(nIntP+0.5_rk))

      loopToExit: do
        p1 = 1.0_rk
        p2 = 0.0_rk
        do j=0 , nIntP-1
          p3 = p2
          p2 = p1
          p1 = ((2.0_rk*j+1.0_rk)*z*p2-j*p3) / (j+1.0_rk)
        end do
        pp = nIntP*(z*p1-p2)/(z*z-1.0_rk)
        z1 = z
        z = z1-p1/pp
        if ( abs(z-z1) < EPS ) EXIT loopToExit
      end do loopToExit

      x(i) = xm-xl*z
      x(nIntP-i+1) = xm+xl*z
      w(i) = 2.0_rk*xl/((1.0_rk-z*z)*pp*pp)
      w(nIntp-i+1) = w(i)

    end do

  end subroutine ply_gaussLegPoints
  ! ------------------------------------------------------------------------ !

end module ply_space_integration_module
