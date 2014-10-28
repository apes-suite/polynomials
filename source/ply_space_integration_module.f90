!> module collecting all datatypes and subroutines related to the space
!! integration in pnpm scheme.
module ply_space_integration_module
  use env_module,         only: rk
  use tem_param_module,   only: PI
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit
  implicit none

  private

  !> Type to specify the space integration of our variational method.
  type space_quadrature_type
    !> Unique characters specifying the space integration type.
    !! If you add an additional space integration method, then you should
    !! add a unique abbreviation for that here. The following relations
    !! hold:
    !! int_type == 'GALEINT' means Gauss-Legendre interpolation
    character(len=32)        :: int_type

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
    !> The dimensions is numVolumePoints.
    real(kind=rk), allocatable :: integration_weights_volume(:)

    !> The integration points to evaluate a surface integral over the
    !! reference element.
    !! The first dimension of this array is numSufacePoints.
    !! The second dimens is 2
    real(kind=rk), allocatable :: integration_points_surface(:,:)

    !> Weights used for quadrature formula of surface integration.
    !> The dimensions is numSurfacePoints.
    real(kind=rk), allocatable :: integration_weights_surface(:)

    !> This vector holds all gauss point coordinates for a one dimensional
    !! quadrature. Additionally the start and end point of the unit interval
    !! (i.e. 0.0 and 1.0) are appended.
    real(kind=rk), allocatable :: integration_points_line(:)

    !> \todo add additional stuff related to the space integration here!

  end type space_quadrature_type

  !> interface to create gauss points for pure gauss quadrature on the
  !! reference element. Depending on the given parameters this will create
  !! quadrature points on the reference volume element or the reference
  !! surface element.
  interface create_volume_gauss_points_cube
    module procedure create_volume_gauss_points_cube
  end interface

  public :: space_quadrature_type,               & 
    &       create_volume_gauss_points_cube,     &
    &       create_volume_gauss_points_cube_2d,  &
    &       create_volume_gauss_points_cube_1d,  &
    &       create_surface_gauss_points_cube,    &
    &       create_surface_gauss_points_cube_2d, &
    &       create_surface_gauss_points_cube_1d, &
    &       gauleg,                              &
    &       create_gauss_points_1d

contains

  !> subroutine to create the quadrature points for a pure gauss quadrature
  !! on the cubic reference volume. The reference cubic element has the space
  !! directions (xi,eta,zeta) (correspond to (x,y,z) in physical space).
  !! \param num_intp_per_direction the number of gauss points per space
  !!        direction
  !! \param points one dimension array which is the ouput of this routine. The
  !!        gauss points are written to this array, so its dimension is:
  !!        num_intp_per_direction * num_intp_per_direction * num_zeta.
  !!        The gauss points are ordered in the following way:
  !!        \f$ points(1,1:3) = (\xi_1 , \eta_1, \zeta_1 ),
  !!            points(2,1:3) = (\xi_2 , \eta_1, \zeta_1 ), \cdots
  !!            (\xi_{numintpperdirection}, \eta_1, \zeta_1 ),
  !!            (\xi_1 , \eta_2, \zeta_1 ), \cdots
  !!            (\xi_{numintpperdirection}, \eta_{numintpperdirection},
  !!             \zeta_1 ), \cdots
  !!            points(numintpperdirection^3,1:3) = (\xi_{numintpperdirection},
  !!                                           \eta_{numintpperdirection},
  !!                                           \zeta_{numintpperdirection} ) \f$
  subroutine create_volume_gauss_points_cube( num_intp_per_direction,      &
    &                                         points, weights, refElemMin, &
    &                                         refElemMax                   )
    !---------------------------------------------------------------------------
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    real(kind=rk), allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !---------------------------------------------------------------------------
    integer :: i, j, k !> loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: numQuadPoints
    !---------------------------------------------------------------------------

    numQuadPoints = num_intp_per_direction**3
    allocate(points(numQuadPoints,3))
    allocate(weights(numQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction   )

    pointNumber = 1
    do k = 1, num_intp_per_direction
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          !> here we build all possible combinations of the one-dimensional
          !! quadrature points to get the three dimensional values.
          points(PointNumber , 1 ) = gaussp1D(i)
          points(PointNumber , 2 ) = gaussp1D(j)
          points(PointNumber , 3 ) = gaussp1D(k)

          weights(PointNumber) = weights1D(i) * weights1D(j) * weights1D(k)

          pointNumber = pointNumber + 1
        end do
      end do
    end do
  end subroutine create_volume_gauss_points_cube


  subroutine create_volume_gauss_points_cube_2d( num_intp_per_direction,      &
    &                                            points, weights, refElemMin, &
    &                                            refElemMax                   )
    !---------------------------------------------------------------------------
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(out) :: points(:,:)
    real(kind=rk),allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !---------------------------------------------------------------------------
    integer :: i, j !> loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: numQuadPoints
    !---------------------------------------------------------------------------

    numQuadPoints = num_intp_per_direction**2
    allocate(points(numQuadPoints,3))
    allocate(weights(numQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction   )

    pointNumber = 1
    do j = 1, num_intp_per_direction
      do i = 1, num_intp_per_direction
        !> here we build all possible combinations of the one-dimensional
        !! quadrature points to get the three dimensional values.
        points(PointNumber , 1 ) = gaussp1D(i)
        points(PointNumber , 2 ) = gaussp1D(j)
        points(PointNumber , 3 ) = 0.0_rk

        weights(PointNumber) = weights1D(i) * weights1D(j) 
        pointNumber = pointNumber + 1
      end do
    end do
  end subroutine create_volume_gauss_points_cube_2d


  subroutine create_volume_gauss_points_cube_1d( num_intp_per_direction,      &
    &                                            points, weights, refElemMin, &
    &                                            refElemMax                   )
    !---------------------------------------------------------------------------
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(out) :: points(:,:)
    real(kind=rk),allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !---------------------------------------------------------------------------
    integer :: i ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    !---------------------------------------------------------------------------

    allocate(points(num_intp_per_direction,3))
    allocate(weights(num_intp_per_direction))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction   )

    pointNumber = 1
    do i = 1, num_intp_per_direction
      !> here we build all possible combinations of the one-dimensional
      !! quadrature points to get the three dimensional values.
      points(PointNumber , 1 ) = gaussp1D(i)
      points(PointNumber , 2 ) = 0.0_rk
      points(PointNumber , 3 ) = 0.0_rk

      weights(PointNumber) = weights1D(i)  
      pointNumber = pointNumber + 1
    end do
  end subroutine create_volume_gauss_points_cube_1d


  subroutine create_surface_gauss_points_cube( num_intp_per_direction,      &
    &                                          points, weights, refElemMin, &
    &                                          refElemMax, dir, align       )
    !---------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------
    integer :: i, j!> loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: nquadPoints
    !---------------------------------------------------------------------------
    nQuadPoints = num_intp_per_direction**2

    allocate(points(nQuadPoints,3))
    allocate(weights(nQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction   )

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
          points(PointNumber , 1 ) = gaussp1D(i)
          points(PointNumber , 2 ) = (-1.0_rk)**align
          points(PointNumber , 3 ) = gaussp1D(j)
          weights(PointNumber) = weights1D(i) * weights1D(j)
          pointNumber = pointNumber + 1
        end do
      end do

    case(3) ! face in z direction, z coord is fixed
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          !> here we build all possible combinations of the one-dimensional
          !! quadrature points to get the three dimensional values.
          points(PointNumber , 1 ) = gaussp1D(i)
          points(PointNumber , 2 ) = gaussp1D(j)
          points(PointNumber , 3 ) = (-1.0_rk)**align
          weights(PointNumber) = weights1D(i) * weights1D(j)
          pointNumber = pointNumber + 1
        end do
      end do

    case default
      write(logUnit(1),*) 'ERROR in create_surface_gauss_points_cube:' &
        &                 // ' unknown face direction, stopping ... '
      call tem_abort()

    end select

  end subroutine create_surface_gauss_points_cube


  subroutine create_surface_gauss_points_cube_2d(num_intp_per_direction,      &
    &                                            points, weights, refElemMin, &
    &                                            refElemMax, dir, align       )
    !---------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------
    integer :: i ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    integer :: nQuadPoints
    !---------------------------------------------------------------------------
    ! The number of quadrature points on the boundary of a 2d volume is the
    ! number of quad points in one direction
    nQuadPoints = num_intp_per_direction

    allocate(points(nQuadPoints,3))
    allocate(weights(nQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction   )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      do i = 1, num_intp_per_direction
        !> here we build all possible combinations of the one-dimensional
        !! quadrature points for 2d case to get the three dimensional values.
        points(PointNumber , 1 ) = (-1.0_rk)**align
        points(PointNumber , 2 ) = gaussp1D(i)
        points(PointNumber , 3 ) = 0.0_rk
        weights(PointNumber) = weights1D(i)
        pointNumber = pointNumber + 1
      end do

    case(2) ! face in y direction, y coord is fixed
      do i = 1, num_intp_per_direction
        !> here we build all possible combinations of the one-dimensional
        !! quadrature points in 2d case to get the three dimensional values.
        points(PointNumber , 1 ) = gaussp1D(i)
        points(PointNumber , 2 ) = (-1.0_rk)**align
        points(PointNumber , 3 ) = 0.0_rk
        weights(PointNumber) = weights1D(i) 
        pointNumber = pointNumber + 1
      end do

    case default
      write(logUnit(1),*) 'ERROR in create_surface_gauss_points_cube_2d:' &
        &                 // ' unknown face direction, stopping ... '
      call tem_abort()

    end select

  end subroutine create_surface_gauss_points_cube_2d


  subroutine create_surface_gauss_points_cube_1d(points, weights, dir, align)
    !---------------------------------------------------------------------------
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    real(kind=rk), allocatable, intent(out) :: weights(:)
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    !! 3 -> z direction
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    !---------------------------------------------------------------------------
    integer :: pointNumber
    integer :: nQuadPoints
    !---------------------------------------------------------------------------
    ! for 1d case, the number of points on the  boundary is 1
    nQuadPoints = 1
    allocate(points(nQuadPoints,3))
    allocate(weights(nQuadPoints))

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
            ! for 1d case, there is no combination of 1d quadrature points
      points(pointNumber , 1 ) = (-1.0_rk)**align
      points(pointNumber , 2 ) = 0.0_rk
      points(pointNumber , 3 ) = 0.0_rk
      ! for 1d case, the weights of the point on the boundary is
      ! automatically 1
      ! weights(pointNumber) = 1.0_rk

    case default
      write(logUnit(1),*) 'ERROR in create_surface_cheb_points_cube_1d:' &
        &                 // ' unknown face direction, stopping ...'
      call tem_abort()

    end select

  end subroutine create_surface_gauss_points_cube_1d


  subroutine create_gauss_points_1d(num_intp_per_direction, numQuadPoints, &
    &                               points, weights, refElemMin, refElemMax)
    !---------------------------------------------------------------------------
    integer, intent(in) :: num_intp_per_direction
    integer, intent(out) :: numQuadPoints
    real(kind=rk), allocatable, intent(out) :: points(:,:)
    real(kind=rk), allocatable, intent(out) :: weights(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !---------------------------------------------------------------------------
    integer :: j!> loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    !---------------------------------------------------------------------------
    numQuadPoints = num_intp_per_direction

    allocate(points(numQuadPoints,1))
    allocate(weights(numQuadPoints))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction)

    pointNumber = 1
    do j = 1, num_intp_per_direction
      !> here we build all possible combinations of the one-dimensional
      !! quadrature points to get the three dimensional values.
      points(PointNumber, 1) = gaussp1D(j)
      weights(PointNumber) = weights1D(j)
      pointNumber = pointNumber + 1
    end do
  end subroutine create_gauss_points_1d


  !> \brief subroutine to create a linearized unique array of surface
  !! integration points with all possible x coordinates for any surface
  !! integration point.
  subroutine create_linearized_line_points( num_intp_per_direction, points, &
    &                                       refElemMin, refElemMax )
    !--------------------------------------------------------------------------
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk), allocatable, intent(out) :: points(:)
    !> Left bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMin
    !> Right bound of the one-dimensional reference element.
    real(kind=rk), intent(in) :: refElemMax
    !--------------------------------------------------------------------------
    real(kind=rk), allocatable :: gaussp1D(:)
    real(kind=rk), allocatable :: weights1D(:)
    !--------------------------------------------------------------------------

    allocate(points(num_intp_per_direction + 2 ))
    allocate(gaussp1D(num_intp_per_direction))
    allocate(weights1D(num_intp_per_direction))

    ! create the one dimensional gauss points first
    call gauleg( x1 = refElemMin, x2 = refElemMax, x = gaussp1D, &
      &          w = weights1D, nIntP = num_intp_per_direction   )
    ! first, copy the one dimensional integration points...
    points(1:num_intp_per_direction) = gaussp1D(:)
    ! second, append the default position on the y-z surface
    points(num_intp_per_direction+1) = refElemMin
    points(num_intp_per_direction+2) = refElemMax

  end subroutine create_linearized_line_points


  !> subroutine to create gauss points and weights for one-dimensional
  !! integration on the interval [x1,x2].
  subroutine gauleg(x1, x2, x, w, nIntP)
    !---------------------------------------------------------------------------
    !> The coordinates of the gauss points on the interval [-1,1].
    !! The array has the length nIntP.
    real(kind=rk), intent(out) :: x(:)
    !> The quadrature weights. The array has the length nIntP.
    real(kind=rk), intent(out) :: w(:)
    !> The number of integration points.
    integer, intent(in) :: nIntP
    !> lower limit of integration interval
    real(kind=rk), intent(in) :: x1
    !> upper limit of integration interval
    real(kind=rk), intent(in) :: x2
    !---------------------------------------------------------------------------
    !> some working variables
    real(kind=rk) :: z1,z,xm,xl,pp,p3,p2,p1;
    !> the relative precision of the points
    real(kind=rk) :: EPS
    integer :: m, i, j
    !---------------------------------------------------------------------------

    EPS= 1.0 / (10.0**(PRECISION(1.0_rk)-2) )
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

  end subroutine gauleg

end module ply_space_integration_module
