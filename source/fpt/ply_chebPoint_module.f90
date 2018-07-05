!> Routines and datatypes related to Chebyshev points
!! are located in this module.
!! \author{Jens Zudrop}
module ply_chebPoint_module
  use env_module,         only: rk
  use tem_param_module,   only: PI
  use tem_aux_module,     only: tem_abort
  use tem_logging_module, only: logUnit

  implicit none

  private

  public :: create_volume_cheb_points_cube
  public :: create_surface_cheb_points_cube
  public :: create_volume_cheb_points_cube_1d
  public :: create_surface_cheb_points_cube_1d
  public :: create_volume_lobattocheb_points_cube_1d
  public :: create_surface_lobattocheb_points_cube_1d
  public :: create_volume_cheb_points_cube_2d
  public :: create_surface_cheb_points_cube_2d
  public :: create_volume_lobattocheb_points_cube_2d
  public :: create_surface_lobattocheb_points_cube_2d
  public :: create_volume_lobattocheb_points_cube
  public :: create_surface_lobattocheb_points_cube


contains


  ! ************************************************************************ !
  !> Generates a given number of Chebyshev points on the unit interval [-1;+1].
  subroutine ply_chebPoint_1D( nPoints, chebPnt1D )
    ! -------------------------------------------------------------------- !
    !> The number of points to generate
    integer, intent(in) :: nPoints
    !> Array with 1D Chebyshev points on the interval [-1,+1]
    real(kind=rk), intent(out), allocatable :: chebPnt1D(:)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iPoint)

    allocate(chebPnt1D(nPoints))
    !$OMP DO
    do iPoint = 1, nPoints
      chebPnt1D(iPoint) = -1.0_rk &
        & * cos( PI / nPoints * ( (iPoint - 1.0_rk) + 1.0_rk / 2.0_rk ) )
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine ply_chebPoint_1D
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Generates a given number of Lobatto-Chebyshev points on the unit interval
  !! [-1;+1].
  !!
  !! The Lobatto-Chebyshev points are ordered from +1 to -1. The fpt
  !! have to be adapted accordingly.
  subroutine ply_lobattoChebPoint_1D( nPoints, chebPnt1D )
    ! -------------------------------------------------------------------- !
    !> The number of points to generate
    integer, intent(in) :: nPoints
    !> Array with 1D Chebyshev points on the interval [-1,+1]
    real(kind=rk), intent(out), allocatable :: chebPnt1D(:)
    ! -------------------------------------------------------------------- !
    integer :: iPoint
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iPoint)

    allocate(chebPnt1D(nPoints))
    !$OMP DO
    do iPoint = 1, nPoints
      chebPnt1D(iPoint) = cos( ( iPoint - 1.0_rk ) * PI / ( nPoints - 1.0_rk ) )
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine ply_lobattoChebPoint_1D
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine create_volume_cheb_points_cube(num_intp_per_direction, points)
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    ! -------------------------------------------------------------------- !
    integer :: i, j, k ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(k, j, i, pointNumber)
    nquadpoints = num_intp_per_direction**3
    allocate(points(nquadpoints,3))

    call ply_chebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1
    !$OMP DO
    do k = 1, num_intp_per_direction
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          ! here we build all possible combinations of the one-dimensional
          ! points to get the three dimensional values.
          points(pointNumber, 1) = chebPnt1D(i)
          points(pointNumber, 2) = chebPnt1D(j)
          points(pointNumber, 3) = chebPnt1D(k)
          pointNumber = pointNumber + 1
        end do
      end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine create_volume_cheb_points_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine create_volume_lobattocheb_points_cube( num_intp_per_direction, &
    &                                               points                  )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    ! -------------------------------------------------------------------- !
    integer :: i, j, k ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(k, j, i, pointNumber)
    nquadpoints = num_intp_per_direction**3
    allocate(points(nquadpoints,3))

    call ply_lobattoChebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1
    !$OMP DO
    do k = 1, num_intp_per_direction
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          ! here we build all possible combinations of the one-dimensional
          ! points to get the three dimensional values.
          points(pointNumber, 1) = chebPnt1D(i)
          points(pointNumber, 2) = chebPnt1D(j)
          points(pointNumber, 3) = chebPnt1D(k)
          pointNumber = pointNumber + 1
        end do
      end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine create_volume_lobattocheb_points_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine create_volume_cheb_points_cube_2d( num_intp_per_direction, points )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    ! -------------------------------------------------------------------- !
    !> loop indices
    integer :: i, j
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(j, i, pointNumber)


    nquadpoints = num_intp_per_direction**2
    allocate(points(nquadpoints,3))

    call ply_chebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1
    !$OMP DO
    do j = 1, num_intp_per_direction
      do i = 1, num_intp_per_direction
        ! here we build all possible combinations of the one-dimensional
        ! points to get the three dimensional values.
        points(pointNumber, 1) = chebPnt1D(i)
        points(pointNumber, 2) = chebPnt1D(j)
        points(pointNumber, 3) = 0.0_rk
        pointNumber = pointNumber + 1
      end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine create_volume_cheb_points_cube_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Creates Lobatto-Chebyshev points (with 3 coordinates) but
  !! for 2D kernels.
  subroutine create_volume_lobattocheb_points_cube_2d(num_intp_per_direction, points)
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    ! -------------------------------------------------------------------- !
    !> loop indices
    integer :: i, j
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(j, i, pointNumber)
    nquadpoints = num_intp_per_direction**2
    allocate(points(nquadpoints,3))

    call ply_lobattoChebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1
    !$OMP DO
    do j = 1, num_intp_per_direction
      do i = 1, num_intp_per_direction
        ! here we build all possible combinations of the one-dimensional
        ! points to get the three dimensional values.
        points(pointNumber, 1) = chebPnt1D(i)
        points(pointNumber, 2) = chebPnt1D(j)
        points(pointNumber, 3) = 0.0_rk
        pointNumber = pointNumber + 1
      end do
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine create_volume_lobattocheb_points_cube_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine create_volume_cheb_points_cube_1d( num_intp_per_direction, points )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk), allocatable, intent(inout) :: points(:,:)
    ! -------------------------------------------------------------------- !
    integer :: i ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i, pointNumber)

    allocate(points(num_intp_per_direction,3))
    call ply_chebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1
    !$OMP DO
    do i = 1, num_intp_per_direction
      ! here we build all possible combinations of the one-dimensional
      ! points to get the three dimensional values.
      points(pointNumber, 1) = chebPnt1D(i)
      points(pointNumber, 2) = 0.0_rk
      points(pointNumber, 3) = 0.0_rk
      pointNumber = pointNumber + 1
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine create_volume_cheb_points_cube_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Creates Lobatto-Chebyshev points (with 3 coordinates) but
  !! for 1D kernels.
  subroutine create_volume_lobattocheb_points_cube_1d( num_intp_per_direction, &
    &                                                  points                  )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    ! -------------------------------------------------------------------- !
    integer :: i ! loop indices
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i, pointNumber)

    allocate(points(num_intp_per_direction,3))
    call ply_lobattoChebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1
    !$OMP DO
    do i = 1, num_intp_per_direction
      ! here we build all possible combinations of the one-dimensional
      ! points to get the three dimensional values.
      points(pointNumber, 1) = chebPnt1D(i)
      points(pointNumber, 2) = 0.0_rk
      points(pointNumber, 3) = 0.0_rk
      pointNumber = pointNumber + 1
    end do
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine create_volume_lobattocheb_points_cube_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Create Chebyshev nodes on a specific face of the reference element.
  subroutine create_surface_cheb_points_cube( num_intp_per_direction, points, &
    &                                         dir, align                      )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
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
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(j, i, pointNumber)

    nquadpoints = num_intp_per_direction**2
    allocate(points(nquadpoints,3))

    call ply_chebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      !$OMP DO
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          points(pointNumber, 1) = (-1.0_rk)**align
          points(pointNumber, 2) = chebPnt1D(i)
          points(pointNumber, 3) = chebPnt1D(j)
          pointNumber = pointNumber + 1
        end do
      end do
      !$OMP END DO
    case(2) ! face in y direction, y coord is fixed
      !$OMP DO
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          points(pointNumber, 1) = chebPnt1D(i)
          points(pointNumber, 2) = (-1.0_rk)**align
          points(pointNumber, 3) = chebPnt1D(j)
          pointNumber = pointNumber + 1
        end do
      end do
      !$OMP END DO
    case(3) ! face in z direction, z coord is fixes
      !$OMP DO
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          points(pointNumber , 1) = chebPnt1D(i)
          points(pointNumber , 2) = chebPnt1D(j)
          points(pointNumber , 3) = (-1.0_rk)**align
          pointNumber = pointNumber + 1
        end do
      end do
      !$OMP END DO
    case default
      call tem_abort( 'ERROR in create_surface_cheb_points_cube: unknown ' &
        & // 'face direction'                                              )
    end select

    !$OMP END PARALLEL

  end subroutine create_surface_cheb_points_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Create Lobatto-Chebyshev nodes on a specific face of the reference element.
  subroutine create_surface_lobattocheb_points_cube( num_intp_per_direction, &
    &                                                points, dir, align      )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
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
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(j, i, pointNumber)

    nquadpoints = num_intp_per_direction**2
    allocate(points(nquadpoints,3))

    call ply_lobattoChebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      !$OMP DO
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          points(pointNumber, 1) = (-1.0_rk)**align
          points(pointNumber, 2) = chebPnt1D(i)
          points(pointNumber, 3) = chebPnt1D(j)
          pointNumber = pointNumber + 1
        end do
      end do
      !$OMP END DO
    case(2) ! face in y direction, y coord is fixed
      !$OMP DO
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          points(pointNumber, 1) = chebPnt1D(i)
          points(pointNumber, 2) = (-1.0_rk)**align
          points(pointNumber, 3) = chebPnt1D(j)
          pointNumber = pointNumber + 1
        end do
      end do
      !$OMP END DO
    case(3) ! face in z direction, z coord is fixes
      !$OMP DO
      do j = 1, num_intp_per_direction
        do i = 1, num_intp_per_direction
          points(pointNumber, 1) = chebPnt1D(i)
          points(pointNumber, 2) = chebPnt1D(j)
          points(pointNumber, 3) = (-1.0_rk)**align
          pointNumber = pointNumber + 1
        end do
      end do
      !$OMP END DO
    case default
      call tem_abort( 'ERROR in create_surface_lobattocheb_points_cube:' &
        & // ' unknown face direction'                                   )
    end select

    !$OMP END PARALLEL

  end subroutine create_surface_lobattocheb_points_cube
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Create Chebyshev nodes on a specific face of the reference element.
  subroutine create_surface_cheb_points_cube_2d( num_intp_per_direction, &
    &                                            points, dir, align      )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: i
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i, pointNumber)

    nquadpoints = num_intp_per_direction
    allocate(points(nquadpoints,3))

    call ply_chebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      !$OMP DO
      do i = 1, num_intp_per_direction
        points(pointNumber, 1) = (-1.0_rk)**align
        points(pointNumber, 2) = chebPnt1D(i)
        points(pointNumber, 3) = 0.0_rk
        pointNumber = pointNumber + 1
      end do
      !$OMP END DO
    case(2) ! face in y direction, y coord is fixed
      !$OMP DO
      do i = 1, num_intp_per_direction
        points(pointNumber, 1) = chebPnt1D(i)
        points(pointNumber, 2) = (-1.0_rk)**align
        points(pointNumber, 3) = 0.0_rk
        pointNumber = pointNumber + 1
      end do
      !$OMP END DO
    case default
      call tem_abort( 'ERROR in create_surface_cheb_points_cube_2d:' &
        & // ' unknown face direction'                               )
    end select

    !$OMP END PARALLEL

  end subroutine create_surface_cheb_points_cube_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Create Lobatto-Chebyshev nodes on a specific face of the reference element.
  subroutine create_surface_lobattocheb_points_cube_2d( &
    & num_intp_per_direction, points, dir, align        )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: num_intp_per_direction
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: i
    integer :: pointNumber
    real(kind=rk), allocatable :: chebPnt1D(:)
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i, pointNumber)

    nquadpoints = num_intp_per_direction
    allocate(points(nquadpoints,3))

    call ply_lobattoChebPoint_1D( num_intp_per_direction, chebPnt1D )

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      !$OMP DO
      do i = 1, num_intp_per_direction
        points(pointNumber, 1) = (-1.0_rk)**align
        points(pointNumber, 2) = chebPnt1D(i)
        points(pointNumber, 3) = 0.0_rk
        pointNumber = pointNumber + 1
      end do
      !$OMP END DO
    case(2) ! face in y direction, y coord is fixed
      !$OMP DO
      do i = 1, num_intp_per_direction
        points(pointNumber, 1) = chebPnt1D(i)
        points(pointNumber, 2) = (-1.0_rk)**align
        points(pointNumber, 3) = 0.0_rk
        pointNumber = pointNumber + 1
      end do
      !$OMP END DO
    case default
      call tem_abort( 'ERROR in create_surface_cheb_points_cube_2d:' &
        & // ' unknown face direction'                               )
    end select

    !$OMP END PARALLEL

  end subroutine create_surface_lobattocheb_points_cube_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Create Chebyshev nodes on a specific face of the reference element.
  subroutine create_surface_cheb_points_cube_1d( points, dir, align )
    ! -------------------------------------------------------------------- !
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: pointNumber
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !

    nquadpoints = 1
    allocate(points(nquadpoints,3))

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      points(pointNumber, 1) = (-1.0_rk)**align
      points(pointNumber, 2) = 0.0_rk
      points(pointNumber, 3) = 0.0_rk
    case default
      call tem_abort( 'ERROR in create_surface_cheb_points_cube_1d:' &
        & // ' unknown face direction'                               )
    end select

  end subroutine create_surface_cheb_points_cube_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Create Lobatto-Chebyshev nodes on a specific face of the reference element.
  subroutine create_surface_lobattocheb_points_cube_1d( points, dir, align )
    ! -------------------------------------------------------------------- !
    real(kind=rk),allocatable, intent(inout) :: points(:,:)
    !> The spatial direction of the face. \n
    !! 1 -> x direction \n
    !! 2 -> y direction \n
    integer :: dir
    !> Left or right face of the reference element
    integer :: align
    ! -------------------------------------------------------------------- !
    integer :: pointNumber
    integer :: nquadpoints
    ! -------------------------------------------------------------------- !

    nquadpoints = 1
    allocate(points(nquadpoints,3))

    pointNumber = 1

    select case(dir)
    case(1) ! face in x direction, x coord is fixed
      points(pointNumber, 1) = (-1.0_rk)**align
      points(pointNumber, 2) = 0.0_rk
      points(pointNumber, 3) = 0.0_rk
    case default
      call tem_abort( 'ERROR in create_surface_cheb_points_cube_1d:' &
        & // ' unknown face direction'                               )
    end select

  end subroutine create_surface_lobattocheb_points_cube_1d
  ! ************************************************************************ !


end module ply_chebPoint_module

