module ply_equidistant_module
  use env_module,         only: rk, labelLen
  use fftw_wrap,          only: fftw_available
  use aotus_module,       only: flu_State, aot_get_val
  use tem_precice_module, only: precice
  implicit none
  private

  !> Datatype to represent facewise nodes
  type ply_equadPoints_type
    real(kind=rk), allocatable :: Points(:,:)
    integer :: nPoints 
  end type ply_equadPoints_type
  public :: create_surface_equidistant_points
  public :: ply_equadPoints_type
  contains

 !****************************************************************************!
!> Routine to generate a list of equidistant points
  subroutine create_surface_equidistant_points(me, nDir, iAlign, nPoly, idir)
    !--------------------------------------------------------------------------!
    type(ply_equadPoints_type), intent(out) :: me
    !> store points to the corresponding spatial directions
    integer, intent(in) :: nDir
    integer, intent(in) :: iAlign
    integer, intent(in) :: idir 
    !> polynomial degree
    integer, intent(in) :: nPoly
    !--------------------------------------------------------------------------!
    !> loop variable for the number of intervals
    integer :: iInterval
    !> loop variable for the number of intervals
    integer :: jInterval
    !>factor to increase number of points, default 1.0 
    integer :: factor
    !--------------------------------------------------------------------------!
    factor = precice%factor_EQ_points
    write(*,*) "Chosen factor for the equidistant points :", factor
    if (nDir>1) then
      me%nPoints = nint((nPoly + 1.0) * factor)
    else
      me%nPoints = 1
    end if
    write(*,*) 'Total number of equidistant Points', me%nPoints
    allocate(me%Points(me%nPoints, 3))

    select case (nDir)
    case(1) !1-dimension
      me%Points(1,1) = (-1.0_rk)**iAlign
      me%Points(1,2) = 0.0_rk
      me%Points(1,3) = 0.0_rk

    case(2) !2-dimension
      select case(idir)
      case(1) ! face in x direction, x coord is fixed
        do iInterval = 1, me%nPoints
          me%Points(iInterval, 1 ) = (-1.0_rk)**iAlign
          me%Points(iInterval, 2 ) = (real((iInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
          me%Points(iInterval, 3 ) = 0.0_rk
        end do !iInterval
      case(2) ! face in y direction, y coord is fixed
        do iInterval = 1, me%nPoints 
          me%Points(iInterval, 1 ) = (real((iInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
          me%Points(iInterval, 2 ) = (-1.0_rk)**iAlign
          me%Points(iInterval, 3 ) = 0.0_rk
        end do  
      end select     
    case(3) ! 3-dimension
      select case(idir)
      case(1) ! face in x direction
        do iInterval = 1, me%nPoints
          do jInterval = 1, me%nPoints
            me%Points(iInterval , 1 ) = (-1.0_rk)**iAlign
            me%Points(iInterval , 2 ) = (real((iInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
            me%Points(jInterval , 3 ) = (real((jInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
        end do
          end do
      case(2) ! face in y direction
        do iInterval = 1, me%nPoints
          do jInterval = 1, me%nPoints
            me%Points(iInterval , 1 ) = (real((iInterval -1) * 2 + 1) / me%nPoints) - 1.0_rk
            me%Points(iInterval , 2 ) = (-1.0_rk)**iAlign
            me%Points(jInterval , 3 ) = (real((jInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
        end do 
          end do
      case(3) ! face in z direction, z coord is fixed
        do iInterval = 1, me%nPoints
          do jInterval = 1, me%nPoints
            me%Points(iInterval , 1 ) = (real((iInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
            me%Points(iInterval , 2 ) = (real((jInterval - 1) * 2 + 1) / me%nPoints) - 1.0_rk
            me%Points(jInterval , 3 ) = (-1.0_rk)**iAlign
          end do
        end do
      end select
    end select
  end subroutine create_surface_equidistant_points
!****************************************************************************!

end module ply_equidistant_module
