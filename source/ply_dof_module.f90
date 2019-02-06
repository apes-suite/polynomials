!> Module provides subroutines, functions and datatypes regarding
!! cell local degrees of freedoms.
module ply_dof_module
  use env_module, only: rk

  implicit none

  private

  public :: ply_dof_nextCoeff, ply_dof_pos, ply_dof_count, ply_dof_2degree, &
    &       ply_degree_2dof

  abstract interface
    subroutine ply_dof_nextCoeff(ansFuncX, ansFuncY, ansFuncZ, maxDegree)
      integer, intent(inout) :: ansFuncX, ansFuncY, ansFuncZ
      integer, intent(in) :: maxdegree
    end subroutine ply_dof_nextCoeff

    function ply_dof_pos(ansFuncX, ansFuncY, ansFuncZ, maxDegree) result(pos)
      integer, intent(in) :: ansFuncX, ansFuncY, ansFuncZ
      integer, intent(in) :: maxdegree
      integer :: pos
    end function ply_dof_pos

    function ply_dof_count(maxPolyDegree) result(DoFs)
      integer, intent(in) :: maxPolyDegree
      integer :: DoFs
    end function ply_dof_count
  end interface

  !> Parameter to identify Q polynomials
  integer, public, parameter :: Q_space = 1
  !> Parameter to identify P polynomials
  integer, public, parameter :: P_space = 2

contains

  elemental function ply_dof_2degree(ndofs, space, ndims) result(deg)
    integer, intent(in) :: ndofs
    integer, intent(in) :: space
    integer, intent(in) :: ndims
    integer :: deg

    integer :: idim, fact, estimate

    select case(space)
    case (Q_space)
      deg = nint(ndofs**(1._rk/real(ndims,kind=rk))) - 1
    case (P_space)
      deg = 0
      do
        estimate = ply_degree_2dof(deg, space, nDims)
        if (estimate >= ndofs) then
          EXIT
        end if
        deg = deg + 1
      end do
    end select

  end function ply_dof_2degree

  elemental function ply_degree_2dof(deg, space, nDims) result(nDofs)
    integer, intent(in) :: deg
    integer, intent(in) :: space
    integer, intent(in) :: nDims
    integer :: nDofs

    select case(space)
    case (Q_space)
      nDofs = (deg+1)**nDims
    case (P_space)
        select case (nDims)
        case(3)
          nDofs = ((deg + 1)  &
            &   *  (deg + 2)  &
            &   *  (deg + 3)) &
            &   / 6
        case(2)
          nDofs = ((deg + 1)  &
            &   *  (deg + 2)) &
            &   / 3
        case(1)
          nDofs = (deg + 1)
        end select
    end select

  end function ply_degree_2dof

end module ply_dof_module
