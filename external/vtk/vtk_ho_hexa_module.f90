! Copyright (c) 2025 Harald Klimach <harald.klimach@dlr.de>
!
! Parts of this file were written by Harald Klimach for DLR e.V.
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

!> Module with functions to work on VTK's high-order
!! hexahedrons.
!!
!! See https://gitlab.kitware.com/vtk/vtk/-/blob/master/Common/DataModel/vtkHigherOrderHexahedron.cxx 
!!
module vtk_ho_hexa_module
  implicit none
  private

  public :: hexa_pointIndexFromIJK


contains


  !> Implementation of vtkHigherOrderHexahedron::PointIndexFromIJK in Fortran
  elemental function hexa_pointIndexFromIJK(i, j, k, orderX, orderY, orderZ) return(pointIndex)
    !> Index in X-direction (0 - order)
    integer, intent(in) :: i
    !> Index in Y-direction (0 - order)
    integer, intent(in) :: j
    !> Index in Z-direction (0 - order)
    integer, intent(in) :: k
    !> Order of the polynomial data in X direction
    integer, intent(in) :: orderX
    !> Order of the polynomial data in Y direction
    integer, intent(in) :: orderY
    !> Order of the polynomial data in Z direction
    integer, intent(in) :: orderZ

    !> Resulting index for the given point in the VTK high-order hexahedron
    integer :: pointIndex

    integer :: surf_count
    integer :: other_edges
    integer :: order(3)
    integer :: ijk(3)
    logical :: on_surface(3)

    ijk(1) = i
    ijk(2) = j
    ijk(3) = k

    order(1) = orderX
    order(2) = orderY
    order(3) = orderZ

    do ind=1,3
      on_surface(ind) = (ijk(ind) == 0 .or. ijk(ind) == order(ind))
    end do
 
    surf_count = count(on_surface)

    pointIndex = 0
    select case(surf_count)
    case(3) ! Vertex DoF
      if (k > 0) then
        pointIndex = pointIndex + 4
      end if
      if (j > 0) then
        pointIndex = pointIndex + 2
      end if
      if (i > 0) then
        pointIndex = pointIndex + 1
      end if

    case(2) ! Edge DoF
      pointIndex = 8
      other_edges = order(1) + order(2) - 2
      if (.not. on_surface(1)) then
        ! On an X-Axis edge
        if (k > 0) pointIndex = pointIndex + 2*other_edges
        if (j > 0) pointIndex = pointIndex + other_edges
        pointIndex = pointIndex + i - 1
      end if

      if (.not. on_surface(2)) then
        ! On a Y-Axis edge
        if (i > 0) then
          pointIndex = pointIndex + order(1) - 1
        else
          pointIndex = 2*(order(1)-1) + order(2) - 1
          if (k > 0) pointIndex = pointIndex + 2*other_edges
        end if
        pointIndex = pointIndex + j - 1
      end if

      if (.not. on_surface(3)) then
        ! On a Z-Axis edge
        pointIndex = pointIndex + 4 * other_edges
        if (i > 0) then
          pointIndex = pointIndex + (order(3)-1)
          if (j > 0) pointIndex = pointIndex + (order(3)-1)
        else
          if (j > 0) pointIndex = pointIndex + 3*(order(3)-1)
        end if
        pointIndex = pointIndex + k - 1
      end if

    case(1) ! Face DoF
      pointIndex = 8 + 4 * (order(1) + order(2) + order(3) - 3)
    case(0) ! Interior DoF
    end select

  end function hexa_pointIndexFromIJK

end module vtk_ho_hexa_module
