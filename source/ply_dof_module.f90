!> Module provides subroutines, functions and datatypes regarding
!! cell local degrees of freedoms.
module ply_dof_module
  use env_module, only: rk

  implicit none

  private

  !! The routine posOfModgCoeffQTens is not available anymore. Use the macro
  !! from ply_dof_module.inc instead.
  public :: nextModgCoeffQTens
  public :: nextModgCoeffPTens
  public :: nextModgCoeffPTens2D
  public :: nextModgCoeffQTens2D
  public :: nextModgCoeffPTens1D
  public :: nextModgCoeffQTens1D

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

  public :: ply_dof_nextCoeff, ply_dof_pos, ply_dof_count

  !> Parameter to identify Q polynomials
  integer, public, parameter :: Q_space = 1
  !> Parameter to identify P polynomials
  integer, public, parameter :: P_space = 2


contains

  !> The x, y and z ansatz degrees are turned into the degrees of the next
  !! ansatz function in the linearized Q tensor
  pure subroutine nextModgCoeffQTens(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncZ
    !> Maximal polynomial degree
    integer, intent(in) :: maxdegree

    integer :: PolyOrd

    PolyOrd = maxdegree + 1

    if (ansFuncX .ne. PolyOrd) then
      ! next X index
      ansFuncX = ansFuncX +1
    elseif (ansFuncY .ne. PolyOrd) then
      ! next y index
      ansFuncX = 1
      ansFuncY = ansFuncY+1
    else
      ! next Z index
      ansFuncX = 1
      ansFuncY = 1
      ansFuncZ = ansFuncZ+1
    end if
  end subroutine nextModgCoeffQTens


  !> The x, y and z ansatz degrees are turned into the degrees of the next
  !! ansatz function in the layered P list
  pure subroutine nextModgCoeffPTens(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncZ
    !> Maximal polynomial degree in each direction
    integer, intent(in) :: maxdegree

    integer :: polyOrd

    polyOrd = maxdegree + 1
    ! - ansatz indices are arranged in layers. Within each layer, the total
    !   degree remains constant.
    ! - within each layer, we have blocks. Within a block, ansFuncZ is
    !   constant, y counts up and x accordingly down.
    ! - within each block, we have items. Each item represents one particular
    !   combination of ansFuncX, -Y, and -Z degrees.

    if (ansFuncX .ne. 1) then
      ! next item
      ansFuncX = ansFuncX -1
      ansFuncY = ansFuncY +1
    elseif (ansFuncY .ne. 1) then
      ! next block
      ansFuncX = ansFuncY -1
      ansFuncY = 1
      ansFuncZ = ansFuncZ +1
    else
      ! next layer
      ansFuncX = ansFuncZ +1
      ansFuncY = 1
      ansFuncZ = 1
    end if
  end subroutine nextModgCoeffPTens


  !> The x and y ansatz degrees are turned into the degrees of the next
  !! ansatz function in the linearized Q tensor
  pure subroutine nextModgCoeffQTens2D(ansFuncX, ansFuncY, maxdegree)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY
    !> Maximal polynomial degree
    integer, intent(in) :: maxdegree

    integer :: PolyOrd

    PolyOrd = maxdegree + 1

    if (ansFuncX .ne. PolyOrd) then
      ! next X index
      ansFuncX = ansFuncX +1
    else
      ! next y index
      ansFuncX = 1
      ansFuncY = ansFuncY+1
    end if
  end subroutine nextModgCoeffQTens2D


  !> The x and y ansatz degrees are turned into the degrees of the next
  !! ansatz function in the layered P list
  pure subroutine nextModgCoeffPTens2D(ansFuncX, ansFuncY)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY

    ! - ansatz indices are arranged in layers. Within each layer, the total
    !   degree remains constant.
    ! - within each layer, we have blocks. Within a block, ansFuncZ is
    !   constant, y counts up and x accordingly down.
    ! - within each block, we have items. Each item represents one particular
    !   combination of ansFuncX, -Y, and -Z degrees.

    if (ansFuncX .ne. 1) then
      ! next item
      ansFuncX = ansFuncX -1
      ansFuncY = ansFuncY +1
    else
      ! next layer
      ansFuncX = ansFuncY + 1
      ansFuncY = 1
    end if

  end subroutine nextModgCoeffPTens2D


  !> The x ansatz degree is turned into the degree of the next
  !! ansatz function in the linearized Q tensor
  pure subroutine nextModgCoeffQTens1D(ansFuncX)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX

    ansFuncX = ansFuncX +1

  end subroutine nextModgCoeffQTens1D


  !> The x ansatz degree is turned into the degree of the next
  !! ansatz function in the layered P list
  pure subroutine nextModgCoeffPTens1D(ansFuncX)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX

    ansFuncX = ansFuncX +1

  end subroutine nextModgCoeffPTens1D


end module ply_dof_module
