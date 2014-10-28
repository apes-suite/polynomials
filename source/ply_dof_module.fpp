?? include "ply_dof_module.inc"
!> Module provides subroutines, functions and datatypes regarding
!! cell local degrees of freedoms.
module ply_dof_module
  use env_module, only: rk

  implicit none

  private

  public :: posOfModgCoeffQTens, nextModgCoeffQTens, getDofsQTens
  public :: posOfModgCoeffPTens, nextModgCoeffPTens, getDofsPTens
  public :: posOfModgCoeffPTens2D, nextModgCoeffPTens2D, getDofsPTens2D
  public :: posOfModgCoeffQTens2D, nextModgCoeffQTens2D, getDofsQTens2D
  public :: posOfModgCoeffPTens1D, nextModgCoeffPTens1D, getDofsPTens1D
  public :: posOfModgCoeffQTens1D, nextModgCoeffQTens1D, getDofsQTens1D

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

  !> Return the position of a given ansatz function combination in the
  !! linearized list of modal coefficients for Q-Tensor product polynomials.
  !! 
  !! There is also a CoCo posOfModgCoeffQTens text available to inline this
  !! function. See ply_dof_module.inc for further details.
  pure function posOfModgCoeffQTens(ansFuncX, ansFuncY, ansFuncZ, &
    &                               maxDegree) result(pos)
    !---------------------------------------------------------------------------
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncZ
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxDegree
    !> The position of the modal coefficient in the list of modal coefficients.
    integer :: pos
    !---------------------------------------------------------------------------

?? copy :: posOfModgCoeffQTens(ansFuncX, ansFuncY, ansFuncZ, maxdegree, pos)

  end function posOfModgCoeffQTens


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


  !> Return the number of degrees of freedom for Q polynomial space
  pure function getDofsQTens(maxPolyDegree) result(dofs)
    !---------------------------------------------------------------------------
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxPolyDegree
    !> The number of degrees of freedom for a Q tensor product polynomial
    integer :: dofs
    !---------------------------------------------------------------------------

    dofs = (maxPolyDegree+1)**3

  end function getDofsQTens


  !> Return the position of a given ansatz function combination in the
  !! linearized list of modal coefficients for P-Tensor product polynomials.
  pure function posOfModgCoeffPTens(ansFuncX, ansFuncY, ansFuncZ, maxdegree) &
    &           result(pos)
    !---------------------------------------------------------------------------
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncZ
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxdegree
    !---------------------------------------------------------------------------
    integer :: pos, degSum, layerStart, blockStart, polyOrd
    !---------------------------------------------------------------------------

    polyOrd = maxdegree + 1

    ! integer divisions are no mistake here.
    degSum = (ansFuncX-1) + (ansFuncY-1) + (ansFuncZ-1)
    layerStart = ((degSum) * (degSum+1) * (degSum+2)) / 6 + 1
    blockStart = (ansFuncZ-1) * (degSum+1) - ((ansFuncZ-2) * (ansFuncZ-1)) / 2
    pos = layerStart + blockStart + (ansFuncY-1)

  end function posOfModgCoeffPTens


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


  !> Return the number of degrees of freedom for broken polynomial space
  pure function getDofsPTens(maxPolyDegree) result(dofs)
    !---------------------------------------------------------------------------
    !> The maximal polynomial degree per spatial direction (for P Tensor product
    !! polynomials this assumed to be the same for each spatial direction).
    integer, intent(in) :: maxPolyDegree
    !> The number of degrees of freedom for a P tensor product polynomial
    integer :: dofs
    !---------------------------------------------------------------------------

    dofs = ((maxPolyDegree+1)*(maxPolyDegree+2)*(maxPolyDegree+3))/6

  end function getDofsPTens


  !> Return the position of a given ansatz function combination in the
  !! linearized list of modal coefficients for Q-Tensor product polynomials.
  pure function posOfModgCoeffQTens2D(ansFuncX, ansFuncY, ansFuncZ, &
    &                                 maxDegree) result(pos)
    !---------------------------------------------------------------------------
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncZ
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxDegree
    !> The position of the modal coefficient in the list of modal coefficients.
    integer :: pos
    !---------------------------------------------------------------------------
    integer :: polyOrd

    polyOrd = maxDegree+1
    pos = ansFuncX + (ansFuncY-1)*polyOrd

  end function posOfModgCoeffQTens2D


  !> The x, y and z ansatz degrees are turned into the degrees of the next
  !! ansatz function in the linearized Q tensor
  pure subroutine nextModgCoeffQTens2D(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
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
    else
      ! next y index
      ansFuncX = 1
      ansFuncY = ansFuncY+1
    end if
  end subroutine nextModgCoeffQTens2D


  !> Return the number of degrees of freedom for Q polynomial space
  pure function getDofsQTens2D(maxPolyDegree) result(dofs)
    !---------------------------------------------------------------------------
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxPolyDegree
    !> The number of degrees of freedom for a Q tensor product polynomial
    integer :: dofs
    !---------------------------------------------------------------------------

    dofs = (maxPolyDegree+1)**2

  end function getDofsQTens2D


  !> Return the position of a given ansatz function combination in the
  !! linearized list of modal coefficients for P-Tensor product polynomials.
  pure function posOfModgCoeffPTens2D(ansFuncX, ansFuncY, ansFuncZ, maxdegree) &
    & result(pos)
    !---------------------------------------------------------------------------
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncY
    !> Ansatz function index in Z direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncZ
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxdegree
    !---------------------------------------------------------------------------
    integer :: pos, degSum, layerStart, polyOrd
    !---------------------------------------------------------------------------

    polyOrd = maxdegree+1

    ! integer divisions are no mistake here.
    degSum = (ansFuncX-1) + (ansFuncY-1)
    layerStart = ((degSum) * (degSum+1)) / 2 + 1
    pos = layerStart + (ansFuncY-1)

  end function posOfModgCoeffPTens2D


  !> The x, y and z ansatz degrees are turned into the degrees of the next
  !! ansatz function in the layered P list
  pure subroutine nextModgCoeffPTens2D(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncZ
    ! Maximal polynomial degree
    integer, intent(in) :: maxDegree

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


  !> Return the number of degrees of freedom for broken polynomial space
  pure function getDofsPTens2D(maxPolyDegree) result(dofs)
    !---------------------------------------------------------------------------
    !> The maximal polynomial degree per spatial direction (for P Tensor product
    !! polynomials this assumed to be the same for each spatial direction).
    integer, intent(in) :: maxPolyDegree
    !> The number of degrees of freedom for a P tensor product polynomial
    integer :: dofs
    !---------------------------------------------------------------------------

    dofs = (maxPolyDegree+1)*(maxPolyDegree+2)/2

  end function getDofsPTens2D

  !> Return the position of a given ansatz function combination in the
  !! linearized list of modal coefficients for Q-Tensor product polynomials.
  pure function posOfModgCoeffQTens1D(ansFuncX, ansFuncY, ansFuncZ, &
    &                                 maxDegree) result(pos)
    !---------------------------------------------------------------------------
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncZ
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxDegree
    !> The position of the modal coefficient in the list of modal coefficients.
    integer :: pos
    !---------------------------------------------------------------------------
    
    pos = ansFuncX 

  end function posOfModgCoeffQTens1D


  !> The x, y and z ansatz degrees are turned into the degrees of the next
  !! ansatz function in the linearized Q tensor
  pure subroutine nextModgCoeffQTens1D(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncZ
    !> Maximal polynomial degree
    integer, intent(in) :: maxdegree

    ansFuncX = ansFuncX +1

  end subroutine nextModgCoeffQTens1D


  !> Return the number of degrees of freedom for Q polynomial space
  pure function getDofsQTens1D(maxPolyDegree) result(dofs)
    !---------------------------------------------------------------------------
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxPolyDegree
    !> The number of degrees of freedom for a Q tensor product polynomial
    integer :: dofs
    !---------------------------------------------------------------------------

    dofs = (maxPolyDegree+1)

  end function getDofsQTens1D


  !> Return the position of a given ansatz function combination in the
  !! linearized list of modal coefficients for P-Tensor product polynomials.
  pure function posOfModgCoeffPTens1D(ansFuncX, ansFuncY, ansFuncZ, maxdegree) &
    & result(pos)
    !---------------------------------------------------------------------------
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncY
    !> Ansatz function index in Z direction. First ansatz function has index 1.
    integer, intent(in) :: ansFuncZ
    !> The maximal polynomial degree per spatial direction
    integer, intent(in) :: maxdegree
    !> The position
    integer :: pos
    !---------------------------------------------------------------------------

    pos = ansFuncX

  end function posOfModgCoeffPTens1D


  !> The x, y and z ansatz degrees are turned into the degrees of the next
  !! ansatz function in the layered P list
  pure subroutine nextModgCoeffPTens1D(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
    !> Ansatz function index in x direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncX
    !> Ansatz function index in y direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncY
    !> Ansatz function index in z direction. First ansatz function has index 1.
    integer, intent(inout) :: ansFuncZ
    ! Maximal polynomial degree
    integer, intent(in) :: maxDegree

    ansFuncX = ansFuncX +1

  end subroutine nextModgCoeffPTens1D


  !> Return the number of degrees of freedom for broken polynomial space
  pure function getDofsPTens1D(maxPolyDegree) result(dofs)
    !---------------------------------------------------------------------------
    !> The maximal polynomial degree per spatial direction (for P Tensor product
    !! polynomials this assumed to be the same for each spatial direction).
    integer, intent(in) :: maxPolyDegree
    !> The number of degrees of freedom for a P tensor product polynomial
    integer :: dofs
    !---------------------------------------------------------------------------

    dofs = (maxPolyDegree+1)

  end function getDofsPTens1D

end module ply_dof_module
