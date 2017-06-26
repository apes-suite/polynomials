!> This module provides the methods to project the polynomial representation in
!! elements onto the representations in their halves in each dimension.
!!
!! To perform the projection for Legendre polynomials we will use the computed
!! coefficients for the Clenshaw algorithm from [[ply_split_legendre_module]].
!! With those the transformation is just a simple triangular matrix
!! multiplication, but we need to take care of the orthogonal degrees of freedom
!! as we want to handle all of them at the same time.
!! Further we want to allow the transformation to be performed for multiple
!! elements at once.
!!
!! In each dimension we need to perform the following coordinate transformation:
!!
!! \[ x = 0.5 \cdot \xi_{left} - 0.5 \]
!! \[ x = 0.5 \cdot \xi_{right} + 0.5 \]
!!
!! Where \(x\) refers to the coordinate in the original (coarse) element, and
!! \(\xi\) to the coordinates in the two (left and right) halves of the element.
module ply_split_element_module
  use env_module, only: rk
  use ply_split_legendre_module, only: ply_split_legendre_matrix

  implicit none

  private

  public :: ply_split_element_1D
  public :: ply_split_element_init

  !> Precomputed matrix to hold the transformation operation to project
  !! Legendre polynomials to its two half intervals.
  !!
  !! This is computed by [[ply_split_legendre_matrix]], see there for details.
  !! There are two triangular matrices stored in this array, one for the
  !! projection to the left half (-1,0) , and one for the projection to the
  !! right half (0,1).
  !!
  !! This is a module variable, as it is only needed to be computed once with
  !! sufficient size. All lower orders are just subarrays out of the larger one.
  real(kind=rk), allocatable :: split_legendre(:,:)


contains


  ! ------------------------------------------------------------------------ !
  !> Initialization of the module.
  !!
  !! This needs to be performed before any call of the actual transformation
  !! [[ply_split_element_1D]].
  !!
  !! The initialization will compute the transformation matrix for Legendre
  !! polynomials with at least nMaxModes. If the initialization was already
  !! called before with the same or larger nMaxModes, the matrix will not be
  !! changed. Thus, calling this routine will only increase the size of the
  !! module variable split_legendre, never decrease it.
  subroutine ply_split_element_init(nMaxModes)
    ! -------------------------------------------------------------------- !
    !> Maximal number of expected modes to perform the splitting for.
    integer, intent(in) :: nMaxModes
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    if (allocated(split_legendre)) then
      if (size(split_legendre, 1) < nMaxModes) deallocate(split_legendre)
    end if

    if (.not. allocated(split_legendre)) then
      allocate(split_legendre(nMaxModes, nMaxModes))
      split_legendre = ply_split_legendre_matrix(nMaxModes)
    end if

  end subroutine ply_split_element_init
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Project a polynomial representation in elements in one dimension.
  !!
  !! For each parent element the projection on the two respective child elements
  !! (half intervals) are computed for one dimension.
  subroutine ply_split_element_1D(parent_data, child_data)
    ! -------------------------------------------------------------------- !
    !> Polynomial representation in the parent elements.
    !!
    !! The first index are the degrees of freedom in elements, the second index
    !! are the elements.
    real(kind=rk), intent(in) :: parent_data(:,:)

    !> Computed projection of the polynomial representation in the child
    !! elements.
    !!
    !! Again, the first index refers to the degrees of freedom, while the
    !! second index are the elements. There need to be twice as many elements
    !! as in the parent_data.
    real(kind=rk), intent(out) :: child_data(:,:)
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    ! Use split_legendre to compute the two child_data elements for each
    ! parent_data element.
    ! We store the left childs in iChild = (iParent*2 - 1), and the right
    ! childs in iParent*2.

  end subroutine ply_split_element_1D
  ! ======================================================================== !

end module ply_split_element_module
