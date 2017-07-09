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
  use ply_modg_basis_module, only: legendre_1d

  implicit none

  private

  public :: ply_split_element_1D
  public :: ply_split_element_init
  public :: ply_split_element_test

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
  !> Project a polynomial representation in elements in one dimension to its
  !! two halves in that direction.
  !!
  !! For each parent element the projection on the two respective child elements
  !! (half intervals) are computed for one dimension.
  !!
  !>@note Preliminary data layout and interface planning.
  !! It might be that we should rather split the index into the direction
  !! in which we perform the operation and all the other directions normal
  !! to that. For a dense matrix this may allow the compiler to detect
  !! the matrix multiply. However, here for the triangular matrix it is not
  !! so sure, whether this would be possible.
  !!@endnote
  !!
  !! As we need to perform this operation in all dimensions, it would be good
  !! to shift the indices around. When doing this, we can stick to the same
  !! implementation for all directions, without the need to put any logic in
  !! here to decide on the current direction.
  !! In 3D we would end up with this chain:
  !! (x,y,z) -> split_element for Z -> (z,x,y)
  !!         -> split_element for Y -> (y,z,x)
  !!         -> split_element for X -> (x,y,z)
  !! Thus, the logic is that we perform the split on the last dimension, and
  !! cycle the indices in the output.
  !!
  !! We can generalize this to arbitrary dimensions.
  !! In 2D it would look like this:
  !! (x,y) -> split_element for Y -> (y,x)
  !!       -> split_element for X -> (x,y)
  !! And in 1D, we just need to perform one transformation:
  !! (x) -> split_element for X -> (x)
  !!
  !! As we allow for a changed number of polynomial degrees in the input and
  !! output, we need to take care of different lengths for each direction.
  !! Thus, we need: the dimensionality, two 1D arrays with the length of this
  !! dimensionality to provide the number of degrees of freedom for each
  !! direction (once for the input, and once for the output).
  !!
  !! We need: nDofs in the direction where the transformation is to be done
  !!          and the nDofs for all normal directions.
  subroutine ply_split_element_1D( nDims, inLen, outLen, parent_data, &
    &                              child_data                         )
    ! -------------------------------------------------------------------- !
    !> Number of dimensions of the polynomial data.
    integer, intent(in) :: nDims

    !> Number degrees of freedom for each direction in parent_Data.
    !!
    !! The first index of parent_data needs to have a length equal to the
    !! product of all inLen components.
    !! The splitting operation will be done in the last dimension.
    integer, intent(in) :: inLen(nDims)

    !> Number degrees of freedom for each direction in child_Data.
    !!
    !! The first index of child_data needs to have a length equal to the
    !! product of all outLen components.
    !! The data will be cyclicly exchanged. Thus, the last dimension in
    !! parent_data corresponds to the first in one in child_data and all
    !! other components are shifted once to the right.
    integer, intent(in) :: outLen(nDims)

    !> Polynomial representation in the parent elements.
    !!
    !! The first index are the degrees of freedom in elements, the second index
    !! are the elements.
    !! In the first index the shape of data has to be in the form
    !! (inLen(1), inLen(2), ... , inLen(nDims)).
    !! The splitting operation is performed on the last dimension in that
    !! data.
    real(kind=rk), intent(in) :: parent_data(:,:)

    !> Computed projection of the polynomial representation in the child
    !! elements.
    !!
    !! Again, the first index refers to the degrees of freedom, while the
    !! second index are the elements. There need to be twice as many elements
    !! as in the parent_data.
    !! Left childs are stored in iChild = (iParent*2 - 1), and the right
    !! childs in iParent*2.
    !!
    !! In the first index the shape of the data has to be in the form
    !! (outLen(1), outLen(2), ... , outLen(nDims)), the data is rotated
    !! in comparison to parent_data and the splitted direction has to be
    !! the first one in child_data (while it was the last in parent_data),
    !! and all other dimensions are shifted by one to the right.
    real(kind=rk), intent(out) :: child_data(:,:)
    ! -------------------------------------------------------------------- !
    integer :: iDir
    integer :: iParent, Lchild, Rchild
    integer :: parentMode, childMode
    integer :: maxrow
    integer :: indep
    integer :: nIndeps
    integer :: nParents
    integer :: parentpos, childpos
    ! -------------------------------------------------------------------- !

    nParents = size(parent_data,2)

    ! Use split_legendre to compute the two child_data elements for each
    ! parent_data element.
    ! We store the left childs in iChild = (iParent*2 - 1), and the right
    ! childs in iParent*2.

    child_data = 0.0_rk

    ! The number of independent modes (in normal directions) is given
    ! by the product of the length in all directions, except the last one.
    nIndeps = 1
    do iDir=1,nDims-1
      nIndeps = nIndeps*inLen(iDir)
    end do

    oldmodes: do parentMode=1,inLen(nDims)
      ! Maximal number modes to compute, as this is a triangular matrix
      ! it is limited by the diagonal (parentMode). However, it may be
      ! that the targe polynomial space in the output is smaller, in this
      ! case we cap the computations and no more than outLen(1) entries
      ! are to be computed.
      maxrow = min(parentMode, outLen(1))

      elemloop: do iParent=1,nParents
        Rchild = iParent*2
        Lchild = Rchild - 1
        newmodes: do childMode=1,maxrow

          do indep=1,nIndeps
            parentpos = indep + nIndeps*(parentMode-1)
            childpos = childmode + (indep-1)*outLen(1)
            child_data(childpos, Lchild) = child_data(childpos, Lchild) &
              &                          + split_legendre( parentmode,  &
              &                                            childmode )  &
              &                            * parent_data(parentpos, iParent)
            child_data(childpos, Rchild) = child_data(childpos, Rchild) &
              &                          + split_legendre( childmode,   &
              &                                            parentmode ) &
              &                            * parent_data(parentpos, iParent)
          end do

        end do newmodes
      end do elemloop

    end do oldmodes

  end subroutine ply_split_element_1D
  ! ======================================================================== !


  ! !!!!!!! !
  ! testing !
  ! !!!!!!! !

  ! ------------------------------------------------------------------------ !
  !> Testing the 1D splitting.
  !!
  subroutine ply_split_element_1D_test(nModes, success)
    ! -------------------------------------------------------------------- !
    !> Number of modes in the (1D) polynomials to use in the check.
    integer, intent(in) :: nModes

    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    integer :: parentModes, childmodes
    integer :: iPoint
    integer :: iElem
    real(kind=rk) :: xi(nModes)
    real(kind=rk) :: x_left(nModes)
    real(kind=rk) :: x_right(nModes)
    real(kind=rk) :: legchild(nModes, nModes)
    real(kind=rk) :: legleft(nModes, nModes)
    real(kind=rk) :: legright(nModes, nModes)
    real(kind=rk) :: rootval(nModes, 2)
    real(kind=rk) :: childval
    real(kind=rk), allocatable :: rootelem(:,:)
    real(kind=rk), allocatable :: childelem(:,:)
    real(kind=rk) :: tolerance
    ! -------------------------------------------------------------------- !

    call ply_split_element_init(nModes)

    tolerance = 8*epsilon(1.0_rk)*nmodes**2
    success = .true.

    ! Some random points to check the resulting child polynomials.
    call random_number(xi)

    legchild = legendre_1D(xi, nModes-1)

    ! The corresponding positions in the left and right half of the root
    ! element.
    x_right = 0.5_rk*xi + 0.5_rk
    x_left  = 0.5_rk*xi - 0.5_rk

    legleft = legendre_1D(x_left, nModes-1)
    legright = legendre_1D(x_right, nModes-1)

    do parentmodes=1,nModes
      allocate(rootelem(parentModes,1))
      call random_number(rootelem)
      do iPoint=1,nModes
        rootval(iPoint,1) = sum( rootelem(:,1)                  &
          &                    * legleft(:parentModes,iPoint) )
        rootval(iPoint,2) = sum( rootelem(:,1)                   &
          &                     * legright(:parentModes,iPoint) )
      end do
      do childmodes=1,parentModes-1
        allocate(childelem(childmodes,2))
        call ply_split_element_1D( nDims       = 1,             &
          &                        inLen       = [parentModes], &
          &                        outLen      = [childModes],  &
          &                        parent_data = rootelem,      &
          &                        child_data  = childelem      )
        success = success                                          &
          &       .and. ( 0.5_rk*(childelem(1,1) + childelem(1,2)) &
          &               - rootelem(1,1) < tolerance              )
        deallocate(childelem)
      end do
      do childmodes=parentmodes,nModes
        allocate(childelem(childmodes,2))
        call ply_split_element_1D( nDims       = 1,             &
          &                        inLen       = [parentModes], &
          &                        outLen      = [childModes],  &
          &                        parent_data = rootelem,      &
          &                        child_data  = childelem      )
        do iElem=1,2
          do iPoint=1,nModes
            childval = sum( childelem(:,iElem)             &
                            * legchild(:childmodes,iPoint) )
            success = success &
              &       .and. ( abs(rootval(iPoint,iElem) - childval) &
              &               < tolerance                           )
          end do
        end do
        deallocate(childelem)
      end do
      deallocate(rootelem)
    end do

  end subroutine ply_split_element_1D_test
  ! ======================================================================== !


  ! ------------------------------------------------------------------------ !
  !> Testing routine for the functions of this module.
  subroutine ply_split_element_test(success)
    ! -------------------------------------------------------------------- !
    !> Indication whether the tests were completed successfully.
    logical, intent(out) :: success
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    call ply_split_element_init(80)

    call ply_split_element_1D_test(nModes = 30, success = success)

    if (.not. success) then
      write(*,*) 'Check for 1D splitting FAILED!'
      RETURN
    end if

  end subroutine ply_split_element_test
  ! ======================================================================== !

end module ply_split_element_module
