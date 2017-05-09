!! This module is used for the projection of QLegendre polynomials from parent
!! to child elements.

module ply_poly_transformation_module
  use env_module,                   only: rk
  use tem_aux_module,               only: tem_abort
  use treelmesh_module,             only: treelmesh_type
  use tem_param_module,             only: childPosition
  use ply_LegPolyProjection_module, only: ply_subsample_type, &
    &                                     ply_array_type 

  implicit none

  private

  public :: ply_Poly_Transformation

contains

  ! ************************************************************************ !
  !> Projection of polynomial data from parent elements to child elements.
  !! The projection is done by a direct transformation of the modal
  !! coeffiecients to another coordinate system with z=ax+b.
  subroutine ply_Poly_Transformation( subsamp, dofReduction, mesh, meshData,   &
    &                                 varDofs, varComps, ndims, refine_tree,   &
    &                                 new_refine_tree, newMeshData, newVarDofs )
    ! -------------------------------------------------------------------- !
    !> Parameters for the subsampling
    type(ply_subsample_type), intent(in) :: subsamp

    !> Factor for reduction of degrees of freedom.
    real(kind=rk), intent(in) :: dofReduction(:)

    !> The mesh related to meshData.
    type(treelmesh_type), intent(in) :: mesh

    !> The data for subsampling.
    type(ply_array_type), intent(in) :: meshData(:)

    !> The number of degrees of freedom for every variable.
    integer, intent(in) :: varDofs(:)

    !> The number of components for every variable.
    integer, intent(in) :: varComps(:)

    !> Number of dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    !> Logical array that marks elements for refinement
    !! of the previous sampling level.
    logical, intent(in) :: refine_tree(:)

    !> Logical array that marks elements for refinement.
    logical, intent(in) :: new_refine_tree(:)

    !> The subsampled data for new_refine_tree.
    type(ply_array_type), allocatable, intent(out) :: newMeshData(:)

    !> The number of dofs for the subsampled data.
    integer, allocatable, intent(out) :: newVarDofs(:)
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: workData(:)
    real(kind=rk), allocatable :: newWorkData(:)
    real(kind=rk), allocatable :: transform_matrix(:,:) 
    integer :: nVars, nDofs, nComponents, nChildDofs, nChilds
    integer :: iVar, max_modes
    ! -------------------------------------------------------------------- !
    nVars = size(varDofs)
    allocate(newVarDofs(nVars))
    allocate(newMeshData(nVars))

    varLoop: do iVar = 1, nVars
      nDofs = varDofs(iVar)
      nComponents = varComps(iVar)

      allocate(workData(size(meshData(iVar)%dat)))
      workData = meshData(iVar)%dat

      nChildDofs = (ceiling(nint(nDofs**(1.0_rk/real(ndims, kind=rk))) &
        &                        * dofReduction(iVar)))**ndims

      if (nChildDofs > nDofs) then
        nChildDofs = nDofs
      elseif (nChildDofs < 1) then
        nChildDofs = 1
      end if

      if(subsamp%sampling_lvl == subsamp%maxsub) then
        nChildDofs = 1
      end if

      max_modes = nint(real(nDofs, kind=rk)**(1.0_rk/real(nDims, kind=rk)))

      ! Get the transformation matrix for the maximal polynomial degree.
      allocate(transform_matrix(max_modes, max_modes))
      call ply_transform_matrix( max_modes = max_modes,       &
        &                        v         = transform_matrix )
!!write(*,*) 'max_modes',max_modes
!!write(*,*) 'transform_matrix',transform_matrix

      call ply_subsampleData( mesh             = mesh,             &
        &                     meshData         = workData,         &
        &                     nDofs            = nDofs,            &
        &                     nChildDofs       = nChildDofs,       &
        &                     nComponents      = nComponents,      &
        &                     transform_matrix = transform_matrix, &
        &                     refine_tree      = refine_tree,      &
        &                     new_refine_tree  = new_refine_tree,  &
        &                     ndims            = ndims,            &
        &                     subsamp          = subsamp,          &
        &                     newMeshData      = newWorkData       )

      allocate(newMeshData(iVar)%dat(size(newWorkData)))
      newMeshData(iVar)%dat = newWorkData

      newVarDofs(iVar) = nChildDofs

      deallocate(workData)
      deallocate(transform_matrix)

    end do varLoop
 
  end subroutine ply_Poly_Transformation
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Coefficients from the recursive formulation of legendre polynomials.
  !! L_n = alpha * x * L_n-1 + beta * L_n-2
  !!
  subroutine ply_transform_matrix(max_modes, v)
    ! -------------------------------------------------------------------- !
    !> The number of modes in a single spatial direction.
    !!
    integer, intent(in) :: max_modes

    !> The transformation matrix.
    !!
    !! Upper triangular matrix is created for shifting and lower triangular
    !! for (-1) * shifting.
    real(kind=rk), allocatable, intent(out) :: v(:,:)
    ! -------------------------------------------------------------------- !
    integer :: m, orig
    real(kind=rk) :: shifting, scaling
    ! -------------------------------------------------------------------- !
    allocate(v(max_modes,max_modes))
    v(:,:) = 0.0

    scaling = 0.5
    shifting = 0.5

    !! Set the first entries of v manually.
    v(1,1) = 1.0
    if (max_modes > 1) then
      v(1,2) = shifting
      v(2,2) = scaling

      if (max_modes > 2) then
        do orig = 3,max_modes
          v(1,orig) = ply_beta(orig-1) * v(1,orig-2)                    &
            &       + ply_alpha(orig-1) * shifting * v(1,orig-1)        &
            &       - scaling * ply_alpha_beta(2,orig-1) * v(2, orig-1)
          do m = 2,orig
            if (m < max_modes) then
              v(m,orig) = ply_beta(orig-1) * v(m,orig-2)                        &
                &       + ply_alpha(orig-1) * shifting * v(m,orig-1)            &
                &       - scaling * ply_alpha_beta(m+1,orig-1) * v(m+1, orig-1) & 
                &       + scaling * ply_alpha_frac(m-1,orig-1) * v(m-1,orig-1)  
            else
              !! Need to skip one summand for v(max_modes,max_modes).
              v(m,orig) = scaling * ply_alpha_frac(m-1,orig-1) * v(m-1,orig-1)
            end if
          end do
        end do
      end if
      !! Fill the lower triangular matrix with help of the entries in upper
      !! triangular matrix.
      do m = 1 , max_modes 
        do orig = 1, m-1
          !! 
          if (mod((m+orig),2) /= 0) then
            v(m ,orig) = - v(orig,m)
          else
            v(m ,orig) = v(orig,m)
          end if
        end do
      end do
    end if

  end subroutine ply_transform_matrix
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Coefficients from the recursive formulation of legendre polynomials.
  !! L_n = alpha * x * L_n-1 + beta * L_n-2
  !!
  subroutine ply_subsampleData(mesh, meshData, nDofs, nChildDofs,   &
    &                          nComponents, transform_matrix,       &
    &                          refine_tree, new_refine_tree, nDims, &
    &                          subsamp, newMeshData                 )
    ! -------------------------------------------------------------------- !
    !> The mesh for the data.
    type(treelmesh_type), intent(in) :: mesh

    !> The data to subsample
    real(kind=rk), intent(in) :: meshData(:)

    !> The number of degrees of freedom.
    integer, intent(in) :: nDofs

    !> The number of degrees of freedom for the child elements.
    integer, intent(in) :: nChildDofs

    !> Number of Components.
    integer, intent(in) :: nComponents

    !> The transformation matrix for the linear transformation z=ax+b
    real(kind=rk), allocatable, intent(in) :: transform_matrix(:,:)

    !> Logical array that marks all elements for refinement for the previous
    !! sampling level.
    logical, intent(in) :: refine_tree(:)

    !> Logical array that marks all elements for refinement for the current
    !! sampling level.
    logical, intent(in) :: new_refine_tree(:)

    !> The number of dimensions in the polynomial representation.
    integer, intent(in) :: nDims

    !> Parameters for subsampling.
    type(ply_subsample_type), intent(in) :: subsamp

    !> The subsampled Data.
    real(kind=rk), allocatable, intent(out) :: newMeshData(:)
    ! -------------------------------------------------------------------- !
    integer :: nChilds, nElems, nElemsToRefine, nElemsNotToRefine
    integer :: nParentElems, childPos, lowElemIndex, upElemIndex
    integer :: iParentElem, iChild, iElem, upChildIndex, lowChildIndex
    integer :: oneDof, noChilds
    real(kind=rk), allocatable :: childData(:)
    ! -------------------------------------------------------------------- !
    nChilds = 2**nDims
    nElems = mesh%nElems
    nElemsToRefine = count(new_refine_tree)
    nElemsNotToRefine = nElems - nElemsToRefine
    nParentElems = size(refine_tree)

    ! Now, we set the correct data for the newMeshData.
    allocate(newMeshData((nElemsToRefine * nChilds * nChildDofs &
      &                  + nElemsNotToRefine) * nComponents))

    newMeshData = 0.0_rk
    upChildIndex = 0
    upElemIndex = 0
    childPos = 0
    
    if (subsamp%sampling_lvl > 1) then

      elementLoop: do iParentElem=1,nParentElems
        ! Check if the parent cell was already refined...
        if (refine_tree(iParentElem)) then
          ! Parent cell was already refined so it contains data with nDofs
          childLoop: do iChild=1,nChilds
            childPos = childPos + 1
            ! Check if the child elems will be refined...
            if (new_refine_tree(childPos)) then
             
              ! Child cell will be refined.
              !
              ! Need to project current elem data to new childs 
              ! with reduced dofs.
              ! Create lower and upper indices for all data of 
              ! iElem in meshData.
              lowElemIndex = upElemIndex + 1 
              upElemIndex = (lowElemIndex-1) + nDofs * nComponents

              ! Project these dofs from the coarse element to the 
              ! finer elements.
              call ply_projDataToChild(                                   &
                &  parentData       = meshData(lowElemIndex:upElemIndex), &
                &  nParentDofs      = nDofs,                              &
                &  nChildDofs       = nChildDofs,                         &
                &  nComponents      = nComponents,                        &
                &  nDimensions      = nDims,                              &
                &  nChilds          = nChilds,                            &
                &  transform_matrix = transform_matrix,                   &
                &  childData        = childData                           )
  
              ! Set the data correctly in newMeshData.
              lowChildIndex = upChildIndex + 1
              upChildIndex = (lowChildIndex-1) + &
                &            nChilds * nChildDofs * nComponents

              newMeshData(lowChildIndex:upChildIndex) = childData
              deallocate(childData)
            else
              ! Child cell won't be refined.
              !
              ! Need projecton from current dofs to 1 dof.
              ! Create lower and upper indices for all data of 
              ! iElem in meshData.
              lowElemIndex = upElemIndex + 1
              upElemIndex = (lowElemIndex-1) + nDofs * nComponents

              ! Projection from nDofs to oneDof (integral mean valuea).
              oneDof = 1
              noChilds = 1
              call ply_projDataToChild(                                   &
                &  parentData       = meshData(lowElemIndex:upElemIndex), &
                &  nParentDofs      = nDofs,                              &
                &  nChildDofs       = oneDof,                             &
                &  nComponents      = nComponents,                        &
                &  nDimensions      = nDims,                              &
                &  nChilds          = noChilds,                           &
                &  transform_matrix = transform_matrix,                   &
                &  childData        = childData                           )
  
              ! Iterate over all childDofs and set the data corectly 
              ! in newMeshData.
              lowChildIndex = upChildIndex + 1
              upChildIndex = (lowChildIndex-1) + nComponents
      
              newMeshData(lowChildIndex:upChildIndex) = childData
              deallocate(childData)
            end if
          end do childLoop

        else
          ! Parent cell wasn't refined so it contains data with only one dof.
          ! Simple copying.
          allocate(childData(nComponents))

          lowElemIndex = upElemIndex + 1
          upElemIndex = (lowElemIndex-1) + nComponents

          childData = meshData(lowElemIndex:upElemIndex)

          lowChildIndex = upChildIndex + 1
          upChildIndex = (lowChildIndex-1) + nComponents

          newMeshData(lowChildIndex:upChildIndex) = childData
          childpos = childpos + 1
          deallocate(childData)

        end if

      end do elementLoop

    else

      elemLoop: do iElem=1,nElems
        if (new_refine_tree(iElem)) then
          ! Create lower and upper indices for all data of iElem in meshData.
          lowElemIndex = upElemIndex + 1 
          upElemIndex = (lowElemIndex-1) + nDofs * nComponents

          ! Project these dofs from the coarse element to the 
          ! finer elements.
          call ply_projDataToChild(                                   &
            &  parentData       = meshData(lowElemIndex:upElemIndex), &
            &  nParentDofs      = nDofs,                              &
            &  nChildDofs       = nChildDofs,                         &
            &  nComponents      = nComponents,                        &
            &  nDimensions      = nDims,                              &
            &  nChilds          = nChilds,                            &
            &  transform_matrix = transform_matrix,                   &
            &  childData        = childData                           )
  
          ! Iterate over all childDofs and set the data corectly in newMeshData
          lowChildIndex = upChildIndex + 1
          upChildIndex = (lowChildIndex-1) + nChilds * nChildDofs * nComponents

          newMeshData(lowChildIndex:upChildIndex) = childData
          deallocate(childData)
        else
          ! Create lower and upper indices for all data of iElem in meshData.
          lowElemIndex = upElemIndex + 1
          upElemIndex = (lowElemIndex-1) + nDofs * nComponents

          ! Projection from nDofs to oneDof (integral mean value).
          oneDof = 1
          noChilds = 1
          call ply_projDataToChild(                                   &
            &  parentData       = meshData(lowElemIndex:upElemIndex), &
            &  nParentDofs      = nDofs,                              &
            &  nChildDofs       = oneDof,                             &
            &  nComponents      = nComponents,                        &
            &  nDimensions      = nDims,                              &
            &  nChilds          = noChilds,                           &
            &  transform_matrix = transform_matrix,                   &
            &  childData        = childData                           )
  
          ! Iterate over all childDofs and set the data corectly in newMeshData
          lowChildIndex = upChildIndex + 1
          upChildIndex = (lowChildIndex-1) + nComponents

          newMeshData(lowChildIndex:upChildIndex) = childData
          deallocate(childData)
        end if
      end do elemLoop
    end if


  end subroutine ply_subsampleData
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Subroutine to project element data from a parent cell to its children.
  subroutine ply_projDataToChild( parentData, nParentDofs, nChildDofs,        &
    &                             nComponents, nDimensions, nChilds,          &
    &                             transform_matrix, childData                 )
    ! -------------------------------------------------------------------- !
    !> The polynomial data for a single parent element.
    real(kind=rk), intent(in) :: parentData(:)

    !> The number of dofs of the parent element.
    integer, intent(in) :: nParentDofs

    !> The total number of dofs for the child cells.
    integer, intent(in) :: nChildDofs

    !> The number of componentns for the given variable.
    integer, intent(in) :: nComponents

    !> The number of dimensions.
    integer, intent(in) :: nDimensions

    !> The number of child elements.
    integer, intent(in) :: nChilds

    !> The transformation matrix for the linear coordinate transformation.
    real(kind=rk), intent(in) :: transform_matrix(:,:)

    !> The new data representation for all child cell of the parent cell.
    real(kind=rk), allocatable, intent(out) :: childData(:)
    ! -------------------------------------------------------------------- !
    integer :: iDimension, iSubElem, iIndep, iComponent, iMode, jMode, i, j
    integer :: nSubElems, nIndeps, lMode, kMode, iChildElem_prev
    integer :: nChildElems_prev, nChildElems_cur
    integer :: child_dofPos, dofPos, childElem
    integer :: parent_modes, child_modes
    integer :: lower_bound, upper_bound, stride
    real(kind=rk), allocatable :: temp_data(:)
    real(kind=rk), allocatable :: temp_childData_prev(:)
    real(kind=rk), allocatable :: temp_childData(:)
    ! -------------------------------------------------------------------- !
    parent_modes = nint(real(nParentDofs,kind=rk)        &
      &                 **(1/real(nDimensions,kind=rk)))

    child_modes = nint(real(nChildDofs,kind=rk)         &
      &                **(1/real(nDimensions,kind=rk)))
!!write(*,*) 'nChildDofs=',nChildDofs
!!write(*,*) 'child_modes=',child_modes
write(*,*) '///////////////////////////////////////////////////////////////////'
write(*,*) 'nChilds=',nChilds
parent_modes = 4
child_modes = 1

    ! transformation matrix looks like this:
    ! [1.0  --  --     shift=0.5   ]
    ! | |  0.5  --                 |
    ! | |   |  0.25             ...|
    ! |             0.125          |
    ! | shift=-0.5        0.0625   |
    ! [    :                    ...]

    allocate(temp_Data(child_modes))
    
    do iDimension = 1,nDimensions
write(*,*) ''
write(*,*) '-------------------------------------------------------------------'
write(*,*) 'iDimension=',iDimension
      if (nChilds > 1) then
        nSubElems = 2
        nChildElems_cur = 2**(iDimension)
        nChildElems_prev = 2**(iDimension-1)
      else
        nChildElems_cur = nChilds
        nChildElems_prev = nChilds
        nSubElems = nChilds
      end if

      stride = child_modes**(iDimension-1) * nComponents
write(*,*) 'stride=',stride

      if (iDimension .eq. 1) then

        nIndeps = parent_modes**(nDimensions-1)

        ! Allocate memory for two childs and set it to zero
        if (iDimension .eq. nDimensions) then
          allocate(childData(nChildElems_cur * parent_modes**(nDimensions-1) * &
            &                child_modes**(nDimensions-2) * nComponents))
          childData = 0.0_rk
        else
          allocate(temp_childData(nChildElems_cur * parent_modes**(nDimensions-1) * &
            &                     child_modes**(nDimensions-2) * nComponents))
          temp_childData = 0.0_rk
        end if

        childElem = 0
  
        do iChildElem_prev = 1, nChildElems_prev
write(*,*) '///////////////////////////////'
write(*,*) 'iChildElem_prev',iChildElem_prev
          do iSubElem = 1, nSubElems
            childElem = childElem + 1
write(*,*) '*******************************'
write(*,*) 'ChildElem=',ChildElem
            do iComponent = 1, nComponents
write(*,*) '-------------------------------'
write(*,*) 'iComponent=',iComponent
              do iIndep = 1, nIndeps
write(*,*) '+++++++++++++++++++++++++++++++'
write(*,*) 'iIndep=',iIndep

                lower_bound = iComponent + &
                  &           (iIndep - 1) * parent_modes * nComponents
write(*,*) 'lower_bound=',lower_bound
                upper_bound = iIndep * parent_modes * nComponents
write(*,*) 'upper_bound=',upper_bound

                temp_data(:) = parentData(lower_bound:upper_bound:stride)

                do iMode = 1, child_modes
                  do jMode = iMode, parent_modes
  
                    child_dofPos = iComponent +                               &
                      &            (iMode-1) * nComponents +                  &
                      &            (iIndep - 1) * child_modes * nComponents + &
                      &            (childElem-1) * nComponents * &
                      &            child_modes**(nDimensions-2) * &
                      &            parent_modes**(nDimensions-1) 
write(*,*) 'child_dofPos=',child_dofPos
  
                    if (iSubElem .eq. 1) then
                      ! upper part of transform_matrix
                      j = iMode
                      i = jMode
                    else
                      ! lower part of transform_matrix
                      j = jMode
                      i = iMode
                    end if
write(*,*) 'j=',j
write(*,*) 'i=',i
                    if (iDimension .eq. nDimensions) then
                      childData(child_dofpos) = childData(child_dofpos) + &
                        &                       temp_Data(jMode) *        &
                        &                       transform_matrix(j,i)
                    else
                      temp_childData(child_dofpos) = temp_childData        &
                        &                            (child_dofpos) +      &
                        &                            temp_Data(jMode) *    &
                        &                            transform_matrix(j,i)
                    end if
  
                  end do 
                end do
              end do
            end do
          end do
        end do
!!call tem_abort()

      elseif (iDimension .eq. 2) then

        ! Copy old temp_childData in temp_childData_prev
        ! This is the data we will work on.
        allocate(temp_childData_prev(size(temp_childData)))
        temp_childData_prev = temp_childData
        deallocate(temp_childData)

        ! Allocate memory for the four childs in y-direction
        if (iDimension .eq. nDimensions) then
          allocate(childData(nChildElems_cur * parent_modes**(nDimensions-2) * &
            &                child_modes**(nDimensions-1) * nComponents) )
          childData = 0.0_rk
        else

          allocate(temp_childData(nChildElems_cur * parent_modes**(nDimensions-2) * &
            &                     child_modes**(nDimensions-1) * nComponents))
          temp_childData = 0.0_rk
        end if

        childElem = 0

        do iChildElem_prev = 1, nChildElems_prev
write(*,*) '///////////////////////////////'
write(*,*) 'iChildElem_prev=',iChildElem_prev
          do iSubElem = 1, nSubElems
            childElem = childElem + 1
write(*,*) '*******************************'
write(*,*) 'ChildElem=',ChildElem
            do iComponent = 1, nComponents
              iIndep = 0
write(*,*) '-------------------------------'
write(*,*) 'iComponent=',iComponent
              do kMode = 1, parent_modes
                do lMode = 1, child_modes
                  iIndep = iIndep + 1
write(*,*) '+++++++++++++++++++++++++++++++'
write(*,*) 'iIndep=',iIndep

                  lower_bound = iComponent + (lMode-1) * nComponents +     &
                    &           (kMode-1) * child_modes*parent_modes * nComponents + &
                    &           (iChildElem_prev-1) * nComponents *               &
                    &           child_modes**(nDimensions-2) * &
                    &           parent_modes**(nDimensions-1)
write(*,*) 'lower_bound=',lower_bound
                  upper_bound = kMode * parent_modes*child_modes * nComponents + &
                    &           (iChildElem_prev-1) * nComponents * &
                    &           child_modes**(nDimensions-2) * &
                    &           parent_modes**(nDimensions-1)
write(*,*) 'upper_bound=',upper_bound

                  temp_data(:) = temp_childData_prev &
                    &            (lower_bound:upper_bound:stride)

                  do iMode = 1, child_modes
                    do jMode = iMode, parent_modes

                      child_dofPos = iComponent +                                &
                        &            (iMode-1) * child_modes * nComponents +     &
                        &            (lMode-1) * nComponents + &
                        &            (kMode-1) * nComponents * child_modes**2 +  &
                        &            (childElem-1) * nComponents * &
                        &            child_modes**(nDimensions-1) * &
                        &            parent_modes**(nDimensions-2)

write(*,*) 'child_dofPos=',child_dofPos
  
                      if (iSubElem .eq. 1) then
                        ! upper part of transform_matrix
                        j = iMode
                        i = jMode
                      else
                        ! lower part of transform_matrix
                        j = jMode
                        i = iMode
                      end if
write(*,*) 'j=',j
write(*,*) 'i=',i

                      if (iDimension .eq. nDimensions) then
                        childData(child_dofpos) = childData(child_dofpos) + &
                          &                       temp_Data(jMode) *        &
                          &                       transform_matrix(j,i)
                      else
                        temp_childData(child_dofpos) = temp_childData        &
                          &                            (child_dofpos) +      &
                          &                            temp_Data(jMode) *    &
                          &                            transform_matrix(j,i)
                      end if
  
                    end do 
                  end do
                end do
              end do
            end do
          end do
        end do
!!call tem_abort()

      else
        deallocate(temp_childData_prev)
        allocate(temp_childData_prev(size(temp_childData)))
        temp_childData_prev = temp_childData

        ! Allocate memory for the eight childs in z-direction
        allocate(childData(nChildElems_cur * child_modes**3 &
          &                 * nComponents))
        childData = 0.0_rk

        childElem = 0

        do iChildElem_prev = 1, nChildElems_prev
write(*,*) '///////////////////////////////'
write(*,*) 'iChildElem_prev=',iChildElem_prev
          do iSubElem = 1, nSubElems
            childElem = childElem + 1
write(*,*) '*******************************'
write(*,*) 'ChildElem=',ChildElem
            do iComponent = 1, nComponents
              iIndep = 0
write(*,*) '*******************************'
write(*,*) 'iComponent=',iComponent
              do kMode = 1,child_modes
                do lMode = 1,child_modes
                  iIndep = iIndep + 1
write(*,*) '-------------------------------'
write(*,*) 'iIndep=',iIndep

                  lower_bound = iComponent + (iIndep - 1) * nComponents + &
                    &           (iChildElem_prev-1) * nComponents * &
                    &           child_modes**(nDimensions-1) * &
                    &           parent_modes**(nDimensions-2)
write(*,*) 'lower_bound=',lower_bound
                  upper_bound = parent_modes * child_modes**(nDimensions-1) * nComponents - &
                    &           (nComponents - iComponent) + &
                    &           (iChildElem_prev-1) * nComponents * &
                    &           child_modes**(nDimensions-1) * &
                    &           parent_modes**(nDimensions-2)
write(*,*) 'upper_bound=',upper_bound

                  temp_data(:) = temp_childData_prev(lower_bound:upper_bound:stride)

                  do iMode = 1, child_modes
                    do jMode = iMode, parent_modes
  
                      child_dofPos = iComponent +                                  &
                        &            (iMode-1) * child_modes**2 * nComponents +    &
                        &            (lMode - 1) * nComponents +                   &
                        &            (kMode - 1) * child_modes * nComponents +  &
                        &            (childElem - 1) * nComponents * &
                        &            child_modes**nDimensions 
write(*,*) 'child_dofPos=',child_dofPos
  
                      if (iSubElem .eq. 1) then
                        ! upper part of transform_matrix
                        j = iMode
                        i = jMode
                      else
                        ! lower part of transform_matrix
                        j = jMode
                        i = iMode
                      end if
write(*,*) 'j=',j
write(*,*) 'i=',i

                      childData(child_dofpos) = childData(child_dofpos) + &
                        &                       temp_Data(jMode) *        &
                        &                       transform_matrix(j,i)

                    end do  
                  end do 
                end do
              end do
            end do
          end do
        end do
      end if
    end do


    deallocate(temp_childData)
    deallocate(temp_childData_prev)
    deallocate(temp_data)
call tem_abort()
    

  end subroutine ply_projDataToChild
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Coefficients from the recursive formulation of legendre polynomials.
  !! L_n = alpha * x * L_n-1 + beta * L_n-2
  function ply_alpha(mode) result(alpha)
    ! -------------------------------------------------------------------- !
    !> The current mode in the polynomial representation.
    integer, intent(in) :: mode

    !> Alpha coefficient from the recursive formulation of legendre
    !! polynomials.
    real(kind=rk) :: alpha
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    if (mode > 0) then
      alpha = real((2 * mode - 1), kind=rk) / real(mode, kind=rk)
    else
      alpha = 0.0
    end if
    return

  end function ply_alpha
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Coefficients from the recursive formulation of legendre polynomials.
  !! L_n = alpha * x * L_n-1 + beta * L_n-2
  !!
  function ply_beta(mode) result(beta)
    ! -------------------------------------------------------------------- !
    !> The current mode in the polynomial representation.
    integer, intent(in) :: mode

    !> Beta coefficient from the recursive formulation of legendre
    !! polynomials.
    real(kind=rk) :: beta
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    if (mode > 0) then
      beta = real((1 - mode), kind=rk) / real(mode, kind=rk)
    else
      beta = 0.0
    end if
    return

  end function ply_beta
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Quotient of two alpha values.
  function ply_alpha_frac(denominator, numerator) result(alpha_frac)
    ! -------------------------------------------------------------------- !
    !> Numerator
    integer, intent(in) :: numerator

    !> Denominator
    integer, intent(in) :: denominator

    !> The quotient of two alpha values.
    real(kind=rk) :: alpha_frac
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    if ( denominator > 0 .and. numerator > 0 ) then
      if ( denominator == numerator ) then
        alpha_frac = 1.0_rk
      else
        alpha_frac = real((2*numerator-1)*denominator, kind=rk)  &
          &            / real((2*denominator-1)*numerator, kind=rk) 
      end if
    else
      alpha_frac = 0.0_rk  
    end if
    return

  end function ply_alpha_frac
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Prodcut of alpha(numerator) * beta(denominator) / alpha(denominator)
  function ply_alpha_beta(denominator, numerator) result(alpha_beta)
    ! -------------------------------------------------------------------- !
    !> Numerator
    integer, intent(in) :: numerator

    !> Denominator
    integer, intent(in) :: denominator

    !> The product of alpha(n) * beta(d) / alpha(d)
    real(kind=rk) :: alpha_beta
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    if (numerator > 0) then
      if ( denominator == numerator ) then
        alpha_beta = real((1-denominator), kind=rk)/ real(denominator, kind=rk)
      else
        alpha_beta = real((2*numerator-1)*(1-denominator), kind=rk) &
        &            / real(numerator*(2*denominator-1), kind=rk)
      end if
    else
      alpha_beta = 0.0_rk
    end if
    return

  end function ply_alpha_beta
  ! ************************************************************************ !

end module ply_poly_transformation_module 
