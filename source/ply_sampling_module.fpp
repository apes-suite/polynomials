?? include "ply_dof_module.inc"
!> This module provides the means to sample polynomial data to break it down
!! into voxels for visualization.
!!
!! Not to be confused with the oversample module!
module ply_sampling_module
  use iso_c_binding, only: c_f_pointer, c_loc
  use env_module, only: labelLen, rk, long_k

  use aotus_module, only: flu_state, aot_get_val, aoterr_Fatal
  use aot_table_module, only: aot_table_open, aot_table_close

  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort
  use tem_bc_prop_module, only: tem_bc_prop_type
  use tem_logging_module, only: logunit
  use tem_refining_module, only: tem_refine_global_subtree
  use tem_subtree_module, only: tem_create_subtree_of, tem_create_tree_from_sub
  use tem_subtree_type_module, only: tem_subtree_type, tem_destroy_subtree
  use tem_time_module, only: tem_time_type
  use tem_tools_module, only: upper_to_lower
  use tem_topology_module, only: tem_coordofid
  use tem_tracking_module, only: tem_tracking_type
  use tem_varsys_module, only: tem_varSys_proc_element, tem_varSys_proc_point, &
    &                          tem_varSys_proc_getParams, &
    &                          tem_varSys_proc_setParams, &
    &                          tem_varSys_proc_setupIndices, &
    &                          tem_varSys_proc_getValOfIndex, &
    &                          tem_varsys_append_statevar, &
    &                          tem_varSys_init, &
    &                          tem_varSys_type, tem_varSys_op_type, &
    &                          tem_varSys_getParams_dummy, &
    &                          tem_varSys_setParams_dummy

  use ply_dof_module, only: q_space, p_space,                           &
    &                       nextModgCoeffQTens, nextModgCoeffPTens,     &
    &                       nextModgCoeffQTens2D, nextModgCoeffPTens2D, &
    &                       nextModgCoeffQTens1D, nextModgCoeffPTens1D
  use ply_modg_basis_module, only: legendre_1D

  implicit none

  private

  !> This is  the data type providing the definitions for the sampling.
  type ply_sampling_type
    !> Maximal number of levels by which any mesh element should be refined.
    !!
    !! A setting of 0 results in no sampling, and the original mesh elements
    !! will be used with the integral mean value (first degree of freedom).
    !! Higher levels provide a limit for the refinement of the mesh.
    !! Note, that even for large settings here, the overall mesh depth is
    !! restricted by the global limit, due to the available space in the
    !! integers representing the treeIDs.
    integer :: max_nlevels = 0

    !> Method to use for the sampling.
    character(len=labelLen) :: method
  end type ply_sampling_type


  public :: ply_sampling_type
  public :: ply_sampling_load
  public :: ply_sample_data
  public :: ply_sampling_free_methodData


  !> Private data type to describe variables with varying polynomial
  !! representation from element to element for each variable.
  type vdata_type
    ! I think we need to make use of growing arrays here.
    ! One for the space (int), nElems
    ! One for the maxdgree (int), nElems
    ! One for the actual data (real kind=rk), nDofs across all elements
  end type vdata_type


  !> Required to make use of c_loc and accessing an array.
  type capsule_array_type
    real(kind=rk), allocatable :: dat(:)
  end type capsule_array_type

contains


  !----------------------------------------------------------------------------!
  !> This subroutine reads the sampling configuration from the Lua script
  !! provided in conf and fills the sampling data in 'me' accordingly.
  !!
  !! The data is expected to be stored in a table under the name 'ply_sampling'.
  subroutine ply_sampling_load(me, conf, parent)
    !> Sampling definition to load.
    type(ply_sampling_type), intent(out) :: me

    !> Configuration to read the sampling settings from.
    type(flu_State), intent(in) :: conf

    !> Parent table in which to look for the sampling settings.
    integer, intent(in), optional :: parent
    !----------------------------------------------------------------------!
    integer :: thandle, iError
    !----------------------------------------------------------------------!

    call aot_table_open( L       = conf,          &
      &                  parent  = parent,        &
      &                  thandle = thandle,       &
      &                  key     = 'ply_sampling' )

    call aot_get_val( L       = conf,           &
      &               thandle = thandle,        &
      &               key     = 'nlevels',      &
      &               val     = me%max_nlevels, &
      &               ErrCode = iError,         &
      &               default = 0               )

    if ( btest(iError, aoterr_Fatal) ) then
      write(logunit(1),*) 'ERROR: nlevels for ply_sampling needs to be an' &
        &                 // ' integer!'
      write(logunit(1),*) 'You provided nlevels but not in a form that could' &
        &                 // ' be interpreted as a number.'
      write(logunit(1),*) 'Aborting the execution, please check your config!'
      call tem_abort()
    end if

    if (me%max_nLevels > 0) then
      !> Actually need to do sampling, read in the rest of configuration
      !! settings.
      call aot_get_val( L       = conf,      &
        &               thandle = thandle,   &
        &               key     = 'method',  &
        &               val     = me%method, &
        &               ErrCode = iError,    &
        &               default = 'fixed'    )

      if ( btest(iError, aoterr_Fatal) ) then
        write(logunit(1),*) 'ERROR: method for ply_sampling needs to be an' &
          &                 // ' string!'
        write(logunit(1),*) 'You provided a method but not in a form that' &
          &                 // ' could be interpreted as a string.'
        write(logunit(1),*) 'Aborting the execution, please check your config!'
        call tem_abort()
      end if

      me%method = adjustl(me%method)
      me%method = upper_to_lower(me%method)

      select case(trim(me%method))
      case('fixed')
        write(logunit(1),'(a,i0,a)') 'Sampling configured to be fixed with ', &
          &                          me%max_nlevels, &
          &                          ' levels.'
      case default
        write(logunit(1),*) 'ERROR: unknown sampling method: ', trim(me%method)
        write(logunit(1),*) '       The sampling method needs to be one of the'
        write(logunit(1),*) '       following:'
        write(logunit(1),*) '       * "fixed" - will refine all elements by a'
        write(logunit(1),*) '                   given number of levels and'
        write(logunit(1),*) '                   provide the polynomial values'
        write(logunit(1),*) '                   at the barycenters of those.'
        write(logunit(1),*) ''
        write(logunit(1),*) 'Stopping!'
        call tem_abort()
      end select
    end if

    call aot_table_close( L       = conf,   &
      &                   thandle = thandle )

  end subroutine ply_sampling_load
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!


  !----------------------------------------------------------------------------!
  !> Sampling polynomial data from a given array and mesh to a new mesh with
  !! a new data array, where just a single degree of freedom per element is
  !! used.
  subroutine ply_sample_data( me, orig_mesh, orig_bcs, varsys, var_degree, &
    &                         var_space, ndims, tracking, time, new_mesh,  &
    &                         resvars                                      )
    !> A ply_sampling_type to describe the sampling method.
    type(ply_sampling_type), intent(in) :: me

    !> The original mesh to be refined.
    type(treelmesh_type), intent(in) :: orig_mesh

    !> Boundary conditions for the original mesh.
    type(tem_BC_prop_type), intent(in) :: orig_bcs

    type(tem_varsys_type), intent(in) :: varsys

    !> Maximal polynomial degree for each variable.
    !!
    !! Needs to be matching the variable definition in the variable system.
    !! @todo Needs to be changed to be an information per element per variable!
    !!       Possibly by defining a variable in the varsys, providing the
    !!       degree.
    integer, intent(in) :: var_degree(:)

    !> Polynomial space for each variable.
    !!
    !! Needs to be matching the variable definition in the variable system.
    integer, intent(in) :: var_space(:)

    !> Number of dimensions in the polynomial representation.
    integer, intent(in) :: ndims

    type(tem_tracking_type), intent(in) :: tracking

    type(tem_time_type), intent(in) :: time

    !> The new mesh with the refined elements.
    type(treelmesh_type), intent(out) :: new_mesh

    !> Resulting system of variables describing the data in the arrays of
    !! subsampled elements.
    type(tem_varsys_type), intent(out) :: resvars
    !----------------------------------------------------------------------!
    type(treelmesh_type) :: tmp_mesh(0:1)
    type(tem_BC_prop_type) :: tmp_bcs(0:1)
    type(tem_subtree_type) :: tmp_subtree
    type(capsule_array_type), pointer :: res
    real(kind=rk), allocatable :: vardat(:)
    real(kind=rk), allocatable :: points(:)
    real(kind=rk), allocatable :: pointval(:,:)
    integer, allocatable :: vardofs(:)
    integer, allocatable :: elempos(:)
    integer :: pointCoord(4)
    integer :: nOrigElems
    integer :: nVars
    integer :: nComponents
    integer :: nChilds
    integer :: cur, prev
    integer :: ans(3)
    integer :: i
    integer :: iLevel
    integer :: iChild
    integer :: iComp
    integer :: ivar
    integer :: iDof
    integer :: iElem
    integer :: iProp
    integer :: varpos
    integer :: childpos, parentpos
    integer :: maxdofs
    integer :: n1D_childs
    integer :: lastdegree
    integer :: ii
    integer :: fak(2)
    integer :: bitlevel
    real(kind=rk) :: legval
    real(kind=rk) :: point_spacing, point_start
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => NULL()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => NULL()
    procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex &
      &                                                  => NULL()
    !----------------------------------------------------------------------!

    ! Procedure to do...
    ! (Internally) create:
    ! Arrays of data, and description of its layout FOR EACH VARIABLE:
    !   Polynomial space and degree for each element
    !   -> Define a datatype for this - see vdata_type

    ! Refine the given mesh according to the configuration in the
    ! ply_sampling_type.
    ! And fill the final data array accordingly.

    select case(trim(me%method))
    case('fixed')
      cur = 0
      call tem_refine_global_subtree( orig_mesh = orig_mesh,        &
        &                             orig_bcs  = orig_bcs,         &
        &                             subtree   = tracking%subtree, &
        &                             ndims     = ndims,            &
        &                             new_mesh  = tmp_mesh(cur),    &
        &                             new_bcs   =  tmp_bcs(cur),    &
        &                             restrict_to_sub = .true.      )
      call tem_create_subTree_of( inTree  = tmp_mesh(0),             &
        &                         bc_prop = tmp_bcs(0),              &
        &                         subtree = tmp_subtree,             &
        &                         inShape = tracking%header%geometry )
      do iLevel=2,me%max_nLevels
        write(logunit(6),*) 'sampling level ', iLevel
        prev = mod(iLevel-2, 2)
        cur = mod(iLevel-1, 2)
        call tem_refine_global_subtree( orig_mesh = tmp_mesh(prev), &
          &                             orig_bcs  = tmp_bcs(prev),  &
          &                             subtree   = tmp_subtree,    &
          &                             ndims     = ndims,          &
          &                             new_mesh  = tmp_mesh(cur),  &
          &                             new_bcs   = tmp_bcs(cur),   &
          &                             restrict_to_sub = .true.    )
        call tem_destroy_subtree(tmp_subtree)
        nullify(tmp_bcs(prev)%property)
        nullify(tmp_bcs(prev)%header)
        if (associated(tmp_mesh(prev)%property)) then
          do iProp=1,size(tmp_mesh(prev)%property)
            deallocate(tmp_mesh(prev)%property(iProp)%ElemID)
          end do
          deallocate(tmp_mesh(prev)%property)
          deallocate(tmp_mesh(prev)%global%property)
        end if
        call tem_create_subTree_of( inTree  = tmp_mesh(cur),           &
          &                         bc_prop = tmp_bcs(cur),            &
          &                         subtree = tmp_subtree,             &
          &                         inShape = tracking%header%geometry )
      end do

      call tem_create_tree_from_sub( intree  = tmp_mesh(cur), &
        &                            subtree = tmp_subtree,   &
        &                            newtree = new_mesh       )

      maxdofs = 0
      nVars = tracking%varmap%varPos%nVals
      call tem_varSys_init( me         = resvars,       &
        &                   systemName = 'sampledVars', &
        &                   length     = nVars          )

      allocate(vardofs(nVars))
      do ivar=1,tracking%varmap%varPos%nVals
        varpos = tracking%varmap%varPos%val(iVar)
        select case(var_space(varpos))
        case (q_space)
          select case(ndims)
          case(3)
?? copy :: getDofsQTens(var_degree(varpos), vardofs(iVar))
          case(2)
?? copy :: getDofsQTens2D(var_degree(varpos), vardofs(iVar))
          case(1)
?? copy :: getDofsQTens1D(var_degree(varpos), vardofs(iVar))
          end select
          maxdofs = max( maxdofs, varsys%method%val(varpos)%nComponents &
            &                     * vardofs(iVar)                       )
        case (p_space)
          select case(ndims)
          case(3)
?? copy :: getDofsPTens(var_degree(varpos), vardofs(iVar))
          case(2)
?? copy :: getDofsPTens2D(var_degree(varpos), vardofs(iVar))
          case(1)
?? copy :: getDofsPTens1D(var_degree(varpos), vardofs(iVar))
          end select
          maxdofs = max( maxdofs, varsys%method%val(varpos)%nComponents &
            &                     * vardofs(iVar)                       )
        end select
      end do

      if (tracking%subtree%useGlobalMesh) then
        ! All elements are refined...
        nOrigElems = orig_mesh%nElems
        allocate( elempos(nOrigElems) )
        elempos = [ (i, i=1,nOrigElems) ]

      else
        nOrigElems = tracking%subtree%nElems
        allocate( elempos(nOrigElems) )
        elempos = tracking%subtree%map2global

      end if

      n1D_childs = 2**me%max_nLevels
      point_spacing = 2.0_rk / real(n1D_childs, kind=rk)
      point_start = 0.5_rk * point_spacing - 1.0_rk

      allocate(vardat(maxdofs*nOrigElems))
      allocate(points(n1D_childs))

      get_element => get_sampled_element
      get_params => tem_varSys_getparams_dummy
      set_params => tem_varSys_setparams_dummy
      nullify(get_point)
      nullify(setup_indices)
      nullify(get_valOfIndex)

      do iChild=1,n1D_childs
        points(iChild) = point_start + ((iChild-1) * point_spacing)
      end do

      varpos = tracking%varmap%varPos%val(1)
      allocate(pointval(var_degree(varpos)+1, n1D_childs))
      pointval = legendre_1D(points = points, degree = var_degree(varpos))
      lastdegree = var_degree(varpos)

      do ivar=1,tracking%varmap%varPos%nVals

        varpos = tracking%varmap%varPos%val(iVar)
        nComponents = varsys%method%val(varpos)%nComponents

        if (var_degree(varpos) /= lastdegree) then
          deallocate(pointval)
          allocate(pointval(var_degree(varpos)+1, n1D_childs))
          pointval = legendre_1D(points = points, degree = var_degree(varpos))
          lastdegree = var_degree(varpos)
        end if

        call varSys%method%val(varpos)                    &
          &        %get_element( varSys  = varSys,        &
          &                      elempos = elempos,       &
          &                      time    = time,          &
          &                      tree    = orig_mesh,     &
          &                      nElems  = nOrigElems,    &
          &                      nDofs   = vardofs(iVar), &
          &                      res     = vardat         )

        if (tracking%subtree%useGlobalMesh) then

          nChilds = (2**ndims)**me%max_nLevels
          allocate(res)
          allocate(res%dat(nComponents*nChilds*nOrigElems))

          do iChild=1,nChilds
            select case(ndims)
            case (3)
              pointCoord = tem_CoordofID( int(iChild, kind=long_k), &
                &                         offset = 1_long_k       ) &
                &        + 1
            case (2)
              ! 2D Bit-sieving coordofid
              bitlevel = 1
              pointCoord = 1
              fak(1) = 1
              fak(2) = 2
              do
                if ((iChild-1) / fak(1) == 0) EXIT
                do ii=1,2
                  pointCoord(ii) = pointCoord(ii) + bitlevel &
                    &                             * mod((iChild-1) / fak(ii), 2)
                end do
                bitlevel = bitlevel*2
                fak = fak*4
              end do
            case (1)
              pointCoord = 1
              pointCoord(1) = iChild
            end select

            iDof = 1
            ans(1) = 1
            ans(2) = 1
            ans(3) = 1
            legval = pointval(ans(1), pointCoord(1))
            do ii=2,ndims
              legval = legval * pointval(ans(ii), pointCoord(ii))
            end do

            do iElem=1,nOrigElems
              parentpos = (iElem-1)*vardofs(iVar)*nComponents
              childpos = (iElem - 1)*nChilds*nComponents &
                &      + (iChild - 1)*nComponents
              do iComp=1,nComponents
                res%dat(childpos+iComp) = legval*vardat(parentpos+iComp)
              end do
            end do

            do iDof=2,vardofs(iVar)
              if (var_space(iVar) == q_space) then
                select case(ndims)
                case (3)
                  call nextModgCoeffQTens( ansFuncX  = ans(1),            &
                    &                      ansFuncY  = ans(2),            &
                    &                      ansFuncZ  = ans(3),            &
                    &                      maxdegree = var_degree(varpos) )
                case (2)
                  call nextModgCoeffQTens2D( ansFuncX  = ans(1),            &
                    &                        ansFuncY  = ans(2),            &
                    &                        ansFuncZ  = ans(3),            &
                    &                        maxdegree = var_degree(varpos) )
                case (1)
                  call nextModgCoeffQTens1D( ansFuncX  = ans(1),            &
                    &                        ansFuncY  = ans(2),            &
                    &                        ansFuncZ  = ans(3),            &
                    &                        maxdegree = var_degree(varpos) )
                end select
              else
                select case(ndims)
                case (3)
                  call nextModgCoeffPTens( ansFuncX  = ans(1),            &
                    &                      ansFuncY  = ans(2),            &
                    &                      ansFuncZ  = ans(3),            &
                    &                      maxdegree = var_degree(varpos) )
                case (2)
                  call nextModgCoeffPTens2D( ansFuncX  = ans(1),            &
                    &                        ansFuncY  = ans(2),            &
                    &                        ansFuncZ  = ans(3),            &
                    &                        maxdegree = var_degree(varpos) )
                case (1)
                  call nextModgCoeffPTens1D( ansFuncX  = ans(1),            &
                    &                        ansFuncY  = ans(2),            &
                    &                        ansFuncZ  = ans(3),            &
                    &                        maxdegree = var_degree(varpos) )
                end select
              end if
              legval = pointval(ans(1), pointCoord(1))
              do ii=2,ndims
                legval = legval * pointval(ans(ii), pointCoord(ii))
              end do
              do iElem=1,nOrigElems
                parentpos = (iElem-1)*vardofs(iVar)*nComponents &
                  &       + (iDof-1)*nComponents
                childpos = (iElem - 1)*nChilds*nComponents &
                  &      + (iChild - 1)*nComponents
                do iComp=1,nComponents
                  res%dat(childpos+iComp) = res%dat(childpos+iComp) &
                    &                     + legval*vardat(parentpos+iComp)
                end do
              end do
            end do

          end do

          ! assign res as method data to the variable in the resvars.
          call tem_varSys_append_stateVar( me          = resvars,            &
            &                              varname     = varsys%varname      &
            &                                                  %val(varpos), &
            &                              nComponents = nComponents,        &
            &                              method_data = c_loc(res),         &
            &                              set_params     = set_params,     &
            &                              get_point      = get_point,      &
            &                              get_element    = get_element,    &
            &                              get_params     = get_params,     &
            &                              setup_indices  = setup_indices,  &
            &                              get_valofindex = get_valofindex  )

          ! now nullify res again, to allow its usage for another allocation:
          nullify(res)

        else
          write(logunit(1),*) 'Not yet supported non-global subtrees!'
          write(logunit(1),*) 'Stopping...'
          call tem_abort()
        end if

      end do

      deallocate(pointval)

    case default
      write(logunit(1),*) 'Not implemented sampling method!'
      write(logunit(1),*) 'Stopping...'
      call tem_abort()
    
    end select

  end subroutine ply_sample_data
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
 

  !----------------------------------------------------------------------------!
  !> Get sampled data.
  !!
  !! This routine provides the get_element function of the variable definition
  !! to access the sampled data array obtained by ply_sample_data.
  subroutine get_sampled_element( fun, varsys, elempos, time, tree, n, &
    &                             nDofs, res                           )
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of elements to obtain for this variable (vectorized access).
    integer, intent(in) :: n

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (nComponents of resulting variable) x (nDegrees of freedom) x (nElems)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    !----------------------------------------------------------------------!
    type(capsule_array_type), pointer :: p
    integer :: datlen(1)
    integer :: iElem
    integer :: nComps
    !----------------------------------------------------------------------!

    nComps = fun%nComponents
    datlen = tree%nElems * nComps

    call c_f_pointer(fun%method_data, p)

    do iElem=1,n
      res(1+(iElem-1)*nComps:iElem*nComps) &
        &  = p%dat(1+(elempos(iElem)-1)*nComps:elempos(iElem)*nComps)
    end do

  end subroutine get_sampled_element
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!


  !----------------------------------------------------------------------------!
  !> Free previously allocated methodData of variable.
  !!
  !! This routine provides a method to free allocated methodData again.
  subroutine ply_sampling_free_methodData(fun)
    !> Description of the method to free the data for.
    class(tem_varSys_op_type), intent(inout) :: fun

    !----------------------------------------------------------------------!
    type(capsule_array_type), pointer :: p
    !----------------------------------------------------------------------!

    call c_f_pointer(fun%method_data, p)
    if (associated(p)) then
      if (allocated(p%dat)) then
        deallocate(p%dat)
      end if
      nullify(p)
      fun%method_data = c_loc(p)
    end if

  end subroutine ply_sampling_free_methodData

end module ply_sampling_module
