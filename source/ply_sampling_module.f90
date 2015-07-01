!> This module provides the means to sample polynomial data to break it down
!! into voxels for visualization.
!!
!! Not to be confused with the oversample module!
module ply_sampling_module
  use env_module, only: labelLen

  use aotus_module, only: flu_state, aot_get_val, aoterr_Fatal
  use aot_table_module, only: aot_table_open, aot_table_close

  use tem_aux_module, only: tem_abort
  use tem_logging_module, only: logunit
  use tem_tools_module, only: upper_to_lower

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


  !> Private data type to describe variables with varying polynomial
  !! representation from element to element for each variable.
  type vdata_type
    ! I think we need to make use of growing arrays here.
    ! One for the space (int), nElems
    ! One for the maxdgree (int), nElems
    ! One for the actual data (real kind=rk), nDofs across all elements
  end type vdata_type


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
  subroutine ply_sample_data()
    ! A ply_sampling_type to describe the sampling method.
    ! Required input: tree/subtree, to describe the original mesh.
    !                 + tracking shape to restrict generated mesh to the tracked
    !                   shape (-> Just use the tracking object?)
    ! Variable system + initial dofs and polyspace for each variable.
    ! Output: New tree
    !         Array of data (1 dof x elements in new mesh for each variable)
    !         -> newelements x nVariables
    ! Possibly new variable system to access the data in that array?

    ! Procedure to do...
    ! (Internally) create:
    ! Arrays of data, and description of its layout FOR EACH VARIABLE:
    !   Polynomial space and degree for each element
    !   -> Define a datatype for this - see vdata_type

    ! Refine the given mesh according to the configuration in the
    ! ply_sampling_type.
    ! And fill the final data array accordingly.
  end subroutine ply_sample_data
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!


end module ply_sampling_module
