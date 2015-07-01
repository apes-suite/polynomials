!> This module provides the means to sample polynomial data to break it down
!! into voxels for visualization.
!!
!! Not to be confused with the oversample module!
module ply_sampling_module
  use aotus_module, only: flu_state, aot_get_val, aoterr_Fatal
  use aot_table_module, only: aot_table_open, aot_table_close

  use tem_aux_module, only: tem_abort
  use tem_logging_module, only: logunit

  implicit none

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
  end type ply_sampling_type


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
        &                 // ' be interpreted for a number.'
      write(logunit(1),*) 'Aborting the execution, please check your config!'
      call tem_abort()
    end if

    call aot_table_close( L       = conf,   &
      &                   thandle = thandle )

  end subroutine ply_sampling_load
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!

end module ply_sampling_module
