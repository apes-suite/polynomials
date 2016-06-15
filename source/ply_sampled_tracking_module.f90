module ply_sampled_tracking_module
  use aotus_module, only: flu_State

  use hvs_output_module, only: hvs_output_init, hvs_output_open, &
    &                          hvs_output_write, hvs_output_close

  use env_module, only: pathLen, labelLen

  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_abort
  use tem_bc_prop_module, only: tem_bc_prop_type
  use tem_comm_env_module, only: tem_comm_env_type
  use tem_logging_module, only: logunit
  use tem_simControl_module, only: tem_simControl_type
  use tem_solveHead_module, only: tem_solveHead_type
  use tem_stencil_module, only: tem_stencilHeader_type
  use tem_time_module, only: tem_time_type, &
    &                        tem_time_reset
  use tem_tracking_module, only: tem_tracking_type,          &
    &                            tem_trackingControl_type,   &
    &                            tem_trackingHeader_type,    &
    &                            tem_tracking_has_triggered, &
    &                            tem_init_tracker_subtree,   &
    &                            tem_load_tracking
  use tem_varMap_module, only: tem_create_varMap
  use tem_varSys_module, only: tem_varsys_type

  use ply_sampling_module, only: ply_sampling_type, &
    &                            ply_sampling_load, &
    &                            ply_sample_data

  implicit none

  private

  type ply_sampled_tracking_type
    type(tem_trackingControl_type) :: trackCtrl
    type(tem_tracking_type), allocatable :: tracking(:)
    type(tem_trackingHeader_type), allocatable :: trackingHeader(:)
    type(treelmesh_type), allocatable :: mesh(:)
    type(tem_varSys_type), allocatable :: varsys(:)
    type(ply_sampling_type) :: sampling
  end type ply_sampled_tracking_type

  public :: ply_sampled_tracking_type
  public :: ply_sampled_tracking_load


contains


  ! ------------------------------------------------------------------------ !
  !> Load the configuration of sampled tracking objects.
  subroutine ply_sampled_tracking_load(me, conf)
    !> Sampled tracking data to load from the config
    type(ply_sampled_tracking_type), intent(out) :: me

    !> Lua config to load the tracking from
    type(flu_State) :: conf
    ! -------------------------------------------------------------------- !
    integer :: iTrack
    ! -------------------------------------------------------------------- !

    call ply_sampling_load( me   = me%sampling, &
      &                     conf = conf         )

    call tem_load_tracking( tCtrl  = me%trackCtrl,      &
      &                     header = me%trackingHeader, &
      &                     conf   = conf               )

    allocate(me%tracking(me%trackCtrl%nTrackings))

    do iTrack=1,me%trackCtrl%nTrackings
      me%tracking(iTrack)%header = me%trackingHeader(iTrack)
    end do

  end subroutine ply_sampled_tracking_load
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! Initialize the sampled tracking entities.
  !
  ! This is necessary to properly setup the tem_tracking data.
  ! It includes building the subtree and the varmap.
  subroutine ply_sampled_track_init( me, mesh, solver, varSys, bc, &
    &                                prefix, stencil               )
    type(ply_sampled_tracking_type), intent(inout) :: me
    type(treelmesh_type), intent(in) :: mesh
    type(tem_solveHead_type), intent(in) :: solver
    type(tem_varSys_type), intent(in) :: varSys
    type(tem_bc_prop_type), intent(in) :: bc
    character(len=labelLen), optional, intent(in) :: prefix
    type(tem_stencilHeader_type), optional, intent(in) :: stencil
    ! -------------------------------------------------------------------- !
    integer :: iTrack
    integer :: iVar
    integer :: nVars
    ! -------------------------------------------------------------------- !

    ! Initialize tracker subTree and remote empty trackers.
    call tem_init_tracker_subTree( me      = me%tracking,  &
      &                            tCtrl   = me%trackCtrl, &
      &                            tree    = mesh,         &
      &                            solver  = solver,       &
      &                            varSys  = varSys,       &
      &                            bc_prop = bc,           &
      &                            stencil = stencil,      &
      &                            prefix  = prefix        )

    do iTrack=1,me%trackCtrl%nTrackings

      ! map variables
      ! create tracking variable position in the global varSys
      call tem_create_varMap( varname = me%tracking(iTrack)%header   &
        &                                                  %varname, &
        &                     varSys  = varSys,                      &
        &                     varMap  = me%tracking(iTrack)%varMap   )

      nVars = me%tracking(iTrack)%varMap%varPos%nVals
      ! Abort if none variables of the variable defined in current
      ! tracking object are found in varSys
      if (nVars==0) then
        write(logUnit(1),*) 'Error: Requested variables: '
        do iVar = 1, size(me%tracking(iTrack)%header%varName)
          write(logUnit(1),*) iVar, &
            &                 trim(me%tracking(iTrack)%header%varName(iVar))
        end do
        write(logUnit(1),*) 'not found in varSys.'
        write(logUnit(1),*) 'Check tracking object: '// &
          &                  trim(me%tracking(iTrack)%header%label)
        call tem_abort()
      end if

    end do

    if (me%trackCtrl%nTrackings > 0) then
      allocate(me%mesh(me%trackCtrl%nTrackings))
      allocate(me%varsys(me%trackCtrl%nTrackings))
    end if

  end subroutine ply_sampled_track_init
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  subroutine ply_sampled_track_output( me, mesh, bc, solver, proc, varSys, &
    &                                  var_degree, var_space, simControl )
    type(ply_sampled_tracking_type), intent(inout)  :: me
    type(treelmesh_type), intent(in)                :: mesh
    type(tem_bc_prop_type), intent(in)              :: bc
    type(tem_solveHead_type), intent(in)            :: solver
    type(tem_comm_env_type), intent(in)             :: proc
    type(tem_varSys_type), intent(in)               :: varSys
    integer, intent(in)                             :: var_degree(:)
    integer, intent(in)                             :: var_space(:)
    type(tem_simControl_type), intent(in), optional :: simControl
    ! -------------------------------------------------------------------- !
    character(len=pathLen) :: basename
    type(tem_time_type)    :: time
    type(tem_varsys_type)  :: sampled_vars
    type(treelmesh_type)   :: sampled_mesh
    integer :: iTrack
    ! -------------------------------------------------------------------- !

    call tem_time_reset(time)
    if (present(simControl)) then
      time = simControl%now
    end if

    do iTrack=1,me%trackCtrl%nTrackings
      if (present(simControl)) then
        ! If a simControl is provided, check each tracking on whether it is to
        ! be written. Without a simControl, we will write all trackings
        ! unconditionally.
        if ( .not. tem_tracking_has_triggered(           &
          &            track      = me%tracking(iTrack), &
          &            simControl = simControl         ) ) CYCLE
      end if

      call ply_sample_data( me         = me%sampling,         &
        &                   orig_mesh  = mesh,                &
        &                   orig_bcs   = bc,                  &
        &                   varsys     = varsys,              &
        &                   var_degree = var_degree,          &
        &                   var_space  = var_space,           &
        &                   tracking   = me%tracking(iTrack), &
        &                   time       = time,                &
        &                   new_mesh   = sampled_mesh,        &
        &                   resvars    = sampled_vars         )

      ! initialize output
      basename = trim(me%tracking(iTrack)%header%prefix) &
        &        // trim(me%tracking(iTrack)%header%label)
      call hvs_output_init(                                                   &
        &             out_file    = me%tracking(iTrack)%output_file,          &
        &             out_config  = me%tracking(iTrack)%header%output_config, &
        &             tree        = sampled_mesh,                             &
        &             varSys      = sampled_vars,                             &
        &             geometry    = me%tracking(iTrack)%header%geometry,      &
        &             basename    = trim(basename),                           &
        &             globProc    = proc,                                     &
        &             solver      = solver                                    )

      call hvs_output_open( out_file   = me%tracking(iTrack)%output_file,    &
        &                   use_iter   = me%tracking(iTrack)%header          &
        &                                  %output_config%vtk%iter_filename, &
        &                   mesh       = sampled_mesh,                       &
        &                   varsys     = sampled_vars                        )

      ! Fill output files with data.
      call hvs_output_write( out_file = me%tracking(iTrack)%output_file, &
        &                    varsys   = sampled_vars,                    &
        &                    mesh     = sampled_mesh                     )

      call hvs_output_close( out_file = me%tracking(iTrack)%output_file, &
        &                    varSys   = sampled_vars,                    &
        &                    mesh     = sampled_mesh                     )

    end do

  end subroutine ply_sampled_track_output


end module ply_sampled_tracking_module
