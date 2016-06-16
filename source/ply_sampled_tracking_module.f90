!> Module that implements tracking with subsampling of polynomials.
module ply_sampled_tracking_module
  use aotus_module, only: flu_State

  use hvs_output_module, only: hvs_output_init, hvs_output_open, &
    &                          hvs_output_write, hvs_output_close

  use env_module, only: pathLen, labelLen

  use treelmesh_module, only: treelmesh_type, unload_treelmesh
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
    &                            tem_tracker,                &
    &                            tem_init_tracker,           &
    &                            tem_trackingControl_type,   &
    &                            tem_trackingHeader_type,    &
    &                            tem_tracking_has_triggered, &
    &                            tem_init_tracker_subtree,   &
    &                            tem_load_tracking
  use tem_varMap_module, only: tem_create_varMap
  use tem_varSys_module, only: tem_varsys_type, tem_empty_varSys

  use ply_sampling_module, only: ply_sampling_type,            &
    &                            ply_sampling_load,            &
    &                            ply_sampling_free_methodData, &
    &                            ply_sample_data

  implicit none

  private

  type ply_sampled_tracking_type
    !> Control for the tracking in total.
    type(tem_trackingControl_type) :: trackCtrl

    !> Individual tracking entities.
    type(tem_tracking_type), allocatable :: tracking(:)

    !> Headers for the configuration of each tracking entity.
    type(tem_trackingHeader_type), allocatable :: trackingHeader(:)

    !> Subsampled mesh for each tracking.
    !!
    !!@todo Actually make use of these, instead of regenerating the mesh
    !!      every time the tracking is written.
    type(treelmesh_type), allocatable :: mesh(:)

    !> Variable system description after subsampling.
    !!
    !!@todo Actuall make use of these, instead of recreating the variable
    !!      system each time a tracking is written.
    type(tem_varSys_type), allocatable :: varsys(:)

    !> Configuration of the subsampling (applied to all trackings).
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
    &                                prefix, stencil, proc, ndofs  )
    !> Sampled tracking variable to initialize. It has to be configured by
    !! [[ply_sampled_tracking_load]] beforehand.
    type(ply_sampled_tracking_type), intent(inout) :: me

    !> The global mesh.
    type(treelmesh_type), intent(in) :: mesh

    !> Information about the solver (used to construct file name strings).
    type(tem_solveHead_type), intent(in) :: solver

    !> Global variable system with description of the data to get the
    !! tracking variables from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Boundary condition properties, used to identify elements close to
    !! the boundary.
    type(tem_bc_prop_type), intent(in) :: bc

    !> A prefix to use in the construction of filenames, can be used to
    !! set a directory to write all tracking data to.
    !!
    !! Defaults to a prefix given by solver%simName and varsys%SystemName
    !! (not a directory).
    character(len=labelLen), optional, intent(in) :: prefix

    !> Description of the stencil in the numerical scheme.
    !!
    !! This is needed to describe elements adjacent to specific boundary
    !! labels.
    type(tem_stencilHeader_type), optional, intent(in) :: stencil

    !> General communication environment
    type(tem_comm_env_type), intent(in) :: proc

    !> Number of degrees of freedom to use in the output.
    integer, intent(in) :: nDofs
    ! -------------------------------------------------------------------- !
    integer :: iTrack
    integer :: iVar
    integer :: nVars
    ! -------------------------------------------------------------------- !

    ! Initialize tracker subTree and remove empty trackers.
    call tem_init_tracker_subTree( me      = me%tracking,  &
      &                            tCtrl   = me%trackCtrl, &
      &                            tree    = mesh,         &
      &                            solver  = solver,       &
      &                            bc_prop = bc,           &
      &                            stencil = stencil,      &
      &                            prefix  = prefix        )

    if (me%sampling%max_nlevels == 0) then
      ! No subsampling to be done, call the general initialization and
      ! exit the routine.
      call tem_init_tracker( me       = me%tracking,  &
        &                    tCtrl    = me%trackCtrl, &
        &                    tree     = mesh,         &
        &                    solver   = solver,       &
        &                    varSys   = varSys,       &
        &                    nDofs    = nDofs,        &
        &                    globProc =  proc         )
      RETURN
    end if

    do iTrack=1,me%trackCtrl%nTrackings

      ! map variables
      ! create tracking variable position in the global varSys
      call tem_create_varMap( varname = me%tracking(iTrack)%header   &
        &                                                  %varname, &
        &                     varSys  = varSys,                      &
        &                     varMap  = me%tracking(iTrack)%varMap   )

      nVars = me%tracking(iTrack)%varMap%varPos%nVals
      ! Abort if none of the variables defined in current
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
  !> Output sampled tracking data.
  !!
  !! Iterates over all tracking instances in the given me variable, checks
  !! whether it should be written at the current point in time (if simControl
  !! is provided), subsamples the data and performs the hvs_output for the
  !! subsampled data.
  !!
  !!@todo Instead of recreating the sampled varsys and mesh everytime the
  !!      tracking is written, store them in the [[ply_sampled_tracking_type]].
  subroutine ply_sampled_track_output( me, mesh, bc, solver, proc, varSys, &
    &                                  var_degree, var_space, simControl )
    !> Sampled tracking instances.
    type(ply_sampled_tracking_type), intent(inout)  :: me

    !> Global mesh, required for the sampling.
    type(treelmesh_type), intent(in)                :: mesh

    !> Boundary properties, needed to inherit boundary information to refined
    !! meshes and allow the extraction of boundary shape geometries.
    type(tem_bc_prop_type), intent(in)              :: bc

    !> Information about the solver, needed for the output file name.
    type(tem_solveHead_type), intent(in)            :: solver

    !> General communication environment
    type(tem_comm_env_type), intent(in)             :: proc

    !> Original variable system
    type(tem_varSys_type), intent(in)               :: varSys

    !> Maximal polynomial degree for each variable
    !!
    !! Needs to match the size of the variable system.
    integer, intent(in)                             :: var_degree(:)

    !> Maximal polynomial space for each variable
    !!
    !! Needs to match the size of the variable system.
    integer, intent(in)                             :: var_space(:)

    !> Simulation control to determine, whether trackings should be written
    !!
    !! If not provided, all trackings will be written unconditionally.
    type(tem_simControl_type), intent(in), optional :: simControl
    ! -------------------------------------------------------------------- !
    character(len=pathLen) :: basename
    type(tem_time_type)    :: time
    type(tem_varsys_type)  :: sampled_vars
    type(treelmesh_type)   :: sampled_mesh
    integer :: iTrack
    integer :: iVar
    ! -------------------------------------------------------------------- !

    call tem_time_reset(time)
    if (present(simControl)) then
      time = simControl%now
      if (me%sampling%max_nlevels == 0) then
        ! No subsampling to be done, call the regular tracker, and leave
        ! the routine.
        call tem_tracker( track      = me%tracking, &
          &               simControl = simControl,  &
          &               varSys     = varsys,      &
          &               tree       = mesh         )
        RETURN
      end if
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

      do ivar=1,sampled_vars%method%nVals
        call ply_sampling_free_methodData(sampled_vars%method%val(iVar))
      end do
      call tem_empty_varSys(sampled_vars)
      call unload_treelmesh(sampled_mesh)

    end do

  end subroutine ply_sampled_track_output


end module ply_sampled_tracking_module
