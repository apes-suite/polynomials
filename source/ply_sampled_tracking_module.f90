!> Module that implements tracking with subsampling of polynomials.
module ply_sampled_tracking_module
  use aotus_module, only: flu_State

  use hvs_output_module, only: hvs_output_init, hvs_output_open,   &
    &                          hvs_output_write, hvs_output_close, &
    &                          hvs_output_finalize

  use env_module, only: pathLen, labelLen

  use treelmesh_module, only: treelmesh_type, free_treelmesh
  use tem_aux_module, only: tem_abort
  use tem_bc_prop_module, only: tem_bc_prop_type
  use tem_comm_env_module, only: tem_comm_env_type
  use tem_logging_module, only: logunit
  use tem_reduction_module, only: tem_reduction_init
  use tem_simControl_module, only: tem_simControl_type
  use tem_solveHead_module, only: tem_solveHead_type
  use tem_stencil_module, only: tem_stencilHeader_type
  use tem_time_module, only: tem_time_type, &
    &                        tem_time_reset
  use tem_tracking_module, only: tem_tracking_type,          &
    &                            tem_tracker,                &
    &                            tem_init_tracker,           &
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
    !> Contains all tracking headers, control and tracking entities
    !! active on local process
    type(tem_tracking_type) :: tracking

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

    !> Dimensionality of the data to sample.
    integer :: ndims
  end type ply_sampled_tracking_type

  public :: ply_sampled_tracking_type
  public :: ply_sampled_tracking_load
  public :: ply_sampled_track_init
  public :: ply_sampled_track_output


contains


  ! ------------------------------------------------------------------------ !
  !> Load the configuration of sampled tracking objects.
  subroutine ply_sampled_tracking_load(me, conf)
    !> Sampled tracking data to load from the config
    type(ply_sampled_tracking_type), intent(out) :: me

    !> Lua config to load the tracking from
    type(flu_State) :: conf
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !

    call ply_sampling_load( me   = me%sampling, &
      &                     conf = conf         )

    call tem_load_tracking( me     = me%tracking, &
      &                     conf   = conf         )

  end subroutine ply_sampled_tracking_load
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! Initialize the sampled tracking entities.
  !
  ! This is necessary to properly setup the tem_tracking data.
  ! It includes building the subtree and the varmap.
  subroutine ply_sampled_track_init( me, mesh, solver, varSys, bc, &
    &                                stencil, proc, ndofs, ndims   )
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

    !> Description of the stencil in the numerical scheme.
    !!
    !! This is needed to describe elements adjacent to specific boundary
    !! labels.
    type(tem_stencilHeader_type), optional, intent(in) :: stencil

    !> General communication environment
    type(tem_comm_env_type), intent(in) :: proc

    !> Number of degrees of freedom to use in the output.
    integer, intent(in) :: nDofs

    !> Number of dimensions in the polynomial representations.
    integer, intent(in) :: nDims
    ! -------------------------------------------------------------------- !
    integer :: iTrack, iConfig
    integer :: iVar
    integer :: nVars
    character(len=pathLen) :: basename
    ! -------------------------------------------------------------------- !

    ! Initialize tracker subTree and remove empty trackers.
    call tem_init_tracker_subTree( me      = me%tracking,  &
      &                            tree    = mesh,         &
      &                            solver  = solver,       &
      &                            bc_prop = bc,           &
      &                            stencil = stencil       )

    me%nDims = nDims

    if (me%sampling%max_nlevels == 0) then
      ! No subsampling to be done, call the general initialization and
      ! exit the routine.
      call tem_init_tracker( me       = me%tracking,  &
        &                    tree     = mesh,         &
        &                    solver   = solver,       &
        &                    varSys   = varSys,       &
        &                    nDofs    = nDofs,        &
        &                    globProc = proc          )
      RETURN
    end if


    do iTrack = 1, me%tracking%control%nActive
      iConfig = me%tracking%instance(iTrack)%pntConfig 

      ! map variables
      ! create tracking variable position in the global varSys
      call tem_create_varMap( varname = me%tracking%config(iConfig)  &
        &                                                  %varname, &
        &                     varSys  = varSys,                      &
        &                     varMap  = me%tracking%instance(iTrack) &
        &                                          %varMap           )

      nVars = me%tracking%instance(iTrack)%varMap%varPos%nVals
      ! Abort if none of the variables defined in current
      ! tracking object are found in varSys
      if (nVars==0) then
        write(logUnit(1),*) 'Error: Requested variables: '
        do iVar = 1, size(me%tracking%config(iConfig)%varName)
          write(logUnit(1),*) iVar, &
            &                 trim(me%tracking%config(iConfig)%varName(iVar))
        end do
        write(logUnit(1),*) 'not found in varSys.'
        write(logUnit(1),*) 'Check tracking object: '// &
          &                  trim(me%tracking%config(iConfig)%label)
        call tem_abort()
      end if

      ! Init spatial reduction
      me%tracking%instance(iTrack)%output_file%ascii%isReduce &
        &  = me%tracking%config(iConfig)%reduction_config%active
      if ( me%tracking%config(iConfig)%reduction_config%active ) then
        call tem_reduction_init( &
          & me               = me%tracking%instance(iTrack)          &
          &                               %output_file%ascii%reduce, &
          & reduction_config = me%tracking%config(iConfig)           &
          &                               %reduction_config,         &
          & varSys           = varSys,                               &
          & varPos           = me%tracking%instance(iTrack)          &
          &                               %varMap%varPos             &
          &                               %val(:nVars)               )
      end if

      if (me%tracking%config(iConfig)%output_config%useGetPoint) then
        ! For point trackings do the initialization here, as no subsampling is
        ! required for them.
        basename = trim(me%tracking%config(iConfig)%prefix) &
          &        // trim(me%tracking%config(iConfig)%label)

        if (me%tracking%instance(iTrack)%subTree%useGlobalMesh) then
          call hvs_output_init(                             &
            &    out_file    = me%tracking%instance(iTrack) &
            &                    %output_file,              &
            &    out_config  = me%tracking%config(iConfig)  &
            &                    %output_config,            &
            &    tree        = mesh,                        &
            &    varSys      = varsys,                      &
            &    varPos      = me%tracking%instance(iTrack) &
            &                               %varMap%varPos  &
            &                               %val(:nVars),   &
            &    geometry    = me%tracking%config(iConfig)  &
            &                    %geometry,                 &
            &    basename    = trim(basename),              &
            &    globProc    = proc,                        &
            &    solver      = solver                       )
        else
          call hvs_output_init(                             &
            &    out_file    = me%tracking%instance(iTrack) &
            &                    %output_file,              &
            &    out_config  = me%tracking%config(iConfig)  &
            &                    %output_config,            &
            &    tree        = mesh,                        &
            &    subtree     = me%tracking%instance(iTrack) &
            &                    %subtree,                  &
            &    varSys      = varsys,                      &
            &    varPos      = me%tracking%instance(iTrack) &
            &                               %varMap%varPos  &
            &                               %val(:nVars),   &
            &    geometry    = me%tracking%config(iConfig)  &
            &                    %geometry,                 &
            &    basename    = trim(basename),              &
            &    globProc    = proc,                        &
            &    solver      = solver                       )
        end if  

      end if

    end do

    if (me%tracking%control%nActive > 0) then
      allocate(me%mesh(me%tracking%control%nActive))
      allocate(me%varsys(me%tracking%control%nActive))
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
    &                                  var_degree, var_space, simControl,  &
    &                                  time                                )
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

    !> Provide a time for the current data set to write in tracking.
    !!
    !! This only is respected if no simControl is provided. If simControl
    !! is present the time information from it will be used instead.
    type(tem_time_type), intent(in), optional       :: time
    ! -------------------------------------------------------------------- !
    character(len=pathLen) :: basename
    type(tem_time_type)    :: loctime
    type(tem_varsys_type)  :: sampled_vars
    type(treelmesh_type)   :: sampled_mesh
    integer :: iTrack, iConfig
    integer :: iVar
    ! -------------------------------------------------------------------- !

    call tem_time_reset(loctime)

    if (present(simControl)) then
      loctime = simControl%now
    else
      loctime = time
    end if

    if (me%sampling%max_nlevels == 0) then
      if (present(simControl)) then
        ! No subsampling to be done, call the regular tracker, and leave
        ! the routine.
        call tem_tracker( track      = me%tracking, &
          &               simControl = simControl,  &
          &               varSys     = varsys,      &
          &               tree       = mesh         )
      else
        do iTrack=1,me%tracking%control%nActive
          iConfig = me%tracking%instance(iTrack)%pntConfig
          call hvs_output_open( out_file   = me%tracking%instance(iTrack) &
            &                                           %output_file,     &
            &                   use_iter   = me%tracking%config(iConfig)  &
            &                                  %output_config%vtk         &
            &                                  %iter_filename,            &
            &                   mesh       = mesh,                        &
            &                   varsys     = varsys,                      &
            &                   time       = time                         )

          call hvs_output_write( out_file = me%tracking%instance(iTrack) &
            &                                          %output_file,     &
            &                    varsys   = varsys,                      &
            &                    mesh     = mesh                         )

          call hvs_output_close( out_file = me%tracking%instance(iTrack) &
            &                                          %output_file,     &
            &                    varSys   = varsys,                      &
            &                    mesh     = mesh                         )
        end do
      end if
      RETURN
    end if

    do iTrack=1,me%tracking%control%nActive
      iConfig = me%tracking%instance(iTrack)%pntConfig
      if (present(simControl)) then
        ! If a simControl is provided, check each tracking on whether it is to
        ! be written. Without a simControl, we will write all trackings
        ! unconditionally.
        if ( .not. tem_tracking_has_triggered(                     &
          &            timeControl = me%tracking%config(iConfig)   &
          &                                     %timeControl,      &
          &            proc        = me%tracking%instance(iTrack)  &
          &                                     %output_file%proc, &
          &            simControl  = simControl         ) ) CYCLE
      end if

      if (.not. me%tracking%config(iConfig)%output_config%useGetPoint) then
        ! Only perform subsampling if not using get_point anyway.
        call ply_sample_data( me          = me%sampling,                  &
          &                   orig_mesh   = mesh,                         &
          &                   orig_bcs    = bc,                           &
          &                   varsys      = varsys,                       &
          &                   var_degree  = var_degree,                   &
          &                   var_space   = var_space,                    &
          &                   ndims       = me%ndims,                     &
          &                   trackInst   = me%tracking%instance(iTrack), &
          &                   trackConfig = me%tracking%config(iConfig),  &
          &                   time        = loctime,                      &
          &                   new_mesh    = sampled_mesh,                 &
          &                   resvars     = sampled_vars                  )

        ! initialize output
        basename = trim(me%tracking%config(iConfig)%prefix) &
          &        // trim(me%tracking%config(iConfig)%label)
        call hvs_output_init(                                               &
          &             out_file    = me%tracking%instance(iTrack)          &
          &                             %output_file,                       &
          &             out_config  = me%tracking%config(iConfig)           &
          &                             %output_config,                     &
          &             tree        = sampled_mesh,                         &
          &             varSys      = sampled_vars,                         &
          &             geometry    = me%tracking%config(iConfig)%geometry, &
          &             basename    = trim(basename),                       &
          &             globProc    = proc,                                 &
          &             solver      = solver                                )

        call hvs_output_open( out_file = me%tracking%instance(iTrack)        &
          &                                         %output_file,            &
          &                   use_iter = me%tracking%config(iConfig)         &
          &                                %output_config%vtk%iter_filename, &
          &                   mesh     = sampled_mesh,                       &
          &                   varsys   = sampled_vars,                       &
          &                   time     = loctime                             )

        ! Fill output files with data.
        call hvs_output_write( out_file = me%tracking%instance(iTrack) &
          &                                          %output_file,     &
          &                    varsys   = sampled_vars,                &
          &                    mesh     = sampled_mesh                 )

        call hvs_output_close( out_file = me%tracking%instance(iTrack) &
          &                                          %output_file,     &
          &                    varSys   = sampled_vars,                &
          &                    mesh     = sampled_mesh                 )

        do ivar=1,sampled_vars%method%nVals
          call ply_sampling_free_methodData(sampled_vars%method%val(iVar))
        end do
        call tem_empty_varSys(sampled_vars)
        call free_treelmesh(sampled_mesh)

      else

        if (me%tracking%instance(iTrack)%subTree%useGlobalMesh) then
          call hvs_output_open( out_file = me%tracking%instance(iTrack) &
            &                                         %output_file,     &
            &                   use_iter = me%tracking%config(iConfig)  &
            &                                %output_config             &
            &                                %vtk%iter_filename,        &
            &                   mesh     = mesh,                        &
            &                   varsys   = varSys,                      &
            &                   time     = loctime                      )

          ! Fill output files with data.
          call hvs_output_write( out_file = me%tracking%instance(iTrack) &
            &                                          %output_file,     &
            &                    varsys   = varSys,                      &
            &                    mesh     = mesh                         )

          call hvs_output_close( out_file = me%tracking%instance(iTrack) &
            &                                          %output_file,     &
            &                    varSys   = varSys,                      &
            &                    mesh     = mesh                         )

        else
          call hvs_output_open( out_file = me%tracking%instance(iTrack) &
            &                                         %output_file,     &
            &                   use_iter = me%tracking%config(iConfig)  &
            &                                %output_config             &
            &                                %vtk%iter_filename,        &
            &                   mesh     = mesh,                        &
            &                   subtree  = me%tracking%instance(iTrack) &
            &                                         %subtree,         &
            &                   varsys   = varSys,                      &
            &                   time     = loctime                      )

          ! Fill output files with data.
          call hvs_output_write( out_file = me%tracking%instance(iTrack) &
            &                                          %output_file,     &
            &                    varsys   = varSys,                      &
            &                    subtree  = me%tracking%instance(iTrack) &
            &                                          %subtree,         &
            &                    mesh     = mesh                         )

          call hvs_output_close( out_file = me%tracking%instance(iTrack) &
            &                                          %output_file,     &
            &                    varSys   = varSys,                      &
            &                    mesh     = mesh,                        &
            &                    subtree  = me%tracking%instance(iTrack) &
            &                                          %subtree          )
        end if
      end if

    end do

  end subroutine ply_sampled_track_output


end module ply_sampled_tracking_module
