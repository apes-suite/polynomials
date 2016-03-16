module ply_subresolution_module
  use env_module, only: pathLen, rk, isLittleEndian, long_k, newunit

  use aotus_module, only: flu_State, aot_get_val, aoterr_Fatal, close_config

  use treelmesh_module, only: treelmesh_type
  use tem_aux_module, only: tem_open_distconf, tem_abort
  use tem_color_prop_module, only: tem_color_prop_type
  use tem_comm_env_module, only: tem_comm_env_type
  use tem_logging_module, only: logunit
  use tem_subres_prop_module, only: tem_subres_prop_type, tem_subres_prop_load
  use tem_tools_module, only: upper_to_lower

  use ply_dof_module, only: P_Space, Q_space, posOfModgCoeffPTens2D, &
    &                       posOfModgCoeffPTens, nextModgCoeffPTens2D, &
    &                       nextModgCoeffPTens

  use ply_transfer_module, only: ply_transfer_P_dim, ply_transfer_dofs

  implicit none

  type ply_subresolution_type
    integer :: polydegree = 0
    integer :: basisType
    type(tem_subres_prop_type) :: subres_prop
  end type ply_subresolution_type


contains


  !****************************************************************************!
  !> Subroutine to load subresolution information for a given tree.
  subroutine ply_subresolution_load(me, tree, proc, coloring)
    !--------------------------------------------------------------------------!
    type(ply_subresolution_type), intent(out) :: me
    type(treelmesh_type), intent(in) :: tree
    type(tem_comm_env_type), intent(in) :: proc
    type(tem_color_prop_type), intent(in) :: coloring
    !--------------------------------------------------------------------------!
    character(len=pathLen) :: configfile
    character :: polyspace
    type(flu_State) :: conf
    integer :: iError
    !--------------------------------------------------------------------------!

    configfile = trim(tree%global%dirname)//'subresolution.lua'

    call tem_subres_prop_load( me       = me%subres_prop, &
      &                        tree     = tree,           &
      &                        coloring = coloring        )

    ! Set the polydegree initially to 0 to ensure a proper setting even if
    ! no subresolution property is present.
    me%polydegree = 0

    ! Only need to do anything, if there is actually a subresolution property.
    ! Checking this via the association status of the property header.
    if (associated(me%subres_prop%header)) then

      call tem_open_distconf( L        = conf,             &
        &                     filename = trim(configfile), &
        &                     proc     = proc              )

      call aot_get_val( L       = conf,          &
        &               key     = 'polydegree',  &
        &               val     = me%polydegree, &
        &               ErrCode = iError         )

      if (btest(iError, aoterr_Fatal)) then
        write(logUnit(1),*) &
          &  'FATAL Error occured, while retrieving subresolution polydegree'
        call tem_abort()
      end if

      call aot_get_val( L       = conf,        &
        &               key     = 'polyspace', &
        &               val     = polyspace,   &
        &               default = 'q',         &
        &               ErrCode = iError       )

      select case(upper_to_lower(trim(polyspace)))
      case('q')
        me%basisType = Q_space

      case('p')
        me%basisType = P_space

      case default
        write(logunit(1),*) 'ERROR in subresolution loading!'
        write(logunit(1),*) 'Unknown polyspace ', trim(polyspace)
        write(logUnit(1),*) 'Supported are:'
        write(logUnit(1),*) '* Q (quadratic with i,j,k <= maxDegree)'
        write(logUnit(1),*) '* P (with i+j+k <= maxDegree)'
        write(logUnit(1),*) 'Stopping....'
        call tem_abort()

      end select

      call close_config(conf)

    end if

  end subroutine ply_subresolution_load
  !****************************************************************************!


  !****************************************************************************!
  !> Get the subresolution data for all elements for a given color and in the
  !! requested format.
  subroutine ply_subres_import_color( me, tree, proc, coloring, iColor,        &
    &                                 target_degree, target_space, target_dim, &
    &                                 subresdat )
    !--------------------------------------------------------------------------!
    type(ply_subresolution_type), intent(in) :: me
    type(treelmesh_type), intent(in) :: tree
    type(tem_comm_env_type), intent(in) :: proc
    type(tem_color_prop_type), intent(in) :: coloring
    integer, intent(in) :: iColor
    integer, intent(in) :: target_degree
    integer, intent(in) :: target_space
    integer, intent(in) :: target_dim
    real(kind=rk), allocatable, intent(out) :: subresdat(:,:)
    !--------------------------------------------------------------------------!
    character(len=pathLen) :: datfile
    integer :: target_Dofs
    integer :: in_dofs
    integer :: read_dofs
    integer :: pre_dofs
    integer :: pdim
    real(kind=rk), allocatable :: indat(:)
    real(kind=rk), allocatable :: predat(:)
    character(len=4) :: datext
    integer :: rl
    integer :: fUnit
    integer :: nElems
    integer(kind=long_k) :: offset
    integer :: iElem
    integer :: iStep, iDof
    integer :: tt_X, tt_Y, tt_Z
    integer :: in_X, in_Y, in_Z
    integer :: tt_pos, in_pos
    integer :: tt_off, tt_zoff
    integer :: in_off, in_zoff
    integer :: minOrd
    integer :: in_dim
    integer :: recs_per_elem
    !--------------------------------------------------------------------------!

    ! Seeder always writes three-dimensional data.
    ! Assume an input dimension of 3:
    in_dim = 3

    ! Only need to do anything, if there is actually a subresolution property.
    ! Checking this via the association status of the property header.
    if (associated(me%subres_prop%header)) then

      if (isLittleEndian) then
        datext = '.lsb'
      else
        datext = '.msb'
      end if

      nElems = me%subres_prop%nElems(iColor)
      offset = me%subres_prop%offset(iColor)

      ! Figure out the target polynomial representation.
      select case(target_space)
      case (Q_Space)
        target_dofs = (target_degree+1)**target_dim

      case (P_Space)
        target_dofs = target_degree+1
        do pdim=2,target_dim
          target_dofs = (target_dofs * (target_degree+pdim) ) / pdim
        end do

      end select

      ! We read the complete data per element in most cases.
      recs_per_elem = 1

      ! Figure out input polynomial representation.
      select case(me%basistype)
      case (Q_Space)
        in_dofs = (me%polydegree+1)**3

        if (target_dim < in_dim) then
          read_dofs = (me%polydegree+1)**target_dim
          recs_per_elem = (me%polydegree+1)**(in_dim-target_dim)
        else
          read_dofs = in_dofs
        end if

      case (P_Space)
        in_dofs = me%polydegree+1
        do pdim=2,in_dim
          in_dofs = (in_dofs * (me%polydegree+pdim) ) / pdim
        end do
        read_dofs = in_dofs

      end select

      ! Allocate arrays accordingly.
      allocate(indat(read_dofs))
      allocate(subresdat(target_dofs, nElems))

      inquire(iolength=rl) indat

      datfile = trim(tree%global%dirname) // 'subresdata_' &
        &       // trim(coloring%color_label(iColor)) // datext

      fUnit = newunit()

      open( unit = fUnit, file = datfile, action = 'read', &
        &   access = 'direct', form = 'unformatted',       &
        &   recl = rl, status = 'old'                      )

      minord = min(target_degree+1, me%polydegree+1)

      subresdat = 0.0_rk

      if  ( (me%basistype == P_Space) .and. (target_dim < in_dim) ) then

        pre_dofs = me%polydegree+1
        do pdim=2,target_dim
          pre_dofs = (pre_dofs * (me%polydegree+pdim) ) / pdim
        end do
        allocate(predat(pre_dofs))

        do iElem=1,nElems

          read(fUnit, rec=offset+iElem) indat

          call ply_transfer_P_dim( indat  = indat,        &
            &                      indim  = in_dim,       &
            &                      outdat = predat,       &
            &                      outdim = target_dim,   &
            &                      degree = me%polydegree )
          call ply_transfer_dofs( indat     = predat,              &
            &                     inspace   = me%basistype,        &
            &                     indegree  = me%polydegree,       &
            &                     outdat    = subresdat(:, iElem), &
            &                     outspace  = target_space,        &
            &                     outdegree = target_degree,       &
            &                     ndims     = target_dim           )
        end do

      else

        do iElem=1,nElems

          read(fUnit, rec=(offset+iElem-1)*recs_per_elem+1) indat

          call ply_transfer_dofs( indat     = indat,               &
            &                     inspace   = me%basistype,        &
            &                     indegree  = me%polydegree,       &
            &                     outdat    = subresdat(:, iElem), &
            &                     outspace  = target_space,        &
            &                     outdegree = target_degree,       &
            &                     ndims     = target_dim           )

        end do

      end if

      close(fUnit)

    else

      ! No subresolution data at all, allocate the array with size 0.
      allocate(subresdat(0,0))

    end if

  end subroutine ply_subres_import_color


end module ply_subresolution_module
