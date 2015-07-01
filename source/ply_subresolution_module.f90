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
    integer :: min_dofs
    integer :: pdim
    real(kind=rk), allocatable :: indat(:)
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
    !--------------------------------------------------------------------------!


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

      ! Figure out input polynomial representation.
      select case(me%basistype)
      case (Q_Space)
        in_dofs = (me%polydegree+1)**3

      case (P_Space)
        in_dofs = me%polydegree+1
        do pdim=2,3
          in_dofs = (in_dofs * (me%polydegree+pdim) ) / pdim
        end do

      end select

      ! Allocate arrays accordingly.
      allocate(indat(in_dofs))
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

      do iElem=1,nElems

        read(fUnit, rec=offset+iElem) indat

        select case(target_space)
        case (Q_Space)
          if (me%basistype == Q_Space) then
            ! Both, target and input are Q Polynomials
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              do tt_Y=0,minord-1
                tt_off = tt_Y*(target_degree+1)
                in_off = tt_Y*(me%polydegree+1)
                subresdat(tt_off+1:tt_off+minord, iElem) &
                  &  = indat(in_off+1:in_off+minord)
              end do

            case (3)
              do tt_Z=0,minord-1
                tt_zoff = tt_Z*(target_degree+1)**2
                in_zoff = tt_Z*(me%polydegree+1)**2
                do tt_Y=0,minord-1
                  tt_off = tt_Y*(target_degree+1) + tt_zoff
                  in_off = tt_Y*(me%polydegree+1) + in_zoff
                  subresdat(tt_off+1:tt_off+minord, iElem) &
                    &  = indat(in_off+1:in_off+minord)
                end do
              end do
            end select

          else

            ! Target is Q, but input is P
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              in_X = 1
              in_Y = 1
              do iDof=1,in_dofs
                in_pos = posOfModgCoeffPTens2D(in_X, in_Y, 0, me%polydegree)
                tt_pos = in_X + (target_degree+1)*(in_Y-1)
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,in_dofs
                  call nextModgCoeffPTens2D(in_X, in_Y, in_Z, me%polydegree)
                  if ((in_X <= minord) .and. (in_Y <= minord)) EXIT
                end do
                if ((in_X > minord) .or. (in_Y > minord)) EXIT
              end do

            case (3)
              in_X = 1
              in_Y = 1
              in_Z = 1
              do iDof=1,in_dofs
                in_pos = posOfModgCoeffPTens(in_X, in_Y, in_Z, me%polydegree)
                tt_pos = in_X + (target_degree+1)*( (in_Y-1) &
                  &                                + (target_degree+1)*(in_Z-1))
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,in_dofs
                  call nextModgCoeffPTens(in_X, in_Y, in_Z, me%polydegree)
                  if ( (in_X <= minord) .and. (in_Y <= minord) &
                    &                   .and. (in_Z <= minord) ) EXIT
                end do
                if ( (in_X > minord) .or. (in_Y > minord) &
                  &                  .or. (in_Z > minord) ) EXIT
              end do
            end select
          end if

        case (P_Space)
          if (me%basistype == Q_Space) then
            ! Target is P, input is Q
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              tt_X = 1
              tt_Y = 1
              do iDof=1,target_dofs
                tt_pos = posOfModgCoeffPTens2D(tt_X, tt_Y, 0, target_degree)
                in_pos = tt_X + (me%polydegree+1)*(tt_Y-1)
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,target_dofs
                  call nextModgCoeffPTens2D(tt_X, tt_Y, tt_Z, target_degree)
                  if ((tt_X <= minord) .and. (tt_Y <= minord)) EXIT
                end do
                if ((tt_X > minord) .or. (tt_Y > minord)) EXIT
              end do

            case (3)
              tt_X = 1
              tt_Y = 1
              tt_Z = 1
              do iDof=1,target_dofs
                tt_pos = posOfModgCoeffPTens(tt_X, tt_Y, tt_Z, target_degree)
                in_pos = tt_X + (me%polydegree+1)*( (tt_Y-1) &
                  &                                + (me%polydegree+1)*(tt_Z-1))
                subresdat(tt_pos, iElem) = indat(in_pos)
                ! Ensure, that next iteration is in the target range
                do istep=iDof,target_dofs
                  call nextModgCoeffPTens(tt_X, tt_Y, tt_Z, target_degree)
                  if ((tt_X <= minord) .and. (tt_Y <= minord) &
                    &                  .and. (tt_Z <= minord) ) EXIT
                end do
                if ((tt_X > minord) .or. (tt_Y > minord) &
                  &                 .or. (tt_Z > minord) ) EXIT
              end do

            end select

          else
            ! Both input and target ar P polynomials
            select case(target_dim)
            case (1)
              subresdat(:minord, iElem) = indat(:minord)

            case (2)
              min_dofs = (minord*(minord+1))/2
              tt_X = 1
              tt_Y = 1
              do iDof=1,min_dofs
                tt_pos = posOfModgCoeffPTens2D(tt_X, tt_Y, 0, target_degree)
                in_pos = posOfModgCoeffPTens2D(tt_X, tt_Y, 0, me%polydegree)
                subresdat(tt_pos, iElem) = indat(in_pos)
                call nextModgCoeffPTens2D(tt_X, tt_Y, tt_Z, minord-1)
              end do

            case (3)
              min_dofs = ( (minord+2)*((minord*(minord+1))/2) ) / 3
              tt_X = 1
              tt_Y = 1
              tt_Z = 1
              do iDof=1,min_dofs
                tt_pos = posOfModgCoeffPTens(tt_X, tt_Y, tt_Z, target_degree)
                in_pos = posOfModgCoeffPTens(tt_X, tt_Y, tt_Z, me%polydegree)
                subresdat(tt_pos, iElem) = indat(in_pos)
                call nextModgCoeffPTens(tt_X, tt_Y, tt_Z, minord-1)
              end do
            end select
          end if

        end select

      end do

      close(fUnit)

    else

      ! No subresolution data at all, allocate the array with size 0.
      allocate(subresdat(0,0))

    end if

  end subroutine ply_subres_import_color


end module ply_subresolution_module
