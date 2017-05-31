!> Module contains the header for the l2p projection method
module ply_l2p_header_module

  use env_module,              only: rk
  use aotus_module,            only: flu_State, aot_get_val
  use aot_out_module,          only: aot_out_type, aot_out_val

  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_float_module,        only: operator(.feq.), operator(.fne.)

  use ply_nodes_header_module, only: ply_nodes_header_type,      &
    &                                assignment(=),              &
    &                                operator(==), operator(>=), &
    &                                operator(/=), operator(<),  &
    &                                operator(<=), operator(>)

  implicit none

  private

  !> l2p projection header type, consisting of the node header which give
  !! information about the type and number of points for the projection
  type ply_l2p_header_type
    type(ply_nodes_header_type) :: nodes_header
    real(kind=rk) :: factor
  end type ply_l2p_header_type

  interface assignment(=)
    module procedure Copy_l2p_header
  end interface

  interface operator(==)
    module procedure isEqual
  end interface

  interface operator(/=)
    module procedure isUnequal
  end interface

  interface operator(<)
    module procedure isSmaller
  end interface

  interface operator(<=)
    module procedure isSmallerOrEqual
  end interface

  interface operator(>)
    module procedure isGreater
  end interface

  interface operator(>=)
    module procedure isGreaterOrEqual
  end interface

  public :: operator(==), operator(/=), operator(<), operator(<=)
  public :: operator(>), operator(>=)
  public :: assignment(=)
  public :: ply_l2p_header_type
  public :: ply_l2p_header_load,  ply_l2p_header_display
  public :: ply_l2p_header_out


contains


  ! ************************************************************************ !
  pure subroutine Copy_l2p_header( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_l2p_header_type), intent(out) :: left
    !> fpt to copy from
    type(ply_l2p_header_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%factor = right%factor
    left%nodes_header = right%nodes_header

  end subroutine Copy_l2p_header
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Load settings to describe a projection method from a Lua table.
  subroutine ply_l2p_header_load( me, conf, thandle )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_header_type), intent(out) :: me
    type(flu_State), intent(inout) :: conf
    integer, intent(in) :: thandle
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    ! for l2p gauss-legendre points are used
    me%nodes_header%nodes_kind = 'gauss-legendre'

    ! fill up l2p header
    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'factor',  &
      &               val     = me%factor, &
      &               default = 1.0_rk,    &
      &               ErrCode = iError     )
    call aot_get_val( L       = conf,                          &
      &               thandle = thandle,                       &
      &               key     = 'lobattoPoints',               &
      &               val     = me%nodes_header%lobattoPoints, &
      &               ErrCode = iError,                        &
      &               default = .false.                        )
    if (me%factor .le. 0) then
      write(logUnit(1),*) 'ERROR in loading projection: factor for ' &
        & // 'projection has to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%factor
      call tem_abort()
    end if

  end subroutine ply_l2p_header_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write L2P settings into a Lua table.
  subroutine ply_l2p_header_out( me, conf )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_header_type), intent(in) :: me
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !

    call aot_out_val( put_conf = conf,     &
      &               vname    = 'factor', &
      &               val      = me%factor )

    call aot_out_val( put_conf = conf,                         &
      &               vname    = 'lobattoPoints',              &
      &               val      = me%nodes_header%lobattoPoints )

  end subroutine ply_l2p_header_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_l2p_header_display( me )
    ! -------------------------------------------------------------------- !
    type(ply_l2p_header_type), intent(in) :: me
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) ' * Kind of projection method = l2p'
    write(logUnit(1),*) ' * Factor to use in projection = ', me%factor
    write(logUnit(1),*) ' * using LobattoPoints =', me%nodes_header%lobattoPoints

    write(logUnit(1),*) ' NOT using fast polynomial transforms for projection!'
    if (me%factor < 2.0_rk) then
      write(logUnit(1),*) ''
      write(logUnit(1),*)                                         &
        & '+-----------------------------------------------------+'
      write(logUnit(1),*)                                         &
        & '| WARNING, the oversampling factor is smaller than 2! |'
      write(logUnit(1),*)                                         &
        & '|        this might lead to bad projections!          |'
      write(logUnit(1),*)                                         &
        & '+-----------------------------------------------------+'
    end if

  end subroutine ply_l2p_header_display
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for equality of two projections.
  !!
  !! Two l2p header are considered to be equal, if their node_header, 
  !! and the factor are equal.
  pure function isEqual( left, right ) result( equality )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !

    equality = ( left%nodes_header == right%nodes_header ) &
      & .and. ( left%factor .feq. right%factor )

  end function isEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for unequality of two projections.
  !!
  !! Two l2p header are considered to be unequal, if their node_header, 
  !! or the factor are not equal.
  pure function isUnequal( left, right ) result( unequality )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !

    unequality = ( left%nodes_header /= right%nodes_header ) &
      & .or. ( left%factor .fne. right%factor )

  end function isUnequal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a < comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header and by the factor.
  pure function isSmaller( left, right ) result( small )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%nodes_header < right%nodes_header) then
      small = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        small = (left%factor < right%factor)
      end if
    end if

  end function isSmaller
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a <= comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header, l2p_blocksize and
  !! last by factor.
  pure function isSmallerOrEqual( left, right ) result( small )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !

    small = .false.
    if (left%nodes_header < right%nodes_header) then
      small = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        small = (left%factor <= right%factor)
      end if
    end if

  end function isSmallerOrEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a > comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header, l2p_blocksize and
  !! last by factor.
  pure function isGreater( left, right ) result( great )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%nodes_header > right%nodes_header) then
      great = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        great = (left%factor > right%factor)
      end if
    end if

  end function isGreater
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a >= comparison of two projections.
  !!
  !! Sorting of l2p header is given by node_header, l2p_blocksize and
  !! last by factor.
  pure function isGreaterOrEqual( left, right ) result( great )
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_l2p_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_l2p_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !

    great = .false.
    if (left%nodes_header > right%nodes_header) then
      great = .true.
    else
      if (left%nodes_header == right%nodes_header) then
        great = (left%factor >= right%factor)
      end if
    end if

  end function isGreaterOrEqual
  ! ************************************************************************ !

end module ply_l2p_header_module
