!> ply_fpt_header_module
!!
!! This module contains all information for the header for the fpt method.

module ply_fpt_header_module

  use env_module,              only: rk
  use aotus_module,            only: flu_State, aot_get_val
  use aot_out_module,          only: aot_out_type, aot_out_val

  use tem_aux_module,          only: tem_abort
  use tem_logging_module,      only: logUnit
  use tem_compileconf_module,  only: vlen
  use tem_float_module
  use ply_nodes_header_module

  implicit none

  private

  !> The recommended minimal blocksize for double precision.
  integer, public, parameter :: ply_fpt_default_blocksize = 64

  !> The default width of the subblocking of the diagonal calculation of the
  !! fpt projection
  integer, public, parameter :: ply_fpt_default_subblockingWidth = 8

  !> Default number of terms to use in FPT blocks. 18 is recommended for
  !! double precision.
  integer, public, parameter :: ply_fpt_default_approx_terms = 18

  !> Type for the fpt header, stores all information needed to initialize the
  !! fpt method later on
  type ply_fpt_header_type
    type(ply_nodes_header_type) :: nodes_header
    !> In case of nonlinear equations, aliasing occurs if the projections
    !! of the nonlinear terms on the testfunctions are not calculated
    !! accurately enough. To avoid these errors it is possible to
    !! extend the transformation vectors of the FPT with zeros. This
    !! factor determines by how many zeros the modal vector is extended
    !! before transformation. This factor has to be chosen properly with
    !! respect of the type of nonlinearity of your equation.
    real(kind=rk) :: factor = 1.0_rk

    !> The blockisze of the fast bases exchange algorithm from
    !! Legendre to Chebyshev polynomials.
    !! A negative number indicates to use the default blocksize of the
    !! algorithm.
    integer :: blocksize = ply_fpt_default_blocksize

    !> The number of approximation terms to use for blocks apart from the
    !! diagonal.
    !!
    !! This defaults to 18, which is recommended for double precision.
    integer :: approx_terms = ply_fpt_default_approx_terms

    !> The striplen, that should be used for vectorized simultaneous
    !! computations of the matrix operation.
    !!
    !! This defaults to the vlen from the TEM_compileconf_module, it might
    !! be set differently here, as we are dealing with a twodimensional
    !! problem here, and the optimal setting might be different from the code
    !! parts.
    integer :: striplen = vlen

    !> The width of the subblocks used during the unrolled base exchange to
    !! ensure a better cache usage.
    !!
    !! The default is a subblocking width of 8.
    integer :: subblockingWidth = ply_fpt_default_subblockingWidth

    !> Should the oversampling factor be adapted to ensure a power of 2
    !! in the oversampled polynomial?
    !!
    !! If this is true, the factor will be increased to ensure
    !! an oversampled representation with a power of 2.
    !! Default is false.
    logical :: adapt_factor_pow2 = .false.
  end type ply_fpt_header_type

  interface assignment(=)
    module procedure Copy_fpt_header
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

  public :: ply_fpt_header_load, ply_fpt_header_display
  public :: ply_fpt_header_type
  public :: ply_fpt_header_out

  public :: assignment(=)


contains


  ! ************************************************************************ !
  pure subroutine Copy_fpt_header(left,right)
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_Fpt_header_type), intent(out) :: left
    !> fpt to copy from
    type(ply_Fpt_header_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%nodes_header = right%nodes_header
    left%factor = right%factor
    left%blocksize = right%blocksize
    left%adapt_factor_pow2 = right%adapt_factor_pow2
    left%approx_terms = right%approx_terms
    left%striplen = right%striplen
    left%subblockingWidth = right%subblockingWidth

  end subroutine Copy_fpt_header
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_fpt_header_load(me, conf, thandle)
    ! -------------------------------------------------------------------- !
    type(ply_fpt_header_type), intent(out) :: me
    type(flu_State), intent(inout) :: conf
    integer, intent(in) :: thandle
    ! -------------------------------------------------------------------- !
    integer :: iError
    logical :: fftMultiThread
    ! -------------------------------------------------------------------- !
    ! check for fpt lib

    ! for fpt chebyshev nodes are used
    me%nodes_header%nodes_kind = 'chebyshev'

    ! fill up the fpt_header
    call aot_get_val( L       = conf,                      &
      &               thandle = thandle,                   &
      &               key     = 'blocksize',               &
      &               val     = me%blocksize ,             &
      &               default = ply_fpt_default_blocksize, &
      &               ErrCode = iError                     )

    if (me%blocksize <= 0) then
      write(logUnit(1),*) 'ERROR in loading projection: blocksize for FPT has' &
        & // ' to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%blocksize
      call tem_abort()
    end if

    call aot_get_val( L       = conf,      &
      &               thandle = thandle,   &
      &               key     = 'factor',  &
      &               val     = me%factor, &
      &               default = 1.0_rk,    &
      &               ErrCode = iError     )

    if (me%factor <= 0) then
      write(logUnit(1),*) 'ERROR in loading projection: factor for projection' &
        & // ' has to be larger than 0!'
      write(logUnit(1),*) 'But it is set to ', me%factor
      call tem_abort()
    end if

    call aot_get_val( L       = conf,            &
      &               thandle = thandle,         &
      &               key     = 'approx_terms',  &
      &               val     = me%approx_terms, &
      &               ErrCode = iError,          &
      &               default = 18               )

    call aot_get_val( L       = conf,        &
      &               thandle = thandle,     &
      &               key     = 'striplen',  &
      &               val     = me%striplen, &
      &               ErrCode = iError,      &
      &               default = vlen         )

    call aot_get_val( L       = conf,                            &
      &               thandle = thandle,                         &
      &               key     = 'subblockingWidth',              &
      &               val     = me%subblockingWidth,             &
      &               ErrCode = iError,                          &
      &               default = ply_fpt_default_subblockingWidth )

    write(logUnit(1), *) 'subblockingWidth = ', me%subblockingWidth

    call aot_get_val( L       = conf,                 &
      &               thandle = thandle,              &
      &               key     = 'adapt_factor_pow2',  &
      &               val     = me%adapt_factor_pow2, &
      &               ErrCode = iError,               &
      &               default = .false.               )

    ! check for lobatto Points
    call aot_get_val( L       = conf,                          &
      &               thandle = thandle,                       &
      &               key     = 'lobattoPoints',               &
      &               val     = me%nodes_header%lobattoPoints, &
      &               ErrCode = iError,                        &
      &               default = .false.                        )

    ! check for the multi-threading of the FFTW version
    call aot_get_val( L       = conf,             &
      &               thandle = thandle,          &
      &               key     = 'fftMultiThread', &
      &               val     = fftMultiThread,   &
      &               ErrCode = iError,           &
      &               default = .false.           )

   write(logUnit(1),*) ' * using fftMultithread = ', fftMultiThread

  end subroutine ply_fpt_header_load
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Write FPT settings into a Lua table.
  subroutine ply_fpt_header_out(me, conf)
    ! -------------------------------------------------------------------- !
    type(ply_fpt_header_type), intent(in) :: me
    type(aot_out_type), intent(inout) :: conf
    ! -------------------------------------------------------------------- !

    ! fill up the fpt_header
    call aot_out_val( put_conf = conf,        &
       &              vname    = 'blocksize', &
       &              val      = me%blocksize )

    call aot_out_val( put_conf = conf,     &
       &              vname    = 'factor', &
       &              val      = me%factor )

    call aot_out_val( put_conf = conf,           &
       &              vname    = 'approx_terms', &
       &              val      = me%approx_terms )

    call aot_out_val( put_conf = conf,       &
       &              vname    = 'striplen', &
       &              val      = me%striplen )

    call aot_out_val( put_conf = conf,               &
       &              vname    = 'subblockingWidth', &
       &              val      = me%subblockingWidth )

    call aot_out_val( put_conf = conf,                &
       &              vname    = 'adapt_factor_pow2', &
       &              val      = me%adapt_factor_pow2 )

    call aot_out_val( put_conf = conf,                         &
       &              vname    = 'lobattoPoints',              &
       &              val      = me%nodes_header%lobattoPoints )

  end subroutine ply_fpt_header_out
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine ply_fpt_header_display (me)
    ! -------------------------------------------------------------------- !
    type(ply_fpt_header_type), intent(in) :: me
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*) ' Using fast polynomial transforms for projection.'
    write(logUnit(1),*)
    write(logUnit(1),*) ' * Kind of projection method = fpt'
    write(logUnit(1),*) ' * Dealising factor to use in projection = ', &
      &                 me%factor
    write(logUnit(1),*) ' * Adapt factor to ensure power of 2 order = ', &
      &                 me%adapt_factor_pow2
    write(logUnit(3),*) '     This setting is only relevant for' &
      &                 //' polynomial degrees < 2*blocksize.'
    write(logUnit(1),*) ' * using LobattoPoints =', &
      &                 me%nodes_header%lobattoPoints
!>TODO VK: writing the parameter approx_terms or not?
    write(logUnit(1),*) ' * Block approximation:'
    write(logUnit(1),*) '   * Blocksize for FPT =', me%blocksize
    write(logUnit(1),*) '   * Number of approximation terms = ', me%approx_terms
    write(logUnit(1),*) '   * Strip length = ', me%striplen
    write(logUnit(1),*) '   * Subblocking width = ', me%subblockingWidth
    write(logUnit(1),*) ''
  end subroutine ply_fpt_header_display
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for equality of two projections.
  !!
  !! Two fpt header are considered to be equal, if their  node_header,
  !! fpt_blocksize or the factor are equal.
  pure function isEqual(left, right) result(equality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is equal??
    logical :: equality
    ! -------------------------------------------------------------------- !

    equality = ( left%nodes_header == right%nodes_header )           &
      & .and. ( left%factor .feq. right%factor )                     &
      & .and. ( left%blocksize == right%blocksize )                  &
      & .and. ( left%approx_terms == right%approx_terms )            &
      & .and. ( left%striplen == right%striplen )                    &
      & .and. ( left%subblockingWidth == right%subblockingWidth )    &
      & .and. ( left%adapt_factor_pow2 .eqv. right%adapt_factor_pow2 )

  end function isEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides the test for unequality of two projections.
  !!
  !! Two fpt header are considered to be unequal, if their  node_header,
  !! fpt_blocksize or the factor are not equal.
  pure function isUnequal(left, right) result(unequality)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    ! -------------------------------------------------------------------- !

    unequality = ( left%nodes_header /= right%nodes_header )         &
      & .or. ( left%factor .fne. right%factor )                      &
      & .or. ( left%blocksize /= right%blocksize )                   &
      & .or. ( left%approx_terms /= right%approx_terms )             &
      & .or. ( left%striplen /= right%striplen )                     &
      & .or. ( left%subblockingWidth /= right%subblockingWidth )     &
      & .or. ( left%adapt_factor_pow2 .neqv. right%adapt_factor_pow2 )

  end function isUnequal
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a < comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isSmaller(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    small = .false.

    if (left%adapt_factor_pow2) then
      left_log = 1
    else
      left_log = 0
    end if
    if (right%adapt_factor_pow2) then
      right_log = 1
    else
      right_log = 0
    end if

    small = left%nodes_header < right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      small = left%factor < right%factor
      if (left%factor .feq. right%factor) then
        small = left%blocksize < right%blocksize
        if (left%blocksize == right%blocksize) then
          small = left%approx_terms < right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            small = left%striplen < right%striplen
            if (left%striplen == right%striplen) then
              small = left%subblockingWidth < right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                small = (left_log < right_log)
              end if
            end if
          end if
        end if
      end if
    end if

  end function isSmaller
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a <= comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isSmallerOrEqual(left, right) result(small)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is smaller??
    logical :: small
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    small = .false.

    if (left%adapt_factor_pow2) then
      left_log=1
    else
      left_log=0
    end if
    if (right%adapt_factor_pow2) then
      right_log=1
    else
      right_log=0
    end if

    small = left%nodes_header < right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      small = left%factor < right%factor
      if (left%factor .feq. right%factor) then
        small = left%blocksize < right%blocksize
        if (left%blocksize == right%blocksize) then
          small = left%approx_terms < right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            small = left%striplen < right%striplen
            if (left%striplen == right%striplen) then
              small = left%subblockingWidth < right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                small = (left_log <= right_log)
              end if
            end if
          end if
        end if
      end if
    end if

  end function isSmallerOrEqual
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a > comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isGreater(left, right) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    great = .false.

    if (left%adapt_factor_pow2) then
      left_log=1
    else
      left_log=0
    end if
    if (right%adapt_factor_pow2) then
      right_log=1
    else
      right_log=0
    end if

    great = left%nodes_header > right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      great = left%factor > right%factor
      if (left%factor .feq. right%factor) then
        great = left%blocksize > right%blocksize
        if (left%blocksize == right%blocksize) then
          great = left%approx_terms > right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            great = left%striplen > right%striplen
            if (left%striplen == right%striplen) then
              great = left%subblockingWidth > right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                great = (left_log > right_log)
              end if
            end if
          end if
        end if
      end if
    end if

  end function isGreater
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> This function provides a >= comparison of two projections.
  !!
  !! Sorting of fpt header is given by node_header, fpt_blocksize and
  !! last by factor.
  pure function isGreaterOrEqual(left, right) result(great)
    ! -------------------------------------------------------------------- !
    !> projection to compare
    type(ply_fpt_header_type), intent(in) :: left
    !> projection to compare against
    type(ply_fpt_header_type), intent(in) :: right
    !> is greater??
    logical :: great
    ! -------------------------------------------------------------------- !
    !> help variables
    integer :: left_log, right_log
    ! -------------------------------------------------------------------- !

    great = .false.

    if (left%adapt_factor_pow2) then
      left_log=1
    else
      left_log=0
    end if
    if (right%adapt_factor_pow2) then
      right_log=1
    else
      right_log=0
    end if

    great = left%nodes_header > right%nodes_header
    if (left%nodes_header == right%nodes_header) then
      great = left%factor > right%factor
      if (left%factor .feq. right%factor) then
        great = left%blocksize > right%blocksize
        if (left%blocksize == right%blocksize) then
          great = left%approx_terms > right%approx_terms
          if ( left%approx_terms ==right%approx_terms) then
            great = left%striplen > right%striplen
            if (left%striplen == right%striplen) then
              great = left%subblockingWidth > right%subblockingWidth
              if (left%subblockingWidth == right%subblockingWidth) then
                great = (left_log >= right_log)
              end if
            end if
          end if
        end if
      end if
    end if
  end function isGreaterOrEqual
  ! ************************************************************************ !

end module ply_fpt_header_module
