?? include 'arrayMacros.inc'
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
module ply_dynarray_project_module

  use env_module, only: long_k, rk, minLength, labelLen, zeroLength
  use aotus_module,              only: flu_State, aot_get_val
  use tem_logging_module,        only: logUnit
  use tem_aux_module,            only: tem_abort

!!  use ply_poly_project_module, only: ply_poly_project_type,      &
!!    &                                assignment(=),              &
!!    &                                operator(==), operator(>=), &
!!    &                                operator(/=), operator(<),  &
!!    &                                operator(<=), operator(>)
  use ply_prj_header_module, only: ply_prj_header_load,    &
    &                                  ply_prj_header_type,    &
    &                                  assignment(=),              &
    &                                  operator(==), operator(>=), &
    &                                  operator(/=), operator(<),  &
    &                                  operator(<=), operator(>)


  implicit none

  private

  !> Projection definition.
  type ply_prj_init_type
    !> Polynomial basis type.
    !!
    !! 3D Monomials have the form x^i * y^j * z^k
    !! - Q_space: quadratic polynomial space (i,j,k) <= maxPolyDegree
    !! - P_space: polynomial space i+j+k <= maxPolyDegree
    integer                       :: basisType
    !> The maximal polynomial degree per spatial direction.
    integer                       :: maxPolyDegree
    !> projection header consits of general information like which kind
    !! of projection is used
    type(ply_prj_header_type) :: header
  end type ply_prj_init_type

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

  interface assignment(=)
    module procedure Copy_ply_prj_init
  end interface


?? copy :: DA_decltxt( projection, type(ply_prj_init_type))

  public :: init, append, dyn_projectionArray_type
  public :: operator(==), operator(/=), operator(<), operator(<=)
  public :: operator(>), operator(>=)
  public :: ply_prj_init_define, ply_prj_init_type
  public :: ply_fill_dynProjectArray


contains

?? copy :: DA_impltxt( projection, type(ply_prj_init_type) , type(ply_prj_init_type))


  !**************************************************************************!
  subroutine Copy_ply_prj_init(left,right)
    !------------------------------------------------------------------------!
    !> fpt to copy to
    type(ply_prj_init_type), intent(out) :: left
    !> fpt to copy from
    type(ply_prj_init_type), intent(in) :: right
    !------------------------------------------------------------------------!

    left%header = right%header

    left%maxPolyDegree = right%maxPolyDegree
    left%basisType = right%basisType

  end subroutine copy_ply_prj_init
  !**************************************************************************!


  !**************************************************************************!
  !> Define a projection, without filling its body.
  subroutine ply_prj_init_define(me, header, maxPolyDegree, basisType)
    !------------------------------------------------------------------------!
    type(ply_prj_init_type), intent(inout) :: me
    type(ply_prj_header_type), intent(in) :: header
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: basisType
    !------------------------------------------------------------------------!

    me%header = header
    me%maxPolyDegree = maxPolyDegree
    me%basisType = basisType

  end subroutine ply_prj_init_define
  !****************************************************************************!


  !****************************************************************************!
  !> Load settings to describe a projection method from a Lua table.
  subroutine ply_fill_dynProjectArray( proj_pos, dyn_projectionArray,basisType, &
    &                                  maxPolyDegree, conf,parent       )
  !----------------------------------------------------------------------------!  
    integer, intent(inout) :: proj_pos
    type(dyn_ProjectionArray_type), intent(inout) :: dyn_projectionArray
    type(flu_State), intent(in)   :: conf
    integer, intent(in)           :: basisType
    integer, intent(in)           :: maxPolyDegree
    integer, intent(in)           :: parent
    !--------------------------------------------------------------------------!
    type(ply_prj_header_type) :: header
    type(ply_prj_init_type) :: proj_init

    !--------------------------------------------------------------------------!  
    ! check if it is general projection (no parent is present) or individuall
    ! projection
    call ply_prj_header_load(me   = header,  &
      &                          conf = conf,    &
      &                          parent = parent )
    ! --> projection table exist, we build up the init_type and check if
    ! it is already in the DA
    ! define the init_type
    call  ply_prj_init_define( me= proj_init,                 &
      &                            header = header,               &
      &                            maxPolyDegree = maxPolyDegree, &
      &                            basisType=basisType            )

    ! Call to the unqiue list -> dynamic array for the projection method
    ! chack if the list needs to be appended...store the position
    write(logUnit(5),*) 'The projection_initialization type is'
    write(logUnit(5),*) ' kind=', proj_init%header%kind 
    write(logUnit(5),*) ' degree=', proj_init%maxPolydegree 
    write(logUnit(5),*) ' basisType=', proj_init%basisType
    write(logUnit(5),*)  'Now append the dynamic list for poly projection'
    call append( me = dyn_projectionArray,  &
      &          val = proj_init,           & 
      &          pos = proj_pos             )

  end subroutine ply_fill_dynProjectArray
  !****************************************************************************!


  !****************************************************************************!
  !> This function provides the test for equality of two projections.
  !!
  !! Two projections are considered to be equal, if their kind, nodes_kind,
  !! maxPolyDegree and oversampling are equal.
  function isEqual(left, right) result(equality)
    !---------------------------------------------------------------------------
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is equal??
    logical :: equality
    !---------------------------------------------------------------------------

    equality = ( left%header == right%header )                      &
      &         .and. ( left%maxPolyDegree == right%maxPolyDegree ) &
      &         .and. ( left%basisType == right%basisType )

  end function isEqual
  !****************************************************************************!


  !****************************************************************************!
  !> This function provides the test for unequality of two projections.
  !!
  !! Two projections are considered to be unequal, if their kind, nodes_kind,
  !! maxpolydegree or factor are not equal.
  function isUnequal(left, right) result(unequality)
    !---------------------------------------------------------------------------
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is unequal??
    logical :: unequality
    !---------------------------------------------------------------------------

    unequality = ( left%header /= right%header )                    &
      &          .or. ( left%maxPolyDegree /= right%maxPolyDegree ) &
      &          .or. ( left%basisType /= right%basisType )

  end function isUnequal
  !****************************************************************************!


  !****************************************************************************!
  !> This function provides a < comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  function isSmaller(left, right) result(small)
    !---------------------------------------------------------------------------
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is smaller??
    logical :: small
    !---------------------------------------------------------------------------

    small = .false.
    if (left%header < right%header) then
      small = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree < right%maxPolyDegree) then
          small = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            small = (left%basisType < right%basisType)
          end if
        end if
      end if
    end if

  end function isSmaller
  !****************************************************************************!


  !****************************************************************************!
  !> This function provides a <= comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  function isSmallerOrEqual(left, right) result(small)
    !---------------------------------------------------------------------------
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is smaller??
    logical :: small
    !---------------------------------------------------------------------------

    small = .false.
    if (left%header < right%header) then
      small = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree < right%maxPolyDegree) then
          small = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            small = (left%basisType <= right%basisType)
          end if
        end if
      end if
    end if

  end function isSmallerOrEqual
  !****************************************************************************!


  !****************************************************************************!
  !> This function provides a > comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  function isGreater(left, right) result(great)
    !---------------------------------------------------------------------------
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is greater??
    logical :: great
    !---------------------------------------------------------------------------

    great = .false.
    if (left%header > right%header) then
      great = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree > right%maxPolyDegree) then
          great = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            great = (left%basisType > right%basisType)
          end if
        end if
      end if
    end if

  end function isGreater
  !****************************************************************************!


  !****************************************************************************!
  !> This function provides a >= comparison of two projections.
  !!
  !! Sorting of projections is given by maxPolyDegree, kind, nodes_kind and
  !! last by factor.
  function isGreaterOrEqual(left, right) result(great)
    !---------------------------------------------------------------------------
    !> projection to compare
    type(ply_prj_init_type), intent(in) :: left
    !> projection to compare against
    type(ply_prj_init_type), intent(in) :: right
    !> is greater??
    logical :: great
    !---------------------------------------------------------------------------

    great = .false.
    if (left%header > right%header) then
      great = .true.
    else
      if (left%header == right%header) then
        if (left%maxPolyDegree > right%maxPolyDegree) then
          great = .true.
        else
          if (left%maxPolyDegree == right%maxPolyDegree) then
            great = (left%basisType >= right%basisType)
          end if
        end if
      end if
    end if

  end function isGreaterOrEqual
  !****************************************************************************!
end module ply_dynarray_project_module

