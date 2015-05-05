module fxt_vecld
  use, intrinsic :: iso_c_binding
  use fxt_types
  IMPLICIT NONE

interface
! /*** create a new vector ***/
! fxt_vecld* fxt_vecld_new(long size);
  type(c_ptr) function fxt_vecld_new(size) bind (c)
    use, intrinsic :: iso_c_binding
    integer(c_long) :: size
  end function fxt_vecld_new


end interface
end module fxt_vecld

! interface
!   type(c_ptr) function fxt_vecl_new(size) bind (c)
!     use, intrinsic :: iso_c_binding
!     integer(c_long) :: size
!   end function fxt_vecl_new
! end interface
! end module fxt_vecld 
