module fxt_types
  use, intrinsic :: iso_c_binding
	IMPLICIT NONE	
  type, bind(c) :: fxt_vecld
    integer(c_long) :: n
    integer(c_long) :: v
  end type fxt_vecld

  type, bind(c) :: fxt_vecl
    integer(c_long) :: n
    integer(c_long) :: v
  end type fxt_vecl

  type, bind(c) :: fxt_flptld
    integer(c_long) :: p
    integer(c_long) :: n
    real(c_double) :: prec
  end type fxt_flptld

  type, bind(c) :: fxt_faltld
    integer(c_long) :: p
    integer(c_long):: n
    type(fxt_vecl) :: mv
    real(c_double) :: prec
    type(fxt_vecld) :: x
    type(fxt_vecld) :: w
    type(falt_alt) :: alt
  end type fxt_faltld

end module fxt_types
