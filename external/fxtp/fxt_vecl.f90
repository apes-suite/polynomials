module fxt_vecl
  use, intrinsic :: iso_c_binding
  use fxt_types
  IMPLICIT NONE 

interface
  type(c_ptr) function fxt_vecl_new(size) bind (c)
    use, intrinsic :: iso_c_binding
    integer(c_long) :: size
  end function fxt_vecl_new
end interface    
end module fxt_vecl
