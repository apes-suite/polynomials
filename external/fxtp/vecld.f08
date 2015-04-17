module fxt_vecld_binding 
implicit none
use, intrinsic :: iso_c_binding
interface 
  subroutine fxt_vecld_del(x) bind(c)
    use, intrinsic :: iso_c_binding
    type(c_ptr), value :: x
  end subroutine fxt_vecld_del
end interface
end module fxt_vecld_binding 
