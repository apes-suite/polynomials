module fxt_fwrap
  use, intrinsic :: iso_c_binding
  implicit none

  type fxtf_flptld
    type(c_ptr) :: flpt
    type(c_ptr) :: w
  end type fxtf_flptld  

  interface
    subroutine fxtf_flptld_evl(v, vn, fplt, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: fplt
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: w
    end subroutine fxtf_flptld_evl

    subroutine fxtf_flptld_exp(u, un, fplt, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: fplt
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: w
    end subroutine fxtf_flptld_exp

    subroutine fxtf_faltld_evl(v, vn, falt, m, u, un, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: falt
      integer(kind=c_long) :: m
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: w
    end subroutine fxtf_faltld_evl

    subroutine fxtf_faltld_exp(u, un, falt, m, v, vn, w) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: u
      integer(kind=c_int) :: un
      type(c_ptr), value :: falt
      integer(kind=c_long) :: m
      type(c_ptr), value :: v
      integer(kind=c_int) :: vn
      type(c_ptr), value :: w
    end subroutine fxtf_faltld_exp

  end interface
 
  contains 
  
  subroutine fxtf_flptld_init(p, n, prec, flpt)
    ! Initialize a flptld data structure
    ! Create a new vector (working vector)
       use, intrinsic :: iso_c_binding
       use :: fxt_fif
       integer, parameter :: rk = selected_real_kind(15)
       integer(8) :: p
       integer(8):: n
       integer(8):: wsize
       real(kind=rk) :: prec
       ! real(kind=rk), dimension(:), allocatable :: w
       ! type(c_ptr), value :: w
       ! type(c_ptr), value :: flpt
       type(fxtf_flptld) :: flpt
       integer :: status
       
       flpt%flpt = fxt_flptld_init(p, n, prec)  
       wsize = fxt_flptld_wsize(flpt%flpt)
       ! allocate(w(0:wsize-1), stat=status)
       flpt%w = fxt_vecld_new(wsize)
  end subroutine fxtf_flptld_init

end module fxt_fwrap
