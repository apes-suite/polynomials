program test_fxtd_n2m2n
  use env_module, only: rk
  use fxt_fwrap, only: fxtf_flptld_type 
  use ply_fxt_module, only:  ply_init_fxt, ply_fxt_type, &
    &                  ply_fxt_n2m_1D, ply_fxt_m2n_1D

  implicit none
  type(ply_fxt_type) :: fxt
  !type(fxtf_flptld_type) :: flpt

  real(kind=rk) :: v_orig(10)
  real(kind=rk), allocatable, target :: u(:,:)
  real(kind=rk), allocatable, target :: v(:,:)
  integer :: nNodes, nModes  

  integer, parameter :: p = 10       ! number of points
  integer, parameter :: n =  9       ! maximal polynomial degree
  real(kind=rk), parameter :: prec = 8*epsilon(1.0_rk) ! Precision for the FMM


  allocate(u(n+1,1))
  allocate(v(p,1))

  call ply_init_fxt(fxt%flpt, p, n, prec)

  call random_number(v_orig)

  write(*,*) 'orig :', v_orig
  v(:,1) = v_orig
  
  ! Test the subroutines m2n and n2m
  nNodes = size(v)
  nModes = size(u)

  ! there ....
  ! transform from physical to wave space
  call ply_fxt_n2m_1D(     fxt        = fxt,      &
    &                      nodal_data = v,        &
    &                      modal_data = u         )

  ! ...and back again
  ! transform from wave to physical space
  call ply_fxt_m2n_1D(     fxt        = fxt,      &
    &                      modal_data = u,        &
    &                      nodal_data = v         )

  write(*,*) 'trafo (after n2m and m2n):', v
  write(*,*) 'Should be the same as orig.'

  if (all(v(:,1) - v_orig < 2.*epsilon(1.0_rk))) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'Data does not match after conversion:'
    write(*,*) 'FAILED'
  end if

end program test_fxtd_n2m2n
