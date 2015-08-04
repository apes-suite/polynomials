program test_fxtd_n2m2n
  use env_module,              only: rk, fin_env
  use tem_logging_module,      only: logUnit
  use ply_dof_module,  only: Q_space
  use ply_dynArray_project_module, only: ply_prj_init_define, &
                                    &    ply_prj_init_type
  use ply_poly_project_module, only: ply_poly_project_fillbody, &
                                &    ply_poly_project_m2n,&
                                &    ply_poly_project_n2m, &
                                &    ply_poly_project_type
  use ply_prj_header_module,        only: ply_prj_header_type
  use tem_general_module,           only: tem_start
  use atl_solver_param_module,      only: solver_param_type
  use fxt_fwrap,                    only: fxtf_flptld_type 
  use ply_fxt_module,               only:  ply_init_fxt, ply_fxt_type, &
    &                                      ply_fxt_n2m_1D, ply_fxt_m2n_1D


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

  real(kind=rk) :: res, newRes
  integer :: power
  type(solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )
  res = 0.0_rk
 !> NA: Commenting it out for the moment
 !NA! allocate(u(n+1,1))
 !NA! allocate(v(p,1))

 !NA! call ply_init_fxt(fxt%flpt, p, n, prec)

 !NA! call random_number(v_orig)

 !NA! write(*,*) 'orig :', v_orig
 !NA! v(:,1) = v_orig
 !NA! 
 !NA! ! Test the subroutines m2n and n2m
 !NA! nNodes = size(v)
 !NA! nModes = size(u)

 !NA! ! there ....
 !NA! ! transform from physical to wave space
 !NA! call ply_fxt_n2m_1D(     fxt        = fxt,      &
 !NA!   &                      nodal_data = v,        &
 !NA!   &                      modal_data = u         )

 !NA! ! ...and back again
 !NA! ! transform from wave to physical space
 !NA! call ply_fxt_m2n_1D(     fxt        = fxt,      &
 !NA!   &                      modal_data = u,        &
 !NA!   &                      nodal_data = v         )

 !NA! write(*,*) 'trafo (after n2m and m2n):', v
 !NA! write(*,*) 'Should be the same as orig.'

 !NA! if (all(v(:,1) - v_orig < 2.*epsilon(1.0_rk))) then
 !NA!   write(*,*) 'PASSED'
 !NA! else
 !NA!   write(*,*) 'Data does not match after conversion:'
 !NA!   write(*,*) 'FAILED'
 !NA! end if



  !Check the 2D fxt projection!
  ! check l2p Q-Space
  do power = 1,7
    write(logUnit(10),*) '---------------------------   CHECKING CHEB->LEG L2P Q-SPACE TRAFO FOR ', 2**power
    call check_fxt_2d(power, newRes)
    if (newRes.gt.res) then
      res = newRes
    end if
    write(logUnit(10),*) '--------------------------- DONE'
  end do

  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if 
  call fin_env()

contains


  subroutine check_fxt_2d(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    !---------------------------------------------
    type(ply_prj_header_type) :: header
    type(ply_poly_project_type)   :: me
    type(ply_prj_init_type)   :: prj_init
    real(kind=rk),allocatable     :: nodal_data(:,:)
    real(kind=rk), allocatable    :: modal_data(:,:)
    real(kind=rk), allocatable    :: oversamp_modal(:,:)
    real(kind=rk), allocatable    :: ref_modes(:)
    type(ply_fxt_type) :: fxt
    !-----------for init
    integer ::basisType, maxdegree, i
    !-----------for oversamp
    integer :: iDegX, iDegY, dof, dofOverSamp

    basisType = Q_space
    maxdegree = power
    header%kind = 'fxt'
    !> todo NA: Check if nodes kind and factor are correctly used for 
    !           fxt here or they should be different
    header%fxt_header%nodes_header%nodes_kind = 'gauss-legendre'
    header%fxt_header%factor = 2.0_rk

    ! define my poly projection type
    call ply_prj_init_define(me= prj_init,             &
      &                          header = header ,         &
      &                          maxPolyDegree= maxdegree, &
      &                          basisType = basistype )
 
    ! fill the projection body according to the header
    call ply_poly_project_fillbody(me, proj_init=prj_init,scheme_dim=2)
  
    allocate(ref_modes(1:me%body_2d%ndofs) ) 
    allocate(modal_data(1:me%body_2d%ndofs,1) ) 
    allocate(oversamp_modal(1:me%body_2d%oversamp_dofs,1) ) 
    allocate(nodal_data(1:me%body_2d%nQuadPoints,1) ) 

    do i=1, me%body_2d%ndofs 
       ref_modes(i)=1.0/real(i, kind=rk)
    end do
    
    modal_data(:,:) = 0.0_rk 
    oversamp_modal(:,:) = 0.0_rk 
    ! oversampling of modes 
    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)
        oversamp_modal(dofOverSamp,1) = ref_modes(dof)
      end do
    end do                  

 
    ! transform from wave to physical space
    call ply_poly_project_m2n(me = me ,                &
      &                       dim = 2 ,                &
      &                       nVars = 1,               &
      &                       nodal_data=nodal_data,   &
      &                       modal_data=oversamp_modal)
  

    
    oversamp_modal(:,:) = 0.0_rk
 
    ! ...and back again
    ! transform from physical to wave space
    call ply_poly_project_n2m(me = me,                  & 
      &                       dim = 2 ,                 &
      &                       nVars = 1,                &
      &                       nodal_data=nodal_data,    &
      &                       modal_data= oversamp_modal)

    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)
        modal_data(dof,:) = oversamp_modal(dofOverSamp,:)
      end do
    end do                  

 
    res= maxval( abs(ref_modes-modal_data(:,1)) ) 

  end subroutine  check_fxt_2d

end program test_fxtd_n2m2n
