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
  use tem_general_module,           only: tem_start, tem_general_type
  use fxt_fwrap,                    only: fxtf_flptld_type, fxtf_flptld_init
  use ply_fxt_module,               only: ply_fxt_type, &
    &                                     ply_fxt_n2m_1D, ply_fxt_m2n_1D
  use ply_oversample_module,        only: ply_convert2oversample, &
    &                                     ply_convertFromOversample


implicit none


  real(kind=rk) :: res1, newRes1
  real(kind=rk) :: res2, newRes2
  real(kind=rk) :: res3, newRes3
  integer :: power
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Polynomials unit test', &
    &            version  = 'utest',                 &
    &            general  = general                  )
  res1 = 0.0_rk
  res2 = 0.0_rk
  res3 = 0.0_rk



  !--------- CHECK FXT 1D -----------------------------------!
  call check_fxt_1d(res1)



  !--------- CHECK FXT 2D -----------------------------------!
  do power = 1,7
    write(logUnit(10),*) '--------------------------- &
    &              CHECKING FXT TRANSFORMATION FOR ', 2**power
    call check_fxt_2d(power, newRes2)
    if (newRes2.gt.res2) then
      res2 = newRes2
    end if
    write(logUnit(10),*) '--------------------------- DONE'
  end do



  !--------- CHECK FXT 3D -----------------------------------!
  do power = 1,7
    write(logUnit(10),*) '---------------------------   CHECKING  &
    &                     FXT TRANSFORMATION  FOR ', 2**power
    call check_fxt_3d (power, newRes3)
    if (newRes3.gt.res3) then
      res3 = newRes3
    end if
    write(logUnit(10),*) '--------------------------- DONE'
  end do


  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  write(*,*) res1, res2, res3
  if(max(res1,res2,res3) .lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
 
  call fin_env()

contains

  subroutine check_fxt_1d(res)
    real(kind=rk), intent(out) :: res
    !---------------------------------------------

    type(ply_fxt_type), target :: fxt
    real(kind=rk) :: v_orig(10)
    real(kind=rk), allocatable :: u(:,:)
    real(kind=rk), allocatable :: v(:,:)
    integer :: nNodes, nModes  
    integer, parameter :: nPoints = 10       ! number of points
    integer, parameter :: maxDegree =  9       ! maximal polynomial degree
    real(kind=rk), parameter :: prec = 8*epsilon(1.0_rk) ! Precision for the FMM

    allocate(u(maxDegree + 1,1))
    allocate(v(nPoints,1))

    call fxtf_flptld_init(flpt    = fxt%flpt,  &
      &                   degree  = maxDegree, &
      &                   nPoints = nPoints,   &
      &                   prec    = prec       )


    call random_number(v_orig)
    
    write(logunit(10),*) 'orig :', v_orig
    v(:,1) = v_orig
    
    ! Test the subroutines m2n and n2m
    nNodes = size(v,1)
    nModes = size(u,1)



    ! there ....
    ! transform from physical to wave space
    call ply_fxt_n2m_1D(  fxt              = fxt,      &
      &                   nodal_data       = v(:,1),   &
      &                   modal_data       = u(:,1),   &
      &                   nNodes           = nNodes,   &
      &                   nModes           = nModes,   &
      &                   oversamp_degree  = maxDegree )

    v = 0.0
    write(logunit(10),*) "modal values = "
    write(logunit(10),*) u 

    ! ...and back again
    ! transform from wave to physical space
    call ply_fxt_m2n_1D(     fxt        = fxt,      &
      &                      modal_data = u(:,1),   &
      &                      nodal_data = v(:,1),   &
      &                      nNodes     = nNodes,   &
      &                      nModes     = nModes,   &
      &                 oversamp_degree = nModes    )

    write(logunit(10),*) 'trafo (after n2m and m2n):', v
    write(logunit(10),*) 'Should be the same as orig.'

    res = maxval(abs(v(:,1) - v_orig) ) 

  end subroutine check_fxt_1d


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
    
    oversamp_modal(:,:) = 0.0_rk 
    modal_data(:,1) = ref_modes
    write(*,*) "modal_data = "
    write(*,*) modal_data
    write(*,*) "------------------------------------------------------"
    ! oversampling of modes 
    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)
        oversamp_modal(dofOverSamp,1) = ref_modes(dof)
      end do
    end do                  
    write(*,*) " modal_data (in oversampled space) = "
    write(*,*) oversamp_modal
    write(*,*) "------------------------------------------------------"

 
    ! transform from wave to physical space
    call ply_poly_project_m2n(me = me ,                &
      &                       dim = 2 ,                &
      &                       nVars = 1,               &
      &                       nodal_data=nodal_data,   &
      &                       modal_data=oversamp_modal)
  
    write(*,*) "converted to nodal_data (in oversampled space) = "
    write(*,*) nodal_data
    write(*,*) "------------------------------------------------------"
    
    oversamp_modal(:,:) = 0.0_rk
 
    ! ...and back again
    ! transform from physical to wave space
    call ply_poly_project_n2m(me = me,                  & 
      &                       dim = 2 ,                 &
      &                       nVars = 1,                &
      &                       nodal_data=nodal_data,    &
      &                       modal_data= oversamp_modal)

    write(*,*) "converting it  back to modal data (in oversampled space) = "
    write(*,*) oversamp_modal
    write(*,*) "------------------------------------------------------"

    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)
        modal_data(dof,:) = oversamp_modal(dofOverSamp,:)
      end do
    end do                  

    write(*,*) "Final  modal_data = "
    write(*,*) modal_data
    write(*,*) "------------------------------------------------------"
 
    res= maxval( abs(ref_modes-modal_data(:,1)) ) 

  end subroutine  check_fxt_2d

  subroutine check_fxt_3D(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    !---------------------------------------------
    type(ply_poly_project_type)   :: me
    type(ply_prj_init_type)   :: prj_init
    type(ply_prj_header_type) :: header
    real(kind=rk),allocatable     :: nodal_data(:,:)
    real(kind=rk), allocatable    :: modal_data(:,:)
    real(kind=rk), allocatable    :: ref_modes(:,:)
    !-----------for init
    integer ::basisType, maxdegree, i
    !-----------for oversamp
    real(kind=rk), allocatable    :: oversamp_modal(:,:)
    integer :: iDegX, iDegY, iDegZ,dof, dofOverSamp

    basisType = Q_space
    maxdegree = power
    header%kind = 'fxt'
    header%fxt_header%nodes_header%nodes_kind = 'gauss-legendre'
    header%fxt_header%factor = 1.0_rk

    ! define my poly projection type
    call ply_prj_init_define(me= prj_init,             &
      &                          header = header ,         &
      &                          maxPolyDegree= maxdegree, &
      &                          basisType = basistype ) 
    ! fill the projection body according to the header
    call ply_poly_project_fillbody(me = me, proj_init=prj_init, scheme_dim=3)
  
    allocate(ref_modes(1:me%body_3d%ndofs,1)) 
    allocate(modal_data(1:me%body_3d%ndofs,1) ) 
    allocate(oversamp_modal(1:me%body_3d%oversamp_dofs,1) ) 
    allocate(nodal_data(1:me%body_3d%nQuadPoints,1) ) 

    do i=1, me%body_3d%ndofs
       ref_modes(i,1)=1.0/real(i, kind=rk)
    end do
  
    modal_data = ref_modes
    write(*,*) "modal_data = "
    write(*,*) modal_data
    write(*,*) "------------------------------------------------------"
    call ply_convert2oversample(state       = modal_data,    &
      &                         ndim        = 3,             &
      &                         poly_proj   = me,            &
      &                         modalCoeffs = oversamp_modal )

    modal_data = 0.0
    write(*,*) " modal_data (in oversampled space) = "
    write(*,*) oversamp_modal
    write(*,*) "------------------------------------------------------"
    call ply_poly_project_m2n(me = me ,                &
      &                       dim = 3 ,                &
      &                       nVars = 1,               &
      &                       nodal_data=nodal_data,   &
      &                       modal_data=oversamp_modal)
    write(*,*) "converted to nodal_data (in oversampled space) = "
    write(*,*) nodal_data
    write(*,*) "------------------------------------------------------"
    call ply_poly_project_n2m(me = me,                  & 
      &                       dim = 3 ,                 &
      &                       nVars = 1,                &
      &                       nodal_data=nodal_data,    &
      &                       modal_data=oversamp_modal )

    write(*,*) "converting it  back to modal data (in oversampled space) = "
    write(*,*) oversamp_modal
    write(*,*) "------------------------------------------------------"

    call ply_convertFromOversample(state       = modal_data,    &
      &                            ndim        = 3,             &
      &                            poly_proj   = me,            &
      &                            modalCoeffs = oversamp_modal )

    write(*,*) "Final  modal_data = "
    write(*,*) modal_data
    write(*,*) "------------------------------------------------------"
    
    res= maxval( abs(ref_modes - modal_data) ) 
    write(*,*) 'power=', power
    write(*,*) 'res=', res
  end subroutine  check_fxt_3D
end program test_fxtd_n2m2n
