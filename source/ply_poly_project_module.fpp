?? include "ply_dof_module.inc"
module ply_poly_project_module
  use env_module,                only: rk, labelLen

  use fftw_wrap,                 only: fftw_available

  use aotus_module,              only: flu_State, aot_get_val
  use aot_table_module,          only: aot_table_open, aot_table_close

  use tem_aux_module,            only: tem_abort
  use tem_logging_module,        only: logUnit
  use tem_tools_module,          only: tem_horizontalSpacer

  use ply_modg_basis_module,     only: evalLegendreTensPoly, scalProdLeg
  use ply_dof_module,            only: Q_space, P_space
  use ply_prj_header_module,     only: ply_prj_header_type,        &
    &                                  assignment(=),              &
    &                                  operator(==), operator(>=), &
    &                                  operator(/=), operator(<),  &
    &                                  operator(<=), operator(>)
  use ply_dynArray_project_module, only: dyn_ProjectionArray_type, &
    &                                    ply_fill_dynProjectArray, &
    &                                    ply_prj_init_type

  use ply_LegFpt_module,           only: ply_legFpt_type, &
    &                                    ply_init_legFpt, &
    &                                    ply_legToPnt,    &
    &                                    ply_PntToLeg,    &
    &                                    assignment(=)
  use ply_legFpt_2D_module,        only: ply_pntToLeg_2D,    &
    &                                    ply_legToPnt_2D
  use ply_legFpt_3D_module,        only: ply_pntToLeg_3D,   &
    &                                    ply_legToPnt_3D
  use ply_l2p_module,              only: ply_l2p_type, &
    &                                    ply_init_l2p,     &
    &                                    assignment(=),    &
    &                                    ply_l2p_trafo_3d, &
    &                                    ply_l2p_trafo_2d, &
    &                                    ply_l2p_trafo_1d

  use ply_nodes_module,        only: init_cheb_nodes, init_cheb_nodes_2d,     &
                                   & init_cheb_nodes_1d, init_gauss_nodes,    &
                                   & init_gauss_nodes_2d, init_gauss_nodes_1d,&
                                   & ply_facenodes_type

  use ply_fxt_module, only: ply_fxt_type, ply_init_fxt,                    &
    &                       ply_fxt_m2n_1D,ply_fxt_m2n_3D, ply_fxt_m2n_2D, &
    &                       ply_fxt_n2m_1D,ply_fxt_n2m_3D, ply_fxt_n2m_2D, &
    &                       ply_fxt_type

  implicit none

  private

  !> Additional data, required for the projection.
  type ply_prj_body_type
    !> The fast polynomial transformation which will be used in case
    !! of nonlinear equations. It is used if fpt is choses as projection
    !! method in the lua file
    type(ply_legFpt_type)  :: fpt
    !> The Legendre Polynomial type for the Fast Orthogonal Function
    !! Transform via fxtpack. It is used if 'fxt' is chosen as projection
    !! method in the lua file
    type(ply_fxt_type) :: fxt
    !> Projection method which cam be used for transfoamtion from modal to
    !! nodal space and vice versa. It is used if 'l2p' is chosen as projection
    !! method in the lua file
    type(ply_l2p_type) :: l2p
    !> Volume quadrature points in the reference element
    real(kind=rk), allocatable :: nodes(:,:)
    !> Facial quadrature nodes (reference element) for all three spatial
    !! direction and left and right face.
    !! These points are necessary to transfer boundary conditions given
    !! in physical space to modal space by projection (l2p or fpt)
    type(ply_faceNodes_type), allocatable :: faces(:,:)
    !> quadrature points including oversampling factor
    integer                    :: nQuadPoints
    !> degree of freedom of the scheme depending on maxPolyDegree
    integer                    :: ndofs
    ! the oversamp_dofs are the degrees of freedom for this 'oversampling
    ! degree and equal to (s*(m+1))**d
    integer                    :: oversamp_dofs
    ! minimal number of dofs, used to enable the flexibilty to
    ! set the oversampling factor < 1
    integer                    :: min_dofs
  end type ply_prj_body_type


  !> Projection definition.
  type ply_poly_project_type
    !> Polynomial basis type.
    !!
    !! 3D Monomials have the form x^i * y^j * z^k
    !! - Q_space: quadratic polynomial space (i,j,k) <= maxPolyDegree
    !! - P_space: polynomial space i+j+k <= maxPolyDegree
    integer                       :: basisType

    !> Kind of projection. Currently available:
    !! - 'l2p', L2-Projection
    !! - 'fpt', Fast Polynomial Transformation. Requires the FFTW.
    !! - 'fxt', Fast Polynomial Transformation. uses FXTPACK
    character(len=labelLen)       :: kind

    !> The maximal polynomial degree per spatial direction.
    integer                       :: maxPolyDegree

    !> Using oversampling, the modal space need to be extended according
    ! to the oversampling factor, thus the oversampling degree is (s*m+1)-1
    integer                       :: oversamp_degree

    ! minimal number of dofs, used to enable the flexibilty to
    ! set the oversampling factor < 1
    integer                       :: min_degree

    !> quadrature points including oversampling factor
    ! per spatial direction
    integer                       :: nQuadPointsPerDir

    !> Logical to indicate whether Chebyshev-Lobatto points or simple
    !! Chebyshev points are used
    logical                       :: lobattoPoints = .false.

    !> projection header consits of general information like which kind
    !! of projection is used
!!    type(ply_prj_header_type) :: header
    
    !> In the body datatype, there is for each dimension the main data
    !! for the projection method stored
    type(ply_prj_body_type)   :: body_1d
    type(ply_prj_body_type)   :: body_2d
    type(ply_prj_body_type)   :: body_3d

  end type ply_poly_project_type


  interface assignment(=)
    module procedure Copy_poly_project
    module procedure Copy_poly_project_body
  end interface

  interface ply_poly_project_m2n
    module procedure ply_poly_project_m2n_multiVar
  end interface

  interface ply_poly_project_n2m
    module procedure ply_poly_project_n2m_multiVar
  end interface

  public :: assignment(=)
  public :: ply_fill_project_list
  public :: ply_poly_project_fillbody
  public :: ply_poly_project_m2n
  public :: ply_poly_project_n2m
  public :: ply_poly_project_type
  public :: ply_faceNodes_type
  public :: get_quadpoints_faces
  public :: get_quadpoints_faces_2d
  public :: get_quadpoints_faces_1d
  public :: ply_prj_body_type


contains


  !**************************************************************************!
  subroutine Copy_poly_project(left,right)
    !------------------------------------------------------------------------!
    !> fpt to copy to
    type(ply_poly_project_type), intent(out) :: left
    !> fpt to copy from
    type(ply_poly_project_type), intent(in) :: right
    !------------------------------------------------------------------------!

    left%body_1d = right%body_1d
    left%body_2d = right%body_2d
    left%body_3d = right%body_3d

    left%kind = right%kind
    left%maxPolyDegree = right%maxPolyDegree
    left%basisType = right%basisType
    left%oversamp_degree = right%oversamp_degree
    left%min_degree = right%min_degree
    left%nquadpointsPerDir = right%nquadpointsPerDir
    left%lobattoPoints = right%lobattoPoints

  end subroutine copy_poly_project
  !**************************************************************************!


  !**************************************************************************!
  subroutine Copy_poly_project_body(left,right)
    !------------------------------------------------------------------------!
    ! fpt to copy to
    type(ply_prj_body_type), intent(out) :: left
    ! fpt to copy from
    type(ply_prj_body_type), intent(in) :: right
    !------------------------------------------------------------------------!

    left%fpt = right%fpt
    left%l2p = right%l2p
    left%fxt = right%fxt
    left%nodes = right%nodes
    left%faces = right%faces
    left%nquadpoints = right%nquadpoints
    left%ndofs = right%ndofs
    left%oversamp_dofs = right%oversamp_dofs

  end subroutine copy_poly_project_body
  !**************************************************************************!


  !***************************************************************************!
  !> Fill ups the bodys accroding to the DA.
  subroutine ply_fill_project_list( minLevel, maxLevel, proj_list, &
    &                               dyn_projectArray, scheme_dim   )
  !---------------------------------------------------------------------------!
    integer, intent(in)           :: minLevel
    integer, intent(in)           :: maxLevel
    type(ply_poly_project_type), intent(inout), allocatable :: proj_list(:)
    type(dyn_ProjectionArray_type), intent(in) :: dyn_projectArray
    integer, intent(in) :: scheme_dim
    !-------------------------------------------------------------------------!
    integer :: ipos
    !-------------------------------------------------------------------------!
    call tem_horizontalSpacer(fUnit=logUnit(2))
    write(logUnit(2),*) 'Loading list of projection methods ... '

    ! allocate the poly_proj_list, the maximum value of all positions will give
    ! the number of elements in the DA
    allocate (proj_list(dyn_projectArray%nVals) )
    write(logUnit(5),*) 'the number of elements in proj_pos is =', &
      &                 dyn_projectArray%nVals
    do ipos=1, dyn_projectArray%nVals
      write(logUnit(5),*) 'for pos=', ipos, 'dyn array is', &
        &                 dyn_projectArray%val(ipos)%basisType,&
        &                 dyn_projectArray%val(ipos)%maxpolydegree, &
        &                 dyn_projectArray%val(ipos)%header%kind

      call ply_poly_project_fillbody(                 &
        &    me         = proj_list(ipos),            &
        &    proj_init  = dyn_projectArray%val(ipos), &
        &    scheme_dim = scheme_dim                  )

      write(logUnit(5),*) ' for position', ipos, &
        &                 'of projection list, projection type is'
      write(logUnit(5),*) '  kind = ', proj_list(ipos)%kind
      write(logUnit(5),*) '  maxPolyDegree =', proj_list(ipos)%maxPolyDegree
      write(logUnit(5),*) '  basisType = ', proj_list(ipos)%basisType
      write(logUnit(5),*) '  oversamp_degree = ', &
        &                 proj_list(ipos)%oversamp_degree
      write(logUnit(5),*) '  min_degree = ', proj_list(ipos)%min_degree
      write(logUnit(5),*) '  nquadpointsPerDir = ', &
        &                 proj_list(ipos)%nquadpointsPerDir
      write(logUnit(5),*) '  lobattoPoints = ', proj_list(ipos)%lobattoPoints

    end do

  end subroutine ply_fill_project_list
  !****************************************************************************!


  !****************************************************************************!
  !> Fill the body of the projection with all required data,
  !! ply_poly_project_define has to be used beforehand to set necessary header
  !! information.
  subroutine ply_poly_project_fillbody(me, proj_init, scheme_dim)
    !--------------------------------------------------------------------------!
    type(ply_poly_project_type), intent(inout)  :: me
    type(ply_prj_init_type), intent (in)    :: proj_init
    integer, intent(in) :: scheme_dim
    !--------------------------------------------------------------------------!
    ! the oversampling order, need to get the number of modes in when
    ! oversampling is used
    integer :: oversampling_order
    ! the number of volume quadrature points per spatial direction
    integer :: numQuadPointsPerDir
    ! projection need to be done for each nScalar
    integer :: nvars = 1
    real(kind=rk) :: log_order, rem_log
    real(kind=rk) :: over_factor
    integer :: lb_log
    !-------------------------------------------------------------------------!
    ! set the kind in the final projection type
    me%kind = trim(proj_init%header%kind)
    me%maxPolyDegree = proj_init%maxPolyDegree
    me%basisType = proj_init%basisType

    select case(trim(proj_init%header%kind))
    case('fpt')
      over_factor = proj_init%header%fpt_header%factor
    case('l2p')
      over_factor = proj_init%header%l2p_header%factor
    case('fxt')
      over_factor = proj_init%header%fxt_header%factor
    end select

    ! Find the oversampling order
    oversampling_order = ceiling( over_factor*(me%maxpolyDegree+1) )
    if (trim(proj_init%header%kind) == 'fpt') then
      ! Ensure an appropriate order in the oversampled polynomial
      ! representation.
      if ( proj_init%header%fpt_header%adapt_factor_pow2               &
        &  .and. (iand(oversampling_order, oversampling_order-1) /= 0) &
        &  ) then

        write(logUnit(1),*) '*** NOTE: oversampling order is increased to' &
          &                 //' next power of 2! ***'
        write(logUnit(2),*) '          original oversampling order would' &
          &                 // ' have been: ', oversampling_order
        ! Oversampling_order is not a power of 2, find the next power of 2
        ! and use that one instead...
        log_order = log(real(oversampling_order, kind=rk))/log(2.0_rk)
        lb_log = max(floor(log_order), 0)
        rem_log = log_order - lb_log
        if (rem_log > epsilon(log_order)*lb_log) then
          oversampling_order = 2**(lb_log+1)
        else
          oversampling_order = 2**lb_log
        end if

      end if
    end if
    if (trim(proj_init%header%kind) == 'fxt') then
      ! Ensure an even order for the oversampled polynomial representation.
      ! There is a bug in FXTPACK resulting faulty conversions for odd
      ! numbers of evaluation points.
      ! To avoid this, we increase the oversampling order by 1, if it is odd.
      oversampling_order = oversampling_order + mod(oversampling_order,2)
    end if

    write(logUnit(1),*) 'Using an oversampled order of: ', oversampling_order
    write(logUnit(1),*) 'maxpolydegree is: ', me%maxPolyDegree
    write(logUnit(2),*) '(Actual oversampling factor: ',                   &
      &                 real(oversampling_order)/real(me%maxpolydegree+1), &
      &                 ')'

    numQuadPointsPerDir = oversampling_order

    me%basisType = proj_init%basisType
    me%maxPolyDegree = proj_init%maxPolyDegree
    me%nquadpointsPerDir = numQuadPointsPerDir
    me%oversamp_degree = oversampling_order-1
    me%min_degree = min(me%maxPolyDegree, me%oversamp_degree)

    ! number of dof depending on q_space or p_space
    if (me%basisType == Q_space) then
       me%body_3d%ndofs = (me%maxPolyDegree+1)**3
       me%body_2d%ndofs = (me%maxPolyDegree+1)**2
       me%body_1d%ndofs = me%maxPolyDegree+1
       me%body_3d%min_dofs = (me%min_degree+1)**3
       me%body_2d%min_dofs = (me%min_degree+1)**2
       me%body_1d%min_dofs = me%min_degree+1
    else !p_space
?? copy :: getDofsPTens(me%maxPolyDegree, me%body_3d%ndofs)
?? copy :: getDofsPTens2d(me%maxPolyDegree, me%body_2d%ndofs)
?? copy :: getDofsPTens1d(me%maxPolyDegree, me%body_1d%ndofs)
?? copy :: getDofsPTens(me%min_degree, me%body_3d%min_dofs)
?? copy :: getDofsPTens2d(me%min_degree, me%body_2d%min_dofs)
?? copy :: getDofsPTens1d(me%min_degree, me%body_1d%min_dofs)
    end if

    me%body_3d%nquadpoints = numQuadPointsPerDir**3
    me%body_2d%nquadpoints = numQuadPointsPerDir**2
    me%body_1d%nquadpoints = numQuadPointsPerDir

    me%body_3d%oversamp_dofs = (oversampling_order)**3
    me%body_2d%oversamp_dofs = (oversampling_order)**2
    me%body_1d%oversamp_dofs = oversampling_order

    select case (trim(proj_init%header%kind))
    case('fpt')
      ! Fill fpt datatype

      ! fpt has option for lobattopoints
      me%lobattopoints = proj_init%header%fpt_header%nodes_header%lobattopoints

      !> Initialize the fpt data type
      call ply_init_legfpt(                                                  &
        &    maxPolyDegree    = me%oversamp_degree,                          &
        &    nIndeps          = 1,                                           &
        &    fpt              = me%body_1d%fpt,                              &
        &    lobattoPoints    = me%lobattoPoints,                            &
        &    blocksize        = proj_init%header%fpt_header%blocksize,       &
        &    approx_terms     = proj_init%header%fpt_header%approx_terms,    &
        &    striplen         = proj_init%header%fpt_header%striplen,        &
        &    subblockingWidth = proj_init%header%fpt_header%subblockingWidth )
      !> Initialization/Create  of the volume quadrature  nodes and the
      !! quadrature points on the face
      call init_cheb_nodes_1d(                              &
        &    me = proj_init%header%fpt_header%nodes_header, &
        &    nodes = me%body_1d%nodes,                      &
        &    faces = me%body_1d%faces,                      &
        &    nQuadPointsPerDir = me%nQuadPointsPerDir       )

      if (scheme_dim >= 2) then
        call ply_init_legfpt(                                                  &
          &    maxPolyDegree    = me%oversamp_degree,                          &
          &    nIndeps          = me%oversamp_degree+1,                        &
          &    fpt              = me%body_2d%fpt,                              &
          &    lobattoPoints    = me%lobattoPoints,                            &
          &    blocksize        = proj_init%header%fpt_header%blocksize,       &
          &    approx_terms     = proj_init%header%fpt_header%approx_terms,    &
          &    striplen         = proj_init%header%fpt_header%striplen,        &
          &    subblockingWidth = proj_init%header%fpt_header%subblockingWidth )
        call init_cheb_nodes_2d(                              &
          &    me = proj_init%header%fpt_header%nodes_header, &
          &    nodes = me%body_2d%nodes,                      &
          &    faces = me%body_2d%faces,                      &
          &    nQuadPointsPerDir = me%nQuadPointsPerDir       )
      end if

      if (scheme_dim >= 3) then
        call ply_init_legfpt(                                                  &
          &    maxPolyDegree    = me%oversamp_degree,                          &
          &    nIndeps          = (me%oversamp_degree+1)**2,                   &
          &    fpt              = me%body_3D%fpt,                              &
          &    lobattoPoints    = me%lobattoPoints,                            &
          &    blocksize        = proj_init%header%fpt_header%blocksize,       &
          &    approx_terms     = proj_init%header%fpt_header%approx_terms,    &
          &    striplen         = proj_init%header%fpt_header%striplen,        &
          &    subblockingWidth = proj_init%header%fpt_header%subblockingWidth )
        call init_cheb_nodes(                                 &
          &    me = proj_init%header%fpt_header%nodes_header, &
          &    nodes = me%body_3d%nodes,                      &
          &    faces = me%body_3d%faces,                      &
          &    nQuadPointsPerDir = me%nQuadPointsPerDir       )
      end if

    case('l2p')
      !> Fill the L2 projection datatype
      !! no lobatto points for gauss nodes implemented
      if (scheme_dim >= 3) then
        call ply_init_l2p(l2p    = me%body_3d%l2p,              &
          &               header = proj_init%header%l2p_header, &
          &               degree = me%oversamp_degree,          &
          &               nDims  = 3,                           &
          &               nodes  = me%body_3d%nodes,            &
          &               faces  = me%body_3d%faces             )
      end if

      if (scheme_dim >= 2) then
        call ply_init_l2p(l2p    = me%body_2d%l2p,              &
          &               header = proj_init%header%l2p_header, &
          &               degree = me%oversamp_degree,          &
          &               nDims  = 2,                           &
          &               nodes  = me%body_2d%nodes,            &
          &               faces  = me%body_2d%faces             )
      end if

      call ply_init_l2p(l2p    = me%body_1d%l2p,              &
        &               header = proj_init%header%l2p_header, &
        &               degree = me%oversamp_degree,          &
        &               nDims  = 1,                           &
        &               nodes  = me%body_1d%nodes,            &
        &               faces  = me%body_1d%faces             )

    case ('fxt')
      !> Fill the fxt Legendre Polynomial datatype
      if (scheme_dim >= 3) then
        call ply_init_fxt(fxt    = me%body_3d%fxt,              &
          &               header = proj_init%header%fxt_header, &
          &               degree = me%oversamp_degree,          &
          &               nDims  = 3,                           &
          &               nodes  = me%body_3d%nodes,            &
          &               faces  = me%body_3d%faces             )
      end if

      if (scheme_dim >= 2) then
        call ply_init_fxt(fxt    = me%body_2d%fxt,              &
          &               header = proj_init%header%fxt_header, &
          &               degree = me%oversamp_degree,          &
          &               nDims  = 2,                           &
          &               nodes  = me%body_2d%nodes,            &
          &               faces  = me%body_2d%faces             )
      end if
        call ply_init_fxt(fxt    = me%body_1d%fxt,              &
          &               header = proj_init%header%fxt_header, &
          &               degree = me%oversamp_degree,          &
          &               nDims  = 1,                           &
          &               nodes  = me%body_1d%nodes,            &
          &               faces  = me%body_1d%faces             )

    case default
      write(logUnit(1),*) 'ERROR in initializing projection:'
      write(logUnit(1),*) 'Unknown projection method <', &
        &                 trim(proj_init%header%kind), &
        &                 '>'
      write(logUnit(1),*) 'Stopping....'
      call tem_abort()

    end select

  end subroutine ply_poly_project_fillbody
  !****************************************************************************!


  !****************** MODAL to NODAL ******************************************!
  !> Convert nDoF modes to nodal values.
  subroutine ply_poly_project_m2n_multiVar(me, dim, nVars, modal_data, &
    &                                      nodal_data)
    !--------------------------------------------------------------------------!
    type(ply_poly_project_type), intent(inout) :: me
    integer, intent(in) :: dim
    !> The number of variables to project. If a variable consists of more than
    !! one component, the number of components has to be passed. If there are
    !! more than one variable, the sum of all components has to be passed (e.g.
    !! 6 when there are two three-dimensional vectors).
    integer, intent(in) :: nVars
    real(kind=rk), intent(inout) :: modal_data(:,:)
    real(kind=rk), intent(inout) :: nodal_data(:,:)
    !--------------------------------------------------------------------------!
    integer :: iVar
    !--------------------------------------------------------------------------!

    select case(trim(me%kind))
    case ('l2p')
      ! for the projection modal to nodal, we do not need to distingusih
      ! between the spaces since the modal values results from the computation
      ! and the projection is on evluation of the nodes at the points with
      ! additional summation
      select case(dim)
      case (1)
        do iVar = 1, nVars
          call ply_l2p_trafo_1D( trafo = me%body_1D%l2p%leg2node, &
            &                    projected = nodal_data(:,iVar),  &
            &                    original  = modal_data(:,iVar)   )
        end do
      case (2)
        do iVar = 1, nVars
          call ply_l2p_trafo_2D( trafo = me%body_2D%l2p%leg2node, &
            &                    projected = nodal_data(:,iVar),  &
            &                    original  = modal_data(:,iVar)   )
        end do
      case (3)
        do iVar = 1, nVars
          call ply_l2p_trafo_3D( trafo = me%body_3D%l2p%leg2node, &
            &                    projected = nodal_data(:,iVar),  &
            &                    original  = modal_data(:,iVar)   )
        end do
      end select

    case ('fpt')

      select case (dim)
      case (3)
        call ply_LegToPnt_3D( fpt       = me%body_3d%fpt, &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data,     &
          &                   nVars     = nVars           )
      case (2)
        call ply_LegToPnt_3D( fpt       = me%body_2d%fpt, &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data,     &
          &                   nVars     = nVars           )
      case (1)
        do iVar = 1,nVars
          call ply_LegToPnt( fpt       = me%body_1d%fpt,     &
            &                pntVal    = nodal_data(:,iVar), &
            &                legCoeffs = modal_data(:,iVar), &
            &                nIndeps   = 1                   )
        end do
      end select

    case ('fxt')
      select case (dim)
      case (3)
        do iVar = 1,nVars
          call ply_fxt_m2n_3D( fxt = me%body_3d%fxt,             &
            &               modal_data = modal_data(:,iVar),     &
            &               nodal_data = nodal_data(:,iVar),     &
            &              oversamp_degree = me%oversamp_degree  )
        end do
      case (2)
        do iVar = 1,nVars
          call ply_fxt_m2n_2D( fxt = me%body_2d%fxt,             &
            &               modal_data = modal_data(:,iVar),     &
            &               nodal_data = nodal_data(:,iVar),     &
            &              oversamp_degree = me%oversamp_degree  )
        end do

      case (1)
        do iVar = 1,nVars
          call ply_fxt_m2n_1D( fxt = me%body_1d%fxt,             &
            &               modal_data = modal_data(:,iVar),     &
            &               nodal_data = nodal_data(:,iVar),     &
            &              oversamp_degree = me%oversamp_degree  )
        end do
      end select
    end select

  end subroutine ply_poly_project_m2n_multivar
  !****************************************************************************!



  !***************** NODAL to MODAL *******************************************!
  !> Convert nodal values to nDoFs modes.
  subroutine ply_poly_project_n2m_multiVar(me, dim, nVars, nodal_data, &
    &                                      modal_data)
    !--------------------------------------------------------------------------!
    type(ply_poly_project_type), intent(inout) :: me
    integer, intent(in) :: dim
    integer, intent(in) :: nVars
    real(kind=rk), intent(inout) :: nodal_data(:,:)
    real(kind=rk), intent(inout) :: modal_data(:,:)
    !--------------------------------------------------------------------------!
    integer :: iVar
    !--------------------------------------------------------------------------!

    select case(trim(me%kind))
    case ('l2p')
      select case (dim)
      case (1)
        do iVar = 1, nVars
          call ply_l2p_trafo_1D( trafo = me%body_1D%l2p%node2leg, &
            &                    projected = modal_data(:,iVar),  &
            &                    original  = nodal_data(:,iVar)   )
        end do
      case (2)
        do iVar = 1, nVars
          call ply_l2p_trafo_2D( trafo = me%body_2D%l2p%node2leg, &
            &                    projected = modal_data(:,iVar),  &
            &                    original  = nodal_data(:,iVar)   )
        end do
      case (3)
        do iVar = 1, nVars
          call ply_l2p_trafo_3D( trafo = me%body_3D%l2p%node2leg, &
            &                    projected = modal_data(:,iVar),  &
            &                    original  = nodal_data(:,iVar)   )
        end do
      end select

    case ('fpt')
      !projection via fpt
      select case (dim)
      case (3)
        call ply_pntToLeg_3D( fpt       = me%body_3d%fpt, &
          &                   nVars     = nVars,          &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data      )
      case (2)
        call ply_pntToLeg_2D( fpt       = me%body_2d%fpt, &
          &                   nVars     = nVars,          &
          &                   pntVal    = nodal_data,     &
          &                   legCoeffs = modal_data      )
      case (1)
        do iVar = 1,nVars
          call ply_pntToLeg( fpt       = me%body_1d%fpt,     &
            &                nIndeps   = 1,                  &
            &                pntVal    = nodal_data(:,iVar), &
            &                legCoeffs = modal_data(:,iVar)  )
        end do
      end select

    case ('fxt')
      select case (dim)
      case (3)
        do iVar = 1, nVars 
          call ply_fxt_n2m_3D(                                  &
            &         fxt              = me%body_3d%fxt,        &
            &         nodal_data       = nodal_data(:,iVar),    &
            &         modal_data       = modal_data(:,iVar),    &
            &         oversamp_degree  = me%oversamp_degree     )
        end do

      case (2)
        do iVar = 1, nVars 
          call ply_fxt_n2m_2D(                                  &
            &         fxt              = me%body_2d%fxt,        &
            &         nodal_data       = nodal_data(:,iVar),    &
            &         modal_data       = modal_data(:,iVar),    &
            &         oversamp_degree  = me%oversamp_degree     )
        end do

      case (1)
        do iVar = 1, nVars 
          call ply_fxt_n2m_1D(                                  &
            &         fxt              = me%body_1d%fxt,        &
            &         nodal_data       = nodal_data(:,iVar),    &
            &         modal_data       = modal_data(:,iVar),    &
            &         oversamp_degree  = me%oversamp_degree     )
        end do
      end select

    case default
       write(logUnit(1),*) 'ERROR in projection nodal to modal'
    end select

  end subroutine ply_poly_project_n2m_multivar
  !***************************************************************************!


  !****************************************************************************!
  !> function to provide the coordinates from the quadrature points on the faces
  ! idir and ialign are inputs which identify which face is needed
  ! faces is allocated as face(dir,align)
  subroutine get_quadpoints_faces(poly_proj, idir, ialign, points)
    !--------------------------------------------------------------------------!
    type(ply_poly_project_type), intent(in) :: poly_proj
    integer, intent(in)                     :: idir
    integer, intent(in)                     :: ialign
    real(kind=rk), allocatable, intent(inout) :: points (:,:)
    !--------------------------------------------------------------------------!

     allocate (points(poly_proj%body_3d%faces(idir,iAlign)%nquadpoints,3))
     points = poly_proj%body_3d%faces(idir,iAlign)%points

  end subroutine get_quadpoints_faces
  !****************************************************************************!


  !****************************************************************************!
  subroutine get_quadpoints_faces_2d(poly_proj, idir, ialign, points)
    !--------------------------------------------------------------------------!
    type(ply_poly_project_type), intent(in) :: poly_proj
    integer, intent(in)                     :: idir
    integer, intent(in)                     :: ialign
    real(kind=rk), allocatable, intent(out) :: points (:,:)
    !--------------------------------------------------------------------------!

     allocate (points(poly_proj%body_2d%faces(idir,iAlign)%nquadpoints,3))
     points = poly_proj%body_2d%faces(idir,iAlign)%points

  end subroutine get_quadpoints_faces_2d
  !****************************************************************************!


  !****************************************************************************!
  subroutine get_quadpoints_faces_1d(poly_proj, idir, ialign, points)
    !--------------------------------------------------------------------------!
    type(ply_poly_project_type), intent(in) :: poly_proj
    integer, intent(in)                     :: idir
    integer, intent(in)                     :: ialign
    real(kind=rk), allocatable, intent(out) :: points (:,:)
    !--------------------------------------------------------------------------!
     allocate (points(poly_proj%body_1d%faces(idir,iAlign)%nquadpoints,3))
     points = poly_proj%body_1d%faces(idir,iAlign)%points
  end subroutine get_quadpoints_faces_1d
  !****************************************************************************!

end module ply_poly_project_module
