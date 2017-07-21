module ply_nodes_module
  use env_module,        only: rk, labelLen
  use fftw_wrap,         only: fftw_available
  use aotus_module,      only: flu_State, aot_get_val

  use tem_aux_module,    only: tem_abort

  use ply_space_integration_module, only: ply_create_surface_gauss_points_cube,   & 
   &                                      ply_create_surface_gauss_points_cube_2d,& 
   &                                      ply_create_surface_gauss_points_cube_1d,& 
   &                                      ply_create_volume_gauss_points_cube,    & 
   &                                      ply_create_volume_gauss_points_cube_2d, & 
   &                                      ply_create_volume_gauss_points_cube_1d 
  use ply_equidistant_module,       only: ply_equadPoints_type,             &
   &                                      create_surface_equidistant_points
  use ply_chebPoint_module,         only: create_volume_cheb_points_cube,           &
   &                                      create_volume_cheb_points_cube_2d,        &
   &                                      create_volume_cheb_points_cube_1d,        &
   &                                      create_volume_lobattocheb_points_cube,    &
   &                                      create_volume_lobattocheb_points_cube_2d, &
   &                                      create_volume_lobattocheb_points_cube_1d, &
   &                                      create_surface_cheb_points_cube,          &
   &                                      create_surface_lobattocheb_points_cube,   &
   &                                      create_surface_cheb_points_cube_2d,       &
   &                                      create_surface_lobattocheb_points_cube_2d,&
   &                                      create_surface_cheb_points_cube_1d,       &
   &                                      create_surface_lobattocheb_points_cube_1d
  use ply_nodes_header_module,      only: ply_nodes_header_type, assignment(=)

  implicit none
  private

  !> Datatype to represent facewise nodes
  type ply_faceNodes_type
    !> The number of face nodes
    integer :: nquadPoints
    !> The face nodes.
    !! First index goes from 1 to nPoints and second index
    !! from 1 to 3 for the 3 spatial coordinates.
    real(kind=rk), allocatable :: points(:,:)
    type (ply_equadPoints_type) :: eQuads
  end type ply_faceNodes_type

  public :: init_gauss_nodes, init_gauss_nodes_2d, init_gauss_nodes_1d
  public :: init_cheb_nodes, init_cheb_nodes_2d, init_cheb_nodes_1d
  public :: init_equi_nodes
  public :: ply_faceNodes_type

  contains

  !****************************************************************************!
  !> routine to intilize quadrature points, 3d
  subroutine init_gauss_nodes(nodes, faces, weights, nQuadPointsPerDir)
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable, intent (inout)  :: nodes(:,:)
    type(ply_faceNodes_type), allocatable, intent (inout)  :: faces(:,:)
    integer, intent (in)           :: nQuadPointsPerDir
    real(kind=rk), allocatable,intent (inout)   :: weights(:)
    integer :: iDir, iAlign
    real(kind=rk), allocatable :: tmp_weights(:)
    !--------------------------------------------------------------------------!

    call ply_create_volume_gauss_points_cube(       &
      & num_intp_per_direction = nQuadPointsPerDir, &
      & points                 = nodes,             &
      & weights                = weights,           &
      & refElemMin             = -1.0_rk,           &
      & refElemMax             =  1.0_rk            )

    ! Build the Gauss-Legendre nodes on the reference faces
    allocate( faces(3,2) )
    do iDir = 1,3
      do iAlign = 1,2
        faces(iDir,iAlign)%nquadpoints = nQuadPointsPerDir**2
        call ply_create_surface_gauss_points_cube(              &
          & num_intp_per_direction = nQuadPointsPerDir,         &
          & points                 = faces(iDir,iAlign)%points, &
          & weights                = tmp_weights,               &
          & refElemMin             = -1.0_rk,                   &
          & refElemMax             = 1.0_rk,                    &
          & dir                    = idir,                      &
          & align                  = iAlign                     )
      deallocate(tmp_weights)
    end do
  end do
end subroutine init_gauss_nodes
!****************************************************************************!


!****************************************************************************!
subroutine init_cheb_nodes(me, nodes, faces, nQuadPointsPerDir )
  !--------------------------------------------------------------------------!
  type(ply_nodes_header_type), intent(in)     :: me
  real(kind=rk), allocatable, intent (inout)  :: nodes(:,:)
  type(ply_faceNodes_type), allocatable, intent (inout)  :: faces(:,:)
  integer, intent (in)           :: nQuadPointsPerDir
  integer :: iDir, iAlign
  !--------------------------------------------------------------------------!

  ! Build the Chebyshev nodes on the reference element (volume)
  if (me%lobattoPoints) then
    call create_volume_lobattocheb_points_cube(     &
      & num_intp_per_direction = nQuadPointsPerDir, &
      & points                 = nodes              )
  else
    call create_volume_cheb_points_cube(            &
      & num_intp_per_direction = nQuadPointsPerDir, &
      & points                 = nodes              )
  end if

  ! Build the Chebyshev nodes on the reference faces
  allocate( faces(3,2) )
  do iDir = 1,3
    do iAlign = 1,2

      faces(iDir,iAlign)%nquadpoints = nQuadPointsPerDir**2
      ! Check which points to generate, either Chebyshev or Lobatto-Chebyshev
      ! points
      if(me%lobattoPoints) then
        call create_surface_lobattocheb_points_cube(             &
          & points                 = faces(iDir, iAlign)%points, &
          & num_intp_per_direction = nQuadPointsPerDir,          &
          & dir                    = iDir,                       &
          & align                  = iAlign                      )
      else
        call create_surface_cheb_points_cube(                    &
          & points                 = faces(iDir, iAlign)%points, &
          & num_intp_per_direction = nQuadPointsPerDir,          &
          & dir                    = iDir,                       &
          & align                  = iAlign                      )
      end if

    end do
  end do
end subroutine init_cheb_nodes
  !****************************************************************************!

 !****************************************************************************!
 subroutine init_equi_nodes(nPoly,faces, nDir,nQuadPointsPerDir)
   !--------------------------------------------------------------------------
   !> polynomial degree
   integer, intent(in) :: nPoly
   type(ply_faceNodes_type), intent (inout) :: faces(:,:)
   integer, intent (in) :: nQuadPointsPerDir
   !> Spatial dimension
   integer, intent(in) :: nDir
   !--------------------------------------------------------------------------!
   integer :: idir, iAlign
   !--------------------------------------------------------------------------!
   ! Build the equidistant points on the reference faces
    do idir = 1,nDir
      do iAlign = 1,2
        call create_surface_equidistant_points( &
          & me     = faces(idir,iAlign)%eQuads, & 
          & nDir   = nDir,                      &     
          & iAlign = iAlign,                    &
          & nPoly  = nPoly,                     &
          & iDir   = idir                       )
      end do
    end do
  end subroutine init_equi_nodes
  !****************************************************************************!


  !****************************************************************************!
  subroutine init_gauss_nodes_2d(nodes, faces, weights, nQuadPointsPerDir)
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable, intent (inout)  :: nodes(:,:)
    type(ply_faceNodes_type), allocatable, intent (inout)  :: faces(:,:)
    integer, intent (in)           :: nQuadPointsPerDir
    real(kind=rk), allocatable,intent (inout)   :: weights(:)
    integer :: iDir, iAlign
    real(kind=rk), allocatable :: tmp_weights(:)
    ! -------------------------------------------------------------------- !

    ! Build Gauss-Legendre Points in the volume, 2d
    call ply_create_volume_gauss_points_cube_2d(    &
      & num_intp_per_direction = nQuadPointsPerDir, &
      & points                 = nodes,             &
      & weights                = weights,           &
      & refElemMin             = -1.0_rk,           &
      & refElemMax             = 1.0_rk             )

    ! Build the Gauss-Legendre nodes on the reference faces
    allocate( faces(2,2) )
    do iDir = 1,2
      do iAlign = 1,2

        faces(iDir,iAlign)%nquadpoints = nQuadPointsPerDir
        call ply_create_surface_gauss_points_cube_2d(           &
          & num_intp_per_direction = nQuadPointsPerDir,         &
          & points                 = faces(iDir,iAlign)%points, &
          & weights                = tmp_weights,               &
          & refElemMin             = -1.0_rk,                   &
          & refElemMax             = 1.0_rk,                    &
          & Dir                    = iDir,                      &
          & Align                  = iAlign                     )
        deallocate(tmp_weights)
      end do
    end do

  end subroutine init_gauss_nodes_2d
  !****************************************************************************!


  !****************************************************************************!
  subroutine init_cheb_nodes_2d(me, nodes, faces, nQuadPointsPerDir)
    !--------------------------------------------------------------------------!
    type(ply_nodes_header_type), intent(in)     :: me
    real(kind=rk), allocatable, intent (inout)  :: nodes(:,:)
    type(ply_faceNodes_type), allocatable, intent (inout)  :: faces(:,:)
    integer, intent (in)           :: nQuadPointsPerDir
    integer :: iDir, iAlign
    !--------------------------------------------------------------------------!

    ! Build the Chebyshev nodes on the reference element (volume)
    if (me%lobattoPoints) then
      call create_volume_lobattocheb_points_cube_2d(  &
        & num_intp_per_direction = nQuadPointsPerDir, &
        & points                 = nodes              )
    else
      call create_volume_cheb_points_cube_2d(         &
        & num_intp_per_direction = nQuadPointsPerDir, &
        & points                 = nodes              )
    end if

    ! Build the Chebyshev nodes on the reference faces
    allocate( faces(2,2) )
    do iDir = 1,2
      do iAlign = 1,2

        faces(iDir,iAlign)%nquadpoints = nQuadPointsPerDir
        ! Check which points to generate, either Chebyshev or
        ! Lobatto-Chebyshev points
        if(me%lobattoPoints) then
          call create_surface_lobattocheb_points_cube_2d(          &
            & points                 = faces(iDir, iAlign)%points, &
            & num_intp_per_direction = nQuadPointsPerDir,          &
            & dir                    = iDir,                       &
            & align                  = iAlign                      )
        else
          call create_surface_cheb_points_cube_2d(                 &
            & points                 = faces(iDir, iAlign)%points, &
            & num_intp_per_direction = nQuadPointsPerDir,          &
            & dir                    = iDir,                       &
            & align                  = iAlign                      )
        end if

      end do
    end do
  end subroutine init_cheb_nodes_2d
  !****************************************************************************!


  !****************************************************************************!
  subroutine init_gauss_nodes_1d(nodes, faces, weights, nQuadPointsPerDir)
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable, intent (inout)  :: nodes(:,:)
    real(kind=rk), allocatable, intent (inout)  :: weights(:)
    type(ply_faceNodes_type), allocatable, intent (inout)  :: faces(:,:)
    integer, intent (in)           :: nQuadPointsPerDir
    integer :: iDir, iAlign
    real(kind=rk), allocatable :: tmp_weights(:)
    !--------------------------------------------------------------------------!

    ! Build the Gauss nodes on the reference element (volume)
    call ply_create_volume_gauss_points_cube_1d(    &
      & num_intp_per_direction = nQuadPointsPerDir, &
      & points                 = nodes,             &
      & weights                = weights,           &
      & refElemMin             = -1.0_rk,           &
      & refElemMax             =  1.0_rk            )

    ! Build the Gauss-Legendre nodes on the reference faces
    allocate( faces(1,2) )
    do iDir = 1,1
      do iAlign = 1,2

        faces(iDir,iAlign)%nquadpoints = 1
        call ply_create_surface_gauss_points_cube_1d( &
          & points  = faces(iDir,iAlign)%points,      &
          & weights = tmp_weights,                    &
          & Dir     = iDir,                           &
          & Align   = iAlign                          )
      deallocate(tmp_weights)
      end do
    end do

  end subroutine init_gauss_nodes_1d
  !****************************************************************************!


  !****************************************************************************!
  subroutine init_cheb_nodes_1d(me, nodes, faces, nQuadPointsPerDir)
    !--------------------------------------------------------------------------!
    type(ply_nodes_header_type), intent(in)     :: me
    real(kind=rk), allocatable, intent (inout)  :: nodes(:,:)
    type(ply_faceNodes_type), allocatable, intent (inout)  :: faces(:,:)
    integer, intent (in)           :: nQuadPointsPerDir
    integer :: iDir, iAlign
    !--------------------------------------------------------------------------!

    ! Build the Chebyshev nodes on the reference element (volume)
    if (me%lobattoPoints) then
      call create_volume_lobattocheb_points_cube_1d(  &
        & num_intp_per_direction = nQuadPointsPerDir, &
        & points                 = nodes              )
    else
      call create_volume_cheb_points_cube_1d(         &
        & num_intp_per_direction = nQuadPointsPerDir, &
        & points                 = nodes              )
    end if

    ! Build the Chebyshev nodes on the reference faces
    allocate( faces(1,2) )
    do iDir = 1,1
      do iAlign = 1,2

        faces(iDir,iAlign)%nquadpoints = 1

        ! Check which points to generate, either Chebyshev or Lobatto-Chebyshev points
        if(me%lobattoPoints) then
          call create_surface_lobattocheb_points_cube_1d( &
            & points = faces(iDir, iAlign)%points,        &
            & dir    = iDir,                              &
            & align  = iAlign                             )
        else
          call create_surface_cheb_points_cube_1d( &
            & points = faces(iDir, iAlign)%points, &
            & dir    = iDir,                       &
            & align  = iAlign                      )
        end if

      end do
    end do
  end subroutine init_cheb_nodes_1d
  !****************************************************************************!

end module ply_nodes_module
