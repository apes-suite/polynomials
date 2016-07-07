!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_module
  use, intrinsic :: iso_c_binding
  use env_module, only: rk
  use tem_param_module, only: PI
  use tem_aux_module, only: tem_abort
  use ply_polyBaseExc_module, only: ply_trafo_params_type, &
                                  & ply_fpt_init, &
                                  & ply_fpt_exec_striped, &
                                  & ply_fpt_exec, &
                                  & ply_legToCheb_param, ply_chebToLeg_param,&
                                  & assignment(=)
  use ply_nodes_module,        only: ply_faceNodes_type
  use fftw_wrap

  implicit none

  private

  !> Datatype for parameters of the FPT used for 1d, 2d and 3d.
  !!
  !! Stores of the parameters for a fast conversion of a modal
  !! Legendre expansion to point values (located at Chebyshev nodes)
  !! and vice versa. \n
  !! The FPT covers: \n
  !! - Transformation from Legendre expansion to point values 
  !!   at Chebyshev nodes \n
  !! - Transformation from point values (Chebyshev nodes) to 
  !!   modal Legendre expansion \n
  type ply_legFpt_type
    !> FPT params for the fast base exchange from Legendre to 
    !! Chebyshev expansion.
    type(ply_trafo_params_type) :: legToChebParams
    !> FPT params for the fast base exchange from Chebyshev to 
    !! Legendre expansion.
    type(ply_trafo_params_type) :: chebToLegParams
    !> FFTW plan for DCT from Chebyshev coefficients to point values.
    type(C_PTR) :: planChebToPnt
    !> FFTW plan for DCT from point values to Chebyshev coefficients.
    type(C_PTR) :: planPntToCheb
  end type ply_legFpt_type

  interface assignment(=)
    module procedure Copy_fpt
  end interface

  public :: ply_legFpt_type, ply_init_legFpt, ply_legToPnt, ply_pntToLeg
  public :: assignment(=)


contains


  !****************************************************************************
  subroutine Copy_fpt(left,right)
    !---------------------------------------------------------------------------
    !> fpt to copy to
    type(ply_legFpt_type), intent(out) :: left
    !> fpt to copy from
    type(ply_legFpt_type), intent(in) :: right
    !---------------------------------------------------------------------------

    left%legToChebParams = right%legToChebParams
    left%chebToLegParams = right%chebToLegParams 

    left%planChebToPnt = right%planChebToPnt
    left%planPntToCheb = right%planPntToCheb

  end subroutine Copy_fpt
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to initialize the fast polynomial transformation
  !! for Legendre expansion.
  !! VK: change maxPolyDegree to number of quadpoints to include the dealising factor 
  subroutine ply_init_legFpt(maxPolyDegree, fpt, blocksize, lobattoPoints, &
    &                        subblockingWidth                              )
    !---------------------------------------------------------------------------
    integer, intent(in) :: maxPolyDegree
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The blocksize of the underlying fast Legendre to Chebyshev 
    !! transformation.
    integer, intent(in), optional :: blocksize
    !> Use Chebyshev-Lobtatto Points (true) or simple Chebyshev points (false)
    logical, intent(in), optional :: lobattoPoints
    !> The width of the subblocks used during the unrolled base exchange to 
    !! ensure a better cache usage.
    integer, optional, intent(in) :: subblockingWidth
    !---------------------------------------------------------------------------
   real(kind=rk), allocatable :: tmpOut(:,:), tmpIn(:,:)
   logical :: lob
   integer :: n
   integer :: nVars=1
   !---------------------------------------------------------------------------
  
   lob = .false.
   if(present(lobattoPoints)) then
     if (lobattoPoints) then
       lob = .true.
     end if
   end if

   ! init the fast Legendre to Chebyshev transformation.
   call ply_fpt_init(n                = maxPolyDegree+1,     &
     &               params           = fpt%legToChebParams, &
     &               trafo            = ply_legToCheb_param, &
     &               blocksize        = blocksize,           &
     &               subblockingWidth = subblockingWidth     )
   call ply_fpt_init(n                = maxPolyDegree+1,     &
     &               params           = fpt%chebToLegParams, &
     &               trafo            = ply_chebToLeg_param, &
     &               blocksize        = blocksize,           &
     &               subblockingWidth = subblockingWidth     )

 
   ! Create the buffers for the intermediate arrays
   n = fpt%legToChebParams%n**2

   ! Temprorary arrays to initialize FFTW real->real transformations
   allocate( tmpIn(n, nVars) )
   allocate( tmpOut(n, nVars) )
     !allocate( fpt%pntVal( fpt%legToChebParams%n  ))
   
   if(.not.lob) then
    ! Init the DCT III ( Leg -> Point values )
    fpt%planChebToPnt = fftw_plan_r2r_1d(fpt%legToChebParams%n,tmpIn,&
                             & tmpOut,FFTW_REDFT01,FFTW_ESTIMATE )

    ! Init the DCT II ( Point values -> Leg )
    fpt%planPntToCheb = fftw_plan_r2r_1d(fpt%legToChebParams%n,tmpIn,&
                             & tmpOut,FFTW_REDFT10,FFTW_ESTIMATE )
   else
     ! Init the DCT I  (Leg -> nodal): To be used with a normalization factor for trafo ...
     fpt%planChebToPnt = fftw_plan_r2r_1d( n = fpt%legToChebParams%n,&
                                        &  in = tmpIn,&   !legCoeffs
                                        &  out = tmpOut, & !pntVal
                                        &  kind = FFTW_REDFT00, &
                                        &  flags = FFTW_ESTIMATE )

     ! Init the DCT I  (nodal -> Leg): To be used with a normalization factor for trafo ...
     fpt%planPntToCheb = fftw_plan_r2r_1d( n = fpt%legToChebParams%n,&
                                        &  in = tmpIn,&   !pntVal
                                        &  out = tmpOut, & !legCoeffs
                                        &  kind = FFTW_REDFT00, &
                                        &  flags = FFTW_ESTIMATE )
   end if
   deallocate( tmpIn, tmpOut )

  end subroutine ply_init_legFpt
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt( fpt, legCoeffs, pntVal, lobattoPoints )
   !---------------------------------------------------------------------------
   real(kind=rk), intent(inout) :: legCoeffs(:) 
   type(ply_legFpt_type), intent(inout) :: fpt
  ! type(ply_legFpt_type), intent(inout) :: fpt
   real(kind=rk), intent(inout) :: pntVal(:)
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iFunc
   !---------------------------------------------------------------------------

   ! ply_fpt_exec on temp (no memory transpose)
   call ply_fpt_exec( alph = legCoeffs,              &
    &                 gam = pntVal,                  &
    &                 nIndeps = 1,                   &
    &                 plan = fpt%planChebToPnt,      &
    &                 lobattoPoints = lobattoPoints, &
    &                 params = fpt%legToChebParams   )
  

  end subroutine ply_legToPnt
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg( fpt, pntVal, legCoeffs, lobattoPoints )
   !---------------------------------------------------------------------------
   type(ply_legFpt_type), intent(inout) :: fpt
   real(kind=rk), intent(inout) :: pntVal(:)
   real(kind=rk), intent(inout) :: legCoeffs(:) 
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iFunc
   !---------------------------------------------------------------------------

   ! ply_fpt_exec on temp (no memory transpose)
   call ply_fpt_exec( alph = pntVal,                 &
    &                 gam = legCoeffs,               &
    &                 nIndeps = 1,                   &
    &                 plan = fpt%planPntToCheb,      &
    &                 lobattoPoints = lobattoPoints, &
    &                 params = fpt%chebToLegParams   )

  end subroutine ply_pntToLeg

end module ply_legFpt_module

