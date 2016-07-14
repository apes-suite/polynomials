!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_module
  use, intrinsic :: iso_c_binding
  use env_module, only: rk
  use tem_param_module, only: PI
  use tem_compileconf_module, only: vlen
  use tem_aux_module, only: tem_abort
  use ply_polyBaseExc_module, only: ply_trafo_params_type, &
    &                               ply_fpt_init,          &
    &                               ply_fpt_exec_striped,  &
    &                               ply_fpt_exec,          &
    &                               ply_legToCheb_param,   &
    &                               ply_chebToLeg_param,   &
    &                               assignment(=)
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


  ! ------------------------------------------------------------------------ !
  !> Subroutine to initialize the fast polynomial transformation
  !! for Legendre expansion.
  subroutine ply_init_legFpt( maxPolyDegree, nIndeps, fpt,               &
    &                         blocksize, approx_terms, striplen,         &
    &                         lobattoPoints, subblockingWidth, fft_flags )
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: maxPolyDegree

    !> Number of independent values that can be computed simultaneously.
    integer, intent(in) :: nIndeps
    type(ply_legFpt_type), intent(inout) :: fpt

    !> Smallest block that is approximated by approx_terms coefficients.
    !!
    !! Please note, that this has to be larger than 2*approx_terms to result
    !! in a reduced number of operations. Default is 64.
    integer, optional, intent(in) :: blocksize

    !> Number of approximation terms used to compute off-diagonal products.
    !!
    !! Defaults to 18, which is the suggested accuracy for double precision.
    integer, optional, intent(in) :: approx_terms

    !> Length to use in vectorization, this is the number of independent
    !! matrix multiplications that are to be done simultaneously.
    !!
    !! Defaults to vlen, which may be set at compile time.
    integer, optional, intent(in) :: striplen

    !> Use Chebyshev-Lobatto Points (true) or simple Chebyshev points (false)
    !!
    !! Default is false.
    logical, intent(in), optional  :: lobattoPoints

    !> The width of the subblocks used during the unrolled base exchange to
    !! ensure a better cache usage.
    integer, optional, intent(in) :: subblockingWidth

    !> Planning flags for the FFT.
    !!
    !! Configuration to how much time to spend on finding an optimal FFT
    !! implementation in the FFTW.
    !! See: http://www.fftw.org/doc/Planner-Flags.html#Planner-Flags
    !!
    !! Defaults to FFTW_MEASURE.
    integer, optional, intent(in) :: fft_flags
    ! -------------------------------------------------------------------- !
    real(kind=rk), allocatable :: tmpOut(:), tmpIn(:)
    logical :: lob
    integer :: n
    integer :: maxstriplen
    integer :: planning_flags
    ! -------------------------------------------------------------------- !

    if (present(striplen)) then
      maxstriplen = min(striplen, nIndeps)
    else
      maxstriplen = min(vlen, nIndeps)
    end if

    if (present(fft_flags)) then
      planning_flags = fft_flags
    else
      planning_flags = FFTW_MEASURE
    end if

    lob = .false.
    if (present(lobattoPoints)) then
      lob = lobattoPoints
    end if

    ! Init the fast Legendre to Chebyshev transformation.
    call ply_fpt_init( n                = maxPolyDegree+1,     &
      &                params           = fpt%legToChebParams, &
      &                trafo            = ply_legToCheb_param, &
      &                blocksize        = blocksize,           &
      &                approx_terms     = approx_terms,        &
      &                striplen         = maxstriplen,         &
      &                subblockingWidth = subblockingWidth     )

    ! Init the fast Chebyshev to Legendre transformation.
    call ply_fpt_init( n                = maxPolyDegree+1,     &
      &                params           = fpt%chebToLegParams, &
      &                trafo            = ply_chebToLeg_param, &
      &                blocksize        = blocksize,           &
      &                approx_terms     = approx_terms,        &
      &                striplen         = maxstriplen,         &
      &                subblockingWidth = subblockingWidth     )

    ! Create the buffers for the intermediate arrays
    n = fpt%legToChebParams%n

    ! Temporary arrays to initialize FFTW real->real transformations
    allocate( tmpIn(n) )
    allocate( tmpOut(n) )

    if (.not.lob) then
      ! Init the DCT III ( Leg -> Point values )
      fpt%planChebToPnt = fftw_plan_r2r_1d( n     = n,             &
        &                                   in    = tmpIn,         &
        &                                   out   = tmpOut,        &
        &                                   kind  = FFTW_REDFT01,  &
        &                                   flags = planning_flags )

      ! Init the DCT II ( Point values -> Leg )
      fpt%planPntToCheb = fftw_plan_r2r_1d( n     = n,             &
        &                                   in    = tmpIn,         &
        &                                   out   = tmpOut,        &
        &                                   kind  = FFTW_REDFT10,  &
        &                                   flags = planning_flags )

    else

      ! Init the DCT I  (Leg -> nodal):
      !   To be used with a normalization factor for trafo ...
      fpt%planChebToPnt = fftw_plan_r2r_1d( n     = n,             &
        &                                   in    = tmpIn,         &
        &                                   out   = tmpOut,        &
        &                                   kind  = FFTW_REDFT00,  &
        &                                   flags = planning_flags )

      ! Init the DCT I  (nodal -> Leg):
      !   To be used with a normalization factor for trafo ...
      fpt%planPntToCheb = fpt%planChebToPnt

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

