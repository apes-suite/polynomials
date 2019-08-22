! Copyright (c) 2012-2014, 2016, 2018 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2016 Langhammer Kay <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop and Harald Klimach
! for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Verena Krupp, Harald Klimach, Peter Vitt
! and Kay Langhammer for University of Siegen.
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !

!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_module
  use, intrinsic :: iso_c_binding
  use env_module,             only: rk
  use tem_compileconf_module, only: vlen
  use ply_polyBaseExc_module, only: ply_trafo_params_type, &
    &                               ply_fpt_init,          &
    &                               ply_fpt_exec,          &
    &                               ply_fpt_single,        &
    &                               ply_legToCheb_param,   &
    &                               ply_chebToLeg_param,   &
    &                               assignment(=)
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

    !> Flag whether to use Lobatto points (include boundary points)
    logical :: use_lobatto_points
  end type ply_legFpt_type

  interface assignment(=)
    module procedure Copy_fpt
  end interface

  public :: ply_legFpt_type, ply_init_legFpt, ply_legToPnt, ply_pntToLeg
  public :: assignment(=)


contains


  ! ************************************************************************ !
  subroutine Copy_fpt( left, right )
    ! -------------------------------------------------------------------- !
    !> fpt to copy to
    type(ply_legFpt_type), intent(out) :: left
    !> fpt to copy from
    type(ply_legFpt_type), intent(in) :: right
    ! -------------------------------------------------------------------- !

    left%legToChebParams = right%legToChebParams
    left%chebToLegParams = right%chebToLegParams

    left%planChebToPnt = right%planChebToPnt
    left%planPntToCheb = right%planPntToCheb

  end subroutine Copy_fpt
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to initialize the fast polynomial transformation
  !! for Legendre expansion.
  subroutine ply_init_legFpt( maxPolyDegree, nIndeps, fpt, blocksize, &
    &                         approx_terms, striplen, lobattoPoints,  &
    &                         subblockingWidth, fft_flags             )
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
    logical, intent(in), optional :: lobattoPoints

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

    fpt%use_lobatto_points = lob

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
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt( fpt, legCoeffs, pntVal, nIndeps )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(inout) :: legCoeffs(:)
    type(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n)
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    if (.not. fpt%use_lobatto_points) then

      do iDof = 1, nIndeps*n, n
        call ply_fpt_single( alph   = legCoeffs(iDof:iDof+n-1), &
          &                  gam    = cheb,                     &
          &                  params = fpt%legToChebParams       )

        ! Normalize the coefficients of the Chebyshev polynomials due
        ! to the unnormalized version of DCT in the FFTW.
        cheb(2:n:2) = -0.5_rk * cheb(2:n:2)
        cheb(3:n:2) =  0.5_rk * cheb(3:n:2)

        call fftw_execute_r2r( fpt%planChebToPnt,    &
          &                    cheb,                 &
          &                    pntVal(iDof:iDof+n-1) )
      end do

    else

      do iDof = 1, nIndeps*n, n
        call ply_fpt_single( alph   = legCoeffs(iDof:iDof+n-1), &
          &                  gam    = cheb,                     &
          &                  params = fpt%legToChebParams       )

        ! Normalize the coefficients of the Chebyshev polynomials due
        ! to the unnormalized version of DCT in the FFTW.
        cheb(2:n-1) = 0.5_rk * cheb(2:n-1)

        call fftw_execute_r2r( fpt%planChebToPnt,    &
          &                    cheb,                 &
          &                    pntVal(iDof:iDof+n-1) )
      end do

    end if ! lobattoPoints

  end subroutine ply_legToPnt
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg( fpt, pntVal, legCoeffs, nIndeps )
    ! -------------------------------------------------------------------- !
    type(ply_legFpt_type), intent(inout) :: fpt
    real(kind=rk), intent(inout) :: pntVal(:)
    real(kind=rk), intent(inout) :: legCoeffs(:)
    integer, intent(in) :: nIndeps
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: cheb(fpt%legToChebParams%n)
    real(kind=rk) :: normFactor
    integer :: iDof
    integer :: n
    ! -------------------------------------------------------------------- !

    n = fpt%legToChebParams%n

    if (.not. fpt%use_lobatto_Points) then

      normFactor = 1.0_rk / real(n,kind=rk)
      do iDof = 1, nIndeps*n, n
        call fftw_execute_r2r( fpt%planPntToCheb,     &
          &                    pntVal(iDof:iDof+n-1), &
          &                    cheb                   )
        ! Normalize the coefficients of the Chebyshev polynomials due
        ! to the unnormalized version of DCT in the FFTW.
        cheb(1) = cheb(1) * 0.5_rk * normfactor
        cheb(2:n:2) = -normFactor * cheb(2:n:2)
        cheb(3:n:2) =  normFactor * cheb(3:n:2)

        call ply_fpt_single( gam    = legCoeffs(iDof:iDof+n-1), &
          &                  alph   = cheb,                     &
          &                  params = fpt%ChebToLegParams       )
      end do

    else

      normFactor = 0.5_rk / real(n-1,kind=rk)
      do iDof = 1, nIndeps*n, n
        call fftw_execute_r2r( fpt%planPntToCheb,     &
          &                    pntVal(iDof:iDof+n-1), &
          &                    cheb                   )
        ! Normalize the coefficients of the Chebyshev polynomials due
        ! to the unnormalized version of DCT in the FFTW.
        cheb(1) = cheb(1) * normFactor
        cheb(2:n-1) = 2.0_rk * normFactor * cheb(2:n-1)
        cheb(n) = cheb(n) * normFactor

        call ply_fpt_single( gam    = legCoeffs(iDof:iDof+n-1), &
          &                  alph   = cheb,                     &
          &                  params = fpt%ChebToLegParams       )
      end do

    end if ! lobattoPoints

  end subroutine ply_pntToLeg
  ! ************************************************************************ !

end module ply_legFpt_module
