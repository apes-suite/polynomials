! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014, 2016, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Simon Zimny and Harald Klimach
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

?? include "ply_dof_module.inc"
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_2D_module
  use, intrinsic :: iso_c_binding


  use fftw_wrap

  use env_module,         only: rk
  use ply_legFpt_module,  only: ply_legFpt_type, &
    &                           ply_legToPnt, &
    &                           ply_pntToLeg, &
    &                           assignment(=)

  implicit none

  private

  interface ply_legToPnt_2D
    module procedure ply_legToPnt_2D_singVar
    module procedure ply_legToPnt_2D_multVar
  end interface ply_legToPnt_2D

  interface ply_pntToLeg_2D
    module procedure ply_pntToLeg_2D_singVar
    module procedure ply_pntToLeg_2D_multVar
  end interface ply_pntToLeg_2D

  public :: ply_legToPnt_2D, ply_pntToLeg_2D


contains


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_2D_singVar( fpt, legCoeffs, pntVal )
   ! --------------------------------------------------------------------- !
   !> The FPT parameters.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:)
   !> The resulting point values (Chebyshev nodes).
   real(kind=rk), intent(inout) :: pntVal(:)
   ! --------------------------------------------------------------------- !
   integer :: striplen, iStrip, n, iAlph, nIndeps
   real(kind=rk), dimension(:), allocatable :: alph
   real(kind=rk), dimension(:), allocatable :: gam
   ! --------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iStrip, iAlph, nIndeps)

   striplen = fpt%legToChebParams%striplen
   n = fpt%legToChebParams%n

   allocate(alph(min(striplen, n)*n))
   allocate(gam(min(striplen, n)*n))
   ! original layout (n = 3):
   !  1  2  3
   !  4  5  6
   !  7  8  9

   ! layout after y-trafo:
   !  1  4  7
   !  2  5  8
   !  3  6  9
   !$OMP DO
   yStripLoop: do iStrip = 1, n, striplen
     ! iAlph is the index of the first element in a line for the transformation
     ! in y-direction.
     do iAlph = iStrip, min(iStrip+striplen-1, n) !y-Trafo
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = legCoeffs(iAlph::n) !ytrafo
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n-iStrip+1)

     call ply_legToPnt( fpt       = fpt,     &
       &                nIndeps   = nIndeps, &
       &                legCoeffs = alph,    &
       &                pntVal    = gam      )

     ! Write gam to pntVal array
     pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

   end do yStripLoop
   !$OMP END DO

   ! x-direction
   !$OMP DO
   xStripLoop: do iStrip = 1, n, striplen
     do iAlph = iStrip, min(iStrip+striplen-1, n)
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = pntVal(iAlph::n) !ztrafo
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n-iStrip+1)

     call ply_legToPnt( fpt       = fpt,     &
       &                nIndeps   = nIndeps, &
       &                legCoeffs = alph,    &
       &                pntVal    = gam      )

     pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

   end do xStripLoop
   !$OMP END DO

   !$OMP END PARALLEL

  end subroutine ply_legToPnt_2D_singVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_2D_multVar( fpt, legCoeffs, pntVal, nVars )
   ! --------------------------------------------------------------------- !
   !> The FPT parameters.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:,:)
   !> The resulting point values (Chebyshev nodes).
   real(kind=rk), intent(inout) :: pntVal(:,:)
   !> The number of scalar variables to transform.
   integer, intent(in) :: nVars
   ! --------------------------------------------------------------------- !
   integer :: iVar
   ! --------------------------------------------------------------------- !

   do iVar = 1, nVars
     call ply_legToPnt_2D(fpt, legCoeffs(:,iVar), pntVal(:,iVar))
   end do

  end subroutine ply_legToPnt_2D_multVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_2D_singVar( fpt, pntVal, legCoeffs )
   ! --------------------------------------------------------------------- !
    !> Parameters of the Fast Polynomial transformation.
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The point values to transform to 2D modal Legendre expansion.
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying DCT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: pntVal(:)
    !> The Legendre coefficients.
    real(kind=rk), intent(inout) :: legCoeffs(:)
   ! --------------------------------------------------------------------- !
    integer :: iStrip, striplen, nIndeps, iAlph, n, n_squared
    real(kind=rk), dimension(:), allocatable :: alph
    real(kind=rk), dimension(:), allocatable :: gam
   ! --------------------------------------------------------------------- !
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(nIndeps, iStrip, iAlph)

    striplen = fpt%chebToLegParams%striplen
    n = fpt%legToChebParams%n
    n_squared = n**2

    allocate(alph(min(striplen, n)*n))
    allocate(gam(min(striplen, n)*n))

    !$OMP DO
    yStripLoop: do iStrip = 1, n, striplen
      do iAlph = iStrip, min(iStrip+striplen-1, n) !y_Trafo
        alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = pntVal(iAlph::n)
      end do

      ! At the end of the array the number of computed strips might be smaller
      nIndeps = min(striplen, n-iStrip+1)

      call ply_pntToLeg( fpt       = fpt,     &
        &                nIndeps   = nIndeps, &
        &                legCoeffs = gam,     &
        &                pntVal    = alph     )

      ! temp -> pntVal (stride-1 writing)
      legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)

    end do yStripLoop ! iStrip
    !$OMP END DO

    ! x-direction
    !$OMP DO
    xStripLoop: do iStrip = 1,n,striplen
      do iAlph = iStrip, min(iStrip+striplen-1, n)
        !ztrafo
        alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = legCoeffs(iAlph::n)
      end do

      ! At the end of the array the number of computed strips might be smaller
      nIndeps = min(striplen, n-iStrip+1)

      call ply_pntToLeg( fpt       = fpt,     &
        &                nIndeps   = nIndeps, &
        &                legCoeffs = gam,     &
        &                pntVal    = alph     )

      legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)

    end do xStripLoop
    !$OMP END DO

    !$OMP END PARALLEL

  end subroutine ply_pntToLeg_2D_singVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_2D_multVar( fpt, pntVal, legCoeffs, nVars )
   ! --------------------------------------------------------------------- !
   !> Parameters of the Fast Polynomial transformation.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The point values to transform to 2D modal Legendre expansion.
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying DCT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: pntVal(:,:)
   !> The Legendre coefficients.
   real(kind=rk), intent(inout) :: legCoeffs(:,:)
   !> The number of scalar variables to transform.
   integer, intent(in) :: nVars
   ! --------------------------------------------------------------------- !
   integer :: iVar
   ! --------------------------------------------------------------------- !

   do iVar = 1, nVars
     call ply_pntToLeg_2D(fpt, pntVal(:,iVar), legCoeffs(:,iVar))
   end do

  end subroutine ply_pntToLeg_2D_multVar
  ! ************************************************************************ !

end module ply_legFpt_2D_module
