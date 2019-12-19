! Copyright (c) 2012-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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

?? include "ply_dof_module.inc"
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_3D_module
  use, intrinsic :: iso_c_binding
  use fftw_wrap
  use env_module,         only: rk
  use ply_legFpt_module,  only: ply_legFpt_type, &
    &                           ply_legToPnt, &
    &                           ply_PntToLeg

  implicit none

  private

  interface ply_LegTopnt_3D
    module procedure ply_LegTopnt_3D_multVar
    module procedure ply_LegTopnt_3D_singVar
  end interface ply_LegToPnt_3D

  interface ply_pntToLeg_3D
    module procedure ply_pntToLeg_3D_multVar
    module procedure ply_pntToLeg_3D_singVar
  end interface ply_pntToLeg_3D

  public :: ply_legToPnt_3D, ply_pntToLeg_3D


contains


  ! ************************************************************************ !
  subroutine ply_legToPnt_3D_singvar( fpt, legCoeffs, pntVal )
   ! --------------------------------------------------------------------- !
   !> The FPT parameters.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:)
   real(kind=rk), intent(inout) :: pntVal(:)
   ! --------------------------------------------------------------------- !
   integer :: n, n_squared, n_cubed
   integer :: nIndeps
   ! The normfactor has to be positive for odd sum of exponents and negative
   ! for even sums of exponents.
   ! (original definition: normfactor = (-1)**(i-1) * (-1)**(j-1) * (-1)**(k-1)
   !  can be written as (-1)**(i+j+k-3) = (-1)**(i+j+k-1) for i+j+k>=3
   !  the power can be avoided completely by normfactor = 1 - 2*mod(i+j+k-1,2)
   !  and the real multiplication can be avoided by putting this into a lookup
   !  table.)
   integer :: iStrip
   integer :: striplen
   integer :: iAlph
   real(kind=rk), allocatable :: alph(:)
   real(kind=rk), allocatable :: gam(:)
   ! --------------------------------------------------------------------- !

   striplen = fpt%legToChebParams%striplen
   n = fpt%legToChebParams%n
   n_squared = fpt%legToChebParams%n**2
   n_cubed = n_squared * fpt%legToChebParams%n

   ! number of strips to execute the fpt on in one call
   ! (usually n_squared, but a smaller value might be assigned at the end of
   ! the array)
   nIndeps = n_squared

   allocate(alph(min(striplen,n_squared)*n))
   allocate(gam(min(striplen,n_squared)*n))


  ! z-direction
     ! original layout (n = 3):
     !  1  2  3   10 11 12   19 20 21
     !  4  5  6   13 14 15   22 23 24
     !  7  8  9   16 17 18   25 26 27
     !
     ! example for striplen = 2:
     ! iStrip           1                   3
     ! iAlph            1         2         3         4
     ! index legCoeffs  1 10 19   2 11 20   3 12 21   4 13 22
     ! index alph/gam   1  2  3   4  5  6   1  2  3   4  5  6
     ! index PntVal     1  2  3   4  5  6   7  8  9  10 11 12
     !
     ! layout after z-trafo: (n = 3)
     !  1 10 19    4 13 22    7 16 25
     !  2 11 20    5 14 23    8 17 26
     !  3 12 21    6 15 24    9 18 27
     !
     ! layout after y-trafo:
     !  1  4  7    2  5  8   3   6  9
     ! 10 13 16   11 14 17   12 15 18
     ! 19 22 25   20 23 26   21 24 27
     !
     ! original layout after x trafo:
     !  1  2  3   10 11 12   19 20 21
     !  4  5  6   13 14 15   22 23 24
     !  7  8  9   16 17 18   25 26 27

   ! zStripLoop: Loop over all strips in z-direction
   zStripLoop: do iStrip = 1, n_squared, striplen
     ! iAlph is the index of the first element in a line for the transformation
     ! in z-direction.
     do iAlph = iStrip, min(iStrip+striplen-1, n_squared)
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = legCoeffs(iAlph::n_squared)
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n_squared-iStrip+1)

     call ply_legToPnt( fpt       = fpt,     &
       &                nIndeps   = nIndeps, &
       &                legCoeffs = alph,    &
       &                pntVal    = gam      )

     pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

   end do zStripLoop

  ! y-direction

   yStripLoop: do iStrip = 1,n_squared,striplen
     do iAlph = iStrip, min(iStrip+striplen-1, n_squared)
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = pntVal(iAlph::n_squared)
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n_squared-iStrip+1)

     call ply_legToPnt( fpt       = fpt,     &
       &                nIndeps   = nIndeps, &
       &                legCoeffs = alph,    &
       &                pntVal    = gam      )

     legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

   end do yStripLoop ! iStrip

  ! x-direction
   xStripLoop: do iStrip = 1,n_squared,striplen
     do iAlph = iStrip, min(iStrip+striplen-1, n_squared)
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = legCoeffs(iAlph::n_squared)
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n_squared-iStrip+1)

     call ply_legToPnt( fpt       = fpt,     &
       &                nIndeps   = nIndeps, &
       &                legCoeffs = alph,    &
       &                pntVal    = gam      )

     pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

   end do xStripLoop


  end subroutine ply_legToPnt_3D_singVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  !!VK: no multivar fashion of this routine is used anymore
  subroutine ply_legToPnt_3D_multVar( fpt, legCoeffs, pntVal, nVars )
   ! --------------------------------------------------------------------- !
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:,:)
   type(ply_legFpt_type), intent(inout) :: fpt
   real(kind=rk), intent(inout) :: pntVal(:,:)
   integer, intent(in) :: nVars
   ! --------------------------------------------------------------------- !
   integer :: iVar
   ! --------------------------------------------------------------------- !

   do iVar = 1, nVars
    call ply_legToPnt_3D( fpt, legCoeffs(:,iVar), pntVal(:,iVar) )
   end do

  end subroutine ply_legToPnt_3D_multVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_3D_multVar( fpt, pntVal, legCoeffs, nVars )
   ! --------------------------------------------------------------------- !
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The point values to transform to 3D modal Legendre expansion.
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying DCT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: pntVal(:,:)
   real(kind=rk), intent(inout) :: legCoeffs(:,:)
   integer, intent(in) :: nVars
   ! --------------------------------------------------------------------- !
   integer :: iVar
   ! --------------------------------------------------------------------- !

   do iVar = 1, nVars
     call ply_pntToLeg_3D( fpt, pntVal(:,iVar), legCoeffs(:,iVar) )
   end do

  end subroutine ply_pntToLeg_3D_multVar
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_3D_singVar( fpt, pntVal, legCoeffs )
   ! --------------------------------------------------------------------- !
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The point values to transform to 3D modal Legendre expansion.
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying DCT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: pntVal(:)
    real(kind=rk), intent(inout) :: legCoeffs(:)
   ! --------------------------------------------------------------------- !
    integer :: striplen
    integer :: iStrip
    integer :: iAlph
    integer :: n
    integer :: n_squared
    integer :: n_cubed
    integer :: nIndeps
    real(kind=rk), dimension(:), allocatable :: alph
    real(kind=rk), dimension(:), allocatable :: gam
   ! --------------------------------------------------------------------- !

    striplen = fpt%legToChebParams%striplen
    n = fpt%legToChebParams%n
    n_squared = fpt%legToChebParams%n**2
    n_cubed = n_squared * fpt%legToChebParams%n

    ! number of strips to execute the fpt on in one call
    ! (usually n_squared, but a smaller value might be assigned at the end of
    ! the array)
    nIndeps = n_squared

    allocate(alph(min(striplen,n_squared)*n))
    allocate(gam(min(striplen,n_squared)*n))

    ! zStripLoop: Loop over all strips in z-direction
    zStripLoop: do iStrip = 1, n_squared, striplen
      ! iAlph is the index of the first element in a line for the transformation in
      ! z-direction.
      do iAlph = iStrip, min(iStrip+striplen-1, n_squared)
        alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = pntVal(iAlph::n_squared)
      end do

      ! At the end of the array the number of computed strips might be smaller
      nIndeps = min(striplen, n_squared-iStrip+1)

      call ply_pntToLeg( fpt       = fpt,     &
        &                nIndeps   = nIndeps, &
        &                legCoeffs = gam,     &
        &                pntVal    = alph     )

      legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)

      ! todo: fft on temp
      ! temp -> pntVal (stride-1 writing)

    end do zStripLoop

    ! y-direction
    yStripLoop: do iStrip = 1,n_squared,striplen
      do iAlph = iStrip, min(iStrip+striplen-1, n_squared)
        alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) &
          & = legCoeffs(iAlph::n_squared)
      end do

      ! At the end of the array the number of computed strips might be smaller
      nIndeps = min(striplen, n_squared-iStrip+1)

      call ply_pntToLeg( fpt       = fpt,     &
        &                nIndeps   = nIndeps, &
        &                legCoeffs = gam,     &
        &                pntVal    = alph     )

        ! todo: fft on temp
        ! temp -> pntVal (stride-1 writing)
      pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

    end do yStripLoop ! iStrip

    ! x-direction
    xStripLoop: do iStrip = 1,n_squared,striplen
      do iAlph = iStrip, min(iStrip+striplen-1, n_squared)
        alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = pntVal(iAlph::n_squared)
      end do

      ! At the end of the array the number of computed strips might be smaller
      nIndeps = min(striplen, n_squared-iStrip+1)

      call ply_pntToLeg( fpt       = fpt,     &
        &                nIndeps   = nIndeps, &
        &                legCoeffs = gam,     &
        &                pntVal    = alph     )

      ! todo: fft on temp
      ! temp -> pntVal (stride-1 writing)

      legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n) = gam(1:nIndeps*n)

    end do xStripLoop


  end subroutine ply_pntToLeg_3D_singVar
  ! ************************************************************************ !

end module ply_legFpt_3D_module
