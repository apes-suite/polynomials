! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2014, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Nikhil Anand, Timo Stentenbach,
! Harald Klimach and Peter Vitt for University of Siegen.
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
!> ply_leg_diff_module
!!
!! This module contains the subroutine for differentiation of the legendre
!! Polynomials in 1D, 2D and 3D.
module ply_leg_diff_module
  use env_module, only: rk

  implicit none

  private

  public :: calcDiff_leg
  public :: calcDiff_leg_2d
  public :: calcDiff_leg_1d
  public :: calcDiff_leg_normal
  public :: calcDiff_leg_2d_normal

contains


  ! ************************************************************************ !
  subroutine calcDiff_leg_normal( legCoeffs, legCoeffsDiff, mPd, nVars, &
    &                             elemLength, iDir, dirVec              )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients. \n
    !! First index is the number of modal coefficients. \n
    !! Second index is the number of velocity components \n
    !! Third index is the number of partial derivatives, i.e. 3 in 3D.
    !real(kind=rk), intent(inout) :: legCoeffsDiff(:,:,:)
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: mPd
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in):: elemLength
    !> The direction to differentiate
    integer, intent(in) :: iDir
    !> The direction vector for the rotation
    integer, optional :: dirVec(3)
    ! -------------------------------------------------------------------- !
    integer :: iVar
    integer :: dofPos, dofPosPrev, dofPos2Prev
    integer :: leg(3), iDeg, iDeg1, iDeg2, iDeg3, DV(3)
    ! -------------------------------------------------------------------- !

    if (present(dirVec)) then
      DV = dirvec
    else
      if (iDir == 1) then
        DV = [3,1,2]
      elseif (iDir ==2) then
        DV = [1,3,2]
      elseif (iDir ==3) then
        DV = [1,2,3]
      endif
    endif

    do iDeg = 1, (mpd+1)**2
      iDeg1 = (iDeg-1)/(mpd+1) + 1      !! do IDeg1 = 1, mPd+1
      iDeg2 = iDeg - (iDeg1-1)*(mpd+1)  !! do IDeg2 = 1, mPd=1   !! iDeg2 = mod(iDeg-1,mpd+1)+1
      iDeg3 = mPd+1
      leg = (/iDeg1, iDeg2, iDeg3/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos )
      ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
      !                              leg(dirVec(2)), &
      !                              leg(dirVec(3)), &
      !                              maxPolyDegree   )

      !legCoeffsDiff(dofPos,:,iDir) = 0.0_rk
      legCoeffsDiff(dofPos,:) = 0.0_rk
      dofPosPrev = dofPos
      leg = (/iDeg1, iDeg2, iDeg3-1/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos )
      ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
      !                              leg(dirVec(2)), &
      !                              leg(dirVec(3)), &
      !                              maxPolyDegree   )

      !legCoeffsDiff(dofPos,:,iDir) = legCoeffs(dofPosPrev,:)
      legCoeffsDiff(dofPos,:) = legCoeffs(dofPosPrev,:)

      do iDeg3 = mPd-1, 1, -1
        leg = (/iDeg1, iDeg2, iDeg3/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos )
        ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
        !                              leg(dirVec(2)), &
        !                              leg(dirVec(3)), &
        !                              mPd   )

        leg = (/iDeg1, iDeg2, iDeg3+1/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofposPrev )
        ! dofposPrev = posOfModgCoeffQTens(leg(dirVec(1)), &
        !                                  leg(dirVec(2)), &
        !                                  leg(dirVec(3)), &
        !                                  mPd   )

        leg = (/iDeg1, iDeg2, iDeg3+2/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos2Prev )
        ! dofpos2Prev = posOfModgCoeffQTens(leg(dirVec(1)), &
        !                                   leg(dirVec(2)), &
        !                                   leg(dirVec(3)), &
        !                                   mPd   )

        do iVar = 1, nVars
          !legCoeffsDiff(dofPos, iVar, iDir) = legCoeffsDiff(dofPos2Prev,iVar,iDir) &
          legCoeffsDiff(dofPos, iVar) = legCoeffsDiff(dofPos2Prev,iVar) &
          &                               + legCoeffs(dofPosPrev, iVar)
        end do
      end do
    end do

    ! Scale the results due to the Jacobians of the mappings
    do dofpos=1,(mpd+1)**3
      ideg3 = (dofpos-1)/(mpd+1)**2 + 1
      iDeg = dofpos - (ideg3-1)*(mpd+1)**2
      iDeg2 = (iDeg-1)/(mpd+1) + 1
      iDeg1 = mod(dofpos-1, mpd+1)  + 1
      leg = (/iDeg1, iDeg2, iDeg3/)
      legCoeffsDiff(dofPos,:) = legCoeffsDiff(dofPos,:)         &
        &                         * (2.0_rk/elemLength)         &
        &                         * (2.0_rk*leg(iDir) - 1.0_rk)
    end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Uncollapsed version of the scaling !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!   do iDeg1 = 1, mPd+1
!!     do iDeg2 = 1, mPd+1
!!       do iDeg3 = 1, mPd+1
!!         leg = (/iDeg1, iDeg2, iDeg3/)
!!         dofPos = 1 + (iDeg1-1)                                   &
!!         &      + (iDeg2-1)*(mPd+1)                               &
!!         &      + (iDeg3-1)*(mPd+1)*(mPd+1)
!!         legCoeffsDiff(dofPos,:,iDir) = legCoeffsDiff(dofPos,:,iDir)       &
!!         &                       * (2.0_rk/elemLength)                     &
!!         &                       * (2.0_rk*leg(iDir) - 1.0_rk)
!!       end do
!!     end do
!!   end do

  end subroutine calcDiff_leg_normal
  ! ************************************************************************ !

! ************************************************************************ !
  subroutine calcDiff_leg_normal_vec( legCoeffs, legCoeffsDiff, mPd, nVars, &
    &                             elemLength, iDir, dirVec              )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients. \n
    !! First index is the number of modal coefficients. \n
    !! Second index is the number of velocity components \n
    !! Third index is the number of partial derivatives, i.e. 3 in 3D.
    !real(kind=rk), intent(inout) :: legCoeffsDiff(:,:,:)
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: mPd
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in):: elemLength
    !> The direction to differentiate
    integer, intent(in) :: iDir
    !> The direction vector for the rotation
    integer, optional :: dirVec(3)
    ! -------------------------------------------------------------------- !
    integer :: iVar
    integer :: dofPos, dofPosPrev, dofPos2Prev
    integer :: leg(3), iDeg,iDeg1, iDeg2, iDeg3, DV(3)
    integer :: legprev(3)
    ! -------------------------------------------------------------------- !

    if (present(dirVec)) then
      DV = dirvec
    else
      if (iDir == 1) then
        DV = [3,1,2]
      elseif (iDir == 2) then
        DV = [1,3,2]
      elseif (iDir == 3) then
        DV = [1,2,3]
      endif
    endif

    varloop: do iVar = 1, nVars

      iDeg3 = mpd+1
      do iDeg = 1, (mpd+1)**2
        iDeg1 = (iDeg-1)/(mpd+1) + 1      !! do IDeg1 = 1, mPd+1
        iDeg2 = iDeg - (iDeg1-1)*(mpd+1)  !! do IDeg2 = 1, mPd=1   !! iDeg2 = mod(iDeg-1,mpd+1)+1

        leg = (/iDeg1, iDeg2, iDeg3/)
?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos )
        ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
        !                              leg(dirVec(2)), &
        !                              leg(dirVec(3)), &
        !                              maxPolyDegree   )
        legCoeffsDiff(dofPos,iVar) = 0.0_rk

        legprev = (/iDeg1, iDeg2, iDeg3-1/)
?? copy :: posOfModgCoeffQTens(legprev(DV(1)), legprev(DV(2)), legprev(DV(3)), mPd, dofposPrev )
        ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
        !                              leg(dirVec(2)), &
        !                              leg(dirVec(3)), &
        !                              maxPolyDegree   )

        legCoeffsDiff(dofPosPrev,iVar) = legCoeffs(dofPos,iVar)

      end do

      do iDeg3 = mPd-1, 1, -1

        !NEC$ ivdep
        do iDeg = 1, (mpd+1)**2
          iDeg1 = (iDeg-1)/(mpd+1) + 1      !! do IDeg1 = 1, mPd+1
          iDeg2 = iDeg - (iDeg1-1)*(mpd+1)  !! do IDeg2 = 1, mPd=1   !! iDeg2 = mod(iDeg-1,mpd+1)+1

          leg = (/iDeg1, iDeg2, iDeg3/)
?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos )
          ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
          !                              leg(dirVec(2)), &
          !                              leg(dirVec(3)), &
          !                              mPd   )

          leg = (/iDeg1, iDeg2, iDeg3+1/)
?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofposPrev )
          ! dofposPrev = posOfModgCoeffQTens(leg(dirVec(1)), &
          !                                  leg(dirVec(2)), &
          !                                  leg(dirVec(3)), &
          !                                  mPd   )

          leg = (/iDeg1, iDeg2, iDeg3+2/)
?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), leg(DV(3)), mPd, dofpos2Prev )
          ! dofpos2Prev = posOfModgCoeffQTens(leg(dirVec(1)), &
          !                                   leg(dirVec(2)), &
          !                                   leg(dirVec(3)), &
          !                                   mPd   )

          legCoeffsDiff(dofPos, iVar) = legCoeffsDiff(dofPos2Prev,iVar) &
            &                         + legCoeffs(dofPosPrev, iVar)
        end do
      end do

      ! Scale the results due to the Jacobians of the mappings
      do dofpos=1,(mpd+1)**3
        ideg3 = (dofpos-1)/(mpd+1)**2 + 1
        iDeg = dofpos - (ideg3-1)*(mpd+1)**2
        iDeg2 = (iDeg-1)/(mpd+1) + 1
        iDeg1 = mod(dofpos-1, mpd+1)  + 1
        leg = (/iDeg1, iDeg2, iDeg3/)
        legCoeffsDiff(dofPos,iVar) = legCoeffsDiff(dofPos,iVar) &
          &                         * (2.0_rk/elemLength)       &
          &                         * (2.0_rk*leg(iDir) - 1.0_rk)
      end do

    end do varloop


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Uncollapsed version of the scaling !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!   do iDeg1 = 1, mPd+1
!!     do iDeg2 = 1, mPd+1
!!       do iDeg3 = 1, mPd+1
!!         leg = (/iDeg1, iDeg2, iDeg3/)
!!         dofPos = 1 + (iDeg1-1)                                   &
!!         &      + (iDeg2-1)*(mPd+1)                               &
!!         &      + (iDeg3-1)*(mPd+1)*(mPd+1)
!!         legCoeffsDiff(dofPos,:,iDir) = legCoeffsDiff(dofPos,:,iDir)       &
!!         &                       * (2.0_rk/elemLength)                     &
!!         &                       * (2.0_rk*leg(iDir) - 1.0_rk)
!!       end do
!!     end do
!!   end do

  end subroutine calcDiff_leg_normal_vec
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Compute the derivative in X direction for 3D Legendre polynomial.
  subroutine calcDiff_leg_x_vec( legCoeffs, legCoeffsDiff, mPd, nVars, &
    &                             elemLength )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients.
    !! * First index is the number of modal coefficients.
    !! * Second index is the number of variable components
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: mPd
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in):: elemLength
    ! -------------------------------------------------------------------- !
    integer :: iVar
    integer :: dofPos, dofPosPrev
    integer :: iDeg, iDeg1, iDeg2, iDeg3
    real(kind=rk) :: derivative((mpd+1)**2,mpd+1)
    ! -------------------------------------------------------------------- !

    varloop: do iVar = 1, nVars

      derivative(:,mpd+1) = 0.0_rk

      iDeg1 = mpd

      do iDeg = 1, (mpd+1)**2
        iDeg3 = (iDeg-1)/(mpd+1) + 1
        iDeg2 = iDeg - (iDeg3-1)*(mpd+1)

?? copy :: posOfModgCoeffQTens(iDeg1+1, iDeg2, iDeg3, mPd, dofpos )
        derivative(iDeg, mpd) = legCoeffs(dofpos,iVar)
      end do

      do iDeg1 = mPd-1, 1, -1

        !NEC$ ivdep
        do iDeg = 1, (mpd+1)**2
          iDeg3 = (iDeg-1)/(mpd+1) + 1
          iDeg2 = iDeg - (iDeg3-1)*(mpd+1)

?? copy :: posOfModgCoeffQTens(iDeg1+1, iDeg2, iDeg3, mPd, dofposPrev)

          derivative(iDeg, iDeg1) = derivative(iDeg, iDeg1+2) &
            &                       + legCoeffs(dofposprev, iVar)
        end do
      end do

      ! Scale the results due to the Jacobians of the mappings
      do dofpos=1,(mpd+1)**3
        ideg3 = (dofpos-1)/(mpd+1)**2 + 1
        iDeg = dofpos - (ideg3-1)*(mpd+1)**2
        iDeg2 = (iDeg-1)/(mpd+1) + 1
        iDeg1 = mod(dofpos-1, mpd+1)  + 1
        legCoeffsDiff(dofPos,iVar) = derivative(iDeg2+(iDeg3-1)*(mpd+1),iDeg1) &
          &                          * (2.0_rk/elemLength)  &
          &                          * (2.0_rk*iDeg1 - 1.0_rk)
      end do

    end do varloop

  end subroutine calcDiff_leg_x_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compute the derivative in Y direction for 3D Legendre polynomial.
  subroutine calcDiff_leg_y_vec( legCoeffs, legCoeffsDiff, mPd, nVars, &
    &                             elemLength )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients.
    !! * First index is the number of modal coefficients.
    !! * Second index is the number of variable components
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: mPd
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in):: elemLength
    ! -------------------------------------------------------------------- !
    integer :: iVar
    integer :: dofPos, dofPosPrev
    integer :: iDeg, iDeg1, iDeg2, iDeg3
    real(kind=rk) :: derivative((mpd+1)**2,mpd+1)
    ! -------------------------------------------------------------------- !

    varloop: do iVar = 1, nVars

      derivative(:,mpd+1) = 0.0_rk

      iDeg2 = mpd

      do iDeg = 1, (mpd+1)**2
        iDeg3 = (iDeg-1)/(mpd+1) + 1
        iDeg1 = iDeg - (iDeg3-1)*(mpd+1)

?? copy :: posOfModgCoeffQTens(iDeg1, iDeg2+1, iDeg3, mPd, dofpos )
        derivative(iDeg, mpd) = legCoeffs(dofpos,iVar)
      end do

      do iDeg2 = mPd-1, 1, -1

        !NEC$ ivdep
        do iDeg = 1, (mpd+1)**2
          iDeg3 = (iDeg-1)/(mpd+1) + 1
          iDeg1 = iDeg - (iDeg3-1)*(mpd+1)

?? copy :: posOfModgCoeffQTens(iDeg1, iDeg2+1, iDeg3, mPd, dofposPrev)

          derivative(iDeg, iDeg2) = derivative(iDeg, iDeg2+2) &
            &                       + legCoeffs(dofposprev, iVar)
        end do
      end do

      ! Scale the results due to the Jacobians of the mappings
      do dofpos=1,(mpd+1)**3
        ideg3 = (dofpos-1)/(mpd+1)**2 + 1
        iDeg = dofpos - (ideg3-1)*(mpd+1)**2
        iDeg2 = (iDeg-1)/(mpd+1) + 1
        iDeg1 = mod(dofpos-1, mpd+1)  + 1
        legCoeffsDiff(dofPos,iVar) = derivative(iDeg1+(iDeg3-1)*(mpd+1),iDeg2) &
          &                          * (2.0_rk/elemLength)  &
          &                          * (2.0_rk*iDeg2 - 1.0_rk)
      end do

    end do varloop

  end subroutine calcDiff_leg_y_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Compute the derivative in Y direction for 3D Legendre polynomial.
  subroutine calcDiff_leg_z_vec( legCoeffs, legCoeffsDiff, mPd, nVars, &
    &                             elemLength )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients.
    !! * First index is the number of modal coefficients.
    !! * Second index is the number of variable components
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: mPd
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in):: elemLength
    ! -------------------------------------------------------------------- !
    integer :: iVar
    integer :: dofPos, dofPosPrev
    integer :: iDeg, iDeg1, iDeg2, iDeg3
    real(kind=rk) :: derivative((mpd+1)**2,mpd+1)
    ! -------------------------------------------------------------------- !

    varloop: do iVar = 1, nVars

      derivative(:,mpd+1) = 0.0_rk

      iDeg3 = mpd

      do iDeg = 1, (mpd+1)**2
        iDeg2 = (iDeg-1)/(mpd+1) + 1
        iDeg1 = iDeg - (iDeg2-1)*(mpd+1)

?? copy :: posOfModgCoeffQTens(iDeg1, iDeg2, iDeg3+1, mPd, dofpos )
        derivative(iDeg, mpd) = legCoeffs(dofpos,iVar)
      end do

      do iDeg3 = mPd-1, 1, -1

        !NEC$ ivdep
        do iDeg = 1, (mpd+1)**2
          iDeg2 = (iDeg-1)/(mpd+1) + 1
          iDeg1 = iDeg - (iDeg2-1)*(mpd+1)

?? copy :: posOfModgCoeffQTens(iDeg1, iDeg2, iDeg3+1, mPd, dofposPrev)

          derivative(iDeg, iDeg3) = derivative(iDeg, iDeg3+2) &
            &                       + legCoeffs(dofposprev, iVar)
        end do
      end do

      ! Scale the results due to the Jacobians of the mappings
      do dofpos=1,(mpd+1)**3
        ideg3 = (dofpos-1)/(mpd+1)**2 + 1
        ideg = dofpos - (ideg3-1)*(mpd+1)**2
        legCoeffsDiff(dofPos,iVar) = derivative(iDeg,iDeg3) &
          &                          * (2.0_rk/elemLength)  &
          &                          * (2.0_rk*iDeg3 - 1.0_rk)
      end do

    end do varloop

  end subroutine calcDiff_leg_z_vec
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine calcDiff_leg_2d_normal( legCoeffs, legCoeffsDiff, mPd, nVars, &
    &                                elemLength, iDir, dirVec              )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients. \n
    !! First index is the number of modal coefficients. \n
    !! Second index is the number of velocity components \n
    !! Third index is the number of partial derivatives, i.e. 3 in 3D.
    !real(kind=rk), intent(inout) :: legCoeffsDiff(:,:,:)
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: mPd
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in):: elemLength
    !> The direction to differentiate
    integer, intent(in) :: iDir
    !> The direction vector for the rotation
    integer, optional :: dirVec(2)
    ! -------------------------------------------------------------------- !
    integer :: iVar
    integer :: dofPos, dofPosPrev, dofPos2Prev
    integer :: leg(2), iDeg1, iDeg2, DV(2)
    ! -------------------------------------------------------------------- !

    if (present(dirVec)) then
      DV = dirvec
    else
      if (iDir == 1) then
        DV = [2,1]
      elseif (iDir ==2) then
        DV = [1,2]
      endif
    endif

    do iDeg1 = 1, mPd+1
      iDeg2 =  mPd+1
      leg = (/iDeg1, iDeg2/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)),1,  mPd, dofpos )
       ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
       !                              leg(dirVec(2)), &
       !                              leg(dirVec(3)), &
       !                              maxPolyDegree   )

      !legCoeffsDiff(dofPos,:,iDir) = 0.0_rk
      legCoeffsDiff(dofPos,:) = 0.0_rk
      dofPosPrev = dofPos
      leg = (/iDeg1, iDeg2-1/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), 1, mPd, dofpos )
       ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
       !                              leg(dirVec(2)), &
       !                              leg(dirVec(3)), &
       !                              maxPolyDegree   )

      !legCoeffsDiff(dofPos,:,iDir) = legCoeffs(dofPosPrev,:)
      legCoeffsDiff(dofPos,:) = legCoeffs(dofPosPrev,:)

      do iDeg2 = mPd-1, 1, -1
        leg = (/iDeg1, iDeg2/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), 1, mPd, dofpos )
         ! dofpos = posOfModgCoeffQTens(leg(dirVec(1)), &
         !                              leg(dirVec(2)), &
         !                              leg(dirVec(3)), &
         !                              mPd   )

        leg = (/iDeg1, iDeg2+1/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), 1, mPd, dofposPrev )
         ! dofposPrev = posOfModgCoeffQTens(leg(dirVec(1)), &
         !                                  leg(dirVec(2)), &
         !                                  leg(dirVec(3)), &
         !                                  mPd   )

        leg = (/iDeg1, iDeg2+2/)

?? copy :: posOfModgCoeffQTens(leg(DV(1)), leg(DV(2)), 1, mPd, dofpos2Prev )
         ! dofpos2Prev = posOfModgCoeffQTens(leg(dirVec(1)), &
         !                                   leg(dirVec(2)), &
         !                                   leg(dirVec(3)), &
         !                                   mPd   )

        do iVar = 1, nVars
          !legCoeffsDiff(dofPos, iVar,iDir) = legCoeffsDiff(dofPos2Prev,iVar,iDir) &
          legCoeffsDiff(dofPos, iVar) = legCoeffsDiff(dofPos2Prev,iVar) &
            &                         + legCoeffs(dofPosPrev, iVar)
        end do
      end do
    end do

    ! Scale the results due to the Jacobians of the mappings
    do iDeg1 = 1, mPd+1
      do iDeg2 = 1, mPd+1
        leg = (/iDeg1, iDeg2/)
        dofPos = 1 + (iDeg1-1) + (iDeg2-1)*(mPd+1)
        !legCoeffsDiff(dofPos,:,iDir) = legCoeffsDiff(dofPos,:,iDir)         &
        legCoeffsDiff(dofPos,:) = legCoeffsDiff(dofPos,:)         &
          &                         * (2.0_rk/elemLength)         &
          &                         * (2.0_rk*leg(iDir) - 1.0_rk)
      end do
    end do

  end subroutine calcDiff_leg_2d_normal
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine calcDiff_leg( legCoeffs, legCoeffsDiff, maxPolyDegree, nVars, &
    &                      elemLength                                      )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients. \n
    !! First index is the number of modal coefficients. \n
    !! Second index is the number of velocity components \n
    !! Third index is the number of partial derivatives, i.e. 3 in 3D.
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:,:)
    integer, intent(in) :: maxPolyDegree
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk),intent(in) :: elemLength
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    if (maxpolydegree > 0) then
      ! Loop over Directions
      call calcDiff_leg_x_vec( legCoeffs     = legCoeffs,            &
        &                      legCoeffsDiff = legCoeffsDiff(:,:,1), &
        &                      mPD           = maxPolyDegree,        &
        &                      nVars         = nVars,                &
        &                      elemLength    = elemLength            )
      call calcDiff_leg_y_vec( legCoeffs     = legCoeffs,            &
        &                      legCoeffsDiff = legCoeffsDiff(:,:,2), &
        &                      mPD           = maxPolyDegree,        &
        &                      nVars         = nVars,                &
        &                      elemLength    = elemLength            )
      call calcDiff_leg_z_vec( legCoeffs     = legCoeffs,            &
        &                      legCoeffsDiff = legCoeffsDiff(:,:,3), &
        &                      mPD           = maxPolyDegree,        &
        &                      nVars         = nVars,                &
        &                      elemLength    = elemLength            )
    else
      legCoeffsDiff = 0.0_rk
    end if

  end subroutine calcDiff_leg
  ! ************************************************************************ !
  

  ! ************************************************************************ !
  subroutine calcDiff_leg_2d( legCoeffs, legCoeffsDiff, maxPolyDegree, nVars, &
    &                         elemLength                                      )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients. \n
    !! First index is the number of modal coefficients. \n
    !! Second index is the number of velocity components \n
    !! Third index is the number of partial derivatives, i.e. 3 in 3D.
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:,:)
    integer, intent(in) :: maxPolyDegree
    !> The number of varibales to differentiate
    integer, intent(in) :: nVars
    !> The physical length of the element to build the derivatives for.
    real(kind=rk), intent(in) :: elemLength
    ! -------------------------------------------------------------------- !
    integer :: dirvec(2,2), iDir
    ! -------------------------------------------------------------------- !

    if (maxpolydegree > 0) then
      dirvec(:,1) = [2, 1]
      dirvec(:,2) = [1, 2]
      ! Loop over Directions
      do iDir = 1,2
        ! Calculate the differentiation for the particular direction
        call calcDiff_leg_2d_normal( legCoeffs = legCoeffs,                   &
          &                          legCoeffsDiff = legCoeffsDiff(:,:,iDir), &
          &                          mPD           = maxPolyDegree,           &
          &                          nVars         = nVars,                   &
          &                          elemLength    = elemLength,              &
          &                          dirvec        = dirvec(:,iDir),          &
          &                          iDir          = iDir                     )
      enddo
    endif

  end subroutine calcDiff_leg_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine calcDiff_leg_1d( legCoeffs, legCoeffsDiff, maxPolyDegree, &
    &                         elemLength                               )
    ! -------------------------------------------------------------------- !
    real(kind=rk), intent(in) :: legCoeffs(:,:)
    !> Modal expansion of the derivative of legCoeffs in terms of Legendre
    !! modal coefficients. \n
    !! First index is the number of modal coefficients. \n
    !! Second index is the number of var components \n
    real(kind=rk), intent(inout) :: legCoeffsDiff(:,:)
    integer, intent(in) :: maxPolyDegree
    !> The physical length of the element to build the derivatives for.
    real(kind=rk),intent(in) :: elemLength
    ! -------------------------------------------------------------------- !
    integer :: iDegX
    integer :: dofPos, dofPosPrev, dofPos2Prev
    ! -------------------------------------------------------------------- !

    ! Build the derivative in x direction
    dofPos = 1 + maxPolyDegree
    legCoeffsDiff(dofPos,:) = 0.0_rk
    if (maxpolydegree > 0) then
      dofPosPrev = dofPos
      dofPos = 1 + (maxPolyDegree-1)
      legCoeffsDiff(dofPos,:) = legCoeffs(dofPosPrev,:)
      do iDegX = maxPolyDegree-1, 1, -1
        dofPos = 1 + (iDegX-1)
        dofPosPrev = 1 + (iDegX)
        dofPos2Prev = 1 + (iDegX+1)
        legCoeffsDiff(dofPos,:) = legCoeffsDiff(dofPos2Prev,:) &
          &                         + legCoeffs(dofPosPrev,:)
      end do
    end if

    do iDegX = 1, maxPolyDegree+1
      dofPos = 1 + (iDegX-1)
      legCoeffsDiff(dofPos,:) = legCoeffsDiff(dofPos,:)     &
        &                         * (2.0_rk/elemLength)     &
        &                         * (2.0_rk*iDegX - 1.0_rk)
    end do

  end subroutine calcDiff_leg_1d
  ! ************************************************************************ !


end module ply_leg_diff_module
