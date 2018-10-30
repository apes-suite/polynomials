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
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iDeg1, iDeg2, iDeg3, iDeg, iVar, leg, dofPos)

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

    !$OMP DO
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
    !$OMP END DO

    ! Scale the results due to the Jacobians of the mappings
    !$OMP DO
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
    !$OMP END DO


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
    !$OMP END PARALLEL

  end subroutine calcDiff_leg_normal
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
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iDeg1, iDeg2, leg, dofPos, dofPosPrev, iVar)

    if (present(dirVec)) then
      DV = dirvec
    else
      if (iDir == 1) then
        DV = [2,1]
      elseif (iDir ==2) then
        DV = [1,2]
      endif
    endif

    !$OMP DO
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
    !$OMP END DO


    ! Scale the results due to the Jacobians of the mappings
    !$OMP DO
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
    !$OMP END DO
    !$OMP END PARALLEL


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
    !> The direction to evaluate the differentiation
    integer :: iDir
    ! -------------------------------------------------------------------- !
    integer :: dirvec(3,3)
    ! -------------------------------------------------------------------- !
    if (maxpolydegree > 0) then
      dirvec(:,1) = [3, 1, 2]
      dirvec(:,2) = [1, 3, 2]
      dirvec(:,3) = [1, 2, 3]

      ! Loop over Directions
      do iDir = 1,3
        ! Calculate the differentiation for the particular direction
        call calcDiff_leg_normal( legCoeffs     = legCoeffs,               &
          &                       legCoeffsDiff = legCoeffsDiff(:,:,iDir), &
          &                       mPD           = maxPolyDegree,           &
          &                       nVars         = nVars,                   &
          &                       elemLength    = elemLength,              &
          &                       dirvec        = dirvec(:,iDir),          &
          &                       iDir          = iDir                     )
      enddo
    endif
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
      end do
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
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(dofPos, iDegX, dofPosPrev, dofPos2Prev)

    ! Build the derivative in x direction
    dofPos = 1 + maxPolyDegree
    legCoeffsDiff(dofPos,:) = 0.0_rk
    if (maxpolydegree > 0) then
      dofPosPrev = dofPos
      dofPos = 1 + (maxPolyDegree-1)
      legCoeffsDiff(dofPos,:) = legCoeffs(dofPosPrev,:)
      !$OMP DO
      do iDegX = maxPolyDegree-1, 1, -1
        dofPos = 1 + (iDegX-1)
        dofPosPrev = 1 + (iDegX)
        dofPos2Prev = 1 + (iDegX+1)
        legCoeffsDiff(dofPos,:) = legCoeffsDiff(dofPos2Prev,:) &
          &                         + legCoeffs(dofPosPrev,:)
      end do
      !$OMP END DO
    end if

    !$OMP DO
    do iDegX = 1, maxPolyDegree+1
      dofPos = 1 + (iDegX-1)
      legCoeffsDiff(dofPos,:) = legCoeffsDiff(dofPos,:)     &
        &                         * (2.0_rk/elemLength)     &
        &                         * (2.0_rk*iDegX - 1.0_rk)
    end do
    !$OMP END DO
    !$OMP END PARALLEL


  end subroutine calcDiff_leg_1d
  ! ************************************************************************ !


end module ply_leg_diff_module
