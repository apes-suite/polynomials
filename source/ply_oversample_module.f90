!> This module provides functions to transfer polynomials from and to the
!! oversampled representation for nodal treatments.
module ply_oversample_module

  use env_module,               only: rk, labelLen
  use ply_dof_module,           only: posOfModgCoeffPTens,   &
    &                                 posOfModgCoeffPTens2D, &
    &                                 nextModgCoeffPTens,    &
    &                                 nextModgCoeffPTens2D,  &
    &                                 Q_space, P_space
  use ply_poly_project_module,  only: ply_poly_project_type

  implicit none

  private

  public :: ply_convert2oversample
  public :: ply_convertFromOversample

contains


  !****************************************************************************!
  !> Copy a single element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  !! Please note, that the oversampled representation is always a Q-Space
  !! polynomial.
  subroutine ply_convert2oversample(state, poly_proj, nDim, modalCoeffs, &
    &                               nScalars                             )
    !--------------------------------------------------------------------------
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> The number of dimensions to determine the correct oversampling routine.
    integer, intent(in) :: nDim

    !> Oversampled array for modal coefficients from state.
    !!
    !! These are always Q-Polynomial representations, as projections only work
    !! with those.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    !--------------------------------------------------------------------------
    select case(nDim)
    case(1)
      call ply_convert2oversample_1d(state, poly_proj, modalCoeffs, nScalars)
    case(2)
      call ply_convert2oversample_2d(state, poly_proj, modalCoeffs, nScalars)
    case(3)
      call ply_convert2oversample_3d(state, poly_proj, modalCoeffs)
    end select
  end subroutine ply_convert2oversample

  !****************************************************************************!
  !> Truncating an oversampled polynomial representation back to the original
  !! representation.
  !!
  !! Note, that the oversampled array, which is given as an input here is
  !! is always a Q-Polynomial representation, while the target truncated state
  !! might also be a P-Polynomial.
  subroutine ply_convertFromOversample(modalCoeffs, poly_proj, nDim, state, &
      &                                nScalars                             )
    !--------------------------------------------------------------------------
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> The number of dimensions to determine the correct oversampling routine.
    integer, intent(in) :: nDim

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(inout) :: state(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    !--------------------------------------------------------------------------
    select case(nDim)
    case(1)
      call ply_convertFromOversample_1d(modalCoeffs, poly_proj, state, nScalars)
    case(2)
      call ply_convertFromOversample_2d(modalCoeffs, poly_proj, state, nScalars)
    case(3)
      call ply_convertFromOversample_3d(modalCoeffs, poly_proj, state)
    end select
  end subroutine ply_convertFromOversample

  !****************************************************************************!
  !> Copy a single element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  !! Please note, that the oversampled representation is always a Q-Space
  !! polynomial.
  subroutine ply_convert2oversample_3d(state, poly_proj, modalCoeffs)
    !--------------------------------------------------------------------------
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> Oversampled array for modal coefficients from state.
    !!
    !! These are always Q-Polynomial representations, as projections only work
    !! with those.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)
    !--------------------------------------------------------------------------
    integer :: pos
    integer :: oversamp_degree
    integer :: mpd1, mpd1_square, mpd1_cube
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp
    !--------------------------------------------------------------------------
    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2
    mpd1_cube = mpd1**3

    !$OMP WORKSHARE
    modalCoeffs(:,:) = 0.0_rk
    !$OMP END WORKSHARE

    if (poly_proj%basisType == Q_Space) then
      !$OMP DO
      do dof = 1, mpd1_cube
        iDegZ = (dof-1)/mpd1_square + 1
        iDegY = (dof-1-(iDegZ-1)*mpd1_square)/mpd1+1
        iDegX = mod(dof-1,mpd1)+1
        dofOverSamp = iDegX + ( iDegY-1  &
          &                     + (iDegZ-1)*(oversamp_degree+1) &
          &                   ) * (oversamp_degree+1)
        modalCoeffs(dofOverSamp,:) = state(dof,:)
      end do
      !$OMP END DO

    else !P_Space

      iDegX = 1
      iDegY = 1
      iDegZ = 1
      !$OMP SINGLE
      do idof = 1, poly_proj%body_3d%min_dofs
        dof = posOfModgCoeffPTens(iDegX, iDegY, iDegZ, poly_proj%maxPolyDegree)
        dofOverSamp = iDegX + ( iDegY-1  &
          &                     + (iDegZ-1)*(oversamp_degree+1) &
          &                   ) * (oversamp_degree+1)
        modalCoeffs(dofOverSamp,:) = state(dof,:)
        call nextModgCoeffPTens(iDegX, iDegY, iDegZ, poly_proj%min_degree)
      end do
      !$OMP END SINGLE
    end if

  end subroutine ply_convert2oversample_3d
  !****************************************************************************!


  !****************************************************************************!
  !> Truncating an oversampled polynomial representation back to the original
  !! representation.
  !!
  !! Note, that the oversampled array, which is given as an input here is
  !! is always a Q-Polynomial representation, while the target truncated state
  !! might also be a P-Polynomial.
  subroutine ply_convertFromOversample_3d(modalCoeffs, poly_proj, state)
    !--------------------------------------------------------------------------
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(inout) :: state(:,:)
    !--------------------------------------------------------------------------
    integer :: pos
    integer :: oversamp_degree
    integer :: mpd1, mpd1_square, mpd1_cube
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp
    !--------------------------------------------------------------------------

    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2
    mpd1_cube = mpd1**3

    if (poly_proj%basisType == Q_Space) then
      !$OMP DO
      do dof = 1, mpd1_cube
        iDegZ = (dof-1)/mpd1_square + 1
        iDegY = (dof-1-(iDegZ-1)*mpd1_square)/mpd1+1
        iDegX = mod(dof-1,mpd1)+1
        dofOverSamp = iDegX + ( iDegY-1  &
          &                     + (iDegZ-1)*(oversamp_degree+1) &
          &                   ) * (oversamp_degree+1)
        state(dof,:) = modalCoeffs(dofOverSamp,:)
      end do
      !$OMP END DO

    else !P_Space

      iDegX = 1
      iDegY = 1
      iDegZ = 1
      !$OMP SINGLE
      do idof = 1, poly_proj%body_3d%min_dofs
        dof = posOfModgCoeffPTens(iDegX, iDegY, iDegZ, poly_proj%maxPolyDegree)
        dofOverSamp = iDegX + ( iDegY-1  &
          &                     + (iDegZ-1)*(oversamp_degree+1) &
          &                   ) * (oversamp_degree+1)
        state(dof,:) = modalCoeffs(dofOverSamp,:)
        call nextModgCoeffPTens(iDegX, iDegY, iDegZ, poly_proj%min_degree)
      end do
      !$OMP END SINGLE
    end if

  end subroutine ply_convertFromoversample_3d
  !****************************************************************************!


  !****************************************************************************!
  !> Copy a single 2D element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  !! Please note, that the oversampled representation is always a Q-Space
  !! polynomial.
  subroutine ply_convert2oversample_2d(state, poly_proj, modalCoeffs, nScalars)
    !--------------------------------------------------------------------------
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> Oversampled array for modal coefficients from state.
    !!
    !! These are always Q-Polynomial representations, as projections only work
    !! with those.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    !--------------------------------------------------------------------------
    integer :: pos
    integer :: oversamp_degree
    integer :: mpd1, mpd1_square
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp, nPVars
    !--------------------------------------------------------------------------
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = nScalars
    else
      nPVars = size(state,2)
    end if

    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2

    ! Initialize oversampled space correct to 0
    !$OMP WORKSHARE
    modalCoeffs(:,:) = 0.0_rk
    !$OMP END WORKSHARE

    if (poly_proj%basisType == Q_Space) then
      !$OMP DO
      do dof = 1, mpd1_square
        iDegX = mod(dof-1,mpd1)+1
        iDegY = (dof-1)/mpd1+1
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(oversamp_degree+1)
        modalCoeffs(dofOverSamp,1:nPVars) = state(dof,1:nPVars)
      end do
      !$OMP END DO

    else !P_Space

      !$OMP SINGLE
      iDegX = 1
      iDegY = 1
      iDegZ = 0 ! not used in posOfModgCoeffPTens_2D, nextModgCoeffPTens
      do idof = 1, poly_proj%body_2d%min_dofs
        dof = posOfModgCoeffPTens2D(iDegX, iDegY, iDegZ, &
          &                         poly_proj%maxPolydegree)
        dofOverSamp = iDegX + (iDegY-1)*(oversamp_degree+1)
        modalCoeffs(dofOverSamp,1:nPVars) = state(dof,1:nPVars)
        call nextModgCoeffPTens2D(iDegX, iDegY, iDegZ, poly_proj%min_degree)
      end do
      !$OMP END SINGLE

    end if

  end subroutine ply_convert2oversample_2d
  !****************************************************************************!


  !****************************************************************************!
  !> Truncating an oversampled 2D polynomial representation back to the original
  !! representation.
  !!
  !! Note, that the oversampled array, which is given as an input here is
  !! is always a Q-Polynomial representation, while the target truncated state
  !! might also be a P-Polynomial.
  subroutine ply_convertFromOversample_2d(modalCoeffs, poly_proj, state, nScalars)
    !--------------------------------------------------------------------------
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(inout) :: state(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    !--------------------------------------------------------------------------
    integer :: pos
    integer :: oversamp_degree
    integer :: mpd1, mpd1_square
    integer :: iDegX, iDegY, iDegZ, idof, dof, dofOverSamp, nPVars
    !--------------------------------------------------------------------------
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = nScalars
    else
      nPVars = size(state,2)
    end if

    ! Information for the oversampling loop
    oversamp_degree = poly_proj%oversamp_degree
    mpd1 = poly_proj%min_degree + 1
    mpd1_square = mpd1**2

    if (poly_proj%basisType == Q_Space) then

      !$OMP DO
      do dof = 1, mpd1_square
       iDegX = mod(dof-1,mpd1)+1
       iDegY = (dof-1)/mpd1+1
       dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(oversamp_degree+1)
       state(dof,1:nPVars) = modalCoeffs(dofOverSamp,1:nPVars)
      end do
      !$OMP END DO

    else !P_Space

      !$OMP SINGLE
      iDegX = 1
      iDegY = 1
      iDegZ = 0 ! not used in posOfModgCoeffPTens_2D, nextModgCoeffPTens
      do idof = 1, poly_proj%body_2d%min_dofs
        dof = posOfModgCoeffPTens2D(iDegX, iDegY, iDegZ, &
          &                         poly_proj%maxPolydegree)
        dofOverSamp = iDegX + (iDegY-1)*(oversamp_degree+1)
        state(dof,1:nPVars) = modalCoeffs(dofOverSamp,1:nPVars)
        call nextModgCoeffPTens2D(iDegX, iDegY, iDegZ, poly_proj%min_degree)
      end do
      !$OMP END SINGLE

    end if

  end subroutine ply_convertFromoversample_2d
  !****************************************************************************!


  !****************************************************************************!
  !> Copy a single 1D element state into a larger array and pad it with zeroes.
  !!
  !! The additional modes might be necessary for proper dealing with
  !! nonlinear operations.
  subroutine ply_convert2oversample_1d(state, poly_proj, modalCoeffs, nScalars)
    !--------------------------------------------------------------------------
    !> Description of the projection method.
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> State in a single element, to be oversampled.
    real(kind=rk), intent(in) :: state(:,:)

    !> Oversampled array for modal coefficients from state.
    real(kind=rk), intent(inout) :: modalCoeffs(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    !--------------------------------------------------------------------------
    integer :: iVar, iPoint, iVP, nPVars
    !--------------------------------------------------------------------------
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = (poly_proj%maxPolyDegree+1)*nScalars
    else
      nPVars = (poly_proj%maxPolyDegree+1)*size(state,2)
    end if

    ! Initialize oversampled space correct to 0
    !$OMP WORKSHARE
    ModalCoeffs(:,:) = 0.0_rk
    !$OMP END WORKSHARE

    !$OMP DO
    do iVP = 1,nPVars
      iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
      iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
      ModalCoeffs(iPoint,iVar) = state(iPoint,iVar)
    end do
    !$OMP END DO

  end subroutine ply_convert2oversample_1d
  !****************************************************************************!


  !****************************************************************************!
  !> Truncating an oversampled 1D polynomial representation back to the original
  !! representation.
  subroutine ply_convertFromOversample_1d(modalCoeffs, poly_proj, state, nScalars)
    !--------------------------------------------------------------------------
    !> Data of the projection method
    type(ply_poly_project_type), intent(in) :: poly_proj

    !> Oversampled modal array for one element
    real(kind=rk), intent(in) :: modalCoeffs(:,:)

    !> Truncated state for one element obtained from the modalCoeffs
    real(kind=rk), intent(inout) :: state(:,:)

    !> The number of scalar variables to convert. If nScalars is not passed
    !! to this subroutine, all variables of argument state will be considered
    !! by this routine.
    integer, intent(in), optional :: nScalars
    !--------------------------------------------------------------------------
    integer :: iVar, iPoint, iVP, nPVars
    !--------------------------------------------------------------------------
    ! Information for the oversampling loop
    if(present(nScalars)) then
      nPVars = (poly_proj%maxPolyDegree+1)*nScalars
    else
      nPVars = (poly_proj%maxPolyDegree+1)*size(state,2)
    end if

    !$OMP DO
    do iVP = 1,nPVars
      iVar = (iVP-1)/(poly_proj%min_degree+1) + 1
      iPoint = iVP - (iVar-1)*(poly_proj%min_degree+1)
      state(iPoint,iVar) = modalCoeffs(iPoint,iVar)
    end do
    !$OMP END DO

  end subroutine ply_convertFromoversample_1d
  !****************************************************************************!

end module ply_oversample_module
