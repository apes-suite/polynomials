!> This module provides routines to transfer degrees of freedom from one
!! polynomial representation to another.
!!
!! These routines allow the copying of data from one order and polynomial space
!! into another.
module ply_transfer_module
  use env_module, only: rk

  use ply_dof_module, only: P_Space, Q_space, posOfModgCoeffPTens2D, &
    &                       posOfModgCoeffPTens, nextModgCoeffPTens2D, &
    &                       nextModgCoeffPTens


  implicit none

  private

  public :: ply_transfer_dofs_1D
  public :: ply_transfer_dofs_2D
  public :: ply_transfer_dofs_3D
  public :: ply_transfer_P_dim


contains


  ! ************************************************************************ !
  !> Transfer of degrees of freedom from one polynomial to another in 1D.
  !!
  !! If the indat is larger than outdat, the higher modes are truncated.
  !! If outdat is larger, higher modes are padded with zeros.
  subroutine ply_transfer_dofs_1D( indat, indegree, outdat, outdegree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Degree of the input polynomial. There are indegree+1 modes expected
    !! in indat.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Degree of the output polynomial. There are outdegree+1 modes expected
    !! in outdat.
    integer, intent(in) :: outdegree
    ! -------------------------------------------------------------------- !
    integer :: minOrd
    ! -------------------------------------------------------------------- !

    minord = min(outdegree+1, indegree+1)

    outdat(:minord) = indat(:minord)
    outdat(minord+1:) = 0.0_rk

  end subroutine ply_transfer_dofs_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transfer of degrees of freedom from one polynomial to another in 2D.
  !!
  !! If the indat is larger than outdat, the higher modes are truncated.
  !! If outdat is larger, higher modes are padded with zeros.
  !!
  !! When different multidimensional polynomial layouts are used, the modes
  !! are copied to to the corresponding locations.
  subroutine ply_transfer_dofs_2D( indat, inspace, indegree,   &
    &                              outdat, outspace, outdegree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Multi-dimensional polynomial layout of the input data.
    !!
    !! Has to be either [[Q_Space]] or [[P_Space]].
    integer, intent(in) :: inspace

    !> Maximal polynomial degree in the input data.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Multi-dimensional polynomial layout of the output data.
    !!
    !! Has to be either [[Q_Space]] or [[P_Space]].
    integer, intent(in) :: outspace

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: outdegree
    ! -------------------------------------------------------------------- !
    integer :: outdofs
    integer :: indofs
    integer :: min_dofs
    integer :: iStep, iDof
    integer :: out_X, out_Y, out_Z
    integer :: in_X, in_Y, in_Z
    integer :: out_pos, in_pos
    integer :: out_off, out_zoff
    integer :: in_off, in_zoff
    integer :: minOrd
    ! -------------------------------------------------------------------- !

    minord = min(outdegree+1, indegree+1)

    select case(inspace)
    case (Q_Space)
      indofs = (indegree+1)**2

    case (P_Space)
      indofs = ((indegree+1)*(indegree+2)) / 2

    end select

    outdat = 0.0_rk

    ospace: select case(outspace)
    case (Q_Space) ospace
      outdofs = (outdegree+1)**2

      ispace_oq: if (inspace == Q_Space) then

        ! Both, output and input are Q Polynomials
        do out_Y=0,minord-1
          out_off = out_Y*(outdegree+1)
          in_off = out_Y*(indegree+1)
          outdat(out_off+1:out_off+minord) = indat(in_off+1:in_off+minord)
        end do

      else ispace_oq

        ! Output is Q, but input is P
        in_X = 1
        in_Y = 1
        do iDof=1,indofs
          in_pos = posOfModgCoeffPTens2D(in_X, in_Y, 0, indegree)
          out_pos = in_X + (outdegree+1)*(in_Y-1)
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the minord range
          do istep=iDof,indofs
            call nextModgCoeffPTens2D(in_X, in_Y, in_Z, indegree)
            if ((in_X <= minord) .and. (in_Y <= minord)) EXIT
          end do
          if ( (in_X > minord) .or. (in_Y > minord) &
             & .or. (in_X+in_Y-2 > indegree)        ) EXIT
        end do

      end if ispace_oq


    case (P_Space) ospace
      outdofs = ((outdegree+1)*(outdegree+2)) / 2

      ispace_op: if (inspace == Q_Space) then

        ! Output is P, input is Q
        out_X = 1
        out_Y = 1
        do iDof=1,outdofs
          out_pos = posOfModgCoeffPTens2D(out_X, out_Y, 0, outdegree)
          in_pos = out_X + (indegree+1)*(out_Y-1)
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the target range
          do istep=iDof,outdofs
            call nextModgCoeffPTens2D(out_X, out_Y, out_Z, outdegree)
            if ((out_X <= minord) .and. (out_Y <= minord)) EXIT
          end do
          if ( (out_X > minord) .or. (out_Y > minord) &
            &  .or. (out_X+out_Y-2 > outdegree)       ) EXIT
        end do

      else ispace_op

        ! Both input and output are P polynomials
        min_dofs = (minord*(minord+1))/2
        outdat(:min_dofs) = indat(:min_dofs)

      end if ispace_op

    end select ospace

  end subroutine ply_transfer_dofs_2d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transfer of degrees of freedom from one polynomial to another in 3D.
  !!
  !! If the indat is larger than outdat, the higher modes are truncated.
  !! If outdat is larger, higher modes are padded with zeros.
  !!
  !! When different multidimensional polynomial layouts are used, the modes
  !! are copied to to the corresponding locations.
  subroutine ply_transfer_dofs_3D( indat, inspace, indegree,   &
    &                              outdat, outspace, outdegree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Multi-dimensional polynomial layout of the input data.
    !!
    !! Has to be either [[Q_Space]] or [[P_Space]].
    integer, intent(in) :: inspace

    !> Maximal polynomial degree in the input data.
    integer, intent(in) :: indegree

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Multi-dimensional polynomial layout of the output data.
    !!
    !! Has to be either [[Q_Space]] or [[P_Space]].
    integer, intent(in) :: outspace

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: outdegree
    ! -------------------------------------------------------------------- !
    integer :: outdofs
    integer :: indofs
    integer :: min_dofs
    integer :: iStep, iDof
    integer :: out_X, out_Y, out_Z
    integer :: in_X, in_Y, in_Z
    integer :: out_pos, in_pos
    integer :: out_off, out_zoff
    integer :: in_off, in_zoff
    integer :: minOrd
    ! -------------------------------------------------------------------- !

    minord = min(outdegree+1, indegree+1)

    select case(inspace)
    case (Q_Space)
      indofs = (indegree+1)**3

    case (P_Space)
      indofs = ((indegree+1)*(indegree+2)*(indegree+3))/6

    end select

    outdat = 0.0_rk

    ospace: select case(outspace)
    case (Q_Space) ospace
      outdofs = (outdegree+1)**3
      ispace_oq: if (inspace == Q_Space) then

        ! Both, output and input are Q Polynomials
        do out_Z=0,minord-1
          out_zoff = out_Z*(outdegree+1)**2
          in_zoff = out_Z*(indegree+1)**2
          do out_Y=0,minord-1
            out_off = out_Y*(outdegree+1) + out_zoff
            in_off = out_Y*(indegree+1) + in_zoff
            outdat(out_off+1:out_off+minord) &
              &  = indat(in_off+1:in_off+minord)
          end do
        end do

      else ispace_oq

        ! Output is Q, but input is P
        in_X = 1
        in_Y = 1
        in_Z = 1
        do iDof=1,indofs
          in_pos = posOfModgCoeffPTens(in_X, in_Y, in_Z, indegree)
          out_pos = in_X + (outdegree+1)*( (in_Y-1) &
            &                             + (outdegree+1)*(in_Z-1))
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the target range
          do istep=iDof,indofs
            call nextModgCoeffPTens(in_X, in_Y, in_Z, indegree)
            if ( (in_X <= minord) .and. (in_Y <= minord) &
              &                   .and. (in_Z <= minord) ) EXIT
          end do
          if ( (in_X > minord) .or. (in_Y > minord) &
            &                  .or. (in_Z > minord) &
            &  .or. (in_X+in_Y+in_Z-3 > indegree)   ) EXIT
        end do

      end if ispace_oq

    case (P_Space) ospace
      outdofs = ((outdegree+1)*(outdegree+2)*(outdegree+3))/6

      ispace_op: if (inspace == Q_Space) then
        ! Output is P, input is Q
        out_X = 1
        out_Y = 1
        out_Z = 1
        do iDof=1,outdofs
          out_pos = posOfModgCoeffPTens(out_X, out_Y, out_Z, outdegree)
          in_pos = out_X + (indegree+1)*( (out_Y-1) &
            &                                + (indegree+1)*(out_Z-1))
          outdat(out_pos) = indat(in_pos)
          ! Ensure, that next iteration is in the target range
          do istep=iDof,outdofs
            call nextModgCoeffPTens(out_X, out_Y, out_Z, outdegree)
            if ((out_X <= minord) .and. (out_Y <= minord) &
              &                   .and. (out_Z <= minord) ) EXIT
          end do
          if ( (out_X > minord) .or. (out_Y > minord) &
            &                   .or. (out_Z > minord) &
            &  .or. (out_X+out_Y+out_Z-3 > outdegree) ) EXIT
        end do

      else ispace_op

        ! Both input and output are P polynomials
        min_dofs = ( (minord+2)*((minord*(minord+1))/2) ) / 3
        outdat(:min_dofs) = indat(:min_dofs)

      end if ispace_op

    end select ospace

  end subroutine ply_transfer_dofs_3d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Transfer the polynomial in P representation from on dimension to
  !! another one.
  subroutine ply_transfer_P_dim( indat, indim, outdat, outdim, degree )
    ! -------------------------------------------------------------------- !
    !> Input data to transfer to output data.
    real(kind=rk), intent(in) :: indat(:)

    !> Dimension of the input polynomial.
    integer, intent(in) :: indim

    !> Output data to fill with input data.
    real(kind=rk), intent(out) :: outdat(:)

    !> Dimension of the output polynomial.
    integer, intent(in) :: outdim

    !> Maximal polynomial degree in the output data.
    integer, intent(in) :: degree
    ! -------------------------------------------------------------------- !
    integer :: nInDofs, nOutDofs
    integer :: iMode
    integer :: iPos
    integer :: iX, iY, iZ
    ! -------------------------------------------------------------------- !

    in_d: select case(indim)
    case(1) in_d

      i1_out_d: select case(outdim)
      case(1) i1_out_d
        outdat = indat

      case(2) i1_out_d
        outdat = 0.0_rk
        do iMode=1,degree+1
          iPos = posOfModgCoeffPTens2D( ansFuncX  = iMode, &
            &                           ansFuncY  = 1,     &
            &                           ansFuncZ  = 1,     &
            &                           maxDegree = degree )
          outdat(iPos) = indat(iMode)
        end do

      case(3) i1_out_d
        outdat = 0.0_rk
        do iMode=1,degree+1
          iPos = posOfModgCoeffPTens( ansFuncX  = iMode, &
            &                         ansFuncY  = 1,     &
            &                         ansFuncZ  = 1,     &
            &                         maxDegree = degree )
          outdat(iPos) = indat(iMode)
        end do

      end select i1_out_d


    case(2) in_d

      i2_out_d: select case(outdim)
      case(1) i2_out_d
        do iMode=1,degree+1
          iPos = posOfModgCoeffPTens2D( ansFuncX  = iMode, &
            &                           ansFuncY  = 1,     &
            &                           ansFuncZ  = 1,     &
            &                           maxDegree = degree )
          outdat(iMode) = indat(iPos)
        end do

      case(2) i2_out_d
        outdat = indat

      case(3) i2_out_d
        outdat = 0.0_rk
        nInDofs = ((degree+1)*(degree+2))/2
        iX = 1
        iY = 1
        iZ = 1
        do iMode=1,nInDofs
          iPos = posOfModgCoeffPTens( ansFuncX  = iX,    &
            &                         ansFuncY  = iY,    &
            &                         ansFuncZ  = iZ,    &
            &                         maxDegree = degree )
          outdat(iPos) = indat(iMode)
          call nextModgCoeffPTens2D( ansFuncX  = iX,    &
            &                        ansFuncY  = iY,    &
            &                        ansFuncZ  = iZ,    &
            &                        maxDegree = degree )
        end do

      end select i2_out_d


    case(3) in_d

      i3_out_d: select case(outdim)
      case(1) i3_out_d
        do iMode=1,degree+1
          iPos = posOfModgCoeffPTens( ansFuncX  = iMode, &
            &                         ansFuncY  = 1,     &
            &                         ansFuncZ  = 1,     &
            &                         maxDegree = degree )
          outdat(iMode) = indat(iPos)
        end do

      case(2) i3_out_d
        nOutDofs = ((degree+1)*(degree+2))/2
        iX = 1
        iY = 1
        iZ = 1
        do iMode=1,nOutDofs
          iPos = posOfModgCoeffPTens( ansFuncX  = iX,    &
            &                         ansFuncY  = iY,    &
            &                         ansFuncZ  = iZ,    &
            &                         maxDegree = degree )
          outdat(iMode) = indat(iPos)
          call nextModgCoeffPTens2D( ansFuncX  = iX,    &
            &                        ansFuncY  = iY,    &
            &                        ansFuncZ  = iZ,    &
            &                        maxDegree = degree )
        end do

      case(3) i3_out_d
        outdat = indat

      end select i3_out_d


    end select in_d

  end subroutine ply_transfer_P_dim
  ! ************************************************************************ !


end module ply_transfer_module
