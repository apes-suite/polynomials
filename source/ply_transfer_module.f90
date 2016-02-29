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
          if ((in_X > minord) .or. (in_Y > minord)) EXIT
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
          if ((out_X > minord) .or. (out_Y > minord)) EXIT
        end do

      else ispace_op

        ! Both input and output are P polynomials
        min_dofs = (minord*(minord+1))/2
        out_X = 1
        out_Y = 1
        do iDof=1,min_dofs
          out_pos = posOfModgCoeffPTens2D(out_X, out_Y, 0, outdegree)
          in_pos = posOfModgCoeffPTens2D(out_X, out_Y, 0, indegree)
          outdat(out_pos) = indat(in_pos)
          call nextModgCoeffPTens2D(out_X, out_Y, out_Z, minord-1)
        end do

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
            &                  .or. (in_Z > minord) ) EXIT
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
              &                  .and. (out_Z <= minord) ) EXIT
          end do
          if ((out_X > minord) .or. (out_Y > minord) &
            &                 .or. (out_Z > minord) ) EXIT
        end do

      else ispace_op

        ! Both input and output are P polynomials
        min_dofs = ( (minord+2)*((minord*(minord+1))/2) ) / 3
        out_X = 1
        out_Y = 1
        out_Z = 1
        do iDof=1,min_dofs
          out_pos = posOfModgCoeffPTens(out_X, out_Y, out_Z, outdegree)
          in_pos = posOfModgCoeffPTens(out_X, out_Y, out_Z, indegree)
          outdat(out_pos) = indat(in_pos)
          call nextModgCoeffPTens(out_X, out_Y, out_Z, minord-1)
        end do

      end if ispace_op

    end select ospace

  end subroutine ply_transfer_dofs_3d
  ! ************************************************************************ !


  ! ************************************************************************ !
  ! Some testing routines....                                                !
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Test for ply_transfer_dofs_1d.
  subroutine ply_test_transfer_1d(success, lu)
    ! -------------------------------------------------------------------- !
    !> Indicator whether the check was successful.
    !!
    !! True if all tests pass successfully, otherwise false.
    logical, intent(out) :: success

    !> A logunit to write messages to.
    integer, intent(in) :: lu
    ! -------------------------------------------------------------------- !

    real(kind=rk) :: smallpoly(3)
    real(kind=rk) :: largepoly(5)
    real(kind=rk) :: tmp(22)
    real(kind=rk) :: eps
    logical :: test_ok

    success = .true.

    eps = epsilon(smallpoly(1))

    ! Test small to large
    write(lu,*) '1D check transfer small to large:'
    smallpoly = [16.0_rk, 8.0_rk, 4.0_rk]
    call ply_transfer_dofs_1d( indat     = smallpoly, &
      &                        indegree  = 2,         &
      &                        outdat    = largepoly, &
      &                        outdegree = 4          )

    test_ok = ( all(abs(largepoly(:3) - smallpoly) < eps) &
      &         .and. all(abs(largepoly(4:)) < eps)       )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallpoly
      write(lu,*) '   large:', largepoly
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized)
    tmp = -5.0_rk
    write(lu,*) '1D check transfer to equal sized but section:'
    call ply_transfer_dofs_1d( indat     = smallpoly, &
      &                        indegree  = 2,         &
      &                        outdat    = tmp(6:8),  &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(6:8) - smallpoly) < eps)  &
      &         .and. all(abs(tmp(:5)+5.0_rk) < eps)  &
      &         .and. all(abs(tmp(10:)+5.0_rk) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallpoly
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small
    write(lu,*) '1D check transfer large to small:'
    largepoly = [16.0_rk, -8.0_rk, 4.0_rk, -2.0_rk, 1.0_rk]
    call ply_transfer_dofs_1d( indat     = largepoly, &
      &                        indegree  = 4,         &
      &                        outdat    = smallpoly, &
      &                        outdegree = 2          )

    test_ok = ( all(abs(largepoly(:3) - smallpoly) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largepoly
      write(lu,*) '   small:', smallpoly
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    if (success) then
      write(lu,*) 'Successfully passed checks for ply_transfer_dofs_1d.'
    else
      write(lu,*) 'FAILED checks for ply_transfer_dofs_1d!'
    end if

  end subroutine ply_test_transfer_1d
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Test for ply_transfer_dofs_2d.
  subroutine ply_test_transfer_2d(success, lu)
    ! -------------------------------------------------------------------- !
    !> Indicator whether the check was successful.
    !!
    !! True if all tests pass successfully, otherwise false.
    logical, intent(out) :: success

    !> A logunit to write messages to.
    integer, intent(in) :: lu
    ! -------------------------------------------------------------------- !

    real(kind=rk) :: smallq(9)
    real(kind=rk) :: largeq(25)
    real(kind=rk) :: smallp(6)
    real(kind=rk) :: largep(15)
    real(kind=rk) :: tmp(22)
    real(kind=rk) :: eps
    logical :: test_ok

    success = .true.

    eps = epsilon(smallq(1))

    ! Test small to large Q-Q
    write(lu,*) '2D Q check transfer small to large:'
    smallq = [ 16.0_rk,  8.0_rk, 4.0_rk, &
      &        27.0_rk,  9.0_rk, 3.0_rk, &
      &       125.0_rk, 25.0_rk, 5.0_rk  ]
    call ply_transfer_dofs_2d( indat     = smallq,  &
      &                        indegree  = 2,       &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = largeq,  &
      &                        outdegree = 4        )

    test_ok = ( all(abs(largeq(:3) - smallq(:3)) < eps) &
      &         .and. all(abs(largeq(4:5)) < eps)       )
    test_ok = test_ok &
      &       .and. ( all(abs(largeq(6:8) - smallq(4:6)) < eps) &
      &               .and. all(abs(largeq(9:10)) < eps)        )
    test_ok = test_ok &
      &       .and. ( all(abs(largeq(11:13) - smallq(7:9)) < eps) &
      &               .and. all(abs(largeq(14:)) < eps)           )
    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '   large:', largeq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test section (equally sized)
    tmp = -5.0_rk
    write(lu,*) '2D Q check transfer to equal sized but section:'
    call ply_transfer_dofs_2d( indat     = smallq,    &
      &                        indegree  = 2,         &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = tmp(6:14), &
      &                        outdegree = 2          )
    test_ok = ( all(abs(tmp(6:14) - smallq) < eps) &
      &         .and. all(abs(tmp(:5)+5.0_rk) < eps)      &
      &         .and. all(abs(tmp(15:)+5.0_rk) < eps)     )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   small:', smallq
      write(lu,*) '     tmp:', tmp
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    ! Test large to small
    write(lu,*) '2D Q check transfer large to small:'
    largeq = [ 16.0_rk, -8.0_rk,  4.0_rk, -2.0_rk, 1.0_rk, &
      &        -8.0_rk,  4.0_rk, -2.0_rk,  1.0_rk, 0.5_rk, &
      &         4.0_rk, -2.0_rk,  1.0_rk, -0.5_rk, 0.3_rk, &
      &        -2.0_rk,  1.0_rk, -0.5_rk,  0.3_rk, 0.1_rk, &
      &         1.0_rk, -0.5_rk,  0.3_rk, -0.1_rk, 0.0_rk  ]
    call ply_transfer_dofs_2d( indat     = largeq, &
      &                        indegree  = 4,      &
      &                        inspace   = Q_Space, &
      &                        outspace  = Q_Space, &
      &                        outdat    = smallq, &
      &                        outdegree = 2       )

    test_ok = ( all(abs(largeq(:3) - smallq(:3)) < eps) )
    test_ok = test_ok &
      &       .and. ( all(abs(largeq(6:8) - smallq(4:6)) < eps) )
    test_ok = test_ok &
      &       .and. ( all(abs(largeq(11:13) - smallq(7:9)) < eps) )

    if (test_ok) then
      write(lu,*) '  OK'
    else
      write(lu,*) '  FAILED!'
      write(lu,*) '   large:', largeq
      write(lu,*) '   small:', smallq
    end if

    success = (success .and. test_ok)

    write(lu,*) '----------------------------------------------------------'

    if (success) then
      write(lu,*) 'Successfully passed checks for ply_transfer_dofs_1d.'
    else
      write(lu,*) 'FAILED checks for ply_transfer_dofs_1d!'
    end if

  end subroutine ply_test_transfer_2d
  ! ************************************************************************ !


end module ply_transfer_module
