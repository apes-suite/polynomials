!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_fpt_ifpt_2D_singVar_lobattoNodes_test
  use env_module,               only: rk, fin_env
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_general_type, tem_start
  use ply_legFpt_module,        only: ply_legFpt_type, ply_init_legFPT
  use ply_legFpt_2D_module,     only: ply_legToPnt_2D,    &
    &                                 ply_pntToLeg_2D

  implicit none

  integer :: iPower
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = general             )

  res = 0.0_rk
  do iPower = 1,6
    call ply_check_legToPnt_2D(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if (res < 1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if 

  call fin_env()


contains


  subroutine ply_check_legToPnt_2D(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    integer :: maxPolyDegree, maxErr
    real(kind=rk), allocatable :: legCoeffs(:), legCoeffsIn(:)
    real(kind=rk), allocatable :: pntVal(:), legVal(:)
    type(ply_legFpt_type) :: fpt

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (per dim): ', maxPolyDegree+1
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (total): ',(maxPolyDegree+1)**2

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs((maxPolyDegree+1)**2))
    allocate(legCoeffsIn((maxPolyDegree+1)**2))
    legCoeffs(:) = real(1,rk)

    ! Init the FPT
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree,   &
      &                   fpt           = fpt,             &
      &                   nIndeps       = maxPolyDegree+1, &
      &                   lobattoPoints = .true.           )

    ! now transform to the Chebyshev nodes
    allocate(pntVal( (maxPolyDegree+1)**2))
    legCoeffsIn = legCoeffs ! Duplicate input vector to make sure that it is not modified in the trafo
    write(logUnit(10),*) 'Calculating FPT ...'
    !$OMP PARALLEL &
    !$OMP DEFAULT(shared)
    call ply_legToPnt_2D( fpt = fpt, legCoeffs = legCoeffsIn, pntVal = pntVal, &
      &                   lobattoPoints = .true. )
    !$OMP END PARALLEL
    write(logUnit(10),*) 'Finished'

    ! now transform back to Legendre coefficients
    allocate(legVal( (maxPolyDegree+1)**2 ))
    write(logUnit(10),*) 'Calculating inverse FPT ...'
    !$OMP PARALLEL &
    !$OMP DEFAULT(shared)
    call ply_pntToLeg_2D( fpt = fpt, pntVal = pntVal, legCoeffs = legVal, &
      &                   lobattoPoints = .true. )
    write(logUnit(10),*) 'Finished'
    !$OMP END PARALLEL

    !!do iDof = 1, (maxPolyDegree+1)**2
    !!  write(*,*) 'Leg coeff ', iDof, ' has error: ', legVal(iDof) - legCoeffs(iDof)
    !!end do

    ! Write out the coefficient with the largest absolute error
    write(*,*) 'For variable ', 1, &
             & ' Leg-Coeff ',maxloc(abs(legVal(:) - legCoeffs(:))),  &
             & ' has largest error of: ' ,maxval(abs(legVal(:) - legCoeffs(:)))
    maxErr = maxloc(abs(legVal(:) - legCoeffs(:)), 1)
    write(*,*) 'Ref. sol ', legCoeffs(maxErr), ' alg delivers: ', legVal(maxErr)

    res = maxval(abs(legVal(:) - legCoeffs(:)))

  end subroutine

end program ply_fpt_ifpt_2D_singVar_lobattoNodes_test
