!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_ifpt_test
  use env_module,               only: rk, fin_env
  use tem_param_module,         only: PI
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_general_type, tem_start
  use ply_legFpt_module,        only: ply_init_legFpt, &
    &                                 ply_legFpt_type, &
    &                                 ply_pntToLeg, ply_legFpt_bu_type
  use ply_modg_basis_module,    only: legendre_1D

  implicit none

  integer :: iPower
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = general             )

  res = 0.0_rk
  do iPower = 1,4
    call ply_check_pntToLeg(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
  call fin_env()

contains

  subroutine ply_check_pntToLeg(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    integer :: maxPolyDegree, iPoint, iPoly
    real(kind=rk), allocatable :: legCoeffs(:)
    real(kind=rk), allocatable :: pntVal(:), legVal(:)
    real(kind=rk), allocatable :: chebPnt(:)
    real(kind=rk), allocatable :: legValChebPnt(:,:)
    type(ply_legFpt_type) :: fpt
    type(ply_legFpt_bu_type) :: bu

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    write(logUnit(10),*) '------- Number of Legendre coefficients: ', maxPolyDegree+1

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs(1:maxPolyDegree+1))
    legCoeffs(:) = 1.0_rk

    ! Create the Chebyshev nodes on the interval [-1,+1]
    allocate(chebPnt(maxPolyDegree+1))
    do iPoint = 1, maxPolyDegree+1
      chebPnt(iPoint) = (-1.0_rk) * cos(PI/(maxPolyDegree+1)*((iPoint-1.0_rk)+1.0_rk/2.0_rk))
      !write(*,*) 'Cehbyshev point', iPoint, ' is at: ', chebPnt(iPoint)
    end do

    ! define the point values (Chebyshev nodes)
    allocate( legValChebPnt(maxPolyDegree+1,maxPolyDegree+1) )
    legValChebPnt(:,:) = legendre_1D(chebPnt, maxPolyDegree)
    allocate(pntVal(maxPolyDegree+1))
    pntVal(:) = 0.0_rk
    write(logUnit(10),*) 'Calculating point values (input) ...'
    do iPoly = 1, maxPolyDegree+1
      pntVal(:) = pntVal(:) + legValChebPnt(iPoly,:) * legCoeffs(iPoly)
    end do
    write(logUnit(10),*) 'Finished'

    ! Init the FPT
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree, &
      &                   nIndeps       = 1,             &
      &                   fpt           = fpt,           &
      &                   bu            = bu             )

    ! now transform to the Legendre coefficients
    allocate(legVal(1:maxPolyDegree+1))
    write(logUnit(10),*) 'Calculating inverse FPT ...'
    call ply_pntToLeg( fpt = fpt, pntVal = pntVal, legCoeffs = legVal, &
        &              nIndeps=1, bu = bu )
    write(logUnit(10),*) 'Finished'

    !!do iPoly = 1, maxPolyDegree+1
    !!  write(*,*) 'Poly degree: ', iPoly, &
    !!           & ' iFPT: ', legVal(iPoly), &
    !!           & ' Ref.: ', legCoeffs(iPoly), &
    !!           & ' error: ', legVal(iPoly)-legCoeffs(iPoly)
    !!end do

    ! Write out the polynomial coefficient with the largest absolute error
    write(*,*) 'Leg-Coefficient ',maxloc(abs(legVal(:) - legCoeffs(:))),  &
              & ' has largest error of: ' ,maxval(abs(legVal(:) - legCoeffs(:)))

    res = maxval(abs(legVal(:) - legCoeffs(:)))

  end subroutine

end program ply_ifpt_test
