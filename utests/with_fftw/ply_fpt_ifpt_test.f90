!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_fpt_ifpt_test
  use env_module,               only: rk, fin_env
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use ply_legFpt_module,        only: ply_init_legFpt, ply_legFpt_type, &
    &                                 ply_legToPnt, ply_pntToLeg
  use ply_modg_basis_module,    only: legendre_1D
  use tem_general_module,       only: tem_general_type, tem_start

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
    call check_fwd_bwd(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if 

  call fin_env()

contains

  subroutine check_fwd_bwd(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    integer :: maxPolyDegree
    real(kind=rk), allocatable :: legCoeffs(:), pntVal(:), legVal(:)
    type(ply_legFpt_type) :: fpt


    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    write(logUnit(10),*) '------- Number of Legendre coefficients: ', maxPolyDegree+1
  
    ! Create the Legendre expansion coefficients
    allocate(legCoeffs(1:maxPolyDegree+1)) 
    allocate(legVal(1:maxPolyDegree+1)) 
    legCoeffs(:) = 1.0_rk
    legVal = legCoeffs

    ! Init the FPT 
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree, &
      &                   nIndeps       = 1,             &
      &                   fpt           = fpt            )

    ! now transform to the Chebyshev nodes
    allocate(pntVal(1:maxPolyDegree+1)) 
    write(logUnit(10),*) 'Calculating FPT ...'
    call ply_legToPnt( fpt = fpt, legCoeffs = legVal, pntVal = pntVal, &
      &                lobattoPoints = .false. ) 
    write(logUnit(10),*) 'Finished'

    ! now transform to the Legendre coefficients
    write(logUnit(10),*) 'Calculating inverse FPT ...'
    call ply_pntToLeg( fpt = fpt, pntVal = pntVal, legCoeffs = legVal, & 
      &                lobattoPoints = .false. ) 
    write(logUnit(10),*) 'Finished'

    ! Write out the polynomial coefficient with the largest absolute error
    write(*,*) 'Leg-Coefficient ',maxloc(abs(legVal(:) - legCoeffs(:))),  &
              & ' has largest error of: ' ,maxval(abs(legVal(:) - legCoeffs(:)))

    res = maxval(abs(legVal(:) - legCoeffs(:)))

  end subroutine check_fwd_bwd

end program ply_fpt_ifpt_test
