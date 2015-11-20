!> Unit test to check functionallity of fast polynomial transformations
!! to Lobatto-Chebyshev-Nodes.
!! \author{Jens Zudrop}
program ply_fpt_lobattoNodes_test
  use env_module,               only: rk, fin_env
  use tem_param_module,         only: PI
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use ply_legFpt_module,        only: ply_init_legFpt, ply_legFpt_type, &
    &                                 ply_legToPnt
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
    call ply_check_legToPnt(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if 
  call fin_env()

contains

  subroutine ply_check_legToPnt(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    integer :: maxPolyDegree, iPoint, iPoly
    real(kind=rk), allocatable :: legCoeffs(:)
    real(kind=rk), allocatable :: pntVal(:), refVal(:)
    real(kind=rk), allocatable :: chebPnt(:)
    real(kind=rk), allocatable :: legValChebPnt(:,:)
    type(ply_legFpt_type) :: fpt
  
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
      chebPnt(iPoint) = cos((iPoint-1.0_rk)*PI/maxPolyDegree);
      !write(*,*) 'Lobatto-Chebyshev-Point', iPoint, chebPnt(iPoint)
    end do
  
    ! define the reference results for the point values (Chebyshev nodes)
    allocate( legValChebPnt(maxPolyDegree+1,maxPolyDegree+1) )
    legValChebPnt(:,:) = legendre_1D(chebPnt, maxPolyDegree)
    allocate(refVal(maxPolyDegree+1))
    refVal(:) = 0.0_rk
    write(logUnit(10),*) 'Calculating reference results ...'
    do iPoly = 1, maxPolyDegree+1
      refVal(:) = refVal(:) + legValChebPnt(iPoly,:) * legCoeffs(iPoly)
    end do
    write(logUnit(10),*) 'Finished'
  
    ! Init the FPT 
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree, fpt = fpt, &
                        & lobattoPoints = .true. ) 
  
    ! now transform to the Chebyshev nodes
    allocate(pntVal(1:maxPolyDegree+1)) 
    write(logUnit(10),*) 'Calculating FPT ...'
    call ply_legToPnt( fpt = fpt, legCoeffs = legCoeffs, pntVal = pntVal , &
      &                lobattoPoints = .true. )
    write(logUnit(10),*) 'Finished'
  
    !!do iPoint = 1, maxPolyDegree+1
    !!  write(*,*) 'Point: ', chebPnt(iPoint), &
    !!           & ' FPT: ', pntVal(iPoint), & 
    !!           & ' Ref.: ', refVal(iPoint), &
    !!           & ' error: ', pntVal(iPoint)-refVal(iPoint)
    !!end do

    ! Write out the point with the largest absolute error
    write(*,*) 'Lobatto-Cheb-Point ',maxloc(abs(pntVal(:) - refVal(:))),  &
             & ' has largest error of: ' ,maxval(abs(pntVal(:) - refVal(:)))

    res = maxval(abs(pntVal(:) - refVal(:)))

  end subroutine

end program ply_fpt_lobattoNodes_test
