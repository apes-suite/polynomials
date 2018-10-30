?? include "ply_dof_module.inc"
!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_fpt_2D_singVar_lobattoNodes_test
  use env_module,               only: rk, fin_env
  use tem_param_module,         only: PI
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_general_type, tem_start
  use ply_legFpt_module,        only: ply_legFpt_type, ply_init_legFPT, &
    &                                 ply_legFpt_bu_type
  use ply_legFpt_2D_module,     only: ply_legToPnt_2D
  use ply_modg_basis_module,    only: evalLegendreTensPoly
  use ply_dof_module,           only: Q_space

  implicit none

  integer :: iPower
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = general      )

  res = 0.0_rk
  do iPower = 1, 4
    call ply_check_legToPnt_2D(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if

  call fin_env()

contains

  subroutine ply_check_legToPnt_2D(power, res)
    integer, intent(in) :: power
    real(kind=rk) :: res
    integer :: maxPolyDegree, iPoint, iPointX, iPointY, iDof
    integer :: pointIndex, funcIndex, iPolyX, iPolyY
    real(kind=rk), allocatable :: legCoeffs(:)
    real(kind=rk), allocatable :: pntVal(:), refVal(:)
    real(kind=rk), allocatable :: chebPnt(:,:), chebPnt1D(:)
    real(kind=rk), allocatable :: legValChebPnt(:,:)
    real(kind=rk) :: rfac
    type(ply_legFpt_type) :: fpt
    type(ply_legFpt_bu_type) :: bu
    integer, allocatable :: rand_seed(:)
    integer :: nSeeds

    ! Init the random number generator
    call random_seed(size=nSeeds)
    allocate(rand_seed(nSeeds))
    rand_seed = 0
    rand_seed(1) = 8345
    call random_seed(put=rand_seed)

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (per dim): ', maxPolyDegree+1
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (total): ',(maxPolyDegree+1)**2

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs((maxPolyDegree+1)**2))
    do iDof = 1, (maxPolyDegree+1)**2
      call random_number(rfac)
      legCoeffs(iDof) = real(1, rk) * rfac
    end do

    ! Create the Chebyshev nodes on the interval [-1,+1]
    write(logUnit(10),*) 'Creating Chebyshev nodes for ref result ...'
    allocate(chebPnt1D(maxPolyDegree+1))
    do iPoint = 1, maxPolyDegree+1
      chebPnt1D(iPoint) = cos((iPoint-1.0_rk)*PI/maxPolyDegree);
    end do
    ! Now, create the 2D Chebyshev nodes in [-1,+1]^3
    allocate(chebPnt( (maxPolyDegree+1)**2,3 ))
    do iPointX = 1, maxPolyDegree+1
      do iPointY = 1, maxPolyDegree+1
        pointIndex = 1 + (iPointX-1) + (iPointY-1)*(maxPolyDegree+1)
        chebPnt(pointIndex,1) = chebPnt1D(iPointX)
        chebPnt(pointIndex,2) = chebPnt1D(iPointY)
        chebPnt(pointIndex,3) = 1.0_rk
      end do
    end do

    ! define the reference results for the point values (Chebyshev nodes)
    !allocate( legValChebPnt((maxPolyDegree+1)**3,(maxPolyDegree+1)**3) )
    !legValChebPnt(:,:) = legendre_1D(chebPnt, maxPolyDegree)
    call evalLegendreTensPoly( coords = chebPnt , nCoords = (maxPolyDegree+1)**2, &
                             & maxPolyDegree = maxPolyDegree, &
                             & basisType = Q_space, &
                             & polyVal = legValChebPnt )
    allocate( refVal((maxPolyDegree+1)**2) )
    refVal(:) = 0.0_rk
    write(logUnit(10),*) 'Calculating reference results ...'
    do iPolyX = 1, maxPolyDegree+1
      do iPolyY = 1, maxPolyDegree+1
?? copy :: posOfModgCoeffQTens(iPolyX, iPolyY,1, maxPolyDegree, funcIndex)
        refVal(:) = refVal(:) + &
                 & legValChebPnt(funcIndex,:) * legCoeffs(funcIndex)
      end do
    end do
    write(logUnit(10),*) 'Finished'

    ! Init the FPT
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree,   &
      &                   nIndeps       = maxpolydegree+1, &
      &                   fpt           = fpt,             &
      &                   lobattoPoints = .true.,          &
      &                   bu            = bu               )

    ! now transform to the Chebyshev nodes
    allocate(pntVal( (maxPolyDegree+1)**2 ))
    write(logUnit(10),*) 'Calculating FPT ...'
    call ply_legToPnt_2D( fpt = fpt, legCoeffs = legCoeffs, pntVal = pntVal, &
        &                 bu = bu )
    write(logUnit(10),*) 'Finished'

    !!do iPoint = 1, (maxPolyDegree+1)**3
    !!  write(*,*) 'Point: ', iPoint, &
    !!           & ' FPT: ', pntVal(iPoint,3), &
    !!           & ' Ref.: ', refVal(iPoint,3), &
    !!           & ' error: ', pntVal(iPoint,3)-refVal(iPoint,3)
    !!end do

    ! Write out the point with the largest absolute error
    write(*,*) 'for variable ', 1, &
           & ' Cheb-Point ',maxloc(abs(pntVal(:) - refVal(:))),  &
           & ' has largest error of: ' ,maxval(abs(pntVal(:) - refVal(:)))

    res = maxval(abs(pntVal(:) - refVal(:)))

  end subroutine

end program ply_fpt_2D_singVar_lobattoNodes_test
