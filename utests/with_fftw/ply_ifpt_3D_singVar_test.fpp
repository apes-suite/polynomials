?? include "ply_dof_module.inc"
!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_ifpt_3D_singVar_test
  use env_module,               only: rk, fin_env
  use tem_param_module,         only: PI
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use tem_general_module,       only: tem_general_type, tem_start
  use ply_legFpt_module,        only: ply_legFpt_type, ply_init_legFPT
  use ply_legFpt_3D_module,     only: ply_pntToLeg_3D
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
  do iPower = 0, 4
    call ply_check_legToPnt_3D(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
  call fin_env()

contains

  subroutine ply_check_legToPnt_3D(power, res)
    integer, intent(in) :: power
    real(kind=rk) :: res
    integer :: maxPolyDegree, iPoint, iPointX, iPointY, iPointZ, iDof
    integer :: pointIndex, funcIndex, iPolyX, iPolyY, iPolyZ
    real(kind=rk), allocatable :: legCoeffs(:), legCoeffsRef(:)
    real(kind=rk), allocatable :: pntVal(:)
    real(kind=rk), allocatable :: chebPnt1D(:)
    real(kind=rk), allocatable :: legValChebPnt(:,:)
    real(kind=rk) :: rfac
    type(ply_legFpt_type) :: fpt
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
    write(logUnit(10),*) '------------------------------------' &
      & // ' Number of Legendre coefficients (per dim): ', maxPolyDegree+1
    write(logUnit(10),*) '------------------------------------' &
      & // ' Number of Legendre coefficients (total): ',(maxPolyDegree+1)**3

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs((maxPolyDegree+1)**3))
    allocate(legCoeffsRef((maxPolyDegree+1)**3))
    do iDof = 1, (maxPolyDegree+1)**3
      call random_number(rfac)
      legCoeffsRef(iDof) = real(1, rk) * rfac
    end do

    ! Create the Chebyshev nodes on the interval [-1,+1]
    write(logUnit(10),*) 'Creating Chebyshev nodes for ref result ...'
    allocate(chebPnt1D(maxPolyDegree+1))
    do iPoint = 1, maxPolyDegree+1
      chebPnt1D(iPoint) = (-1.0_rk) * &
                & cos(PI/(maxPolyDegree+1)*((iPoint-1.0_rk)+1.0_rk/2.0_rk))
    end do

    ! define the reference results for the point values (Chebyshev nodes)
    allocate( legValChebPnt((maxPolyDegree+1),(maxPolyDegree+1)) )
    legValChebPnt(:,:) = legendre_1D(chebPnt1D, maxPolyDegree)
    allocate(pntVal( (maxPolyDegree+1)**3 ))
    pntVal(:) = 0.0_rk
    write(logUnit(10),*) 'Calculating reference results ...'
    do iPolyX = 1, maxPolyDegree+1
      do iPolyY = 1, maxPolyDegree+1
        do iPolyZ = 1, maxPolyDegree+1
?? copy :: posOfModgCoeffQTens(iPolyX, iPolyY, iPolyZ, maxPolyDegree, funcIndex)
          do iPointX = 1, maxPolyDegree+1
            do iPointY = 1, maxPolyDegree+1
              do iPointZ = 1, maxPolyDegree+1
                pointIndex = 1 + (iPointX-1) + (iPointY-1)*(maxPolyDegree+1) &
                  &        + (iPointZ-1)*((maxPolyDegree+1)**2)
                pntVal(pointIndex) = pntVal(pointIndex)               &
                  &                + legValChebPnt(iPolyX, iPointX)   &
                  &                  * legValChebPnt(iPolyY, iPointY) &
                  &                  * legValChebPnt(iPolyZ, iPointZ) &
                  &                  * legCoeffsRef(funcIndex)
              end do
            end do
          end do
        end do
      end do
    end do
    write(logUnit(10),*) 'Finished'

    ! Init the FPT
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree,        &
      &                   nIndeps       = (maxPolyDegree+1)**2, &
      &                   fpt           = fpt                   )

    ! now transform to the Chebyshev nodes
    write(logUnit(10),*) 'Calculating FPT ...'
    call ply_pntToLeg_3D( fpt = fpt, pntVal = pntVal, legCoeffs = legCoeffs )
    write(logUnit(10),*) 'Finished'

    !!do iPoint = 1, (maxPolyDegree+1)**3
    !!  write(*,*) 'Point: ', iPoint, &
    !!           & ' FPT: ', pntVal(iPoint,3), &
    !!           & ' Ref.: ', refVal(iPoint,3), &
    !!           & ' error: ', pntVal(iPoint,3)-refVal(iPoint,3)
    !!end do

    ! Write out the point with the largest absolute error
    write(*,*) 'for variable ', 1, &
             & ' Cheb-Point ',maxloc(abs(legCoeffs(:) - legCoeffsRef(:))),  &
             & ' has largest error of: ' ,maxval(abs(legCoeffs(:) - legCoeffsRef(:)))

    res = maxval(abs(legCoeffs(:) - legCoeffsRef(:)))

  end subroutine

end program ply_ifpt_3D_singVar_test
