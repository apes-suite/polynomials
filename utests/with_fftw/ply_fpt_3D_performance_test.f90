!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_fpt_3D_performance_test
  use mpi, only: mpi_wtime
  use env_module,               only: rk, fin_env
  use tem_logging_module,       only: logUnit, tem_logging_init_primary
  use tem_general_module,       only: tem_general_type, tem_start
  use ply_legFpt_module,        only: ply_legFpt_type, ply_init_legFPT
  use ply_legFpt_3D_module,     only: ply_legToPnt_3D, &
    &                                 ply_pntToLeg_3D

  implicit none

  integer :: iPower
  integer, parameter :: maxpower = 8
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'FPT 3D Performance Test', &
    &            version  = '1',                       &
    &            general  = general                    )
  call tem_logging_init_primary( level = 1,               &
    &                            rank = general%proc%rank )

  res = 0.0_rk
  do iPower = 1,maxpower
    call ply_check_legToPnt_3D(iPower, newRes)
    if (newRes.gt.res) then
      res = newRes
    end if
  end do

  write(*,*) 'Maximal deviation:', res
  if (res <  1.e-08) then
    write(logUnit(1),*) 'PASSED'
  end if 

  call fin_env()

contains

  subroutine ply_check_legToPnt_3D(power,res)
    integer, intent(in) :: power
    real(kind=rk) :: res
    integer :: maxPolyDegree, iVar, nVars
    real(kind=rk), allocatable :: legCoeffs(:,:), legCoeffsIn(:,:)
    real(kind=rk), allocatable :: pntVal(:,:), legVal(:,:)
    type(ply_legFpt_type) :: fpt
    real(kind=rk) :: starttime, stoptime

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree =  2**power-1  ! maxPolyDegree+1 has to be a power of 2
    nVars = 3
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (per dim): ', maxPolyDegree+1
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (total): ',(maxPolyDegree+1)**3
  
    ! Create the Legendre expansion coefficients
    allocate(legCoeffs((maxPolyDegree+1)**3,nVars)) 
    allocate(legCoeffsIn((maxPolyDegree+1)**3,nVars)) 
    do iVar = 1, nVars
      legCoeffs(:,iVar) = real(iVar, rk)
    end do
  
    ! Init the FPT 
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree,        &
      &                   nIndeps       = (maxpolydegree+1)**2, &
      &                   fpt           = fpt                   )
  
    ! now transform to the Chebyshev nodes
    allocate(pntVal( (maxPolyDegree+1)**3, nVars )) 
    legCoeffsIn = legCoeffs
    starttime = MPI_Wtime()
    call ply_legToPnt_3D( fpt       = fpt,         &
      &                   legCoeffs = legCoeffsIn, &
      &                   pntVal    = pntVal,      &
      &                   nVars     = nVars        ) 
    stoptime = MPI_Wtime()
    write(*,*) 'Time for degree ', maxpolydegree, ' trafo:   ', stoptime - starttime

    ! now transform back to Legendre coefficients
    allocate(legVal( (maxPolyDegree+1)**3,nVars )) 
    starttime = MPI_Wtime()
    call ply_pntToLeg_3D( fpt       = fpt,    &
      &                   pntVal    = pntVal, &
      &                   legCoeffs = legVal, &
      &                   nVars     = nVars   )
    stoptime = MPI_Wtime()
    write(*,*) 'Time for degree ', maxpolydegree, ' inverse: ', stoptime - starttime
  
    res = maxval(abs(legVal(:,:) - legCoeffs(:,:)))

  end subroutine ply_check_legToPnt_3D

end program ply_fpt_3D_performance_test
