! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012-2016, 2018-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop and Harald Klimach for
! German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Peter Vitt, Verena Krupp,
! and Nikhil Anand for University of Siegen.
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !

!> Unit test to check performance of l2p transformations.
program ply_l2p_3D_performance_test
  use mpi, only: mpi_wtime
  use env_module,            only: rk, fin_env
  use tem_logging_module,    only: logUnit, tem_logging_init_primary
  use tem_general_module,    only: tem_general_type, tem_start
  use ply_l2p_header_module, only: ply_l2p_header_type, ply_l2p_header_define
  use ply_l2p_module,        only: ply_l2p_type, ply_init_l2p, ply_l2p_trafo_3D

  !mpi!nprocs = 1

  implicit none

  integer :: iPower
  integer, parameter :: maxpower = 7
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'L2P 3D Performance Test', &
    &            version  = '1',                       &
    &            general  = general                    )
  call tem_logging_init_primary( level = 1,               &
    &                            rank = general%proc%rank )

  res = 0.0_rk
  do iPower = 1,maxpower
    call ply_check_legToPnt_3D(iPower, newRes)
    write(*,*) 'deviation:', newres
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
    type(ply_l2p_header_type) :: header
    type(ply_l2p_type) :: trafo
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

    ! Init the L2 Projection
    call ply_l2p_header_define( me         = header,     &
      &                         nodes_kind = 'chebyshev' )
    call ply_init_l2p( degree = maxPolyDegree, &
      &                l2p    = trafo,         &
      &                header = header         )

    ! now transform to the Chebyshev nodes
    allocate(pntVal( (maxPolyDegree+1)**3, nVars ))
    legCoeffsIn = legCoeffs
    starttime = MPI_Wtime()
    do iVar = 1, nVars
      call ply_l2p_trafo_3D( trafo     = trafo%leg2node,      &
        &                    original  = legCoeffsIn(:,iVar), &
        &                    projected = pntVal(:,iVar)       )
    end do
    stoptime = MPI_Wtime()
    write(*,*) 'Time for degree ', maxpolydegree, ' trafo:   ', &
      &        stoptime - starttime

    ! now transform back to Legendre coefficients
    allocate(legVal( (maxPolyDegree+1)**3,nVars ))
    starttime = MPI_Wtime()
    do iVar = 1, nVars
      call ply_l2p_trafo_3D( trafo     = trafo%node2leg, &
        &                    projected = legVal(:,iVar), &
        &                    original  = pntVal(:,iVar)  )
    end do
    stoptime = MPI_Wtime()
    write(*,*) 'Time for degree ', maxpolydegree, ' inverse: ', &
      &        stoptime - starttime
  
    res = maxval(abs(legVal(:,:) - legCoeffs(:,:)))

  end subroutine ply_check_legToPnt_3D

end program ply_l2p_3D_performance_test
