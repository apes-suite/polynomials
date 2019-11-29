! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2012-2016, 2018-2019 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013-2014 Verena Krupp <v.krupp@grs-sim.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Kay Langhammer <kay.langhammer@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop for German Research School
! for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Peter Vitt, Verena Krupp,
! Nikhil Anand and Kay Langhammer for University of Siegen.
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

?? include "ply_dof_module.inc"
!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program ply_fpt_3D_test
  use env_module,               only: rk, fin_env
  use tem_param_module,         only: PI
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_general_type, tem_start
  use ply_legFpt_module,        only: ply_legFpt_type, ply_init_legFPT
  use ply_legFpt_3D_module,     only: ply_legToPnt_3D
  use ply_modg_basis_module,    only: evalLegendreTensPoly
  use ply_dof_module,           only: Q_space

  !mpi!nprocs = 1

  implicit none

  integer :: iPower
  real(kind=rk) :: res, newRes
  type(tem_general_type) :: general

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = general             )

  res = 0.0_rk
!do iPower = 0, 4
  do iPower = 2,2
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
    integer :: pointIndex, funcIndex, iPolyX, iPolyY, iPolyZ, iVar, nVars
    real(kind=rk), allocatable :: legCoeffs(:,:)
    real(kind=rk), allocatable :: pntVal(:,:), refVal(:,:)
    real(kind=rk), allocatable :: chebPnt(:,:), chebPnt1D(:)
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
    nVars = 3
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (per dim): ', maxPolyDegree+1
    write(logUnit(10),*) '------------------------------------' // &
      & ' Number of Legendre coefficients (total): ',(maxPolyDegree+1)**3

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs((maxPolyDegree+1)**3, nVars))
    do iVar = 1, nVars
      do iDof = 1, (maxPolyDegree+1)**3
        call random_number(rfac)
        legCoeffs(iDof, iVar) = real(iVar, rk) * rfac
      end do
    end do

    ! Create the Chebyshev nodes on the interval [-1,+1]
    write(logUnit(10),*) 'Creating Chebyshev nodes for ref result ...'
    allocate(chebPnt1D(maxPolyDegree+1))
    do iPoint = 1, maxPolyDegree+1
      chebPnt1D(iPoint) = (-1.0_rk) * &
                & cos(PI/(maxPolyDegree+1)*((iPoint-1.0_rk)+1.0_rk/2.0_rk))
    end do
    ! Now, create the 3D Chebyshev nodes in [-1,+1]^3
    allocate(chebPnt( (maxPolyDegree+1)**3,3 ))
    do iPointX = 1, maxPolyDegree+1
      do iPointY = 1, maxPolyDegree+1
        do iPointZ = 1, maxPolyDegree+1
          pointIndex = 1 + (iPointX-1) + (iPointY-1)*(maxPolyDegree+1) &
                     & + (iPointZ-1)*((maxPolyDegree+1)**2)
          chebPnt(pointIndex,1) = chebPnt1D(iPointX)
          chebPnt(pointIndex,2) = chebPnt1D(iPointY)
          chebPnt(pointIndex,3) = chebPnt1D(iPointZ)
        end do
      end do
    end do

    ! define the reference results for the point values (Chebyshev nodes)
    !allocate( legValChebPnt((maxPolyDegree+1)**3,(maxPolyDegree+1)**3) )
    !legValChebPnt(:,:) = legendre_1D(chebPnt, maxPolyDegree)
    call evalLegendreTensPoly( coords = chebPnt , nCoords = (maxPolyDegree+1)**3, &
                             & maxPolyDegree = maxPolyDegree, &
                             & basisType = Q_space, &
                             & polyVal = legValChebPnt )
    allocate( refVal((maxPolyDegree+1)**3, nVars) )
    refVal(:,:) = 0.0_rk
    write(logUnit(10),*) 'Calculating reference results ...'
    do iPolyX = 1, maxPolyDegree+1
      do iPolyY = 1, maxPolyDegree+1
        do iPolyZ = 1, maxPolyDegree+1
?? copy :: posOfModgCoeffQTens(iPolyX, iPolyY, iPolyZ, maxPolyDegree, funcIndex)
          do iVar = 1, nVars
            refVal(:,iVar) = refVal(:,iVar) + &
                   & legValChebPnt(funcIndex,:) * legCoeffs(funcIndex, iVar)
          end do
        end do
      end do
    end do
    write(logUnit(10),*) 'Finished'

    ! Init the FPT
    call ply_init_legFpt( maxPolyDegree = maxPolyDegree,        &
      &                   nIndeps       = (maxpolydegree+1)**2, &
      &                   fpt           = fpt                   )

    ! now transform to the Chebyshev nodes
    allocate(pntVal( (maxPolyDegree+1)**3, nVars ))
    write(logUnit(10),*) 'Calculating FPT ...'
    call ply_legToPnt_3D( fpt = fpt, legCoeffs = legCoeffs, pntVal = pntVal, &
      &                   nVars = nVars )
    write(logUnit(10),*) 'Finished'

    !!do iPoint = 1, (maxPolyDegree+1)**3
    !!  write(*,*) 'Point: ', iPoint, &
    !!           & ' FPT: ', pntVal(iPoint,3), &
    !!           & ' Ref.: ', refVal(iPoint,3), &
    !!           & ' error: ', pntVal(iPoint,3)-refVal(iPoint,3)
    !!end do

    ! Write out the point with the largest absolute error
    do iVar = 1, nVars
      write(*,*) 'for variable ', iVar, &
             & ' Cheb-Point ',maxloc(abs(pntVal(:, iVar) - refVal(:, iVar))),  &
             & ' has largest error of: ' ,maxval(abs(pntVal(:,iVar) - refVal(:,iVar)))
    end do

    res = maxval(abs(pntVal(:,:) - refVal(:,:)))

  end subroutine

end program ply_fpt_3D_test
