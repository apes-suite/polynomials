! Copyright (c) 2012-2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2012 Jan Hueckelheim <j.hueckelheim@grs-sim.de>
! Copyright (c) 2013-2014,2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2018 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
!
! Parts of this file were written by Jens Zudrop, Jan Hueckelheim, Melven
! Zoellner and Harald Klimach for German Research School for Simulation
! Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Verena Krupp, Peter Vitt,
! Daniel Fleischer and Tobias Girresser for University of Siegen.
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

!> Module provides subroutines, functions and datatypes regarding
!! cell local degrees of freedoms.
module ply_dof_module
  use env_module, only: rk

  implicit none

  private

  public :: ply_dof_nextCoeff, ply_dof_pos, ply_dof_count, ply_dof_2degree

  abstract interface
    subroutine ply_dof_nextCoeff(ansFuncX, ansFuncY, ansFuncZ, maxDegree)
      integer, intent(inout) :: ansFuncX, ansFuncY, ansFuncZ
      integer, intent(in) :: maxdegree
    end subroutine ply_dof_nextCoeff

    function ply_dof_pos(ansFuncX, ansFuncY, ansFuncZ, maxDegree) result(pos)
      integer, intent(in) :: ansFuncX, ansFuncY, ansFuncZ
      integer, intent(in) :: maxdegree
      integer :: pos
    end function ply_dof_pos

    function ply_dof_count(maxPolyDegree) result(DoFs)
      integer, intent(in) :: maxPolyDegree
      integer :: DoFs
    end function ply_dof_count
  end interface

  !> Parameter to identify Q polynomials
  integer, public, parameter :: Q_space = 1
  !> Parameter to identify P polynomials
  integer, public, parameter :: P_space = 2

contains

  elemental function ply_dof_2degree(ndofs, space, ndims) result(deg)
    integer, intent(in) :: ndofs
    integer, intent(in) :: space
    integer, intent(in) :: ndims
    integer :: deg

    integer :: idim, fact, estimate

    select case(space)
    case (Q_space)
      deg = nint(ndofs**(1._rk/real(ndims,kind=rk))) - 1
    case (P_space)
      deg = 0
      fact = 1
      do idim=1,ndims
        fact = fact*idim
      end do
      do
        estimate = 1
        do idim=1,ndims
          estimate = estimate * (deg + idim)
        end do
        estimate = estimate/fact
        if (estimate >= ndofs) then
          deg = estimate
          EXIT
        end if
      end do
    end select

  end function ply_dof_2degree

end module ply_dof_module
