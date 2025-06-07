! Copyright (c) 2025 Harald Klimach <harald.klimach@dlr.de>
!
! Parts of this file were written by Harald Klimach for DLR e.V.
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
program vtk_ho_hexa_test
  use vtk_ho_hexa_module, only: hexa_pointIndexFromIJK
  use sample_vtkHOhex_module, only: hexorders, indexMap, fillMap

  implicit none

  integer :: i,j,k
  integer :: pind
  logical :: success = .TRUE.

  call fillMap()

  write(*,*) 'Testing the VTK High-Order Hexahedron module...'
  do k=0,hexorders(3)
    do j=0,hexorders(2)
      do i=0,hexorders(1)
        pind = hexa_pointIndexFromIJK(i, j, k, hexorders(1),     &
          &                           hexorders(2), hexorders(3) &
          &                                                      )
        if ( pind /= indexMap(i, j, k) ) then
          success = .FALSE.
        end if
      end do
    end do
  end do
  if (success) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'unexpected index ordering!'
    write(*,*) 'FAILED'
  end if

end program vtk_ho_hexa_test
