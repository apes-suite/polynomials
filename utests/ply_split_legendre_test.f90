program ply_split_legendre_test_prog
  use ply_split_legendre_module, only: ply_split_legendre_test
  implicit none

  logical :: success

  call ply_split_legendre_test(success)

  if (success) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program ply_split_legendre_test_prog
