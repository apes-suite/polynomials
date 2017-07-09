program ply_split_element_test_prog
  use env_module, only: rk
  use ply_split_element_module, only: ply_split_element_test

  implicit none


  logical :: success

  call ply_split_element_test(success)

  if (success) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program ply_split_element_test_prog
