program ply_transfer_test
  use env_module, only: stdOutUnit
  use ply_transfer_test_module

  implicit none


  logical :: test_ok
  logical :: passed_tests

  call ply_test_transfer_1d( success = test_ok,   &
    &                        lu      = stdOutUnit )

  passed_tests = test_ok

  call ply_test_transfer_2d( success = test_ok,   &
    &                        lu      = stdOutUnit )

  passed_tests = passed_tests .and. test_ok

  if (passed_tests) write(*,*) 'PASSED'

end program ply_transfer_test
