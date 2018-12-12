!> Testing functions from the polynomial base exchange module.
program ply_polyBaseExc_test
  use env_module, only: rk
  use ply_polyBaseExc_module, only: ply_lambda

  implicit none

  real(kind=rk), parameter :: max_deviation = 5.0E-15_rk
  real(kind=rk) :: testnums(10) = [  5._rk, 10._rk, 15._rk, 17._rk, 20._rk, &
    &                               25._rk, 40._rk, 50._rk, 75._rk, 90._rk  ]
  real(kind=rk) :: lam_res(10)
  real(kind=rk) :: ref_lam(10) = [ 0.43618981487127939_rk, &
    &                              0.31230114333906128_rk, &
    &                              0.25605656734380250_rk, &
    &                              0.24075907021388787_rk, &
    &                              0.22221375806199573_rk, &
    &                              0.19900256214091094_rk, &
    &                              0.15762056118567075_rk, &
    &                              0.14106825029753825_rk, &
    &                              0.11527776545731630_rk, &
    &                              0.10526295596826729_rk  ]

  integer :: inum

  lam_res = ply_lambda(testnums)
  do inum=1,10
    write(*,*) 'lambda(', testnums(inum), ') =', lam_res(inum)
  end do
  write(*,*) ''
  if ( maxval(abs(lam_res - ref_lam)) <= max_deviation ) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program ply_polyBaseExc_test
