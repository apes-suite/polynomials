module fftw_wrap
  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  logical, parameter :: fftw_available = .true.

end module fftw_wrap
