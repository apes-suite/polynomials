!> This is a small utility to approximate a discontinuity at a given
!! location in the interval [-1,1].
!!
!! The jump is located at a given point jota, and the target function is
!! 0 for x <= jota and 1 for x > jota.
!! The interval is subdivided by level bisections towards jota.
!! Intervals for which the lower bound is smaller than jota are assumed to
!! be 0, all others are 1.
!! This interval values are used to determine the values at Chebyshev nodes
!! and a FPT is done to find the Legendre modes.
!!
!! Required settings are therefore:
!! - jota (location of the jump)
!! - level (number of bisections towards jota)
!! - target polynomial degree
!! - oversampling factor
!!
!! These have to be provided as command line parameters in this order.
!!
!! Resulting output is:
!! - Interval subdivision
!! - Approximated jump location
!! - Chebyshev nodes
!! - Legendre modes
!! - L2-Error
program approximate_1D_jump
  use env_module, only: rk
  use ply_dof_module, only: q_space
  use ply_dynArray_project_module, only: ply_prj_init_type
  use ply_poly_project_module, only: ply_poly_project_type, &
    &                                ply_poly_project_fillbody, &
    &                                ply_poly_project_n2m

  implicit none

  character(len=32) :: argstring
  real(kind=rk) :: jota
  integer :: level
  integer :: polydegree

  integer :: ibis
  integer :: ival
  integer :: ipoint
  integer :: imode
  integer :: bis_lb
  real(kind=rk) :: interval_jump
  real(kind=rk) :: ofact
  real(kind=rk) :: bis_half
  real(kind=rk) :: x_point
  real(kind=rk), allocatable :: bisect(:)
  real(kind=rk), allocatable :: nodal_data(:,:)
  real(kind=rk), allocatable :: modal_data(:,:)
  type(ply_prj_init_type) :: project
  type(ply_poly_project_type) :: polypro

  call get_command_argument(1,argstring)
  read(argstring,*) jota
  call get_command_argument(2,argstring)
  read(argstring,*) level
  call get_command_argument(3,argstring)
  read(argstring,*) polydegree
  call get_command_argument(4,argstring)
  read(argstring,*) ofact

  if ( (jota <= -1.0_rk) .or. (jota >= 1.0_rk) ) then
    write(*,*) 'First parameter (jota) needs to be a real value in'
    write(*,*) 'the interval (-1,1), but is:', jota
    write(*,*) 'Stopping!'
    STOP
  end if

  if (level < 0) then
    write(*,*) 'Second parameter (level) needs to be non-negative!'
    write(*,*) 'But it is:', level
    write(*,*) 'Stopping!'
    STOP
  end if

  if (polydegree < 0) then
    write(*,*) 'Third parameter (polydegree) needs to be non-negative!'
    write(*,*) 'But it is:', polydegree
    write(*,*) 'Stopping!'
    STOP
  end if

  if (ofact < 0.0_rk) then
    write(*,*) 'Fourth parameter (oversampling) needs to be non-negative!'
    write(*,*) 'But it is:', ofact
    write(*,*) 'Stopping!'
    STOP
  end if

  allocate(bisect(0:level+1))

  bisect = 1.0_rk
  bisect(0) = -1.0_rk
  bis_lb = 0
  do ibis=1,level

    bis_half = 0.5_rk * ( bisect(bis_lb+1) - bisect(bis_lb) )
    do ival=level,bis_lb+1,-1
      bisect(ival+1) = bisect(ival)
    end do
    bisect(bis_lb+1) = bisect(bis_lb) + bis_half
    if (bisect(bis_lb+1) < jota) then
      bis_lb = bis_lb+1
    end if

  end do

  write(*,*) 'Bisection points:'
  do ibis=0,level+1
    write(*,*) bisect(ibis)
  end do
  write(*,*) '-----------------'

  ival = 0
  do
    if (bisect(iVal) > jota) EXIT
    iVal = iVal + 1
  end do
  interval_jump = bisect(ival)

  write(*,*) 'Approximated jump location:', interval_jump

  project%basisType = q_space
  project%maxPolyDegree = polydegree
  project%header%kind = 'fpt'
  project%header%fpt_header%factor = ofact
  project%header%fpt_header%approx_terms = 18
  project%header%fpt_header%striplen = 128
  project%header%fpt_header%adapt_factor_pow2 = .false.
  project%header%fpt_header%nodes_header%lobattopoints = .false.
  project%header%fpt_header%nodes_header%nodes_kind = 'chebyshev'
  call ply_poly_project_fillbody( me         = polypro, &
    &                             proj_init  = project, &
    &                             scheme_dim = 1        )

  allocate(nodal_data(polypro%nQuadPointsPerDir,1))
  allocate(modal_data(polypro%nQuadPointsPerDir,1))

  write(*,*)
  write(*,*) 'Chebyshev points:'
  do ipoint=1,polypro%nQuadPointsPerDir
    x_point = polypro%body_1D%nodes(ipoint,1)
    write(*,*) x_point
    if (x_point >= interval_jump) then
      nodal_data(ipoint,1) = 1.0_rk
    else
      nodal_data(ipoint,1) = 0.0_rk
    end if
  end do
  write(*,*) '-----------------'

  call ply_poly_project_n2m( me         = polypro,    &
    &                        dim        = 1,          &
    &                        nVars      = 1,          &
    &                        nodal_data = nodal_data, &
    &                        modal_data = modal_data  )

  write(*,*) 'Legendre modes:'
  do imode=1,polydegree+1
    write(*,*) modal_data(imode,1)
  end do
  write(*,*) '-----------------'

end program approximate_1D_jump
