?? include "ply_dof_module.inc"
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_3D_module
  use, intrinsic :: iso_c_binding
  !$ use omp_lib
  use env_module, only: rk
  use tem_aux_module, only: tem_abort
  use tem_logging_module, only: logUnit
  use ply_polyBaseExc_module, only: ply_trafo_params_type, &
                                  & ply_fpt_init, &
                                  & ply_fpt_exec_striped, &
                                  & ply_fpt_exec, & 
                                  & ply_legToCheb_param, ply_chebToLeg_param, &
                                  & assignment(=)
  use fftw_wrap
  use ply_nodes_module,        only: ply_faceNodes_type
  use ply_legFpt_module,       only: ply_legFpt_type

  implicit none

  private

  interface ply_LegTopnt_3D
    module procedure ply_LegTopnt_3D_multVar
    module procedure ply_LegTopnt_3D_singVar
  end interface ply_LegToPnt_3D

  interface ply_pntToLeg_3D
    module procedure ply_pntToLeg_3D_multVar
    module procedure ply_pntToLeg_3D_singVar
  end interface ply_pntToLeg_3D

  public :: ply_init_legFpt_3D, ply_legToPnt_3D, ply_pntToLeg_3D


contains


  !****************************************************************************
  !> Subroutine to initialize the fast polynomial transformation
  !! for Legendre expansion.
  !! \todo VK: change maxPolyDegree to number of quadpoints to include the 
  !! dealising factor 
  subroutine ply_init_legFpt_3D( maxPolyDegree, nVars, fpt, &
    &                            blocksize, approx_terms, striplen, &
    &                            lobattoPoints, subblockingWidth )
    !---------------------------------------------------------------------------
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: nVars
    type(ply_legFpt_type), intent(inout) :: fpt

    !> Smallest block that is approximated by approx_terms coefficients.
    !!
    !! Please note, that this has to be larger than 2*approx_terms to result
    !! in a reduced number of operations. Default is 64.
    integer, optional, intent(in) :: blocksize

    !> Number of approximation terms used to compute off-diagonal products.
    !!
    !! Defaults to 18, which is the suggested accuracy for double precision.
    integer, optional, intent(in) :: approx_terms

    !> Length to use in vectorization, this is the number of independent
    !! matrix multiplications that are to be done simultaneously.
    !!
    !! Defaults to 512, but the optimal setting is platform specific.
    integer, optional, intent(in) :: striplen

    !> Use Chebyshev-Lobtatto Points (true) or simple Chebyshev points (false)
    logical, intent(in), optional  :: lobattoPoints
    !> The width of the subblocks used during the unrolled base exchange to 
    !! ensure a better cache usage.
    integer, optional, intent(in) :: subblockingWidth
    !---------------------------------------------------------------------------
    integer :: n_cubed
    real(kind=rk), allocatable :: tmpOut(:,:), tmpIn(:,:)
    logical :: lob
    !$ integer :: fftwMultThread
    !---------------------------------------------------------------------------

    lob = .false.
    if (present(lobattoPoints)) then
      lob = lobattoPoints
    end if

    ! init the fast Legendre to Chebyshev transformation.
    call ply_fpt_init( n                = maxPolyDegree+1,     &
      &                params           = fpt%legToChebParams, &
      &                trafo            = ply_legToCheb_param, &
      &                blocksize        = blocksize,           &
      &                approx_terms     = approx_terms,        &
      &                striplen         = striplen,            &
      &                subblockingWidth = subblockingWidth     )

    call ply_fpt_init( n                = maxPolyDegree+1,     &
      &                params           = fpt%chebToLegParams, &
      &                trafo            = ply_chebToLeg_param, &
      &                blocksize        = blocksize,           &
      &                approx_terms     = approx_terms,        &
      &                striplen         = striplen,            &
      &                subblockingWidth = subblockingWidth     )

    ! Create the buffers for the intermediate arrays
    n_cubed = fpt%legToChebParams%n**3

    ! Temporary arrays to initialize FFTW real->real transformations
    allocate( tmpIn(n_cubed, nVars) )
    allocate( tmpOut(n_cubed, nVars) )
 
    ! If we have OpenMP parallelism, we have to init FFTW with the
    ! corrseponding function call
    !$ fftwMultThread = fftw_init_threads() 
    !$ if(fftwMultThread.eq.0) then
    !$   write(logUnit(1),*) 'ERROR in ply_init_legFpt_3D: ' &
    !$     &                 //'Not able to init OpenMP parallel FFTW'
    !$   write(logUnit(1),*) 'Stopping!'
    !$   call tem_abort()    
    !$ end if

    ! Tell FFTW how many threads at max we want to use: at max
    ! we use the number of OMP threads in the current team.
    !$ call fftw_plan_with_nthreads( omp_get_max_threads() )

    if (.not.lob) then

      ! Init the DCT III ( Leg -> Point values )
      fpt%planChebToPnt = fftw_plan_r2r_3d( n0    = fpt%legToChebParams%n, &
        &                                   n1    = fpt%legToChebParams%n, &
        &                                   n2    = fpt%legToChebParams%n, &
        &                                   in    = tmpIn,                 &
        &                                   out   = tmpOut,                &
        &                                   kind0 = FFTW_REDFT01,          &
        &                                   kind1 = FFTW_REDFT01,          &
        &                                   kind2 = FFTW_REDFT01,          &
        &                                   flags = FFTW_ESTIMATE          )

      ! Init the DCT II ( Point values -> Leg )
      fpt%planPntToCheb = fftw_plan_r2r_3d( n0    = fpt%legToChebParams%n, &
        &                                   n1    = fpt%legToChebParams%n, &
        &                                   n2    = fpt%legToChebParams%n, &
        &                                   in    = tmpIn,                 &
        &                                   out   = tmpOut,                &
        &                                   kind0 = FFTW_REDFT10,          &
        &                                   kind1 = FFTW_REDFT10,          &
        &                                   kind2 = FFTW_REDFT10,          &
        &                                   flags = FFTW_ESTIMATE          )

    else

      ! Init the DCT I ( Leg -> Point values )
      fpt%planChebToPnt = fftw_plan_r2r_3d( n0    = fpt%legToChebParams%n, &
        &                                   n1    = fpt%legToChebParams%n, &
        &                                   n2    = fpt%legToChebParams%n, &
        &                                   in    = tmpIn,                 &
        &                                   out   = tmpOut,                &
        &                                   kind0 = FFTW_REDFT00,          &
        &                                   kind1 = FFTW_REDFT00,          &
        &                                   kind2 = FFTW_REDFT00,          &
        &                                   flags = FFTW_ESTIMATE          )

      ! Init the DCT I ( Point values -> Leg )
      fpt%planPntToCheb = fftw_plan_r2r_3d( n0    = fpt%legToChebParams%n, &
        &                                   n1    = fpt%legToChebParams%n, &
        &                                   n2    = fpt%legToChebParams%n, &
        &                                   in    = tmpIn,                 &
        &                                   out   = tmpOut,                &
        &                                   kind0 = FFTW_REDFT00,          &
        &                                   kind1 = FFTW_REDFT00,          &
        &                                   kind2 = FFTW_REDFT00,          &
        &                                   flags = FFTW_ESTIMATE          )

    end if

    deallocate( tmpIn, tmpOut )

  end subroutine ply_init_legFpt_3D
  !****************************************************************************


  !****************************************************************************
  subroutine ply_legToPnt_3D_singvar( fpt, legCoeffs, pntVal, lobattoPoints )
   !---------------------------------------------------------------------------
   !> The FPT parameters.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified 
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:) 
   real(kind=rk), intent(inout) :: pntVal(:)
   logical, intent(in)  :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iFuncX, iFuncY, iFuncZ, funcIndex
   integer :: iFunc, iDof

   integer :: n, n_squared, n_cubed
   integer :: nIndeps
   ! The normfactor has to be positive for odd sum of exponents and negative
   ! for even sums of exponents.
   ! (original definition: normfactor = (-1)**(i-1) * (-1)**(j-1) * (-1)**(k-1)
   !  can be written as (-1)**(i+j+k-3) = (-1)**(i+j+k-1) for i+j+k>=3
   !  the power can be avoided completely by normfactor = 1 - 2*mod(i+j+k-1,2)
   !  and the real multiplication can be avoided by putting this into a lookup
   !  table.)
   real(kind=rk), parameter :: normFactor(0:1) = [ -0.125_rk, 0.125_rk ]
   integer :: iStrip
   integer :: strip_ub
   integer :: striplen
   integer :: stride
   integer :: alph_lb, alph_ub, iAlph, jAlph
   integer :: iDim
   integer :: itest
   integer :: linesPerStrip
   integer :: jPerStrip
   integer :: linesInAlph
   integer :: iIndex
!   real(kind=rk), dimension (fpt%legToChebParams%striplen) :: alph
!   real(kind=rk), dimension (fpt%legToChebParams%striplen) :: gam
   real(kind=rk), dimension(:), allocatable :: alph
   real(kind=rk), dimension(:), allocatable :: gam
   integer, dimension(:), allocatable :: alphIndices
   !---------------------------------------------------------------------------

   striplen = fpt%legToChebParams%striplen
   n = fpt%legToChebParams%n
   n_squared = fpt%legToChebParams%n**2
   n_cubed = n_squared * fpt%legToChebParams%n
   
   ! number of strips to execute the fpt on in one call
   ! (usually n_squared, but a smaller value might be assigned at the end of 
   ! the array)
   nIndeps = n_squared !'

!   allocate(alph(min(striplen, fpt%legToChebParams%n)))
!   allocate(gam(min(striplen, fpt%legToChebParams%n)))
   allocate(alph(min(striplen,n_squared)*n))
   allocate(gam(min(striplen,n_squared)*n))
   allocate(alphIndices(n_squared))
   ! Dimension-by-dimension transform Legendre expansion to Chebyshev expansion
   ! ... transformation in X direction (Leg->Cheb)
!'   call ply_fpt_exec_striped( nIndeps = n_squared,          &
!'     &                        alph    = legCoeffs,          &
!'     &                        gam     = pntVal,             &
!'     &                        params  = fpt%legToChebParams )
!'   ! ... transformation in Y direction (Leg->Cheb)
!'   call ply_fpt_exec_striped( nIndeps = n_squared,          &
!'     &                        alph    = pntVal,             &
!'     &                        gam     = legCoeffs,          &
!'     &                        params  = fpt%legToChebParams )
!'   ! ... transformation in Z direction (Leg->Cheb)
!'   call ply_fpt_exec_striped( nIndeps = n_squared,          &
!'     &                        alph    = legCoeffs,          &
!'     &                        gam     = pntVal,             &
!'     &                        params  = fpt%legToChebParams )

! stride for reading
  !\ testloop
  write(*,*)'striplen', striplen
  write(*,*)'initial values'
  do iTest = 1,n_cubed
    write(*,*)'iTest=', iTest, 'pntVal=', pntVal(itest),'legCoeff=', legCoeffs(itest)
  end do

  ! z-direction
!   stride = n_squared    
     ! original layout (n = 3):
     !  1  2  3   10 11 12   19 20 21
     !  4  5  6   13 14 15   22 23 24
     !  7  8  9   16 17 18   25 26 27
     ! 
     ! example for striplen = 2:
     ! iStrip           1                   3
     ! iAlph            1         2         3         4
     ! index legCoeffs  1 10 19   2 11 20   3 12 21   4 13 22  
     ! index alph/gam   1  2  3   4  5  6   1  2  3   4  5  6
     ! index PntVal     1  2  3   4  5  6   7  8  9  10 11 12
     !
     ! layout after z-trafo: (n = 3)
     !  1 10 19    4 13 22    7 16 25 
     !  2 11 20    5 14 23    8 17 26
     !  3 12 21    6 15 24    9 18 27

   ! zStripLoop: Loop over all strips in z-direction
   zStripLoop: do iStrip = 1, n_squared, striplen
!'     write(*,*)'iStrip', iStrip
     ! iAlph is the index of the first element in a line for the transformation in 
     ! z-direction. 
     do iAlph = iStrip, min(iStrip+striplen-1, n_squared)  !z_Trafo
       
       ! alph_ub is the index of the last element in a line for the transformation 
       ! in z-direction. 
!       write(*,*)'n_cubed', n_cubed
       alph_ub = iAlph + n_cubed - n_squared  !z-trafo

       ! legCoeffs -> copy data for 1D-transformation into alph (strided reading)
!       write(*,*)'iAlph', iAlph
!       write(*,*)'alph Index', iAlph-striplen*(iStrip-1),iAlph+n-striplen*(iStrip-1)-1
!       write(*,*)'legCoeffs Index',iAlph,alph_ub,n_squared
!       alph((iAlph-1)*n+1-striplen*(iStrip-1):(iAlph-1)*n+n-striplen*(iStrip-1)) = &
!'       write(*,*)'(iAlph-iStrip)*n+1, (iAlph-iStrip+1)*n', (iAlph-iStrip)*n+1, (iAlph-iStrip+1)*n
!'       write(*,*) 'iAlph, alph_ub, n_squared', iAlph, alph_ub, n_squared
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = &
           & legCoeffs(iAlph:alph_ub:n_squared) !ztrafo

        !\ testloop
!        write(*,*)'after ialph', iAlph
!        do iTest = 1,n_cubed
!          write(*,*)'iTest=', iTest, 'alph=', alph(itest),'legCoeff=', legCoeffs(itest)
!        end do

     end do

  !\ testloop before
  write(*,*)'values before fpt'
  do iTest = 1,n_cubed
    write(*,*)'iTest=', iTest, 'alph=', alph(itest),'gam=', gam(itest)
  end do

       ! At the end of the array the number of computed strips might be smaller
       nIndeps = min(striplen, n_squared-iStrip+1)

       ! ply_fpt_exec on temp (no memory transpose)
       call ply_fpt_exec( alph = alph,                  &
        &                 gam = gam,                    &
        &                 nIndeps = nIndeps,            &
        &                 params = fpt%legToChebParams  )
!'         gam = alph !\
  !\ testloop after
  write(*,*)'values after fpt'
  do iTest = 1,n_cubed
    write(*,*)'iTest=', iTest, 'alph=', alph(itest),'gam=', gam(itest)
  end do

!         pntVal((iStrip-1)*striplen*n+1:min(iStrip*striplen*n, n_cubed))  = gam(:)
         pntVal((iStrip-1)*n+1:min((iStrip+striplen-1)*n, n_cubed))  = gam(:)

       ! todo: fft on temp
       ! temp -> pntVal (stride-1 writing)
!      write(*,*) 'iAlph', iAlph
!      write(*,*) 'alph_ub', alph_ub
!      write(*,*) 'size of alph', size(alph)
!      write(*,*) 'size of pntval dim1', size(pntVal, 1)
!      write(*,*) 'shape of pntval', shape(pntVal)
!      write(*,*) 'size of gam', size(gam)
!      write(*,*) 'stride', stride
!      write(*,*) 'iStrip', iStrip
!      write(*,*) 'linesPerStrip', linesPerStrip
!      write(*,*) 'min(iStrip+linesPerStrip, n_squared)', &
!                  min(iStrip+linesPerStrip, n_squared)
!      write(*,*) 'pntVal(iAlph:alph_ub:stride)', & 
!                  pntVal(iAlph:alph_ub:stride)
!      write(*,*) 'alph', alph(:)
!      write(*,*) 'gam', gam(:)

     ! Write the temporary output into pntval
!'     do iAlph = iStrip, min(iStrip+linesPerStrip, n_squared)
!        write(*,*) 'ialph', iAlph
!       write(*,*) 'min((iStrip-1)+linesPerStrip*n+1, n_cubed-n_squared)', & 
!        &  min((iStrip-1)*linesPerStrip*n+1, n_cubed-n_squared) 
!       write(*,*) 'iStrip*linesPerStrip*n', iStrip*linesPerStrip*n

!       do iTEST = 1,size(pntVal)
!         write(*,*) 'pntVal', iTEST, pntval(iTEST)
!       end do

 !      write(*,*) 'pntVal(min((iStrip-1)+linesPerStrip*n+1, n_cubed-n_squared) :  &
 !        &      min(iStrip*linesPerStrip*n, n_cubed)', &
 !       & pntVal((iStrip-1)*linesPerStrip*n+1 : min(iStrip*linesPerStrip*n, n_cubed))
!        pntVal(iAlph:alph_ub:stride) = gam((iAlph-1)*stride+1:iAlph*stride)
!'     end do

   end do zStripLoop

  !\ testloop
  write(*,*)'after z-trafo'
  do iTest = 1,n_cubed
    write(*,*)'iTest=', iTest, 'pntVal=', pntVal(itest),'legCoeff=', legCoeffs(itest)
  end do

  ! y-direction

     ! alphIndices: Array which contains the values for iAlpha.
     j_yLoop: do jAlph = 1, n
       i_yLoop: do iAlph = jAlph, jAlph+n_squared-n, n
          alphIndices((jAlph-1)*n + (iAlph-jAlph)/n +1) = iAlph
       end do i_yLoop
     end do j_yLoop

     !  .  .  .    .  .  .    .  .  .
     ! original layout (n = 3):
     !  1  2  3   10 11 12   19 20 21
     !  4  5  6   13 14 15   22 23 24
     !  7  8  9   16 17 18   25 26 27
     !
     ! layout after z-trafo: (n = 3)
     !  1 10 19    4 13 22    7 16 25 
     !  2 11 20    5 14 23    8 17 26
     !  3 12 21    6 15 24    9 18 27

     ! todo For y-trafo:
     !  val  =  1  4  7    2  5  8    3  6  9   10 13 16    
     !  iStrip= 1                     3                  for striplen = 2
     !  jAlph = 1                                2         {1, 2, 3}
     !  iIndex= 1          2          3          4   
     !  iAlph = 1          4          7          2        
     !  index = 1 10 19    4 13 22    7 16 25    2 11 16  (pntVal)
     !  index = 1  2  3    4  5  6    1  2  3    4  5  6  (alph, gam)
     !  index = 1  2  3    4  5  6    7  8  9   10 11 12  (legCoeffs)
     !  
     ! layout after y-trafo:
     !  1  4  7   10 13 16   19 22 25
     !  2  5  8   11 14 17   20 23 26
     !  3  6  9   12 15 18   21 24 27
     !
     ! original layout after x trafo: 
     !  1  2  3   10 11 12   19 20 21
     !  4  5  6   13 14 15   22 23 24
     !  7  8  9   16 17 18   25 26 27
     
     ! y-direction
!     alph_lb = 1


   yStripLoop: do iStrip = 1,n_squared,striplen
      do iIndex = iStrip,min(iStrip+striplen, n_squared)
        iAlph = alphIndices(iIndex)
        alph((iIndex-iStrip)*n+1 : min((iIndex-iStrip)*n+n, n_cubed)) = &
           &  pntVal(iAlph:iAlph+n_cubed-n_squared:n_squared)
      end do

       ! At the end of the array the number of computed strips might be smaller
       nIndeps = min(striplen, n_squared-iStrip+1)

         ! ply_fpt_exec on temp (no memory transpose)
         call ply_fpt_exec( alph = alph,                  &
           &                gam = gam,                    &
           &                nIndeps = nIndeps,            &
        &                 params = fpt%legToChebParams    )
 
         ! todo: fft on temp
         ! temp -> pntVal (stride-1 writing)

        legCoeffs((iStrip-1)*n+1 : min((iStrip-1+striplen)*n, n_cubed)) &
          &        = gam(:)

   end do yStripLoop ! iStrip

!"2  y-direction for leg_fpt_2D
!"2     jLoop: do jAlph = (iStrip-1)*jPerStrip, iStrip*jPerStrip
!"2       iLoop: do iAlph = (jAlph-1)*n_squared+1, (jAlph-1)*n_squared+n
!"2         alph(iAlph) = pntVal(iAlph:iAlph+n_cubed-n_squared:n_squared)
!"2       end do iLoop
!"2     end do jLoop

     !"1 There are two loops for the transformation in y-direction:
     !"1 e.g. for n = 3
     !"1   |  1,  4,  7;  2,  5,  8;  3,  6,  9;|  j = 1
     !"1   | 10, 13, 16; 11, 14, 17; 12, 15, 18;|  j = 2 
     !"1   | 19, 22, 25; 20, 23, 26; 21, 24, 27 |  j = 3
     !"1
     !"1 i = 1, 2, 3, 10, 11, 12, 19, 20, 21
     !"1   The outerLoop(jLoop) has 3 iterations        
     
!"1     jLoop: do jAlph = (iStrip-1)*jPerStrip+1, iStrip*jPerStrip
!"1       iLoop: do iAlph = (jAlph-1)*n_squared+1, (jAlph-1)*n_squared+n
!"1         alph((iLine-1)*n+1:iLine*n) = pntVal(iAlph:iAlph+(n_squared-n):n)
!"1         iLine = iLine+1
!"1       end do iLoop
!"1     end do jLoop
!     jLoop: do jAlph = 1, n_cubed - n + iAlph, n_squared ! y-trafo  
!       iLoop: do iAlph = jAlph, jAlph + n - 1            ! y-trafo  
!  !      alph_ub = iAlph + stride*(n-1)
!         alph_ub = iAlph + n_squared - n 
!  
!         ! legCoeffs -> copy data for 1D-transformation into alph (strided reading)
!         alph(iStrip:iStrip+n-1) = legCoeffs(iAlph:alph_ub:stride) 
!  
!         ! Exit Loop if array is full
!         linesInAlph = linesInAlph+1
!         if (linesInAlph.ge.linesPerStrip) exit jLoop
!       end do iLoop
!     end do jLoop
    

!'         legCoeffs(iAlph:alph_ub:stride) = gam(:)
!'
!'  
!'         ! exit strip-loop if when y-trafo is completed
!'     do iStrip = 1,ceiling(real(n_cubed)/linesPerSt<F2>rip)    
!'         if (iStrip.ge.ceiling(real(n_cubed)/linesPerStrip)) exit

  ! x-direction
     ! layout after y-trafo:
     !  1  4  7   10 13 16   19 22 25
     !  2  5  8   11 14 17   20 23 26
     !  3  6  9   12 15 18   21 24 27
     !
     ! x-trafo:
     ! iStrip= 1                   3                     5
     ! jAlph = 1                             10         
     ! iAlph = 1         2         3         10         11
     ! iIndex= 1         2         3          4          5 
     ! index = 1  4  7   2  5  8   3  6  9   10 13 16   11 14 17 (pntVal)
     ! index = 1  2  3   4  5  6   1  2  3    4  5  6    1  2  3 (alph/gam)
     ! index = 1  2  3   4  5  6   7  8  9   10 11 12   13 14 15 (legCoeffs)

     ! original layout after x trafo: 
     !  1  2  3   10 11 12   19 20 21
     !  4  5  6   13 14 15   22 23 24
     !  7  8  9   16 17 18   25 26 27

  !\ testloop
  write(*,*)'after y-trafo'
  do iTest = 1,n_cubed
    write(*,*)'iTest=', iTest, 'pntVal=', pntVal(itest),'legCoeff=', legCoeffs(itest)
  end do
   ! x-trafo
   !  (also y-direction for leg_fpt_2D)

   ! alphIndices: Array which contains the values for iAlpha.
   j_xLoop: do jAlph = 1, n_cubed-n_squared+1, n_squared
!     write(*,*)'jAlph', jAlph
     i_xLoop: do iAlph = jAlph, jAlph+n-1
!       write(*,*)'iAlph', iAlph
        iIndex = (jAlph-1)/n + iAlph-jAlph+1
!        write(*,*)'iIndex', iIndex
        alphIndices(iIndex) = iAlph
!'   write(*,*)'iIndex', iIndex, 'alphIndices', alphIndices
     end do i_xLoop
   end do j_xLoop


   xStripLoop: do iStrip = 1,n_squared,striplen
      do iIndex = iStrip, min(iStrip+striplen, n_squared)
        iAlph = alphIndices(iIndex)
        alph((iIndex-iStrip)*n+1 : min((iIndex-iStrip)*n+n, n_cubed)) = &
           &  legCoeffs(iAlph:iAlph+n_squared-n:n)
      end do

       ! At the end of the array the number of computed strips might be smaller
       nIndeps = min(striplen, n_squared-iStrip+1)

!   xStripLoop: do iStrip = 1, n_squared, striplen
!     alph(iStrip*) = legCoeffs(alphIndices(iIndex))  !
       ! ply_fpt_exec on temp (no memory transpose)
       call ply_fpt_exec( alph = alph,                  &
         &                gam = gam,                    &
         &                nIndeps = nIndeps,            &
         &                 params = fpt%legToChebParams )

     ! todo: fft on temp
     ! temp -> pntVal (stride-1 writing)
     pntVal(iStrip : min(iStrip+striplen*n, n_cubed)) = gam(:)
   end do xStripLoop

  !\ testloop
  write(*,*)'after x-trafo'
  do iTest = 1,n_cubed
    write(*,*)'iTest=', iTest, 'pntVal=', pntVal(itest),'legCoeff=', legCoeffs(itest)
  end do

!!     alph_lb = (iStrip-1)*jPerStrip*n_squared+1
!     j_xLoop: do jAlph = (iStrip-1)*jPerStrip+1, iStrip*jPerStrip
!!         write(*,*)'jAlph', jAlph
!       i_xLoop: do iAlph = (jAlph-1)*n_squared+1, (jAlph-1)*n_squared+n
!!         write(*,*)'iAlph', iAlph
!!         write(*,*)'alph_lb+1:alph_lb+n', alph_lb+1, alph_lb+n
!!         write(*,*)'iAlph: iAlph+n_squared-n: n', iAlph, iAlph+n_squared-n, n
!         alph(alph_lb:alph_lb+n-1) = legCoeffs(iAlph: iAlph+n_squared-n: n)
!         alph_lb = alph_lb+n
!       end do i_xLoop
!     end do j_xLoop

!"2   x-trafo
!"2   do iStrip = 1,ceiling(real(n_cubed)/(linesPerStrip*n))    

!'     do iAlph = (iStrip-1)*n+1, min(iStrip*striplen, n_cubed) !x-trafo 
  !    alph_ub = iAlph + stride*(n-1)
!       alph_ub = iAlph + n - 1 ! x-trafo
  
       ! legCoeffs -> copy data for 1D-transformation into alph (stride=1 reading)
!"2     do iAlph = (iStrip-1)*linesPerStrip*n+1, n_cubed-n+1
!"2           alph((iStrip-1)*linesPerStrip*n+1:min(iStrip*linesPerStrip*n, n_cubed)) = &
!"2           & legCoeffs((iStrip-1)*linesPerStrip*n+1:n_cubed:n)
!"2     end do
!"1       ! legCoeffs -> copy data for 1D-transformation into alph (stride=1 reading)
!"1       alph(iAlph:iAlph+n) = legCoeffs(iAlph:iAlph+n) 
     
!       ! ply_fpt_exec on temp (no memory transpose)
!       call ply_fpt_exec( alph = alph,                  &
!         &                gam = gam,                    &
!         &                params = fpt%legToChebParams  )
! 
!     ! todo: fft on temp
!     ! temp -> pntVal (stride-1 writing)
!     pntVal((iStrip-1)*striplen+1:min(iStrip*striplen, n_cubed)) = gam(:)
!   end do xStripLoop

!    write(*,*) 'n', fpt%legToChebParams%n
!    write(*,*) 'striplen', striplen
!    write(*,*) 'alph_lb', alph_lb
!    write(*,*) 'alph_ub', alph_ub
!    write(*,*) 'istrip', istrip
!    write(*,*) 'strip_ub', strip_ub
!    write(*,*) 'stride', stride
!    write(*,*) 'size of alph', size(alph)
!    write(*,*) 'size of legcoeffs', size(legCoeffs)
!    write(*,*) 'size of legcoeff slices', size(legCoeffs(alph_lb:alph_ub:stride))


    !'gam -> legCoeffs (stride-1 writing)
    ! (just for the test, double memory transpose)
    ! pntVal(iStrip:strip_ub) = gam

!' ! y-direction
!' stride = fpt%legToChebParams%n
!' 
!' do iStrip = 0,nIndeps-1,striplen 
!'   ! Calculate the upper bound of the current strip
!'   strip_ub = min(iStrip + striplen, nIndeps)
!' 
!'   ! pntVal -> copy data for Y-transformation into temp (strided reading)
!'   temp = pntVal(iStrip:fpt%legToChebParams%n:n_cubed)
!' 
!'   ! ply_fpt_exec on temp (no memory transpose)
!'   call ply_fpt_exec( nIndeps = nIndeps,            & 
!'                      alph = temp,                  &
!'                      gam = legCoeffs,              &
!'                      params = fpt%legToChebParams, &
!'                      strip_lb = iStrip,            &
!'                      strip_ub = strip_ub           )
!'   ! todo: fft on temp
!'   ! temp -> legCoeffs (stride-1 writing)
!'   legCoeffs(iStrip:strip_ub) = temp
!' end do
!' 
!' ! x-direction
!' do iStrip = 0,nIndeps-1,striplen 
!'   ! legCoeffs -> copy data for X-transformation into temp (strided reading)
!'   strip_ub = min(iStrip + striplen, nIndeps)
!'   temp = legCoeffs(iStrip:strip_ub)
!'   ! ply_fpt_exec on temp (no memory transpose)
!'   call ply_fpt_exec( nIndeps = nIndeps,            & 
!'                      alph = temp,                  &
!'                      gam = pntVal,                 &
!'                      params = fpt%legToChebParams, &
!'                      strip_lb = iStrip,            &
!'                      strip_ub = strip_ub           )
!'   ! todo: fft on temp
!'   ! temp -> pntVal (stride-1 writing)
!'   pntVal(iStrip:strip_ub) = temp
!' end do

   if (.not. lobattoPoints) then

     ! Normalize the coefficients of the Chebyshev polynomials due
     ! to the unnormalized version of DCT-III in the FFTW.
     !$OMP DO
     do funcIndex = 1, n_cubed
       iFuncZ = (funcIndex-1)/n_squared + 1
       iDof = funcIndex - (iFuncZ-1)*n_squared
       iFuncY = (iDof-1)/fpt%legToChebParams%n + 1
       iFuncX = mod(iDof-1, fpt%legToChebParams%n) + 1
       legCoeffs(funcIndex) = normFactor(mod(iFuncX+iFuncY+iFuncZ,2)) &
         &                    * pntVal(funcIndex)
     end do
     !$OMP END DO
     !$OMP DO
     do FuncIndex = 1, n_squared
       legCoeffs(funcIndex) = legCoeffs(funcIndex)*2.0_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFunc = 1, n_squared
       iFuncZ = (iFunc-1)/fpt%legToChebParams%n + 1
       iFuncX = mod(iFunc-1, fpt%legToChebParams%n) + 1
       funcIndex = iFuncX + (iFuncZ-1) * n_squared
       legCoeffs(funcIndex) = legCoeffs(funcIndex)*2.0_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFunc = 1, n_squared
       funcIndex = 1 + (iFunc-1)*fpt%legToChebParams%n
       legCoeffs(funcIndex) = legCoeffs(funcIndex)*2.0_rk
     end do
     !$OMP END DO
  
     ! Transform Chebyshev expansion to point values at Chebyshev nodes
     ! by DCT III
     !$OMP SINGLE
     call fftw_execute_r2r( fpt%planChebToPnt, legCoeffs(:), pntVal(:) )
     !$OMP END SINGLE
   
   else
     legCoeffs(1) = pntVal(1)
     iFunc = fpt%legToChebParams%n
     legCoeffs(iFunc) = pntVal(iFunc)

     iFunc = 1 + (fpt%legToChebParams%n-1)*fpt%legToChebParams%n
     legCoeffs(iFunc) = pntVal(iFunc)

     iFunc = n_squared
     legCoeffs(iFunc) = pntVal(iFunc)

     iFunc = 1 + (fpt%legToChebParams%n-1)*n_squared
     legCoeffs(iFunc) = pntVal(iFunc)
     iFunc = fpt%legToChebParams%n + (fpt%legToChebParams%n-1)*n_squared
     legCoeffs(iFunc) = pntVal(iFunc)

     iFunc = 1 + (fpt%legToChebParams%n-1)*(fpt%legToChebParams%n &
       &                                    + n_squared)
     legCoeffs(iFunc) = pntVal(iFunc)

     iFunc = n_cubed
     legCoeffs(iFunc) = pntVal(iFunc)

     ! Normalize the coefficients of the Chebyshev polynomials 
     !$OMP DO
     do iFunc = 1, (fpt%legToChebParams%n-2)**3
       iFuncZ = (iFunc-1)/((fpt%legToChebParams%n-2)**2) + 2
       iDof = iFunc - (iFuncZ-2) * (fpt%legToChebParams%n-2)**2
       iFuncY = (iDof-1)/(fpt%legToChebParams%n-2) + 2
       iFuncX = mod(iDof-1, fpt%legToChebParams%n-2) + 2
?? copy :: posOfModgCoeffQTens(iFuncX, iFuncY, iFuncZ, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.125_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFunc = 1, (fpt%legToChebParams%n-2)**2
       iFuncY = (iFunc-1)/(fpt%legToChebParams%n-2) + 2
       iFuncX = mod(iFunc-1, fpt%legToChebParams%n-2) + 2
       ! first combination ...
?? copy :: posOfModgCoeffQTens(iFuncX, iFuncY, 1, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.25_rk
?? copy :: posOfModgCoeffQTens(iFuncX, iFuncY, fpt%legToChebParams%n, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.25_rk
       ! second combination ...
?? copy :: posOfModgCoeffQTens(iFuncX, 1, iFuncY, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.25_rk
?? copy :: posOfModgCoeffQTens(iFuncX, fpt%legToChebParams%n, iFuncY, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.25_rk
       ! third combination ...
?? copy :: posOfModgCoeffQTens(1, iFuncX, iFuncY, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.25_rk
?? copy :: posOfModgCoeffQTens(fpt%legToChebParams%n, iFuncX, iFuncY, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.25_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFuncX = 2, fpt%legToChebParams%n-1
       ! first combination ...
?? copy :: posOfModgCoeffQTens(iFuncX, 1, 1, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(iFuncX, 1, fpt%legToChebParams%n, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(iFuncX, fpt%legToChebParams%n, 1, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(iFuncX, fpt%legToChebParams%n, fpt%legToChebParams%n, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
       ! second combination ...
?? copy :: posOfModgCoeffQTens(1, iFuncX, 1, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(1, iFuncX, fpt%legToChebParams%n, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(fpt%legToChebParams%n, iFuncX, 1, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(fpt%legToChebParams%n, iFuncX, fpt%legToChebParams%n, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
       ! third combination ...
?? copy :: posOfModgCoeffQTens(1, 1, iFuncX, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(1, fpt%legToChebParams%n, iFuncX, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(fpt%legToChebParams%n, 1, iFuncX, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
?? copy :: posOfModgCoeffQTens(fpt%legToChebParams%n, fpt%legToChebParams%n, iFuncX, fpt%legToChebParams%n-1, funcIndex)
       legCoeffs(funcIndex) = pntVal(funcIndex) * 0.5_rk
     end do
     !$OMP END DO

     ! Transform Chebyshev expansion to point values at Chebyshev nodes by DCT I
     !$OMP SINGLE
!\     call fftw_execute_r2r( fpt%planChebToPnt, legCoeffs, pntVal )
     !$OMP END SINGLE
 
   end if

  end subroutine ply_legToPnt_3D_singVar
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  !!VK: no multivar fashion of this routine is used anymore
  subroutine ply_legToPnt_3D_multVar(fpt,legCoeffs,pntVal,nVars,lobattoPoints )
   !---------------------------------------------------------------------------
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified 
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:,:) 
   type(ply_legFpt_type), intent(inout) :: fpt
   real(kind=rk), intent(inout) :: pntVal(:,:)
   integer, intent(in) :: nVars
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iVar
   !---------------------------------------------------------------------------

   do iVar = 1, nVars
    call ply_legToPnt_3D(fpt, legCoeffs(:,iVar), pntVal(:,iVar), lobattoPoints) 
   end do

  end subroutine ply_legToPnt_3D_multVar
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_3D_multVar(fpt,pntVal,legCoeffs,nVars,lobattoPoints )
   !---------------------------------------------------------------------------
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The point values to transform to 3D modal Legendre expansion.
   !! \attention Although this array serves as input only, it is modified 
   !! inside of this routine by the underlying DCT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: pntVal(:,:)
   real(kind=rk), intent(inout) :: legCoeffs(:,:) 
   integer, intent(in) :: nVars
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iVar
   !---------------------------------------------------------------------------

   do iVar = 1, nVars
     call ply_pntToLeg_3D(fpt, pntVal(:,iVar), legCoeffs(:,iVar), lobattoPoints)
   end do

  end subroutine ply_pntToLeg_3D_multVar
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_3D_singVar( fpt, pntVal, legCoeffs, lobattoPoints )
   !---------------------------------------------------------------------------
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The point values to transform to 3D modal Legendre expansion.
   !! \attention Although this array serves as input only, it is modified 
   !! inside of this routine by the underlying DCT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: pntVal(:)
   real(kind=rk), intent(inout) :: legCoeffs(:) 
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iFunc, iDof, iFuncX, iFuncY, iFuncZ, funcIndex
   integer :: n_squared
   integer :: n_cubed
   real(kind=rk) :: normFactor(0:1)
   !---------------------------------------------------------------------------

   n_squared = fpt%chebToLegParams%n**2
   n_cubed = n_squared * fpt%chebToLegParams%n

   if ( .not. lobattoPoints ) then
     ! Transform the point values to Chebyshev polynomials by DCT II
     !$OMP SINGLE
     call fftw_execute_r2r( fpt%planPntToCheb, pntVal(:), legCoeffs(:) )
     !$OMP END SINGLE

     normfactor(1) = 1.0_rk / real(n_cubed, kind=rk)
     normfactor(0) = - normfactor(1)

     ! Apply normalization factors of the DCT II
     !$OMP DO
     do iFunc = 1, n_cubed
       iFuncX = (iFunc-1)/n_squared + 1
       iDof = iFunc - (iFuncX-1)*n_squared
       iFuncY = (iDof-1)/fpt%chebToLegParams%n + 1
       iFuncZ = iDof - (iFuncY-1)*(fpt%chebToLegParams%n)
       pntVal(iFunc) = legCoeffs(iFunc) &
         &             * normFactor(mod(iFuncX+iFuncY+iFuncZ,2))
     end do
     !$OMP END DO
     ! ... take care of the normalization factors for the constants per spatial direction.
     !$OMP DO
     do iFunc = 1, n_squared
       iFuncY = (iFunc-1)/fpt%chebToLegParams%n + 1
       iFuncZ = iFunc - (iFuncY-1)*fpt%chebToLegParams%n
?? copy :: posOfModgCoeffQTens(1,iFuncY,iFuncZ,fpt%chebToLegParams%n-1, funcIndex)
       pntVal(funcIndex) = pntVal(funcIndex) * 0.5_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFunc = 1, n_squared
       iFuncX = (iFunc-1)/fpt%chebToLegParams%n + 1
       iFuncZ = iFunc - (iFuncX-1)*fpt%chebToLegParams%n
?? copy :: posOfModgCoeffQTens(iFuncX,1,iFuncZ,fpt%chebToLegParams%n-1, funcIndex)
       pntVal(funcIndex) = pntVal(funcIndex) * 0.5_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFunc = 1, n_squared
       iFuncX = (iFunc-1)/fpt%chebToLegParams%n + 1
       iFuncY = iFunc - (iFuncX-1)*fpt%chebToLegParams%n
?? copy :: posOfModgCoeffQTens(iFuncX,iFuncY,1,fpt%chebToLegParams%n-1, funcIndex)
       pntVal(funcIndex) = pntVal(funcIndex) * 0.5_rk
     end do
     !$OMP END DO

   else

     ! Transform the point values to Chebyshev polynomials by DCT I
     !$OMP SINGLE
     call fftw_execute_r2r( fpt%planPntToCheb, pntVal(:), legCoeffs(:) )
     !$OMP END SINGLE

     ! Apply normalization factors
     !$OMP DO
     do iFunc = 1,(fpt%chebToLegParams%n-2)**3
       iFuncX = (iFunc-1)/((fpt%chebToLegParams%n-3+1)**2) + 2
       iDof = iFunc - (iFuncX-2)*((fpt%chebToLegParams%n-3+1)**2)
       iFuncY = (iDof-1)/(fpt%chebToLegParams%n-3+1) + 2
       iFuncZ = mod(iFunc-1,fpt%chebToLegParams%n-3+1) + 2
?? copy :: posOfModgCoeffQTens(iFuncX, iFuncY, iFuncZ, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 8.0_rk
     end do
     !$OMP END DO
     !$OMP DO
     do iFunc = 1,(fpt%chebToLegParams%n-2)**2
       iFuncX = (iFunc-1)/(fpt%chebToLegParams%n-3+1) + 2
       iFuncY = mod(iFunc-1,fpt%chebToLegParams%n-3+1) + 2
       ! ... first combination
?? copy :: posOfModgCoeffQTens(iFuncX, iFuncY, 1, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 4.0_rk
?? copy :: posOfModgCoeffQTens(iFuncX, iFuncY, fpt%chebToLegParams%n, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 4.0_rk
       ! ... second combination
?? copy :: posOfModgCoeffQTens(iFuncX, 1, iFuncY, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 4.0_rk
?? copy :: posOfModgCoeffQTens(iFuncX, fpt%chebToLegParams%n, iFuncY, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 4.0_rk
       ! ... third combination
?? copy :: posOfModgCoeffQTens(1, iFuncX, iFuncY, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 4.0_rk
?? copy :: posOfModgCoeffQTens(fpt%chebToLegParams%n, iFuncX, iFuncY, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 4.0_rk
     end do
     !$OMP END DO
     !$OMP DO 
     do iFuncX = 2, fpt%chebToLegParams%n-1
       ! ... first combination
?? copy :: posOfModgCoeffQTens(iFuncX, 1, 1, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(iFuncX, 1, fpt%chebToLegParams%n, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(iFuncX, fpt%chebToLegParams%n, 1, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(iFuncX, fpt%chebToLegParams%n, fpt%chebToLegParams%n, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
       ! ... second combination
?? copy :: posOfModgCoeffQTens(1, iFuncX, 1, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(1, iFuncX, fpt%chebToLegParams%n, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(fpt%chebToLegParams%n, iFuncX, 1, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(fpt%chebToLegParams%n, iFuncX, fpt%chebToLegParams%n,  fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
       ! ... third combination
?? copy :: posOfModgCoeffQTens(1, 1, iFuncX, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(1, fpt%chebToLegParams%n, iFuncX, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(fpt%chebToLegParams%n, 1, iFuncX, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
?? copy :: posOfModgCoeffQTens(fpt%chebToLegParams%n, fpt%chebToLegParams%n, iFuncX, fpt%chebToLegParams%n-1, funcIndex)
       legCoeffs(funcIndex) = legCoeffs(funcIndex) * 2.0_rk
     end do
     !$OMP END DO

     normfactor(0) = 1.0_rk / ((2.0_rk*(fpt%chebToLegParams%n-1))**3)

     !$OMP WORKSHARE
     pntVal(:) = legCoeffs(:) * normfactor(0)
     !$OMP END WORKSHARE

   end if

   ! Dimension-by-dimension transform Chebyshev polynomials to Legendre polynomial
   ! ... transformation in X direction (Cheb->Leg)
   call ply_fpt_exec_striped(nIndeps = n_squared,          &
     &                       alph    = pntVal,             &
     &                       gam     = legCoeffs,          &
     &                       params  = fpt%chebToLegParams )
   ! ... transformation in Y direction (Cheb->Leg)
   call ply_fpt_exec_striped(nIndeps = n_squared,          &
     &                       alph    = legCoeffs,          &
     &                       gam     = pntVal,             &
     &                       params  = fpt%chebToLegParams )
   ! ... transformation in Z direction (Cheb->Leg)
   call ply_fpt_exec_striped(nIndeps = n_squared,          &
     &                       alph    = pntVal,             &
     &                       gam     = legCoeffs,          &
     &                       params  = fpt%chebToLegParams )

  end subroutine ply_pntToLeg_3D_singVar

end module ply_legFpt_3D_module

