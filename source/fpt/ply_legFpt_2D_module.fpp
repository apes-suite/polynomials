?? include "ply_dof_module.inc"
!> Module providing datatypes and routines for a fast
!! transformation of Legendre expansion to point values.
!! \author{Jens Zudrop}
module ply_legFpt_2D_module
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
  use ply_nodes_module,       only: ply_faceNodes_type
  use fftw_wrap
  use ply_legFpt_module,      only: ply_legFpt_type, &
                                  & assignment(=)

  implicit none

  private

  interface ply_legToPnt_2D
    module procedure ply_legToPnt_2D_singVar
    module procedure ply_legToPnt_2D_multVar
  end interface ply_legToPnt_2D

  interface ply_pntToLeg_2D
    module procedure ply_pntToLeg_2D_singVar
    module procedure ply_pntToLeg_2D_multVar
  end interface ply_pntToLeg_2D

  public :: ply_init_legFpt_2D, ply_legToPnt_2D, ply_pntToLeg_2D


contains


  !****************************************************************************
  !> Subroutine to initialize the fast polynomial transformation
  !! for Legendre expansion.
  !! VK: change maxPolyDegree to number of quadpoints to include the dealising factor
  subroutine ply_init_legFpt_2D( maxPolyDegree, nVars, fpt, blocksize, &
    &                            striplen, lobattoPoints, subblockingWidth )
    !---------------------------------------------------------------------------
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: nVars
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The blocksize of the underlying fast Legendre to Chebyshev transformation.
    integer, intent(in), optional :: blocksize

    !> Length to use in vectorization, this is the number of independent
    !! matrix multiplications that are to be done simultaneously.
    !!
    !! Defaults to 512, but the optimal setting is platform specific.
    integer, optional, intent(in) :: striplen
    !> Use Chebyshev-Lobtatto Points (true) or simple Chebyshev points (false)
    logical, intent(in), optional :: lobattoPoints
    !> The width of the subblocks used during the unrolled base exchange to
    !! ensure a better cache usage.
    integer, optional, intent(in) :: subblockingWidth
    !---------------------------------------------------------------------------
    integer :: n
    real(kind=rk), allocatable :: tmpOut(:,:), tmpIn(:,:)
    logical :: lob
    !$ integer :: fftwMultThread
    !---------------------------------------------------------------------------

    lob = .false.
    if(present(lobattoPoints)) then
      if (lobattoPoints) then
        lob = .true.
      end if
    end if

    ! init the fast Legendre to Chebyshev transformation.
    call ply_fpt_init( n                = maxPolyDegree + 1,   &
      &                params           = fpt%legToChebParams, &
      &                trafo            = ply_legToCheb_param, &
      &                blocksize        = blocksize,           &
      &                striplen         = striplen,            &
      &                subblockingWidth = subblockingWidth     )
    call ply_fpt_init( n                = maxPolyDegree + 1,   &
      &                params           = fpt%chebToLegParams, &
      &                trafo            = ply_chebToLeg_param, &
      &                blocksize        = blocksize,           &
      &                striplen         = striplen,            &
      &                subblockingWidth = subblockingWidth     )

    ! Create the buffers for the intermediate arrays
!   n = fpt%legToChebParams%n**2
    n = fpt%legToChebParams%n

    ! Temprorary arrays to initialize FFTW real->real transformations
    allocate( tmpIn(n, nVars) )
    allocate( tmpOut(n, nVars) )

    ! If we have OpenMP parallelism, we have to init FFTW with the
    ! corrseponding function call
    !$ fftwMultThread = fftw_init_threads()
    !$ if(fftwMultThread.eq.0) then
    !$   write(logUnit(1),*) 'ERROR in ply_init_legFpt_2D: Not able to init OpenMP parallel FFTW, stopping...'
    !$   call tem_abort()
    !$ end if

    ! Tell FFTW how many threads at max we want to use: at max
    ! we use the number of OMP threads in the current team.
    !$ call fftw_plan_with_nthreads( omp_get_max_threads() )

    if(.not.lob) then
      ! Init the DCT III ( Leg -> Point values )
      fpt%planChebToPnt = fftw_plan_r2r_1d( n = fpt%legToChebParams%n,&
                                          & in = tmpIn, &!legCoeffs, &
                                          & out = tmpOut, &!fpt%pntVal, &
                                          & kind = FFTW_REDFT01, &
                                          & flags = FFTW_ESTIMATE )

      ! Init the DCT II ( Point values -> Leg )
      fpt%planPntToCheb = fftw_plan_r2r_1d( n = fpt%legToChebParams%n, &
                                          & in = tmpIn, & !fpt%pntVal,&
                                          & out = tmpOut, & !legCoeffs, &
                                          & kind = FFTW_REDFT10, &
                                          & flags = FFTW_ESTIMATE )

      deallocate( tmpIn, tmpOut )
    else
      ! Init the DCT I  (Leg -> nodal): To be used with a normalization factor for trafo ...
      fpt%planChebToPnt = fftw_plan_r2r_1d( n = fpt%legToChebParams%n,&
                                          & in = tmpIn, &!legCoeffs, &
                                          & out = tmpOut, &!fpt%pntVal, &
                                          & kind = FFTW_REDFT00, &
                                          & flags = FFTW_ESTIMATE )

      ! Init the DCT I  (nodal -> Leg): To be used with a normalization factor for trafo ...
      fpt%planPntToCheb = fftw_plan_r2r_1d( n = fpt%legToChebParams%n, &
                                          & in = tmpIn, & !fpt%pntVal,&
                                          & out = tmpOut, & !legCoeffs, &
                                          & kind = FFTW_REDFT00, &
                                          & flags = FFTW_ESTIMATE )

      deallocate( tmpIn, tmpOut )
    end if

  end subroutine ply_init_legFpt_2D
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_2D_singVar( fpt, legCoeffs, pntVal, lobattoPoints )
   !---------------------------------------------------------------------------
   !> The FPT parameters.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:)
   !> The resulting point values (Chebyshev nodes).
   real(kind=rk), intent(inout) :: pntVal(:)
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iFunc, iFuncX, iFuncY, funcIndex
   integer :: striplen, iStrip, n, iAlph, nIndeps
   integer :: iDof
   real(kind=rk), dimension(:), allocatable :: alph
   real(kind=rk), dimension(:), allocatable :: gam
   real(kind=rk) :: normFactor
   !---------------------------------------------------------------------------

   striplen = fpt%legToChebParams%striplen
   n = fpt%legToChebParams%n

   allocate(alph(min(striplen, n)*n))
   allocate(gam(min(striplen, n)*n))
   ! original layout (n = 3):
   !  1  2  3
   !  4  5  6
   !  7  8  9

   ! layout after y-trafo:
   !  1  4  7
   !  2  5  8
   !  3  6  9
write(*,*) '////////////////'
write(*,*) 'y-Direction'
   yStripLoop: do iStrip = 1, n, striplen
     ! iAlph is the index of the first element in a line for the transformation in 
     ! y-direction. 
     do iAlph = iStrip, min(iStrip+striplen-1, n)  !y-Trafo
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = &
           & legCoeffs(iAlph::n) !ytrafo
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n-iStrip+1)
 
     ! ply_fpt_exec on temp (no memory transpose)
     call ply_fpt_exec( alph = alph,                   &
       &                gam = gam,                     &
       &                nIndeps = nIndeps,             &
       &                plan = fpt%planChebToPnt,      &
       &                lobattoPoints = lobattoPoints, &
       &                params = fpt%legToChebParams   )

     ! Write gam to pntVal array
     pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)

   end do yStripLoop

  ! x-direction
   xStripLoop: do iStrip = 1, n, striplen
     do iAlph = iStrip, min(iStrip+striplen-1, n)  
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = &
           & pntVal(iAlph::n) !ztrafo
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n-iStrip+1)
write(*,*) '////////////////'
write(*,*) 'x-Direction'
     ! ply_fpt_exec on temp (no memory transpose)
     call ply_fpt_exec( alph = alph,                   &
       &                gam = gam,                     &
       &                nIndeps = nIndeps,             &
       &                plan = fpt%planChebToPnt,      &
       &                lobattoPoints = lobattoPoints, &
       &                params = fpt%legToChebParams   )

     pntVal((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)

   end do xStripLoop

  end subroutine ply_legToPnt_2D_singVar
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_legToPnt_2D_multVar(fpt,legCoeffs, pntVal, nVars,lobattoPoints )
   !---------------------------------------------------------------------------
   !> The FPT parameters.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The Legendre coefficients to convert to point values (Chebyshev nodes).
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying FPT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: legCoeffs(:,:)
   !> The resulting point values (Chebyshev nodes).
   real(kind=rk), intent(inout) :: pntVal(:,:)
   !> The number of scalar variables to transform.
   integer, intent(in) :: nVars
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iVar
   !---------------------------------------------------------------------------

   do iVar = 1, nVars
     call ply_legToPnt_2D(fpt, legCoeffs(:,iVar), pntVal(:,iVar), lobattoPoints)
   end do

  end subroutine ply_legToPnt_2D_multVar
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_2D_singVar( fpt, pntVal, legCoeffs, lobattoPoints )
    !---------------------------------------------------------------------------
    !> Parameters of the Fast Polynomial transformation.
    type(ply_legFpt_type), intent(inout) :: fpt
    !> The point values to transform to 2D modal Legendre expansion.
    !! \attention Although this array serves as input only, it is modified
    !! inside of this routine by the underlying DCT algorithm. So, when
    !! this routine returns from its call the original values of pntVal will
    !! be modified.
    real(kind=rk), intent(inout) :: pntVal(:)
    !> The Legendre coefficients.
    real(kind=rk), intent(inout) :: legCoeffs(:)
    logical, intent(in) :: lobattoPoints
    !---------------------------------------------------------------------------
    integer :: iFunc, iFuncX, iFuncY, funcIndex, iDof
    integer :: iStrip, striplen, nIndeps, iAlph, n, n_squared
    real(kind=rk), dimension(:), allocatable :: alph
    real(kind=rk), dimension(:), allocatable :: gam
    real(kind=rk) :: normFactor, inv_ndofs
    !---------------------------------------------------------------------------

   striplen = fpt%chebToLegParams%striplen
   n = fpt%legToChebParams%n
   n_squared = n**2

   allocate(alph(min(striplen, n)*n))
   allocate(gam(min(striplen, n)*n))

   yStripLoop: do iStrip = 1, n, striplen
     do iAlph = iStrip, min(iStrip+striplen-1, n)  !y_Trafo
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = &
           & pntVal(iAlph::n) 
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n-iStrip+1)

     ! ply_fpt_exec on temp (no memory transpose)
     call ply_fpt_exec( alph = alph,                   &
       &                gam = gam,                     &
       &                nIndeps = nIndeps,             &
       &                plan = fpt%planPntToCheb,      &
       &                lobattoPoints = lobattoPoints, &
       &                params = fpt%chebToLegParams   )
 
     ! temp -> pntVal (stride-1 writing)
     legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)


   end do yStripLoop ! iStrip

  ! x-direction
   xStripLoop: do iStrip = 1,n,striplen
     do iAlph = iStrip, min(iStrip+striplen-1, n) 
       alph((iAlph-iStrip)*n+1:(iAlph-iStrip+1)*n) = &
           & legCoeffs(iAlph::n) !ztrafo
     end do

     ! At the end of the array the number of computed strips might be smaller
     nIndeps = min(striplen, n-iStrip+1)

     ! ply_fpt_exec on temp (no memory transpose)
     call ply_fpt_exec( alph = alph,                   &
       &                gam = gam,                     &
       &                nIndeps = nIndeps,             &
       &                plan = fpt%planPntToCheb,      &
       &                lobattoPoints = lobattoPoints, &
       &                params = fpt%chebToLegParams   )

     legCoeffs((iStrip-1)*n+1 : (iStrip+nIndeps-1)*n)  = gam(1:nIndeps*n)

   end do xStripLoop

  end subroutine ply_pntToLeg_2D_singVar
  !****************************************************************************


  !****************************************************************************
  !> Subroutine to transform Legendre expansion to point values
  !! at Chebyshev nodes.
  subroutine ply_pntToLeg_2D_multVar(fpt,pntVal,legCoeffs,nVars,lobattoPoints )
   !---------------------------------------------------------------------------
   !> Parameters of the Fast Polynomial transformation.
   type(ply_legFpt_type), intent(inout) :: fpt
   !> The point values to transform to 2D modal Legendre expansion.
   !! \attention Although this array serves as input only, it is modified
   !! inside of this routine by the underlying DCT algorithm. So, when
   !! this routine returns from its call the original values of pntVal will
   !! be modified.
   real(kind=rk), intent(inout) :: pntVal(:,:)
   !> The Legendre coefficients.
   real(kind=rk), intent(inout) :: legCoeffs(:,:)
   !> The number of scalar variables to transform.
   integer, intent(in) :: nVars
   logical, intent(in) :: lobattoPoints
   !---------------------------------------------------------------------------
   integer :: iVar
   !---------------------------------------------------------------------------

   do iVar = 1, nVars
     call ply_pntToLeg_2D(fpt, pntVal(:,iVar), legCoeffs(:,iVar), lobattoPoints)
   end do

  end subroutine ply_pntToLeg_2D_multVar
  !****************************************************************************

end module ply_legFpt_2D_module

