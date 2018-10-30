!> Managing the variable system description for sampled data.
module ply_sampling_varsys_module
  use env_module, only: rk

  use treelmesh_module, only: treelmesh_type
  use tem_time_module, only: tem_time_type
  use tem_tracking_module, only: tem_tracking_instance_type
  use tem_varsys_module, only: tem_varsys_type, tem_varSys_init

  implicit none

  private

  public :: ply_sampling_var_type
  public :: ply_sampling_varsys_for_track
  public :: ply_sampling_var_allocate
  public :: ply_sampling_var_move
  public :: ply_sampling_var_compute_elemdev

  !> Small helping type to allow arrays of arrays for the variable data.
  type ply_sampling_var_type
    integer :: nDeviating
    integer, allocatable :: degree(:)
    integer, allocatable :: first(:)
    real(kind=rk), pointer :: dat(:)
    logical, allocatable :: deviates(:)
  end type ply_sampling_var_type


contains


  ! ------------------------------------------------------------------------ !
  !> Create a variable system for the given tracking instance.
  subroutine ply_sampling_varsys_for_track( varsys, trackInst, mesh, nDims, &
    &                                       var_degree, sample_varsys, var, &
    &                                       time  )
    ! -------------------------------------------------------------------- !
    !> Variable system describing the access to the original data to sample.
    type(tem_varsys_type), intent(in) :: varsys

    !> The tracking object that should be sampled.
    type(tem_tracking_instance_type), intent(in) :: trackInst

    !> Original mesh describing the spatial organisation of the data to
    !! sample.
    type(treelmesh_type), intent(in) :: mesh

    !> Dimensionality of the data to sample.
    integer, intent(in) :: nDims

    !> Maximal polynomial degree for each variable.
    !!
    !! Needs to be matching the variable definition in the variable system.
    !! @todo Needs to be changed to be an information per element per variable!
    !!       Possibly by defining a variable in the varsys, providing the
    !!       degree.
    integer, intent(in) :: var_degree(:)

    !> Variable system for the sampled data.
    type(tem_varsys_type), intent(out) :: sample_varsys

    !> Extracted data for all the variables requested in the given tracking
    !! instance.
    type(ply_sampling_var_type), pointer :: var(:)

    !> Point in time to get the data for.
    type(tem_time_type), intent(in) :: time
    ! -------------------------------------------------------------------- !
    integer :: nVars
    integer :: nElems
    integer :: nScalars
    integer :: varpos
    integer :: nComponents
    integer :: nDofs
    integer :: i
    integer :: iElem
    integer :: iVar
    integer :: iComponent
    integer :: iScalar
    integer :: iTotComp
    integer :: lastScalar
    integer, allocatable :: elempos(:)
    real(kind=rk), allocatable :: elemdat(:)
    ! -------------------------------------------------------------------- !

    nVars = trackInst%varmap%varPos%nVals

    nullify(var)

    nScalars = sum( varsys%method%val(trackInst%varmap%varPos%val(:nVars)) &
      &                   %nComponents                                     )

    call tem_varSys_init( me         = sample_varsys, &
      &                   systemName = 'sampledVars', &
      &                   length     = nVars          )

    if (trackInst%subtree%useGlobalMesh) then
      nElems = mesh%nElems
      allocate( elempos(nElems) )
      elempos = [ (i, i=1,nElems) ]
    else
      nElems = trackInst%subtree%nElems
      allocate( elempos(nElems) )
      elempos = trackInst%subtree%map2global
    end if

    allocate(var(nScalars))

    iScalar = 1
    variables: do iVar=1,nVars

      varpos = trackInst%varmap%varPos%val(iVar)
      nComponents = varsys%method%val(varpos)%nComponents
      nDofs = (var_degree(iVar)+1)**nDims

      do iComponent=1,nComponents
        iTotComp = iScalar+iComponent-1
        call ply_sampling_var_allocate( var     = var(iTotComp), &
          &                             nElems  = nElems,        &
          &                             datalen = nDofs*nElems   )
      end do

      isScalar: if (nComponents == 1) then

        ! If there is only one component, we can directly copy the data
        ! into vardat.
        call varSys%method%val(varpos)%get_element( &
          &    varSys  = varSys,                    &
          &    elempos = elempos,                   &
          &    time    = time,                      &
          &    tree    = mesh,                      &
          &    nElems  = nElems,                    &
          &    nDofs   = nDofs,                     &
          &    res     = var(iScalar)%dat           )

      else isScalar

        ! If the variable is not scalar, we need to take care of the
        ! component data and separate them for `var%dat`.
        ! To avoid overly large memory consumption, we do this element by
        ! element.
        allocate(elemDat(nComponents*ndofs))

        do iElem=1,nElems
          call varSys%method%val(varpos)%get_element( &
            &    varSys  = varSys,                    &
            &    elempos = elempos(iElem:iElem),      &
            &    time    = time,                      &
            &    tree    = mesh,                      &
            &    nElems  = 1,                         &
            &    nDofs   = ndofs,                     &
            &    res     = elemdat                    )
          do iComponent=1,nComponents
            iTotComp = iScalar+iComponent-1
            var(iTotComp)%dat(1+nDofs*(iElem-1):nDofs*iElem) &
              & = elemdat(iComponent::nComponents)
          end do
        end do

        deallocate(elemDat)

      end if isScalar

      ! For each element we store the position of its first degree of freedom
      ! (the last is given by the first of the next element), and the
      ! polynomial degree in each element.
      ! This is redundant, but convenient.
      lastScalar = iScalar+nComponents-1
      do iComponent=0,nComponents-1
        var(iScalar+iComponent)%first(1) = 1
        do iElem=1,nElems
          var(iScalar+iComponent)%first(iElem+1)     &
            & = var(iScalar+iComponent)%first(iElem) + nDofs
          var(iScalar+iComponent)%degree(iElem) = var_degree(iVar)
        end do
      end do

      iScalar = iScalar + nComponents

    end do variables


  end subroutine ply_sampling_varsys_for_track
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Allocate memory for a sampled variable.
  subroutine ply_sampling_var_allocate(var, nElems, datalen)
    ! -------------------------------------------------------------------- !
    !> The variable to allocate the space for.
    type(ply_sampling_var_type), intent(inout) :: var

    !> Number of elements the data lives in.
    integer, intent(in) :: nElems

    !> Size of the container to use for representation of the polynomial
    !! data across all elements.
    integer, intent(in) :: datalen
    ! -------------------------------------------------------------------- !
    ! -------------------------------------------------------------------- !
    allocate(var%degree(nelems))
    allocate(var%deviates(nelems))
    allocate(var%first(nelems+1))
    allocate(var%dat(datalen))
  end subroutine ply_sampling_var_allocate
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> Move the variable data from source to destination.
  !!
  !! If there is data in destination, it will be discarded.
  !! The data in source becomes accessible via destination and source
  !! itself gets nullified.
  subroutine ply_sampling_var_move(source, destination)
    ! -------------------------------------------------------------------- !
    !> Variable data to move (and make accessible via destination).
    !!
    !! Source itself will be null after moving.
    type(ply_sampling_var_type), pointer :: source(:)

    !> Pointer to refer to the data in source. If destination already
    !! contains data, this data will be discarded.
    type(ply_sampling_var_type), pointer :: destination(:)
    ! -------------------------------------------------------------------- !
    integer :: nScalars
    integer :: iScalar
    ! -------------------------------------------------------------------- !

    if (associated(destination)) then
      nScalars = size(destination)
      do iScalar=1,nScalars
        if (allocated(destination(iScalar)%degree)) &
          & deallocate(destination(iScalar)%degree)
        if (allocated(destination(iScalar)%first)) &
          & deallocate(destination(iScalar)%first)
        if (associated(destination(iScalar)%dat)) &
          & deallocate(destination(iScalar)%dat)
        if (allocated(destination(iScalar)%deviates)) &
          & deallocate(destination(iScalar)%deviates)
      end do
    end if
    destination => source
    nullify(source)

  end subroutine ply_sampling_var_move
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> This routine computes for each element whether the solution in it is
  !! considered to be deviating from the mean above the given threshold or
  !! not. The logical result is stored in `var%deviates` for each element.
  !!
  !! The total number of deviating elements is stored in `var%nDeviating`.
  !!
  !! The variation is computed by the sum of the absolute values of all higher
  !! modes divided by the first mode.
  !! As we are using series of Legendre polynomials this also is a
  !! bounding estimation for the maximal (relative) deviation from the mean in
  !! this element.
  !!
  !! A variation of 0.01 for example would imply that the state in the element
  !! is guaranteed to nowhere deviate from the mean by more than 1 percent.
  !! However, this is a very rough estimation and the actual maximal deviation
  !! from the mean is probably much lower (at least for sufficiently high
  !! polynomial degrees).
  !!
  !! If the mean is too close to 0, we use epsilon for the normalization
  !! instead of the actual mean.
  !!
  !! The computation is done for the current data found in `var%dat`, any
  !! previous computations of these flags will be discarded by this routine.
  subroutine ply_sampling_var_compute_elemdev(var, threshold, min_mean)
    ! -------------------------------------------------------------------- !
    !> Variable data to compute the deviation for.
    type(ply_sampling_var_type), intent(inout) :: var

    !> Relative threshold to use as decision whether an element has a high
    !! deviation or not.
    !!
    !! If the absolute value of higher modes sums to a larger value than
    !! threshold times the first mode (integral mean), the element is marked
    !! as deviating.
    real(kind=rk), intent(in) :: threshold

    !> A minimal mean value to use as comparison (to cut off changes that are
    !! too close to 0).
    !!
    !! This should be small but has to be larger than 0.
    real(kind=rk), intent(in) :: min_mean
    ! -------------------------------------------------------------------- !
    real(kind=rk) :: variation
    real(kind=rk) :: absmean
    integer :: iElem
    integer :: nElems
    integer :: ndofs
    ! -------------------------------------------------------------------- !

    if (allocated(var%deviates)) deallocate(var%deviates)
    var%nDeviating = 0

    if (allocated(var%first)) then
      nElems = size(var%first)-1
      allocate(var%deviates(nElems))

      do iElem=1,nElems
        ndofs = var%first(iElem+1) - var%first(iElem) - 1
        absmean = max( abs(var%dat(var%first(iElem))), min_mean )
        variation = sum( abs(var%dat(var%first(iElem)+1:var%first(iElem+1)-1)) )
        var%deviates(iElem) = (variation > threshold*absmean)
        if (var%deviates(iElem)) var%nDeviating = var%nDeviating + 1
      end do
    end if


  end subroutine ply_sampling_var_compute_elemdev
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

end module ply_sampling_varsys_module
