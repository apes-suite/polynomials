! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of 
!   modal coefficients in.
?? text :: posOfModgCoeffQTens(ansFuncX, ansFuncY, ansFuncZ, maxdegree, result)
  ?result? = ?ansFuncX?                                      &
    &      + ( ( ?ansFuncY?-1)                             &
    &      + (?ansFuncZ?-1)*(?maxdegree?+1))*(?maxdegree?+1)
?? end text posOfModgCoeffQTens


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of 
!   modal coefficients in.
?? text :: posOfModgCoeffQTens2D(ansFuncX, ansFuncY, maxdegree, result)
  ?result? = ?ansFuncX? + (?ansFuncY?-1)*(?maxdegree?+1)
?? end text posOfModgCoeffQTens2D


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for Q-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of 
!   modal coefficients in.
?? text :: posOfModgCoeffQTens1D(ansFuncX, result)
  ?result? = ?ansFuncX?
?? end text posOfModgCoeffQTens1D


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * Ansatzfunction index in z direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of 
!   modal coefficients in.
?? text :: posOfModgCoeffPTens(ansFuncX, ansFuncY, ansFuncZ, result)
  ! integer divisions are no mistake here.
  ?result? = (((?ansFuncX? + ?ansFuncY? + ?ansFuncZ? - 3) &
    &     * (?ansFuncX? + ?ansFuncY? + ?ansFuncZ? - 2) &
    &     * (?ansFuncX? + ?ansFuncY? + ?ansFuncZ? - 1)) &
    &   / 6 + 1)             &
    & + ((?ansFuncZ?-1) * (?ansFuncX? + ?ansFuncY? + ?ansFuncZ? -2) &
    &   - ((?ansFuncZ?-2) * (?ansFuncZ?-1)) / 2) &
    & + (?ansFuncY?-1)
?? end text posOfModgCoeffPTens 


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * Ansatzfunction index in y direction. Index starts with 1.
! * The maximal polynomial degree per spatial direction.
! * The variable to store the position of the modal coefficient in the list of 
!   modal coefficients in.
?? text :: posOfModgCoeffPTens2D(ansFuncX, ansFuncY, result)
  ! integer divisions are no mistake here.
  ?result? = ((((?ansFuncX? - 1) + (?ansFuncY? - 1))            &
    &   * (((?ansFuncX? - 1) + (?ansFuncY? - 1)) + 1)) / 2 + 1) &
    & + (?ansFuncY? - 1)
?? end text posOfModgCoeffPTens2D


! Return the position of a given ansatz function combination in the
! linearized list of modal coefficients for P-Tensor product polynomials.
! You must provide
! * Ansatzfunction index in x direction. Index starts with 1.
! * The variable to store the position of the modal coefficient in the list of 
!   modal coefficients in.
?? text :: posOfModgCoeffPTens1D(ansFuncX, result)
  ! integer divisions are no mistake here.
  ?result? = ?ansFuncX?
?? end text posOfModgCoeffPTens1D


! Return the number of degrees of freedom for Q polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial
?? text :: getDofsQTens(maxPolyDegree, result)
  ?result? = ((?maxPolyDegree?)+1)**3
?? end text getDofsQTens


! Return the number of degrees of freedom for broken polynomial space
! Your must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial
?? text :: getDofsPTens(maxPolyDegree, result)
  ?result? = (((?maxPolyDegree?) + 1) &
    &   * ((?maxPolyDegree?) + 2) &
    &   * ((?maxPolyDegree?) + 3)) &
    & / 6
??  end text getDofsPTens


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial
?? text :: getDofsQTens2D(maxPolyDegree, result)
  ?result? = ((?maxPolyDegree?)+1)**2
?? end text getDofsQTens2D


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * A variable to store the number of degrees of freedom for a P tensor product
!   polynomial
?? text :: getDofsPTens2D(maxPolyDegree, result)
  ?result? = ((?maxPolyDegree?)+1)*((?maxPolyDegree?)+2)/2
?? end text getDofsPTens2D


! Return the number of degrees of freedom for Q polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction
! * The variable to store the number of degrees of freedom for a Q tensor
!   product polynomial
?? text :: getDofsQTens1D(maxPolyDegree, result)
  ?result? = ((?maxPolyDegree?)+1)
?? end text getDofsQTens1D


! Return the number of degrees of freedom for broken polynomial space
! You must provide:
! * The maximal polynomial degree per spatial direction (for P Tensor product
!   polynomials this assumed to be the same for each spatial direction).
! * The variable to store the number of degrees of freedom for a P tensor
!   product polynomial
?? text :: getDofsPTens1D(maxPolyDegree, result)
  ?result? = ((?maxPolyDegree?)+1)
?? end text getDofsPTens1D

