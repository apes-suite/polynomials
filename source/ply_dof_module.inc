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
