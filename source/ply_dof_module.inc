! Copyright (c) 2014,2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt and Harald Klimach for
! University of Siegen.
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !
!
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

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
?? text :: nextModgCoeffPTens(ansFuncX, ansFuncY, ansFuncZ)
  ! - ansatz indices are arranged in layers. Within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. Within a block, ansFuncZ is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. Each item represents one particular
  !   combination of ansFuncX, -Y, and -Z degrees.

  if (?ansFuncX? .ne. 1) then
    ! next item
    ?ansFuncX? = ?ansFuncX? - 1
    ?ansFuncY? = ?ansFuncY? + 1
  elseif (?ansFuncY? .ne. 1) then
    ! next block
    ?ansFuncX? = ?ansFuncY? - 1
    ?ansFuncY? = 1
    ?ansFuncZ? = ?ansFuncZ? + 1
  else
    ! next layer
    ?ansFuncX? = ?ansFuncZ? + 1
    ?ansFuncY? = 1
    ?ansFuncZ? = 1
  end if
?? end text nextModgCoeffPTens

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
?? text :: nextModgCoeffPTens2D(ansFuncX, ansFuncY)
  ! - ansatz indices are arranged in layers. Within each layer, the total
  !   degree remains constant.
  ! - within each layer, we have blocks. Within a block, ansFuncZ is
  !   constant, y counts up and x accordingly down.
  ! - within each block, we have items. Each item represents one particular
  !   combination of ansFuncX, -Y, and -Z degrees.
  if (?ansFuncX? .ne. 1) then
    ! next item
    ?ansFuncX? = ?ansFuncX? - 1
    ?ansFuncY? = ?ansFuncY? + 1
  else
    ! next layer
    ?ansFuncX? = ?ansFuncY? + 1
    ?ansFuncY? = 1
  end if
?? end text nextModgCoeffPTens2D

! The x ansatz degree is turned into the degree of the next
! ansatz function in the layered P list
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
?? text :: nextModgCoeffPTens1D(ansFuncX)
  ?ansFuncX? = ?ansFuncX? + 1
?? end text nextModgCoeffPTens1D

! The x, y and z ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Ansatz function index in z direction. First ansatz function has index 1.
! * Maximal polynomial degree
?? text :: nextModgCoeffQTens(ansFuncX, ansFuncY, ansFuncZ, maxdegree)
  if (?ansFuncX? .ne. (?maxdegree? + 1)) then
    ! next X index
    ?ansFuncX? = ?ansFuncX? + 1
  elseif (?ansFuncY? .ne. (?maxdegree? + 1)) then
    ! next y index
    ?ansFuncX? = 1
    ?ansFuncY? = ?ansFuncY? + 1
  else
    ! next Z index
    ?ansFuncX? = 1
    ?ansFuncY? = 1
    ?ansFuncZ? = ?ansFuncZ? + 1
  end if
?? end text nextModgCoeffQTens

! The x and y ansatz degrees are turned into the degrees of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
! * Ansatz function index in y direction. First ansatz function has index 1.
! * Maximal polynomial degree
?? text :: nextModgCoeffQTens2D(ansFuncX, ansFuncY, maxdegree)
  if (?ansFuncX? .ne. (?maxdegree? + 1)) then
    ! next X index
    ?ansFuncX? = ?ansFuncX? + 1
  else
    ! next y index
    ?ansFuncX? = 1
    ?ansFuncY? = ?ansFuncY? + 1
  end if
?? end text nextModgCoeffQTens2D

! The x ansatz degree is turned into the degree of the next
! ansatz function in the linearized Q tensor
! You must provide:
! * Ansatz function index in x direction. First ansatz function has index 1.
?? text :: nextModgCoeffQTens1D(ansFuncX)
  ?ansFuncX? = ?ansFuncX? + 1
?? end text nextModgCoeffQTens1D
