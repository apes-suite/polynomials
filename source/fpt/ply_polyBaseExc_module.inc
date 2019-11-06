! Copyright (c) 2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Langhammer Kay <kay.langhammer@student.uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt, Harald Klimach and
! Kay Langhammer for University of Siegen.
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

?? text :: EightTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                  &
              & + matrix(m, iDiag)     * alph(iFun + col_off)      &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off +  2) &
              & + matrix(m, iDiag + 2) * alph(iFun + col_off +  4) &
              & + matrix(m, iDiag + 3) * alph(iFun + col_off +  6) &
              & + matrix(m, iDiag + 4) * alph(iFun + col_off +  8) &
              & + matrix(m, iDiag + 5) * alph(iFun + col_off + 10) &
              & + matrix(m, iDiag + 6) * alph(iFun + col_off + 12) &
              & + matrix(m, iDiag + 7) * alph(iFun + col_off + 14)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: SevenTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                  &
              & + matrix(m, iDiag)     * alph(iFun + col_off)      &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off +  2) &
              & + matrix(m, iDiag + 2) * alph(iFun + col_off +  4) &
              & + matrix(m, iDiag + 3) * alph(iFun + col_off +  6) &
              & + matrix(m, iDiag + 4) * alph(iFun + col_off +  8) &
              & + matrix(m, iDiag + 5) * alph(iFun + col_off + 10) &
              & + matrix(m, iDiag + 6) * alph(iFun + col_off + 12)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: SixTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                  &
              & + matrix(m, iDiag)     * alph(iFun + col_off)      &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off +  2) &
              & + matrix(m, iDiag + 2) * alph(iFun + col_off +  4) &
              & + matrix(m, iDiag + 3) * alph(iFun + col_off +  6) &
              & + matrix(m, iDiag + 4) * alph(iFun + col_off +  8) &
              & + matrix(m, iDiag + 5) * alph(iFun + col_off + 10)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: FiveTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                 &
              & + matrix(m, iDiag)     * alph(iFun + col_off)     &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off + 2) &
              & + matrix(m, iDiag + 2) * alph(iFun + col_off + 4) &
              & + matrix(m, iDiag + 3) * alph(iFun + col_off + 6) &
              & + matrix(m, iDiag + 4) * alph(iFun + col_off + 8)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: FourTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                 &
              & + matrix(m, iDiag)     * alph(iFun + col_off)     &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off + 2) &
              & + matrix(m, iDiag + 2) * alph(iFun + col_off + 4) &
              & + matrix(m, iDiag + 3) * alph(iFun + col_off + 6)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: ThreeTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                 &
              & + matrix(m, iDiag)     * alph(iFun + col_off)     &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off + 2) &
              & + matrix(m, iDiag + 2) * alph(iFun + col_off + 4)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: TwoTimesUnrolledLoop(m_lowerbound, m_upperbound, m_next)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal)                                 &
              & + matrix(m, iDiag)     * alph(iFun + col_off)     &
              & + matrix(m, iDiag + 1) * alph(iFun + col_off + 2)
          end do
        end do
      end do
      ?m_next? = m
?? end text

?? text :: BaseLoop(m_lowerbound, m_upperbound)
      do m_blocking = ?m_lowerbound?, ?m_upperbound?, subblockingWidth
        do m = m_blocking,                                       &
          & min(m_blocking + subblockingWidth - 1, ?m_upperbound?)
          row_off = m + block_offset
          col_off = row_off + remainder + diag_off

          !NEC$ ivdep
          do indep = strip_lb + 1, strip_ub
            iFun = (indep - 1) * n
            iVal = iFun + row_off
            gam(iVal) = gam(iVal) + matrix(m, iDiag) * alph(iFun + col_off)
          end do
        end do
      end do
?? end text
