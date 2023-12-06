!
!  Dalton, a molecular electronic structure program
!  Copyright (C) The Dalton Authors (see AUTHORS file for details).
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!

C     Search certain word in a specific file.
C     The file must be open before calling to this subroutine.
C     J. Aucar 2023
      subroutine SearchWord(file_unit, word_to_search, found)
#include "priunit.h"
      integer, intent(in) :: file_unit
      character(len=*), intent(in) :: word_to_search
      logical, intent(out) :: found
      character(len=100) :: line
      integer :: i

!     Initialize the control variable
      found = .false.

!     Read line by line from the file
      do
         read(file_unit, '(A)', iostat=i) line
c      write(lupri,*) "AUCAR", word_to_search
         if (i /= 0) exit ! End of file
!     Search for the word in the line
         if (index(line, word_to_search) /= 0) then
            found = .true.
            exit
         end if
      end do

      end subroutine SearchWord

