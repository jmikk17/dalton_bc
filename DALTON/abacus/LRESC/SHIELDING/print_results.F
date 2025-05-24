
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
C

C     /* Deck print_lrshield */
      SUBROUTINE PRINT_LRSHIELD_CORR(names,tensor)
C Purpose:
C     Print information about LRESC-SHIELDING correction
C Author:
C     J. J. Aucar 2022
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"

      REAL*8, intent(in)       :: tensor(3,3)
      REAL*8                   :: iso
      CHARACTER*8, intent(in)  :: names

      CALL QENTER('PRINT_LRSHIELD_CORR')

C     Information about calculation of LRESC
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(9X,A,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &      names, ':',tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'

      CALL QEXIT('PRINT_LRSHIELD_CORR')
      RETURN
      END


C     /* Deck print_lrshield_SUM */
      SUBROUTINE PRINT_LRSHIELD_SUM(phrase,tensor)
C Purpose:
C     Print information about sum of LRESC-SHIELDING corrections
C Author:
C     J. J. Aucar 2022
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"

      REAL*8, intent(in)       :: tensor(3,3)
      REAL*8                   :: iso
      CHARACTER*21, intent(in)  :: phrase
c      character(len=:), allocatable :: phrase

      CALL QENTER('PRINT_LRSHIELD_SUM')

C     Information about calculation of LRESC
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(4X,A,2X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &      phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'

      CALL QEXIT('PRINT_LRSHIELD_SUM')
      RETURN
      END

C     /* Deck print_lrshield_TOT */
      SUBROUTINE PRINT_LRSHIELD_TOT(phrase,tensor)
C Purpose:
C     Print information about LRESC-SHIELDING overall contributions
C Author:
C     J. J. Aucar 2022
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"

      REAL*8, intent(in)        :: tensor(3,3)
      REAL*8                    :: iso
      CHARACTER*14, intent(in)  :: phrase

      CALL QENTER('PRINT_LRSHIELD_TOT')

C     Information about calculation of LRESC
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(5X,A,8X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &      phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'

      CALL QEXIT('PRINT_LRSHIELD_TOT')
      RETURN
      END
