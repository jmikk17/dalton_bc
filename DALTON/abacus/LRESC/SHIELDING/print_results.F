
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
      SUBROUTINE PRINT_LRSHIELD
C Purpose:
C     Print LRESC-SHIELDING results section
C Author:
C     J. J. Aucar 2022
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"

      CALL QENTER('PRINT_LRSHIELD')

C     Paramagnetic 1st singlet or triplet
C
      IF ((SIGMAP1S).OR.(SIGMAP1T)) THEN
         CALL HEADER('Paramagnetic Corrections',-1)
         IF(SIGMAP1S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA1S]'
               WRITE(LUPRI,*)'       ---------------            xx  '//
     &         '          yy             zz          iso'
               call PRINT_LRSHIELD_CORR('Lkin    ',LRLKIN)
               call PRINT_LRSHIELD_CORR('PsoKin  ',LRPSOK)
            ENDIF
         ENDIF

         IF (SIGMAP1T) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA1T]'
               WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
               call PRINT_LRSHIELD_CORR('FcKin   ',LRFCZK)
               call PRINT_LRSHIELD_CORR('SdKin   ',LRSDZK)
               call PRINT_LRSHIELD_CORR('FcBso   ',LRFCBS)
               call PRINT_LRSHIELD_CORR('SdBso   ',LRSDBS)
            ENDIF
         ENDIF

         IF (SIGMAP3S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA3S]'
               WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
               call PRINT_LRSHIELD_CORR('PSOOZMv ',QRPSOOZMV)
               call PRINT_LRSHIELD_CORR('PSOOZDw ',QRPSOOZDW)
            ENDIF
         ENDIF

         IF (SIGMAP3T) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA3T]'
               WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
               call PRINT_LRSHIELD_CORR('FCOZSO  ',QROZFCSO)
               call PRINT_LRSHIELD_CORR('SDOZSO  ',QROZSDSO)
            ENDIF
         ENDIF

         WRITE(LUPRI,'(/,7X,A)') '---------------             xx  '//
     &         '          yy             zz          iso'
         IF (SIGMAP1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                      SGP1S)
         IF (SIGMAP1T) call PRINT_LRSHIELD_SUM('2nd Order Triplets:  ',
     &                      SGP1T)
         IF (SIGMAP3S) call PRINT_LRSHIELD_SUM('3nd Order Singlets:  ',
     &                      SGP3S)
         IF (SIGMAP3T) call PRINT_LRSHIELD_SUM('3nd Order Triplets:  ',
     &                      SGP3T)
      ENDIF
C
C     diamagnetic 0th or 1st order singlet
C
      IF ((SIGMAD0S).OR.(SIGMAD1S)) THEN
         WRITE(LUPRI,*)
         CALL HEADER(' Diamagnetic Corrections :',-1)
         IF(SIGMAD0S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.DIAM0S]'
               WRITE(LUPRI,*)'       ---------------            xx  '//
     &         '          yy             zz          iso'
               call PRINT_LRSHIELD_CORR('FC      ',LRFCAV)
               call PRINT_LRSHIELD_CORR('AngPso  ',LRANGP)
               call PRINT_LRSHIELD_CORR('DiaKin  ',LRDIAK)

            ENDIF
         ENDIF

         IF(SIGMAD1S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.DIAM1S]'
               WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
               call PRINT_LRSHIELD_CORR('DiaMv   ',LRDIAM)
               call PRINT_LRSHIELD_CORR('DiaDw   ',LRDIAD)
            ENDIF
         ENDIF

         WRITE(LUPRI,'(/,A)')'       ---------------             xx  '//
     &         '          yy             zz          iso'
         IF (SIGMAD0S) call PRINT_LRSHIELD_SUM('1st Order Singlets:  ',
     &                                         SGD0S)
         IF (SIGMAD1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                                         SGD1S)
      ENDIF
C
C     All of them in detail
C
      IF (LRESCALL) THEN
         CALL HEADER('Total Relativistic Corrections :',-1)
         WRITE(LUPRI,'(A)')'                                   xx  '//
     &         '          yy             zz          iso'
         IF((SIGMAP1S.AND.SIGMAP1T).OR.LRESCALL)
     &      call PRINT_LRSHIELD_TOT('Paramagnetic :',
     &      SGP1S+SGP1T+SGP3S+SGP3T)
         IF((SIGMAD1S.AND.SIGMAD0S).OR.LRESCALL)
     &      call PRINT_LRSHIELD_TOT('Diamagnetic  :',SGD0S+SGD1S)
      call PRINT_LRSHIELD_TOT('Total LRESC  :',SLRESC)
      END IF

      WRITE(LUPRI,*) 

      CALL QEXIT('PRINT_LRSHIELD')
      RETURN
      END

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
