
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
#include "mxcent.h"
#include "nuclei.h"

      CHARACTER*15  :: slash

      slash = '---------------'

      CALL QENTER('PRINT_LRSHIELD')

C     Paramagnetic 1st singlet or triplet
C
      IF ((SIGMAP1S).OR.(SIGMAP1T).OR.(SIGMAP3S).OR.(SIGMAP3T)) THEN
         IF (PRTALL1) THEN
         CALL HEADER('Paramagnetic Corrections for '
     &                //NAMN(LRATOM)//' :',-1)
         END IF
         IF(SIGMAP1S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA1S]'
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A,10X,A)') 
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('OZK     ',LROZK,'iso')
               call PRINT_LRSHIELD_CORR('PsoKin  ',LRPSOK,'iso')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A)') 
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('OZK     ',LROZK,'ani')
               call PRINT_LRSHIELD_CORR('PsoKin  ',LRPSOK,'ani')
               END IF
            ENDIF
         ENDIF

         IF (SIGMAP1T) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA1T]'
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A,10X,A)') 
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('FcKin   ',LRFCZK,'iso')
               call PRINT_LRSHIELD_CORR('FcDer   ',LRFCDE,'iso')
               call PRINT_LRSHIELD_CORR('SdKin   ',LRSDK,'iso')
c               call PRINT_LRSHIELD_CORR('SdKin-2 ',LRSDK2,'iso')
               call PRINT_LRSHIELD_CORR('SdDer   ',LRSDD,'iso')
c               call PRINT_LRSHIELD_CORR('SdDer-2 ',LRSDD2,'iso')
               call PRINT_LRSHIELD_CORR('FcBso   ',LRFCBS,'iso')
               call PRINT_LRSHIELD_CORR('SdBso   ',LRSDBS,'iso')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A)') 
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('FcKin   ',LRFCZK,'ani')
               call PRINT_LRSHIELD_CORR('FcDer   ',LRFCDE,'ani')
               call PRINT_LRSHIELD_CORR('SdKin   ',LRSDK,'ani')
c               call PRINT_LRSHIELD_CORR('SdKin-2 ',LRSDK2,'ani')
               call PRINT_LRSHIELD_CORR('SdDer   ',LRSDD,'ani')
c               call PRINT_LRSHIELD_CORR('SdDer-2 ',LRSDD2,'ani')
               call PRINT_LRSHIELD_CORR('FcBso   ',LRFCBS,'ani')
               call PRINT_LRSHIELD_CORR('SdBso   ',LRSDBS,'ani')
               END IF
            ENDIF
         ENDIF

         IF (SIGMAP3S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA3S]'
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A,10X,A)') 
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('PSOOZMv ',QRPSOOZMV,'iso')
               call PRINT_LRSHIELD_CORR('PSOOZDw ',QRPSOOZDW,'iso')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A)') 
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('PSOOZMv ',QRPSOOZMV,'ani')
               call PRINT_LRSHIELD_CORR('PSOOZDw ',QRPSOOZDW,'ani')
               END IF
            ENDIF
         ENDIF

         IF (SIGMAP3T) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA3T]'
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A,10X,A)') 
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('FCOZSO  ',QROZFCSO,'iso')
               call PRINT_LRSHIELD_CORR('SDOZSO  ',QROZSDSO,'iso')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A)') 
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('FCOZSO  ',QROZFCSO,'ani')
               call PRINT_LRSHIELD_CORR('SDOZSO  ',QROZSDSO,'ani')
               END IF
            ENDIF
         ENDIF

         IF (PRTALL1) THEN
         WRITE(LUPRI,'(/,7X,A,13X,A,12X,A,13X,A,10X,A)') 
     &    slash,'xx','yy','zz','iso'
         IF (SIGMAP1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                      SGP1S,'iso')
         IF (SIGMAP1T) call PRINT_LRSHIELD_SUM('2nd Order Triplets:  ',
     &                      SGP1T,'iso')
         IF (SIGMAP3S) call PRINT_LRSHIELD_SUM('3nd Order Singlets:  ',
     &                      SGP3S,'iso')
         IF (SIGMAP3T) call PRINT_LRSHIELD_SUM('3nd Order Triplets:  ',
     &                      SGP3T,'iso')
         IF (LRANISO) THEN
         WRITE(LUPRI,'(7X,A,13X,A,12X,A,13X,A)') 
     &    slash,'xy','xz','yz'
         IF (SIGMAP1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                      SGP1S,'ani')
         IF (SIGMAP1T) call PRINT_LRSHIELD_SUM('2nd Order Triplets:  ',
     &                      SGP1T,'ani')
         IF (SIGMAP3S) call PRINT_LRSHIELD_SUM('3nd Order Singlets:  ',
     &                      SGP3S,'ani')
         IF (SIGMAP3T) call PRINT_LRSHIELD_SUM('3nd Order Triplets:  ',
     &                      SGP3T,'ani')
         END IF
         END IF
      ENDIF
C
C     diamagnetic 0th or 1st order singlet
C
      IF ((SIGMAD0S).OR.(SIGMAD1S)) THEN
         IF (PRTALL1) THEN
            WRITE(LUPRI,*)
            CALL HEADER(' Diamagnetic Corrections for '
     &                //NAMN(LRATOM)//' :',-1)
         END IF
         IF(SIGMAD0S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.DIAM0S]'
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A,10X,A)') 
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('FC      ',LRFCAV,'iso')
               call PRINT_LRSHIELD_CORR('AngPso  ',LRANGP,'iso')
               call PRINT_LRSHIELD_CORR('DiaKin  ',LRDIAK,'iso')
               call PRINT_LRSHIELD_CORR('SD      ',LRSDAV,'iso')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A)') 
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('FC      ',LRFCAV,'ani')
               call PRINT_LRSHIELD_CORR('AngPso  ',LRANGP,'ani')
               call PRINT_LRSHIELD_CORR('DiaKin  ',LRDIAK,'ani')
               call PRINT_LRSHIELD_CORR('SD      ',LRSDAV,'ani')
               END IF
            ENDIF
         ENDIF

         IF(SIGMAD1S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.DIAM1S]'
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A,10X,A)') 
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('DiaMv   ',LRDIAM,'iso')
               call PRINT_LRSHIELD_CORR('DiaDw   ',LRDIAD,'iso')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,12X,A,12X,A,13X,A)') 
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('DiaMv   ',LRDIAM,'ani')
               call PRINT_LRSHIELD_CORR('DiaDw   ',LRDIAD,'ani')
               END IF
            ENDIF
         ENDIF

         IF (PRTALL1) THEN
         WRITE(LUPRI,'(/,7X,A,13X,A,12X,A,13X,A,10X,A)') 
     &    slash,'xx','yy','zz','iso'
         IF (SIGMAD0S) call PRINT_LRSHIELD_SUM('1st Order Singlets:  ',
     &                                         SGD0S,'iso')
         IF (SIGMAD1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                                         SGD1S,'iso')
         IF (LRANISO) THEN
         WRITE(LUPRI,'(7X,A,13X,A,12X,A,13X,A)') 
     &    slash,'xy','xz','yz'
         IF (SIGMAD0S) call PRINT_LRSHIELD_SUM('1st Order Singlets:  ',
     &                                         SGD0S,'ani')
         IF (SIGMAD1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                                         SGD1S,'ani')
         END IF
         END IF
      ENDIF
C
C     All of them in detail
C
c      IF (LRESCALL) THEN
         CALL HEADER('Total Relativistic Corrections for '
     &                //NAMN(LRATOM)//' :',-1)
         WRITE(LUPRI,'(35X,A,12X,A,13X,A,10X,A)') 'xx','yy','zz','iso'
         IF(LRESCALL) call PRINT_LRSHIELD_TOT('Paramagnetic :',
     &                     SGP1S+SGP1T+SGP3S+SGP3T,'iso')
         IF(LRESCALL) call PRINT_LRSHIELD_TOT('Diamagnetic  :',
     &                     SGD0S+SGD1S,'iso')
      call PRINT_LRSHIELD_TOT('Total LRESC  :',SLRESC,'iso')
c      END IF

      IF (LRANISO) THEN
         WRITE(LUPRI,'(/,35X,A,12X,A,13X,A)') 'xy','xz','yz'
         IF(LRESCALL) call PRINT_LRSHIELD_TOT('Paramagnetic :',
     &                     SGP1S+SGP1T+SGP3S+SGP3T,'ani')
         IF(LRESCALL) call PRINT_LRSHIELD_TOT('Diamagnetic  :',
     &                     SGD0S+SGD1S,'ani')
      call PRINT_LRSHIELD_TOT('Total LRESC  :',SLRESC,'ani')

      END IF

c      END IF

      WRITE(LUPRI,*) 

      CALL QEXIT('PRINT_LRSHIELD')
      RETURN
      END

C     /* Deck print_lrshield */
      SUBROUTINE PRINT_LRSHIELD_CORR(names,tensor,elements)
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
      CHARACTER*3, intent(in)  :: elements

      CALL QENTER('PRINT_LRSHIELD_CORR')

C     Information about calculation of LRESC
      IF (elements == 'iso') THEN
         iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
         WRITE (LUPRI,'(9X,A,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &   names, ':',tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'
      END IF
      IF (elements == 'ani') THEN
         WRITE (LUPRI,'(9X,A,A,9X,F13.5,X,F13.5,X,F13.5,13X,A)')
     &   names, ':',tensor(1,2),tensor(1,3),tensor(2,3),' [ppm]'
      END IF

      CALL QEXIT('PRINT_LRSHIELD_CORR')
      RETURN
      END


C     /* Deck print_lrshield_SUM */
      SUBROUTINE PRINT_LRSHIELD_SUM(phrase,tensor,elements)
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
      CHARACTER*3, intent(in)  :: elements

      CALL QENTER('PRINT_LRSHIELD_SUM')

C     Information about calculation of LRESC
      IF (elements == 'iso') THEN
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(4X,A,2X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &      phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'
      END IF
      IF (elements == 'ani') THEN
      WRITE (LUPRI,'(4X,A,2X,F13.5,X,F13.5,X,F13.5,13X,A)')
     &      phrase,tensor(1,2),tensor(1,3),tensor(2,3),' [ppm]'
      END IF

      CALL QEXIT('PRINT_LRSHIELD_SUM')
      RETURN
      END

C     /* Deck print_lrshield_TOT */
      SUBROUTINE PRINT_LRSHIELD_TOT(phrase,tensor,elements)
C Purpose:
C     Print information about LRESC-SHIELDING overall contributions
C Author:
C     J. J. Aucar 2022
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"

      CHARACTER*3, intent(in)  :: elements
      REAL*8, intent(in)        :: tensor(3,3)
      REAL*8                    :: iso
      CHARACTER*14, intent(in)  :: phrase

      CALL QENTER('PRINT_LRSHIELD_TOT')

C     Information about calculation of LRESC
      IF (elements == 'iso') THEN
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(5X,A,8X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &      phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'
      END IF
      IF (elements == 'ani') THEN
      WRITE (LUPRI,'(5X,A,8X,F13.5,X,F13.5,X,F13.5,14X,A)')
     &      phrase,tensor(1,2),tensor(1,3),tensor(2,3),'[ppm]'
      END IF

      CALL QEXIT('PRINT_LRSHIELD_TOT')
      RETURN
      END
