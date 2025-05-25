
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
      USE lrescinf_module
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "mxcent.h"
#include "nuclei.h"

      CHARACTER*15  :: slash
      CHARACTER*9  :: aba,lib,avg

      slash = '---------------'


      aba = '(abacus) '
      lib = '(libcint)'
      avg = '(avg)    '

      CALL QENTER('PRINT_LRSHIELD')
C     Paramagnetic 1st singlet or triplet
C
      IF ((SIGMAP1S).OR.(SIGMAP1T).OR.(SIGMAP3S).OR.
     &   (SIGMAP3T)) THEN
         IF (PRTALL1) THEN
         CALL HEADER('Paramagnetic Corrections for '
     &                //NAMN(LRATOM)//' :',-1)
         END IF

         IF(SIGMAP1S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA1S]'
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A,10X,A)')
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('OZK     ',aba,LROZK,'iso ')
               call PRINT_LRSHIELD_CORR('PSOK    ',avg,LRPSOK,'iso ')
               IF (EXTRAINFO) call PRINT_LRSHIELD_CORR
     &            ('PSOK    ',aba,LRPSOKaba,'iso ')
               IF (LRANISO) THEN
                  WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &                  slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('OZK     ',aba,LROZK,'ani1')
               call PRINT_LRSHIELD_CORR('PSOK    ',avg,LRPSOK,'ani1')
               IF (EXTRAINFO) call PRINT_LRSHIELD_CORR
     &            ('PSOK    ',aba,LRPSOKaba,'ani1')
                  WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &                  slash,'yx','zx','zy'
                  call PRINT_LRSHIELD_CORR('OZK     ',aba,LROZK,'ani2')
                  call PRINT_LRSHIELD_CORR('PSOK    ',avg,LRPSOK,'ani2')
                  IF (EXTRAINFO) call PRINT_LRSHIELD_CORR
     &               ('PSOK    ',aba,LRPSOKaba,'ani2')
               END IF
            ENDIF
         ENDIF

         IF (SIGMAP1T) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA1T]'
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A,10X,A)')
     &               slash,'xx','yy','zz','iso'
               IF (EXTRAINFO) THEN
                  call PRINT_LRSHIELD_CORR('Fc-Kin  ',aba,LRFCZK,'iso ')
                  call PRINT_LRSHIELD_CORR('Fc-Der  ',aba,LRFCDE,'iso ')
                  call PRINT_LRSHIELD_CORR('Sd-Kin  ',aba,LRSDK,'iso ')
                  call PRINT_LRSHIELD_CORR('Sd-Der  ',aba,LRSDD,'iso ')
               ELSE
                  call PRINT_LRSHIELD_CORR
     &            ('SZK(Fc)   ',aba,LRFCZK+LRFCDE,'iso ')
                  call PRINT_LRSHIELD_CORR
     &            ('SZK(Sd) ',aba,LRSDK+LRSDD,'iso ')
               END IF
               call PRINT_LRSHIELD_CORR('BSO(Fc) ',aba,LRFCBS,'iso ')
               call PRINT_LRSHIELD_CORR('BSO(Sd) ',aba,LRSDBS,'iso ')
               IF (LRANISO) THEN
                  WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &                  slash,'xy','xz','yz'
               IF (EXTRAINFO) THEN
                  call PRINT_LRSHIELD_CORR('Fc-Kin  ',aba,LRFCZK,'ani1')
                  call PRINT_LRSHIELD_CORR('Fc-Der  ',aba,LRFCDE,'ani1')
                  call PRINT_LRSHIELD_CORR('Sd-Kin  ',aba,LRSDK,'ani1')
                  call PRINT_LRSHIELD_CORR('Sd-Der  ',aba,LRSDD,'ani1')
               ELSE
                  call PRINT_LRSHIELD_CORR
     &            ('SZK(Fc) ',aba,LRFCZK+LRFCDE,'ani1')
                  call PRINT_LRSHIELD_CORR
     &            ('SZK(Sd) ',aba,LRSDK+LRSDD,'ani1')
               END IF
               call PRINT_LRSHIELD_CORR('BSO(Fc) ',aba,LRFCBS,'ani1')
               call PRINT_LRSHIELD_CORR('BSO(Sd) ',aba,LRSDBS,'ani1')
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'yx','zx','zy'
               IF (EXTRAINFO) THEN
                  call PRINT_LRSHIELD_CORR('Fc-Kin  ',aba,LRFCZK,'ani2')
                  call PRINT_LRSHIELD_CORR('Fc-Der  ',aba,LRFCDE,'ani2')
                  call PRINT_LRSHIELD_CORR('Sd-Kin  ',aba,LRSDK,'ani2')
                  call PRINT_LRSHIELD_CORR('Sd-Der  ',aba,LRSDD,'ani2')
               ELSE
                  call PRINT_LRSHIELD_CORR
     &            ('SZK(Fc) ',aba,LRFCZK+LRFCDE,'ani1')
                  call PRINT_LRSHIELD_CORR
     &            ('SZK(Sd) ',aba,LRSDK+LRSDD,'ani2')
               END IF
               call PRINT_LRSHIELD_CORR('BSO(Fc) ',aba,LRFCBS,'ani2')
               call PRINT_LRSHIELD_CORR('BSO(Sd) ',aba,LRSDBS,'ani2')
               END IF
            ENDIF
         ENDIF

         IF (SIGMAP3S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA3S]'
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A,10X,A)')
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('Mv/p    ',aba,QRPSOOZMV,'iso ')
               call PRINT_LRSHIELD_CORR('Dw/p    ',aba,QRPSOOZDW,'iso ')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('Mv/p    ',aba,QRPSOOZMV,'ani1')
               call PRINT_LRSHIELD_CORR('Dw/p    ',aba,QRPSOOZDW,'ani1')
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'yx','zx','zy'
               call PRINT_LRSHIELD_CORR('Mv/p    ',aba,QRPSOOZMV,'ani2')
               call PRINT_LRSHIELD_CORR('Dw/p    ',aba,QRPSOOZDW,'ani2')
               END IF
            ENDIF
         ENDIF

         IF (SIGMAP3T) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.PARA3T]'
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A,10X,A)')
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('SO(Fc)  ',aba,QROZFCSO,'iso ')
               call PRINT_LRSHIELD_CORR('SO(Sd)  ',aba,QROZSDSO,'iso ')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('SO(Fc)  ',aba,QROZFCSO,'ani1')
               call PRINT_LRSHIELD_CORR('SO(Sd)  ',aba,QROZSDSO,'ani1')
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'yx','zx','zy'
               call PRINT_LRSHIELD_CORR('SO(Fc)  ',aba,QROZFCSO,'ani2')
               call PRINT_LRSHIELD_CORR('SO(Sd)  ',aba,QROZSDSO,'ani2')
               END IF
            ENDIF
         ENDIF

         IF (PRTALL1) THEN
         WRITE(LUPRI,'(/,7X,A,17X,A,12X,A,13X,A,10X,A)') 
     &    slash,'xx','yy','zz','iso'
         IF (SIGMAP1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                      SGP1S,'iso ')
         IF (SIGMAP1T) call PRINT_LRSHIELD_SUM('2nd Order Triplets:  ',
     &                      SGP1T,'iso ')
         IF (SIGMAP3S) call PRINT_LRSHIELD_SUM('3nd Order Singlets:  ',
     &                      SGP3S,'iso ')
         IF (SIGMAP3T) call PRINT_LRSHIELD_SUM('3nd Order Triplets:  ',
     &                      SGP3T,'iso ')
         IF (LRANISO) THEN
         WRITE(LUPRI,'(7X,A,17X,A,12X,A,13X,A)') 
     &    slash,'xy','xz','yz'
         IF (SIGMAP1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                      SGP1S,'ani1')
         IF (SIGMAP1T) call PRINT_LRSHIELD_SUM('2nd Order Triplets:  ',
     &                      SGP1T,'ani1')
         IF (SIGMAP3S) call PRINT_LRSHIELD_SUM('3nd Order Singlets:  ',
     &                      SGP3S,'ani1')
         IF (SIGMAP3T) call PRINT_LRSHIELD_SUM('3nd Order Triplets:  ',
     &                      SGP3T,'ani1')
         WRITE(LUPRI,'(7X,A,17X,A,12X,A,13X,A)') 
     &    slash,'yx','zx','zy'
         IF (SIGMAP1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                      SGP1S,'ani2')
         IF (SIGMAP1T) call PRINT_LRSHIELD_SUM('2nd Order Triplets:  ',
     &                      SGP1T,'ani2')
         IF (SIGMAP3S) call PRINT_LRSHIELD_SUM('3nd Order Singlets:  ',
     &                      SGP3S,'ani2')
         IF (SIGMAP3T) call PRINT_LRSHIELD_SUM('3nd Order Triplets:  ',
     &                      SGP3T,'ani2')
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
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A,10X,A)')
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('FC/d    ',aba,LRFCAV,'iso ')
               if (PYINTS) then
                  call PRINT_LRSHIELD_CORR
     &            ('FC/d    ',lib,LRFCAVpy,'iso ')
               end if
               call PRINT_LRSHIELD_CORR('SD/d    ',aba,LRSDAV,'iso ')
               call PRINT_LRSHIELD_CORR
     &         ('PSO-OZ  ',avg,LRANGPavg,'iso ')
               if (EXTRAINFO) then
                  call PRINT_LRSHIELD_CORR('PSO-OZ  ',aba,LRANGP,'iso ')
               end if
               call PRINT_LRSHIELD_CORR('DIAK    ',aba,LRDIAK,'iso ')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('FC/d    ',aba,LRFCAV,'ani1')
               if (PYINTS) then
               call PRINT_LRSHIELD_CORR('FC/d    ',lib,LRFCAVpy,'ani1')
               end if
               call PRINT_LRSHIELD_CORR('SD/d    ',aba,LRSDAV,'ani1')
               call PRINT_LRSHIELD_CORR
     &         ('PSO-OZ  ',avg,LRANGPavg,'ani1')
               if (EXTRAINFO) then
                  call PRINT_LRSHIELD_CORR('PSO-OZ  ',aba,LRANGP,'ani1')
               end if
               call PRINT_LRSHIELD_CORR('DIAK    ',aba,LRDIAK,'ani1')
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'yx','zx','zy'
               call PRINT_LRSHIELD_CORR('FC/d    ',aba,LRFCAV,'ani2')
               call PRINT_LRSHIELD_CORR('SD/d    ',aba,LRSDAV,'ani2')
               if (PYINTS) then
               call PRINT_LRSHIELD_CORR('FC/d    ',lib,LRFCAVpy,'ani2')
               end if
               call PRINT_LRSHIELD_CORR
     &         ('PSO-OZ  ',avg,LRANGPavg,'ani2')
               if (EXTRAINFO) then
                  call PRINT_LRSHIELD_CORR('PSO-OZ  ',aba,LRANGP,'ani2')
               end if
               call PRINT_LRSHIELD_CORR('DIAK    ',aba,LRDIAK,'ani2')
               END IF
            ENDIF
         ENDIF

         IF(SIGMAD1S) THEN
            IF(PRTALL1) THEN
               WRITE(LUPRI,'(/,7X,A)') 'Detailed info : [.DIAM1S]'
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A,10X,A)')
     &               slash,'xx','yy','zz','iso'
               call PRINT_LRSHIELD_CORR('Mv/d    ',aba,LRDIAM,'iso ')
               call PRINT_LRSHIELD_CORR('Dw/d    ',aba,LRDIAD,'iso ')
               IF (LRANISO) THEN
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'xy','xz','yz'
               call PRINT_LRSHIELD_CORR('Mv/d    ',aba,LRDIAM,'ani1')
               call PRINT_LRSHIELD_CORR('Dw/d    ',aba,LRDIAD,'ani1')
               WRITE(LUPRI,'(8X,A,16X,A,12X,A,13X,A)')
     &               slash,'yx','zx','zy'
               call PRINT_LRSHIELD_CORR('Mv/d    ',aba,LRDIAM,'ani2')
               call PRINT_LRSHIELD_CORR('Dw/d    ',aba,LRDIAD,'ani2')
               END IF
            ENDIF
         ENDIF

         IF (PRTALL1) THEN
         WRITE(LUPRI,'(/,7X,A,17X,A,12X,A,13X,A,10X,A)') 
     &    slash,'xx','yy','zz','iso'
         IF (SIGMAD0S) call PRINT_LRSHIELD_SUM('1st Order Singlets:  ',
     &                                         SGD0S,'iso ')
         IF (SIGMAD1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                                         SGD1S,'iso ')
         IF (LRANISO) THEN
         WRITE(LUPRI,'(7X,A,17X,A,12X,A,13X,A)') 
     &    slash,'xy','xz','yz'
         IF (SIGMAD0S) call PRINT_LRSHIELD_SUM('1st Order Singlets:  ',
     &                                         SGD0S,'ani1')
         IF (SIGMAD1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                                         SGD1S,'ani1')
         WRITE(LUPRI,'(7X,A,17X,A,12X,A,13X,A)') 
     &    slash,'yx','zx','zy'
         IF (SIGMAD0S) call PRINT_LRSHIELD_SUM('1st Order Singlets:  ',
     &                                         SGD0S,'ani2')
         IF (SIGMAD1S) call PRINT_LRSHIELD_SUM('2nd Order Singlets:  ',
     &                                         SGD1S,'ani2')
         END IF
         END IF
      ENDIF

C     Non-rel contributions
C
         CALL HEADER('Total Non-relativistic Contributions for '
     &                //NAMN(LRATOM)//' :',-1)

      IF (.NOT.LRANISO) THEN
         WRITE(LUPRI,'(35X,A,12X,A,13X,A,10X,A,9X,A,9X,A)') 
     &               'xx','yy','zz','iso','anis','asym'

         IF(LRESCALL) call PRINT_LRSHIELD_TOT('Paramagnetic :',
     &                     SGPNR,'iso ')
         IF(LRESCALL) call PRINT_LRSHIELD_TOT('Diamagnetic  :',
     &                     SGDNR,'iso ')
      call PRINT_LRSHIELD_TOT('Total non-rel  :',SGPNR+SGDNR,'ist ')
      END IF

c     Print full (non-rel) matrices
      IF (LRANISO) THEN
        CALL HEADER('Diamagnetic contributions (ppm):',1)
        CALL PRINT_TENSOR(SGDNR,'tensor')

        CALL HEADER('Paramagnetic contributions (ppm):',1)
        CALL PRINT_TENSOR(SGPNR,'tensor')

        CALL HEADER('Total shielding tensor (ppm):',1)
        CALL PRINT_TENSOR(SGPNR+SGDNR,'tensor')

        CALL HEADER('Principal values and axes:',1)
        CALL PRINT_TENSOR(SGPNR+SGDNR,'eigens')

        CALL HEADER('Anisotropic parameters:',1)
        CALL PRINT_TENSOR(SGPNR+SGDNR,'aniso ')
      END IF
C
C     Relativistic corrections
      CALL HEADER('Total Relativistic Corrections for '
     &                //NAMN(LRATOM)//' :',-1)

      IF (.NOT.LRANISO) THEN
        WRITE(LUPRI,'(35X,A,12X,A,13X,A,10X,A)') 'xx','yy','zz','iso'
        IF(LRESCALL) call PRINT_LRSHIELD_TOT('Paramagnetic :',
     &                       SGP1S+SGP1T+SGP3S+SGP3T,'iso ')
        IF(LRESCALL) call PRINT_LRSHIELD_TOT('Diamagnetic  :',
     &                       SGD0S+SGD1S,'iso ')
        call PRINT_LRSHIELD_TOT('Total LRESC  :',SLRESC,'iso ')
      END IF

c     Print full (rel corrections) matrices
      IF (LRANISO) THEN
        CALL HEADER('Diamagnetic corrections (ppm):',1)
        CALL PRINT_TENSOR(SGD0S+SGD1S,'tensor')

        CALL HEADER('Paramagnetic corrections (ppm):',1)
        CALL PRINT_TENSOR(SGP1S+SGP1T+SGP3S+SGP3T,'tensor')

        CALL HEADER('Total corrections shielding tensor (ppm):',1)
        CALL PRINT_TENSOR(SLRESC,'tensor')
      END IF

C
c     Total relativistic values
      CALL HEADER('Total Relativistic Values for '
     &                //NAMN(LRATOM)//' :',-1)

      IF (.NOT.LRANISO) THEN
        WRITE(LUPRI,'(35X,A,12X,A,13X,A,10X,A,9X,A,9X,A)')
     &                    'xx','yy','zz','iso','anis','asym'

        IF(LRESCALL) call PRINT_LRSHIELD_TOT('Paramagnetic :',
     &                       SGPNR+SGP1S+SGP1T+SGP3S+SGP3T,'iso ')
        IF(LRESCALL) call PRINT_LRSHIELD_TOT('Diamagnetic  :',
     &                       SGDNR+SGD0S+SGD1S,'iso ')
        call PRINT_LRSHIELD_TOT('Total LRESC  :',
     &                          SGPNR+SGDNR+SLRESC,'ist ')
      END IF

c     Print full (relativistic) matrices
      IF (LRANISO) THEN
        CALL HEADER('Diamagnetic contributions (ppm):',1)
        CALL PRINT_TENSOR(SGDNR+SGD0S+SGD1S,'tensor')

        CALL HEADER('Paramagnetic contributions (ppm):',1)
        CALL PRINT_TENSOR(SGPNR+SGP1S+SGP1T+SGP3S+SGP3T,'tensor')

        CALL HEADER('Total shielding tensor (ppm):',1)
        CALL PRINT_TENSOR(SGPNR+SGDNR+SLRESC,'tensor')

        CALL HEADER('Principal values and axes:',1)
        CALL PRINT_TENSOR(SGPNR+SGDNR+SLRESC,'eigens')

        CALL HEADER('Anisotropic parameters:',1)
        CALL PRINT_TENSOR(SGPNR+SGDNR+SLRESC,'aniso ')
      END IF

      WRITE(LUPRI,*) 

      CALL QEXIT('PRINT_LRSHIELD')
      RETURN
      END

C     /* Deck print_lrshield */
      SUBROUTINE PRINT_LRSHIELD_CORR(names,integrals,tensor,elements)
C Purpose:
C     Print information about LRESC-SHIELDING correction
C Author:
C     J. J. Aucar 2022
      USE lrescinf_module
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"

      REAL*8, intent(in)       :: tensor(3,3)
      REAL*8                   :: iso
      CHARACTER*8, intent(in)  :: names
      CHARACTER*9, intent(in)  :: integrals
      CHARACTER*9              :: loc_integrals
      CHARACTER*4, intent(in)  :: elements

      CALL QENTER('PRINT_LRSHIELD_CORR')
      IF (EXTRAINFO) THEN
         loc_integrals = integrals
      ELSE
         loc_integrals = '         '
      END IF


C     Information about calculation of LRESC
      IF (elements == 'iso ') THEN
         iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
         WRITE (LUPRI,'(9X,2A,4X,A,F12.4,2X,F12.4,3X,2(F12.4),A)')
     &   names, loc_integrals,':',
     &   tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'
      ELSEIF (elements == 'ani1') THEN
         WRITE (LUPRI,'(9X,2A,4X,A,F12.4,2X,F12.4,3X,F12.4,13X,A)')
     &   names, loc_integrals, ':',
     &   tensor(1,2),tensor(1,3),tensor(2,3),'[ppm]'
      ELSEIF (elements == 'ani2') THEN
         WRITE (LUPRI,'(9X,2A,4X,A,F12.4,2X,F12.4,3X,F12.4,13X,A)')
     &   names, loc_integrals, ':',
     &   tensor(2,1),tensor(3,1),tensor(3,2),'[ppm]'
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
      USE lrescinf_module
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"

      REAL*8, intent(in)       :: tensor(3,3)
      REAL*8                   :: iso
      CHARACTER*21, intent(in)  :: phrase
      CHARACTER*4, intent(in)  :: elements

      CALL QENTER('PRINT_LRSHIELD_SUM')

C     Information about calculation of LRESC
      IF (elements == 'iso ') THEN
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(4X,A,6X,F12.4,2X,F12.4,3X,F12.4,F12.4,A)')
     &      phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'
c      END IF
c      IF (elements == 'ani1') THEN
      ELSEIF (elements == 'ani1') THEN
      WRITE (LUPRI,'(4X,A,6X,F12.4,2X,F12.4,3X,F12.4,13X,A)')
     &      phrase,tensor(1,2),tensor(1,3),tensor(2,3),'[ppm]'
      ELSEIF (elements == 'ani2') THEN
      WRITE (LUPRI,'(4X,A,6X,F12.4,2X,F12.4,3X,F12.4,13X,A)')
     &      phrase,tensor(2,1),tensor(3,1),tensor(3,2),'[ppm]'
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
      USE lrescinf_module
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"

      CHARACTER*4, intent(in)  :: elements
      REAL*8, intent(in)        :: tensor(3,3)
      REAL*8                    :: iso
      CHARACTER*14, intent(in)  :: phrase
      DOUBLE PRECISION :: parameters(4)
      DOUBLE PRECISION :: principal_axes(3,3)
      DOUBLE PRECISION :: eigenvalues(3)

      CALL QENTER('PRINT_LRSHIELD_TOT')

C     Information about calculation of LRESC
c     For total (NR and REL) also ANIS and ASMY are printed
      IF (elements == 'ist ') THEN
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
C     Initialize parameters
      DATA parameters / 0.d0, 0.d0, 0.d0, 0.d0 /
      call shielddiag(tensor,3,parameters,principal_axes,eigenvalues)
      WRITE(LUPRI,'(5X,A,8X,F12.4,2X,F12.4,
     &              3X,F12.4,F12.4,F12.4,X,F12.4,X,A)')
     & phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,parameters(3),
     & parameters(4),' [ppm]'
      END IF

      IF (elements == 'iso ') THEN
      iso =(tensor(1,1)+tensor(2,2)+tensor(3,3))/3.0
      WRITE (LUPRI,'(5X,A,8X,F12.4,2X,F12.4,
     &               3X,F12.4,F12.4,26X,A)')
     &      phrase,tensor(1,1),tensor(2,2),tensor(3,3),iso,' [ppm]'
      END IF
      IF (elements == 'ani1') THEN
      WRITE (LUPRI,'(5X,A,8X,F12.4,2X,F12.4,
     &               3X,F12.4,39X,A)')
     &      phrase,tensor(1,2),tensor(1,3),tensor(2,3),'[ppm]'
      END IF
      IF (elements == 'ani2') THEN
      WRITE (LUPRI,'(5X,A,8X,F12.4,2X,F12.4,
     &               3X,F12.4,39X,A)')
     &      phrase,tensor(2,1),tensor(3,1),tensor(3,2),'[ppm]'
      END IF

      CALL QEXIT('PRINT_LRSHIELD_TOT')
      RETURN
      END

C     /* Deck print_TENSOR */
      SUBROUTINE PRINT_TENSOR(tensor,options)
C Purpose:
C     Print LRESC-SHIELDING tensors and associated quantities.
C Author:
C     J. J. Aucar 2025
      USE lrescinf_module
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"

      REAL*8, intent(in)  :: tensor(3,3)
      DOUBLE PRECISION    :: parameters(4)
      character*1         :: coordinate(3)
      character*6         :: options
      DOUBLE PRECISION :: principal_axis(3,3)
      DOUBLE PRECISION :: eigenvalues(3)

      CALL QENTER('PRINT_TENSOR')

      coordinate(1) = 'x'
      coordinate(2) = 'y'
      coordinate(3) = 'z'

C     Initialize parameters
      DATA parameters / 0.d0, 0.d0, 0.d0, 0.d0 /
      DATA eigenvalues / 0.d0, 0.d0, 0.d0 /

      call shielddiag(tensor,3,parameters,principal_axis,eigenvalues)

C     Print full matrix
      IF (options == 'tensor') THEN
         write(lupri,'(20X,A)') 'Bx            By             Bz'
         do i = 1,3
           WRITE(LUPRI,'(8X,A,3X,F12.4,2X,F12.4,
     &                   3X,F12.4)')
     &     coordinate(i),tensor(i,1),tensor(i,2),tensor(i,3)
         end do
         WRITE(LUPRI,'(/,8X,A,1X,F12.4/)') 'iso',
     &        (tensor(1,1)+tensor(2,2)+tensor(3,3)) / 3.0D0
C     Print eigenvalues and eigenstates
      ELSE IF (options == 'eigens') THEN
         write(lupri,'(20X,A)') 'Ax            Ay             Az'
         do i = 1, 3
            write(lupri,'(F12.4,F12.4,2X,F12.4,3X,F12.4)' )
     &       eigenvalues(i),
     &       principal_axis(i,1),principal_axis(i,2),principal_axis(i,3)
         end do
C     Print anisotropic parameters
      ELSE IF (options == 'aniso ') THEN
         write(lupri,'(20X,A)') 'anis          asym'
         write(lupri,'(12X,F12.4,2X,F12.4)' ) parameters(3),
     &         parameters(4)
      END IF

      CALL QEXIT('PRINT_TENSOR')
      RETURN
      END
