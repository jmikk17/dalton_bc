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
C  /* Deck lrscinp */
C This module read nesserary input for LRESC calculations.
C jim-gesc : May-2012
C
C
      SUBROUTINE LRSCINP_EFG(WORD)
C
C Purpose:
C     Options initialization of LRESC modules to EFG.
C     This routine is called by ABAINP : abacus/abadrv.F
C Author:
c	Juan J. Aucar, 2021
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "lrescinf.h"
#include "nuclei.h"
#include "abainf.h"
      PARAMETER (NTABLE = 6, D0 = 0.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "spnout.h"
C                      1         2         3         4
      DATA TABLE /'.PRINT ','.PRTALL','.ORBCON','.PQPKIN',
C                      5         6
     &             '.C2    ','.C4    '/
C
c      write(lupri,*) ' LRSCOPTS :  called with : ' , WORD

      NEWDEF = (WORD .EQ. '*LROPTE')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
c           write(lupri,*) ' @LRSCINP : reading : ' , WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                  GOTO (1,2,3,4,5,6), I
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                  CALL PRTAB(NTABLE,TABLE,WORD1//
     &             ' input keywords',LUPRI)
                  GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &               '" not recognized in LRESC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//
     &         ' input keywords',LUPRI)
               CALL QUIT('Input keywords in LRESC Input, LRSCINP')



 1             CONTINUE
C   SET PRINT LEVEL
                  READ (LUCMD,*) JJAPRT
c                 write(lupri,*)' ..Reading print level :', JJAPRT
                  ICHANG = ICHANG + 1
                  GOTO 100
 2             CONTINUE
                  PRTALL2 = .TRUE.
                  GOTO 100
 3            CONTINUE
                 ORBCON=.TRUE. !Orbital contributions por LRESC corrections to EFG expectation values
               GOTO 100
 4             CONTINUE
                  PRTALL2 = .TRUE.
                  PQPKINLRESC= .TRUE.
            GOTO 100
 5             CONTINUE
               EFGC2flag = .TRUE.
               EFGC4flag = .FALSE.
               GOTO 100
 6             CONTINUE
               EFGC2flag = .FALSE.
               EFGC4flag = .TRUE.
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &                             '" not recognized in LRESC'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in LRESC')
            END IF
         END IF
 300  CONTINUE

      RETURN
      END
C -----------------------------------------------------------------------
C
C  /* Deck lrscini */
      SUBROUTINE LRSCINI_EFG
C
C Purpose:
C     Initialize the variables common LRESCINI: include/lrescinf.h.
C     This routine is call by ABAINIALL : abacus/abadrv.F
C Author:
C     Juan J. Aucar       April 2021 (EFG)
C
#include "implicit.h"
#include "lrescinf.h"
#include "mxcent.h"
#include "cbiqr.h"
#include "abainf.h"
#include "inforb.h"              
C
      PARAMETER(D0 = 0.0D+0)

      JJAPRT    = 0
      EFGC2 = 0.0D0
      EFGC4 = 0.0D0

      ORBCON=.FALSE.
      PQPKINLRESC= .FALSE.
      EFGC2FLAG = .TRUE.
      EFGC4FLAG = .TRUE.
C
C     Initialize /CUADRA/
C
      IPRINT = IPRDEF
      IPRQR  = IPRINT
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      THRESH = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
      NABAPP = 0
      LBFREQ = 1
      LCFREQ = 1
      CALL DZERO (QBFREQ,NFMAX)
      CALL DZERO (QCFREQ,NFMAX)
C
      RETURN
      END
c    ---------------------------------------
C  /* Deck lrscdrv */
      SUBROUTINE LRSCDRV_EFG(WORK,LWORK)
C
C Purpose:
C     Main routine of LRESC module.
C     This routine is call by ABACTL : abacus/abadrv.F
C Author:
C     Juan J. Aucar, 2021. Revisited in 2024
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
#include "nuclei.h"
C
C
      CALL QENTER('LRSCDRV_EFG')

C     Print information about LRESC metodology
      CALL HEADER_LREFG

c c^-2
      CALL AVELRSC_EFG('EFG ',WORK,LWORK) !Contains NR value
      IF (EFGC2FLAG) THEN
            IF (PQPKINLRESC) CALL AVELRSC_EFG('EFG2',WORK,LWORK)
            CALL LINEARLR_EFG('J-Dw',WORK,LWORK) 
            CALL LINEARLR_EFG('J-Mv',WORK,LWORK)
      END IF

c c^-4
      IF (EFGC4FLAG) THEN
c            call generate_xyz() !Generate xyz file. For use with pySCF
c            CALL LINEARLR_EFG('Mvc4',WORK,LWORK) ! From SS contribution
c            CALL LINEARLR_EFG('Dwc4',WORK,LWORK) ! From SS contribution
c           Triplet responses must be calculated at the end            
c            CALL LINEARLR_EFG('JSO2',WORK,LWORK) !wrong results...

C           Import generated PySCF matrix to AOPROPER for SO-S
            call import_matrix(3*NUCIND+3,"efg")
            CALL LINEARLR_EFG('SOSx',WORK,LWORK)
            CALL LINEARLR_EFG('SOSy',WORK,LWORK)
            CALL LINEARLR_EFG('SOSz',WORK,LWORK)

c           SO-L correction
            CALL QRLREFG('EFG-SO-SO',WORK,LWORK)
      END IF


C     First and Second order Total Contributions
      DO I=1,NUCIND
         EFGC2(I,6)=EFGC2(I,1)+EFGC2(I,2)+EFGC2(I,3)+EFGC2(I,4)
         EFGC4(I,7)=EFGC4(I,1)+EFGC4(I,4)+EFGC4(I,5)+EFGC4(I,6)
      END DO

c     Print Results
      call PRINT_LREFG

      CALL QEXIT('LRSCDRV_EFG')
      RETURN
      END
c    ---------------------------------------
c    ---------------------------------------