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
      PARAMETER (NTABLE = 4, D0 = 0.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "spnout.h"
C                      1         2         3		4
      DATA TABLE /'.PRINT ','.PRTALL','.ORBCON','.PQPKIN'/
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
                  GOTO (1,2,3,4), I
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
c jim-dbg ICHANG > 0 default variables change
c jim-dbg ICHANG < 0 no changes in default variables


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

      ORBCON=.FALSE.
      PQPKINLRESC= .FALSE.
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
C     Juan J. Aucar, 2021
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
C this is for print section
#include "nuclei.h"
C
C
      LOGICAL PING
      PING=.FALSE.
cv      write(lupri,*)'@ LRSCDRV : '
      CALL QENTER('LRSCDRV_EFG')



c c^-2
            CALL AVELRSC_EFG('EFG ',WORK,LWORK)
            CALL AVELRSC_EFG('EFG2',WORK,LWORK)
            CALL LINEARLR_EFG('J-Dw',WORK,LWORK) 
            CALL LINEARLR_EFG('J-Mv',WORK,LWORK)


      CALL TITLER('EFG(ZZ) LRESC RELATIVISTIC CORRECTIONS','*',124)
      WRITE(lupri,*)
      WRITE(lupri,*)'   First order relativistic '// 
     & 'corrections to the zz component of the electronic EFG.'
      WRITE(LUPRI,*)
       WRITE(LUPRI,*)'   -  Original LRESC papers at :'
       WRITE(LUPRI,'(33X,A)')' J. Chem. Phys.'//
     &               ' 118, 2 (2003), doi: 10.1063/1.1525808' 
       WRITE(LUPRI,'(33X,A)')' Int. J. Quantum Chem. (2019), '//
     & 'doi: 10.1002/qua.25935'
       WRITE(LUPRI,'(33X,A)')' Int. J. Quantum Chem. (2021), '//
     & 'doi: 10.1002/qua.26769'
   

      WRITE(LUPRI,*)
      WRITE(LUPRI,*)'   -  Corrections calculated in all symmetry '// 
     & 'independent centers'

      IF (PRTALL2) THEN 
      WRITE(LUPRI,*)'   -  To identify LRESC corrections '// 
     & 'see related papers'
      END IF
        
      IF (.NOT. PRTALL2) THEN 
         WRITE(LUPRI,'(14X,A)')'*  For LRESC corrections, add '//
     &   '.PRTALL in your input file under *LROPTE'
      ENDIF
    
C
C NR value
C
      IF (PRTALL2) THEN
            WRITE(LUPRI,*)
            CALL HEADER('Non-relativistic Values (NR) [au]',-1)


            WRITE(LUPRI,'(10X,A,19X,A)') 'Atom','NR'
            WRITE(LUPRI,'(10X,A)') REPEAT('-',28)
            IJ = 0
            DO I=1,NUCIND
                  WRITE (LUPRI,'(10X,A,9X,F15.5)')
     &            NAMN(I),EFGC0(I)
                  IJ = IJ + 1
            END DO
       ENDIF




C First order corrections
C
      IF (PRTALL2) THEN
            WRITE(LUPRI,*)
            CALL HEADER('First Order Corrections (1st Order) [au]',-1)
            WRITE(LUPRI,*)'         Atom           Mv  '//
     &              '            Dw              Lap            Total'
            WRITE(LUPRI,'(10X,A)') REPEAT('-',67)
            IJ = 0
            DO I=1,NUCIND
               relall=EFGC2(I,1)+EFGC2(I,2)+EFGC2(I,5)
               WRITE (LUPRI,'(10X,A,F15.5,X,F15.5,X,F15.5,X,F15.5)')
     &         NAMN(I),EFGC2(I,1),EFGC2(I,2),EFGC2(I,5),relall
               IJ = IJ + 1
            END DO


            IF (PQPKINLRESC) THEN
                  WRITE(LUPRI,*)' '
                  WRITE(LUPRI,*)' '
                  WRITE(LUPRI,'(10X,A)') 'Extra info :' 
                  WRITE(LUPRI,'(10X,A)') REPEAT('-',12)
                  WRITE(LUPRI,'(10X,A)')' pqpkin should approach'//
     &            '  Lap in the basis set limit'
                  WRITE(LUPRI,*)' '
                  WRITE(LUPRI,*)'         Atom           pqp'//
     &              '     +       kin    =        pqpkin'
                  WRITE(LUPRI,'(10X,A)') REPEAT('-',53)
                  IJ = 0
                  DO I=1,NUCIND
                        WRITE (LUPRI,'(10X,A,F15.5,X,F15.5,X,F15.5)')
     &               NAMN(I),EFGC2(I,3),EFGC2(I,4),EFGC2(I,3)+EFGC2(I,4)
                        IJ = IJ + 1
                  END DO
            END IF

      ENDIF
C


      WRITE(LUPRI,*) 
      WRITE(LUPRI,*) 
      WRITE(LUPRI,*) 
      WRITE(LUPRI,'(27X,A)') REPEAT('=',29)
      WRITE(LUPRI,'(27X,A)') 'Total EFG (zz component) [au]'
      WRITE(LUPRI,'(27X,A)') REPEAT('=',29)
      WRITE(LUPRI,*) 
       
      WRITE(LUPRI,*) '         Atom           '//
     &         'NR           1st Order          Total'

      WRITE(LUPRI,'(10X,A)') REPEAT('-',52)
      IJ = 0
      DO I=1,NUCIND
         relall=EFGC2(I,1)+EFGC2(I,2)+EFGC2(I,5)
         total=EFGC0(I)+relall
         WRITE (LUPRI,'(10X,A,A,F15.5,X,F15.5,X,F15.5)')
     &       '@',NAMN(I),EFGC0(I),relall,total
         IJ = IJ + 1
      END DO
      WRITE(LUPRI,*) 


      CALL QEXIT('LRSCDRV_EFG')
      RETURN
      END
c    ---------------------------------------
c    ---------------------------------------

