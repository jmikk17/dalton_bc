
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
C     /* Deck print_lrefg */
      SUBROUTINE PRINT_LREFG
C Purpose:
C     Print LRESC-EFG results section
C Author:
C     J. J. Aucar 2024
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "mxcent.h"
#include "nuclei.h"

      REAL*8 :: REL
     
      CALL QENTER('PRINT_LREFG')    


C     NR value
      IF (PRTALL2) THEN
         WRITE(LUPRI,*)
         CALL HEADER('Non-relativistic Values (NR) [au]',-1)
         WRITE(LUPRI,'(10X,A,19X,A)') 'Atom','NR'
         WRITE(LUPRI,'(10X,A)') REPEAT('-',28)
         IJ = 0
         DO I=1,NUCIND
            WRITE (LUPRI,'(10X,A,9X,F15.5)')
     &                    NAMN(I),EFGC0(I)
            IJ = IJ + 1
         END DO
      ENDIF



C     First order corrections
      IF (PRTALL2.AND.EFGC2FLAG) THEN
         WRITE(LUPRI,*)
         CALL HEADER('First Order Corrections (1st Order) [au]',-1)
      WRITE(LUPRI,'(10X,A)') 'Atom           Mv              Dw      '//
     &               '       pqp             kin             Total'


         WRITE(LUPRI,'(10X,A)') REPEAT('-',83)
         IJ = 0
         DO I=1,NUCIND
            WRITE(LUPRI,'(10X,A,5(F15.5,X))')
     &      NAMN(I),EFGC2(I,1),EFGC2(I,2),
     &      EFGC2(I,3),EFGC2(I,4),EFGC2(I,6)
            IJ = IJ + 1
         END DO

         IF (PQPKINLRESC.AND.EFGC2FLAG) THEN
            WRITE(LUPRI,'(2/,10X,A)') 'Extra info :' 
            WRITE(LUPRI,'(10X,A)') REPEAT('-',12)
            WRITE(LUPRI,'(10X,A)')'Lap should approach'//
     &               '  pqpkin in the basis set limit'
            WRITE(LUPRI,*)' '
            WRITE(LUPRI,*)'         Atom           pqp'//
     &              '     +       kin    =       pqpkin'//
     %              '          Lap'
            WRITE(LUPRI,'(10X,A)') REPEAT('-',67)
            IJ = 0
            DO I=1,NUCIND
               WRITE (LUPRI,'(10X,A,4(F15.5,X))')
     &      NAMN(I),EFGC2(I,3),EFGC2(I,4),
     &      EFGC2(I,3)+EFGC2(I,4),+EFGC2(I,5)
               IJ = IJ + 1
            END DO
         END IF

      ENDIF


C     Second order corrections
      IF (EFGC4FLAG) THEN
      WRITE(LUPRI,*) 
      CALL HEADER('Second Order Corrections (2nd Order) [au]',-1)

      WRITE(LUPRI,*) 
      WRITE(LUPRI,*) '         Atom           '//
     &         'SO-L            SO-S'
      WRITE(LUPRI,'(10X,A)') REPEAT('-',36)
      DO I=1,NUCIND
         WRITE (LUPRI,'(10X,A,A,F15.5,X,F15.5)')
     &   '@',NAMN(I),EFGC4(I,1),EFGC4(I,4)+EFGC4(I,5)+EFGC4(I,6)
      END DO
      WRITE(LUPRI,*) 

c      WRITE(LUPRI,*) '         Atom           '//
c     &         'MV-pqp          Dw-pqp'
c      WRITE(LUPRI,'(10X,A)') REPEAT('-',52)
c      DO I=1,NUCIND
c         WRITE (LUPRI,'(10X,A,A,F15.5,X,F15.5)')
c     &       '@',NAMN(I),EFGC4(I,2),EFGC4(I,3)
c      END DO
c      WRITE(LUPRI,*) 

      END IF
      
C     Sum of contributions
      IF (EFGC2FLAG) THEN
         WRITE(LUPRI,'(3/,27X,A)') REPEAT('=',29)
         WRITE(LUPRI,'(27X,A)') 'Total EFG (zz component) [au]'
         WRITE(LUPRI,'(27X,A)') REPEAT('=',29)
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) '          Atom           NR           '//
     &                  '1st Order      2nd Order          Total'

         WRITE(LUPRI,'(10X,A)') REPEAT('-',68)
         IJ = 0
         REL = 0
         DO I=1,NUCIND
            REL = EFGC2(I,6) + EFGC4(I,7)
            WRITE (LUPRI,'(10X,A,A,F15.5,X,F15.5,X,F15.5,X,F15.5)')
     &                   '@',NAMN(I),EFGC0(I),EFGC2(I,6),
     &                    EFGC4(I,7),REL+EFGC0(I)
            IJ = IJ + 1
         END DO
         WRITE(LUPRI,*) 
      END IF


      CALL QEXIT('PRINT_LREFG')
      RETURN
      END

