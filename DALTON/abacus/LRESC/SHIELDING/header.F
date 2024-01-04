
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

C     /* Deck header_lrshield */
      SUBROUTINE HEADER_LRSHIELD
C Purpose:
C     Print information about LRESC-SHIELDING module              
C Author:
C     J. J. Aucar 2022 (moved from workplace)
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
#include "nuclei.h"

      CALL QENTER('LRSCDRV_SHI')
C
C     Information about calculation of LRESC
C
C ---- Default
      CALL TITLER
     &('LRESC RELATIVISTIC CORRECTIONS TO SHIELDING','*',124)

      WRITE(LUPRI,*)'   -  Original LRESC papers at :'
      write(lupri,'(7X,A)') 'J. Chem. Phys. 118, 2 (2003),'//
     &            'doi:10.1063/1.1525808'
c      WRITE(LUPRI,'(33X,A)')' Mol. Phys. 101, 20, 3103-3109 (2003) '
c     &                       ' DOI: 10.1080/00268970310001617784'
c      WRITE(LUPRI,'(33X,A)')' J. Chem. Phys. 125, 064107 (2006), '//
c     & 'doi:10.1063/1.2244572'
      WRITE(LUPRI,'(7X,A)')'Int. J. Quantum Chem. 118, '//
     &    'e25487 (2017), doi.org/10.1002/qua.25487'

c      IF (LRANISO) THEN
      WRITE(LUPRI,'(/,4X,A)') '-  Definitions according to : '
      WRITE(LUPRI,'(7X,A)') 'J. Mason, Solid state Nuc. Magn. '//
     &                         'Res. 2 (1993), 285 '
      WRITE(LUPRI,'(/,4X,A)') '-  Tensor treatment according to : '
      WRITE(LUPRI,'(7X,A)') 'F. A. Anet et. alter, Concepts Magn. '//
     &                         'Reson. (1991), 3 '
c      END IF

      IF (LRESCALL) THEN
            WRITE(lupri,'(/,4X,A)')
     &   '-  All LRESC corrections were activated'
      ELSE
            WRITE(lupri,'(/,4X,A)') 
     &   '-  Specific LRESC corrections were activated'
      ENDIF
      WRITE(LUPRI,*)

c     For sublist of atoms (not implemented yet) this should be modified
c      IF (LRATOM_changed) THEN
c         WRITE(LUPRI,'(4X,A,I2,A,A)') '-  Selected Atom    : #',LRATOM,
c     &         ": ", NAMN(LRATOM)
c         WRITE(LUPRI,'(4X,A,3F10.4)') '-  Atom coordinates :',
c     &   CORD(1,LRATOM), CORD(2,LRATOM),CORD(3,LRATOM)
c      ELSE
c      IF (.not.LRATOM_changed) THEN
      IF (NLRATOMS.EQ.NUCIND) THEN
              WRITE(LUPRI,'(4X,A)') '-  All Atoms selected'
      END IF

      IF (.NOT.LRATOM_changed) THEN
         WRITE(LUPRI,'(/,4X,A)')'-  Remember : '
         WRITE(LUPRI,'(7X,A)')'GAUGEO must be placed at desired '//
     &   'nucleus.'
      END IF
      IF (.NOT. PRTALL1) THEN
         WRITE(LUPRI,'(7X,A)')'For printing LRESC corrections,'//
     &   ' add .PRTALL in your input file under *LROPTS'

      ENDIF

C ---Not Default
      IF ((GAUCHANG).OR.(JIMPRT.NE.0).OR.(.NOT. LRESCALL))
     & CALL HEADER('Changes of defaults for LRESC :',3)

      IF (JIMPRT.NE.0) THEN
         WRITE (LUPRI,'(A,I5)')
     &      '--- Print level : ', JIMPRT
      END IF

      IF (.NOT. LRESCALL) THEN
         WRITE (LUPRI,'(4X,A)')
     &      '-  Not all corrections were activated. ' //
     &      'Only the following will be done : '
         IF (SIGMAD0S) WRITE (LUPRI,'(7X,A)')
     &   'Diamagnetic zero order singlet corrections'
         IF (SIGMAD1S) WRITE (LUPRI,'(7X,A)')
     &   'Diamagnetic first order singlet corrections'
         IF (SIGMAP1S)  WRITE (LUPRI,'(7X,A)')
     &   'Paramagnetic first order singlet corrections'
         IF (SIGMAP1T) WRITE (LUPRI,'(7X,A)')
     &   'Paramagnetic first order triplet corrections'
         IF (SIGMAP3S) WRITE (LUPRI,'(7X,A)')
     &   'Paramagnetic third order singlet corrections'
         IF (SIGMAP3T) WRITE (LUPRI,'(7X,A)')
     &   'Paramagnetic third order triplet corrections'
      ENDIF
      WRITE(LUPRI,*)' '


      CALL QEXIT('LRSCDRV_SHI')
      RETURN
      END
