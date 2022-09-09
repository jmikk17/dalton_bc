
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
     &('ABACUS - LRESC RELATIVISTIC CORRECTIONS TO SHIELDING','*',124)

      WRITE(lupri,*)
      WRITE(lupri,*)' First and Second order relativistic '//
     & 'corrections to nuclear shielding constant.'
      WRITE(LUPRI,*)

      IF (LRESCALL) THEN
            WRITE(lupri,*)
     &   '* LRESC corrections to paramagnetic component were activated'
      ELSE
            WRITE(lupri,*) '* LRESC corrections were activated'
      ENDIF
      WRITE(LUPRI,*)

c      WRITE(LUPRI,*)'   -  For definitions on terms named ABC,'//
c     &              ' refered to relativistic correction '
c      WRITE(LUPRI,*)'      of sigma, see JCP. 137, 214319(2012). '
      WRITE(LUPRI,*)'   -  Original LRESC papers at : J. Chem. Phys.'//
     &               ' 118, 2 (2003), doi:10.1063/1.1525808'
      WRITE(LUPRI,'(33X,A)')' Mol. Phys. 101, 20, 3103-3109 (2003) '
c     &                       ' DOI: 10.1080/00268970310001617784'
      WRITE(LUPRI,'(33X,A)')' J. Chem. Phys. 125, 064107 (2006), '//
     & 'doi:10.1063/1.2244572'

      WRITE(LUPRI,*)
      WRITE(LUPRI,'(4X,A,I2,A,A)') '-  Selected Atom    : #',LRATOM,
     &         ": ", NAMN(LRATOM)
      WRITE(LUPRI,'(4X,A,3F10.4)') '-  Atom coordinates :',
     &  CORD(1,LRATOM), CORD(2,LRATOM),CORD(3,LRATOM)
      WRITE(LUPRI,'(4X,A)')'-  Remember : '
      WRITE(LUPRI,'(14X,A)')'*  GAUGEO should be placed at selected '//
     &  'nucleus, on INTEGRALS SECTION, in a.u'
      IF (.NOT. PRTALL1) THEN
         WRITE(LUPRI,'(14X,A)')'*  For individual contributions, add '//
     &   '.PRTALL in your input file under *LROPTS'
      ENDIF
      WRITE(LUPRI,'(14X,A)')'*  Third Order Singlet and Triplet '//
     & 'Corrections must be run in order to get'
      WRITE(LUPRI,'(17X,A)')'all relativistic corrections. Not '//
     &                      ' yet implemented, but can be done via'
      WRITE(LUPRI,'(17X,A)')'RESPONSE calculation (see Manual for'//
     &                      ' details)'

C ---Not Default
      IF ((GAUCHANG).OR.(JIMPRT.NE.0).OR.(.NOT. LRESCALL))
     & CALL HEADER('Changes of defaults for LRESC :',0)

      IF (GAUCHANG) THEN
         WRITE(lupri,*)'*** WARNING:'
         WRITE(lupri,*)'GAUGEO must be placed on selected nucleus at'
         WRITE(lupri,*)'a.u. due to Hermite section (1 body integrals)'
      END IF

      IF (JIMPRT.NE.0) THEN
         WRITE (LUPRI,'(A,I5)')
     &      '--- Print level : ', JIMPRT
      END IF

      IF (.NOT. LRESCALL) THEN
         WRITE (LUPRI,'(A)')
     &      '--- Not all corrections were activated. ' //
     &      'Only the following will be done : '
         IF (SIGMAP1S)  WRITE (LUPRI,'(A)')'  LRESC : '//
     &   '   Paramagnetic first order singlet corrections to shielding'
         IF (SIGMAP1T) WRITE (LUPRI,'(A)')'  LRESC : '//
     &   '   Paramagnetic first order triplet corrections to shielding'
         IF (SIGMAP3S) WRITE (LUPRI,'(A)')'  LRESC : '//
     &   '   Paramagnetic third order singlet corrections to shielding'
         IF (SIGMAP3T) WRITE (LUPRI,'(A)')'  LRESC : '//
     &   '   Paramagnetic third order triplet corrections to shielding'
            WRITE (LUPRI,*) 'LRESC : '
         IF (SIGMAD0S) WRITE (LUPRI,*)
     &   '    Diamagnetic zero order singlet corrections to shielding'
         IF (SIGMAD1S) WRITE (LUPRI,*)
     &   '   Diamagnetic first order singlet corrections to shielding'
         IF (NUCSPI .GT. 0)
     &      WRITE (LUPRI,'(A)')' XXXXXXXX NUCSPI GT 0 ??  '//
     &           '    '
      ENDIF
      WRITE(LUPRI,*)' '

      CALL QEXIT('LRSCDRV_SHI')
      RETURN
      END
