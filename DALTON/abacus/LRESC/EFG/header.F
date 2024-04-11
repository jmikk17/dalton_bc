
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
      SUBROUTINE HEADER_LREFG
C Purpose:
C     Print information about LRESC-EFG module              
C Author:
C     J. J. Aucar 2024 (moved from workplace)
c#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
#include "nuclei.h"

      CALL QENTER('HEADER_LREFG')
C
C     Information about calculation of LRESC
C
C ---- Default
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


      CALL QEXIT('HEADER_LREFG')
      RETURN
      END
