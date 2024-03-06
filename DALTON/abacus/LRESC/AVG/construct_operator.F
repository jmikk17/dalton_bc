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
C
      SUBROUTINE construct_operator(KEY)
C
C     Purpose:
C     Construct operators matrices as productor of at 
C     least two known ones.
C     Author: Juan J. Aucar         2024

#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inforb.h"
#include "nuclei.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "chrxyz.h"
#include "chrnos.h"
#include "orgcom.h"
      double precision cmos(NORBT,NBAST)
      LOGICAL FOUND
      CHARACTER*4 KEY
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      double precision, dimension(NBAST,NBAST) :: AOPERATOR,BOPERATOR
      double precision, dimension(NBAST,NBAST) :: aux0,aux1,aux2,aux3
      double precision, dimension(NBAST*NORBT) :: xmos
      CHARACTER*8 RTNLBL(2)

c     NBAST: Total number of (atomic) basis functions
c     NORBT: Total number of molecular orbitals
C
C
C   KEY = 'PSOK':   2.0(<KINENERG> . I . <PSO> + <PSO> . I . <KINENERG>)
C
      CALL QENTER('CONSTRUCT OP')
      CALL TIMER('START ',TIMEIN,TIMOUT)
C
      npos1 = 3*LRATOM-2
cb      write(lupri,*)' selected atom is :', LRATOM
      do i = 1, 3
         ICOORD = i !(x,y,z)
         ISCOOR = LRATOM*3 + (ICOORD - 3)
         IFIRST=ISCOOR/100
         ISECND = MOD(ISCOOR,100)/10
         ITHIRD = MOD(MOD(ISCOOR,100),10)
         nuc_coord(i)=CHRNOS(IFIRST)//CHRNOS(ISECND)//CHRNOS(ITHIRD)
      end do

      nuc_number=CHRNOS(LRATOM/100)//CHRNOS(MOD(LRATOM,100)/10)//
     &              CHRNOS(MOD(LRATOM,10))

C      cmos(i,j)   :  Molecular orbitals coeficients matrix elements (from SIRIUS.RST)
       CALL RD_SIRIFC('CMO',FOUND,xmos)
       k = 1
       cmos=0.0
       do i=1,NORBT                !i corresponds to a MO
          do j=1, NBAST           !j corresponds to j coef. of the ith MO
             cmos(i,j) = xmos(k)  ! Transforms vector into matrix
             k = k + 1
          end do
       end do

      aux0=matmul(transpose(cmos),cmos)
C ---------------------------------------------------------------------
      IF(KEY.EQ."PSOK") THEN
         CALL GETDAT(RTNLBL(1),RTNLBL(2))
         RTNLBL(1)='1 A  21 '
C        Replace time information with symmetry information
         RTNLBL(2)='SQUARE  '

         call READAOP(AOPERATOR, NBAST, 'KINENERG') !gets KINENERG matrix
         DO I= 1, 3 !x, y, z
            call READAOP(BOPERATOR, NBAST, 'PSO '//nuc_coord(I)//' ') !gets PSO_i matrix
c           Construct A.I.B + B.I.A
            aux1 = matmul(matmul(AOPERATOR,aux0),BOPERATOR)    !KINENERG . I . PSO
            aux2 = matmul(matmul(BOPERATOR,aux0),AOPERATOR)    !PSO . I . KINENERG
            aux3 = 0.D0 !Initialize/reset matrix
            aux3 = 2.0*(aux1+aux2)
c           aux3 = 2.0*(0.5*aux2) !This corresponds to the built-in PSO-K integrals

            call WRITEAOPJUAN(aux3,NBAST,'PSOKJ'//nuc_coord(I),RTNLBL) !Write Matrix in AOPROPER

c           For debugging purposes:
c           write(lupri,*) "IMPRIMIENDO EL OPERADOR: "
c     &     //'PSOKJ'//nuc_coord(I)
c           call OUTPUT(aux3,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)  ! prints aux3
         END DO
      ENDIF
C ---------------------------------------------------------------------

      CALL TIMER ('CONSTRUCT OP',TIMEIN,TIMOUT)
      CALL QEXIT('CONSTRUCT OP')
      RETURN
      END
