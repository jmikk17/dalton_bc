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

      SUBROUTINE AVGchecks(PRT,CMOS,OVERLP)
C
C       Linux version, Melo, azua, giribet, alejandrito, 
C       Aucar, Danian, romero, homero,bart y lisa
C
C       Revisited - J. Aucar 2020
#include "inftap.h"
#include "inforb.h"
#include "priunit.h"

      integer, intent(in) :: PRT
      double precision, dimension(NORBT,NBAST), intent(in) :: CMOS
      double precision, dimension(NBAST,NBAST), intent(in) :: OVERLP
      double precision naos, nmos, Av, uno
      integer*8 k
      double precision Aux(NBAST,NBAST), P(NBAST,NBAST)
      
       naos = NBAST ! NBAST: Total number of (atomic) basis functions
       nmos= NORBT !NORBT: Total number of molecular orbitals
C ----------------------------------------------------------
C         BASE CHECKING
C ----------------------------------------------------------
      IF (PRT.GT.0) THEN 
         write(lupri,*) '*******************************************' 
         write(lupri,*) '********** Basis set check    *************'
         write(lupri,*) '*******************************************' 
         IF (PRT.GT.3) then
            write(lupri,*) ' **** INTEGRALS OUTPUT   DETECTED ! *******'
            write(lupri,*) ' ************** CMOS care !! **************'
            write(lupri,*) ' ***** DALTON.OUT has the transposed **'
            write(lupri,*) ' ******      CMOS(NOCCT,naos)       ******'
            write(lupri,*) ' *****************************************' 


         write(lupri,*) 'CMOS (checks.F):' 
         call OUTPUT(CMOS,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)  ! prints CMOS
         write(lupri,*) 'OVERLAP (checks.F):' 
         call OUTPUT(OVERLP,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints OVERLP
         ENDIF
C        ---------------------------------------------
C          CHECK  CMOS has norm = 1 in atomic basis !!
C        ---------------------------------------------
         uno = 0.0
         Av = 0.0
         write(lupri,*) 
         write(lupri,*) 'Norm of  CMOS by rows in atomic basis :' 
         do i = 1,NOCCT
            do j = 1,naos
               uno = uno + CMOS(i,j)*CMOS(i,j)
               uno = sqrt(uno)
            enddo
            write(lupri,*) '  Orbital', i, ' : ',uno
            Av = Av + uno 
         enddo
         write(lupri,*) '** Sum of row norms :', Av
         write(lupri,*) 
C     ------------------------
C     CHEQUEO Tr(PS) Mulliken 
C     ------------------------
         Av = 0.0
         P = 0.0
         Aux = 0.0
         do i = 1,naos
            do j = 1,naos
               do k=1,NOCCT
                  P(i,j) = P(i,j)+CMOS(k,i)*CMOS(k,j)
               enddo
            enddo
         enddo
         do i = 1,naos
            do j = 1,naos
               do k=1,naos
               Aux(i,j) = Aux(i,j) + P(i,k)*OVERLP(k,j)
               enddo
            enddo
         enddo
         do k=1,naos
            Av = Av + Aux(k,k)
         enddo
         write(lupri,*) '  Tr(P.S) = ',Av
      ENDIF
      RETURN
      END SUBROUTINE
