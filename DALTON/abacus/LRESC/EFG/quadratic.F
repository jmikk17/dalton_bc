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
      SUBROUTINE QRLREFG(KEY,WORK,LWORK)
      USE lrescinf_module

#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "nuclei.h"
#include "chrnos.h"
#include "chrxyz.h"

      DIMENSION WORK(LWORK)
      CHARACTER*9 :: KEY
      CHARACTER*8 :: operators_trio(3)
      CHARACTER*1 :: spin_trio(3)

      Character*8, dimension(:,:), allocatable :: A_operators
      Character*8, dimension(:,:), allocatable :: B_operators
      Character*8, dimension(:,:), allocatable :: C_operators
      character*8 LISTA1(4*MXCOOR+9)
      
      REAL*8, dimension(:), allocatable :: QR(:),TOTALQR(:)
c      REAL*8 QR(:),TOTALQR(:)
c      REAL*8 QRvector(3)


      INTEGER                  :: NopA, NopB, NopC



      
      CALL QENTER('QRLREFG')
      CALL GETTIM(TIMQR,WALQR)

      allocate(QR(NUCIND))
      allocate(TOTALQR(NUCIND))
      TOTALQR = 0.d0
      QR = 0.d0


      IJ = 1
      DO I=1,NUCIND
         LISTA1(IJ)= 'ZZEFG'//CHRNOS(I/10)//CHRNOS(MOD(I,10))//
     &      CHRNOS(1)
         IJ = IJ + 1
      END DO

      IF (KEY.EQ.'EFG-SO-SO') THEN
         NopA = 1
         NopB = 3
         NopC = 3
         allocate(A_operators(NopA,NUCIND))
         allocate(B_operators(NopB,NUCIND))
         allocate(C_operators(NopC,NUCIND))
         spin_trio(1)='0'
         spin_trio(2)='1'
         spin_trio(3)='1'


         A_operators(1, :) = reshape(LISTA1, (/NUCIND/))

         B_operators(1,:) = 'X1SPNORB'
         B_operators(2,:) = 'Y1SPNORB'
         B_operators(3,:) = 'Z1SPNORB'
         C_operators(1,:) = 'X1SPNORB'
         C_operators(2,:) = 'Y1SPNORB'
         C_operators(3,:) = 'Z1SPNORB'

      END IF


      do j = 1, NUCIND
      do i = 1, 3
         operators_trio(1)=A_operators(1,j)
         operators_trio(2)=B_operators(i,j)
         operators_trio(3)=C_operators(i,j)
         call create_QRinput(operators_trio,spin_trio)
         LUPRI = 6232
         OPEN(UNIT=LUPRI, FILE='LRESC-EFG-QR',status='unknown',
     &              position='append')
         call RSPDRV(WORK,LWORK)
         CLOSE(LUPRI)
         LUPRI  = 2
         CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
         CALL GPCLOSE(LURSP3,'DELETE')
         CALL GPCLOSE(LURSP4,'DELETE')
         CALL GPCLOSE(LURSP5,'DELETE')
         call get_QRvalue(QR(i),operators_trio)
         TOTALQR(j) = TOTALQR(j) + QR(i)
      end do
      end do


      
      IF (KEY.EQ.'EFG-SO-SO') THEN
         EFGC4(:,1) = -0.5D0 * C4EFGSOL * TOTALQR(1:NUCIND)
      END IF


      deallocate(A_operators)
      deallocate(B_operators)
      deallocate(C_operators)
      deallocate(QR)
      deallocate(TOTALQR)


      CALL GETTIM(TEND,WEND)
      TIMQR = TEND - TIMQR
      WALQR = WEND - WALQR
      IF (JIMPRT .GE. 2) THEN
         CALL TIMTXT('    Total CPU  time used in LRESC-EFG-QR('
     &                     //KEY//'):',TIMQR,LUPRI)
         CALL TIMTXT('    Total wall time used in LRESC-EFG-QR('
     &                     //KEY//'):',WALQR,LUPRI)
      END IF
      CALL QEXIT('QRLREFG')

      RETURN

      END
