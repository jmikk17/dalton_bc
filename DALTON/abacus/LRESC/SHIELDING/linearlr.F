
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
C  /* Deck get_lr */
      SUBROUTINE GET_LR(KEY,TRIPLETE)

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

      REAL*8 LR(3,3),TOTALLR(3,3)
      REAL*8 LRtensor(3,3)
      CHARACTER*4 :: KEY
      Character*8, dimension(:,:,:), allocatable :: A_operators
      Character*8, dimension(:,:,:), allocatable :: B_operators

      INTEGER                  :: NopA, NopB
      CHARACTER*8, dimension(2) :: operators_duo
      INTEGER :: SUMA_B
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      LOGICAL :: TRIPLETE,TRPLET_FLAG

      CALL QENTER('LRSLINEAR')

c     NopA: number of A operators for each shield tensor element
c     NopB: number of B operators for each shield tensor element

C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'DIAM'  DIAMAGNETIC <<DIA,MV>>
C         'DIAD'  DIAMAGNETIC <<DIA,DW>>
C         'PSOK'  PARAMAGNETIC with built-in integrals and AVG integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'OZK '  PARAMAGNETIC
C         'FCBS' 'SDBS'
C         'FCZK' 'FCZK'
C         'FCDE' 'FCDE'
      IF (JIMPRT .GE. 1) THEN
         IF (KEY.EQ.'DIAM') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic Second Order Singlet (DIAM)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'DIAD') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic Second Order Singlet (DIAD)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'NRPA') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic non-relativistic contribution',
     &   ' (NRPA) is under calculation'
         IF (KEY.EQ.'PSOK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (PSOK)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'PSKI') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (PSKI)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'OZK ') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (OZK)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'FCZK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCZK)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'FCDE') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCDE)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'SDK ') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDK)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'SDD ') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDD)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'FCBS') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCBS)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'SDBS') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDBS)',
     &   ' correction is under calculation'
C         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
      END IF
      TOTALLR = 0.d0
      LR = 0.d0


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

      IF((KEY.EQ."DIAM").OR.(KEY.EQ."DIAD")) THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                    ij
         A_operators(1,1,1) = nuc_coord(1)//'NSCO'//CHRXYZ(1) !11
         A_operators(1,2,1) = nuc_coord(1)//'NSCO'//CHRXYZ(2) !12
         A_operators(1,3,1) = nuc_coord(1)//'NSCO'//CHRXYZ(3) !13
         A_operators(2,1,1) = nuc_coord(2)//'NSCO'//CHRXYZ(1) !21
         A_operators(2,2,1) = nuc_coord(2)//'NSCO'//CHRXYZ(2) !22
         A_operators(2,3,1) = nuc_coord(2)//'NSCO'//CHRXYZ(3) !23
         A_operators(3,1,1) = nuc_coord(3)//'NSCO'//CHRXYZ(1) !31
         A_operators(3,2,1) = nuc_coord(3)//'NSCO'//CHRXYZ(2) !32
         A_operators(3,3,1) = nuc_coord(3)//'NSCO'//CHRXYZ(3) !33
C        -----------------------------------
         allocate(B_operators(3,3,NopB))
         IF (KEY.EQ."DIAD") B_operators(:,:,1) ='DARWIN  '
         IF (KEY.EQ."DIAM") B_operators(:,:,1) ='MASSVELO'
      END IF


      IF(KEY.EQ."NRPA") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
C        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                      ij
         A_operators(1,1,1) = 'PSO '//nuc_coord(1)//' '       !11
         A_operators(1,2,1) = 'PSO '//nuc_coord(1)//' '       !12
         A_operators(1,3,1) = 'PSO '//nuc_coord(1)//' '       !13
         A_operators(2,1,1) = 'PSO '//nuc_coord(2)//' '       !21
         A_operators(2,2,1) = 'PSO '//nuc_coord(2)//' '       !22
         A_operators(2,3,1) = 'PSO '//nuc_coord(2)//' '       !23
         A_operators(3,1,1) = 'PSO '//nuc_coord(3)//' '       !31
         A_operators(3,2,1) = 'PSO '//nuc_coord(3)//' '       !32
         A_operators(3,3,1) = 'PSO '//nuc_coord(3)//' '       !33
C        -----------------------------------
         allocate(B_operators(3,3,NopB))
         B_operators(1,1,1) = 'XANGMOM '                      !11
         B_operators(1,2,1) = 'YANGMOM '                      !12
         B_operators(1,3,1) = 'ZANGMOM '                      !13
         B_operators(2,1,1) = 'XANGMOM '                      !21
         B_operators(2,2,1) = 'YANGMOM '                      !22
         B_operators(2,3,1) = 'ZANGMOM '                      !23
         B_operators(3,1,1) = 'XANGMOM '                      !31
         B_operators(3,2,1) = 'YANGMOM '                      !32
         B_operators(3,3,1) = 'ZANGMOM '                      !33
      END IF


      IF(KEY.EQ."PSKK") THEN !PSOK with abacus integrals
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
c        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                    ij
         A_operators(1,1,1) = 'PSOKE'//nuc_coord(1)           !11
         A_operators(1,2,1) = 'PSOKE'//nuc_coord(1)           !12
         A_operators(1,3,1) = 'PSOKE'//nuc_coord(1)           !13
         A_operators(2,1,1) = 'PSOKE'//nuc_coord(2)           !21
         A_operators(2,2,1) = 'PSOKE'//nuc_coord(2)           !22
         A_operators(2,3,1) = 'PSOKE'//nuc_coord(2)           !23
         A_operators(3,1,1) = 'PSOKE'//nuc_coord(3)           !31
         A_operators(3,2,1) = 'PSOKE'//nuc_coord(3)           !32
         A_operators(3,3,1) = 'PSOKE'//nuc_coord(3)           !33
c        -----------------------------------
         allocate(B_operators(3,3,NopB))
         B_operators(1,1,1) = 'XANGMOM '                      !11
         B_operators(1,2,1) = 'YANGMOM '                      !12
         B_operators(1,3,1) = 'ZANGMOM '                      !13
         B_operators(2,1,1) = 'XANGMOM '                      !21
         B_operators(2,2,1) = 'YANGMOM '                      !22
         B_operators(2,3,1) = 'ZANGMOM '                      !23
         B_operators(3,1,1) = 'XANGMOM '                      !31
         B_operators(3,2,1) = 'YANGMOM '                      !32
         B_operators(3,3,1) = 'ZANGMOM '                      !33
      END IF

      IF(KEY.EQ."PSOK") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                    ij
         A_operators(1,1,1) = 'PSOKJ'//nuc_coord(1)           !11
         A_operators(1,2,1) = 'PSOKJ'//nuc_coord(1)           !12
         A_operators(1,3,1) = 'PSOKJ'//nuc_coord(1)           !13
         A_operators(2,1,1) = 'PSOKJ'//nuc_coord(2)           !21
         A_operators(2,2,1) = 'PSOKJ'//nuc_coord(2)           !22
         A_operators(2,3,1) = 'PSOKJ'//nuc_coord(2)           !23
         A_operators(3,1,1) = 'PSOKJ'//nuc_coord(3)           !31
         A_operators(3,2,1) = 'PSOKJ'//nuc_coord(3)           !32
         A_operators(3,3,1) = 'PSOKJ'//nuc_coord(3)           !33
C        -----------------------------------
         allocate(B_operators(3,3,NopB))
         B_operators(1,1,1) = 'XANGMOM '                      !11
         B_operators(1,2,1) = 'YANGMOM '                      !12
         B_operators(1,3,1) = 'ZANGMOM '                      !13
         B_operators(2,1,1) = 'XANGMOM '                      !21
         B_operators(2,2,1) = 'YANGMOM '                      !22
         B_operators(2,3,1) = 'ZANGMOM '                      !23
         B_operators(3,1,1) = 'XANGMOM '                      !31
         B_operators(3,2,1) = 'YANGMOM '                      !32
         B_operators(3,3,1) = 'ZANGMOM '                      !33
      END IF


      IF(KEY.EQ."OZK ") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                    ij
         A_operators(1,1,1) = 'PSO '//nuc_coord(1)//' '       !11
         A_operators(1,2,1) = 'PSO '//nuc_coord(1)//' '       !12
         A_operators(1,3,1) = 'PSO '//nuc_coord(1)//' '       !13
         A_operators(2,1,1) = 'PSO '//nuc_coord(2)//' '       !21
         A_operators(2,2,1) = 'PSO '//nuc_coord(2)//' '       !22
         A_operators(2,3,1) = 'PSO '//nuc_coord(2)//' '       !23
         A_operators(3,1,1) = 'PSO '//nuc_coord(3)//' '       !31
         A_operators(3,2,1) = 'PSO '//nuc_coord(3)//' '       !32
         A_operators(3,3,1) = 'PSO '//nuc_coord(3)//' '       !33
C        -----------------------------------
         allocate(B_operators(3,3,NopB))
         B_operators(1,1,1) = CHRXYZ(1)//'OZKE   '            !11
         B_operators(1,2,1) = CHRXYZ(2)//'OZKE   '            !12
         B_operators(1,3,1) = CHRXYZ(3)//'OZKE   '            !13
         B_operators(2,1,1) = CHRXYZ(1)//'OZKE   '            !21
         B_operators(2,2,1) = CHRXYZ(2)//'OZKE   '            !22
         B_operators(2,3,1) = CHRXYZ(3)//'OZKE   '            !23
         B_operators(3,1,1) = CHRXYZ(1)//'OZKE   '            !31
         B_operators(3,2,1) = CHRXYZ(2)//'OZKE   '            !32
         B_operators(3,3,1) = CHRXYZ(3)//'OZKE   '            !33
      END IF

      IF (KEY.EQ."FCZK") THEN
         NopA = 1
         TRPLET_FLAG = .TRUE.
C        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                          ij
         A_operators(:,:,:) = ' '
         A_operators(1,1,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number  !11
         A_operators(2,2,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number  !22
         A_operators(3,3,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number  !33
C        -----------------------------------
         NopB = 1
         allocate(B_operators(3,3,NopB))
         IF(KEY.EQ."FCZK") THEN
            B_operators(:,:,:) = ' '
            B_operators(1,1,1) = 'KINENERG'                      !11
            B_operators(2,2,1) = 'KINENERG'                      !22
            B_operators(3,3,1) = 'KINENERG'                      !33
         END IF

      END IF


      IF(KEY.EQ."FCDE") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .TRUE.
C        -----------------------------------
         allocate(A_operators(3,3,NopA))
c                  3*i+j-3                                          ij
         A_operators(:,:,:) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
C        -----------------------------------
         allocate(B_operators(3,3,NopB))
         B_operators(1,1,1) = 'dd/dxdx '
         B_operators(1,2,1) = 'dd/dxdy '
         B_operators(1,3,1) = 'dd/dxdz '
C        ----------------------------------
         B_operators(2,1,1) = 'dd/dxdy '
         B_operators(2,2,1) = 'dd/dydy '
         B_operators(2,3,1) = 'dd/dydz '
C        ----------------------------------
         B_operators(3,1,1) = 'dd/dxdz '
         B_operators(3,2,1) = 'dd/dydz '
         B_operators(3,3,1) = 'dd/dzdz '
      END IF

      IF (KEY.EQ."SDK ") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .TRUE.
         allocate(A_operators(3,3,NopA))
         allocate(B_operators(3,3,NopB))
C        -----------------------------------
         A_operators(1,1,1) = 'SD '//nuc_coord(1)//' x'
         A_operators(1,2,1) = 'SD '//nuc_coord(1)//' y'
         A_operators(1,3,1) = 'SD '//nuc_coord(1)//' z'
         A_operators(2,1,1) = 'SD '//nuc_coord(2)//' x'
         A_operators(2,2,1) = 'SD '//nuc_coord(2)//' y'
         A_operators(2,3,1) = 'SD '//nuc_coord(2)//' z'
         A_operators(3,1,1) = 'SD '//nuc_coord(3)//' x'
         A_operators(3,2,1) = 'SD '//nuc_coord(3)//' y'
         A_operators(3,3,1) = 'SD '//nuc_coord(3)//' z'
C        -----------------------------------
         B_operators(:,:,:) = 'KINENERG'
C        -----------------------------------
      END IF


      IF (KEY.EQ."SDDE") THEN  !OLD - never called - just in case
         NopA = 3
         NopB = 3
         TRPLET_FLAG = .TRUE.
         allocate(A_operators(3,3,NopA))
         allocate(B_operators(3,3,NopB))
C        ----------------------------------- 
c                      3*i+j-3                                       ij
         DO K=1, 3
            A_operators(1,1,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !11
            A_operators(1,2,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !12
            A_operators(1,3,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !13
            A_operators(2,1,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !21
            A_operators(2,2,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !22
            A_operators(2,3,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !23
            A_operators(3,1,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !31
            A_operators(3,2,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !32
            A_operators(3,3,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !33
         END DO
C        -----------------------------------
c        K=1
         B_operators(1,1,1) = 'dd/dxdx ' !11
         B_operators(2,1,1) = 'dd/dxdx ' !21
         B_operators(3,1,1) = 'dd/dxdx ' !31
         B_operators(1,2,1) = 'dd/dxdy ' !12
         B_operators(2,2,1) = 'dd/dxdy ' !22
         B_operators(3,2,1) = 'dd/dxdy ' !32
         B_operators(1,3,1) = 'dd/dxdz ' !13
         B_operators(2,3,1) = 'dd/dxdz ' !23
         B_operators(3,3,1) = 'dd/dxdz ' !33
C        -----------------------------------
c        K=2
         B_operators(1,1,2) = 'dd/dxdy ' !11
         B_operators(2,1,2) = 'dd/dxdy ' !21
         B_operators(3,1,2) = 'dd/dxdy ' !31
         B_operators(1,2,2) = 'dd/dydy ' !12
         B_operators(2,2,2) = 'dd/dydy ' !22
         B_operators(3,2,2) = 'dd/dydy ' !32
         B_operators(1,3,2) = 'dd/dydz ' !13
         B_operators(2,3,2) = 'dd/dydz ' !23
         B_operators(3,3,2) = 'dd/dydz ' !33
C        -----------------------------------
c        K=3
         B_operators(1,1,3) = 'dd/dxdz ' !11
         B_operators(2,1,3) = 'dd/dxdz ' !21
         B_operators(3,1,3) = 'dd/dxdz ' !31
         B_operators(1,2,3) = 'dd/dydz ' !12
         B_operators(2,2,3) = 'dd/dydz ' !22
         B_operators(3,2,3) = 'dd/dydz ' !32
         B_operators(1,3,3) = 'dd/dzdz ' !13
         B_operators(2,3,3) = 'dd/dzdz ' !23
         B_operators(3,3,3) = 'dd/dzdz ' !33
C        ----------------------------------
      END IF

      IF (KEY.EQ."SDD ") THEN
         NopA = 3
         NopB = 3
         TRPLET_FLAG = .TRUE.
         allocate(A_operators(3,3,NopA))
         allocate(B_operators(3,3,NopB))
C        ----------------------------------- 
         A_operators(1,1,1) = 'dd/dxdx '
         A_operators(1,2,1) = 'dd/dxdy '
         A_operators(1,3,1) = 'dd/dxdz '
         A_operators(2,1,1) = 'dd/dxdx '
         A_operators(2,2,1) = 'dd/dxdy '
         A_operators(2,3,1) = 'dd/dxdz '
         A_operators(3,1,1) = 'dd/dxdx '
         A_operators(3,2,1) = 'dd/dxdy '
         A_operators(3,3,1) = 'dd/dxdz '

         A_operators(1,1,2) = 'dd/dxdy '
         A_operators(1,2,2) = 'dd/dydy '
         A_operators(1,3,2) = 'dd/dydz '
         A_operators(2,1,2) = 'dd/dxdy '
         A_operators(2,2,2) = 'dd/dydy '
         A_operators(2,3,2) = 'dd/dydz '
         A_operators(3,1,2) = 'dd/dxdy '
         A_operators(3,2,2) = 'dd/dydy '
         A_operators(3,3,2) = 'dd/dydz '

         A_operators(1,1,3) = 'dd/dxdz '
         A_operators(1,2,3) = 'dd/dydz '
         A_operators(1,3,3) = 'dd/dzdz '
         A_operators(2,1,3) = 'dd/dxdz '
         A_operators(2,2,3) = 'dd/dydz '
         A_operators(2,3,3) = 'dd/dzdz '
         A_operators(3,1,3) = 'dd/dxdz '
         A_operators(3,2,3) = 'dd/dydz '
         A_operators(3,3,3) = 'dd/dzdz '
C        -----------------------------------
         DO K=1, 3
            B_operators(1,1,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !11
            B_operators(1,2,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !12
            B_operators(1,3,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !13
            B_operators(2,1,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !21
            B_operators(2,2,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !22
            B_operators(2,3,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !23
            B_operators(3,1,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !31
            B_operators(3,2,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !32
            B_operators(3,3,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !33
         END DO
C        -----------------------------------
      END IF

      IF (KEY.EQ."FCBS") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .TRUE.
         allocate(A_operators(3,3,NopA))
         allocate(B_operators(3,3,NopB))
C        -----------------------------------
c                  3*i+j-3                               ij
C        -----------------------------------
         A_operators(:,:,:) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
C        -----------------------------------
         B_operators(1,1,1) = 'SOMF  XX' !11
         B_operators(1,2,1) = 'SOMF  XY' !12
         B_operators(1,3,1) = 'SOMF  XZ' !13
         B_operators(2,1,1) = 'SOMF  YX' !21
         B_operators(2,2,1) = 'SOMF  YY' !22
         B_operators(2,3,1) = 'SOMF  YZ' !23
         B_operators(3,1,1) = 'SOMF  ZX' !31
         B_operators(3,2,1) = 'SOMF  ZY' !32
         B_operators(3,3,1) = 'SOMF  ZZ' !33
C        -----------------------------------
      END IF

         IF(KEY.EQ."SDBS") THEN
         NopA = 3
         NopB = 3
         TRPLET_FLAG = .TRUE.
         allocate(A_operators(3,3,NopA))
         allocate(B_operators(3,3,NopB))
C        -----------------------------------
         DO K=1, 3
            A_operators(1,1,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !11
            A_operators(1,2,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !12
            A_operators(1,3,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !13
            A_operators(2,1,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !21
            A_operators(2,2,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !22
            A_operators(2,3,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !23
            A_operators(3,1,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !31
            A_operators(3,2,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !32
            A_operators(3,3,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !33
         END DO
C        -----------------------------------
         DO K=1, 3
            B_operators(1,1,K) = 'SOMF  '//CHRXYZ(K)//'X' !11
            B_operators(1,2,K) = 'SOMF  '//CHRXYZ(K)//'Y' !12
            B_operators(1,3,K) = 'SOMF  '//CHRXYZ(K)//'Z' !13
            B_operators(2,1,K) = 'SOMF  '//CHRXYZ(K)//'X' !21
            B_operators(2,2,K) = 'SOMF  '//CHRXYZ(K)//'Y' !22
            B_operators(2,3,K) = 'SOMF  '//CHRXYZ(K)//'Z' !23
            B_operators(3,1,K) = 'SOMF  '//CHRXYZ(K)//'X' !31
            B_operators(3,2,K) = 'SOMF  '//CHRXYZ(K)//'Y' !32
            B_operators(3,3,K) = 'SOMF  '//CHRXYZ(K)//'Z' !33
         END DO
C        -----------------------------------
         END IF



c     Elemento ij del tensor de shielding
      do i = 1, 3 
      do j = 1, 3
         do k = 1, NopA
           IF ((NopA.Eq.1).or.(NopB.Eq.1)) THEN !Combine all A with all B
            do l = 1, NopB
               operators_duo(1) = A_operators(i,j,k)
               operators_duo(2) = B_operators(i,j,l)

                  call get_LRvalue(LR(i,j),operators_duo,TRIPLETE) !LR from here has an extra minus sign
               TOTALLR(i,j) = TOTALLR(i,j) - LR(i,j) !LR sign is corrected here
c               IF (JIMPRT .GE. 3) write(lupri,*) operators_duo(1),
c     &         ' ',operators_duo(2)," LR: ",LR(i,j)
            end do
         else
            l = k
            operators_duo(1) = A_operators(i,j,k)
            operators_duo(2) = B_operators(i,j,l)

            call get_LRvalue(LR(i,j),operators_duo,TRIPLETE) !LR from here has an extra minus sign
            TOTALLR(i,j) = TOTALLR(i,j) - LR(i,j) !LR sign is corrected here
c         IF (JIMPRT .GE. 3) write(lupri,*) 'Labels:', operators_duo(1),
c     &      ' ',operators_duo(2)," LR:",LR(i,j)
         end if
         end do

         LRtensor(i,j) = TOTALLR(i,j)
         IF (JIMPRT .GE. 3) THEN
c            write(lupri,*) KEY,"(",i,",",j,"):",LRtensor(i,j)
c            write(lupri,*)
            IF (JIMPRT .GE. 3) write(lupri,*) '    Labels: ', 
     &       operators_duo(1),' ',operators_duo(2),". LR:",LR(i,j)
         END IF

      end do
      end do

      IF (KEY.EQ.'DIAM') LRDIAM = LRtensor
      IF (KEY.EQ.'DIAD') LRDIAD = LRtensor
      IF (KEY.EQ.'PSOK') LRPSOK = LRtensor
      IF (KEY.EQ.'PSKK') LRPSOKaba = LRtensor
      IF (KEY.EQ.'NRPA') LRPNR = LRtensor
      IF (KEY.EQ.'OZK ') LROZK  = LRtensor
      IF (KEY.EQ.'FCZK') LRFCZK = LRtensor
      IF (KEY.EQ.'FCDE') LRFCDE = LRtensor
      IF (KEY.EQ.'FCBS') LRFCBS = LRtensor
      IF (KEY.EQ.'SDBS') LRSDBS = LRtensor
      IF (KEY.EQ.'SDK ') LRSDK = LRtensor
      IF (KEY.EQ.'SDD ') LRSDD = LRtensor

      deallocate(A_operators)
      deallocate(B_operators)

      CALL QEXIT('LRSLINEAR')

      RETURN

      END 


C  /* Deck newlrinput_all */
      SUBROUTINE LRINPUT_ALL(LRATOMS,NNLRATOMS,TRIPLETE)

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

      CHARACTER*8, allocatable :: operators(:)
      INTEGER                  :: NopA, NopB, NLRoperators
      INTEGER :: SUMA_B, NNLRATOMS
      INTEGER, intent(in)  :: LRATOMS(NNLRATOMS)
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      character*10 :: namefile
      INTEGER :: fileID
      LOGICAL :: TRIPLETE

      CALL QENTER('NEWLRINPUT_SINGLETS_ALL')

c     NopA: number of A operators for each shield tensor element
c     NopB: number of B operators for each shield tensor element

C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'DIAM'  DIAMAGNETIC <<DIA,MV>>
C         'DIAD'  DIAMAGNETIC <<DIA,DW>>
C         'PSOK'  PARAMAGNETIC with built-in integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'OZK '  PARAMAGNETIC OZ-K
C         'FCBSO' 'SDBSO'
C         'FCZK' 'FCZK'
C         'FCDE' 'FCDE'

      namefile='DALTON.INP'
      fileID = -1
      CALL GPOPEN(fileID,namefile,'UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      write(fileID,'(a)') '**DALTON INPUT'
      write(fileID,'(a)') '**RESPONSE'
      write(fileID,'(a)') '*LINEAR'
      IF (TRIPLETE) write(fileID,'(a)') '.TRIPLET'

      DO J = 1, size(LRATOMS)
      LRATOM = LRATOMS(J)


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

      NLRoperators = 0

      IF (.not.(TRIPLETE)) THEN
         IF (SIGMAPNR) NLRoperators = NLRoperators + 6
         IF (SIGMAD1S) NLRoperators = NLRoperators + 11
         IF (SIGMAP1S) NLRoperators = NLRoperators + 12
         IF ((SIGMAP1S).AND.(EXTRAINFO)) NLRoperators = NLRoperators + 3
      END IF
      IF ((TRIPLETE).AND.(SIGMAP1T)) NLRoperators = NLRoperators + 26

      allocate(operators(NLRoperators))
      el = 1


      IF (.not.(TRIPLETE)) THEN
      IF (SIGMAPNR) THEN ! Non-relativistic paramagnetic contribution
C        ---------------------------------------------------
         operators(el)   = 'XANGMOM '
         operators(el+1) = 'YANGMOM '
         operators(el+2) = 'ZANGMOM '
C        ---------------------------------------------------
         operators(el+3) = 'PSO '//nuc_coord(1)//' '
         operators(el+4) = 'PSO '//nuc_coord(2)//' '
         operators(el+5) = 'PSO '//nuc_coord(3)//' '
C        ---------------------------------------------------
         el = el + 6
      END IF

      IF (SIGMAD1S) THEN !SIGMAD1S: DIAM and DIAD
C        ---------------------------------------------------
         operators(el)   = nuc_coord(1)//'NSCO'//CHRXYZ(1)
         operators(el+1) = nuc_coord(1)//'NSCO'//CHRXYZ(2)
         operators(el+2) = nuc_coord(1)//'NSCO'//CHRXYZ(3)
         operators(el+3) = nuc_coord(2)//'NSCO'//CHRXYZ(1)
         operators(el+4) = nuc_coord(2)//'NSCO'//CHRXYZ(2)
         operators(el+5) = nuc_coord(2)//'NSCO'//CHRXYZ(3)
         operators(el+6) = nuc_coord(3)//'NSCO'//CHRXYZ(1)
         operators(el+7) = nuc_coord(3)//'NSCO'//CHRXYZ(2)
         operators(el+8) = nuc_coord(3)//'NSCO'//CHRXYZ(3)
C        ---------------------------------------------------
C        ---------------------------------------------------
C        ---------------------------------------------------
         operators(el+9)  = 'DARWIN  '
C        ---------------------------------------------------
         operators(el+10) = 'MASSVELO'
C        ---------------------------------------------------
         el = el + 11
      END IF

      IF (SIGMAP1S) THEN !SIGMAP1S: PSOK and OZK
C        ---------------------------------------------------
         operators(el) = 'XANGMOM '
         operators(el+1) = 'YANGMOM '
         operators(el+2) = 'ZANGMOM '
C        ---------------------------------------------------
C        ---------------------------------------------------
         operators(el+3) = 'PSO '//nuc_coord(1)//' '
         operators(el+4) = 'PSO '//nuc_coord(2)//' '
         operators(el+5) = 'PSO '//nuc_coord(3)//' '
C        ---------------------------------------------------
         operators(el+6)  = CHRXYZ(1)//'OZKE   '
         operators(el+7) = CHRXYZ(2)//'OZKE   '
         operators(el+8) = CHRXYZ(3)//'OZKE   '
C        ---------------------------------------------------
         operators(el+9)   = 'PSOKJ'//nuc_coord(1)
         operators(el+10) = 'PSOKJ'//nuc_coord(2)
         operators(el+11) = 'PSOKJ'//nuc_coord(3)
         el = el + 12
C        ---------------------------------------------------
         IF (EXTRAINFO) THEN
            operators(el)   = 'PSOKE'//nuc_coord(1)
            operators(el+1) = 'PSOKE'//nuc_coord(2)
            operators(el+2) = 'PSOKE'//nuc_coord(3)
            el = el + 3
         END IF
C        ---------------------------------------------------
      END IF
      END IF

      IF (TRIPLETE) THEN
      IF (SIGMAP1T) THEN !FCZK, FCDE, FCBS, SDK1, SDD and SDBS
C        ---------------------------------------------------
         operators(el)   = 'FC '//NAMN(LRATOM)(1:2)//nuc_number !FCZK & FCBS
C        ---------------------------------------------------
         operators(el+1) = 'SD '//nuc_coord(1)//' x'            !SDK1 & SDBS
         operators(el+2) = 'SD '//nuc_coord(1)//' y'
         operators(el+3) = 'SD '//nuc_coord(1)//' z'
         operators(el+4) = 'SD '//nuc_coord(2)//' x'
         operators(el+5) = 'SD '//nuc_coord(2)//' y'
         operators(el+6) = 'SD '//nuc_coord(2)//' z'
         operators(el+7) = 'SD '//nuc_coord(3)//' x'
         operators(el+8) = 'SD '//nuc_coord(3)//' y'
         operators(el+9) = 'SD '//nuc_coord(3)//' z'
C        ---------------------------------------------------
C        ---------------------------------------------------
         operators(el+10) = 'KINENERG'                          !FCZK
C        ---------------------------------------------------
         operators(el+11) = 'dd/dxdx '                          !FCDE,SDK1 & SDD
         operators(el+12) = 'dd/dxdy '
         operators(el+13) = 'dd/dxdz '
         operators(el+14) = 'dd/dydy '
         operators(el+15) = 'dd/dydz '
         operators(el+16) = 'dd/dzdz '
C        ---------------------------------------------------
         operators(el+17) = 'SOMF  XX'                          !FCBS & SDBS
         operators(el+18) = 'SOMF  XY'
         operators(el+19) = 'SOMF  XZ'
         operators(el+20) = 'SOMF  YX'
         operators(el+21) = 'SOMF  YY'
         operators(el+22) = 'SOMF  YZ'
         operators(el+23) = 'SOMF  ZX'
         operators(el+24) = 'SOMF  ZY'
         operators(el+25) = 'SOMF  ZZ'
C        --------------------------------------------------- !NEW SDZK &SDD1
c         operators(el+26) = 'XYEFG'//CHRNOS(LRATOM/10)//
c     &   CHRNOS(MOD(LRATOM,10))//CHRNOS(1)
c         operators(el+27) = 'XZEFG'//CHRNOS(LRATOM/10)//
c     &   CHRNOS(MOD(LRATOM,10))//CHRNOS(1)
c         operators(el+28) = 'YZEFG'//CHRNOS(LRATOM/10)//
c     &   CHRNOS(MOD(LRATOM,10))//CHRNOS(1)
C        ---------------------------------------------------
         el = el + 26
      END IF
      END IF

      do i = 1, SIZE(operators,DIM=1)
         write(fileID,'(a)') '.PROPRT'
         write(fileID,'(a)') operators(i)
      end do


      deallocate(operators)
      END DO

      write(fileID,'(a)') '**END OF DALTON INPUT'
      CALL GPCLOSE(fileID,'KEEP')
      



c      CALL GETTIM(TEND,WEND)
c      TIMHER = TEND - TIMHER
c      WALHER = WEND - WALHER
c      CALL TIMTXT(' Total CPU  time used in '//KEY//' '//NAMN(LRATOM)//
c     &            ' LR_addopnames:',TIMHER,LUPRI)
c      CALL TIMTXT(' Total wall time used in '//KEY//' '//NAMN(LRATOM)//
c     &            ' LR_addopnames:',WALHER,LUPRI)
      CALL QEXIT('NEWLRINPUT_SINGLETS_ALL')

      RETURN

      END 



      subroutine get_LRvalue(LRVAL,labels_in,TRIPLETE)
      USE lrescinf_module
      logical :: file_exists, file_open, TRIPLETE
      INTEGER :: fileID
      CHARACTER*164 STRING, TMPSTRING
      CHARACTER*8, intent(in) :: labels_in(2)
      CHARACTER*8:: labels_read(2)
      CHARACTER*10 :: LRfile
      REAL*8, intent(out) :: LRVAL

c

#include "priunit.h"
#include "chrnos.h"      

C Based on get_QRvalue
      IF (TRIPLETE) THEN
         LRfile = 'LRESC-LR-T'
      ELSE
         LRfile = 'LRESC-LR-S'
      END IF

      IF (.not.(labels_in(1).EQ.' ')) THEN
      IF (.not.(labels_in(2).EQ.' ')) THEN

      inquire(file = LRfile, exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            fileID = -1
            CALL GPOPEN(fileID,LRfile,'OLD',' ',
     &                  'FORMATTED',IDUMMY,.FALSE.)
         end if
         rewind(fileID)
      else
         call quit(LRfile//' file not found')
      end if

 1    READ (fileID,'(A164)',END=3) STRING
c     Getting LR when its called for the first time
      IF (STRING(1:30) .EQ. '@ FREQUENCY INDEPENDENT SECOND') THEN
         READ (fileID,'(A164)') TMPSTRING !Empty line
 4       READ (fileID,'(A164)') TMPSTRING
         IF ( (TMPSTRING(7:14) .EQ. labels_in(1) .AND.
     &       TMPSTRING(18:25) .EQ. labels_in(2)) .OR.
     &       (TMPSTRING(7:14) .EQ. labels_in(2) .AND.
     &       TMPSTRING(18:25) .EQ. labels_in(1)) ) THEN
                BACKSPACE(fileID)
                READ (fileID,'(31X,F19.12)') LRVAL
         ELSE
            GOTO 4
         END IF
      ELSE
         GOTO 2  
      END IF

c     Getting LR when its not called for the first time
 2    IF (STRING(2:25) .EQ. 'Linear response function') THEN
          READ (fileID,'(36X,A8)') labels_read(1)
          READ (fileID,'(36X,A8)') labels_read(2)
      IF (labels_read(1) .EQ. labels_in(1) .AND.
     &    labels_read(2) .EQ. labels_in(2) ) then
         READ (fileID,'(/,48X,F15.12)') LRVAL
      ELSE
         GOTO 1
      END IF

      ELSE
         GOTO 1
      END IF

 3           CONTINUE

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

      ELSE
         LRVAL = 0.d0
      END IF
      END IF


      RETURN
      END
