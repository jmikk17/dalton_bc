
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

#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
#include "nuclei.h"
#include "chrnos.h"
#include "chrxyz.h"

      REAL*8 LR(3,3),TOTALLR(3,3)
      REAL*8 LRtensor(3,3)
      CHARACTER*4 :: KEY
      CHARACTER*8, allocatable :: operadores_A(:,:)
      CHARACTER*8, allocatable :: operadores_B(:,:)
      INTEGER                  :: NopA, NopB
      CHARACTER*8, dimension(2) :: operadores_duo
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
C         'PSOK'  PARAMAGNETIC with built-in integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'LKIN'  PARAMAGNETIC
C         'FCBS' 'SDBS'
C         'FCZK' 'FCZK'
      IF (JIMPRT .GE. 1) THEN
         IF (KEY.EQ.'DIAM') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic Second Order Singlet (DIAM)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'DIAD') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic Second Order Singlet (DIAD)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'PSOK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (PSOK)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'PSOKI') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (PSOKI)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'LKIN') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (LKIN)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'FCZK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCZK)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'SDZK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDZK)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'FCBS') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCBS)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'SDBS') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDBS)',
     &   ' correction is under calculation...'
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
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                    ij
         operadores_A(1,1) = nuc_coord(1)//'NSCO'//CHRXYZ(1) !11
         operadores_A(2,1) = nuc_coord(1)//'NSCO'//CHRXYZ(2) !12
         operadores_A(3,1) = nuc_coord(1)//'NSCO'//CHRXYZ(3) !13
         operadores_A(4,1) = nuc_coord(2)//'NSCO'//CHRXYZ(1) !21
         operadores_A(5,1) = nuc_coord(2)//'NSCO'//CHRXYZ(2) !22
         operadores_A(6,1) = nuc_coord(2)//'NSCO'//CHRXYZ(3) !23
         operadores_A(7,1) = nuc_coord(3)//'NSCO'//CHRXYZ(1) !31
         operadores_A(8,1) = nuc_coord(3)//'NSCO'//CHRXYZ(2) !32
         operadores_A(9,1) = nuc_coord(3)//'NSCO'//CHRXYZ(3) !33
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         IF (KEY.EQ."DIAD") operadores_B(:,1) ='DARWIN  '
         IF (KEY.EQ."DIAM") operadores_B(:,1) ='MASSVELO'
      END IF


      IF(KEY.EQ."PSOK") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                    ij
         operadores_A(1,1) = 'PSOKE'//nuc_coord(1)           !11 
         operadores_A(2,1) = 'PSOKE'//nuc_coord(1)           !12 editar
         operadores_A(3,1) = 'PSOKE'//nuc_coord(1)           !13 editar
         operadores_A(4,1) = 'PSOKE'//nuc_coord(1)           !21 editar
         operadores_A(5,1) = 'PSOKE'//nuc_coord(2)           !22
         operadores_A(6,1) = 'PSOKE'//nuc_coord(1)           !23 editar
         operadores_A(7,1) = 'PSOKE'//nuc_coord(1)           !31 editar
         operadores_A(8,1) = 'PSOKE'//nuc_coord(1)           !32 editar
         operadores_A(9,1) = 'PSOKE'//nuc_coord(3)           !33
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         operadores_B(1,1) = 'XANGMOM '                      !11 
         operadores_B(2,1) = 'XANGMOM '                      !12 
         operadores_B(3,1) = 'XANGMOM '                      !13 
         operadores_B(4,1) = 'XANGMOM '                      !21 
         operadores_B(5,1) = 'YANGMOM '                      !22 
         operadores_B(6,1) = 'XANGMOM '                      !23 
         operadores_B(7,1) = 'XANGMOM '                      !31 
         operadores_B(8,1) = 'XANGMOM '                      !32 
         operadores_B(9,1) = 'ZANGMOM '                      !33 
      END IF

      IF(KEY.EQ."LKIN") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                    ij
         operadores_A(1,1) = 'PSO '//nuc_coord(1)//' '       !11 
         operadores_A(2,1) = 'PSO '//nuc_coord(1)//' '       !12
         operadores_A(3,1) = 'PSO '//nuc_coord(1)//' '       !13
         operadores_A(4,1) = 'PSO '//nuc_coord(2)//' '       !21
         operadores_A(5,1) = 'PSO '//nuc_coord(2)//' '       !22
         operadores_A(6,1) = 'PSO '//nuc_coord(2)//' '       !23
         operadores_A(7,1) = 'PSO '//nuc_coord(3)//' '       !31
         operadores_A(8,1) = 'PSO '//nuc_coord(3)//' '       !32
         operadores_A(9,1) = 'PSO '//nuc_coord(3)//' '       !33
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         operadores_B(1,1) = CHRXYZ(1)//'OZKE   '            !11 
         operadores_B(2,1) = CHRXYZ(2)//'OZKE   '            !12 
         operadores_B(3,1) = CHRXYZ(3)//'OZKE   '            !13 
         operadores_B(4,1) = CHRXYZ(1)//'OZKE   '            !21 
         operadores_B(5,1) = CHRXYZ(2)//'OZKE   '            !22 
         operadores_B(6,1) = CHRXYZ(3)//'OZKE   '            !23 
         operadores_B(7,1) = CHRXYZ(1)//'OZKE   '            !31 
         operadores_B(8,1) = CHRXYZ(2)//'OZKE   '            !32 
         operadores_B(9,1) = CHRXYZ(3)//'OZKE   '            !33 
      END IF

      IF((KEY.EQ."FCZK").OR.(KEY.EQ."FCBS")) THEN
         NopA = 1
         TRPLET_FLAG = .TRUE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                          ij
         operadores_A(:,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number  !all
C        -----------------------------------
         NopB = 1
         allocate(operadores_B(9,NopB))
         IF(KEY.EQ."FCZK") THEN
            operadores_B(:,1) = 'KINENERG'                            !all
         END IF

         IF(KEY.EQ."FCBS") THEN
            operadores_B(1,1) = 'SOMF  XX'                      !11 
            operadores_B(2,1) = 'SOMF  YX'                      !12 
            operadores_B(3,1) = 'SOMF  ZX'                      !13 
            operadores_B(4,1) = 'SOMF  YX'                      !21 
            operadores_B(5,1) = 'SOMF  YY'                      !22 
            operadores_B(6,1) = 'SOMF  ZY'                      !23 
            operadores_B(7,1) = 'SOMF  ZX'                      !31 
            operadores_B(8,1) = 'SOMF  ZY'                      !32 
            operadores_B(9,1) = 'SOMF  ZZ'                      !33 
         END IF
      END IF


c     Under construction
      IF((KEY.EQ."SDZK").OR.(KEY.EQ."SDBS")) THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .TRUE.
         allocate(operadores_A(9,NopA))
         allocate(operadores_B(9,NopB))
C        -----------------------------------
c                  3*i+j-3                               ij
         operadores_A(1,1) = 'SD '//nuc_coord(1)//' x'  !11
         operadores_A(2,1) = 'SD '//nuc_coord(1)//' y'  !12
         operadores_A(3,1) = 'SD '//nuc_coord(1)//' z'  !13
         operadores_A(4,1) = 'SD '//nuc_coord(2)//' x'  !21
         operadores_A(5,1) = 'SD '//nuc_coord(2)//' y'  !22
         operadores_A(6,1) = 'SD '//nuc_coord(2)//' z'  !23
         operadores_A(7,1) = 'SD '//nuc_coord(3)//' x'  !31
         operadores_A(8,1) = 'SD '//nuc_coord(3)//' y'  !32
         operadores_A(9,1) = 'SD '//nuc_coord(3)//' z'  !33
C        -----------------------------------
         IF(KEY.EQ."SDZK") THEN
c                     3*i+j-3               ij
            operadores_B(1,1) = 'dd/dxdx ' !11
            operadores_B(2,1) = 'dd/dxdy ' !12
            operadores_B(3,1) = 'dd/dxdz ' !13
            operadores_B(4,1) = 'dd/dxdy ' !21
            operadores_B(5,1) = 'dd/dydy ' !22
            operadores_B(6,1) = 'dd/dydz ' !23
            operadores_B(7,1) = 'dd/dxdz ' !31
            operadores_B(8,1) = 'dd/dydz ' !32
            operadores_B(9,1) = 'dd/dzdz ' !33
         END IF
C        -----------------------------------
         IF(KEY.EQ."SDBS") THEN
c                     3*i+j-3               ij
            operadores_B(1,1) = 'SOMF  XX' !11
            operadores_B(2,1) = 'SOMF  YX' !12
            operadores_B(3,1) = 'SOMF  ZX' !13
            operadores_B(4,1) = 'SOMF  YX' !21
            operadores_B(5,1) = 'SOMF  YY' !22
            operadores_B(6,1) = 'SOMF  ZY' !23
            operadores_B(7,1) = 'SOMF  ZX' !31
            operadores_B(8,1) = 'SOMF  ZY' !32
            operadores_B(9,1) = 'SOMF  ZZ' !33
         END IF
      END IF

c     Elemento ij del tensor de shielding
      do i = 1, 3 
      do j = 1, 3
c         LR = 0.d0
c         TOTALLR = 0.d0
         do k = 1, NopA
         do l = 1, NopB
            operadores_duo(1) = operadores_A(3*i+j-3,k)
            operadores_duo(2) = operadores_B(3*i+j-3,l)

            call get_LRvalue(LR(i,j),operadores_duo,TRIPLETE) !LR from here has an extra minus sign
c            WRITE(LUPRI,*) "AUCAR-LR",operadores_duo,i,j,
c     &                     calfa*calfa*1.0D6*LR(i,j)
c            WRITE(LUPRI,*) "AUCAR-LR",calfa*calfa*1.0D6*LR
            TOTALLR(i,j) = TOTALLR(i,j) + LR(i,j)
            IF (JIMPRT .GT. 2) write(lupri,*) operadores_duo(1),
     &         ' ',operadores_duo(2)," LR: ",LR(i,j)
         end do
         end do

         LRtensor(i,j) = TOTALLR(i,j)
         IF (JIMPRT .GT. 2) THEN
            write(lupri,*) KEY,"(",i,",",j,"):",LRtensor(i,j)
            write(lupri,*)
         END IF

      end do
      end do

      IF (KEY.EQ.'DIAM') LRDIAM = LRtensor
      IF (KEY.EQ.'DIAD') LRDIAD = LRtensor
      IF (KEY.EQ.'PSOK') LRPSOK = LRtensor
      IF (KEY.EQ.'LKIN') LRLKIN = LRtensor
      IF (KEY.EQ.'FCZK') LRFCZK = LRtensor
      IF (KEY.EQ.'FCBS') LRFCBS = LRtensor
      IF (KEY.EQ.'SDZK') LRSDZK = LRtensor
      IF (KEY.EQ.'SDBS') LRSDBS = LRtensor

      deallocate(operadores_A)
      deallocate(operadores_B)

      CALL QEXIT('LRSLINEAR')

      RETURN

      END 


C  /* Deck newlrinput_all */
      SUBROUTINE LRINPUT_ALL(LRATOMS,NLRATOMS,TRIPLETE)

#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
#include "nuclei.h"
#include "chrnos.h"
#include "chrxyz.h"

      CHARACTER*8, allocatable :: operadores(:)
      INTEGER                  :: NopA, NopB, NLRoperators
      INTEGER :: SUMA_B, NLRATOMS
      INTEGER, intent(in)  :: LRATOMS(NLRATOMS)
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      character*10 :: namefile
      INTEGER :: fileID
      LOGICAL :: TRIPLETE

      CALL QENTER('NEWLRINPUT_SINGLETS_ALL')
c      CALL GETTIM(TIMHER,WALHER)

c     NopA: number of A operators for each shield tensor element
c     NopB: number of B operators for each shield tensor element

C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'DIAM'  DIAMAGNETIC <<DIA,MV>>
C         'DIAD'  DIAMAGNETIC <<DIA,DW>>
C         'PSOK'  PARAMAGNETIC with built-in integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'LKIN'  PARAMAGNETIC
C         'FCBS' 'SDBS'
C         'FCZK' 'FCZK'

      namefile='DALTON.INP'
      fileID = -1
c      CALL GPOPEN(fileID,namefile,'NEW',' ','FORMATTED',
c     &            IDUMMY,.FALSE.)
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
         IF (SIGMAD1S) NLRoperators = NLRoperators + 11
         IF (SIGMAP1S) NLRoperators = NLRoperators + 12
      END IF
      IF ((TRIPLETE).AND.(SIGMAP1T)) NLRoperators = NLRoperators + 23

      allocate(operadores(NLRoperators))
      el = 1

      IF (.not.(TRIPLETE)) THEN
      IF (SIGMAD1S) THEN !SIGMAD1S: DIAM and DIAD
C        ---------------------------------------------------
         operadores(el)   = nuc_coord(1)//'NSCO'//CHRXYZ(1)
         operadores(el+1) = nuc_coord(1)//'NSCO'//CHRXYZ(2)
         operadores(el+2) = nuc_coord(1)//'NSCO'//CHRXYZ(3)
         operadores(el+3) = nuc_coord(2)//'NSCO'//CHRXYZ(1)
         operadores(el+4) = nuc_coord(2)//'NSCO'//CHRXYZ(2)
         operadores(el+5) = nuc_coord(2)//'NSCO'//CHRXYZ(3)
         operadores(el+6) = nuc_coord(3)//'NSCO'//CHRXYZ(1)
         operadores(el+7) = nuc_coord(3)//'NSCO'//CHRXYZ(2)
         operadores(el+8) = nuc_coord(3)//'NSCO'//CHRXYZ(3)
C        ---------------------------------------------------
C        ---------------------------------------------------
C        ---------------------------------------------------
         operadores(el+9)  = 'DARWIN  '
C        ---------------------------------------------------
         operadores(el+10) = 'MASSVELO'
C        ---------------------------------------------------
         el = el + 11
      END IF

      IF (SIGMAP1S) THEN !SIGMAP1S: PSOK and LKIN
C        ---------------------------------------------------
         operadores(el)   = 'PSOKE'//nuc_coord(1)
         operadores(el+1) = 'PSOKE'//nuc_coord(2)
         operadores(el+2) = 'PSOKE'//nuc_coord(3)
C        ---------------------------------------------------
         operadores(el+3) = 'XANGMOM '
         operadores(el+4) = 'YANGMOM '
         operadores(el+5) = 'ZANGMOM '
C        ---------------------------------------------------
C        ---------------------------------------------------
         operadores(el+6) = 'PSO '//nuc_coord(1)//' '
         operadores(el+7) = 'PSO '//nuc_coord(2)//' '
         operadores(el+8) = 'PSO '//nuc_coord(3)//' '
C        ---------------------------------------------------
         operadores(el+9)  = CHRXYZ(1)//'OZKE   '
         operadores(el+10) = CHRXYZ(2)//'OZKE   '
         operadores(el+11) = CHRXYZ(3)//'OZKE   '
C        ---------------------------------------------------
         el = el + 12
      END IF
      END IF

      IF (TRIPLETE) THEN
      IF (SIGMAP1T) THEN !FCZK, FCBS, SDZK and SDBS
C        ---------------------------------------------------
         operadores(el)   = 'FC '//NAMN(LRATOM)(1:2)//nuc_number !FCZK & FCBS
C        ---------------------------------------------------
         operadores(el+1) = 'SD '//nuc_coord(1)//' x'            !SDZK & SDBS
         operadores(el+2) = 'SD '//nuc_coord(1)//' y'
         operadores(el+3) = 'SD '//nuc_coord(1)//' z'
         operadores(el+4) = 'SD '//nuc_coord(2)//' x'
         operadores(el+5) = 'SD '//nuc_coord(2)//' y'
         operadores(el+6) = 'SD '//nuc_coord(2)//' z'
         operadores(el+7) = 'SD '//nuc_coord(3)//' x'
         operadores(el+8) = 'SD '//nuc_coord(3)//' y'
         operadores(el+9) = 'SD '//nuc_coord(3)//' z'
C        ---------------------------------------------------
C        ---------------------------------------------------
         operadores(el+10) = 'KINENERG'                          !FCZK
C        ---------------------------------------------------
         operadores(el+11) = 'dd/dxdx '                          !SDZK
         operadores(el+12) = 'dd/dxdy '
         operadores(el+13) = 'dd/dxdz '
         operadores(el+14) = 'dd/dydy '
         operadores(el+15) = 'dd/dydz '
         operadores(el+16) = 'dd/dzdz '
C        ---------------------------------------------------
         operadores(el+17) = 'SOMF  XX'                          !FCBS & SDBS
         operadores(el+18) = 'SOMF  YX'
         operadores(el+19) = 'SOMF  ZX'
         operadores(el+20) = 'SOMF  YY'
         operadores(el+21) = 'SOMF  ZY'
         operadores(el+22) = 'SOMF  ZZ'
C        ---------------------------------------------------
         el = el + 23
      END IF
      END IF

      do i = 1, SIZE(operadores,DIM=1)
         write(fileID,'(a)') '.PROPRT'
         write(fileID,'(a)') operadores(i)
      end do


      deallocate(operadores)
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
      logical :: file_exists, file_open, TRIPLETE
      INTEGER :: fileID
      CHARACTER*164 STRING, TMPSTRING
      CHARACTER*8, intent(in) :: labels_in(2)
      CHARACTER*8:: labels_read(2)
      CHARACTER*10 :: LRfile
      REAL*8, intent(out) :: LRVAL

c

#include "priunit.h"
#include "lrescinf.h"
#include "chrnos.h"      

C Based on get_QRvalue
      IF (TRIPLETE) THEN
         LRfile = 'LRESC-LR-T'
      ELSE
         LRfile = 'LRESC-LR-S'
      END IF

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

      RETURN
      END
