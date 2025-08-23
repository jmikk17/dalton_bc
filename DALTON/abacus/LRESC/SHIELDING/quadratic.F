
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
C  /* Deck QUADRATICLR */
      SUBROUTINE QUADRATICLR(KEY,WORK,LWORK)

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
      CHARACTER*8 :: operadores_trio(3)
      CHARACTER*1 :: spin_trio(3)

c      double precision, dimension(:), allocatable :: xmos
c      double precision, dimension(:,:), allocatable :: cmos
c      double precision, dimension(:,:), allocatable :: Identity
c      double precision, dimension(:,:), allocatable :: OVERLP
c      integer :: naos,nmos

      CHARACTER*8 RTNLBL(2)

      LOGICAL FOUND, matched
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number

      REAL*8 QR(3,3),TOTALQR(3,3)
      REAL*8 QRtensor(3,3)
      Character*8, dimension(:,:,:), allocatable :: operadores_A
      Character*8, dimension(:,:,:), allocatable :: operadores_B
      Character*8, dimension(:,:,:), allocatable :: operadores_C

      INTEGER                  :: NopA, NopB, NopC
      LOGICAL :: TRIPLETE,TRPLET_FLAG, file_exists


      
      CALL QENTER('LRSCQUAD')
      CALL GETTIM(TIMQR,WALQR)

      matched = .FALSE.

c      naos = NBAST ! NBAST: Total number of (atomic) basis functions
c      nmos= NORBT  ! NORBT: Total number of molecular orbitals

      TOTALQR = 0.d0
      QR = 0.d0

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

      IF ((KEY.EQ.'OZ-PSO-MV').or.(KEY.EQ.'OZ-PSO-Dw')) THEN
         NopA = 1
         NopB = 1
         NopC = 1
         allocate(operadores_A(3,3,NopA))
         allocate(operadores_B(3,3,NopB))
         allocate(operadores_C(3,3,NopC))
         spin_trio(:)='0'
         operadores_A(:,1,1) = 'XANGMOM '
         operadores_A(:,2,1) = 'YANGMOM '
         operadores_A(:,3,1) = 'ZANGMOM '
         operadores_B(1,:,1) = 'PSO '//nuc_coord(1)//' '
         operadores_B(2,:,1) = 'PSO '//nuc_coord(2)//' '
         operadores_B(3,:,1) = 'PSO '//nuc_coord(3)//' '
         IF (KEY.EQ.'OZ-PSO-MV') operadores_C(:,:,:) = 'MASSVELO'
         IF (KEY.EQ.'OZ-PSO-Dw') operadores_C(:,:,:) = 'DARWIN  '
      END IF

      IF (KEY.EQ.'OZ-FC-SO ') THEN
         NopA = 1
         NopB = 1
         NopC = 1
         allocate(operadores_A(3,3,NopA))
         allocate(operadores_B(3,3,NopB))
         allocate(operadores_C(3,3,NopC))
         spin_trio(1)='0'
         spin_trio(2)='1'
         spin_trio(3)='1'
         operadores_A(:,1,1) = 'XANGMOM '
         operadores_A(:,2,1) = 'YANGMOM '
         operadores_A(:,3,1) = 'ZANGMOM '
         operadores_B(:,:,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         operadores_C(1,:,1) = 'X1SPNORB'
         operadores_C(2,:,1) = 'Y1SPNORB'
         operadores_C(3,:,1) = 'Z1SPNORB'
      END IF

      IF (KEY.EQ.'OZ-SD-SO ') THEN
         NopA = 3
         NopB = 3
         NopC = 3
         allocate(operadores_A(3,3,NopA))
         allocate(operadores_B(3,3,NopB))
         allocate(operadores_C(3,3,NopC))
         spin_trio(1)='0'
         spin_trio(2)='1'
         spin_trio(3)='1'
         operadores_A(:,1,:) = 'XANGMOM '
         operadores_A(:,2,:) = 'YANGMOM '
         operadores_A(:,3,:) = 'ZANGMOM '

         DO K=1, 3
            operadores_B(1,1,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !11
            operadores_B(1,2,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !12
            operadores_B(1,3,K) = 'SD '//nuc_coord(1)//' '//CHRXYZ(-K)  !13
            operadores_B(2,1,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !21
            operadores_B(2,2,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !22
            operadores_B(2,3,K) = 'SD '//nuc_coord(2)//' '//CHRXYZ(-K)  !23
            operadores_B(3,1,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !31
            operadores_B(3,2,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !32
            operadores_B(3,3,K) = 'SD '//nuc_coord(3)//' '//CHRXYZ(-K)  !33
         END DO
         operadores_C(:,:,1) = 'X1SPNORB'
         operadores_C(:,:,2) = 'Y1SPNORB'
         operadores_C(:,:,3) = 'Z1SPNORB'
      END IF


c      IF (KEY.EQ.'SZ-PSO-MV') THEN
c        constante = 1000000/(137.036**2)
c         allocate(operadores(3,3))
c         allocate(spin_rank(3,3))
c         operadores(1,1) = 'X1SPNORB'
c         operadores(2,1) = 'Y1SPNORB'
c         operadores(3,1) = 'Z1SPNORB'
c         spin_rank(:,1) = '1'
c         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
c         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
c         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
c         spin_rank(:,2) = '0'
c         operadores(:,3) = 'MASSVELO'
c         spin_rank(:,3) = '0'
c      END IF

c      IF (KEY.EQ.'SZ-PSO-Dw') THEN
c        constante = 1000000/(137.036**2)
c         allocate(operadores(3,3))
c         allocate(spin_rank(3,3))
c         operadores(1,1) = 'X1SPNORB'
c         operadores(2,1) = 'Y1SPNORB'
c         operadores(3,1) = 'Z1SPNORB'
c         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
c         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
c         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
c         operadores(:,3) = 'DARWIN  '
c         spin_rank(:,1) = '1'
c         spin_rank(:,2) = '0'
c         spin_rank(:,3) = '0'
c      END IF




c     Está entrando dos veces a la condición, no?
c      IF (KEY.EQ.'SZ-SD-SO '.OR. KEY.EQ.'SZ-FC-SO ') THEN
c        CALL RD_SIRIFC('CMO',FOUND,xmos)
c        cmos=0.0
c        k = 1
c        do i=1,nmos                  !i corresponds to a MO
c          do j=1,naos                              !j corresponds to j coef. of the ith MO
c            cmos(i,j) = xmos(k)                ! Transforms vector into matrix
c           write(lupri,*) "JUAN-cmos: ",i,j,cmos(i,j)
c           write(lupri,*)
c            k = k + 1
c          end do
c        end do
c
c        CALL GETDAT(RTNLBL(1),RTNLBL(2))
c        RTNLBL(1)='1 A  21 '
C       Replace time information with symmetry information
c        RTNLBL(2)='SQUARE  '
c        RTNLBL(2)='ANTISYMM'
c        call WRITEAOPJUAN(matmul(transpose(cmos),cmos),NBAST,
c     &        'IDENTIDA',RTNLBL) !Write Matrix in AOPROPER
c        write(LUPRI,*) "MATRIZ IDENTIDAD:"
c        call OUTPUT(matmul(transpose(cmos),cmos),1,NORBT,1,NORBT,
c     &         NORBT,NORBT,1,LUPRI) ! print matrix

c        END IF
c
c      IF (KEY.EQ.'SZ-FC-SO ') THEN
c        Constante: (137,036^-2/2)*1000000
c         constante = 1000000/(2*(137.036**2))
c         allocate(operadores(3,3))
c         allocate(spin_rank(3,3))
c         operadores(1,2) = 'X1SPNORB'
c         operadores(2,2) = 'Y1SPNORB'
c         operadores(3,2) = 'Z1SPNORB'
c         operadores(:,2) = 'IDENTIDA'
c         spin_rank(:,2) = '1'
c         operadores(1,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
c         operadores(2,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
c         operadores(3,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
c         spin_rank(:,1) = '1'
c         operadores(1,3) = 'X1SPNORB'
c         operadores(2,3) = 'Y1SPNORB'
c         operadores(3,3) = 'Z1SPNORB'
c         spin_rank(:,3) = '1'
c      END IF

c      IF (KEY.EQ.'SZ-SD-SO ') THEN
cc         CALL GETDAT(RTNLBL(1),RTNLBL(2))
cc         RTNLBL(1)='1 A  21 '
Cc        Replace time information with symmetry information
cc         RTNLBL(2)='SYMMETRI'
cc         Identity = 0                          ! Initialize the array.
cc         forall(j = 1:NBAST) Identity(j,j) = 1 ! Set the diagonal.
cc         call WRITEAOPJUAN(Identity,NBAST,'IDENTITY',RTNLBL) !Write Matrix in AOPROPER

cc         WRITE(lupri,*)' IDENTITY : '
cc         CALL OUTPUT(Identity,1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
cc         do i=1,NBAST
cc         do j=1,NBAST
cc            write(lupri,*) i,j,OVERLP(i,j)
cc         enddo
cc         enddo
cc         WRITE(lupri,*) '   '

cc         WRITE(lupri,*)' OVERLAP : '
cc         call READAOP(OVERLP, NBAST, 'OVERLAP ') !gets OVERLAP matrix
cc         CALL OUTPUT(OVERLP,1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
cc         WRITE(lupri,*)
cc         do i=1,NBAST
cc         do j=1,NBAST
cc            write(lupri,*) i,OVERLP(i,i)
cc         enddo
cc         enddo


cc         CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO))
cc         CALL RD_SIRIFC('CMO',FOUND,xmos)
cc                cmos=0.0
cc       k = 1
cc       do i=1,nmos                  !i corresponds to a MO
cc          do j=1,naos                              !j corresponds to j coef. of the ith MO
cc             cmos(i,j) = xmos(k)                ! Transforms vector into matrix
cc             write(lupri,*) "JUAN-cmos: ",i,j,cmos(i,j)
cc             write(lupri,*)
cc            k = k + 1
cc         end do
cc       end do

cc         do i=1,NBAST
cc         do j=1,NBAST
cc            write(lupri,*) "JUAN-Ident: ",
cc     &                                i,j,matmul(transpose(cmos),cmos)!(i,j)
cc             write(lupri,*)
cc         enddo
cc         enddo


cc         CALL GETDAT(RTNLBL(1),RTNLBL(2))
cc         RTNLBL(1)='1 A  21 '
Cc        Replace time information with symmetry information
cc         RTNLBL(2)='SQUARE  '
cc        RTNLBL(2)='ANTISYMM'
cc      call WRITEAOPJUAN(matmul(transpose(cmos),cmos),NBAST,
cc    &        'IDENTIDA',RTNLBL) !Write Matrix in AOPROPER
cc        Constante: (137,036^-2/2)*1000000
cc         constante = 1000000/(2*(137.036**2))
c         allocate(operadores(9,3))
c         allocate(spin_rank(9,3))
c         operadores(:,1) = 'OVERLAP '
c         operadores(:,1) = 'IDENTIDA'
c         spin_rank(:,1) = '1'
c                              nuc coord
c         operadores(1,2) = 'SD '//nuc_coord(1)//' x'
c         operadores(2,2) = 'SD '//nuc_coord(1)//' y'
c         operadores(3,2) = 'SD '//nuc_coord(1)//' z'
c         operadores(4,2) = 'SD '//nuc_coord(2)//' x'
c         operadores(5,2) = 'SD '//nuc_coord(2)//' y'
c         operadores(6,2) = 'SD '//nuc_coord(2)//' z'
c         operadores(7,2) = 'SD '//nuc_coord(3)//' x'
c         operadores(8,2) = 'SD '//nuc_coord(3)//' y'
c         operadores(9,2) = 'SD '//nuc_coord(3)//' z'
c         spin_rank(:,2) = '1'

c         operadores(1,3) = 'X1SPNORB'
c         operadores(2,3) = 'Y1SPNORB'
c         operadores(3,3) = 'Z1SPNORB'
c         spin_rank(:,3) = '1'

c        END IF

!     Logical for first call to read LRESC-QR file
      inquire(file = 'LRESC-QR', exist = file_exists)

c     ij shielding tensor element
      do i = 1, 3 !
      do j = 1, 3 !
         IF (NopA+NopB+NopC.LT.6) THEN !Combine all A,B,C operators (none of them has 3 components)
            do k = 1, NopA
            do l = 1, NopB
            do m = 1, NopC
               operadores_trio(1)=operadores_A(i,j,k)
               operadores_trio(2)=operadores_B(i,j,l)
               operadores_trio(3)=operadores_C(i,j,m)
               if (file_exists)
     &            call get_QRvalue(QR(i,j),operadores_trio,matched)
               if (.not.matched) then
                 call create_QRinput(operadores_trio,spin_trio)
                 LUPRI = 6232
                 OPEN(UNIT=LUPRI, FILE='LRESC-QR',status='unknown',
     &                position='append')
                 call RSPDRV(WORK,LWORK)
                 CLOSE(LUPRI)
                 LUPRI  = 2
                 CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
                 CALL GPCLOSE(LURSP3,'DELETE')
                 CALL GPCLOSE(LURSP4,'DELETE')
                 CALL GPCLOSE(LURSP5,'DELETE')
                 call get_QRvalue(QR(i,j),operadores_trio,matched)
               END IF
               matched = .FALSE.
               TOTALQR(i,j) = TOTALQR(i,j) + QR(i,j)
            end do
            end do
            end do

         ELSE
            do k = 1, 3
               operadores_trio(1)=operadores_A(i,j,k)
               operadores_trio(2)=operadores_B(i,j,k)
               operadores_trio(3)=operadores_C(i,j,k)
               if (file_exists)
     &            call get_QRvalue(QR(i,j),operadores_trio,matched)
               if (.not.matched) then
                 call create_QRinput(operadores_trio,spin_trio)
                 LUPRI = 6232
                 OPEN(UNIT=LUPRI, FILE='LRESC-QR',status='unknown',
     &                position='append')
                 call RSPDRV(WORK,LWORK)
                 CLOSE(LUPRI)
                 LUPRI  = 2
                 CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
                 CALL GPCLOSE(LURSP3,'DELETE')
                 CALL GPCLOSE(LURSP4,'DELETE')
                 CALL GPCLOSE(LURSP5,'DELETE')
                 call get_QRvalue(QR(i,j),operadores_trio,matched)
               END IF
               matched = .FALSE.
               TOTALQR(i,j) = TOTALQR(i,j) + QR(i,j)
            end do
         END IF
         QRtensor(i,j) = TOTALQR(i,j)

      end do
      end do

      IF (KEY.EQ.'OZ-PSO-MV') QRPSOOZMV = QRtensor
      IF (KEY.EQ.'OZ-PSO-Dw') QRPSOOZDW = QRtensor
      IF (KEY.EQ.'OZ-FC-SO ') QROZFCSO = QRtensor
      IF (KEY.EQ.'OZ-SD-SO ') QROZSDSO = QRtensor

      deallocate(operadores_A)
      deallocate(operadores_B)
      deallocate(operadores_C)


      CALL GETTIM(TEND,WEND)
      TIMQR = TEND - TIMQR
      WALQR = WEND - WALQR
      IF (JIMPRT .GE. 2) THEN
         CALL TIMTXT('    Total CPU  time used in LRESC-QR('
     &                     //KEY//'):',TIMQR,LUPRI)
         CALL TIMTXT('    Total wall time used in LRESC-QR('
     &                     //KEY//'):',WALQR,LUPRI)
      END IF
      CALL QEXIT('LRSCQUAD')

      RETURN

      END 


      subroutine create_QRinput(operadores,spin_rank)
      character*10 :: namefile
      logical            :: file_exists, file_open
      character(100) :: line
      INTEGER :: fileID
      CHARACTER*8, intent(in) :: operadores(3)
      CHARACTER*1, intent(in) :: spin_rank(3)

#include "priunit.h"
      namefile='DALTON.INP'
      fileID = -1
      CALL GPOPEN(fileID,namefile,'NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      write(fileID,'(a)') '**DALTON INPUT'
      write(fileID,'(a)') '**RESPONSE'
      write(fileID,'(a)') '*LINEAR'
      write(fileID,'(a)') '.PRINT'
      write(fileID,'(a)') ' 0'
      write(fileID,'(a)') '*QUADRA'
      write(fileID,'(a)') '.PRINT'
      write(fileID,'(a)') ' 0'
      write(fileID,'(a)') '.ISPABC'
      write(fileID,'(a)') spin_rank(1)//'    '//
     &                    spin_rank(2)//'    '//spin_rank(3)
      write(fileID,'(a)') '.APROP'
      write(fileID,'(a)') operadores(1)//' '//spin_rank(1)
      write(fileID,'(a)') '.BPROP'
      write(fileID,'(a)') operadores(2)//' '//spin_rank(2)
      write(fileID,'(a)') '.CPROP'
      write(fileID,'(a)') operadores(3)//' '//spin_rank(3)
      write(fileID,'(a)') '**END OF DALTON INPUT'

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

      RETURN
      END

      subroutine get_QRvalue(QUAVAL,labels_in,match)
      logical :: file_exists, file_open, match
      INTEGER :: fileID
      CHARACTER STRING*75
      CHARACTER*8, intent(in) :: labels_in(3)
      CHARACTER*8:: labels_read(3)
      REAL*8 FRQB, FRQC
      REAL*8, intent(out) :: QUAVAL
c

#include "priunit.h"

      match = .FALSE.

C Based on some lines of SUBROUTINE BCCHK

      inquire(file = 'LRESC-QR', exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            fileID = -1
            CALL GPOPEN(fileID,'LRESC-QR','OLD',' ','FORMATTED',
     &           IDUMMY,.FALSE.)
         end if
         rewind(fileID)
      else
        call quit('LRESC-QR file not found')
      end if


 346  READ (fileID,'(A75)',END=347) STRING
      IF (STRING(1:8) .EQ. '@ Quadra') THEN
          READ (fileID,'(30X,A8)') labels_read(1)
          READ (fileID,'(30X,A8)') labels_read(2)
          READ (fileID,'(30X,A8)') labels_read(3)

      IF (labels_read(1) .EQ. labels_in(1) .AND.
     &    labels_read(2) .EQ. labels_in(2) .AND.
     &    labels_read(3) .EQ. labels_in(3) ) then
         READ (fileID,'(/,30X,3F20.8)') FRQB, FRQC, QUAVAL
         match = .true.
      ELSE
         GOTO 346
      END IF

      ELSE
         GOTO 346
      END IF

 347         CONTINUE

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

      RETURN
      END
