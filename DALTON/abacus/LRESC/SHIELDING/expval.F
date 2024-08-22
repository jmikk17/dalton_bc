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
      SUBROUTINE AVELRSC(KEY,WORK,LWORK)
C
C     Purpose:
C     Get average values from AOPROPER operators in order
C     to get LRESC
C     Author:
C           J. I. Melo
C     Edited:
C           Juan J. Aucar         April 2020              
C...  Note:
C...  This subroutine was written by Juan Ignacio Melo using
C...  the subroutine ABACTOCD  as a model (2012)
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "trkoor.h"
#include "maxorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "cbilnr.h"
#include "infpri.h"
      integer indice
      LOGICAL FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1,LISTA1(100)
      CHARACTER*6 LABEL2
      CHARACTER*4 KEY
      CHARACTER*3 char
      PARAMETER (D05=0.5D0,D025=0.25)
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      INTEGER naoss, nmos
      real(8), allocatable :: xmos(:)
      real(8), allocatable :: cmos(:,:)
      real(8), allocatable :: AOPERATOR(:,:)
      real(8), allocatable :: BOPERATOR(:,:)
      real(8), allocatable :: aux0(:,:)
      real(8), allocatable :: aux1(:,:)
      real(8), allocatable :: aux2(:,:)
      real(8), allocatable :: aux3(:,:)
      CHARACTER*8 RTNLBL(2)

C
#include "cbiexc.h"
#include "inflin.h"
#include "infvar.h"
#include "infdim.h"
#include "inforb.h"
#include "nuclei.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "maxmom.h"
#include "maxaqn.h"
#include "symmet.h"
#include "abainf.h"
#include "gnrinf.h"

C
#include "lrescinf.h"
#include "chrxyz.h"
#include "chrnos.h"
#include "orgcom.h"
C
      IPRLNR = JIMPRT
      IPRRSP = -1
       naoss = NBAST ! NBAST: Total number of (atomic) basis functions
       nmos= NORBT !NORBT: Total number of molecular orbitals
cxu
C
C
C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'FCAV'  <Fc>
C         'DIAK'  <Dia.Kin>
C         'ANGP'  <L.Pso>
C         'SDAV'  <SD>
C
      CALL QENTER('AVELRSC')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (JIMPRT .GE. 1) THEN
!         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
         IF (KEY.EQ.'DINR') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic non-relativistic (DINR)',
     &   ' contribution is under calculation'
         IF (KEY.EQ.'FCAV') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic First Order Singlet (FC)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'DIAK') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic First Order Singlet (DIAK)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'ANGP') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic First Order Singlet (ANGP)',
     &   ' correction is under calculation'
         IF (KEY.EQ.'SDAV') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic First Order Singlet (SDAV)',
     &   ' correction is under calculation'
!         IF (KEY.EQ.'FCt1') WRITE(LUPRI,'(5X,A,A)')
      END IF
C

C
C
C     Get reference state
C     ===================
C
C     1. Work Allocations:
C
      LUDV   = N2ASHX
      LPVX   = 0
      KFREE  = 1
      LFREE  = LWORK
c      IF (JIMPRT.GE.4) Then
c         WRITE(lupri,'(A,3F12.8)') ' orgcom.h : GAGORG :', GAGORG
c         WRITE(lupri,'(A,3F12.8)') '            ORIGIN :', ORIGIN
c         WRITE(lupri,'(A,3F12.8)') '            CMXYZ  :', CMXYZ
c         WRITE(lupri,*)
c         WRITE(lupri,*) ' alocando 1 : ANTES MEMGET :'
c         WRITE(lupri,*) ' ----------------------------'
c         WRITE(lupri,*) '     KFREE = 1'
c         WRITE(lupri,*) '     LFREE = LWORK :         ', LWORK
c         WRITE(lupri,*) '                             '
C
c         WRITE(lupri,*) ' COMMON VARIABLES on LINEAR '
c         WRITE(lupri,*) ' ----------------------------'
c         WRITE(lupri,*) '     N2ASHX : ni idea ' , N2ASHX
c         WRITE(lupri,*) '     NASHT # Active Orbitals = 0 ? :', NASHT
c         WRITE(lupri,*) '     LISTA1y2 (4*MXCOOR+9) = ', 4*MXCOOR+9
c         WRITE(lupri,*) '     MXCOOR              = ', MXCOOR
c         WRITE(lupri,*) '     NCMOT = NORB * NORB = ', NCMOT
c         WRITE(lupri,*) '   '
c         WRITE(lupri,*) ' memget....  '
c      ENDIF
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK ,KFREE ,LFREE)
C
      KWORK1 = KFREE
      WORK1  = LFREE
c      IF (JIMPRT.GE.4) Then
c         WRITE(lupri,*) '   '
c         WRITE(lupri,*) ' AFTER MEMGET  '
c         WRITE(lupri,*) ' ----------------------------'
c         WRITE(lupri,*) '        KCMO, NCMOT     =  ', KCMO,NCMOT
c         WRITE(lupri,*) '        KUDV, LUDV      =  ', KUDV,LUDV
c         WRITE(lupri,*) '        KPVX, LPVX      =  ', KPVX,LPVX
c         WRITE(lupri,*) '        KPXINDX, LCINDX =  ', KXINDX,KXINDX
c         WRITE(lupri,*) '        KWORK1 = KFREE :   ', KFREE
c         WRITE(lupri,*) '        WORK1  = LFREE :   ', LFREE
c         WRITE(lupri,*) '   '
c      ENDIF

      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO))
 
      IF (.NOT.FOUND) 
     & CALL QUIT('***Error AVELRSC: CMO is not in SIRIFC***')
c      IF (JIMPRT.GT.4) THEN
c         WRITE(lupri,*)' CMOS : '
c         CALL OUTPUT(WORK(KCMO),1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
c         WRITE(lupri,*) '   '
c      ENDIF

      ISYM = 1
C
C     we keep this just in case
c      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
C
C
C     Construct property-integrals and WRITE to LUPROP
C     ================================================
C
C     2. Work Allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
c      IF (JIMPRT.GE.4) Then
c         WRITE(lupri,*) '        '
c         WRITE(lupri,*) ' stills to allocate : '
c         WRITE(lupri,*) ' ----------------------------'
c         WRITE(lupri,*) ' KIDSYM = KWORK1           : ' , KIDSYM
c         WRITE(lupri,*) ' KIDADR = KIDSYM + 9MXCENT : ' , KIDADR
c         WRITE(lupri,*) ' KWORK2 = KIDADR + 9MXCENT : ' , KWORK2
c         WRITE(lupri,*) ' LWORK2 = LWORK - KWORK2   : ' , LWORK2
c      ENDIF
C
C
C  
C ==============================================================================
C
C  Starting Labels stuff
C
C ==============================================================================
      NLAB = 0
      LISTA1='        '
      npos1 = 3*LRATOM-2
!      write(lupri,*)' selected atom is :', LRATOM
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

C===============================================================================
C  SINGLET CALCULATIONS : DNR
C ==============================================================================
C ---------------------------------------------------------------------
C
C      Look for LABELS : NSCO
C     ---------------------------
      IF(KEY.EQ."DINR") THEN
        NLAB=9
        LISTA1(1) = nuc_coord(1)//'NSCO'//CHRXYZ(1) !11
        LISTA1(2) = nuc_coord(1)//'NSCO'//CHRXYZ(2) !12
        LISTA1(3) = nuc_coord(1)//'NSCO'//CHRXYZ(3) !13
        LISTA1(4) = nuc_coord(2)//'NSCO'//CHRXYZ(1) !21
        LISTA1(5) = nuc_coord(2)//'NSCO'//CHRXYZ(2) !22
        LISTA1(6) = nuc_coord(2)//'NSCO'//CHRXYZ(3) !23
        LISTA1(7) = nuc_coord(3)//'NSCO'//CHRXYZ(1) !31
        LISTA1(8) = nuc_coord(3)//'NSCO'//CHRXYZ(2) !32
        LISTA1(9) = nuc_coord(3)//'NSCO'//CHRXYZ(3) !33
      ENDIF
C ---------------------------------------------------------------------
C===============================================================================
C  SINGLET CALCULATIONS : D1S
C ==============================================================================
C
C      Look for LABELS : FERMI
C     ---------------------------
      IF(KEY.EQ."FCAV") THEN
         NLAB=1
         WRITE(LISTA1(1),'(A3,A2,I3.3)') 'FC ',
     &        NAMN(LRATOM)(1:2), LRATOM
      ENDIF
C
C ---------------------------------------------------------------------
C
C      Look for LABELS : DIAKIN
C     ---------------------------
      IF(KEY.EQ."DIAK") THEN
        NLAB=9
        indice = 1
        DO I=npos1,npos1+2
        DO IJ = 1, 3
           LISTA1(indice)= CHRNOS(I/100)//CHRNOS(I/10)//
     &               CHRNOS(MOD(I,10))//'NSKE'// CHRXYZ(IJ)
           indice = indice + 1
        END DO
        END DO
      ENDIF
C
C
C ---------------------------------------------------------------------
C
C      Look for LABELS : SDAV
C     ---------------------------
      IF(KEY.EQ."SDAV") THEN
        NLAB=9
        indice = 1
        DO I=npos1,npos1+2
        DO IJ = 1, 3
           LISTA1(indice)= 'SD '//CHRNOS(I/100)//CHRNOS(I/10)//
     &               CHRNOS(MOD(I,10))//' '//CHRXYZ(-IJ)
           indice = indice + 1
        END DO
        END DO
      ENDIF
C ---------------------------------------------------------------------
C
C      Look for LABELS : ANGPSO
C     ---------------------------
      IF(KEY.EQ."ANGP") THEN
        NLAB=9
        indice = 1
        DO I=npos1,npos1+2
        DO IJ = 1, 3
           LISTA1(indice)= CHRNOS(I/100)//CHRNOS(I/10)//
     &               CHRNOS(MOD(I,10))//'PSOZ'// CHRXYZ(IJ)
           indice = indice + 1
        END DO
        END DO
      ENDIF

C ---------------------------------------------------------------------
C
C      Look for LABELS : FCt0 (Fermy type integrals)
C     ---------------------------
      IF(KEY.EQ."FCt0") THEN
         NLAB=1
         WRITE(LISTA1(1),'(A3,A2,I3.3)') 'FC ',
     &        NAMN(LRATOM)(1:2), LRATOM
      ENDIF


C      Look for LABELS : FCt1
C     ---------------------------
      IF(KEY.EQ."FCt1") THEN
        NLAB=1
        LISTA1(1)='FC1  '//nuc_number
      ENDIF
C
C      Look for LABELS : FCt2
C     ---------------------------
      IF(KEY.EQ."FCt2") THEN
        NLAB=9
        indice = 1
        DO I= 1, 3 !x, y, z
        DO IJ = 1, 3 !x, y ,z
           LISTA1(indice)= 'FC2'//chrxyz(-I)//chrxyz(-IJ)//nuc_number
           indice = indice + 1
        END DO
        END DO
      ENDIF
C
C      Look for LABELS : FCt3
C     ---------------------------
      IF(KEY.EQ."FCt3") THEN
        NLAB=9
        indice = 1
        DO I= 1, 3 !x, y, z
        DO IJ = 1, 3 !x, y, z
           LISTA1(indice)= 'FC3'//chrxyz(-I)//chrxyz(-IJ)//nuc_number
           indice = indice + 1
        END DO
        END DO
      ENDIF
C
C
C      Look for LABELS : NEWJ
C     ---------------------------
      IF(KEY.EQ."NEWJ") THEN
        NLAB=9

!      cmos(i,j)   :  Molecular orbitals coeficients matrix elements (from SIRIUS.RST)
       allocate(xmos(NBAST*NORBT))
       allocate(cmos(NORBT,NBAST))
       CALL RD_SIRIFC('CMO',FOUND,xmos)
       k = 1
       cmos=0.0
       do i=1,nmos                !i corresponds to a MO
          do j=1, naoss           !j corresponds to j coef. of the ith MO
             cmos(i,j) = xmos(k)  ! Transforms vector into matrix
             k = k + 1
          end do
       end do
       deallocate(xmos)

      allocate(aux0(NBAST,NBAST))
      aux0=matmul(transpose(cmos),cmos)
      deallocate(cmos)
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(1)='1 A  21 '
!     Replace time information with symmetry information
      RTNLBL(2)='SQUARE  '

      allocate(AOPERATOR(NBAST,NBAST))
      allocate(BOPERATOR(NBAST,NBAST))
      allocate(aux1(NBAST,NBAST))
      allocate(aux2(NBAST,NBAST))
      allocate(aux3(NBAST,NBAST))

      indice = 1
      DO I= 1, 3 !x, y, z
      DO IJ = 1, 3 !x, y, z
         LISTA1(indice)=nuc_coord(IJ)//'psoz'//chrxyz(I)
         call READAOP(AOPERATOR, naoss, chrxyz(I)//'ANGMOM ') !gets ANGMOM matrix
         call READAOP(BOPERATOR, naoss, 'PSO '//nuc_coord(IJ)//' ') !gets PSO matrix
!!!         call READAOP(BOPERATOR, naoss, chrxyz(I)//'ANGMOM ') !gets ANGMOM matrix !TESTEANDO eRRORRR
!        Construct A.I.B + B.I.A
         aux1 = matmul(matmul(AOPERATOR,aux0),BOPERATOR)    !X/Y/Z ANGMOMG . I . PSO
         aux2 = matmul(matmul(BOPERATOR,aux0),AOPERATOR)    !PSO . I . X/Y/Z ANGMOM
         aux3 = 0.D0 !Initialize/reset matrix
         aux3 = 0.5*(aux1+aux2)

         call WRITEAOPJUAN(aux3,naoss,
     &   nuc_coord(IJ)//'psoz'//chrxyz(I),RTNLBL) !Write Matrix in AOPROPER

         indice = indice + 1
!        For debugging purposes:
!        write(lupri,*) "IMPRIMIENDO EL OPERADOR: "
!     &                 //nuc_coord(IJ)//'psoz'//chrxyz(I)
!        call OUTPUT(aux3,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)  ! prints aux1
      END DO
      END DO

      deallocate(aux0)
      deallocate(AOPERATOR)
      deallocate(BOPERATOR)
      deallocate(aux1)
      deallocate(aux2)
      deallocate(aux3)
      
      ENDIF
C
C ---------------------------------------------------------------------
C  Print LABELS
C -----------------------------------
c      IF(JIMPRT.GE.3) THEN
c         WRITE(lupri,*) '@AVELRSC setting  LABEL :'
c         DO i =1, NLAB
c            WRITE(LUPRI,*)'   LABEL :', LISTA1(I)
c         ENDDO
c      ENDIF
C ---------------------------------------------------------------------
C
C   
C
      KJ = 0
      IDIP = 0
      DO 300 IDIP = 1,NLAB
C
C           3. Work Allocations:
C
         KGD1   = KWORK1
         KWRKG1 = KGD1
         LWRKG1 = LWORK - KWRKG1
         KSLV   = KGD1 + 2*NVARPT
         KLAST  = KSLV + 2*NVARPT
c         WRITE(lupri,*)' NVARPT :', NVARPT
         IF (KLAST.GT.LWORK) 
     &   CALL STOPIT('KLAST GT LWORK on AVELRESC',' ',KLAST,LWORK)
         KWRK = KLAST
         LWRK = LWORK - KLAST + 1
C
         KJ = KJ+1
!           Get average value of operator
            CALL PRP1AVE(LISTA1(IDIP),SNDPRP,WORK(KCMO),WORK(KUDV),
     $      WORK(KWRKG1),LWRKG1,JIMPRT)
            IF (JIMPRT.GE.3) THEN
               WRITE (LUPRI,'(1A,I2,3A,F20.12)')
     &         '     #',KJ,' Expectation Value for Operator : < ',
     &         LISTA1(IDIP), ' > = ', SNDPRP
            ENDIF
C
C       =========================================================
C
C
C              WRITE properties into the various property matrices
C              ===================================================
C
C         <NSCO> : 
C       -----------
         IF (KEY.EQ.'DINR') THEN
            IF(KJ.EQ.1) LRDNR(1,1)=SNDPRP
            IF(KJ.EQ.2) LRDNR(1,2)=SNDPRP
            IF(KJ.EQ.3) LRDNR(1,3)=SNDPRP
            IF(KJ.EQ.4) LRDNR(2,1)=SNDPRP
            IF(KJ.EQ.5) LRDNR(2,2)=SNDPRP
            IF(KJ.EQ.6) LRDNR(2,3)=SNDPRP
            IF(KJ.EQ.7) LRDNR(3,1)=SNDPRP
            IF(KJ.EQ.8) LRDNR(3,2)=SNDPRP
            IF(KJ.EQ.9) LRDNR(3,3)=SNDPRP
         ENDIF
C         <Fc> : 
C       -----------
         IF (KEY.EQ.'FCAV') THEN
            LRFCAV(1,1)= SNDPRP
            LRFCAV(2,2)= SNDPRP
            LRFCAV(3,3)= SNDPRP
         ENDIF
C         <DiaKin> :
C       ------------------
         IF (KEY.EQ.'DIAK') THEN
            IF(KJ.EQ.1) LRDIAK(1,1)=SNDPRP
            IF(KJ.EQ.2) LRDIAK(1,2)=SNDPRP
            IF(KJ.EQ.3) LRDIAK(1,3)=SNDPRP
            IF(KJ.EQ.4) LRDIAK(2,1)=SNDPRP
            IF(KJ.EQ.5) LRDIAK(2,2)=SNDPRP
            IF(KJ.EQ.6) LRDIAK(2,3)=SNDPRP
            IF(KJ.EQ.7) LRDIAK(3,1)=SNDPRP
            IF(KJ.EQ.8) LRDIAK(3,2)=SNDPRP
            IF(KJ.EQ.9) LRDIAK(3,3)=SNDPRP
         ENDIF
C         <Angpso> :
C       ------------------
         IF (KEY.EQ.'ANGP') THEN
            IF(KJ.EQ.1) LRANGP(1,1)=SNDPRP
            IF(KJ.EQ.2) LRANGP(1,2)=SNDPRP
            IF(KJ.EQ.3) LRANGP(1,3)=SNDPRP
            IF(KJ.EQ.4) LRANGP(2,1)=SNDPRP
            IF(KJ.EQ.5) LRANGP(2,2)=SNDPRP
            IF(KJ.EQ.6) LRANGP(2,3)=SNDPRP
            IF(KJ.EQ.7) LRANGP(3,1)=SNDPRP
            IF(KJ.EQ.8) LRANGP(3,2)=SNDPRP
            IF(KJ.EQ.9) LRANGP(3,3)=SNDPRP
         ENDIF
C         < SDAV > :
C       ------------------
         IF (KEY.EQ.'SDAV') THEN
            IF(KJ.EQ.1) LRSDAV(1,1)=SNDPRP
            IF(KJ.EQ.2) LRSDAV(1,2)=SNDPRP
            IF(KJ.EQ.3) LRSDAV(1,3)=SNDPRP
            IF(KJ.EQ.4) LRSDAV(2,1)=SNDPRP
            IF(KJ.EQ.5) LRSDAV(2,2)=SNDPRP
            IF(KJ.EQ.6) LRSDAV(2,3)=SNDPRP
            IF(KJ.EQ.7) LRSDAV(3,1)=SNDPRP
            IF(KJ.EQ.8) LRSDAV(3,2)=SNDPRP
            IF(KJ.EQ.9) LRSDAV(3,3)=SNDPRP
        ENDIF
C         <Fct0> :
C       ------------------
         IF (KEY.EQ.'FCt0') THEN
!            IF(KJ.EQ.1) write(lupri,*) "news:",KEY,SNDPRP*
!     &      calfa*calfa*calfa*calfa*1.0D+6
            LRFCAV0(1,1)=SNDPRP
            LRFCAV0(2,2)=SNDPRP
            LRFCAV0(3,3)=SNDPRP
         ENDIF
C         <Fct1> :
C       ------------------
         IF (KEY.EQ.'FCt1') THEN
!            IF(KJ.EQ.1) write(lupri,*) "news:",KEY,SNDPRP*
!     &      calfa*calfa*calfa*calfa*1.0D+6
!            IF(KJ.EQ.9) LR(3,3)=SNDPRP
            LRFCAV1(1,1)=SNDPRP
            LRFCAV1(2,2)=SNDPRP
            LRFCAV1(3,3)=SNDPRP
         ENDIF
C         <Fct2> :
C       ------------------
         IF (KEY.EQ.'FCt2') THEN
!            IF(KJ.EQ.1) write(lupri,*) "news:",KEY,SNDPRP*
!     &      calfa*calfa*calfa*calfa*1.0D+6
            IF(KJ.EQ.1) LRFCAV2(1,1)=SNDPRP
            IF(KJ.EQ.2) LRFCAV2(1,2)=SNDPRP
            IF(KJ.EQ.3) LRFCAV2(1,3)=SNDPRP
            IF(KJ.EQ.4) LRFCAV2(2,1)=SNDPRP
            IF(KJ.EQ.5) LRFCAV2(2,2)=SNDPRP
            IF(KJ.EQ.6) LRFCAV2(2,3)=SNDPRP
            IF(KJ.EQ.7) LRFCAV2(3,1)=SNDPRP
            IF(KJ.EQ.8) LRFCAV2(3,2)=SNDPRP
            IF(KJ.EQ.9) LRFCAV2(3,3)=SNDPRP
         ENDIF
C         <Fct3> :
C       ------------------
         IF (KEY.EQ.'FCt3') THEN
!            IF(KJ.EQ.1) write(lupri,*) "news:",KEY,SNDPRP*
!     &      calfa*calfa*calfa*calfa*1.0D+6
            IF(KJ.EQ.1) LRFCAV3(1,1)=SNDPRP
            IF(KJ.EQ.2) LRFCAV3(1,2)=SNDPRP
            IF(KJ.EQ.3) LRFCAV3(1,3)=SNDPRP
            IF(KJ.EQ.4) LRFCAV3(2,1)=SNDPRP
            IF(KJ.EQ.5) LRFCAV3(2,2)=SNDPRP
            IF(KJ.EQ.6) LRFCAV3(2,3)=SNDPRP
            IF(KJ.EQ.7) LRFCAV3(3,1)=SNDPRP
            IF(KJ.EQ.8) LRFCAV3(3,2)=SNDPRP
            IF(KJ.EQ.9) LRFCAV3(3,3)=SNDPRP
         ENDIF
C         <NEWJ> :
C       ------------------
         IF (KEY.EQ.'NEWJ') THEN
!             write(lupri,*) 'Corrección ',KJ
!          write(lupri,*) calfa*calfa*calfa*calfa*1.0D+6*CANGP*SNDPRP
            IF(KJ.EQ.1) LRANGPavg(1,1)=SNDPRP
            IF(KJ.EQ.2) LRANGPavg(1,2)=SNDPRP
            IF(KJ.EQ.3) LRANGPavg(1,3)=SNDPRP
            IF(KJ.EQ.4) LRANGPavg(2,1)=SNDPRP
            IF(KJ.EQ.5) LRANGPavg(2,2)=SNDPRP
            IF(KJ.EQ.6) LRANGPavg(2,3)=SNDPRP
            IF(KJ.EQ.7) LRANGPavg(3,1)=SNDPRP
            IF(KJ.EQ.8) LRANGPavg(3,2)=SNDPRP
            IF(KJ.EQ.9) LRANGPavg(3,3)=SNDPRP
         ENDIF
  300 CONTINUE
C


      CALL TIMER ('AVELRSC',TIMEIN,TIMOUT)
C
      CALL QEXIT('AVELRSC')
      RETURN
      END
c
c



C  /* Deck Regenerate_integrals */
      SUBROUTINE REGENERATE_INTEGRALS(COORDENADAS)

#include "implicit.h"
#include "priunit.h"
#include "lrescinf.h"
#include "chrnos.h"
#include "chrxyz.h"

      real(8) :: COORDENADAS(3)
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      character*10 :: namefile
      logical :: file_exists, file_open
      INTEGER :: fileID
      CALL QENTER('REGENERATE_INTEGRALS')
c      CALL GETTIM(TIMHER,WALHER)

      namefile='DALTON.INP'
      fileID = -1
      CALL GPOPEN(fileID,namefile,'UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      write(fileID,'(a)') '**DALTON INPUT'


c      write(fileID,'(a)') '*MOLBAS'
c      write(fileID,'(a)') '.UNCONTRACTED'
c      write(fileID,'(a)') '.NUCMOD'
c      write(fileID,'(a)') ' 1'


      write(fileID,'(a)') '**INTEGRALS'
c     CARE! Some Gauge variable stills unchanged.
      write(fileID,'(a)') '.GAUGEO'
      write(fileID,'(3f11.7)') COORDENADAS
c      write(fileID,'(a)') '.PROPRI'
      write(fileID,'(a)') '.LR-SHI'
c      write(fileID,'(a)') '*LR-SHI'
c      write(fileID,'(a)') '.LRGAUGEO'
c      write(fileID,'(3f11.7)') COORDENADAS
      write(fileID,'(a)') '**END OF DALTON INPUT'

      CALL GPCLOSE(fileID,'KEEP')
      CALL QEXIT('REGENERATE_INTEGRALS')

      RETURN

      END

