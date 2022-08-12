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
c#include "sigma.h"
#include "maxorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "cbilnr.h"
c#include "suscpt.h"
#include "infpri.h"
      double precision cmos(NORBT,NBAST)
      double precision OVERLP(NBAST,NBAST)
      LOGICAL FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1,LISTA1(100)
      CHARACTER*6 LABEL2
      CHARACTER*4 KEY
      CHARACTER*3 char
      DIMENSION SIGMAMV(9),SIGMADW(9)
      PARAMETER (D05=0.5D0,D025=0.25)
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
c#include "infsop.h"

C
#include "lrescinf.h"
#include "chrxyz.h"
#include "chrnos.h"
#include "orgcom.h"
C
      IPRLNR = JIMPRT
      IPRRSP = -1
cxu
C
C
C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'FCAV'  <Fc>
C         'DIAK'  <Dia.Kin>
C         'ANGP'  <L.Pso>
C
      CALL QENTER('AVELRSC')
      CALL TIMER('START ',TIMEIN,TIMOUT)

      IF (JIMPRT .GE. 2) THEN
         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
         WRITE(LUPRI,*)
         IF (KEY.EQ.'FCAV') WRITE(LUPRI,'(A)')
     &      '   Diamagnetic First Order Singlet : FC '
         IF (KEY.EQ.'DIAK') WRITE(LUPRI,'(A)')
     &      '   Diamagnetic First Order Singlet : DIAK'
         IF (KEY.EQ.'ANGP') WRITE(LUPRI,'(A)')
     &      '   Diamagnetic First Order Singlet : ANGP '
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
      IF (JIMPRT.GT.3) Then
         WRITE(lupri,'(A,3F12.8)') ' orgcom.h : GAGORG :', GAGORG
         WRITE(lupri,'(A,3F12.8)') '            ORIGIN :', ORIGIN
         WRITE(lupri,'(A,3F12.8)') '            CMXYZ  :', CMXYZ
         WRITE(lupri,*)
         WRITE(lupri,*) ' alocando 1 : ANTES MEMGET :'
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '     KFREE = 1'
         WRITE(lupri,*) '     LFREE = LWORK :         ', LWORK
         WRITE(lupri,*) '                             '
C
         WRITE(lupri,*) ' COMMON VARIABLES on LINEAR '
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '     N2ASHX : ni idea ' , N2ASHX
         WRITE(lupri,*) '     NASHT # Active Orbitals = 0 ? :', NASHT
         WRITE(lupri,*) '     LISTA1y2 (4*MXCOOR+9) = ', 4*MXCOOR+9
         WRITE(lupri,*) '     MXCOOR              = ', MXCOOR
         WRITE(lupri,*) '     NCMOT = NORB * NORB = ', NCMOT
         WRITE(lupri,*) '   '
         WRITE(lupri,*) ' memget....  '
      ENDIF
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK ,KFREE ,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK ,KFREE ,LFREE)
c                  TYPE, KBASE, LENGTH, WORK, KFREE, LFREE
c            dimensiona work(KCMO, KCMO+NCMOT)
C
      KWORK1 = KFREE
      WORK1  = LFREE
      IF (JIMPRT.GT.3) Then
         WRITE(lupri,*) '   '
         WRITE(lupri,*) ' AFTER MEMGET  '
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '        KCMO, NCMOT     =  ', KCMO,NCMOT
         WRITE(lupri,*) '        KUDV, LUDV      =  ', KUDV,LUDV
         WRITE(lupri,*) '        KPVX, LPVX      =  ', KPVX,LPVX
         WRITE(lupri,*) '        KPXINDX, LCINDX =  ', KXINDX,KXINDX
         WRITE(lupri,*) '        KWORK1 = KFREE :   ', KFREE
         WRITE(lupri,*) '        WORK1  = LFREE :   ', LFREE
         WRITE(lupri,*) '   '
      ENDIF

      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO))
C          RD_SIRIFC( KEY ,FOUND,   AMAT   ,  WRK      ,LWRK)
 
      IF (.NOT.FOUND) 
     & CALL QUIT('***Error AVELRSC: CMO is not in SIRIFC***')
      IF (JIMPRT.GT.5) THEN
         WRITE(lupri,*)' CMOS : '
         CALL OUTPUT(WORK(KCMO),1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
         WRITE(lupri,*) '   '
      ENDIF
cx ACA es para alguna capa activa
cx      IF (NASHT .GT. 0) THEN
cx         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
cx         WRITE(lupri,*)'jim  DV found on RD_SIFC '
cx         IF (.NOT.FOUND)
cx     &      CALL QUIT('ROUTINE error: DV not found on SIRIFC')
cx         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
cx      END IF
C
      ISYM = 1
      IF (JIMPRT.GT.3) Then
         WRITE(lupri,*) '   '
         WRITE(lupri,*) ' about to call LNRVAR'
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '           ISYM   : ', ISYM
         WRITE(lupri,*) '   KWORK1=KFREE   : ', KWORK1
         WRITE(lupri,*) '    WORK1=LFREE   : ', WORK1
         WRITE(lupri,*) '   '
      ENDIF
C not needed for mean value
C      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
C     we keep this just in case
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
C
C     SOPPA :
C
cdx      IF (ABASOP) THEN
C
C        Initialize XINDX
C
cdx         CALL DZERO(WORK(KXINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
cdx         CALL SET2SOPPA(WORK(KXINDX+KABSAD-1),WORK(KXINDX+KABTAD-1),
cdx     *                  WORK(KXINDX+KIJSAD-1),WORK(KXINDX+KIJTAD-1),
cdx     *                  WORK(KXINDX+KIJ1AD-1),WORK(KXINDX+KIJ2AD-1),
cdx     *                  WORK(KXINDX+KIJ3AD-1),WORK(KXINDX+KIADR1-1))
C
C
cdx         REWIND (LUSIFC)
cdx         IF (CCPPA) THEN
cdx            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
cdx         ELSE
cdx            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
cdx         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
cdx         CALL READT (LUSIFC,LPVMAT,WORK(KPVX))
C
cdx         IF (IPRLNR.GT.10) THEN
cdx            IF (CCPPA) THEN
cdx               WRITE(LUPRI,'(/A)')' EXCIT1 : CCSD correlation ',
cdx     &                           'coefficients'
cdx            ELSE
cdx               WRITE(LUPRI,'(/A,A)')' EXCIT1 :',
cdx     &                              ' MP2 correlation coefficients'
cdx            ENDIF
cdx            CALL OUTPUT(WORK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
cdx         END IF
C
C        reads the MP2 or CCSD second order one particle density matrix 
C
cdx         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
cdx         IF (IPRLNR.GT.10) THEN
cdx            IF (CCPPA) THEN
cdx               WRITE(LUPRI,'(/A)')' RSPMC : CCSD density'
cdx            ELSE
cdx               WRITE(LUPRI,'(/A)')' RSPMC : MP2 density'
cdx            END IF
cdx            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
cdx     &                  LUPRI)
cdx         END IF
C
cdx         CALL SOPUDV(WORK(KUDV))
cdx      END IF
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
      IF (JIMPRT.GT.3) Then
         WRITE(lupri,*) '        '
         WRITE(lupri,*) ' stills to allocate : '
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) ' KIDSYM = KWORK1           : ' , KIDSYM
         WRITE(lupri,*) ' KIDADR = KIDSYM + 9MXCENT : ' , KIDADR
         WRITE(lupri,*) ' KWORK2 = KIDADR + 9MXCENT : ' , KWORK2
         WRITE(lupri,*) ' LWORK2 = LWORK - KWORK2   : ' , LWORK2
      ENDIF
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
cb      write(lupri,*)' selected atom is :', LRATOM

C===============================================================================
C  SINGLET CALCULATIONS : D1S and P1S
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
         NLAB=3
         LABEL2='NSKE'
         IJ = 1
         DO I=npos1,npos1+2
            LISTA1(IJ)= CHRNOS(I/100)//CHRNOS(I/10)//
     &                    CHRNOS(MOD(I,10))//'NSKE'//CHRXYZ(IJ)
            IJ = IJ + 1
         END DO
      ENDIF
C
C
C ---------------------------------------------------------------------
C
C      Look for LABELS : ANGPSO
C     ---------------------------
      IF(KEY.EQ."ANGP") THEN
        NLAB=3 
        LABEL2='PSOZ'
        IJ = 1
        DO I=npos1,npos1+2
           LISTA1(IJ)= CHRNOS(I/100)//CHRNOS(I/10)//
     &               CHRNOS(MOD(I,10))//'PSOZ'// CHRXYZ(IJ)
c           WRITE(*,*) LISTA1(IJ)
          IJ = IJ + 1
        END DO
      ENDIF

C
C
C
C   AngPso is done with a unity matrix inside : AngPso = ANGMOM.1.PSO
C
cx      IF(KEY.EQ."ANGP") THEN
cx        NLAB = 3   ! this is to control loop ahead when calling angpso
cxC
cxC        Look for LABELS : ANGPSO
cxC       ---------------------------
cx        LISTA1(1)='XANGMOM'
cx        LISTA1(2)='YANGMOM'
cx        LISTA1(3)='ZANGMOM'
cxC
cxC       Look for LABELS : PSO 
cxC       ----------------------       
cx        LABEL1='PSO'
cx        IJ = 4
cx        DO I=npos1,npos1+2
cx           LISTA1(IJ)= 'PSO '//CHRNOS(I/100)//CHRNOS(I/10)//
cx     &            CHRNOS(MOD(I,10))
cx           IJ = IJ + 1
cx        END DO
cx      ENDIF
C
C ---------------------------------------------------------------------

c        WRITE(lupri,*)' NAOS  ' , NAOS(1)
c        WRITE(lupri,*)' NASHT ' , NASHT
c        WRITE(lupri,*)' NCMOT ' , NCMOT
c        WRITE(lupri,*)'   KCMO  ', KCMO
c        WRITE(lupri,*)'   NORBT ', NORBT
c        WRITE(lupri,*)'   NORB  ', NORB(1)
c        WRITE(lupri,*)'   N2BASX  ', N2BASX 
c        WRITE(lupri,*)'   NNBASX  ', NNBASX 
c        WRITE(lupri,*)'NISHT,NSSHT,NOCCT,NBAST,NRHFT,NVIRT :' ,
c     &                 NISHT,NSSHT,NOCCT,NBAST,NRHFT,NVIRT
c  idea pasarle n2basx a ver si lee bien aoproper.
c        DO j = 1 , NORBT
c           WRITE(lupri,*) j ,  WORK(KCMO+j)
c        ENDDO

C  Print LABELS
C -----------------------------------
      IF(JIMPRT.GE.2) THEN
         WRITE(lupri,*) '@AVELRSC setting  LABEL :'
         DO i =1, NLAB
            WRITE(LUPRI,*)'   LABEL :', LISTA1(I)
         ENDDO
      ENDIF
C
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
C        Starting Calculations depending on angpso or not
C
         KJ = KJ+1
cx         IF (KEY.EQ.'ANGP') THEN
cx            CALL ANGPSO(WORK(KCMO),LWORK,JIMPRT,LISTA1(IDIP),
cx     $      LISTA1(IDIP+3),SNDPRP)
cx            IF (JIMPRT.GE.2) THEN
cx               WRITE (LUPRI,'(1A,I2,5A,F20.12)')
cx     &         '#',KJ,' Expectation Value for Operator : < ',
cx     &         LISTA1(IDIP), '.',LISTA1(IDIP+3),' > = ',SNDPRP
cx            ENDIF
cx         ELSE
!           Get average value of operator
            CALL PRP1AVE(LISTA1(IDIP),SNDPRP,WORK(KCMO),WORK(KUDV),
     $      WORK(KWRKG1),LWRKG1,JIMPRT)
            IF (JIMPRT.GE.2) THEN
               WRITE (LUPRI,'(1A,I2,3A,F20.12)')
     &         '#',KJ,' Expectation Value for Operator : < ',
     &         LISTA1(IDIP), ' > = ',SNDPRP
c            ENDIF
           ENDIF
cx         ENDIF
C
C       =========================================================
C
C
C              WRITE properties into the various property matrices
C              ===================================================
C
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
            IF(KJ.EQ.2) LRDIAK(2,2)=SNDPRP
            IF(KJ.EQ.3) LRDIAK(3,3)=SNDPRP
         ENDIF
C         <Angpso> :
C       ------------------
         IF (KEY.EQ.'ANGP') THEN
            IF(KJ.EQ.1) LRANGP(1,1)=SNDPRP
            IF(KJ.EQ.2) LRANGP(2,2)=SNDPRP
            IF(KJ.EQ.3) LRANGP(3,3)=SNDPRP
         ENDIF
  300 CONTINUE
C


      CALL TIMER ('AVELRSC',TIMEIN,TIMOUT)
C
      CALL QEXIT('AVELRSC')
      RETURN
      END
C...
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
