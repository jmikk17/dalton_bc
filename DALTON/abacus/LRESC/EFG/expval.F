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
      SUBROUTINE AVELRSC_EFG(KEY,WORK,LWORK)
C
C     Purpose:
C     Get average values from AOPROPER operators in order
C     to get LRESC corrections                  
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
      IPRLNR = JJAPRT
      IPRRSP = -1
cxu
C
C
C  LINEAR LRESC SINGLET ROUTINE
C         'EFG ' <ZZEFG.kin>, <p.q.p> and <(p^2 qzz)>
C
      CALL QENTER('AVELRSC')
      CALL TIMER('START ',TIMEIN,TIMOUT)

      IF (JJAPRT .GE. 2) THEN
         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
         WRITE(LUPRI,*)
      IF (KEY.EQ.'EFG ') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Electric Field Gradient : pqp and kin '
         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
      IF (KEY.EQ.'EFG2') WRITE(LUPRI,'(A)')
     &      '   LINEARLR, Electric Field Gradient : nabla^2(qzz) '
      WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
      END IF
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
      IF (JJAPRT.GT.3) Then
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
         WRITE(lupri,*) '     N2ASHX : ' , N2ASHX
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
      IF (JJAPRT.GT.3) Then
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
      IF (JJAPRT.GT.5) THEN
         WRITE(lupri,*)' CMOS : '
         CALL OUTPUT(WORK(KCMO),1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
         WRITE(lupri,*) '   '
      ENDIF

      ISYM = 1
      IF (JJAPRT.GT.3) Then
         WRITE(lupri,*) '   '
         WRITE(lupri,*) ' about to call LNRVAR'
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '           ISYM   : ', ISYM
         WRITE(lupri,*) '   KWORK1=KFREE   : ', KWORK1
         WRITE(lupri,*) '    WORK1=LFREE   : ', WORK1
         WRITE(lupri,*) '   '
      ENDIF


C     we keep this just in case
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
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
      IF (JJAPRT.GT.3) Then
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

C===============================================================================
C  SINGLET CALCULATIONS :
C ==============================================================================
C       EFG - Using AVG
C --------------------------------
      IF(KEY.EQ."EFG ") THEN
         NLAB=NUCIND         
         IJ = 1
         DO I=1,NUCIND
      LISTA1(IJ)= 'ZZEFG'//CHRNOS(I/10)//CHRNOS(MOD(I,10))//
     &            CHRNOS(1)
           IJ = IJ + 1
         END DO
      END IF

C      Look for LABELS : LAPLACIAN qzz
C     ---------------------------
      IF(KEY.EQ."EFG2") THEN
      NLAB=NUCIND 
      IJ = 1
      DO I=1, NUCIND
            LISTA1(IJ)= 'LEFG '//CHRNOS(0)//CHRNOS(0)//CHRNOS(3*I)
      LISTA1(IJ)= 'LEFG '//CHRNOS(3*I/100)//CHRNOS(MOD(3*I,100)/10)
     &   //CHRNOS(MOD(3*I,10))
      IJ = IJ + 1
      END DO
      ENDIF

C  Print LABELS
C -----------------------------------
      IF(JJAPRT.GE.2) THEN
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
         KJ = KJ+1
!           Get average value of operator
            CALL PRP1AVE(LISTA1(IDIP),SNDPRP,WORK(KCMO),WORK(KUDV),
     $      WORK(KWRKG1),LWRKG1,JJAPRT)
            IF (JJAPRT.GE.2) THEN
               WRITE (LUPRI,'(1A,I2,3A,F20.12)')
     &         '#',KJ,' Expectation Value for Operator : < ',
     &         LISTA1(IDIP), ' > = ',SNDPRP
            ENDIF

			IF (KEY.EQ.'EFG ') THEN 
				EFGC0(IDIP)=SNDPRP
			ENDIF
cx         ENDIF
C
C       =========================================================
C
C
C              WRITE properties into the various property matrices
C              ===================================================
C
C         <nabla^2(qzz)> : (EFG)
C       ------------------
      IF (KEY.EQ.'EFG2') THEN
            EFGC2(KJ,5)=CEFGlap*calfa*calfa*SNDPRP
      ENDIF
  300 CONTINUE
C


C ZZ-Electric Field Gradient Corrections
C --------------------------------
      IF(KEY.EQ."EFG ") THEN
         call EFGdrv(WORK(KCMO),JJAPRT,cmos,OVERLP)
         call AVGchecks(JJAPRT,cmos,OVERLP) !Basis and Tr(P.S) Mulliken checks.
      END IF


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
