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
      SUBROUTINE LINEARLR_EFG(KEY,WORK,LWORK)
C
C     Purpose:
C     Get linear response function values in order
C     to obtain LRESC corrections                  
C     Author:
C           Juan J. Aucar         April 2020              
C...
C...  This subroutine was written by Juan Ignacio Melo using
C...  the subroutine ABACTOCD  as a model (2012)
C
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
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL1,LABEL2,LISTA1(4*MXCOOR+9),LISTA2(4*MXCOOR+9)
      CHARACTER*4 KEY
      CHARACTER*3 char
      PARAMETER (D05=0.5D0,D025=0.25)
      LOGICAL TODOINT
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
cxu  sacar lo de abajo
      IPRLNR = JJAPRT
cxu
C
C
c LINEAR LRESC SINGLETS ROUTINE for EFG
c J. Aucar - 2021
C   KEY = 'J-Dw' <<qzz,Dw>>
C         'J-Mv' <<qzz,Mv>>
      CALL QENTER('LRSCLIN')
      CALL TIMER('START ',TIMEIN,TIMOUT)

       IPRRSP = JJAPRT-1

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
      IF (JJAPRT.GT.2) Then 
         WRITE(lupri,'(A,3F12.8)') ' orgcom.h : GAGORG :', GAGORG
         WRITE(lupri,'(A,3F12.8)') '            ORIGIN :', ORIGIN
         WRITE(lupri,'(A,3F12.8)') '            CMXYZ  :', CMXYZ
         WRITE(lupri,*)
         WRITE(lupri,*) ' alocando 1 : BEFORE MEMGET :'
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
      LWORK1 = LFREE
      IF (JJAPRT.GT.5) Then
         WRITE(lupri,*) '   '
         WRITE(lupri,*) ' AFTER MEMGET  '
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '        KCMO, NCMOT     =  ', KCMO,NCMOT
         WRITE(lupri,*) '        KUDV, LUDV      =  ', KUDV,LUDV
         WRITE(lupri,*) '        KPVX, LPVX      =  ', KPVX,LPVX
         WRITE(lupri,*) '        KPXINDX, LCINDX =  ', KXINDX,KXINDX
         WRITE(lupri,*) '        KWORK1 = KFREE :   ', KFREE
         WRITE(lupri,*) '        LWORK1 = LFREE :   ', LFREE
         WRITE(lupri,*) '   '
      ENDIF

      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO))
C          RD_SIRIFC( KEY ,FOUND,   AMAT   ,  WRK      ,LWRK)
      IF (.NOT.FOUND)
     & CALL QUIT('***Error LRSCLIN: CMO is not in SIRIFC')

      ISYM = 1
      IF (JJAPRT.GT.5) Then
         WRITE(lupri,*) '   '
         WRITE(lupri,*) ' about to call LNRVAR'
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) '           ISYM   : ', ISYM
         WRITE(lupri,*) '   KWORK1=KFREE   : ', KWORK1
         WRITE(lupri,*) '   LWORK1=LFREE   : ', LWORK1
         WRITE(lupri,*) '   '
      ENDIF


      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
c tirar esto
      IPRCIX = -1
c
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
      LWORK2 = LWORK1  - KWORK2
      IF (JJAPRT.GT.5) Then
         WRITE(lupri,*) '        '
         WRITE(lupri,*) ' stills allocate : '
         WRITE(lupri,*) ' ----------------------------'
         WRITE(lupri,*) ' KIDSYM = KWORK1           : ' , KIDSYM
         WRITE(lupri,*) ' KIDADR = KIDSYM + 9MXCENT : ' , KIDADR
         WRITE(lupri,*) ' KWORK2 = KIDADR + 9MXCENT : ' , KWORK2
         WRITE(lupri,*) ' LWORK2 = LWORK - KWORK2   : ' , LWORK2
      ENDIF
C

C  blanquear todas las labels antes
       NLAB1 = 0
       NLAB2 = 0 

C ==============================================================================
C
C  Staritng Labels stuff
C
C ==============================================================================

         npos1 = 3*LRATOM-2
C
C===============================================================================
C  EFG CALCULATIONS : Dw, Mv
C ==============================================================================
       IF(KEY.EQ."J-Dw") THEN ! we will use dalton built in integrals
C        -------------------------------------       
C        Look for LABELS : ZZEFG---
         LABEL1='ZZEFG  '
         IJ = 1
         DO I=1,NUCIND
            LISTA1(IJ)= 'ZZEFG'//CHRNOS(I/10)//CHRNOS(MOD(I,10))//
     &      CHRNOS(1)
            IJ = IJ + 1
         END DO
         NLAB1 = NUCIND
C        Look for LABELS : Darwin
C        -------------------------       
         LABEL2='DARWIN  '
           LISTA2(1)='DARWIN  '
         NLAB2 = 1
      END IF
      
      
      
       IF(KEY.EQ."J-Mv") THEN ! we will use dalton built in integrals
C        Look for LABELS : ZZEFG--- 
C        -------------------------------------       
         LABEL1='ZZEFG  '
         IJ = 1
         DO I=1,NUCIND
            LISTA1(IJ)= 'ZZEFG'//CHRNOS(I/10)//CHRNOS(MOD(I,10))//
     &      CHRNOS(1)
            LABSYM(IJ)=1
            IJ = IJ + 1
         END DO
         NLAB1 = NUCIND
C        Look for LABELS : Mass Velocity
C        -------------------------       
         LABEL2='MASSVELO'
           LISTA2(1)='MASSVELO'
         NLAB2 = 1
      END IF



C
C  Print Section for LABELS on LISTAS
C  ----------------------------------------------
C
      IF (JJAPRT.GE.2) THEN
         WRITE(LUPRI,*)' @LinearLR setting  LABEL1 :'
         DO i =1, nlab1
            WRITE(LUPRI,*)'   LABEL1 :', LISTA1(I)
         ENDDO
         WRITE(LUPRI,*)
C 
         WRITE(LUPRI,*)' @LinearLR setting  LABEL2 :'
         DO i =1, nlab2
            WRITE(LUPRI,*)'   LABEL2 :', LISTA2(I)
         ENDDO
      ENDIF
C
C
C ---------------------------------------------------------------------
C
C     Set variables for ABARSP and logicals
C
      CICLC  = .FALSE.      ! TRUE for CI calculations
      HFCLC  = NASHT .LE. 1 ! .T. RHF-closed shell or 1e in one active orbital
      TRIPLE = .FALSE.      ! .T. for triplet perturbation operators
      EXECLC = .FALSE.      ! false for linear response equations
      IF(KEY.EQ."J-Dw".OR.KEY.EQ."J-Mv") THEN
         NABATY = 1      ! = 1 for real operators .. -1 for imm. op.
      END IF

      NABAOP = 1      ! number of right hand sides. dejarlo asi . solo 1
C
C     Zero the property tensors
cdx      IF (MAGSUS) CALL DZERO(SUSDZD,9)

C
C        Loop over the right operators which are the
C        the dipole velocity operators
C        ===========================================
C
      LUSOVE = 456
      LUGDVE = 457
      LUREVE = 458
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      KJ = 0
      DO 300 IDIP = 1,NLAB1
         LABEL1 = LISTA1(IDIP)
         ISYM = 1   ! ISYM deberia ser 1 ...
c        set variables for response module
         IF(JJAPRT.GT.3) THEN
            WRITE(lupri,*) ' about to call LNRVAR'
            WRITE(lupri,*) ' ----------------------------'
         ENDIF
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)

C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            IF (KLAST.GT.LWORK)
     &       CALL STOPIT('KLAST GT LWORK on LINEARLR',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1

cx            WRITE(lupri,*) ' KLAST ',KLAST
cx            WRITE(lupri,*) ' LWORK ',LWORK
cx            WRITE(lupri,*) ' IF KLAST GT LWORK you will get an error '

C
C           Find right hand side for right operator and WRITE to file
C           =========================================================
C
            KSYMOP = ISYM
            TRPLET = .FALSE.
            ANTSYM = 0                   !juan Aucar ANTSYM = -1 if triplet (chequear) !!
C           ANTSYM : matrix symmetry of PRPMO matrix
C          (1: symmetric, -1: antisymmetric, 0: unknown)

!
! Gradient Property Vector
!
            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
            IF (JJAPRT.GE.3) THEN
               WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL1
               CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
            ENDIF
C
C           Calculate eigenvector and WRITE to file
C           =======================================
C
            CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &            FRVAL,NFRVAL,NABATY,NABAOP,LABEL1,LUGDVE,LUSOVE,
     &            LUREVE,THCLNR,MAXITE,IPRRSP,MXRM,MXPHP,
     &            WORK(KWRK),LWRK)
C
C           Loop over the left side  property operators
C           ===========================================
C
            DO 200 IPL = 1, NLAB2
C
C              Find label and symmetry of the left side operator
C
               LABEL2 = LISTA2(IPL)
               KSYM   = 1
               KJ = KJ +1
cb               WRITE(lupri,*) '     else KSYM   = 1      ', KSYM 
C
C              If symmetry of right operator equals symmetry of
C              the left operator, that is if ISYM = KSYM, then
C              ================================================
C              (otherwise 2. order property SNDPRP is zero)
C
cx             IF (KSYM.EQ.ISYM) THEN
               KSYMOP = ISYM
               TRPLET = .FALSE.
C
C                 Find right hand side for left operator
C                 ========================================
C
               CALL GETGPV(LABEL2,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &             WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
C
               IF (JJAPRT.GT.3) THEN
                  WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL2
                  CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
               ENDIF
C
C                 Form second order properties SNDPRP
C                 ===================================
C
               REWIND LUSOVE
               CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
               IF (JJAPRT.GT.3) THEN
                  WRITE (LUPRI,'(2A)') 'Solution Vector, label: ',LABEL1
                  CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,NVARPT,2,1,LUPRI)
               ENDIF
C
               SNDPRP = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
               IF (JJAPRT.GE.2) THEN
               WRITE (LUPRI,'(1A,I2,5A,F50.12)')
     &         '#',KJ,' Response of operators: <<',LABEL2,';',LABEL1,
     &            '>> = ',SNDPRP
               ENDIF
C
C              WRITE properties into the various property matrices
C              ===================================================
C
C         <<Mv,ZZEFGabc>>
C       ------------------
      IF(KEY.EQ."J-Mv") THEN
            EFGC2(KJ,1)=CEFGMV*SNDPRP
      ENDIF

C         <<Dw,ZZEFGabc>>
C       ------------------
      IF(KEY.EQ."J-Dw") THEN
            EFGC2(KJ,2)=CEFGDW*SNDPRP
      ENDIF
  200 CONTINUE
cm         END IF  NFRVAL
  300 CONTINUE

      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
C
C
      IF (JJAPRT.GT.10) THEN
         WRITE(LUPRI,*)
         CALL TIMER ('LRSCLIN',TIMEIN,TIMOUT)
      ENDIF
C
      CALL QEXIT('LRSCLIN')
      RETURN
      END
C...
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
