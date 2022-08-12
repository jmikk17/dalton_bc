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
C
C  /* Deck cprhsn */
      SUBROUTINE CPRHSN(RHST1,RHST2,T1AM,T2AM,
     *                  WORK,LWORK)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHSN written by Henrik Koch 25-Sep-1993
C     to calculate vector for obtaining corrections to CP amplitudes
C     by using the AO-integrals directly from disk.
C
C
C
C     NB! It is assumed that the vectors are allocated in the following
C     order:
C           T1AM(*), OMEGA1(*), OMEGA2(*), T2AM(*),  WORK(*).
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
#include "ccisao.h"
#include "blocks.h"
#include "ccfield.h"
#include "ccsections.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "eribuf.h"
#include "ccnoddy.h"
#include "cbirea.h"
#include "r12int.h"
#include "ccr12int.h"
#include "qm3.h"
!#include "qmmm.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (XMHALF = -0.5D0, XMONE= -1.0D0 )
      PARAMETER (ISYM0 = 1)
C
      LOGICAL CC1BSA
C
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION RHST1(*),RHST2(*),T1AM(*),T2AM(*),WORK(LWORK)
C
      CHARACTER CFIL*6,DFIL*6, FN3SRT*8, FNDELD*6, CDUMMY*8
      CHARACTER FNCKJD*6, FNDKBC*4, FNTOC*8, FN3VI*6, FN3VI2*8
      CHARACTER FNIADJ*8, FNIJDA*8, CPFIL*8, DPFIL*8
      CHARACTER MODEL*10, DFILg*7
C
      PARAMETER (FNIADJ = 'CCXIADJ0', FNIJDA = 'CCXIJDA0')
      PARAMETER (CPFIL  = 'CC_CPR12', DPFIL  = 'CC_DPR12')
C

      INTEGER IGLMRHS(8,8),IGLMVIS(8,8),NGLMDS(8),ICMO(8,8),NCMO(8),
     &        IMAIJM(8,8),NMAIJM(8),
     &        IMATIJM(8,8),NMATIJM(8),NGAMSM(8),IGAMSM(8,8),
     &        IRGIJS(8,8),NRGIJS(8),IR1BASM(8,8),NR1BASM(8),
     &        IR2BASM(8,8),NR2BASM,IR1XBASM(8,8),NR1XBASM(8),
     &        IR2XBASM(8,8),IMATF(8,8),NMATF(8),IMAKLM(8,8),NMAKLM(8)
      INTEGER NADP(8),IADP(8,8),NLAMDX(8),ILAMDX(8,8)
C
      REAL*8, ALLOCATABLE :: DENMAT(:), FOCKMAT(:), FOCKTEMP(:)
C
      CALL QENTER('CCRHSN')
      NEWGAM = .FALSE.
      OMEGSQ = .FALSE.
      OMEGOR = .FALSE.
C
C-----------------------------------------------------------
C     For energy calculation trial vector is totalsymmetric.
C-----------------------------------------------------------
C
      ISYMTR = 1
C
C-----------------------------------------
C     Save CC1B flag and if CC1A set true.
C-----------------------------------------
C
      CC1BSA = CC1B
      IF ( CC1A ) CC1B = .TRUE.
C
      IF ( IPRINT .GT. 10 ) THEN
C
         WRITE(LUPRI,*) ' In cpsd_rhs : '
         WRITE(LUPRI,*) ' CCSD, CC2: ',CCSD,CC2
         WRITE(LUPRI,*) ' CC1A, CC1B, CC3: ', CC1A, CC1B, CC3
C
      ENDIF
C
C----------------
C     Open files.
C----------------
C
      LUC = -1
      LUD = -1
      LUDg = -1
      CFIL = 'PMAT_C'
      DFIL = 'PMAT_D'
      DFILg = 'PMAT_Dg'
C
      IF (DEBUG) WRITE(LUPRI,*) 'DUMPCD = ',DUMPCD
      IF (DUMPCD) THEN
         CALL WOPEN2(LUC,CFIL,64,0)
         CALL WOPEN2(LUD,DFIL,64,0)
         CALL WOPEN2(LUDg,DFILg,64,0)
C
      END IF

C
C----------------------------------
C     Initialize timing parameters.
C----------------------------------
C
      TIMALL  = SECOND()
      TIMA    = 0.0D00
      TIMB    = 0.0D00
      TIMC    = 0.0D00
      TIMD    = 0.0D00
      TIMG    = 0.0D00
      TIMH    = 0.0D00
      TIMGAM  = 0.0D00
      TIMLAM  = 0.0D00
      TIMRDAO = 0.0D00
      TIMHER1 = 0.0D00
      TIMHER2 = 0.0D00
      TIMT2AO = 0.0D00
      TIMFCK  = 0.0D00
      TIMDM   = 0.0D00
      TIMT2TR = 0.0D00
      TIMTRBT = 0.0D00
C
C---------------------------
C     Check inconsistencies.
C---------------------------
C
      IF (NEWGAM) THEN
         IF ((.NOT. DUMPCD) .OR. (.NOT. OMEGOR)) THEN
            WRITE(LUPRI,*) 'NEWGAM requires both DUMPCD and OMEGOR'
            CALL QUIT('ERROR: NEWGAM inconsistency')
         END IF
      END IF
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      KLAMDP = 1
      KLAMIP = KLAMDP + NLAMDT
      KLAMDH = KLAMIP + 1
      KGAMMA = KLAMDH + NLAMDT
      KEND1 = KGAMMA + NGAMMA(ISYMOP)
      KE1PIM = KEND1
C
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHSN')
      ENDIF
C
C------------------------------------
C     Save the CC amplitudes on disk.
C------------------------------------
C
      LURHS1 = -1
      CALL GPOPEN(LURHS1,'CCRHS1','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURHS1)
      WRITE (LURHS1) (T1AM(I), I = 1,NT1AMX)
      WRITE (LURHS1) (T2AM(I), I = 1,NT2AMX)
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
      TIMLAM  = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)
      TIMLAM  = SECOND() - TIMLAM
C
C-------------------------------
C     Prepare the t2-amplitudes.
C-------------------------------
C
      CALL DCOPY(NT2AMX,T2AM,1,RHST2,1)
      CALL CC_T2SQ(RHST2,T2AM,ISYMTR)
C
C-----------------------------------------
C     Construct the transposed amplitudes.
C-----------------------------------------
C
      IF ((.NOT. DIRECT) .AND. T2TCOR) THEN
C
         KT2AMT = KEND1
         KEND1  = KT2AMT + NT2SQ(1)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient core in CCRHSN')
         END IF
C
         JSYM = 1
         CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
         CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND1),LWRK1,JSYM)
C
      END IF
C
C-------------------------------
C     Initialize OMEGA1 & OMEGA2
C-------------------------------
C
      CALL DZERO(RHST1,NT1AM(ISYMOP))
      CALL DZERO(RHST2,NT2AO(ISYMOP))
C
C----------------------
C     Initialize GAMMA.
C----------------------
C
      IF (.NOT. NEWGAM) CALL DZERO(WORK(KGAMMA),NGAMMA(ISYMOP))
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed t2am.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),
     &                        WORK(KODBC1),WORK(KODBC2),
     &                        WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),
     &                        WORK(KRDPP1),WORK(KRDPP2),
     &                        WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CCRHSN')
               END IF
C
               IF (T2TCOR) THEN
                  KT2AMT = KEND1
                  KEND1  = KT2AMT + NT2SQ(1)
                  LWRK1  = LWORK  - KEND1
                  IF (LWRK1 .LT. 0) THEN
                     CALL QUIT('Insufficient core in CCRHSN')
                  END IF
C
                  JSYM = 1
                  CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
                  CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND1),LWRK1,JSYM)
               END IF
C
            ELSE
               NUMDIS = 1
               KRECNR = KENDSV
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               IT2DEL(IDEL) = ICDEL1
               ICDEL1 = ICDEL1 + NT2BCD(ISYDIS)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
C
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
C
C------------------------------------------
C              Work space allocation no. 3.
C------------------------------------------
C
               KSCRM = KEND2
               KEND3 = KSCRM + NT2BCD(ISYMD)
               LWRK3 = LWORK - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
C
C----------------------------------------------------------------
C              Construct the partially transformed T2-amplitudes.
C----------------------------------------------------------------
C
               DTIME   = SECOND()
               ICON = 1
               ISYMLH = 1
               CALL CC_T2AO(T2AM,WORK(KLAMDH),ISYMLH,WORK(KSCRM),
     *                         WORK(KEND3),LWRK3,IDEL,ISYMD,
     *                         ISYMTR,ICON)
               DTIME   = SECOND() - DTIME
               TIMT2AO = TIMT2AO + DTIME
C
C-----------------------------------
C              Calculate the B-term.
C-----------------------------------
C
               DTIME   = SECOND()
               CALL CCRHS_B(WORK(KXINT),RHST2,WORK(KLAMDP),
     *                         WORK(KLAMDH),WORK(KSCRM),WORK(KEND3),
     *                         LWRK3,IDEL,ISYMD)
               DTIME   = SECOND() - DTIME
               TIMB    = TIMB     + DTIME
C
C------------------------------------------
C              Work space allocation no. 4.
C------------------------------------------
C
               KDSRHF = KEND3
               KEND4  = KDSRHF + NDSRHF(ISYMD)
               LWRK4  = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND4,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
C
C--------------------------------------------------------
C              Transform one index in the integral batch.
C--------------------------------------------------------
C
               DTIME   = SECOND()
               ISYMLP  = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                     ISYMLP,WORK(KEND4),LWRK4,ISYDIS)
               DTIME   = SECOND() - DTIME
               TIMTRBT = TIMTRBT + DTIME
C
C-------------------------------------------------------------
C              Calculate the gamma matrix entering the A-term.
C-------------------------------------------------------------
C
               DTIME   = SECOND()
               CALL CPRHS_GAM(WORK(KDSRHF),WORK(KGAMMA),WORK(KLAMDP),
     *                           WORK(KLAMDH),WORK(KSCRM),WORK(KEND4),
     *                           LWRK4,IDEL,ISYMD)
               DTIME   = SECOND() - DTIME
               TIMGAM  = TIMGAM   + DTIME
C
C-----------------------------------
C              Calculate the C-term.
C-----------------------------------
C
               DTIME   = SECOND()
C
               FACTC = XMHALF
C
               ICON = 2
               IV = 1
C
               IOPTR12 = 0
               IOPTE = 0
C
               CALL CPRHS_C(WORK(KXINT),WORK(KDSRHF),RHST2,
     *                         WORK(KT2AMT),ISYMOP,
     *                         WORK(KLAMDP),WORK(KLAMIP),
     *                         WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                         WORK(KLAMDP),ISYMTR,
     *                         WORK(KSCRM),WORK(KE1PIM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                         IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
CTesT
C              WRITE(LUPRI,*) 'E1PIM after CCRHS_C:'
C              WRITE(LUPRI,*) 'Norm^2: ',
C    &                       DDOT(NADP(1),WORK(KE1PIM),1,WORK(KE1PIM),1)
C              DO ISYM = 1,NSYM
C                CALL OUTPUT(WORK(KE1PIM+IADP(ISYM,ISYM)),
C    &                       1,NVIR(ISYM),1,MBAS1(ISYM)+MBAS2(ISYM),
C    &                       NVIR(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
C    &                       1, LUPRI)
C              END DO
C              CALL FLSHFO(LUPRI)
CTesT
C
               DTIME   = SECOND() - DTIME
               TIMC    = TIMC     + DTIME
C
C-----------------------------------------------
C              Transform the cluster amplitudes.
C-----------------------------------------------
C
               CALL CC_MTCME(WORK(KSCRM),WORK(KEND4),LWRK4,
     *                       ISYMD,ISYMTR)
C
C-----------------------------------
C              Calculate the D-term.
C-----------------------------------
C
               DTIME   = SECOND()
C
               FACTD = HALF
C
               ICON   = 2
               IV = 1
C
               IOPTR12 = 0
               IOPTE = 0
C
               DTIME   = SECOND() - DTIME
               CALL CPRHS_D(WORK(KXINT),WORK(KDSRHF),RHST2,T2AM,
     *                         ISYMTR,WORK(KLAMDP),WORK(KLAMIP),
     *                         WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                         WORK(KLAMDH),ISYMTR,
     *                         WORK(KSCRM),WORK(KE1PIM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTD,ICON,IOPTR12,
     *                         IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
               TIMD    = TIMD     + DTIME
               DTIMEg   = SECOND() - DTIMEg
               CALL CPRHS_Dg(WORK(KXINT),WORK(KDSRHF),RHST2,T2AM,
     *                         ISYMTR,WORK(KLAMDP),WORK(KLAMIP),
     *                         WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                         WORK(KLAMDH),ISYMTR,
     *                         WORK(KSCRM),WORK(KE1PIM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTD,ICON,IOPTR12,
     *                         IOPTE,LUDg,DFILg,LUDP,DPFIL,IV)
               TIMDg    = TIMDg     + DTIMEg
CTesT
C              WRITE(LUPRI,*) 'E1PIM after CCRHS_D:'
C                            WRITE(LUPRI,*) 'Norm^2: ',
C    &                       DDOT(NADP(1),WORK(KE1PIM),1,WORK(KE1PIM),1)
C              DO ISYM = 1,NSYM
C                CALL OUTPUT(WORK(KE1PIM+IADP(ISYM,ISYM)),
C    &                       1,NVIR(ISYM),1,MBAS1(ISYM)+MBAS2(ISYM),
C    &                       NVIR(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
C    &                       1, LUPRI)
C              END DO
C              CALL FLSHFO(LUPRI)
CTesT
C
C-----------------------------------
C              Calculate the G-term.
C-----------------------------------
C
               DTIME   = SECOND()
               ISYMP1 = 1
               ISYMH1 = 1
               CALL CCRHS_G(WORK(KDSRHF),RHST1,WORK(KLAMDP),ISYMP1,
     *                      WORK(KLAMDH),ISYMH1,WORK(KSCRM),WORK(KEND4),
     *                      LWRK4,ISYDIS,ISYMD,ISYMTR)
               DTIME   = SECOND() - DTIME
               TIMG    = TIMG     + DTIME
C
C-----------------------------------
C              Calculate the H-term.
C-----------------------------------
C
               DTIME   = SECOND()
               CALL CCRHS_H(WORK(KDSRHF),RHST1,WORK(KLAMDP),
     *                      WORK(KLAMDH),WORK(KSCRM),WORK(KEND4),
     *                      LWRK4,ISYDIS,ISYMD,ISYMTR)
               DTIME   = SECOND() - DTIME
               TIMH    = TIMH     + DTIME
C
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
CTesT
C        WRITE(LUPRI,*) 'E1PIM after CCSDR12AO:'
C        WRITE(LUPRI,*) 'Norm^2: ',
C    &                 DDOT(NADP(1),WORK(KE1PIM),1,WORK(KE1PIM),1)
C        DO ISYM = 1,NSYM
C          CALL OUTPUT(WORK(KE1PIM+IADP(ISYM,ISYM)),
C    &                 1,NVIR(ISYM),1,MBAS1(ISYM)+MBAS2(ISYM),
C    &                 NVIR(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
C    &                 1, LUPRI)
C        END DO
C        CALL FLSHFO(LUPRI)
CTesT
C
C------------------------
C     Recover work space.
C------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('After  Delta Loop: RHST1')
C         CALL CC_PRP(OMEGA1,OMEGA2,1,1,0)
      ENDIF
C
C-------------------------------------------------
C     Transform the Omega2 vector to the MO basis.
C-------------------------------------------------
C
      IF (NT2AM(ISYMOP) .GT. 2*NT2AMX) THEN
         WRITE(LUPRI,*)
     &        'Length of T2AM is smaller than RHST2 in MO basis'
         CALL QUIT('Insufficient space in CC_T2MO')
      ENDIF
C
      IF ( .NOT. CC2 ) THEN
C
C---------------------------------------
C        Save the CC amplitudes on disk.
C---------------------------------------
C
         WRITE (LURHS1) (T2AM(I), I = 1,NT2AM(ISYMOP))
C
C----------------------------------------------------
C        Transform the Omega2 vector to the MO basis.
C----------------------------------------------------
C
         IF (NT2AM(ISYMOP) .GT. 2*NT2AMX) THEN
            WRITE(LUPRI,*)
     *        'Length of T2AM is smaller than RHST2 in AO basis'
            CALL QUIT('Insufficient space in CC_T2MO')
         ENDIF
C
         TIMOME2 = SECOND()
         ISYMBF = ISYMOP
         ICON = 1
         
         CALL CC_T2MO(FAKE,PHONEY,ISYMOP,RHST2,T2AM,WORK(KGAMMA),
     *                WORK(KLAMDP),WORK(KLAMDP),ISYMTR,
     *                WORK(KEND1),LWRK1,ISYMBF,ICON)
         CALL DCOPY(NT2AM(ISYMTR),T2AM,1,RHST2,1)
         TIMOME2 = SECOND() - TIMOME2
C
C--------------------------------------------
C        Restore the CC amplitudes from disk.
C--------------------------------------------
C
         REWIND (LURHS1)
         READ (LURHS1)
         READ (LURHS1)
         READ (LURHS1) (T2AM(I), I = 1,NT2AM(ISYMOP))
C
      ENDIF

C
C---------------------
C     Reallocate T2TP.
C---------------------
C
      IF (DIRECT .AND. T2TCOR) THEN
C
         KT2AMT = KEND1
         KEND2  = KT2AMT + NT2SQ(1)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2. LT. 0) THEN
            CALL QUIT('Insufficient memory in CCSD_RHS')
         END IF
C
      ELSE
C
         KEND2 = KEND1
         LWRK2 = LWRK1
C
      END IF
C
C---------------------
C     Reallocate T2TP.
C---------------------
C
      IF ((DIRECT .AND. T2TCOR) .OR. (CCSDT .AND. T2TCOR)) THEN
C
         KT2AMT = KEND1
         KEND2  = KT2AMT + NT2SQ(1)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2. LT. 0) THEN
            CALL QUIT('Insufficient memory in CCSD_RHS')
         END IF
C
      ELSE
C
         KEND2 = KEND1
         LWRK2 = LWRK1
C
      END IF
C
C----------------------
C     Recalculate T2TP.
C----------------------
C
      IF (T2TCOR) THEN
C
         JSYM = 1
         CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
         CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND2),LWRK2,JSYM)
C
      END IF
C
C----------------------
C     Calculate A-term.
C----------------------
C
      IOPT = 1
      TIMA     = SECOND()
      CALL CCRHS_A(RHST2,T2AM,WORK(KGAMMA),WORK(KEND2),LWRK2,
     *                ISYMTR,ISYMTR,IOPT)
      TIMA     = SECOND() - TIMA
c
C--------------------------------------
C     If (DUMPCD) calculate the C-term.
C--------------------------------------
C
      IF (DUMPCD .AND. (.NOT. CC2)) THEN
C
         ISYVEC = 1
         ISYCIM = 1
         IOPT   = 1
         IVECNR = 1
C
         TIMCIO = SECOND()
         CALL CPRHS_CIO(RHST2,WORK(KT2AMT),WORK(KLAMDH),
     *                     WORK(KEND2),LWRK2,ISYVEC,ISYCIM,
     *                     LUC,CFIL,IVECNR,IOPT)
C
         TIMCIO  = SECOND() - TIMCIO
C
      ENDIF
C
C---------------------
C     Scale T2 to 2T2.
C---------------------
C
      DTIME    = SECOND()
      CALL DSCAL(NT2SQ(1),TWO,T2AM,1)

      DTIME    = SECOND() - DTIME
      TIMT2TR  = TIMT2TR + DTIME
C
C--------------------------------------------------
C     If (DUMPCD) calculate the D-term and Dg-term.
C--------------------------------------------------
C
      IF (DUMPCD .AND. (.NOT. CC2)) THEN
C
         ISYDIM = 1
         ISYVEC = 1
         IOPT = 1
         IVECNR = 1
C
         TIMDIO = SECOND()
         CALL CCRHS_DIO(RHST2,T2AM,WORK(KLAMDH),WORK(KEND2),LWRK2,
     *                  ISYVEC,ISYDIM,LUD,DFIL,IVECNR,IOPT)
         TIMDIO  = SECOND() - TIMDIO
         TIMDIOg = SECOND()
         CALL CCRHS_DIO(RHST2,WORK(KT2AMT),WORK(KLAMDH),WORK(KEND2),
     *                  LWRK2,ISYVEC,ISYDIM,LUDg,DFILg,IVECNR,IOPT)
         TIMDIOg  = SECOND() - TIMDIOg
      END IF
C
C------------------------
C     Scale final result.
C------------------------
C
C     CALL DSCAL(NT1AM,TWO,OMEGA1,1)
C     CALL DSCAL(NT2IND,TWO,OMEGA2,1)
C
      IF (IPRINT .GT. 25) THEN
         CALL AROUND('END OF CPRHS:RHS T1')
         DO 300 ISYM = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
            KOFF = IT1AM(ISYM,ISYM) + 1
C            CALL OUTPUT(OMEGA1(KOFF),1,NVIR(ISYM),1,NRHF(ISYM),
C     *                  NVIR(ISYM),NRHF(ISYM),1,LUPRI)
  300    CONTINUE
         WRITE(LUPRI,*)
         CALL AROUND('END OF CPRHS:RHS T2')
         DO 310 ISYM = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
            KOFF = IT2AM(ISYM,ISYM) + 1
            CALL OUTPAK(RHST2(KOFF),NT1AM(ISYM),1,LUPRI)
  310    CONTINUE
      ENDIF
      TIMALL  = SECOND() - TIMALL
      IF ( IPRINT .GT. 2) THEN
         WRITE(LUPRI,9999) 'RHS - TOTAL', TIMALL
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,9999) 'CCRHS_A    ', TIMA
         WRITE(LUPRI,9999) 'CCRHS_B    ', TIMB
         WRITE(LUPRI,9999) 'CCRHS_C    ', TIMC
         WRITE(LUPRI,9999) 'CCRHS_CIO  ', TIMCIO
         WRITE(LUPRI,9999) 'CCRHS_C-tot', TIMCIO + TIMC
         WRITE(LUPRI,9999) 'CCRHS_D    ', TIMD
         WRITE(LUPRI,9999) 'CCRHS_Dg   ', TIMDg
         WRITE(LUPRI,9999) 'CCRHS_DIO  ', TIMDIO
         WRITE(LUPRI,9999) 'CCRHS_D-tot', TIMDIO + TIMD
         WRITE(LUPRI,9999) 'CCRHS_DIOg ', TIMDIOg
         WRITE(LUPRI,9999) 'CCRHS_D-totg', TIMDIOg + TIMDg
         WRITE(LUPRI,9999) 'CCRHS_G    ', TIMG
         WRITE(LUPRI,9999) 'CCRHS_H    ', TIMH
         WRITE(LUPRI,9999) 'CCRHS_GAM  ', TIMGAM
         WRITE(LUPRI,9999) 'CCRHS_LAM  ', TIMLAM
         WRITE(LUPRI,9999) 'CCRHS_RDAO ', TIMRDAO
         WRITE(LUPRI,9999) 'HERDIS1    ', TIMHER1
         WRITE(LUPRI,9999) 'HERDIS2    ', TIMHER2
         WRITE(LUPRI,9999) 'CC_T2AO    ', TIMT2AO
         WRITE(LUPRI,9999) 'CCRHS_FCK  ', TIMFCK
         WRITE(LUPRI,9999) 'CCRHS_DM   ', TIMDM
         WRITE(LUPRI,9999) 'CCRHS_TRBT ', TIMTRBT
         WRITE(LUPRI,9999) 'CCRHS_T2TR ', TIMT2TR
      ENDIF
9999  FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
C
C-----------------------------------------
C     Restore the CC amplitudes from disk.
C-----------------------------------------
C
      REWIND (LURHS1)
      READ(LURHS1) (T1AM(I), I = 1,NT1AMX)
      READ(LURHS1) (T2AM(I), I = 1,NT2AMX)
      CALL GPCLOSE(LURHS1,'DELETE')
C
C-----------------
C     Close files.
C-----------------
C
      IF (DUMPCD) THEN
         CALL WCLOSE2(LUC,CFIL,'KEEP')
         CALL WCLOSE2(LUD,DFIL,'KEEP')
         CALL WCLOSE2(LUDg,DFILg,'KEEP')
      END IF
C
C
C-----------------------
C     Restore CC1B flag.
C-----------------------
C
      CC1B = CC1BSA
C
      CALL QEXIT('CCRHSN')
C
      RETURN
      END
C
C
C#####################################################################
C                                                                    #
C     The following subroutines are CCRHS subroutines modified       #
C     to give only the contributions required by the CPSD_ENERGY     #
C     subroutine.                                                    #
C                                                                    #
C#####################################################################
C
C
      SUBROUTINE CPRHS_GAM(DSRHF,GAMMA,XLAMDP,XLAMDH,SCRM,
     *                     WORK,LWORK,IDEL,ISYMD)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_GAM written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 21-July-1994
C
C     Purpose: Calculate the gamma intermediate.
C     Adaptation: Remove double amplitude correction from gamma
C                 so it does not contribute in the A-term.
C
#include "implicit.h"
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*)
      DIMENSION WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KLAMDA = 1
      KEND1  = KLAMDA + NRHF(ISYMD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_GAM')
      ENDIF
C
C---------------------------------------
C     Copy XLAMDH vector for given IDEL.
C---------------------------------------
C
      KOFF1 = ILMRHF(ISYMD) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NRHF(ISYMD),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KLAMDA),1)
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMAG = MULD2H(ISYML,ISYDIS)
C
C---------------------------
C        Dynamic allocation.
C---------------------------
C
         KSCR1  = KEND1
         KSCR2  = KSCR1  + N2BST(ISYMAG)
         KSCR3  = KSCR2  + NT1AO(ISYMAG)
         KSCR4  = KSCR3  + NT1AM(ISYMAG)
         KEND2  = KSCR4  + NMATIJ(ISYMAG)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_GAM')
         ENDIF
C
         CALL CPRHS_GAM1(DSRHF,GAMMA,SCRM,WORK(KLAMDA),XLAMDP,XLAMDH,
     *                   WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),
     *                   WORK(KSCR4),WORK(KEND2),LWRK2,ISYMD,ISYML,
     *                   ISYMAG)
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CPRHS_GAM1(DSRHF,GAMMA,SCRM,XLAM,
     *              XLAMDP,XLAMDH,SCR1,SCR2,SCR3,SCR4,WORK,
     *              LWORK,ISYMD,ISYML,ISYMAG)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate the gamma intermediate.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*),XLAM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*),WORK(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMKC = ISYMAG
C
      DO 100 L = 1,NRHF(ISYML)
C
         KOFF1 = IDSRHF(ISYMAG,ISYML) + NNBST(ISYMAG)*(L - 1) + 1
C
         CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,SCR1)
C
         DO 110 ISYMG = 1,NSYM
C
            ISYMA = MULD2H(ISYMG,ISYMAG)
            ISYMK = ISYMA
            ISYMC = ISYMG
            ISYMI = ISYMG
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASG = MAX(NBAS(ISYMG),1)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            KOFF2 = ILMRHF(ISYMK) + 1
            KOFF3 = IAODIS(ISYMA,ISYMG) + 1
            KOFF4 = IT1AOT(ISYMK,ISYMG) + 1
C
            CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMA),
     *                 ONE,XLAMDP(KOFF2),NBASA,SCR1(KOFF3),NBASA,
     *                 ZERO,SCR2(KOFF4),NRHFK)
C
            KOFF5 = ILMVIR(ISYMC) + 1
            KOFF6 = IT1AMT(ISYMK,ISYMC) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),NBAS(ISYMG),
     *                 ZERO,SCR2(KOFF4),NRHFK,XLAMDH(KOFF5),NBASG,
     *                 ZERO,SCR3(KOFF6),NRHFK)
C
            KOFF7 = ILMRHF(ISYMI) + 1
            KOFF8 = IMATIJ(ISYMK,ISYMI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF7),NBASG,
     *                 ZERO,SCR4(KOFF8),NRHFK)
C
  110    CONTINUE
C
         DO 120 ISYMJ = 1,NSYM
C
            ISYMLJ = MULD2H(ISYML,ISYMJ)
            ISYMKI = MULD2H(ISYMLJ,ISYMOP)
            ISYMCI = MULD2H(ISYMJ,ISYMD)
C
            KSCR5 = 1
            KEND1 = KSCR5 + NMATIJ(ISYMKI)
C
            IF (ISYMKI .GT. ISYMLJ) GOTO 120
C
            DO 130 J = 1,NRHF(ISYMJ)
C
               DO 140 ISYMI = 1,NSYM
C
                  ISYMC = MULD2H(ISYMI,ISYMCI)
                  ISYMK = MULD2H(ISYMI,ISYMKI)
C
                  NVIRC = MAX(NVIR(ISYMC),1)
                  NRHFK = MAX(NRHF(ISYMK),1)
C
                  KOFF2 = IT1AMT(ISYMK,ISYMC) + 1
                  KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *                  + NT1AM(ISYMCI)*(J - 1)
     *                  + IT1AM(ISYMC,ISYMI) + 1
                  KOFF4 = KSCR5 + IMATIJ(ISYMK,ISYMI)
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NVIR(ISYMC),ONE,SCR3(KOFF2),NRHFK,
     *                       SCRM(KOFF3),NVIRC,ZERO,WORK(KOFF4),NRHFK)
C
  140          CONTINUE
C
               IF (ISYMJ .EQ. ISYMD) THEN
                  CALL DAXPY(NMATIJ(ISYMKI),XLAM(J),SCR4,1,
     *                       WORK(KSCR5),1)
               ENDIF
C
               NLJ = IMATIJ(ISYML,ISYMJ) + NRHF(ISYML)*(J - 1) + L
C
               IF (ISYMOP .EQ. 1) THEN
                  KKILJ = IGAMMA(ISYMKI,ISYMLJ) + NLJ*(NLJ-1)/2
                  DO 150 NKI = 1,NLJ
C
                     KOFF = KSCR5 + NKI - 1
                     NKILJ = KKILJ + NKI
                     GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  150             CONTINUE
               ELSE
                  KKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                  + NMATIJ(ISYMKI)*(NLJ - 1)
                  DO 160 NKI = 1,NMATIJ(ISYMKI)
C
                     KOFF = KSCR5 + NKI - 1
                     NKILJ = KKILJ + NKI
                     GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  160             CONTINUE
               END IF
C
  130       CONTINUE
  120    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C
C
C
C  /* Deck ccrhs_c */
      SUBROUTINE CPRHS_C(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                   XLAMDP,XLAMIP,XLAMDH,
     *                   XLAMPC,ISYMPC,XLAMHC,ISYMHC,SCRM,E1PIM,
     *                   WORK,LWORK,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                   IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_C written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 27-July-1994
C     Generalisation for CCLR by Ove Christiansen august-september 1995
C     (right transformation) and september 1996 (F-matrix).
C
C     Extended for CCSDR12, C. Neiss spring 2006
C     IOPTR12 = 1 Calculate both conv. C and r12 C' intermediates;
C                 T2-dependent contr. to C' interm. is added with a prefactor
C                 of 2*FACTC
C     IOPTE   = 1 Calculate the T-dependent part of the 
C                 E_{a delta')^1' intermediate (on E1PIM).
C
C     Purpose: Calculate C-term.
C     Adaptation: Remove the Sum(dl)T(al,di)*I(lckd) contribution
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION XINT(*),DSRHF(*),OMEGA2(*),XLAMDH(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMIP(*),SCRM(*),XLAMPC(*),XLAMHC(*)
      DIMENSION T2AM(*),E1PIM(*)
      CHARACTER CFIL*(*),CPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C--------------------------------------
C     Dynamic allocation of work space.
C--------------------------------------
C
      KSCR1 = 1
      KSCR2 = KSCR1 + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3 = KSCR2 + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      IF (IOPTR12.EQ.1) THEN
         KSCR4  = KEND1
         KEND1  = KSCR4  + NT2BCD(ISYAIK)
      END IF
      
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CCRHS_C')
      END IF
C
C--------------------------------------
C     Transpose the cluster amplitudes.
C--------------------------------------
C
      IF (ICON .GE. 2) THEN
         IF (.NOT. T2TCOR) THEN
            CALL CCSD_T2TP(T2AM,WORK(KEND1),LWRK1,ISYMT2)
         ENDIF
         IF (.NOT. DUMPCD) CALL CCSD_T2MTP(SCRM,WORK(KEND1),LWRK1,ISYMD)
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      IF (.NOT. CC2) THEN
         CALL CPRHS_C1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,SCRM,E1PIM,
     *                 WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),WORK(KSCR4),
     *                 XLAMDP,XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                 WORK(KEND1),LWRK1,
     *                 ISYDIS,IDEL,ISYMD,FACTC,ICON,IOPTR12,IOPTE,
     *                 LUC,CFIL,LUCP,CPFIL,IV)
      ENDIF
C
C--------------------------------------
C     Transpose the cluster amplitudes.
C--------------------------------------
C
      IF (ICON .GE. 2) THEN
         IF (.NOT. T2TCOR) THEN
            CALL CCSD_T2TP(T2AM,WORK(KEND1),LWRK1,ISYMT2)
         ENDIF
         IF (.NOT. DUMPCD) CALL CCSD_T2MTP(SCRM,WORK(KEND1),LWRK1,ISYMD)
      ENDIF
C
      RETURN
      END
      SUBROUTINE CPRHS_C1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,SCRM,E1PIM,
     *                    SCR1,SCR2,SCR3,SCR4,XLAMDP,XLAMIP,
     *                    XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,WORK,
     *                    LWORK,ISYDIS,IDEL,ISYDEL,FACTC,ICON,IOPTR12,
     *                    IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_C1 written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 27-July-1994
C
C     modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTC) ( assumes DUMCD )
C     and - calculate intermediates for CCLR.
C
C     modification by Ove Christiansen 17-9-1996 for calculating
C     local C-intermediate for F-matrix transformation.
C
C     Thus:
C
C     Modification to calculate terms in 2C1 right transformation in CCLR:
C
C                         1. if icon = 2 both contributions are calculated,
C                            otherwise if ICON = 
C                            1:only the integral (ki | ac)
C                              = (k i-bar | a c) + (k i | a-bar c)
C
C                         3: (k i-bar | a c) + (k i | a-bar c)
C                              + FACTC*Sum(xT*int)
C                                where xT may be non total symmetric.
C
C                         2. Allow for general transformation matrix for
C                            alpha to a(XLAMPC) and for i (XLAMHC).
C                            (the extra i transformation introduces new
C                             blocks which is only entered when 
C                             icon =1 or 3)
C
C                         3. If icon diff. from 2 (we have linear response)
C                            The C intermediate is stored according to
C                            the number of simultaneous trial vector
C                            given by IV. This is ensured using IT2DLR.
C
C     Thus in energy calc: icon = 2,fact = 1/2
C     For right transformation: 
C         icon=1,fact=anything, iv = current vector being transformed
C     For F-matrix transformation:
C         icon=3,fact=1.0, NB - not implemented several vectors yet.
C
C     extended for CCSDR12, C. Neiss spring 2006
C
C     Purpose: Calculate C-term intermediate.
C     Adaptation: Remove the Sum(dl)T(al,di)*I(lckd) contribution
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER (TWO=2.0D0)
      DIMENSION XINT(*),OMEGA2(*),T2AM(*),DSRHF(*)
      DIMENSION SCRM(*),E1PIM(*)
      DIMENSION SCR1(*), SCR2(*), SCR3(*), SCR4(*)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*),XLAMPC(*),XLAMHC(ISYMHC)
      DIMENSION WORK(LWORK)
      INTEGER   NADP(8),IADP(8,8),IBASX(8)
      CHARACTER CFIL*(*),CPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "r12int.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (IOPTE.EQ.1) THEN
        IF (.NOT.CCR12) CALL QUIT('IOPTE only implemented for CC-R12')
        IBASX(1) = 0
        DO ISYM = 2, NSYM
          IBASX(ISYM) = IBASX(ISYM-1) + MBAS2(ISYM-1)
        END DO
        DO ISYM = 1, NSYM
          NADP(ISYM) = 0
          DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IADP(ISYM1,ISYM2) = NADP(ISYM)
            NADP(ISYM) = NADP(ISYM) + 
     &                   NVIR(ISYM1)*(MBAS1(ISYM2)+MBAS2(ISYM2))
          END DO
        END DO
      END IF  
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      ISAIK2 = MULD2H(ISYDIS,ISYMT2)
      IF (ISYAIK .NE. ISAIK2) CALL QUIT('Symmetry mismatch in CCRHS_C')
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF (ICON .GE. 2) THEN
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_C1')
         ENDIF
C
         DO 200 ISYMK  = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            ISYMDL = MULD2H(ISYDIS,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1 = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2 = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3 = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ZERO,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
C
         !save a copy for first contribution:
         IF (IOPTR12.EQ.1) THEN
           CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR4,1)
         END IF
C
      ENDIF
C
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_C1')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_C1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
      IF (DUMPCD) GOTO 800
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
  800 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C-----------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      IF ( ICON .EQ. 1 ) CALL DZERO(SCR1,NT2BCD(ISYAIK))
C
      DO 810 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 820 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
            ISYMAL= MULD2H(ISYMPC,ISYMA)
            ISYALI= MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 830 I = 1,NRHF(ISYMI)
C
               NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MALI = IT1AO(ISYMAL,ISYMI) + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI- 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTC,SCR1(KOFF3),NRHFK)
C
               IF (IOPTE.EQ.1) THEN
                 IF (ISYMI.EQ.ISYMK) THEN
                   KOFF3 = IT2BCT(ISYMK,ISYMAI) + 
     &                     NRHF(ISYMK)*(NAI - 1) + I
                   IF (IDEL.GT.NBAST) THEN
                     D = IDEL-IBASX(ISYDEL)-NBAST+MBAS1(ISYDEL)
                   ELSE
                     D = IDEL-IBAS(ISYDEL)
                   END IF
C                  WRITE(LUPRI,*)'ISYDEL, IDEL, D:',ISYDEL, IDEL, D
                   KOFFE = IADP(ISYMA,ISYDEL) + 
     &                     NVIR(ISYMA)*(D-1) + 1
                   CALL DAXPY(NVIR(ISYMA),1.5D0,SCR1(KOFF3),NRHF(ISYMK),
     &                        E1PIM(KOFFE),1)
                 END IF
               END IF
C
  830       CONTINUE
  820    CONTINUE
  810 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     I is C1 transformed.
C-----------------------------------------------
C
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 850 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 860 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 870 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  870          CONTINUE
  860       CONTINUE
  850    CONTINUE
C
      ENDIF
C---------------------------------------------------------
C     Dump to disk the new contribution.
C     energy calc icon = 2
C     rsp calc. write to position given by it2dlr(idel,iv)
C---------------------------------------------------------
C
      IF ( ICON .EQ. 2 ) THEN
C
         IOFF = IT2DEL(IDEL) + 1
C
      ELSE
C
         IOFF = IT2DLR(IDEL,IV) + 1
C
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUC,CFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
      IF (IOPTR12.EQ.1) THEN
        CALL DAXPY(NT2BCD(ISYAIK),FACTC,SCR4,1,SCR1,1)
        IF (NT2BCD(ISYAIK) .GT. 0) THEN
          CALL PUTWA2(LUCP,CPFIL,SCR1,IOFF,NT2BCD(ISYAIK))
        END IF
      END IF
C
  999 CONTINUE
C
      RETURN
      END
C
C  /* Deck ccrhs_cio */
      SUBROUTINE CPRHS_CIO(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                     ISYVEC,ISYCIM,LUC,CFIL,IV,IOPT)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_CIO written by:
C
C     asm 17-aug-1994
C
C     Ove Christiansen 30-7-1995: modified to account for general
C                                 non. total symmetric vectors (ISYVEC) and
C                                 intermediates (ISYCIM). LUC and CFIL is
C                                 used to control from which file the
C                                 intermediate is obtained.
C
C                                 if iopt = 1 the C intermediate is assumed
C                                    to be as in energy calc.
C
C                                 if iopt ne. 1 we use the intermediate
C                                    on luc with address given according to
C                                    transformed vector nr iv (iv is not 1
C                                    if several vectors are transformed
C                                    simultaneously.)
C
C                                 in energy calc: iv=1,iopt=1
C
C     PURPOSE:
C             Calculate the C-term making I/O
C     Adaptaion: Call CP_PUTC instead of CC_PUTC
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*)
      DIMENSION WORK(LWORK)
      CHARACTER CFIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*)
     &        'I/O in C-term not implemented for square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_CIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYCIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYCIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for '//
     &           'allocation in CCRHS_CIO-1')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR2 + NDISAI*NT1AO(ISYMDK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for '//
     &              'allocation in CCRHS_CIO-2')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUC,CFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C-----------------------------------------
C              Transform delta index to c.
C-----------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 190 I = 1,NRHF(ISYMI)
C
                  DO 200 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     IF ((NAI .LT. IFSTAI) .OR. (NAI .GT. ILSTAI))
     *                  GOTO 200
C
                     DO 210 ISYMJ = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
C
                        DO 220 J = 1,NRHF(ISYMJ)
C
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J-1) + A
C
                           CALL CP_PUTC(WORK(KSCR2),OMEGA2,ISYMAI,
     *                                  ISYMAJ,ISYMBI,ISYMBJ,ISYMB,
     *                                  ISYMI,ISYMJ,NAI,NAJ,I,J,
     *                                  IFSTAI)
C
  220                   CONTINUE
  210                CONTINUE
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
C
C
C  /* Deck cp_putc */
      SUBROUTINE CP_PUTC(SCR2,OMEGA2,ISYMAI,ISYMAJ,ISYMBI,ISYMBJ,
     *                   ISYMB,ISYMI,ISYMJ,NAI,NAJ,I,J,IFSTAI)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CC_PUTC written by:
C     Ove Christiansen 30-10-1995: Put in C contribution in omega vector
C                                  avoid troble on cray with optimization.
C
C     Adaptation: Remove first contribution to the C-term
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SCR2(*),OMEGA2(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (ISYMAJ .EQ. ISYMBI) THEN
C
         DO 400 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + INDEX(NAJ,NBI)
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
  400    CONTINUE
C
      ENDIF
C
      IF (ISYMAJ .LT. ISYMBI) THEN
C
         DO 500 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + NT1AM(ISYMAJ)*(NBI - 1)
     *            + NAJ
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
  500    CONTINUE
C
      ENDIF
C
      IF (ISYMBI .LT. ISYMAJ) THEN
C
         DO 600 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + NT1AM(ISYMBI)*(NAJ - 1)
     *            + NBI
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
  600    CONTINUE
C
      ENDIF
      END
C
C
C
C
C
C /* Deck cprhs_d */
      SUBROUTINE CPRHS_D(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                   XLAMDP,XLAMIP,XLAMDH,
     *                   XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                   SCRM,E1PIM,WORK,LWORK,IDEL,ISYMD,FACTD,ICON,
     *                   IOPTR12,IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_D written by Henrik Koch 9-Jan-1994
C
C     Generalisation for CCLR by Ove Christiansen august-september 1995
C     (right transformation) and september 1996 (F-matrix).
C
C     adapted for CCSDR12, C. Neiss, spring 2006
C     IOPTR12 = 1 Calculate both conv. D and r12 D' intermediates
C                 T2-dependent contr. to D' interm. is added with a prefactor
C                 of 2*FACTD
C     IOPTE   = 1 Calculate the T-dependent part of the
C                 E_{a delta')^1' intermediate (on E1PIM).
C 
C     Purpose: Calculate D-term.
C     Adaptation: Remove the sum(2*t(ai,dl)-t(di,al))*L(ldkc) contribution
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),DSRHF(*),OMEGA2(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*),SCRM(*)
      DIMENSION XLAMPC(*),XLAMHC(*)
      DIMENSION T2AM(*),E1PIM(*)
      CHARACTER DFIL*(*),DPFIL
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      IF (ISYMT2 .NE. ISYMPC) CALL QUIT('Symmetry Mismatch in CCRHS_D' )
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KSCR2  = KSCR1  + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      IF (IOPTR12.EQ.1) THEN
         KSCR4  = KEND1
         KEND1  = KSCR4  + NT2BCD(ISYAIK)
      END IF

      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_D')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CPRHS_D1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *              SCRM,E1PIM,WORK(KSCR1),
     *              WORK(KSCR2),WORK(KSCR3),WORK(KSCR4),XLAMDP,XLAMIP,
     *              XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *              WORK(KEND1),LWRK1,ISYDIS,IDEL,
     *              ISYMD,FACTD,ICON,IOPTR12,IOPTE,
     *              LUD,DFIL,LUDP,DPFIL,IV)
C
      RETURN
      END
      SUBROUTINE CPRHS_D1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                    SCRM,E1PIM,SCR1,SCR2,SCR3,SCR4,
     *                    XLAMDP,XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,
     *                    ISYMHC,WORK,LWORK,ISYDIS,IDEL,ISYDEL,FACTD,
     *                    ICON,IOPTR12,IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_D written by Henrik Koch 3-Jan-1994
C
C     Modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTD). NB: Assumes DUMCD. 
C     - calculate intermediates for CCLR.
C
C     29-9-1995 (17-9-1996 F-matrix.) Ove Christiansen:  
C     
C                 1. If Icon = 2 both contributions are calculated,
C                    for total sym. case. Otherwise 
C                    a.ICON = 1 only the integral Laikc(del)
C                               = La-bar,i,k,c + La,i-bar,k,c
C                      for Jacobian right transformation
C                    b.ICON = 3 
C                          La-bar,i,k,c + La,i-bar,k,c + Tx*Int
C                      for F-matrix times vector.
C                              
C                 2. Allow for general transformation matrix for
C                    alpha to a(XLAMPC) and for i (XLAMHC).
C                    (the extra i transformation introduces new
C                     blocks which is only entered when icon = 1 or 3)
C
C                 3. If icon diff. from 2 (we have linear response)
C                    The D intermediate is stored according to
C                    the number of simultaneous trial vector
C                    given by IV. This is ensured using IT2DLR.
C
C
C     This to calculate terms in 2C1 right transformation in CCLR.
C
C     adapted for CCSDR12, C. Neiss spring 2006
C
C     Purpose: Calculate D-term.
C     Adaptation: Remove the sum(2*t(ai,dl)-t(di,al))*L(ldkc) contribution.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER(TWO=2.0D0)
      DIMENSION XINT(*),OMEGA2(*),T2AM(*),DSRHF(*)
      DIMENSION SCRM(*),E1PIM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*)
      DIMENSION XLAMPC(*),XLAMHC(*)
      DIMENSION WORK(LWORK)
      INTEGER   NADP(8),IADP(8,8),IBASX(8)
      CHARACTER DFIL*(*),DPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "r12int.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (IOPTE.EQ.1) THEN
        IF (.NOT.CCR12) CALL QUIT('IOPTE only implemented for CC-R12')
        IBASX(1) = 0
        DO ISYM = 2, NSYM
          IBASX(ISYM) = IBASX(ISYM-1) + MBAS2(ISYM-1)
        END DO
        DO ISYM = 1, NSYM
          NADP(ISYM) = 0
          DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IADP(ISYM1,ISYM2) = NADP(ISYM)
            NADP(ISYM) = NADP(ISYM) +
     &                   NVIR(ISYM1)*(MBAS1(ISYM2)+MBAS2(ISYM2))
          END DO
        END DO
      END IF
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF (ICON .GE. 2) THEN
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_D1')
         ENDIF
C
         DO 200 ISYMK = 1,NSYM
C
            ISYMDL = MULD2H(ISYMK,ISYDIS)
            ISYMAI = MULD2H(ISYAIK,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1  = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2  = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3  = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ZERO,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
C
         !save a copy of first contribution:
         IF (IOPTR12.EQ.1) THEN
           CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR4,1)
         END IF
C
      ENDIF
C
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_D1')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_D1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 1) .OR. (ICON.EQ.3)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
      CALL DSCAL(NT2BGD(ISYDIS),-ONE,SCR2,1)
C
      ISALIK = MULD2H(ISYDIS,ISYMHC)
C
      CALL DSCAL(NT2BGD(ISALIK),-ONE,SCR3,1)
C
      DO 340 ISYMK = 1,NSYM
C
         ISYALG = MULD2H(ISYMK,ISYDIS)
         ISYALI = MULD2H(ISYMHC,ISYALG)
         NT1AOM = MAX(NT1AO(ISYALG),NT1AO(ISYALI))
C
         KSCR10 = 1
         KSCR11 = KSCR10 + N2BST(ISYALG)
         KEND1  = KSCR11 + NT1AOM
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_D1')
         END IF
C
         DO 350 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYALG,ISYMK) + NNBST(ISYALG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYALG,WORK(KSCR10))
C
            ISYALI = ISYALG
            CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
C------------------------------
C           Usual contribution.
C------------------------------
C
            DO 360 ISYMI = 1,NSYM
C
               ISYMAL = MULD2H(ISYMI,ISYALI)
               ISYMG  = ISYMI
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF2),NBASAL,XLAMDH(KOFF3),NBASG,
     *                    ZERO,WORK(KOFF4),NBASAL)
C
  360       CONTINUE
C
            NRHFK = MAX(NRHF(ISYMK),1)
            KOFF5 = IT2BGT(ISYMK,ISYALI) + K
            CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,SCR2(KOFF5),
     *                 NRHFK)
C
C----------------------------------------------------
C           In 2C1 linear tronsformation extra  cont.
C----------------------------------------------------
C
            IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
               ISYALI = MULD2H(ISYALG,ISYMHC)
C
               CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
               DO 361 ISYMI = 1,NSYM
C
                  ISYMAL = MULD2H(ISYMI,ISYALI)
                  ISYMG  = MULD2H(ISYMI,ISYMHC)
C
                  NBASAL = MAX(NBAS(ISYMAL),1)
                  NBASG  = MAX(NBAS(ISYMG),1)
C
                  KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
                  KOFF3 = IGLMRH(ISYMG,ISYMI) + 1
                  KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KOFF2),NBASAL,
     *                       XLAMHC(KOFF3),NBASG,
     *                       ZERO,WORK(KOFF4),NBASAL)
C
  361          CONTINUE
C
               NRHFK = MAX(NRHF(ISYMK),1)
               KOFF5 = IT2BGT(ISYMK,ISYALI) + K
               CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,
     *                    SCR3(KOFF5),NRHFK)
C
            ENDIF
C
  350    CONTINUE
C
  340 CONTINUE
C
      IF (DUMPCD) GOTO 700
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
  700 CONTINUE
C
C--------------------------------------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     for 2C1 transformation this means lamdpc is a C1 transformed lambda.
C--------------------------------------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      DO 710 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 720 ISYMI = 1,NSYM
C
            ISYMA  = MULD2H(ISYMI,ISYMAI)
            ISYMAL = MULD2H(ISYMPC,ISYMA)
            ISYALI = MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 730 I = 1,NRHF(ISYMI)
C
               NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
               MALI  = IT1AO(ISYMAL,ISYMI)  + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI - 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTD ,SCR1(KOFF3),NRHFK)
C
               IF (IOPTE.EQ.1) THEN
                 IF (ISYMI.EQ.ISYMK) THEN
                   KOFF3 = IT2BCT(ISYMK,ISYMAI) +
     &                     NRHF(ISYMK)*(NAI - 1) + I
                   IF (IDEL.GT.NBAST) THEN
                     D = IDEL-IBASX(ISYDEL)-NBAST+MBAS1(ISYDEL)
                   ELSE
                     D = IDEL-IBAS(ISYDEL)
                   END IF
                   KOFFE = IADP(ISYMA,ISYDEL) +
     &                     NVIR(ISYMA)*(D-1) + 1
                   CALL DAXPY(NVIR(ISYMA),-0.5D0,SCR1(KOFF3),
     &                        NRHF(ISYMK),E1PIM(KOFFE),1)
                 END IF
               END IF
C
  730       CONTINUE
  720    CONTINUE
  710 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     I is C1 transformed.
C-----------------------------------------------
C
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 750 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 760 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 770 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  770          CONTINUE
  760       CONTINUE
  750    CONTINUE
C
      ENDIF
C
C---------------------------------------
C     Dump to disk the new contribution.
C---------------------------------------
C
C
      IF ( ICON .EQ. 2 ) THEN
         IOFF = IT2DEL(IDEL) + 1
      ELSE
         IOFF = IT2DLR(IDEL,IV) + 1
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUD,DFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
      IF (IOPTR12.EQ.1) THEN
        CALL DAXPY(NT2BCD(ISYAIK),FACTD,SCR4,1,SCR1,1)
        IF (NT2BCD(ISYAIK) .GT. 0) THEN
          CALL PUTWA2(LUDP,DPFIL,SCR1,IOFF,NT2BCD(ISYAIK))
        END IF
      END IF
C
  999 CONTINUE
C
      RETURN
      END
C
C
C
C /* Deck cprhs_dg */
      SUBROUTINE CPRHS_Dg(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                   XLAMDP,XLAMIP,XLAMDH,
     *                   XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                   SCRM,E1PIM,WORK,LWORK,IDEL,ISYMD,FACTD,ICON,
     *                   IOPTR12,IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_D written by Henrik Koch 9-Jan-1994
C
C     Generalisation for CCLR by Ove Christiansen august-september 1995
C     (right transformation) and september 1996 (F-matrix).
C
C     adapted for CCSDR12, C. Neiss, spring 2006
C     IOPTR12 = 1 Calculate both conv. D and r12 D' intermediates
C                 T2-dependent contr. to D' interm. is added with a prefactor
C                 of 2*FACTD
C     IOPTE   = 1 Calculate the T-dependent part of the
C                 E_{a delta')^1' intermediate (on E1PIM).
C 
C     Purpose: Calculate D-term.
C     Adaptation: Remove the sum(2*t(ai,dl)-t(di,al))*L(ldkc) contribution
C                 Calculate g(aick) instead of L(aick)
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),DSRHF(*),OMEGA2(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*),SCRM(*)
      DIMENSION XLAMPC(*),XLAMHC(*)
      DIMENSION T2AM(*),E1PIM(*)
      CHARACTER DFIL*(*),DPFIL
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      IF (ISYMT2 .NE. ISYMPC) CALL QUIT('Symmetry Mismatch in CCRHS_D' )
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KSCR2  = KSCR1  + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      IF (IOPTR12.EQ.1) THEN
         KSCR4  = KEND1
         KEND1  = KSCR4  + NT2BCD(ISYAIK)
      END IF

      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_D')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CPRHS_Dg1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *              SCRM,E1PIM,WORK(KSCR1),
     *              WORK(KSCR2),WORK(KSCR3),WORK(KSCR4),XLAMDP,XLAMIP,
     *              XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *              WORK(KEND1),LWRK1,ISYDIS,IDEL,
     *              ISYMD,FACTD,ICON,IOPTR12,IOPTE,
     *              LUD,DFIL,LUDP,DPFIL,IV)
C
      RETURN
      END
      SUBROUTINE CPRHS_Dg1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                    SCRM,E1PIM,SCR1,SCR2,SCR3,SCR4,
     *                    XLAMDP,XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,
     *                    ISYMHC,WORK,LWORK,ISYDIS,IDEL,ISYDEL,FACTD,
     *                    ICON,IOPTR12,IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
C
C     Written by Andreas Erbs Hillers-Bendtsen, Frederik Ørsted Kjeldal,
C     Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021
C
C     Adaptation of the subroutine CCRHS_D written by Henrik Koch 3-Jan-1994
C
C     Modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTD). NB: Assumes DUMCD. 
C     - calculate intermediates for CCLR.
C
C     29-9-1995 (17-9-1996 F-matrix.) Ove Christiansen:  
C     
C                 1. If Icon = 2 both contributions are calculated,
C                    for total sym. case. Otherwise 
C                    a.ICON = 1 only the integral Laikc(del)
C                               = La-bar,i,k,c + La,i-bar,k,c
C                      for Jacobian right transformation
C                    b.ICON = 3 
C                          La-bar,i,k,c + La,i-bar,k,c + Tx*Int
C                      for F-matrix times vector.
C                              
C                 2. Allow for general transformation matrix for
C                    alpha to a(XLAMPC) and for i (XLAMHC).
C                    (the extra i transformation introduces new
C                     blocks which is only entered when icon = 1 or 3)
C
C                 3. If icon diff. from 2 (we have linear response)
C                    The D intermediate is stored according to
C                    the number of simultaneous trial vector
C                    given by IV. This is ensured using IT2DLR.
C
C
C     This to calculate terms in 2C1 right transformation in CCLR.
C
C     adapted for CCSDR12, C. Neiss spring 2006
C
C     Purpose: Calculate D-term.
C     Adaptation: Remove the sum(2*t(ai,dl)-t(di,al))*L(ldkc) contribution
C                 Calculate g(aick) instead of L(aick)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER(TWO=2.0D0)
      DIMENSION XINT(*),OMEGA2(*),T2AM(*),DSRHF(*)
      DIMENSION SCRM(*),E1PIM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*)
      DIMENSION XLAMPC(*),XLAMHC(*)
      DIMENSION WORK(LWORK)
      INTEGER   NADP(8),IADP(8,8),IBASX(8)
      CHARACTER DFIL*(*),DPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "r12int.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (IOPTE.EQ.1) THEN
        IF (.NOT.CCR12) CALL QUIT('IOPTE only implemented for CC-R12')
        IBASX(1) = 0
        DO ISYM = 2, NSYM
          IBASX(ISYM) = IBASX(ISYM-1) + MBAS2(ISYM-1)
        END DO
        DO ISYM = 1, NSYM
          NADP(ISYM) = 0
          DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IADP(ISYM1,ISYM2) = NADP(ISYM)
            NADP(ISYM) = NADP(ISYM) +
     &                   NVIR(ISYM1)*(MBAS1(ISYM2)+MBAS2(ISYM2))
          END DO
        END DO
      END IF
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF (ICON .GE. 2) THEN
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_D1')
         ENDIF
C
         DO 200 ISYMK = 1,NSYM
C
            ISYMDL = MULD2H(ISYMK,ISYDIS)
            ISYMAI = MULD2H(ISYAIK,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1  = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2  = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3  = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ZERO,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
C
         !save a copy of first contribution:
         IF (IOPTR12.EQ.1) THEN
           CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR4,1)
         END IF
C
      ENDIF
C
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_D1')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_D1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 1) .OR. (ICON.EQ.3)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
      CALL DSCAL(NT2BGD(ISYDIS),-ONE,SCR2,1)
C
      ISALIK = MULD2H(ISYDIS,ISYMHC)
C
      CALL DSCAL(NT2BGD(ISALIK),-ONE,SCR3,1)
C
      DO 340 ISYMK = 1,NSYM
C
         ISYALG = MULD2H(ISYMK,ISYDIS)
         ISYALI = MULD2H(ISYMHC,ISYALG)
         NT1AOM = MAX(NT1AO(ISYALG),NT1AO(ISYALI))
C
         KSCR10 = 1
         KSCR11 = KSCR10 + N2BST(ISYALG)
         KEND1  = KSCR11 + NT1AOM
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_D1')
         END IF
C
         DO 350 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYALG,ISYMK) + NNBST(ISYALG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYALG,WORK(KSCR10))
C
            ISYALI = ISYALG
            CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
C------------------------------
C           Usual contribution.
C------------------------------
C
            DO 360 ISYMI = 1,NSYM
C
               ISYMAL = MULD2H(ISYMI,ISYALI)
               ISYMG  = ISYMI
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF2),NBASAL,XLAMDH(KOFF3),NBASG,
     *                    ZERO,WORK(KOFF4),NBASAL)
C
  360       CONTINUE
C
            NRHFK = MAX(NRHF(ISYMK),1)
            KOFF5 = IT2BGT(ISYMK,ISYALI) + K
C           g term
            CALL DCOPY(NT1AO(ISYALI),WORK(KSCR11),1,SCR2(KOFF5),
     *                 NRHFK)
            CALL DSCAL(NT1AO(ISYALI),-TWO,SCR2(KOFF5),NRHFK)
C
C----------------------------------------------------
C           In 2C1 linear tronsformation extra  cont.
C----------------------------------------------------
C
            IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
               ISYALI = MULD2H(ISYALG,ISYMHC)
C
               CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
               DO 361 ISYMI = 1,NSYM
C
                  ISYMAL = MULD2H(ISYMI,ISYALI)
                  ISYMG  = MULD2H(ISYMI,ISYMHC)
C
                  NBASAL = MAX(NBAS(ISYMAL),1)
                  NBASG  = MAX(NBAS(ISYMG),1)
C
                  KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
                  KOFF3 = IGLMRH(ISYMG,ISYMI) + 1
                  KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KOFF2),NBASAL,
     *                       XLAMHC(KOFF3),NBASG,
     *                       ZERO,WORK(KOFF4),NBASAL)
C
  361          CONTINUE
C
               NRHFK = MAX(NRHF(ISYMK),1)
               KOFF5 = IT2BGT(ISYMK,ISYALI) + K
               CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,
     *                    SCR3(KOFF5),NRHFK)
C
            ENDIF
C
  350    CONTINUE
C
  340 CONTINUE
C
      IF (DUMPCD) GOTO 700
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
  700 CONTINUE
C
C--------------------------------------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     for 2C1 transformation this means lamdpc is a C1 transformed lambda.
C--------------------------------------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      DO 710 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 720 ISYMI = 1,NSYM
C
            ISYMA  = MULD2H(ISYMI,ISYMAI)
            ISYMAL = MULD2H(ISYMPC,ISYMA)
            ISYALI = MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 730 I = 1,NRHF(ISYMI)
C
               NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
               MALI  = IT1AO(ISYMAL,ISYMI)  + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI - 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTD ,SCR1(KOFF3),NRHFK)
C
               IF (IOPTE.EQ.1) THEN
                 IF (ISYMI.EQ.ISYMK) THEN
                   KOFF3 = IT2BCT(ISYMK,ISYMAI) +
     &                     NRHF(ISYMK)*(NAI - 1) + I
                   IF (IDEL.GT.NBAST) THEN
                     D = IDEL-IBASX(ISYDEL)-NBAST+MBAS1(ISYDEL)
                   ELSE
                     D = IDEL-IBAS(ISYDEL)
                   END IF
                   KOFFE = IADP(ISYMA,ISYDEL) +
     &                     NVIR(ISYMA)*(D-1) + 1
                   CALL DAXPY(NVIR(ISYMA),-0.5D0,SCR1(KOFF3),
     &                        NRHF(ISYMK),E1PIM(KOFFE),1)
                 END IF
               END IF
C
  730       CONTINUE
  720    CONTINUE
  710 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     I is C1 transformed.
C-----------------------------------------------
C
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 750 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 760 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 770 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  770          CONTINUE
  760       CONTINUE
  750    CONTINUE
C
      ENDIF
C
C---------------------------------------
C     Dump to disk the new contribution.
C---------------------------------------
C
C
      IF ( ICON .EQ. 2 ) THEN
         IOFF = IT2DEL(IDEL) + 1
      ELSE
         IOFF = IT2DLR(IDEL,IV) + 1
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUD,DFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
      IF (IOPTR12.EQ.1) THEN
        CALL DAXPY(NT2BCD(ISYAIK),FACTD,SCR4,1,SCR1,1)
        IF (NT2BCD(ISYAIK) .GT. 0) THEN
          CALL PUTWA2(LUDP,DPFIL,SCR1,IOFF,NT2BCD(ISYAIK))
        END IF
      END IF
C
  999 CONTINUE
C
      RETURN
      END