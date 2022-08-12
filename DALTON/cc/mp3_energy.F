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

!00000000000000000000000000000000000000000000000000
!00000000000000000000000000000000000000000000000000
!00000000000000000000000000000000000000000000000000
!00000000000000000000000111110000000000000000000000
!00000000000000000000001111111000000000000000000000
!00000000000000000000011111111110000000000000000000
!00000000000000000000111111110110000000000000000000
!00000000000000000000111111000110000000000000000000
!00000000000000000001111110011111000000000000000000
!00000000000000000001111111111101000000000000000000
!00000000000000000001111111011101100000000000000000
!00000000000000000001111111000001100000000000000000
!00000000000000000001111111100001100000000000000000
!00000000000000000001111111110001000000000000000000
!00000000000000000000111111100000000000000000000000
!00000000000000000000011111111000000000000000000000
!00000000000000000000011111110000000000000000000000
!00000000000000000000001111100000000000000000000000
!00000000000000000000001111111000000000000000000000
!00000000000000000000001111111000000000000000000000
!00000000000000000000001111110000000000000000000000
!00000000000000000000001111110000000000000000000000
!00000000000000000000000111100000000000000000000000
!00000000000000000000010011100000000000000000000000
!00000000000000000000000001100000000000000000000000
!00000000000000000000000000000000000000000000000000
!00000000000000000000000000000000000000000000000000
!00000000000000000000000000000000000000000000000000
!00000000001110000000000000000000000000000000000000
!00000000001110000000000000000000000000000000000000
!00000000001110000000000000000000000000000000000000
!00000000001110000000000000000000000000010000000000
!00000000001110000000000000000000000000010000000000
!00000000001110000000000000000000000000010000000000
!00000000001110000000000100000000000001110000000000
!00000000001110000000000100000000000001110000000000
!00000000011110000000000100010000000011110000000000
!00000000000000000000000000000000000000000000000000

      SUBROUTINE MP3_ENERGY(WORK,LWORK,ECPSD2,ECPSD3)
!######################################################################################
!                                                                                     #
!  Subroutine for the calculation of MP3 energy corrections to a HF energy in a CP    #
!  formalism                                                                          #
!                                                                                     #
!  Written by Frederik Ørsted Kjeldal, Andreas Erbs Hillers-Bendtsen,                 #
!  Nicolai Machholdt Høyer, and Kurt V. Mikkelsen in Jan 2021                         #
!                                                                                     #
!######################################################################################
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccinftap.h"
#include "inftap.h"
#include "ccnoddy.h"
#include "ccfield.h"
C
      DIMENSION WORK(LWORK)
      DOUBLE PRECISION ECPSD2, ECPSD3
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER(XMONE = -1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      WRITE (LUPRI,'(A,/)') '  '
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                           STARTING MP3 CALCULATION     '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A,/)')
     *'*********************************************************'//
     *'**********************'
C      
C
C
      IF (CPORDER.EQ.0) THEN
            WRITE(LUPRI,'(//10x,A)')
     *      'The CPS(D-0) solution is the parent CCS solution'
            WRITE(LUPRI,'(//10x,A)') 
     *      'The CP program will not calculate any corrections'
      ENDIF
C
      IF (CPORDER.EQ.1) THEN
            WRITE(LUPRI,'(//10x,A)')
     *      'The First order correction to the energy is zero'
            WRITE(LUPRI,'(//10x,A)')
     *      'So the CPS(D-1) solution is the parent CCS solution'
            WRITE(LUPRI,'(//10x,A)') 
     *      'The CP program will not calculate any corrections'
      ENDIF
C
C-------------------------------------------------
C Initialize second order correction to the energy
C-------------------------------------------------
C
      IF (CPORDER.GE.2) THEN
            WRITE (LUPRI,'(//10x,A)')
     *      'Initializing calculation of first order amplitude'//
     *      ' corrections'
            CALL FLUSH(LUPRI)
C
C           DYNAMIC ALLOCATION OF MEMORY
C
            KFOCKD = 1
            KIAJB = KFOCKD + NORBTS
            KT1AM = KIAJB + NT2AMX
            KRHST1 = KT1AM + NT1AMX
            KRHST2 = KRHST1 + NT1AM(ISYMOP)
            KT2P1 = KRHST2  +
     *          MAX(NTAMP,NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
            KT1P2 = KT2P1 +
     *          MAX(NT2SQ(ISYMOP),(NT2AMX+NTAMR12),NT2R12(1),NTG2SQ(1))     
            KT2P2 = KT1P2 + NT1AMX
            KT2P3 = KT2P2 + 
     *          MAX(NT2SQ(ISYMOP),(NT2AMX+NTAMR12),NT2R12(1),NTG2SQ(1))
            KEND1 = KT2P3 + NT2AMX
            LWRK1 = LWORK - KEND1
C
            IF (LWRK1.LT.0) THEN
              WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
              CALL QUIT('Insufficient space in CPSD')
            END IF  
C
C           GET CANONICAL ORBITAL ENERGIES
C
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            REWIND LUSIFC
C
            CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
            READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *                  LSYM,MS2
            READ (LUSIFC) NISHT
            ESCF = EMCSCF
C
            CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
            CALL GPCLOSE(LUSIFC,'KEEP')
            IF (FREEZE) THEN
                  CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
            END IF
C
            CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C           GET ( IA | JB ) INTEGRALS
C
            REWIND(LUIAJB)
            READ(LUIAJB) (WORK(KIAJB+ I-1), I = 1,NT2AMX)
C
C           Calculate 1. order corection to T2
C           Start loop over IAJB
C
            ECPSD2 = 0.0
            ISYOPE = 1
C
            TIMEE2 = SECOND()
C
            DO 100 ISYMIJ = 1,NSYM
              ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
              DO 110 ISYMJ = 1,NSYM
                ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
                DO 120 ISYMB = 1,NSYM
                  ISYMA = MULD2H(ISYMB,ISYMAB)
C
                  ISYMAI = MULD2H(ISYMA,ISYMI)
                  ISYMBJ = MULD2H(ISYMB,ISYMJ)
                  ISYMBI = MULD2H(ISYMB,ISYMI)
                  ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
                  !Loop over occupied
                  DO I=1,NRHF(ISYMI)
                    INDEXI = IRHF(ISYMI) + I
C
                    !Loop over occupied
                    DO J=1,NRHF(ISYMJ)
                      INDEXJ = IRHF(ISYMJ) + J
C
                      !Loop over virtual
                      DO A=1,NVIR(ISYMA)
                        INDEXA = IVIR(ISYMA) + A
C
                        !Loop over virtual
                        DO B=1,NVIR(ISYMB)
                          INDEXB = IVIR(ISYMB) + B
C
                          ! Indexing for epsilon
                          NBJ = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J - 1) + B
                          NAI = IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I - 1) + A
                          NBI = IT1AM(ISYMB,ISYMI)
     *                        + NVIR(ISYMB)*(I - 1) + B
                          NAJ = IT1AM(ISYMA,ISYMJ)
     *                        + NVIR(ISYMA)*(J - 1) + A
C
                          IF (ISYMAI.EQ.ISYMBJ) THEN
                            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + INDEX(NAI,NBJ)
                          ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1)+NAI
                          ELSE IF (ISYMBJ.LT.ISYMAI) THEN
                            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                          END IF
C
                          IF (ISYMAJ.EQ.ISYMBI) THEN
                            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                            + INDEX(NAJ,NBI)
                          ELSE IF (ISYMAJ.LT.ISYMBI) THEN
                            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                            + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                          ELSE IF (ISYMBI.LT.ISYMAJ) THEN
                            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                            + NT1AM(ISYMBI)*(NAJ-1)+NBI
                          END IF
C
C                         Calc. can. orb. diff.
                          epsi_iajb =(WORK(KFOCKD+INDEXA-1)
     *                               +WORK(KFOCKD+INDEXB-1)
     *                               -WORK(KFOCKD+INDEXI-1)
     *                               -WORK(KFOCKD+INDEXJ-1))
C
                          !Save first order aux amp. corr.
                          WORK(KT2P1+NAIBJ-1) = -WORK(KIAJB+NAIBJ-1)
     *                                          /(epsi_iajb)
C
                          !Add energy correction
                          ECPSD2 = ECPSD2 + (TWO*WORK(KIAJB+NAIBJ-1)
     *                                    - WORK(KIAJB+NAJBI-1))
     *                                    * WORK(KT2P1+NAIBJ-1)
                        END DO
                      END DO
                    END DO
                  END DO
  120           CONTINUE
  110         CONTINUE                
  100       CONTINUE
C
            WRITE(LUPRI,'(//10x,A)') 
     *      'A wonder has occurred: We got the first order'//
     *      ' amplitude corrections!'
C
            WRITE(LUPRI,'(//10x,A, F14.10)')
     *      'The second order correction to the energy is: ', ECPSD2
C
            TIMEE2  = SECOND() - TIMEE2
C
            CALL FLUSH(LUPRI)
      ENDIF
C
C------------------------------------------------
C Initialize third order correction to the energy
C------------------------------------------------
C
      IF (CPORDER.GE.3) THEN
          WRITE (LUPRI,'(//10x,A)')
     *      'Initializing calculation of second order amplitude'//
     *      ' corrections'
          CALL FLUSH(LUPRI)
C
C         ZERO FIRST ORDER AMPLTIUDES TO GET UNTRANSFORMED INTEGRALS
C
          CALL DZERO(WORK(KT1AM),NT1AMX)
C
C Calculate RHS used to obtain 2. order correction to T2
C
          CALL DZERO(WORK(KRHST1),NT1AMX)
          CALL DZERO(WORK(KRHST2),MAX(NTAMP,NT2AO(ISYMOP),
     *               2*NT2ORT(ISYMOP)))
C
          TIMERHS1  = SECOND()
C
          CALL CPRHSN(WORK(KRHST1),WORK(KRHST2),
     *     WORK(KT1AM),WORK(KT2P1),WORK(KEND1),LWRK1)
C
          TIMERHS1  = SECOND() - TIMERHS1
C
C  Calculate 2. order correction to T2
C      - If fourth order energy correction also Construct
C        Jacobian used to obtain 2. order correction to T1
C  Start loop over IAJB
C
C
          TIMEE3 = SECOND()
C
          ECPSD3 = 0.0
          DO 200 ISYMIJ = 1,NSYM
            ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
            DO 210 ISYMJ = 1,NSYM
              ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
              DO 220 ISYMB = 1,NSYM
                ISYMA = MULD2H(ISYMB,ISYMAB)
C
                ISYMAI = MULD2H(ISYMA,ISYMI)
                ISYMBJ = MULD2H(ISYMB,ISYMJ)
                ISYMBI = MULD2H(ISYMB,ISYMI)
                ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
                !Loop over occupied
                DO I=1,NRHF(ISYMI)
                  INDEXI = IRHF(ISYMI) + I
C
                  !Loop over occupied
                  DO A=1,NVIR(ISYMA)
                    INDEXA = IVIR(ISYMA) + A   
C
                    !Loop over virtual
                    DO J=1,NRHF(ISYMJ)
                      INDEXJ = IRHF(ISYMJ) + J
C
                      !Loop over virtual
                      DO B=1,NVIR(ISYMB)
                        INDEXB = IVIR(ISYMB) + B
C
                        NAI = IT1AM(ISYMA,ISYMI)
     *                              + NVIR(ISYMA)*(I-1) + A
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                                    + NVIR(ISYMB)*(J-1) + B
                        NBI = IT1AM(ISYMB,ISYMI)
     *                                    + NVIR(ISYMB)*(I-1) + B
                        NAJ = IT1AM(ISYMA,ISYMJ)
     *                                    + NVIR(ISYMA)*(J-1) + A
C
                        IF (ISYMAI.EQ.ISYMBJ) THEN
                          NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                          + INDEX(NAI,NBJ)
                        ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                          NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                          + NT1AM(ISYMAI)*(NBJ-1)+NAI
                        ELSE IF (ISYMBJ.LT.ISYMAI) THEN
                          NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                          + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                        END IF
C
                        IF (ISYMAJ.EQ.ISYMBI) THEN
                          NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                          + INDEX(NAJ,NBI)
                        ELSE IF (ISYMAJ.LT.ISYMBI) THEN
                          NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                          + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                        ELSE IF (ISYMBI.LT.ISYMAJ) THEN
                          NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                          + NT1AM(ISYMBI)*(NAJ-1)+NBI
                        END IF
C
C                       Calc. can. orb. diff.
                        epsi_iajb = (WORK(KFOCKD+INDEXA-1)
     *                              +WORK(KFOCKD+INDEXB-1)
     *                              -WORK(KFOCKD+INDEXI-1)
     *                              -WORK(KFOCKD+INDEXJ-1))
C
C                       Construct 2. order correction to T2
C
                        WORK(KT2P2+NAIBJ-1)= - (WORK(KRHST2+NAIBJ-1))
     *                                       /epsi_iajb
C
C                       Calc. 3. order energy correction
C
                        ECPSD3 = ECPSD3 + (TWO*WORK(KIAJB+NAIBJ-1)
     *                                  - WORK(KIAJB+NAJBI-1))
     *                                  * WORK(KT2P2+NAIBJ-1)
                      END DO
                    END DO
                  END DO
                END DO
  220         CONTINUE
  210       CONTINUE                
  200     CONTINUE
C
          WRITE(LUPRI,'(//10x,A)')
     *'A wonder has occurred: We got the second order'//
     *' amplitude corrections!'
C
          WRITE(LUPRI,'(//10x,A,F14.10)')
     *'The third order correction to the energy is: ', ECPSD3
          TIMEE3  = SECOND() - TIMEE3
          CALL FLUSH(LUPRI)    
      END IF
C
      IF (IPRINT .GT. 9) THEN
C
        IF (CPORDER.GE.2) THEN
          WRITE(LUPRI,9999) 'ECPS(D-2)      ', TIMEE2
        END IF
C
        IF (CPORDER.GE.3) THEN
          WRITE(LUPRI,9999) 'CPRHS CPS(D-3) ', TIMERHS1
          WRITE(LUPRI,9999) 'ECPS(D-3)      ', TIMEE3
        END IF
C
      END IF
C
 9999  FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
C
      WRITE (LUPRI,'(A,/)') '  '
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                            FINSHED MP3 CALCULATION     '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A,/)')
     *'*********************************************************'//
     *'**********************'
      END SUBROUTINE
