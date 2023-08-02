C
C  /* Deck so_dens3 */
      SUBROUTINE SO_DENS3(DENS3IJ,LDENSIJ,DENS3AB,LDENSAB,
     &                    T2AM,LT2AM,T2AM3,LT2AM3,WORK,LWORK)
C
C     This routine is part of the atomic integral direct TOPPA program.
C
C     Juliane Fuglsbjerg, May 2022
C     
C     PURPOSE: Calculate third-order TOPPA one-electron correction to the density
C              matrices Dij and Dab from the T2 and MP3 T2 amplitudes.
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION DENS3IJ(LDENSIJ), DENS3AB(LDENSAB)
      DIMENSION T2AM(LT2AM), T2AM3(LT2AM3)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_DENS3')
C
C---------------------------------
C     Initialize density matrices.
C---------------------------------
C
      CALL DZERO(DENS3IJ,LDENSIJ)
      CALL DZERO(DENS3AB,LDENSAB)
C
      DO 100 ISYMB = 1,NSYM
C
         DO 200 ISYMK = 1,NSYM
C
            ISYMBK = MULD2H(ISYMB,ISYMK)
            ISYMAJ = ISYMBK
C
            DO 300 ISYMJ = 1,NSYM
C
               ISYMA  = MULD2H(ISYMJ,ISYMAJ)
               ISYMAK = MULD2H(ISYMA,ISYMK)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
               DO 400 B = 1,NVIR(ISYMB)
C
C--------------------------------------------------------------
C                 The necessary space in WORK is allocated for
C                 the current symmetry with a seperate space for
C                 the different amplitudes.
C
C                 The amplitudes for all K, J, A and the 
C                 given B are collected in WORK.
C                 
C                 First order coloumb-minus-exchange like amplitudes are 
C                 collected in KT2CE1 and KT2CE2             
C                 First order amplitudes in KT2AM1 and KT2AM2
C                 Second order coloumb-minus-exchange like amplitudes in 
C                 KT2CE31 and KT2CE32
C                 Second order amplitudes in KT2AM31 and KT2AM32
C--------------------------------------------------------------
C
                  LT2 = NVIR(ISYMA)*NRHF(ISYMK)*NRHF(ISYMJ)
C
                  KT2CE1 = 1
                  KT2CE2 = KT2CE1 + LT2
                  KT2AM1 = KT2CE2 + LT2
                  KT2AM2 = KT2AM1 + LT2
                  KT2CE31 = KT2AM2 + LT2
                  KT2CE32 = KT2CE31 + LT2
                  KT2AM31 = KT2CE32 + LT2
                  KT2AM32 = KT2AM31 + LT2
                  KEND1   = KT2AM32 + LT2
                  LWORK1  = LWORK   - KEND1
C
                  CALL SO_MEMMAX ('SO_DENS3.1',LWORK1)
                  IF (LWORK1 .LT. 0)
     &                CALL STOPIT('SO_DENS3.1',' ',KEND1,LWORK)
C
                  DO 500 K = 1,NRHF(ISYMK)
C
                     DO 600 J = 1,NRHF(ISYMJ)
C
                        DO 700 A = 1,NVIR(ISYMA)
C
                           NAJ   = IT1AM(ISYMA,ISYMJ)
     &                           + NVIR(ISYMA)*(J-1) + A
                           NBK   = IT1AM(ISYMB,ISYMK)
     &                           + NVIR(ISYMB)*(K-1) + B
                           NAJBK = IT2AM(ISYMAJ,ISYMBK)
     &                           + INDEX(NAJ,NBK)
C
                           NAK   = IT1AM(ISYMA,ISYMK)
     &                           + NVIR(ISYMA)*(K-1) + A
                           NBJ   = IT1AM(ISYMB,ISYMJ)
     &                           + NVIR(ISYMB)*(J-1) + B
                           NAKBJ = IT2AM(ISYMAK,ISYMBJ)
     &                           + INDEX(NAK,NBJ)
C
                           NSQKJ = NRHF(ISYMK)*(J-1) + K
                           NSQJK = NRHF(ISYMJ)*(K-1) + J
                           NAKJ  = NVIR(ISYMA)*(NSQKJ-1) + A
                           NAJK  = NVIR(ISYMA)*(NSQJK-1) + A
C
                           NSQAK = NVIR(ISYMA)*(K-1) + A
                           NSQKA = NRHF(ISYMK)*(A-1) + K
                           NJAK  = NRHF(ISYMJ)*(NSQAK-1) + J
                           NJKA  = NRHF(ISYMJ)*(NSQKA-1) + J
C
                           T2MP  = TWO * T2AM(NAJBK) - T2AM(NAKBJ)
C
                           WORK(KT2CE1 + NAKJ - 1) = T2MP
                           WORK(KT2CE2 + NAJK - 1) = T2MP
C
                           WORK(KT2AM1 + NJAK - 1) = T2AM(NAJBK)
                           WORK(KT2AM2 + NJKA - 1) = T2AM(NAJBK)
C                           
                           T2MP3 = TWO * T2AM3(NAJBK) - T2AM3(NAKBJ)
C
                           WORK(KT2CE31 + NAKJ - 1) = T2MP3
                           WORK(KT2CE32 + NAJK - 1) = T2MP3
C
                           WORK(KT2AM31 + NJAK - 1) = T2AM3(NAJBK)
                           WORK(KT2AM32 + NJKA - 1) = T2AM3(NAJBK)
C
  700                   CONTINUE
C
  600                CONTINUE
C
  500             CONTINUE
C
C---------------------------------------------------------------------
C                 Calculate Dij by adding in contributions for each B.
C---------------------------------------------------------------------
C
                  NTOTJ  = MAX(NRHF(ISYMJ),1)
                  NTOTAK = MAX(NVIR(ISYMA)*NRHF(ISYMK),1)
C
                  KOFF1  = IIJDEN(ISYMJ,ISYMJ) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMJ),NRHF(ISYMJ),
     &                       NVIR(ISYMA)*NRHF(ISYMK),-ONE,
     &                       WORK(KT2AM1),NTOTJ,WORK(KT2CE31),NTOTAK,
     &                       ONE,DENS3IJ(KOFF1),NTOTJ)
C
                  CALL DGEMM('N','N',NRHF(ISYMJ),NRHF(ISYMJ),
     &                       NVIR(ISYMA)*NRHF(ISYMK),-ONE,
     &                       WORK(KT2AM31),NTOTJ,WORK(KT2CE1),NTOTAK,
     &                       ONE,DENS3IJ(KOFF1),NTOTJ)
C
C---------------------------------------------------------------------
C                 Calculate Dab by adding in contributions for each B.
C---------------------------------------------------------------------
C
                  NTOTA  = MAX(NVIR(ISYMA),1)
                  NTOTJK = MAX(NRHF(ISYMJ)*NRHF(ISYMK),1)
C
                  KOFF2  = IABDEN(ISYMA,ISYMA) + 1
C
                  CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMA),
     &                       NRHF(ISYMJ)*NRHF(ISYMK),ONE,
     &                       WORK(KT2CE32),NTOTA,WORK(KT2AM2),NTOTJK,
     &                       ONE,DENS3AB(KOFF2),NTOTA)
C     
                  CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMA),
     &                       NRHF(ISYMJ)*NRHF(ISYMK),ONE,
     &                       WORK(KT2CE2),NTOTA,WORK(KT2AM32),NTOTJK,
     &                       ONE,DENS3AB(KOFF2),NTOTA)
C
  400          CONTINUE
C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_DENS3')
C
      RETURN
      END
