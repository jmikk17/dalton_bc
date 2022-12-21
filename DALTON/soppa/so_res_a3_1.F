C  /* Deck so_res_a3_1 */
      SUBROUTINE SO_RES_A3_1(RES1,LRES1,BDENSIJ,LBDENSIJ,BDENSAB,
     &                      LBDENSAB,DSRHF,CMO,LCMO,IDEL,ISDEL,
     &                      ISYDIS,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Kasper F. Schaltz November 2022
C
C     PURPOSE: Calculate all terms belonging to the A'^(3)(1) matrix, 
C     which requires two backtransforms on the first and second index
C     of the two-particle integrals
C     This sub-routine combined with so_res_a3_2.F calculates all the 
C     terms in the A'^(3)(1)
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (QUARTER = 0.25D0)
C
      DIMENSION RES1(LRES1), DSRHF(*)
      DIMENSION BDENSIJ(LBDENSIJ), BDENSAB(LBDENSAB)
      DIMENSION CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_A3_1')
C
C     Term: del_ij sum_dl (d b | l l ) rho_da  (line 3 first term)
C
      ISYMD = ISDEL
C      
      LCDB  = NVIR(ISYMD)
      IF (LCDB.LE.0) GOTO 999
C
      DO 070 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMTR)
         ISALBE = MULD2H(ISYMI,ISYDIS)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.15',LWORK2)
         IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_1.15',' ',KEND2,LWORK)
C
         DO 071 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) * (I - 1) + 1
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
C   
C-----------------------------------------------------------------
C                               ~
C              Generate ( l l | i delta ) for
C              given i and delta in KSCR2.
C-----------------------------------------------------------------
C
            DO 072 ISYML = 1,NSYM
C
               ISALFA = ISYML
               ISBETA = ISYML
C
               LSCR2  = NRHF(ISYML)*NRHF(ISYML)
               LSCR3  = NBAS(ISBETA)*NRHF(ISYML)
C
               KSCR2  = KEND2
               KSCR3  = KSCR2 + LSCR2
               KEND3  = KSCR3 + LSCR3
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_A3_1.16',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.16',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
C               NTOTC  = MAX(NVIR(ISYMC),1)
               NTOTA  = MAX(NRHF(ISYML),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYML) + 1
C
C              ( alpha, beta |..) * C(alpha, l) => ( beta, l | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYML),
     &                    NBAS(ISALFA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTBE)
C              ( beta, l | .. ) * C(beta, l) => ( l, l | ..)
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYML),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR3),NTOTBE,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTA)
C
C
               KOFF5  = IT1AM(ISYMA,ISYML) + 1
               KOFFBDENS1 = ILMVIR(ISYMA) - ILMVIR(1) + 1
C
               COUNT = 0
               DO A = 1, NVIR(ISYMA)
                  DO L = 1,NRHF(ISYML)
                    RES1(KOFF5+COUNT+A-1) = RES1(KOFF5+COUNT+A-1) +
     &              WORK(KSCR2 + NRHF(ISYML)*(L-1) + L) *
     &              BDENSAB(KOFFBDENS1 + A)
                  END DO
                  COUNT = COUNT + NRHF(ISYML)
               END DO
C 
C
  072       CONTINUE
C
  071    CONTINUE
C
  070 CONTINUE
C

C     Term: 1/4 del_kj del_ca sum_kc (k i | c b ) rho_ca 
C      (line 11 3rd term)
C
      ISYMB = ISDEL
C
      LCDB  = NVIR(ISYMB)
      IF (LCDB.LE.0) GOTO 999
C
      DO 380 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMKI = ISALBE
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.77',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.77',' ',KEND2,LWORK)
C
         DO 381 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 382 ISYMK = 1,NSYM
C
               ISYMI  = MULD2H(ISYMK,ISYMKI)
               ISYMA  = MULD2H(ISYMI,ISYMTR)

C-----------------------------------------------------------------
C                               ~
C              Generate ( j i | j delta) for
C              given j and delta in KSCR2.
C-----------------------------------------------------------------
C
               ISALFA = ISYMI
               ISBETA = ISYMK
C
               LSCR2  = NRHF(ISYMK)*NRHF(ISYMI)
               LSCR3  = NBAS(ISBETA)*NRHF(ISYMI)
C
               KSCR2  = KEND2
               KSCR3  = KSCR2 + LSCR2
               KEND3  = KSCR3 + LSCR3
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_A3_1.78',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.78',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYMI) + 1
               KOFF4  = ILMRHF(ISYMK) + 1
C
C              ( alpha, beta | .. ) * C(alpha,i) => ( beta, i | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYMI),
     &                    NBAS(ISALFA),QUARTER,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,ZERO,WORK(KSCR3),NTOTBE)
C                                                          ~
C              C( beta, j) * ( beta, i | ..)  = > ( j, i | j delta)
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,CMO(KOFF4),NTOTBE,
     &                    WORK(KSCR3),NTOTBE,ZERO,WORK(KSCR2),NTOTK)
C
               KOFFBDENS1 = ILMVIR(ISYMA) - ILMVIR(1) + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF5  = IT1AM(ISYMA,ISYMI) + 1
C                     ~
C              (j,i | j,delta) * rho_{a delta} => U(a,i)
               COUNT = 0
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA) 
                    RES1(KOFF5+COUNT+A-1) = RES1(KOFF5+COUNT+A-1) +
     &              WORK(KSCR2 + NRHF(ISYMK)*(I-1) + J) *
     &              BDENSAB(KOFFBDENS1 + A)
                  END DO
                  COUNT = COUNT + NVIR(ISYMA)
               END DO               

  382       CONTINUE
C
  381    CONTINUE
C
  380 CONTINUE
C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_1')
C
      RETURN
      END
