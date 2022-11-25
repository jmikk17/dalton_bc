C  /* Deck so_res_a3_1 */
      SUBROUTINE SO_RES_A3_1(RES1,LRES1,BDENSIJ,LBDENSIJ,BDENSAB,
     &                      LBDENSAB,DSRHF,CMO,LCMO,IDEL,ISDEL,
     &                      ISYDIS,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Kasper F. Schaltz November 2022
C
C     PURPOSE: Calculate all terms belonging to the A'^(3)(1) matrix.
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
      ISYMB = ISDEL
C      
      LCDB  = NRHF(ISYMD)
      IF (LCDB.LE.0) GOTO 999
C
      DO 070 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMLL = ISALBE
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.1',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_RES_A3_1.1',' ',KEND2,LWORK)
C
         DO 071 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) * (I - 1) + 1
C
C-----------------------------------------------------------------------
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 072 ISYMK = 1,NSYM
C
               ISKBJ  = MULD2H(ISYMK,ISYMBJ)
               ISYML  = ISKBJ
               ISYMLK = MULD2H(ISYML,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMLJ = MULD2H(ISYML,ISYMJ)
               ISYML  = MULD2H(ISYML,ISYMLL)               
C
C-----------------------------------------------------------------
C              Generate ( l l | i delta ) for
C              given i and delta in KSCR2.
C-----------------------------------------------------------------
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
               CALL SO_MEMMAX ('SO_RES_A3_1.2',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.2',' ',KEND3,LWORK)
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
C               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYML),
C     &                    NBAS(ISALFA),ONE,
C     &                    WORK(KOFF2),NTOTAL,
C     &                    CMO(KOFF3),NTOTAL,
C     &                    ZERO,WORK(KSCR3),NTOTBE)
C              ( beta, l | .. ) * C(beta, l) => ( l, l | ..)
C               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYML),
C     &                    NBAS(ISBETA),ONE,
C     &                    WORK(KSCR3),NTOTBE,
C     &                    CMO(KOFF3),NTOTBE,
C     &                    ZERO,WORK(KSCR2),NTOTA)
C
C     This part hasn't been completely changed yet
C
               NTOTA  = MAX(NRHF(ISYML),1)
               NTOTL  = MAX(NRHF(ISYML),1)
C
               KOFF4  = ILMVIR(ISYMA) - ILMVIR(1) + 1
C                        ~                                        ~
C               ( l, l | j, delta ) * rho_{a delta} = > x ( l l , j a)
C               CALL DGEMM('N','T',NRHF(ISYML),NRHF(ISYML),NRHF(ISYML),
C     &                    -ONE,WORK(KSCR2),NTOTA,BDENSAB(KOFF4),NTOTA,
C     &                    ONE,RES1(KOFF3),NTOTC)
C                        ~
C               ( a, c | j, delta ) * T2( a,k, j, delta ) = > x ( c, k)
C               CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMA),
C     &                    -ONE,WORK(KSCR2),NTOTA,T2M1(KOFFX2M1),NTOTA,
C     &                    ONE,RES1(KOFF3),NTOTC)     
C                                     ~
C              T2(a,k,j,delta)*(k,i | j,delta) => x( a, i)
C               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
C     &                    ONE,T2M1(KOFFX2M1),NTOTA,WORK(KSCR3),
C     &                    NTOTK,ONE,RES1(KOFF4),NTOTA)
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
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMKI = ISALBE
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.3',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_RES_A3_1.3',' ',KEND2,LWORK)
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
               ISKBJ  = MULD2H(ISYMK,ISYMBJ)
               ISYMA  = ISKBJ
               ISYMI  = MULD2H(ISYMK,ISYMKI)
               ISYMAK = MULD2H(ISYMA,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)

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
               CALL SO_MEMMAX ('SO_RES_A3_1.4',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.4',' ',KEND3,LWORK)
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
C              C( beta, k) * ( beta, i | ..)  = > ( j, i | j delta)
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,CMO(KOFF4),NTOTBE,
     &                    WORK(KSCR3),NTOTBE,ZERO,WORK(KSCR2),NTOTK)
C
               KOFFBDENS1 = ILMVIR(ISYMA) - ILMVIR(1) + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF4  = IT1AM(ISYMA,ISYMI) + 1
C                     ~
C              (j,i | j,delta) * rho_{a delta} => U(a,i)
               COUNT = 0
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA) 
                    RES1(KOFF4+COUNT+A-1) = RES1(KOFF4+COUNT+A-1) +
     &              WORK(KSCR2 + NRHF(ISYMK)*(I-1)+J) *
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
