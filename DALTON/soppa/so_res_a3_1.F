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
C     This sub-routine combined with so_res_a3_2.F and so_res_a3_3.F
C     calculates all the terms in the A'^(3)(1)
C     The equation for the A'^(3)(1) matrix that has been used can be
C     found in my master thesis (currently unpublished)

#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1(LRES1), DSRHF(*)
      DIMENSION BDENSIJ(LBDENSIJ), BDENSAB(LBDENSAB)
      DIMENSION CMO(LCMO)
      DIMENSION WORK(LWORK)

      INTEGER COUNT
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_A3_1')

      I_ACT_TERM = 0

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 4)) then
C
C     Term: - 1/2 del_ij sum_cd (a c | d b ) rho_cd (line 1 fourth term)
C     Implemented as (c a | b d)
C     Calculated by setting d = delta and and transforming the integrals
C     with the trial vector (before this routine is called) 
      ISYMD = ISDEL
C      
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 040 ISYMI = 1,NSYM
C
         ISYMA  = MULD2H(ISYMI,ISYMTR)
         ISALBE = MULD2H(ISYMI,ISYDIS)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.040',LWORK2)
         IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_1.040',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
         DO 041 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) * (I - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
C   
C-----------------------------------------------------------------
C                               ~
C              Generate ( a k | i delta ) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 042 ISYMC = 1,NSYM
            
               ISALFA = ISYMC
               ISBETA = ISYMA
C
               LSCR2  = NBAS(ISBETA)*NVIR(ISYMC)
               LSCR3  = NVIR(ISYMA)*NVIR(ISYMC)
C
               KSCR2  = KEND2
               KSCR3  = KSCR2 + LSCR2
               KEND3  = KSCR3 + LSCR3
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_A3_1.041',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.041',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTC  = MAX(NVIR(ISYMC),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMVIR(ISYMC) + 1
               KOFF4  = ILMVIR(ISYMA) + 1
C
C              ( beta, alpha |..) * C(alpha, c) => ( beta, c | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NVIR(ISYMC),
     &                    NBAS(ISALFA),-HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),NTOTBE)
C              ( c, beta | .. ) * C(beta, a) => ( c, a | ..)
               CALL DGEMM('T','N',NVIR(ISYMC),NVIR(ISYMA),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),NTOTBE,
     &                    CMO(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),NTOTC)
C
C                     ~
C              (c,a | i,delta) * rho_{c delta} => U(a,i)
C
               KOFF5  = IT1AM(ISYMA,ISYMI) + (I-1)*NVIR(ISYMA)
               KOFF6  = ILMVIR(ISYMC)
C
               DO IA = 1, NVIR(ISYMA)
                  DO IC = 1,NVIR(ISYMC)
                    RES1(KOFF5+IA) = RES1(KOFF5+IA) +
     &              WORK(KSCR3 - 1 + (IA-1)*NVIR(ISYMC) + IC) *
     &              BDENSAB(KOFF6 + (IDEL-1)*NVIR(ISYMC) + IC)
                  END DO
               END DO
C 
  042       CONTINUE
C
  041    CONTINUE
C
  040 CONTINUE

      END IF


      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 8)) then
C
C     Term: - 1/2 del_ij sum_kl (a k | l b ) rho_kl (line 2 fourth term)
C     Implemented as (k a | b l)
C     Calculated by setting l = delta and and transforming the integrals
C     with the trial vector (before this routine is called) 
      ISYML = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 080 ISYMI = 1,NSYM
C
         ISYMA  = MULD2H(ISYMI,ISYMTR)
         ISALBE = MULD2H(ISYMI,ISYDIS)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.080',LWORK2)
         IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_1.080',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
         DO 081 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) * (I - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
C   
C-----------------------------------------------------------------
C                               ~
C              Generate ( a k | i delta ) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 082 ISYMK = 1,NSYM
            
               ISALFA = ISYMK
               ISBETA = ISYMA
C
               LSCR2  = NBAS(ISBETA)*NRHF(ISYMK)
               LSCR3  = NVIR(ISYMA)*NRHF(ISYMK)
C
               KSCR2  = KEND2
               KSCR3  = KSCR2 + LSCR2
               KEND3  = KSCR3 + LSCR3
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_A3_1.081',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.081',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYMK) + 1
               KOFF4  = ILMVIR(ISYMA) + 1
C
C              ( beta, alpha |..) * C(alpha, k) => ( beta, k | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYMK),
     &                    NBAS(ISALFA),-HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),NTOTBE)
C              ( k, beta | .. ) * C(beta, a) => ( k, a | ..)
               CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMA),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),NTOTBE,
     &                    CMO(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),NTOTK)
C
C                     ~
C              (k,a | i,delta) * rho_{k delta} => U(a,i)
C
               KOFF5  = IT1AM(ISYMA,ISYMI) + (I-1)*NVIR(ISYMA)
               KOFF6  = ILMRHF(ISYMK)
C
               DO IA = 1, NVIR(ISYMA)
                  DO K = 1,NRHF(ISYMK)
                    RES1(KOFF5+IA) = RES1(KOFF5+IA) +
     &              WORK(KSCR3 - 1 + (IA-1)*NRHF(ISYMK) + K) *
     &              BDENSIJ(KOFF6 + (IDEL-1)*NRHF(ISYMK) + K)
                  END DO
               END DO
C 
  082       CONTINUE
C
  081    CONTINUE
C
  080 CONTINUE

      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 12)) then
C     Term: 1/2 sum_c (j i| c b) rho_ca   (line 3 fourth term)
C     Implemented as (j i| b c)
C     Calculated by setting c = delta and transforming the integrals
C     with the trial vector (before this routine is called) 
C
      ISYMC = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 120 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMI  = MULD2H(ISALBE,ISYMJ)
         ISYMA  = MULD2H(ISYMI,ISYMTR)
         ISALFA = ISYMJ
         ISBETA = ISYMI
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.120',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.120',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
         DO 121 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               ~
C              Generate ( j i | j delta) for
C              given j and delta in KSCR3.
C-----------------------------------------------------------------
C
               LSCR2  = NBAS(ISBETA)
               LSCR3  = NRHF(ISYMI)

C
               KSCR2  = KEND2
               KSCR3  = KSCR2 + LSCR2
               KEND3  = KSCR3 + LSCR3
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_A3_1.121',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.121',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTJ  = MAX(NRHF(ISYMJ),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYMJ) + 1
               KOFF4  = ILMRHF(ISYMI) + 1
C
C              C(j,alpha) * ( alpha, beta | .. ) => ( j, beta | ..)
               CALL DGEMM('T','N',1,NBAS(ISBETA),
     &                    NBAS(ISALFA),HALF,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,ZERO,WORK(KSCR2),1)
C                                                          ~
C              ( j, beta | ..) * C( beta, i)  = > ( j, i | j delta)
               CALL DGEMM('N','N',1,NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,WORK(KSCR2),1,
     &                    CMO(KOFF4),NTOTBE,ZERO,WORK(KSCR3),1)
C
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTK  = MAX(NRHF(ISYMK),1)
C
               KOFF5  = IT1AM(ISYMA,ISYMI)
               KOFF6  = ILMVIR(ISYMA) - ILMVIR(1)
C                     ~
C              (j,i | j,delta) * rho_{a delta} => U(a,i)
               COUNT = 0                             
               DO I = 1,NRHF(ISYMI)
                  DO IA = 1,NVIR(ISYMA)
                    RES1(KOFF5+COUNT+IA) = RES1(KOFF5+COUNT+IA) +
     &              WORK(KSCR3 - 1 + I) *
     &              BDENSAB(KOFF6 + (IDEL-1)*NVIR(ISYMA) + IA)
                  END DO
                  COUNT = COUNT + NVIR(ISYMA)
               END DO               
C
  121    CONTINUE
C
  120 CONTINUE

      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 13)) then

C     Term: sum_l (a i | b l ) rho_jl  (line 4 first term)
C     Implemented as (a i | b l)
C     Calculated by setting l = delta and transforming the integrals
C     with the trial vector (before this routine is called) 
C
      ISYML = ISDEL
C
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 130 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMAI = ISALBE
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.130',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.130',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
         DO 131 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               ~
C              Generate ( a i | j delta) for
C              given j and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 132 ISYMA = 1,NSYM

                  ISALFA = ISYMA
                  ISBETA = MULD2H(ISALFA, ISALBE)
                  ISYMI  = ISBETA

                  LSCR2  = NBAS(ISALFA)*NRHF(ISYMI)
                  LSCR3  = NRHF(ISYMI)*NVIR(ISYMA)
C
                  KSCR2  = KEND2
                  KSCR3  = KSCR2 + LSCR2
                  KEND3  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.131',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.131',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYMI) + 1
                  KOFF4  = ILMVIR(ISYMA) + 1
C
C              ( alpha, beta |..) * C(beta, i)   => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C              b(a,alpha) * ( alpha, i | .. ) => ( a, i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)

C                     ~
C              (a,i | j,delta) * rho_{j delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = ILMRHF(ISYMI)
C
                  COUNT  = 0
                  DO I = 1, NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+COUNT+IA) =RES1(KOFF5+COUNT+IA)
     &                        + WORK(KSCR3 - 1 + COUNT + IA) *
     &                        BDENSIJ(KOFF6-1+J+ (IDEL-1)*NRHF(ISYMJ))
                        END DO
                        COUNT = COUNT + NVIR(ISYMA)
                  END DO
C             
  132       CONTINUE
C
  131    CONTINUE
C
  130 CONTINUE

      END IF

C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_1')
C
      RETURN
      END
