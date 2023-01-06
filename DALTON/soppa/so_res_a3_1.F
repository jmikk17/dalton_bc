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
C     The equation for the A'^(3)(1) matrix that has been used can be
C     found in my master thesis (currently unpublished)

#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (QUARTER = 0.25D0)
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
C

C     Term: 1/2 sum_c (j i| c b) rho_ca   (line 2 third term)
C     Implemented as (j i| b c)
C     Calculated by setting c = delta and transforming the integrals
C     with the trial vector (before this routine is called) 
C
      ISYMC = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 060 ISYMJ = 1,NSYM
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
         CALL SO_MEMMAX ('SO_RES_A3_1.60',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.60',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
         DO 061 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               ~
C              Generate ( j i | j delta) for
C              given j and delta in KSCR2.
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
               CALL SO_MEMMAX ('SO_RES_A3_1.61',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.61',' ',KEND3,LWORK)
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
  061    CONTINUE
C
  060 CONTINUE

C     Term: 1/2 sum_l (a i | l b ) rho_jl  (line 3 fourth term)
C     Implemented as (a i | b l)
C     Calculated by setting l = delta and transforming the integrals
C     with the trial vector (before this routine is called) 
C
      ISYML = ISDEL
C
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 080 ISYMJ = 1,NSYM
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
         CALL SO_MEMMAX ('SO_RES_A3_1.080',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.080',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
         DO 081 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               ~
C              Generate ( a i | j delta) for
C              given j and delta in KSCR2.
C-----------------------------------------------------------------
C
            DO 082 ISYMA = 1,NSYM

                  ISALFA = ISYMA
                  ISBETA = MULD2H(ISALFA, ISALBE)
                  ISYMI  = ISBETA

                  LSCR2  = NBAS(ISALFA)*NRHF(ISYMI)
                  LSCR3  = NRHF(ISYMI)*NVIR(ISYMA)
C
                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.081',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.081',' ',KEND4,LWORK)
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
     &                    NBAS(ISBETA),HALF,
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
  082       CONTINUE
C
  081    CONTINUE
C
  080 CONTINUE


C     Term: 1/2 sum_l (a j | l b ) rho_li  (line 3 fourth term)
C     Implemented as (a j | b l)
C     Calculated by setting l = delta and transforming the integrals
C     with the trial vector (before this routine is called) 
C
      ISYML = ISDEL
C
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 110 ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISBETA = ISYMJ
         ISALFA = MULD2H(ISBETA, ISALBE)
         ISYMA  = ISALFA
         ISYMI  = MULD2H(ISYMA, ISYMTR)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.110',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.110',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
         DO 111 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) * (J - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               ~
C              Generate ( a j | j delta) for
C              given j and delta in KSCR2.
C-----------------------------------------------------------------
C
                  LSCR2  = NBAS(ISALFA)
                  LSCR3  = NRHF(ISYMI)
C
                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.111',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.111',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYMJ) + (J-1)*NBAS(ISYMJ)+1
                  KOFF4  = ILMRHF(ISYMA) + 1
C
C              ( alpha, beta |..) * C(beta, j)   => ( alpha, j | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),1,
     &                    NBAS(ISBETA),HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                                                   
C              b(a,alpha) * ( alpha, j | .. ) => ( a, j | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),1,
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)

C                     ~
C              (a,j | j,delta) * rho_{i delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = ILMRHF(ISYMI)
C
                  COUNT  = 0
                  DO I = 1, NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+COUNT+IA) =RES1(KOFF5+COUNT+IA)
     &                        + WORK(KSCR3 - 1 + IA) *
     &                        BDENSIJ(KOFF6 + (IDEL-1)*NRHF(ISYMI) + I)
                        END DO
                        COUNT = COUNT + NVIR(ISYMA)
                  END DO
C             

C
  111    CONTINUE
C
  110 CONTINUE


C     Term: 3/2 del_ij sum_cl (c b | l l ) rho_ca  (line 4 first term)
C     Implemented as (l l | b c)
C     Calculated by setting c = delta and and transforming the integrals
C     with the trial vector (before this routine is called) 
      ISYMC = ISDEL
C      
      LCDB  = NVIR(ISYMC)
      IF (LCDB.LE.0) GOTO 999
C
      DO 120 ISYMI = 1,NSYM
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
         CALL SO_MEMMAX ('SO_RES_A3_1.120',LWORK2)
         IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_1.120',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
         DO 121 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) * (I - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
C   
C-----------------------------------------------------------------
C                               ~
C              Generate ( l l | i delta ) for
C              given i and delta in KSCR2.
C-----------------------------------------------------------------
C
            DO 122 ISYML = 1,NSYM
C
               ISALFA = ISYML
               ISBETA = ISYML
C
               LSCR2  = NBAS(ISBETA)*NRHF(ISYML)
               LSCR3  = NRHF(ISYML)*NRHF(ISYML)
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
               NTOTL  = MAX(NRHF(ISYML),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYML) + 1
C
C              ( beta, alpha |..) * C(alpha, l) => ( beta, l | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYML),
     &                    NBAS(ISALFA),ONE+HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),NTOTBE)
C              ( l, beta | .. ) * C(beta, l) => ( l, l | ..)
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYML),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),NTOTBE,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR3),NTOTL)
C
C                     ~
C              (l,l | i,delta) * rho_{a delta} => U(a,i)
C
               KOFF4  = ILMVIR(ISYMA) - ILMVIR(1)
               KOFF5  = IT1AM(ISYMA,ISYMI) + (I-1)*NVIR(ISYMA)
C
               DO IA = 1, NVIR(ISYMA)
                  DO L = 1,NRHF(ISYML)
                    RES1(KOFF5+IA) = RES1(KOFF5+IA) +
     &              WORK(KSCR3 - 1 + NRHF(ISYML)*(L-1) + L) *
     &              BDENSAB(KOFF4 + (IDEL-1)*NVIR(ISYMA) + IA)
                  END DO
               END DO
C 
  122       CONTINUE
C
  121    CONTINUE
C
  120 CONTINUE
C
C     Term: del_ij sum_kl (a k | l b ) rho_kl  (line 4 first term)
C     Implemented as (a k | b l)
C     Calculated by setting l = delta and and transforming the integrals
C     with the trial vector (before this routine is called) 
      ISYML = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
      DO 220 ISYMI = 1,NSYM
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
         CALL SO_MEMMAX ('SO_RES_A3_1.220',LWORK2)
         IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_1.220',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              ~
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
         DO 221 I = 1,NRHF(ISYMI)
C
            KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) * (I - 1) + 1

            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
C   
C-----------------------------------------------------------------
C                               ~
C              Generate ( a k | i delta ) for
C              given i and delta in KSCR2.
C-----------------------------------------------------------------
C
            DO 222 ISYMK = 1,NSYM
            
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
               CALL SO_MEMMAX ('SO_RES_A3_1.221',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_A3_1.221',' ',KEND3,LWORK)
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
     &                    NBAS(ISALFA),ONE,
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
                  DO K = 1,NRHF(ISYML)
                    RES1(KOFF5+IA) = RES1(KOFF5+IA) +
     &              WORK(KSCR3 - 1 + (IA-1)*NRHF(ISYMK) + K) *
     &              BDENSIJ(KOFF6 + (IDEL-1)*NRHF(ISYMK) + K)
                  END DO
               END DO
C 
  222       CONTINUE
C
  221    CONTINUE
C
  220 CONTINUE

C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_1')
C
      RETURN
      END
