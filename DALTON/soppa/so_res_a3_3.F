C  /* Deck so_res_a3_3 */
      SUBROUTINE SO_RES_A3_3(RES1,LRES1,B2DENSIJ,LB2DENSIJ,B2DENSAB,
     &                      LB2DENSAB,XINT,CMO,LCMO,IDEL,ISDEL,ISYDIS,
     &                      ISYMTR,OBTRIAL,LOBTRIAL,VBTRIAL,LVBTRIAL,
     &                      WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Kasper F. Schaltz November 2022
C
C     PURPOSE: Calculate all terms belonging to the A'^(3)(1) matrix, 
C     which utilizes a double backtransformed density matrix.
C     This sub-routine combined with so_res_a3_1.F and so_res_a3_2.F 
C     calculates all the terms in the A'^(3)(1) matrix.
C     The equation for the A'^(3)(1) matrix that has been used can be
C     found in my master thesis (currently unpublished)

#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES1(LRES1), XINT(*)
      DIMENSION B2DENSIJ(LB2DENSIJ), B2DENSAB(LB2DENSAB)
      DIMENSION CMO(LCMO)
      DIMENSION WORK(LWORK)
      DIMENSION OBTRIAL(LOBTRIAL), VBTRIAL(LVBTRIAL)

      INTEGER COUNT, ALPHA, GAMMA
C OBTRIAL = Trialvector with occupied index not transformed
C VBTRIAL = Trialvector with virtual index not transformed
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
      CALL QENTER('SO_RES_A3_3')

      I_ACT_TERM = 0

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 1)) then
C     Term: - del_ab sum_cd (i j | c d ) rho_cd  (line 1 first term)
C     Implemented as (i j | c d)
C
      ISYMD = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
      DO 010 ISYMGA = 1,NSYM
C
         ISALBE = MULD2H(ISYMGA,ISYDIS)
         ISYMAI = ISALBE
         
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_3.10',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_3.10',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | gamma delta ) for given 
C           gamma and delta.
C-----------------------------------------------------------------------
C
         DO 011 GAMMA = 1,NBAS(ISYMGA)
C
            KOFF1 = IDSBAS(ISALBE,ISYMGA) + NNBST(ISALBE) 
     &                  * (GAMMA - 1) + 1

            CALL CCSD_SYMSQ(XINT(KOFF1),ISALBE,WORK(KSCR1))

C-----------------------------------------------------------------
C                               
C              Generate ( a i | gamma delta) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 012 ISYMI = 1,NSYM

                  ISYMA  = MULD2H(ISYMI,ISYMTR)
                  ISALFA = ISYMA
                  ISBETA = ISYMI

                  LSCR2  = NVIR(ISYMA)*NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)*NVIR(ISYMA)
C
                  KSCR2  = KEND2
                  KSCR3  = KSCR2 + LSCR2
                  KEND3  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.011',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.011',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYMI) + 1
                  KOFF4  = IMATAV(ISYMI,ISYMA) + 1
C
C              b(a,alpha) * ( alpha, beta | .. ) => ( a, beta | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NBAS(ISBETA),
     &                    NBAS(ISALFA),-HALF,
     &                    VBTRIAL(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)
C
C              ( a, beta |..) * C(beta, i)   => ( a, i | ..)
                  CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C                     
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISALFA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR2 - 1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSAB(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  012       CONTINUE
C
  011    CONTINUE
C
  010 CONTINUE

      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 2)) then
C     Term: del_ab 1/2 sum_cd (j c | d i ) rho_cd  (line 1 second term)
C     Implemented as (j c | i d)
C
      ISYMD = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | i delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND1  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_3.20',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_3.20',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND1),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over i
C-----------------------------------------------------------------------
C
      DO 020 ISYMI = 1,NSYM
C
         ISALBE = MULD2H(ISYMI,ISYDIS)
         ISYMA  = MULD2H(ISYMI,ISYMTR)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_3.21',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_3.21',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
         DO 021 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) 
     &                  * (I - 1) + KDSRHF

            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               
C              Generate ( c j | i delta) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 022 ISYMC = 1,NSYM

                  ISALFA = ISYMC
                  ISBETA = ISYMA

                  LSCR2  = NBAS(ISALFA)*NVIR(ISYMA)
C
                  KSCR2  = KEND2
                  KEND3  = KSCR2 + LSCR2
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.022',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.022',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = IMATAV(ISYMI,ISYMA) + 1
C
C              ( alpha, beta |..) * b(beta, a)  => ( alpha, a | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NVIR(ISYMA),
     &                    NBAS(ISBETA),HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    VBTRIAL(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C                     
C              (alpha, a | i,delta) * rho_{alpha delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISALFA,ISDEL)
C
                  DO IA = 1,NVIR(ISYMA)
                        DO ALPHA = 1,NBAS(ISALFA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR2 -(IA-1)*NBAS(ISALFA)+ALPHA)*
     &                      B2DENSAB(KOFF6 +(IDEL-1)*NBAS(ISALFA)+ALPHA)
                        END DO
                  END DO
C             
  022       CONTINUE
C
  021    CONTINUE
C
  020 CONTINUE

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 3)) then
C     Term: - del_ab sum_cd (a b | c d ) rho_cd  (line 1 third term)
C     Implemented as (a b | c d)
C
      ISYMD = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
      DO 030 ISYMGA = 1,NSYM
C
         ISALBE = MULD2H(ISYMGA,ISYDIS)
         ISYMAI = ISALBE
         
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_3.30',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_3.30',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | gamma delta ) for given 
C           gamma and delta.
C-----------------------------------------------------------------------
C
         DO 031 GAMMA = 1,NBAS(ISYMGA)
C
            KOFF1 = IDSBAS(ISALBE,ISYMGA) + NNBST(ISALBE) 
     &                  * (GAMMA - 1) + 1

            CALL CCSD_SYMSQ(XINT(KOFF1),ISALBE,WORK(KSCR1))

C-----------------------------------------------------------------
C                               
C              Generate ( a i | gamma delta) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 032 ISYMI = 1,NSYM

                  ISYMA  = MULD2H(ISYMI,ISYMTR)
                  ISALFA = ISYMI
                  ISBETA = ISYMA

                  LSCR2  = NVIR(ISYMA)*NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)*NVIR(ISYMA)
C
                  KSCR2  = KEND2
                  KSCR3  = KSCR2 + LSCR2
                  KEND3  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.031',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.031',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYMI) + 1
                  KOFF4  = IT1AO(ISYMI,ISYMA) + 1
C
C               ( alpha beta | .. ) * b(beta,i) => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),NTOTAL,
     &                    OBTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),NTOTAL)
C
C              C(a,alpha) * ( alpha, i |..)  => ( a, i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF3),NTOTAL,
     &                    WORK(KOFF2),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTA)
C
C                     
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISALFA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR2 - 1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSAB(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  032       CONTINUE
C
  031    CONTINUE
C
  030 CONTINUE

      END IF


      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 5)) then
C     Term: - del_ab sum_kl (i j | k l ) rho_kl  (line 2 first term)
C     Implemented as (i j | k l)
C
      ISYML = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
      DO 050 ISYMGA = 1,NSYM
C
         ISALBE = MULD2H(ISYMGA,ISYDIS)
         ISYMAI = ISALBE
         
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_3.50',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_3.50',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | gamma delta ) for given 
C           gamma and delta.
C-----------------------------------------------------------------------
C
         DO 051 GAMMA = 1,NBAS(ISYMGA)
C
            KOFF1 = IDSBAS(ISALBE,ISYMGA) + NNBST(ISALBE) 
     &                  * (GAMMA - 1) + 1

            CALL CCSD_SYMSQ(XINT(KOFF1),ISALBE,WORK(KSCR1))

C-----------------------------------------------------------------
C                               
C              Generate ( a i | gamma delta) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 052 ISYMI = 1,NSYM

                  ISYMA  = MULD2H(ISYMI,ISYMTR)
                  ISALFA = ISYMA
                  ISBETA = ISYMI

                  LSCR2  = NVIR(ISYMA)*NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)*NVIR(ISYMA)
C
                  KSCR2  = KEND2
                  KSCR3  = KSCR2 + LSCR2
                  KEND3  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.051',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.051',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYMI) + 1
                  KOFF4  = IMATAV(ISYMI,ISYMA) + 1
C
C              b(a,alpha) * ( alpha, beta | .. ) => ( a, beta | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NBAS(ISBETA),
     &                    NBAS(ISALFA),-HALF,
     &                    VBTRIAL(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)
C
C              ( a, beta |..) * C(beta, i)   => ( a, i | ..)
                  CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C                     
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISALFA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR2 - 1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSIJ(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  052       CONTINUE
C
  051    CONTINUE
C
  050 CONTINUE

      END IF


      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 6)) then
C     Term: del_ab 1/2 sum_kl (j k | l i ) rho_kl  (line 2 second term)
C     Implemented as (j k | i l)
C
      ISYML = ISDEL
C
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | i delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_3.60',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_3.60',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over i
C-----------------------------------------------------------------------
C
      DO 060 ISYMI = 1,NSYM
C
         ISALBE = MULD2H(ISYMI,ISYDIS)
         ISYMA  = MULD2H(ISYMI,ISYMTR)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_3.61',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_3.61',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | i delta ) for given i and
C           delta.
C-----------------------------------------------------------------------
C
         DO 061 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) 
     &                  * (I - 1) + KDSRHF

            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               
C              Generate ( k j | i delta) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 062 ISYMK = 1,NSYM

                  ISALFA = ISYMK
                  ISBETA = ISYMA

                  LSCR2  = NBAS(ISALFA)*NVIR(ISYMA)
C
                  KSCR2  = KEND2
                  KEND3  = KSCR2 + LSCR2
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.062',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.062',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = IMATAV(ISYMI,ISYMA) + 1
C
C              ( alpha, beta |..) * b(beta, a)  => ( alpha, a | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NVIR(ISYMA),
     &                    NBAS(ISBETA),HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    VBTRIAL(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C                     
C              (alpha, a | i,delta) * rho_{alpha delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISALFA,ISDEL)
C
                  DO IA = 1,NVIR(ISYMA)
                        DO ALPHA = 1,NBAS(ISALFA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR2 -(IA-1)*NBAS(ISALFA)+ALPHA)*
     &                      B2DENSIJ(KOFF6 +(IDEL-1)*NBAS(ISALFA)+ALPHA)
                        END DO
                  END DO
C             
  062       CONTINUE
C
  061    CONTINUE
C
  060 CONTINUE

      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 7)) then
C     Term: -1/2 sum_kl (a b | k l ) rho_bc  (line 2 third term)
C     Implemented as (k l | a b)
C
      ISYMB = ISDEL
C
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | a delta) 
C-----------------------------------------------------------------------
C
      LDSVIR = NDSVIR(ISDEL)
C
      KDSVIR  = 1
      KEND2  = KDSVIR + LDSVIR
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.70',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.70',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over a
C-----------------------------------------------------------------------
C
      DO 070 ISYMA = 1,NSYM
C
         ISALBE = MULD2H(ISYMA,ISYDIS)
         ISYMI  = MULD2H(ISYMA,ISYMTR)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.071',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.071',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | a delta ) for given a and
C           delta.
C-----------------------------------------------------------------------
C
         DO 071 IA = 1,NVIR(ISYMA)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMA) + NNBST(ISALBE) 
     &                  * (IA - 1) + KDSVIR

            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               
C              Generate (  | a delta) by multiplying with a double 
C              backtransformed density matrix for given a and delta 
C              in KSCR3.
C-----------------------------------------------------------------
C
            DO 072 ISALFA = 1,NSYM

                  ISBETA  = ISALFA
C                  ISBETA = MULD2H(ISALFA, ISALBE)

                  LSCR2  = NBAS(ISALFA)*NBAS(ISBETA)
                  LSCR3  = NBAS(ISALFA)*NBAS(ISBETA)
C
                  KSCR2  = KEND2
                  KSCR3  = KSCR2 + LSCR2
                  KEND3  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.072',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.072',' ',KEND3,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = IAODIS(ISALFA,ISBETA) + 1
C
C        (alpha, beta | a delta ) * rho(beta, alpha) => (alpha |a delta)
                  CALL DGEMM('N','T',NBAS(ISALFA),NBAS(ISBETA),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    B2DENSIJ(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C        (alpha, beta | a delta ) * rho(alpha, beta) => (alpha |a delta)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISALFA,ISBETA)
C
C                  DO BETA = 1,NBAS(ISBETA)
C                        DO ALPHA = 1,NBAS(ISALFA)
C                             WORK(KSCR2+(BETA-1)*NBAS(ISALFA)+ALPHA) = 
C     &                       WORK(KOFF2+(BETA-1)*NBAS(ISALFA)+ALPHA) *
C     &                       B2DENSIJ(KOFF6+(BETA-1)*NBAS(ISALFA)+ALPHA)
C                        END DO
C                  END DO       
C              ( alpha alpha | a, delta) * b_{delta i} => U(a,i)
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IMATAO(ISDEL,ISYMI)

                  DO I = 1,NRHF(ISYMI)
                        DO ALPHA = 1,NBAS(ISALFA)
                             RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) =
     &                       RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                       WORK(KSCR2-1+(ALPHA-1)*NBAS(ISALFA)+ALPHA)*
     &                       OBTRIAL(KOFF6 + (I-1)*NBAS(ISDEL) + IDEL)
                        END DO
                  END DO
 
  072       CONTINUE

  071    CONTINUE

  070 CONTINUE

      END IF
C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_3')
C
      RETURN
      END