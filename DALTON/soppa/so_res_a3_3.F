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

      I_ACT_TERM = 7
      
      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 1)) then
C     Term: - del_ab sum_cd (i j | c d ) rho_cd  (line 1 first term)
C     Implemented as (i j | c d)
C
      ISYMD = ISDEL
C
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
         ISYMGA = ISDEL
C
         ISALBE = MULD2H(ISYMGA,ISYDIS)
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
C                  ISYMA  = MULD2H(ISYMI,ISALBE)
                  ISALFA = ISYMA
                  ISBETA = ISYMI

                  LSCR2  = NBAS(ISALFA)*NRHF(ISYMI)
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
                  KOFF4  = IMATAV(ISYMA,ISYMI) + 1
C
C              ( alpha, beta |..) * C(beta, i)   => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),-ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C              b(a,alpha) * ( alpha, i | .. ) => ( a, i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    VBTRIAL(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)
C
C                  
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISYMGA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR3 -1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSAB(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  012       CONTINUE
C
  011    CONTINUE
C
      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 2)) then
C     Term: del_ab 1/2 sum_cd (j c | d i ) rho_cd  (line 1 second term)
C     Implemented as (c j | i d)
C
      ISYMD = ISDEL
C
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

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
      CALL SO_MEMMAX ('SO_RES_A3_3.20',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_3.20',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

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
         KSCR1  = KEND2
         KEND3  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND3
C
         CALL SO_MEMMAX ('SO_RES_A3_3.21',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_3.21',' ',KEND3,LWORK)
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
     &                  * (I - 1) + 1

            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               
C              Generate ( alpha j | i delta) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
                  ISYMC = ISDEL

                  ISALFA = ISYMC
                  ISBETA = ISYMA

                  LSCR2  = NBAS(ISALFA)*NVIR(ISYMA)
C
                  KSCR2  = KEND3
                  KEND4  = KSCR2 + LSCR2
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.022',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.022',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = IMATAV(ISYMA,ISYMI) + 1
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
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                       + WORK(KSCR2-1+(IA-1)*NBAS(ISALFA)+ALPHA)*
     &                      B2DENSAB(KOFF6 +(IDEL-1)*NBAS(ISALFA)+ALPHA)
                        END DO
                  END DO             
C
  021    CONTINUE
C
  020 CONTINUE
      
      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 3)) then
C     Term: del_ij sum_cd (a b | c d ) rho_cd  (line 1 third term)
C     Implemented as (a b | c d)
C
      ISYMD = ISDEL
C
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
         ISYMGA = ISDEL
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
C              given gamma and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 032 ISYMI = 1,NSYM

                  ISYMA  = MULD2H(ISYMI,ISYMTR)
C                  ISYMA  = MULD2H(ISYMI,ISALBE)
                  ISALFA = ISYMA
                  ISBETA = ISYMI

                  LSCR2  = NRHF(ISYMI)*NBAS(ISALFA)
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
                  KOFF3  = ILMVIR(ISYMA) + 1
                  KOFF4  = IT1AO(ISYMA,ISYMI) + 1
C
C               ( alpha, beta |..) * b(beta,i)  => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    OBTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                C(a,alpha) * ( alpha i | .. ) => ( a, i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF3),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)
C                     
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISYMGA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR3 -1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSAB(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  032       CONTINUE
C
  031 CONTINUE
C

      END IF


      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 5)) then
C     Term: - del_ab sum_kl (i j | k l ) rho_kl  (line 2 first term)
C     Implemented as (j i | k l)
C
      ISYML = ISDEL
C
C      LCDB  = NRHF(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
         ISYMGA = ISDEL
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

                  LSCR2  = NBAS(ISALFA)*NRHF(ISYMI)
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
                  KOFF4  = IMATAV(ISALFA,ISYMA) + 1
C
C              ( alpha, beta |..) * C(beta, i)   => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),-ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C              b(a,alpha) * ( alpha, i | .. ) => ( a, i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    VBTRIAL(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)
C                     
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISYMGA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR3 - 1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSIJ(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  052       CONTINUE
C
  051    CONTINUE
C
      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 6)) then
C     Term: del_ab 1/2 sum_kl (j k | l i ) rho_kl  (line 2 second term)
C     Implemented as (j k | i l)
C
      ISYML = ISDEL
C
C      LCDB  = NRHF(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

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
         KSCR1  = KEND2
         KEND3  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND3
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
                  ISYMK = ISDEL

                  ISALFA = ISYMK
                  ISBETA = ISYMA

                  LSCR2  = NBAS(ISALFA)*NVIR(ISYMA)
C
                  KSCR2  = KEND3
                  KEND4  = KSCR2 + LSCR2
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_3.062',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_3.062',' ',KEND4,LWORK)
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
     &                        +WORK(KSCR2-1+(IA-1)*NBAS(ISALFA)+ALPHA)*
     &                      B2DENSIJ(KOFF6 +(IDEL-1)*NBAS(ISALFA)+ALPHA)
                        END DO
                  END DO
C             
  061    CONTINUE
C
  060 CONTINUE

      END IF

            IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 7)) then
C     Term: del_ij sum_kl (a b | k l ) rho_kl  (line 2 third term)
C     Implemented as (a b | k l)
C
      ISYML = ISDEL
C
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Start loop over gamma
C-----------------------------------------------------------------------
C
         ISYMGA = ISDEL
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
         DO 071 GAMMA = 1,NBAS(ISYMGA)
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
            DO 072 ISYMI = 1,NSYM

                  ISYMA  = MULD2H(ISYMI,ISYMTR)
                  ISALFA = ISYMA
                  ISBETA = ISYMI

                  LSCR2  = NRHF(ISYMI)*NBAS(ISALFA)
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
                  KOFF3  = ILMVIR(ISYMA) + 1
                  KOFF4  = IT1AO(ISBETA,ISYMI) + 1
C
C               ( alpha, beta |..) * b(beta,i)  => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    OBTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                C(a,alpha) * ( alpha i | .. ) => ( a, i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF3),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)
C                     
C              (a i | gamma,delta) * rho_{gamma delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IAODIS(ISYMGA,ISDEL)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF5+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR3 -1 +(I-1)*NVIR(ISYMA)+IA) *
     &                      B2DENSIJ(KOFF6 +(IDEL-1)*NBAS(ISYMGA)+GAMMA)
                        END DO
                  END DO
C             
  072       CONTINUE
C
  071 CONTINUE

      END IF
C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_3')
C
      RETURN
      END