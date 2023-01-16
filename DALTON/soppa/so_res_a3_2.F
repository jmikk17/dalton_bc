C  /* Deck so_res_a3_1 */
      SUBROUTINE SO_RES_A3_2(RES1,LRES1,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                      BDENSIJ,BDENSIJT,LBDENSIJ,BDENSAB,BDENSABT,
     &                      LBDENSAB,XINT,CMO,LCMO,IDEL,ISDEL,ISYDIS,
     &                      ISYMTR,OBTRIAL,LOBTRIAL,VBTRIAL,LVBTRIAL,
     &                      WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Kasper F. Schaltz November 2022
C
C     PURPOSE: Calculate all terms belonging to the A'^(3)(1) matrix, 
C     which either requires three backtransformations or is calculated
C     by using the trialvector to change the first or second index.
C     This sub-routine combined with so_res_a3_1.F calculates all the 
C     terms in the A'^(3)(1)
C     The equation for the A'^(3)(1) matrix that has been used can be
C     found in my master thesis (currently unpublished)

#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (QUARTER = 0.25D0)
C
      DIMENSION RES1(LRES1), XINT(*)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION BDENSIJ(LBDENSIJ), BDENSAB(LBDENSAB)
      DIMENSION BDENSIJT(LBDENSIJ), BDENSABT(LBDENSAB)
      DIMENSION CMO(LCMO)
      DIMENSION WORK(LWORK)
      DIMENSION OBTRIAL(LOBTRIAL), VBTRIAL(LVBTRIAL)

      INTEGER COUNT
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
      CALL QENTER('SO_RES_A3_2')

C     Term:  - 2 sum_l (j i | l l ) rho_ba (line 1 first term)
C     Implemented as (l l | i j)
C
      ISYMJ = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Do  (alpha beta | gamma delta) * C(gamma,i) =  
C           (alpha beta | i delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.010',LWORK2)
      IF (LWORK2 .LT. 0) 
     &      CALL STOPIT('SO_RES_A3_2.010',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)
C-----------------------------------------------------------------------
C                                  
C           Start loop over l
C-----------------------------------------------------------------------            
      DO 010 ISYML = 1,NSYM
       
            ISALBE = MULD2H(ISYML,ISYDIS)
            ISALFA = ISYML
            ISBETA = ISYML      
            ISYMI  = MULD2H(ISYMTR, ISDEL)
            ISYMA  = MULD2H(ISYMI, ISYMTR)
            ISYMB  = MULD2H(ISYMTR, ISDEL)
 
            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.011',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.011',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | i delta ) for given i 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 011 I = 1,NRHF(ISYMI)
C
                 KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) 
     &                  * (I - 1) + KDSRHF

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( l l | i delta ) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C 
C
                  LSCR2  = NBAS(ISALFA)*NRHF(ISALFA)
                  LSCR3  = NRHF(ISBETA)*NRHF(ISBETA)
                  LSCR4  = NVIR(ISYMA)

                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KSCR4  = KSCR3 + LSCR3
                  KEND4  = KSCR4 + LSCR4
                  LWORK3 = LWORK - KEND4

                  CALL SO_MEMMAX ('SO_RES_A3_2.012',LWORK3)
                  IF (LWORK3 .LT. 0)
     &               CALL STOPIT('SO_RES_A3_2.012',' ',KEND4,LWORK)
     
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTL  = MAX(NRHF(ISYML),1)

                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYML) + 1

C              ( beta, alpha |..) * C(alpha, l) => ( beta, l | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYML),
     &                    NBAS(ISALFA),-TWO,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),NTOTBE)
C              ( l, beta | .. ) * C(beta, l) => ( l, l | ..)
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYML),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),NTOTBE,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR3),NTOTL)

C                                               ~
C               b_delta b * rho_ba => rho_delta a
C
                  KOFF5  = IABDEN(ISYMA,ISYMB) + 1
                  KOFF6  = IMATAV(ISDEL,ISYMB) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','N',1,NVIR(ISYMA),
     &                    NVIR(ISYMB),ONE,
     &                    DENSAB(KOFF5),1,
     &                    VBTRIAL(KOFF6),NTOTB,
     &                    ZERO,WORK(KSCR4),1)
C                                   ~          
C              (l l, i delta) * rho_a delta => U(a,i)
C
               KOFF7  = IT1AM(ISYMA,ISYMI) + (I-1)*NVIR(ISYMA)

               COUNT = 0
               DO IA = 1, NVIR(ISYMA)
                  DO L = 1,NRHF(ISYML)
                    RES1(KOFF7 + IA) = RES1(KOFF7+IA) +
     &              WORK(KSCR3 - 1 + NRHF(ISYML)*(L-1) + L) *
     &              WORK(KSCR4 - 1 + IA)
                  END DO
                  COUNT = COUNT + NVIR(ISYMA)
               END DO    

C
  011    CONTINUE
C
  010 CONTINUE

C     Term:  1/2 sum_l (j l | l i ) rho_ba (line 1 second term)
C     Implemented as (i l | l j)
C
      ISYMJ = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Do  (alpha beta | gamma delta) * C(gamma,l) =  
C           (alpha beta | l delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.020',LWORK2)
      IF (LWORK2 .LT. 0) 
     &      CALL STOPIT('SO_RES_A3_2.020',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)
C-----------------------------------------------------------------------
C                                  
C           Start loop over l
C-----------------------------------------------------------------------            
      DO 020 ISYML = 1,NSYM
       
            ISALBE = MULD2H(ISYML,ISYDIS)
            ISBETA = ISYML
            ISALFA = MULD2H(ISBETA,ISALBE)      
            ISYMI  = ISALFA
            ISYMA  = MULD2H(ISYMI, ISYMTR)
            ISYMB  = MULD2H(ISYMTR, ISDEL)
 
            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.021',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.021',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | l delta ) for given l 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 021 L = 1,NRHF(ISYML)
C
                 KOFF1 = IDSRHF(ISALBE,ISYML) + NNBST(ISALBE) 
     &                  * (L - 1) + KDSRHF

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( i l | l delta ) for
C              given l and delta in KSCR3.
C-----------------------------------------------------------------
C 
C
                  LSCR2  = NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)
                  LSCR4  = NVIR(ISYMA)

                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KSCR4  = KSCR3 + LSCR3
                  KEND4  = KSCR4 + LSCR4
                  LWORK3 = LWORK - KEND4

                  CALL SO_MEMMAX ('SO_RES_A3_2.022',LWORK3)
                  IF (LWORK3 .LT. 0)
     &               CALL STOPIT('SO_RES_A3_2.022',' ',KEND4,LWORK)
     
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTL  = MAX(NRHF(ISYML),1)

                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYML) + 1
                  KOFF4  = ILMRHF(ISYMI) + 1

C              ( beta, alpha |..) * C(alpha, l) => ( beta, l | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),1,
     &                    NBAS(ISALFA),HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),NTOTBE)
C              C(i, beta) * ( beta, i | .. )  => ( i, l | ..)
               CALL DGEMM('T','N',NRHF(ISYMI),1,
     &                    NBAS(ISBETA),ONE,
     &                    CMO(KOFF4),NTOTBE,
     &                    WORK(KSCR2),NTOTBE,
     &                    ZERO,WORK(KSCR3),1)

C                                               ~
C               b_delta b * rho_ba => rho_delta a
C
                  KOFF5  = IABDEN(ISYMA,ISYMB) + 1
                  KOFF6  = IMATAV(ISDEL,ISYMB) + 1

                  NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','N',1,NVIR(ISYMA),
     &                    NVIR(ISYMB),ONE,
     &                    DENSAB(KOFF5),1,
     &                    VBTRIAL(KOFF6),NTOTB,
     &                    ZERO,WORK(KSCR4),1)
C                                   ~          
C              (i l, l delta) * rho_a delta => U(a,i)
C
               KOFF7  = IT1AM(ISYMA,ISYMI)

               COUNT = 0
               DO I = 1, NRHF(ISYMI)
                  DO IA = 1,NVIR(ISYMA)
                    RES1(KOFF7+COUNT+IA) = RES1(KOFF7+COUNT+IA) +
     &              WORK(KSCR3 - 1 + I) *
     &              WORK(KSCR4 - 1 + IA)
                  END DO
                  COUNT = COUNT + NVIR(ISYMA)
               END DO    

C
  021    CONTINUE
C
  020 CONTINUE

C     Term: - 1/2 sum_c (a c | c b ) rho_ij (line 1 third term)
C     Implemented as (a c | c b)
C
      ISYMB = ISDEL
C      
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Do  (alpha beta | gamma delta) * C(gamma,c) =  
C           (alpha beta | c delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSVIR(ISDEL)
C
      KDSVIR  = 1
      KEND2  = KDSVIR + LDSVIR
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.030',LWORK2)
      IF (LWORK2 .LT. 0) 
     &      CALL STOPIT('SO_RES_A3_2.030',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)
C-----------------------------------------------------------------------
C                                  
C           Start loop over c
C-----------------------------------------------------------------------            
      DO 030 ISYMC = 1,NSYM
       
            ISALBE = MULD2H(ISYMC,ISYDIS)
            ISBETA = ISYMC
            ISALFA = MULD2H(ISBETA,ISALBE)      
            ISYMI  = ISALFA
            ISYMJ  = ISYMI
            ISYMA  = MULD2H(ISYMI, ISYMTR)
 
            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.031',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.031',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | c delta ) for given c 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 031 C = 1,NVIR(ISYMC)
C
                 KOFF1 = IDSVIR(ISALBE,ISYMC) + NNBST(ISALBE) 
     &                  * (C - 1) + KDSVIR

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( a c | c delta ) for
C              given c and delta in KSCR3.
C-----------------------------------------------------------------
C 
C
                  LSCR2  = NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)
                  LSCR4  = NVIR(ISYMA)

                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KSCR4  = KSCR3 + LSCR3
                  KEND4  = KSCR4 + LSCR4
                  LWORK3 = LWORK - KEND4

                  CALL SO_MEMMAX ('SO_RES_A3_2.032',LWORK3)
                  IF (LWORK3 .LT. 0)
     &               CALL STOPIT('SO_RES_A3_2.032',' ',KEND4,LWORK)
     
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTC  = MAX(NVIR(ISYMC),1)

                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMVIR(ISYMC) + 1
                  KOFF4  = ILMVIR(ISYMA) + 1

C              ( beta, alpha |..) * C(alpha, c) => ( beta, c | ..)
               CALL DGEMM('T','N',NBAS(ISBETA),1,
     &                    NBAS(ISALFA),-HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),NTOTBE)
C              C(a, beta) * ( beta, a | .. )  => ( a, c | ..)
               CALL DGEMM('T','N',NVIR(ISYMA),1,
     &                    NBAS(ISBETA),ONE,
     &                    CMO(KOFF4),NTOTBE,
     &                    WORK(KSCR2),NTOTBE,
     &                    ZERO,WORK(KSCR3),1)

C                                               ~
C               b_delta j * rho_ij => rho_delta i
C
                  KOFF5  = IIJDEN(ISYMI,ISYMJ) + 1
                  KOFF6  = IMATAO(ISDEL,ISYMJ) + 1

                  NTOTJ = MAX(NRHF(ISYMJ),1)

               CALL DGEMM('N','N',1,NRHF(ISYMI),
     &                    NRHF(ISYMJ),ONE,
     &                    DENSIJ(KOFF5),1,
     &                    OBTRIAL(KOFF6),NTOTJ,
     &                    ZERO,WORK(KSCR4),1)
C                                   ~          
C              (a c, c delta) * rho_i delta => U(a,i)
C
               KOFF7  = IT1AM(ISYMA,ISYMI)

               COUNT = 0
               DO I = 1, NRHF(ISYMI)
                  DO IA = 1,NVIR(ISYMA)
                    RES1(KOFF7+COUNT+IA) = RES1(KOFF7+COUNT+IA) +
     &              WORK(KSCR3 - 1 + I) *
     &              WORK(KSCR4 - 1 + IA)
                  END DO
                  COUNT = COUNT + NVIR(ISYMA)
               END DO    

C
  031    CONTINUE
C
  030 CONTINUE


C     Term: 1/2 sum_c (j i | a c ) rho_bc  (line 2 first term)
C     Implemented as (j i | a c)
C     Calculated by setting c = delta 
      ISYMC = ISDEL
C      
      LCDB  = NVIR(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
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
      CALL SO_MEMMAX ('SO_RES_A3_2.40',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.40',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over a
C-----------------------------------------------------------------------
      DO 040 ISYMA = 1,NSYM
       
            ISALBE = MULD2H(ISYMA,ISYDIS)
            ISYMJI = ISALBE

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.41',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.41',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | a delta ) for given a 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 041 IA = 1,NVIR(ISYMA)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMA) + NNBST(ISALBE) 
     &                  * (IA - 1) + KDSVIR

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C              Generate ( j i | a delta ) for
C              given a and delta in KSCR3.
C-----------------------------------------------------------------
C
                  DO 042 ISYMI = 1,NSYM

                        ISBETA = ISYMI
                        ISALFA = MULD2H(ISBETA, ISALBE)
                        ISYMB  = ISALFA

                        LSCR2  = NRHF(ISYMI)*NBAS(ISBETA)
                        LSCR3  = NRHF(ISYMI)
C
                        KSCR2  = KEND3
                        KSCR3  = KSCR2 + LSCR2
                        KEND4  = KSCR3 + LSCR3
                        LWORK3 = LWORK - KEND4
C
                        CALL SO_MEMMAX ('SO_RES_A3_2.42',LWORK3)
                        IF (LWORK3 .LT. 0)
     &                    CALL STOPIT('SO_RES_A3_2.42',' ',KEND4,LWORK)
C
                        NTOTAL = MAX(NBAS(ISALFA),1)
                        NTOTBE = MAX(NBAS(ISBETA),1)
                        NTOTB  = MAX(NVIR(ISYMB),1)
C
                        KOFF2  = ILMRHF(ISYMI) + 1
                        KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                        KOFF4  = IT1AO(ISYMA,ISYMI) + 1
C
C               ( alpha, beta |..) * C(beta, i)  => ( alpha, i | ..)
                        CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),HALF,
     &                    WORK(KOFF3),NTOTAL,
     &                    CMO(KOFF2),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                                                   ~
C               b(b,alpha) * ( alpha, i | .. ) => ( b, i | ..)
                        CALL DGEMM('T','N',NVIR(ISYMB),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),NTOTAL,
     &                    VBTRIAL(KOFF4),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTB)

C               ~
C              (b i | a delta) * rho_{b delta} => U(a,i)
C
                        KOFF5  = IT1AM(ISYMA,ISYMI)
                        KOFF6  = ILMVIR(ISYMB) - ILMVIR(1)
C
                        COUNT = 0   
                        DO I = 1, NRHF(ISYMI)
                          DO IB = 1, NVIR(ISYMB)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) = 
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                        WORK(KSCR3 + COUNT + IB)
     &                        * BDENSAB(KOFF6 +(1-IDEL)*NVIR(ISYMB)+ IB)
                          END DO
                          COUNT = COUNT + NVIR(ISYMB)
                        END DO
C
  042       CONTINUE
C
  041    CONTINUE
C
  040 CONTINUE


C     Term: 1/2 sum_c (a i | j c ) rho_bc  (line 2 second term)
C     Implemented as (a i | c j)
C
      ISYMJ = ISDEL
C
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | b delta) by multiplying with rho_gamma b
C-----------------------------------------------------------------------
C
      LDSVIR = NDSVIR(ISDEL)
C
      KDSVIR  = 1
      KEND2  = KDSVIR + LDSVIR
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.50',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.50',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),BDENSABT,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over b
C-----------------------------------------------------------------------
C
      DO 050 ISYMB = 1,NSYM
C
         ISALBE = MULD2H(ISYMB,ISYDIS)
         ISYMAI = ISALBE
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = 1
         KEND2  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.051',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.051',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | b delta ) for given b and
C           delta.
C-----------------------------------------------------------------------
C
         DO 051 B = 1,NVIR(ISYMB)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMB) + NNBST(ISALBE) 
     &                  * (B - 1) + KDSVIR

            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C

C-----------------------------------------------------------------
C                               
C              Generate ( a i | b delta) for
C              given b and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 052 ISYMA = 1,NSYM

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
                  CALL SO_MEMMAX ('SO_RES_A3_2.052',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.052',' ',KEND3,LWORK)
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

C                     
C              (a,i | j,delta) * b_{delta b} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IMATAV(ISDEL,ISYMB)
C
                  COUNT  = 0
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+COUNT+IA) =RES1(KOFF5+COUNT+IA)
     &                        + WORK(KSCR3 - 1 + COUNT + IA) *
     &                        VBTRIAL(KOFF6 + B)
                        END DO
                        COUNT = COUNT + NVIR(ISYMA)
                  END DO
C             
  052       CONTINUE
C
  051    CONTINUE
C
  050 CONTINUE


C     Term: 1/4 del_ij sum_cl (c l | l b ) rho_ca  (line 4 second term)
C     Implemented as (l b | l c)
C     Calculated by setting c = delta
      ISYMC  = ISDEL
C      
      LCDB   = NVIR(ISYMC)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | l delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.130',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.130',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over l
C-----------------------------------------------------------------------
      DO 130 ISYML = 1,NSYM
       
            ISALBE = MULD2H(ISYML,ISYDIS)
            ISALFA = ISYML
            ISBETA = MULD2H(ISALFA, ISALBE)
            ISYMI  = MULD2H(ISYMTR, ISBETA)
            ISYMA  = MULD2H(ISYMI, ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.131',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.131',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | l delta ) for given l 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 131 L = 1,NRHF(ISYML)
C
                  KOFF1 = IDSRHF(ISALBE,ISYML) + NNBST(ISALBE) 
     &                  * (L - 1) + KDSRHF

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           ~
C              Generate ( l i | l delta ) for
C              given l and delta in KSCR3.
C-----------------------------------------------------------------
C 
C
                  LSCR2  = NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)
C
                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.132',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.132',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
C
                  KOFF2  = ILMRHF(ISYML) + (L-1)*NBAS(ISYML)+1
                  KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF4  = IT1AO(ISYMA,ISYMI) + 1
C
C              C(l, alpha) *( alpha, beta |..)  => ( l, beta | ..)
                  CALL DGEMM('T','N',1,NBAS(ISBETA),
     &                    NBAS(ISALFA),QUARTER,
     &                    CMO(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),1)

C                                                   ~
C              ( l, beta | .. ) * b(beta,i) => ( l, i | ..)
                  CALL DGEMM('N','N',1,NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),1,
     &                    OBTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),1)

C                 ~
C              (l,i | l,delta) * rho_{a delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6 = ILMVIR(ISYMA) - ILMVIR(1)
C
                  COUNT = 0
                  DO I = 1, NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+COUNT+IA) =RES1(KOFF5+COUNT+IA)
     &                        + WORK(KSCR3 - 1 +(I-1)*NRHF(ISYMI)+ L) *
     &                        BDENSAB(KOFF6 + (IDEL-1)*NVIR(ISYMA)+ IA)
                        END DO
                        COUNT = COUNT + NVIR(ISYMA)
                  END DO
C
  131    CONTINUE
C
  130 CONTINUE


C     Term: 1/2 del_ab sum_kl (j k | l l ) rho_ik  (line 5 first term)
C     Implemented as (l l | k j)
C     Calculated by setting j = delta and multiplying with rho_i gamma 
C     first
      ISYMJ = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Do (alpha beta | gamma delta) * rho_i gamma = 
C           (alpha beta | i delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.140',LWORK2)
      IF (LWORK2 .LT. 0) 
     &      CALL STOPIT('SO_RES_A3_2.140',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),BDENSIJT,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)
C-----------------------------------------------------------------------
C                                  
C           Start loop over i
C-----------------------------------------------------------------------            
      DO 140 ISYMI = 1,NSYM
       
            ISALBE = MULD2H(ISYMI,ISYDIS)     
            ISYMA  = MULD2H(ISYMI,ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.141',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.141',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | i delta ) for given i 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 141 I = 1,NRHF(ISYMI)
C
                 KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) 
     &                  * (I - 1) + KDSRHF
                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( l l | i delta ) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C 
                  DO 142 ISYML = 1,NSYM
C
                  ISALFA = ISYML
                  ISBETA = ISYML

                  LSCR2  = NBAS(ISBETA)*NRHF(ISBETA)
                  LSCR3  = NRHF(ISBETA)*NRHF(ISBETA)

                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4

                  CALL SO_MEMMAX ('SO_RES_A3_2.142',LWORK3)
                  IF (LWORK3 .LT. 0)
     &               CALL STOPIT('SO_RES_A3_2.142',' ',KEND4,LWORK)
     
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTL  = MAX(NRHF(ISYML),1)

                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMRHF(ISYML) + 1

C              ( alpha, beta |..) * C(beta, l) => ( alpha, l | ..)
               CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYML),
     &                    NBAS(ISBETA),HALF,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C              ( l, alpha | .. ) * C(alpha, l) => ( l, l | ..)
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYML),
     &                    NBAS(ISALFA),ONE,
     &                    WORK(KSCR2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTL)

C                 
C              (l,l | i,delta) * b_delta a => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI) + (I-1)*NVIR(ISYMA)
                  KOFF6  = IMATAV(ISYMI,ISYMA)
C
                  DO IA = 1, NVIR(ISYMA)
                        DO L = 1,NRHF(ISYML)
                              RES1(KOFF5 + IA) =
     &                         RES1(KOFF5 + IA)
     &                        + WORK(KSCR3-1 + (L-1)*NRHF(ISYML) + L) *
     &                        VBTRIAL(KOFF6 + NBAS(ISDEL)*(IA-1)+ IDEL)
                        END DO 
                  END DO
C
  142       CONTINUE
C  
  141    CONTINUE
C
  140 CONTINUE

C     Term: -1/2 del_ab sum_kl (j l | l k ) rho_ik  (line 5 second term)
C     Implemented as (l k | l j)
C     Calculated by setting j = delta
      ISYMJ = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | l delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.150',LWORK2)
      IF (LWORK2 .LT. 0)
     &    CALL STOPIT('SO_RES_A3_2.150',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over l
C-----------------------------------------------------------------------
      DO 150 ISYML = 1,NSYM
       
            ISALBE = MULD2H(ISYML,ISYDIS)
            ISALFA = ISYML
            ISBETA = MULD2H(ISALFA, ISALBE)
            ISYMI  = ISBETA
            ISYMA  = MULD2H(ISYMI, ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.151',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.151',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | l delta ) for given l 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 151 L = 1,NRHF(ISYML)
C
                  KOFF1 = IDSRHF(ISALBE,ISYML) + NNBST(ISALBE) 
     &                  * (L - 1) + KDSRHF

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( l i | l delta ) for
C              given l and delta in KSCR3.
C-----------------------------------------------------------------
C 
                  LSCR2  = NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)
C
                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.152',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.152',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTI  = MAX(NRHF(ISYMI),1)
C
                  KOFF2  = ILMRHF(ISYML) + (L-1)*NBAS(ISYML) + 1
                  KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF4  = ILMRHF(ISYMI) + 1
C
C              C(l, alpha) *( alpha, beta |..)  => ( l, beta | ..)
                  CALL DGEMM('T','N',1,NBAS(ISBETA),
     &                    NBAS(ISALFA),-HALF,
     &                    CMO(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),1)

C                                                  
C              ( l, beta | .. ) * rho(beta,i) => ( l, i | ..)
                  CALL DGEMM('N','T',1,NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),1,
     &                    BDENSIJ(KOFF4),NTOTI,
     &                    ZERO,WORK(KSCR3),1)

C                 ~
C              (l,i | l,delta) * b_{delta a} => U(a,i)
C 
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IMATAV(ISYMA,ISYMI)
C
                  COUNT  = 0
                  DO I = 1, NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+COUNT+IA) =RES1(KOFF5+COUNT+IA)
     &                        + WORK(KSCR3 - 1 + I) *
     &                        VBTRIAL(KOFF6 + NBAS(ISDEL)*(IA-1)+ IDEL)
                        END DO
                        COUNT = COUNT + NVIR(ISYMA)
                  END DO
C
  151    CONTINUE
C
  150 CONTINUE

C     Term: 1/4 del_ab sum_kc (j c | c k ) rho_ik  (line 5 third term)
C     Implemented as (c k | c j)
C     Calculated by setting j = delta 
      ISYMJ = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | c delta)
C-----------------------------------------------------------------------
C
      LDSVIR = NDSVIR(ISDEL)
C
      KDSVIR  = 1
      KEND2  = KDSVIR + LDSVIR
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.160',LWORK2)
      IF (LWORK2 .LT. 0)
     &    CALL STOPIT('SO_RES_A3_2.160',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over c
C-----------------------------------------------------------------------
      DO 160 ISYMC = 1,NSYM
       
            ISALBE = MULD2H(ISYMC,ISYDIS)
            ISALFA = ISYMC
            ISBETA = MULD2H(ISALFA, ISALBE)
            ISYMI  = ISBETA
            ISYMA  = MULD2H(ISYMI, ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.161',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.161',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | c delta ) for given c 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 161 IC = 1,NVIR(ISYMC)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMC) + NNBST(ISALBE) 
     &                  * (IC - 1) + KDSVIR

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( c i | c delta ) for
C              given c and delta in KSCR3.
C-----------------------------------------------------------------
C 
                  LSCR2  = NBAS(ISBETA)
                  LSCR3  = NRHF(ISYMI)
C
                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.162',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.162',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
C
                  KOFF2  = ILMVIR(ISYMC) + (C-1)*NBAS(ISYMC) + 1
                  KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF4  = ILMRHF(ISYMI) + 1
C
C              C(alpha, c) *( alpha, beta |..)  => ( c, beta | ..)
                  CALL DGEMM('T','N',1,NBAS(ISBETA),
     &                    NBAS(ISALFA),QUARTER,
     &                    CMO(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),1)

C                                                     
C              ( c, beta | .. ) * rho(beta,i) => ( c, i | ..)
                  CALL DGEMM('N','T',1,NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),1,
     &                    BDENSIJ(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),1)

C                 
C              (c,i | c,delta) * b_{delta_a} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFF6  = IMATAV(ISYMA,ISYMI)
C
                  COUNT  = 0
                  DO I = 1, NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF5+COUNT+IA) =RES1(KOFF5+COUNT+IA)
     &                        + WORK(KSCR3 - 1 + I) *
     &                        VBTRIAL(KOFF6 + (IA-1)*NBAS(ISDEL) + IDEL)
                        END DO
                        COUNT = COUNT + NVIR(ISYMA)
                  END DO
  161    CONTINUE
C
  160 CONTINUE

C     Term:  1/4 del_ab sum_kc (c c | j k ) rho_ik (line 5 fourth term)
C     Implemented as (c c | k j)
C     Calculated by setting j = delta and starting with multiplying 
C     with rho_i gamma

      ISYMJ = ISDEL
C      
      LCDB  = NRHF(ISDEL)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Do (alpha beta | gamma delta) * rho_i gamma = 
C           (alpha beta | i delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISDEL)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.170',LWORK2)
      IF (LWORK2 .LT. 0) 
     &      CALL STOPIT('SO_RES_A3_2.170',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),BDENSIJ,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)
C-----------------------------------------------------------------------
C                                  
C           Start loop over c
C-----------------------------------------------------------------------            
      DO 170 ISYMC = 1,NSYM
       
            ISALBE = MULD2H(ISYMC,ISYDIS)
            ISALFA = ISYMC
            ISBETA = ISYMC      
            ISYMI  = MULD2H(ISYMTR, ISDEL)
            ISYMA  = MULD2H(ISYMI, ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.171',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.171',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | i delta ) for given i 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 171 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) 
     &                  * (I - 1) + KDSRHF
                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( c c | i delta ) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C 
C
                  LSCR2  = NBAS(ISBETA)*NVIR(ISBETA)
                  LSCR3  = NVIR(ISBETA)*NVIR(ISBETA)

                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4

                  CALL SO_MEMMAX ('SO_RES_A3_2.172',LWORK3)
                  IF (LWORK3 .LT. 0)
     &               CALL STOPIT('SO_RES_A3_2.172',' ',KEND4,LWORK)
     
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTC  = MAX(NVIR(ISYMC),1)

                  KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF3  = ILMVIR(ISYMC) + 1

C              ( alpha, beta |..) * C(beta, c) => ( alpha, c | ..)
               CALL DGEMM('N','N',NBAS(ISALFA),NVIR(ISYMC),
     &                    NBAS(ISBETA),QUARTER,
     &                    WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C              ( c, alpha | .. ) * C(alpha, c) => ( c, c | ..)
               CALL DGEMM('T','N',NRHF(ISYML),NRHF(ISYML),
     &                    NBAS(ISALFA),ONE,
     &                    WORK(KSCR2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTC)

C                 
C              (c,c | i,delta) * b_delta a => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI) + (I-1)*NVIR(ISYMA)
                  KOFF6  = IMATAV(ISYMI,ISYMA)
C
                  DO IA = 1, NVIR(ISYMA)
                        DO IC = 1,NVIR(ISYMC)
                              RES1(KOFF5+IA) =RES1(KOFF5+IA)
     &                        + WORK(KSCR3-1 + (IC-1)*NVIR(ISYMC) + IC)*
     &                        VBTRIAL(KOFF6 + (IA-1)*NBAS(ISDEL)+IDEL)
                        END DO
                  END DO
C
  171    CONTINUE
C
  170 CONTINUE

C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_2')
C
      RETURN
      END