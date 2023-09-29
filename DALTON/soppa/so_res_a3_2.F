C  /* Deck so_res_a3_3 */
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
C     This sub-routine combined with so_res_a3_1.F and so_res_a3_3.F 
C     calculates all the terms in the A'^(3)(1)
C     The equation for the A'^(3)(1) matrix that has been used can be
C     found in my master thesis (currently unpublished)

#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
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

      I_ACT_TERM = 0

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 9)) then
C     Term: - sum_c (a i | j c ) rho_bc  (line 3 first term)
C     Implemented as (a i | j c)
C
      ISYMC = ISDEL
C
C      LCDB  = NRHF(ISDEL)
C      IF (LCDB.LE.0) GOTO 999

C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | b delta) by multiplying with b_gamma_b
C-----------------------------------------------------------------------
C
      LDSVIR = NDSVIR(ISDEL)
C
      KDSVIR  = 1
      KEND2  = KDSVIR + LDSVIR
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.90',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.90',' ',KEND2,LWORK)
C
      KOFF1 = 1
      ISYMLP = 1
      DO ISYMG = 1,NSYM
C
         ISYMC  = MULD2H(ISYMLP,ISYMG)
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         NNBSAB = MAX(NNBST(ISYMAB),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         KOFF2  = 1 + IGLMVI(ISYMG,ISYMC) - IGLMVI(1,1)
         KOFF3  = 1 + IDSVIR(ISYMAB,ISYMC)
C
         CALL DGEMM('N','N',NNBST(ISYMAB),NVIR(ISYMC),NBAS(ISYMG),
     *              ONE,XINT(KOFF1),NNBSAB,VBTRIAL(KOFF2),NBASG,
     *              ZERO,WORK(KOFF3),NNBSAB)
C
         KOFF1 = KOFF1 + NNBST(ISYMAB)*NBAS(ISYMG)
      END DO
C-----------------------------------------------------------------------
C                                  
C           Start loop over b
C-----------------------------------------------------------------------
C
      DO 090 ISYMB = 1,NSYM
C
         ISALBE = MULD2H(ISYMB,ISYDIS)
C
         LSCR1  = N2BST(ISALBE)
C
         KSCR1  = KEND2
         KEND3  = KSCR1 + LSCR1
         LWORK2 = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_RES_A3_1.91',LWORK2)
         IF (LWORK2 .LT. 0)
     &     CALL STOPIT('SO_RES_A3_1.91',' ',KEND2,LWORK)
C
C-----------------------------------------------------------------------
C                                              
C           Get a squared set of ( alfa beta | b delta ) for given b and
C           delta.
C-----------------------------------------------------------------------
C
         DO 091 IB = 1,NVIR(ISYMB)
C
                  KOFF4 = IDSVIR(ISALBE,ISYMB) + NNBST(ISALBE)
     &                  * (IB - 1) + 1

            CALL CCSD_SYMSQ(WORK(KOFF4),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                               
C              Generate ( a i | b delta) for
C              given b and delta in KSCR3.
C-----------------------------------------------------------------
C
            DO 092 ISYMI = 1,NSYM

                  ISBETA = ISYMI
                  ISYMA  = MULD2H(ISYMI,ISYMTR)
                  ISALFA = ISYMA

                  LSCR2  = NBAS(ISALFA)*NRHF(ISYMI)
                  LSCR3  = NVIR(ISYMA)*NRHF(ISYMI)
C
                  KSCR2  = KEND3
                  KSCR3  = KSCR2 + LSCR2
                  KEND4  = KSCR3 + LSCR3
                  LWORK3 = LWORK - KEND4
C
                  CALL SO_MEMMAX ('SO_RES_A3_2.92',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.92',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  KOFF5  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF6  = ILMRHF(ISYMI) + 1
                  KOFF7  = ILMVIR(ISYMA) + 1
C
C              ( alpha, beta |..) * C(beta, i)   => ( alpha, i | ..)
                  CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),-ONE,
     &                    WORK(KOFF5),NTOTAL,
     &                    CMO(KOFF6),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C
C              C(a,alpha) * ( alpha, i | .. ) => ( a i | ..)
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF7),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)

C                     
C              (a,i | b,delta) * rho{b_delta} => U(a,i)
C
                  KOFF8  = IT1AM(ISYMA,ISYMI)
                  KOFF9  = ILMVIR(ISYMB) - ILMVIR(1)
C
                  DO I = 1,NRHF(ISYMI)
                        DO IA = 1,NVIR(ISYMA)
                              RES1(KOFF8+(I-1)*NVIR(ISYMA)+IA) =
     &                         RES1(KOFF8+(I-1)*NVIR(ISYMA)+IA)
     &                        + WORK(KSCR3 -1 + (I-1)*NVIR(ISYMA)+IA) *
     &                        BDENSAB(KOFF9 + (IDEL-1)*NVIR(ISYMB)+IB)
                        END DO
                  END DO
C             
  092       CONTINUE
C
  091    CONTINUE
C
  090 CONTINUE
                  
      END IF


      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 10)) then
C     Term: 1/2 sum_c (j i | c a ) rho_bc  (line 3 second term)
C     Implemented as (j i | a c)
C     Calculated by setting c = delta 
      ISYMC = ISDEL
C      
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999
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
      CALL SO_MEMMAX ('SO_RES_A3_2.100',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.100',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK,LWORK,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over a
C-----------------------------------------------------------------------
      DO 100 ISYMA = 1,NSYM
       
            ISALBE = MULD2H(ISYMA,ISYDIS)
            ISYMI  = MULD2H(ISYMA,ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.101',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.101',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | a delta ) for given a 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 101 IA = 1,NVIR(ISYMA)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMA) + NNBST(ISALBE) 
     &                  * (IA - 1) + 1

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C              Generate ( j i | a delta ) for
C              given a and delta in KSCR3.
C-----------------------------------------------------------------
C
C                  DO 102 ISYMB = 1,NSYM
                        ISYMB  = MULD2H(ISYMI,ISALBE)

                        ISALFA = ISYMB
                        ISBETA = ISYMI

                        LSCR2  = NRHF(ISYMI)*NBAS(ISBETA)
                        LSCR3  = NVIR(ISYMB)*NRHF(ISYMI)
C
                        KSCR2  = KEND3
                        KSCR3  = KSCR2 + LSCR2
                        KEND4  = KSCR3 + LSCR3
                        LWORK3 = LWORK - KEND4
C
                        CALL SO_MEMMAX ('SO_RES_A3_2.102',LWORK3)
                        IF (LWORK3 .LT. 0)
     &                    CALL STOPIT('SO_RES_A3_2.102',' ',KEND4,LWORK)
C
                        NTOTAL = MAX(NBAS(ISALFA),1)
                        NTOTBE = MAX(NBAS(ISBETA),1)
                        NTOTB  = MAX(NVIR(ISYMB),1)
C
                        KOFF2  = ILMRHF(ISYMI) + 1
                        KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                        KOFF4  = IMATAV(ISYMB,ISALFA) + 1
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
     &                    NBAS(ISALFA),ONE,
     &                    VBTRIAL(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTB)

C               ~
C              (b i | a delta) * rho_{b delta} => U(a,i)
C
                        KOFF5  = IT1AM(ISYMA,ISYMI)
                        KOFF6  = ILMVIR(ISYMB) - ILMVIR(1)
C
                        DO I = 1, NRHF(ISYMI)
                          DO IB = 1, NVIR(ISYMB)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) = 
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                        WORK(KSCR3 - 1 +(I-1)*NVIR(ISYMB) + IB)
     &                        * BDENSAB(KOFF6 +(IDEL-1)*NVIR(ISYMB)+ IB)
                          END DO
                        END DO
C
C  102       CONTINUE
C
  101    CONTINUE
C
  100 CONTINUE

      END IF
      
C      I_ACT_TERM = 11
      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 11)) then
C     Term: - sum_c (b j | i c ) rho_ac  (line 3 third term)
C     Implemented as (b j | i c)
C     Calculated by setting c = delta 
      ISYMC = ISDEL
C      
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999
C
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
      CALL SO_MEMMAX ('SO_RES_A3_2.110',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.110',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over i
C-----------------------------------------------------------------------
      DO 110 ISYMI = 1,NSYM
       
            ISALBE = MULD2H(ISYMI,ISYDIS)
            ISYMBJ = ISALBE
            ISYMA  = MULD2H(ISYMTR, ISYMI)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_2.111',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_2.111',' ',KEND3,LWORK)
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | i delta ) for given i 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 111 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE) 
     &                  * (I - 1) + 1

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           ~
C              Generate ( j j | i delta ) for
C              given i and delta in KSCR3.
C-----------------------------------------------------------------
C
                  DO 112 ISYMJ = 1,NSYM

                        ISBETA = ISYMJ
                        ISALFA = ISYMJ

                        LSCR2  = NRHF(ISYMJ)*NBAS(ISBETA)
                        LSCR3  = NRHF(ISYMJ)*NRHF(ISYMJ)
C
                        KSCR2  = KEND3
                        KSCR3  = KSCR2 + LSCR2
                        KEND4  = KSCR3 + LSCR3
                        LWORK3 = LWORK - KEND4
C
                        CALL SO_MEMMAX ('SO_RES_A3_2.112',LWORK3)
                        IF (LWORK3 .LT. 0)
     &                    CALL STOPIT('SO_RES_A3_2.112',' ',KEND4,LWORK)
C
                        NTOTAL = MAX(NBAS(ISALFA),1)
                        NTOTBE = MAX(NBAS(ISBETA),1)
                        NTOTJ  = MAX(NRHF(ISYMJ),1)
C
                        KOFF2  = ILMRHF(ISYMJ) + 1
                        KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                        KOFF4  = IT1AO(ISBETA,ISYMJ) + 1
C                                                            ~
C               ( alpha, beta |..) * b(beta, j)  => ( alpha, j | ..)
                        CALL DGEMM('N','N',NBAS(ISBETA),NRHF(ISYMJ),
     &                    NBAS(ISBETA),-ONE,
     &                    WORK(KOFF3),NTOTAL,
     &                    OBTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                                                      ~
C               C(j,alpha) * ( alpha, i | .. ) => ( j, j | ..)
                        CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMJ),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF2),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTJ)

C               ~         
C              (j j | i delta) * rho_{a delta} => U(a,i)
C
                        KOFF5  = IT1AM(ISYMA,ISYMI)
                        KOFF6  = ILMVIR(ISYMA) - ILMVIR(1)
C
                        DO IA = 1, NVIR(ISYMA)
                          DO J = 1, NRHF(ISYMJ)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) = 
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                        WORK(KSCR3 - 1 + (J-1)*NRHF(ISYMJ) + J)
     &                        * BDENSAB(KOFF6 +(IDEL-1)*NVIR(ISYMA)+ IA)
                          END DO
                        END DO
C
  112       CONTINUE
C
  111    CONTINUE
C
  110 CONTINUE

      END IF


      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 14)) then
C     Term: - 1/2 sum_l (a b | i l ) rho_jl  (line 4 second term)
C     Implemented as (i l | a b)
C     Calculated by setting b = delta 
      ISYMB = ISDEL
C      
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999
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
      CALL SO_MEMMAX ('SO_RES_A3_2.140',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.140',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK,LWORK,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over a
C-----------------------------------------------------------------------
      DO 140 ISYMA = 1,NSYM
       
            ISALBE = MULD2H(ISYMA,ISYDIS)
            ISYMI  = MULD2H(ISYMTR, ISYMA)

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
C           Get a squared set of ( alpha beta | a delta ) for given a 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 141 IA = 1,NVIR(ISYMA)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMA) + NNBST(ISALBE) 
     &                  * (IA - 1) + 1

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( i j | a delta ) for
C              given a and delta in KSCR3.
C-----------------------------------------------------------------
C
C                  DO 142 ISYMJ = 1,NSYM
                        ISYMJ  = MULD2H(ISYMI,ISALBE)

                        ISBETA = ISYMJ
                        ISALFA = ISYMI

                        LSCR2  = NRHF(ISYMJ)*NBAS(ISBETA)
                        LSCR3  = NRHF(ISYMI)*NRHF(ISYMJ)
C
                        KSCR2  = KEND3
                        KSCR3  = KSCR2 + LSCR2
                        KEND4  = KSCR3 + LSCR3
                        LWORK3 = LWORK - KEND4
C
                        CALL SO_MEMMAX ('SO_RES_A3_2.142',LWORK3)
                        IF (LWORK3 .LT. 0)
     &                    CALL STOPIT('SO_RES_A3_2.142',' ',KEND4,LWORK)
C
                        NTOTAL = MAX(NBAS(ISALFA),1)
                        NTOTBE = MAX(NBAS(ISBETA),1)
                        NTOTI  = MAX(NRHF(ISYMI),1)
                        NTOTJ  = MAX(NRHF(ISYMJ),1)
C
                        KOFF2  = ILMRHF(ISYMJ) + 1
                        KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                        KOFF4  = ILMRHF(ISYMI) + 1
C                                                           
C               ( alpha, beta |..) * rho(beta, j)  => ( alpha, j | ..)
                        CALL DGEMM('N','T',NBAS(ISALFA),NRHF(ISYMJ),
     &                    NBAS(ISBETA),-HALF,
     &                    WORK(KOFF3),NTOTAL,
     &                    BDENSIJ(KOFF2),NTOTJ,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                                                      
C               C(i,alpha) * ( alpha, j | .. ) => ( i, j | ..)
                        CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTI)

C                        
C              (i j | a delta) * b_{delta j} => U(a,i)
C
                        KOFF5  = IT1AM(ISYMA,ISYMI)
                        KOFF6  = IT1AO(ISDEL,ISYMJ)
C
                        DO I = 1, NRHF(ISYMI)
                          DO J = 1, NRHF(ISYMJ)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) = 
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                        WORK(KSCR3 - 1 + (J-1)*NRHF(ISYMI) + I)
     &                        * OBTRIAL(KOFF6 +(J-1)*NBAS(ISDEL)+IDEL)
                          END DO
                        END DO
C
C  142       CONTINUE
C
  141    CONTINUE
C
  140 CONTINUE

      END IF


      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 15)) then
C     Term: sum_l (b j | a l ) rho_il  (line 4 third term)
C     Implemented as (a l | j b)
C     Calculated by setting b = delta 
      ISYMB = ISDEL
C     
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | j delta)
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
C           Start loop over j
C-----------------------------------------------------------------------
      DO 150 ISYMJ = 1,NSYM
       
            ISALBE = MULD2H(ISYMJ,ISYDIS)
            ISYMAI = ISALBE

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
C           Get a squared set of ( alpha beta | j delta ) for given j 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 151 J = 1,NRHF(ISYMJ)
C
                  KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE) 
     &                  * (J - 1) + 1

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( a i | j delta ) for
C              given j and delta in KSCR3.
C-----------------------------------------------------------------
C
                  DO 152 ISYMI = 1,NSYM

                        ISBETA = ISYMI
                        ISYMA  = MULD2H(ISYMI,ISYMTR)
                        ISALFA = ISYMA

                        LSCR2  = NRHF(ISYMI)*NBAS(ISBETA)
                        LSCR3  = NVIR(ISYMA)*NRHF(ISYMI)
C
                        KSCR2  = KEND3
                        KSCR3  = KSCR2 + LSCR2
                        KEND4  = KSCR3 + LSCR3
                        LWORK4 = LWORK - KEND4
C
                        CALL SO_MEMMAX ('SO_RES_A3_2.152',LWORK4)
                        IF (LWORK4 .LT. 0)
     &                    CALL STOPIT('SO_RES_A3_2.152',' ',KEND4,LWORK)
C
                        NTOTAL = MAX(NBAS(ISALFA),1)
                        NTOTBE = MAX(NBAS(ISBETA),1)
                        NTOTA  = MAX(NVIR(ISYMA),1)
                        NTOTI  = MAX(NRHF(ISYMI),1)
C
                        KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
                        KOFF3  = ILMRHF(ISYMI) + 1
                        KOFF4  = ILMVIR(ISYMA) + 1
C                                                           
C               ( alpha, beta |..) * rho(beta, i)  => ( alpha, i | ..)
                        CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    BDENSIJT(KOFF3),NTOTBE,
     &                    ZERO,WORK(KSCR2),NTOTAL)
C                                                      
C               C(a,alpha) * ( alpha, i | .. ) => ( a, i | ..)
                        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTA)

C               
C              (a i | j delta) * b_{delta j} => U(a,i)
C
                        KOFF5  = IT1AM(ISYMA,ISYMI)
                        KOFF6  = IT1AO(ISDEL,ISYMJ)
C
                        DO I = 1, NRHF(ISYMI)
                          DO IA = 1, NVIR(ISYMA)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) = 
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                        WORK(KSCR3 - 1 + (I-1)*NVIR(ISYMA) + IA)
     &                        * OBTRIAL(KOFF6 +(J-1)*NBAS(ISDEL)+IDEL)
                          END DO
                        END DO
C
  152       CONTINUE
C
  151    CONTINUE
C
  150 CONTINUE

      END IF

      IF ((I_ACT_TERM .EQ. 0) .or. (I_ACT_TERM .EQ. 16)) then
C     Term: - 1/2 sum_c (a b | j l ) rho_il  (line 4 fourth term)
C     Implemented as (j l | a b)
C     Calculated by setting b = delta
      ISYMB = ISDEL
C      
C      LCDB  = NVIR(ISDEL)
C      IF (LCDB.LE.0) GOTO 999
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
      CALL SO_MEMMAX ('SO_RES_A3_2.160',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.160',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK,LWORK,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over a
C-----------------------------------------------------------------------
      DO 160 ISYMA = 1,NSYM
       
            ISALBE = MULD2H(ISYMA,ISYDIS)
            ISYMI  = MULD2H(ISYMTR, ISYMA)

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
C           Get a squared set of ( alpha beta | a delta ) for given a 
C           and delta.
C-----------------------------------------------------------------------
C
            DO 161 IA = 1,NVIR(ISYMA)
C
                  KOFF1 = IDSVIR(ISALBE,ISYMA) + NNBST(ISALBE) 
     &                  * (IA - 1) + 1

                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           
C              Generate ( j i | a delta ) for
C              given a and delta in KSCR3.
C-----------------------------------------------------------------
C
C                  DO 162 ISYMJ = 1,NSYM
                        ISYMJ  = MULD2H(ISYMI,ISALBE)

                        ISALFA = ISYMJ
                        ISBETA = ISYMI

                        LSCR2  = NRHF(ISYMI)*NBAS(ISALFA)
                        LSCR3  = NRHF(ISYMI)*NRHF(ISYMJ)
C
                        KSCR2  = KEND3
                        KSCR3  = KSCR2 + LSCR2
                        KEND4  = KSCR3 + LSCR3
                        LWORK3 = LWORK - KEND4
C
                        CALL SO_MEMMAX ('SO_RES_A3_2.162',LWORK3)
                        IF (LWORK3 .LT. 0)
     &                    CALL STOPIT('SO_RES_A3_2.162',' ',KEND4,LWORK)
C
                        NTOTAL = MAX(NBAS(ISALFA),1)
                        NTOTBE = MAX(NBAS(ISBETA),1)
                        NTOTI  = MAX(NRHF(ISYMI),1)
                        NTOTJ  = MAX(NRHF(ISYMJ),1)
C
                        KOFF2  = ILMRHF(ISYMI) + 1
                        KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                        KOFF4  = ILMRHF(ISYMJ) + 1
C                                                           
C               ( alpha, beta |..) * rho(beta, i)  => ( alpha, i | ..)
                        CALL DGEMM('N','T',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),-HALF,
     &                    WORK(KOFF3),NTOTAL,
     &                    BDENSIJ(KOFF2),NTOTI,
     &                    ZERO,WORK(KSCR2),NTOTAL)

C                                                      
C               C(j,alpha) * ( alpha, i | .. ) => ( j, i | ..)
                        CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR2),NTOTAL,
     &                    ZERO,WORK(KSCR3),NTOTJ)

C                        
C              (j i | a delta) * b_{delta j} => U(a,i)
C
                        KOFF5  = IT1AM(ISYMA,ISYMI)
                        KOFF6  = IT1AO(ISDEL,ISYMJ)
C
                        DO I = 1, NRHF(ISYMI)
                          DO J = 1, NRHF(ISYMJ)
                              RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) = 
     &                        RES1(KOFF5+(I-1)*NVIR(ISYMA) + IA) + 
     &                        WORK(KSCR3 - 1 + (I-1)*NRHF(ISYMJ) + J)
     &                        * OBTRIAL(KOFF6 +(J-1)*NBAS(ISDEL)+IDEL)
                          END DO
                        END DO
C
C  162       CONTINUE
C
  161    CONTINUE
C
  160 CONTINUE

      END IF
      

C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_2')
C
      RETURN
      END
