C  /* Deck so_res_a3_1 */
      SUBROUTINE SO_RES_A3_2(RES1,LRES1,BDENSIJ,LBDENSIJ,BDENSAB,
     &                      LBDENSAB,XINT,CMO,LCMO,IDEL,ISDEL,
     &                      ISYDIS,ISYMTR,BTRIAL,LBTRIAL,WORK,LWORK)
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

#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (QUARTER = 0.25D0)
C
      DIMENSION RES1(LRES1), XINT(*)
      DIMENSION BDENSIJ(LBDENSIJ), BDENSAB(LBDENSAB)
      DIMENSION CMO(LCMO)
      DIMENSION WORK(LWORK)
      DIMENSION BTRIAL(LBTRIAL)
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

C     Term: -1/4 del_ij sum_dl (d l | l b ) rho_da  (line 3 second term)
C
      ISYMD = ISDEL
C      
      LCDB  = NVIR(ISYMD)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | l delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSRHF(ISYMD)
C
      KDSRHF  = 1
      KEND2  = KDSRHF + LDSRHF
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.80',LWORK2)
      IF (LWORK2 .LT. 0) 
     &    CALL STOPIT('SO_RES_A3_2.80',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT(XINT,WORK(KDSRHF),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)

C-----------------------------------------------------------------------
C                                  
C           Start loop over l
C-----------------------------------------------------------------------
      DO 080 ISYML = 1,NSYM
       
            ISALBE = MULD2H(ISYML,ISYDIS)
            ISALFA = ISYML
            ISBETA = MULD2H(ISALFA, ISALBE)
            ISYMI  = MULD2H(ISYMTR, ISBETA)
            ISYMA  = MULD2H(ISYMI, ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_1.81',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_1.81',' ',KEND3,LWORK)
C
            DO 081 L = 1,NRHF(ISYML)
C
                  KOFF1 = IDSRHF(ISALBE,ISYML) + NNBST(ISALBE) 
     &                  * (L - 1) + KDSRHF
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | l delta ) for given l 
C           and delta.
C-----------------------------------------------------------------------
C
                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           ~
C              Generate ( l i | l delta ) for
C              given l and delta in KSCR2.
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
                  CALL SO_MEMMAX ('SO_RES_A3_2.82',LWORK3)
                  IF (LWORK3 .LT. 0)
     &              CALL STOPIT('SO_RES_A3_2.82',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
C
                  KOFF2  = ILMRHF(ISYML) + (L-1)*NBAS(ISYML)+1
                  KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF4  = IT1AO(ISYMA,ISYMI) + 1
C
C              C(alpha, l) *( alpha, beta |..)  => ( l, beta | ..)
                  CALL DGEMM('T','N',1,NBAS(ISBETA),
     &                    NBAS(ISALFA),-QUARTER,
     &                    CMO(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),1)

C                                                   ~
C              ( l, beta | .. ) * b(beta,i) => ( l, i | ..)
                  CALL DGEMM('N','N',1,NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),1,
     &                    BTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),1)

C                 ~
C              (l,i | l,delta) * rho_{a delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFFBDENS1 = ILMVIR(ISYMA) - ILMVIR(1)
C
                  COUNT = 0
                  DO A = 1, NVIR(ISYMA)
                        DO I = 1,NRHF(ISYMI)
                              RES1(KOFF5+COUNT+A) = RES1(KOFF5+COUNT+A)
     &                        + WORK(KSCR3 - 1 + I) *
     &                        BDENSAB(KOFFBDENS1 + A)
                        END DO
                        COUNT = COUNT + NRHF(ISYMI)
                  END DO
C
  081    CONTINUE
C
  080 CONTINUE


C     Term: 1/4 del_ij sum_dc (d c | c b ) rho_da  (line 3 third term)
C
      ISYMD = ISDEL
C      
      LCDB  = NVIR(ISYMD)
      IF (LCDB.LE.0) GOTO 999
C
C-----------------------------------------------------------------------
C                                  
C           Transform (alpha beta | gamma delta) to  
C           (alpha beta | c delta)
C-----------------------------------------------------------------------
C
      LDSRHF = NDSVIR(ISYMD)
C
      KDSVIR  = 1
      KEND2  = KDSVIR + LDSVIR
      LWORK2 = LWORK - KEND2
C
      CALL SO_MEMMAX ('SO_RES_A3_2.90',LWORK2)
      IF (LWORK2 .LT. 0) 
     &      CALL STOPIT('SO_RES_A3_2.90',' ',KEND2,LWORK)
C
      ISYMLP = 1
      CALL CCTRBT_VIR(XINT,WORK(KDSVIR),CMO,
     &            ISYMLP,WORK(KEND2),LWORK2,ISYDIS)
C-----------------------------------------------------------------------
C                                  
C           Start loop over c
C-----------------------------------------------------------------------            
      DO 090 ISYMC = 1,NSYM
       
            ISALBE = MULD2H(ISYMC,ISYDIS)
            ISALFA = ISYMC
            ISBETA = MULD2H(ISALFA, ISALBE)
            ISYMI  = MULD2H(ISYMTR, ISBETA)
            ISYMA  = MULD2H(ISYMI, ISYMTR)

            LSCR1  = N2BST(ISALBE)
            KSCR1  = KEND2
            KEND3  = KSCR1 + LSCR1
            LWORK3 = LWORK - KEND3

            CALL SO_MEMMAX ('SO_RES_A3_1.91',LWORK3)
            IF (LWORK3 .LT. 0)
     &          CALL STOPIT('SO_RES_A3_1.91',' ',KEND3,LWORK)
C
            DO 091 C = 1,NVIR(ISYMC)
C
                 KOFF1 = IDSVIR(ISALBE,ISYML) + NNBST(ISALBE) 
     &                  * (L - 1) + KDSVIR
C
C-----------------------------------------------------------------------
C                                  
C           Get a squared set of ( alpha beta | c delta ) for given c 
C           and delta.
C-----------------------------------------------------------------------
C
                  CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR1))
C
C-----------------------------------------------------------------
C                           ~
C              Generate ( c i | c delta ) for
C              given c and delta in KSCR2.
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
                  CALL SO_MEMMAX ('SO_RES_A3_2.92',LWORK3)
                  IF (LWORK3 .LT. 0)
     &               CALL STOPIT('SO_RES_A3_2.92',' ',KEND4,LWORK)
C
                  NTOTAL = MAX(NBAS(ISALFA),1)
                  NTOTBE = MAX(NBAS(ISBETA),1)
C
                  KOFF2  = ILMVIR(ISYMC) + (C-1)*NBAS(ISYMC)+1
                  KOFF3  = KSCR1 + IAODIS(ISALFA,ISBETA)
                  KOFF4  = IT1AO(ISYMA,ISYMI) + 1
C
C              C(alpha, c) *( alpha, beta |..)  => ( c, beta | ..)
                  CALL DGEMM('T','N',1,NBAS(ISBETA),
     &                    NBAS(ISALFA),-QUARTER,
     &                    CMO(KOFF2),NTOTAL,
     &                    WORK(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR2),1)

C                                                   ~
C              ( c, beta | .. ) * b(beta,i) => ( c, i | ..)
                  CALL DGEMM('N','N',1,NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,
     &                    WORK(KSCR2),1,
     &                    BTRIAL(KOFF4),NTOTBE,
     &                    ZERO,WORK(KSCR3),1)

C                 ~
C              (c,i | c,delta) * rho_{a delta} => U(a,i)
C
                  KOFF5  = IT1AM(ISYMA,ISYMI)
                  KOFFBDENS1 = ILMVIR(ISYMA) - ILMVIR(1)
C
                  COUNT = 0
                  DO A = 1, NVIR(ISYMA)
                        DO I = 1,NRHF(ISYMI)
                              RES1(KOFF5+COUNT+A) = RES1(KOFF5+COUNT+A)
     &                        + WORK(KSCR3 - 1 + I) *
     &                        BDENSAB(KOFFBDENS1 + A)
                        END DO
                        COUNT = COUNT + NRHF(ISYMI)
                  END DO
C
  091    CONTINUE
C
  090 CONTINUE

C-----------------------
C     Remove from trace.
C-----------------------
C
  999 CALL QEXIT('SO_RES_A3_1')
C
      RETURN
      END