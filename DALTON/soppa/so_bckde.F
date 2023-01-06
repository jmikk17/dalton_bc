C
C  /* Deck so_bckde */
      SUBROUTINE SO_BCKDE(DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                   BDENSIJ,BDENSIJT,LBDENSIJ,BDENSAB,BDENSABT,
     &                   LBDENSAB,CMO,LCMO)

C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Kasper F. Schaltz October 2022
C     
C     PURPOSE: Back transform one index of the one-particle density 
C              matrices
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION BDENSIJ(LBDENSIJ), BDENSIJT(LBDENSIJ)
      DIMENSION BDENSAB(LBDENSAB), BDENSABT(LBDENSAB)
      DIMENSION CMO(LCMO)
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
      CALL QENTER('SO_BCKDE')

C     Transform virtual index
C     rho_ab to rho_adelta
C
      DO ISYMA = 1,NSYM
C
            NTOTAL = MAX(NBAS(ISYMA),1)
            NTOTA  = MAX(NVIR(ISYMA),1)       
C
            KOFF1  = IABDEN(ISYMA,ISYMA) + 1
            KOFF2  = ILMVIR(ISYMA) + 1
            KOFF3  = ILMVIR(ISYMA) - ILMVIR(1) + 1
C
            CALL DGEMM('N','T',NVIR(ISYMA),NBAS(ISYMA),
     &         NVIR(ISYMA),ONE,DENSAB(KOFF1),NTOTA,
     &         CMO(KOFF2),NTOTAL,ZERO,BDENSAB(KOFF3),NTOTA)
C     Calculate transpose of BDENSAB aswell  -  rho_delta a
            CALL DGEMM('N','N',NBAS(ISYMA),NVIR(ISYMA),
     &         NVIR(ISYMA),ONE,CMO(KOFF2),NTOTAL,
     &         DENSAB(KOFF1),NTOTA,ZERO,BDENSABT(KOFF3),NTOTAL)
      END DO
C      
C     Transform occupied index
C     rho_ij to rho_idelta
      DO ISYMJ = 1,NSYM
C
            NTOTAL = MAX(NBAS(ISYMJ),1)
            NTOTJ  = MAX(NRHF(ISYMJ),1) 
C
            KOFF1  = IIJDEN(ISYMJ,ISYMJ) + 1
            KOFF2  = ILMRHF(ISYMJ) + 1 
C
            CALL DGEMM('N','T',NRHF(ISYMJ),NBAS(ISYMJ),
     &           NRHF(ISYMJ),ONE,DENSIJ(KOFF1),NTOTJ,
     &           CMO(KOFF2),NTOTAL,ZERO,BDENSIJ(KOFF2),NTOTJ)
C     Calculate transpose of BDENSIJ aswell  -  rho_delta j
            CALL DGEMM('N','N',NBAS(ISYMJ),NRHF(ISYMJ),
     &           NRHF(ISYMJ),ONE,CMO(KOFF2),NTOTAL,
     &           DENSIJ(KOFF1),NTOTJ,ZERO,BDENSIJ(KOFF2),NTOTAL)
      END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_BCKDE')
C
      RETURN
      END
