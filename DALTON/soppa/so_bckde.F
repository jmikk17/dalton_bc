C
C  /* Deck so_bckde */
      SUBROUTINE SO_BCKDE(DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                   BDENSIJ,LBDENSIJ,BDENSAB,LBDENSAB,
     &                   CMO, LCMO)

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
      DIMENSION BDENSIJ(LBDENSIJ), BDENSAB(LBDENSAB)
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
            NTOTAL = MAX(NVIR(ISYMA),1)          
C
            KOFF1  = IABDEN(ISYMA,ISYMA)
            KOFF2  = ILMVIR(ISYMA) + 1
            KOFF3  = ILMVIR(ISYMA) - ILMVIR(1) + 1
C
            CALL DGEMM('N','T',NVIR(ISYMA),NORB(ISYMA),
     &         NVIR(ISYMA),ONE,DENSAB(KOFF1),NTOTAL,
     &         CMO(KOFF2),NORB,ZERO,BDENSAB(KOFF3),NTOTAL)
        END DO
C      
C     Transform occupied index
C     rho_ij to rho_idelta
        DO ISYMJ = 1,NSYM
C
               NTOTAL  = MAX(NRHF(ISYMJ),1)
C
               KOFF1  = IIJDEN(ISYMJ,ISYMJ) + 1
               KOFF2  = ILMRHF(ISYMJ) + 1 
C
C               WRITE(LUPRI,'(I5,A)') NRHF(ISYMJ),'NRHF(ISYMJ)'
C               WRITE(LUPRI,'(I5,A)') NORB(ISYMJ),'NORB(ISYMJ)'
C               WRITE(LUPRI,'(F14.8,A)') KOFF1,'KOFF1'
C               DO k = 1,NRHF(ISYMJ)*NRHF(ISYMJ)
C                  WRITE(LUPRI,'(I5,F14.8,A)')
C     &                 k,DENSIJ(KOFF1+k-1),'DENSIJ(KOFF1)'
C               END DO
C               WRITE(LUPRI,'(F14.8,A)') KOFF2,'KOFF2'
C               WRITE(LUPRI,'(F14.8,A)') CMO(KOFF2),'CMO(KOFF2)'
C               
C
               CALL DGEMM('N','T',NRHF(ISYMJ),NORB(ISYMJ),
     &           NRHF(ISYMJ),ONE,DENSIJ(KOFF1),NRHF(ISYMJ),
     &           CMO(KOFF2),NORB(ISYMJ),ZERO,BDENSIJ(KOFF2),NTOTAL)
c
C               WRITE(LUPRI,'(F14.8,A)') BDENSIJ(KOFF2),'BDENSIJ(KOFF2)'
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
