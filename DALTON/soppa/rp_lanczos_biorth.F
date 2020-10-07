C
C  /* Deck rp_lanczos_biorth */
      SUBROUTINE RP_LANCZOS_BIORTH(NOLDTR,LNEWLAN,NEWLAN1E,
     &                       NEWLAN1D,SQNORM,NORM,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Nov 2019
C
C     PURPOSE: Bi-orthogonalize the new Lanczos vector against all previous
C              Lanczos vectors (both left/right, excitation/de-excitation) 
C              and normalize. Return the norm and the orthonormalized
C              new Lanczos vector.
C          
C     NOLDTR      NOLDTR+1 is a number of previous Lanczos vectors 
C     LNEWLAN     Length of Lanczos vectors
C     NEWLAN1E    The upper (exc.) part of the new Lanczos vector    
C     NEWLAN1D    The lower (de-exc.) part of the new Lanczos vector 
C     SQNORM      The square of the norm of the new Lanczos vector    
C     NORM        The norm of the new Lanczos vector  
C          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          
         use so_info, only: sop_dp, sop_lanc_chain_len
C          
      IMPLICIT NONE
C            
#include "priunit.h"
C
C#include "ccsdsym.h"
C#include "ccorb.h"
#include "soppinf.h"
C
C----------------
C     Parameters
C----------------
C      
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D+00, ONE = 1.0D+00
      REAL(SOP_DP), PARAMETER :: THRLDP = 1.0D-20, THROUND = 1.0D-4,
     &                           T1MIN = 1.0D-8 
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      REAL(SOP_DP), INTENT(INOUT) :: WORK(LWORK)
      INTEGER, INTENT(IN) :: LWORK, NOLDTR, LNEWLAN
      REAL(SOP_DP), INTENT(INOUT), DIMENSION(LNEWLAN) :: NEWLAN1E,
     &                                                   NEWLAN1D 
      REAL(SOP_DP), INTENT(OUT) :: SQNORM, NORM
C
C---------------------
C     Local variables
C---------------------
C
      INTEGER :: KLAN1E, KLAN1D, LWORK1, KEND1
      INTEGER :: ITURN, IPRVTR, I
      REAL(SOP_DP) :: DOTP, INVNORM, OLD_NORM, NORM_RATIO
C
C-------------------------
C     BLAS functions used
C-------------------------
C
      REAL(SOP_DP) :: DDOT 
C      
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LANCZOS_BIORTH')
C
C-----------------------------------------------------------
C     Allocation of work space for previous Lanczos vectors.
C-----------------------------------------------------------
C
      KLAN1E   = 1
      KLAN1D   = KLAN1E  + LNEWLAN
C
      KEND1   = KLAN1D + LNEWLAN
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_LANCZOS_BIORTH',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_LANCZOS_BIORTH',' ',
     &                                           KEND1,LWORK)
C
C-------------------------------------------------------------
C     Check the norm of the new Lanczos vector.
C     Print the new raw Lanczos vector.
C     Check for serious break-down/exhausted invariant space.
C-------------------------------------------------------------
C
      OLD_NORM = DDOT(LNEWLAN,NEWLAN1E,1,NEWLAN1E,1)
     &         - DDOT(LNEWLAN,NEWLAN1D,1,NEWLAN1D,1)
C      
      IF ( IPRSOP .GE. 7 ) THEN
C
         CALL AROUND('Raw new Lanczos vector in RP_LANCZOS_BIORTH')
C
         WRITE(LUPRI,'(I8,1X,F16.10,5X,F16.10)')
     &           (I,NEWLAN1E(I),NEWLAN1D(I),I=1,LNEWLAN)
C
      END IF
C
      IF (DABS(OLD_NORM) .LE. T1MIN) THEN
          WRITE(LUPRI,'(A,2X,I8,2X,A,D29.15)')
     &   'WARNING: Norm of the new Lanczos vector  at it. = ',NOLDTR+1,
     &   'is close to zero (1.0D-08), Square Norm = ', OLD_NORM
         WRITE(LUPRI,'(A)') 'Might be a serious break-down' 
      END IF
C
C
C--------------------------------------------
C     Bi-orthogonalize twice.
C         Loop over previous Lanczos vectors.
C--------------------------------------------
C
      DO ITURN = 1, 2
          DO IPRVTR = 1,NOLDTR+1
C
C---------------------------------------
C            Read previous Lanczos vectors.
C---------------------------------------
C
             CALL SO_READ(WORK(KLAN1E),LNEWLAN,LUTR1E,FNTR1E,IPRVTR)
             CALL SO_READ(WORK(KLAN1D),LNEWLAN,LUTR1D,FNTR1D,IPRVTR)
C
             IF ( IPRSOP .GE. 9 ) THEN
C
                CALL AROUND('Previous trial vector in
     &                            RP_LANCZOS_BIORTH')
C
                WRITE(LUPRI,'(I8,1X,F16.10,5X,F16.10)')
     &               (I,WORK(KLAN1E+I-1),WORK(KLAN1D+I-1),I=1,LNEWLAN)
C
             END IF
C
C------------------------------------------------------------------------
C            Orthogonalize new Lanczos vector against previous Lanczos
C            vectors and their de-excitation/paired partners.
C------------------------------------------------------------------------
C

             DOTP = DDOT(LNEWLAN,NEWLAN1E,1,WORK(KLAN1E),1)
     &            - DDOT(LNEWLAN,NEWLAN1D,1,WORK(KLAN1D),1)
C
             CALL DAXPY(LNEWLAN,-DOTP,WORK(KLAN1E),1,NEWLAN1E,1)
             CALL DAXPY(LNEWLAN,-DOTP,WORK(KLAN1D),1,NEWLAN1D,1)
C
C
             DOTP = DDOT(LNEWLAN,NEWLAN1E,1,WORK(KLAN1D),1)
     &            - DDOT(LNEWLAN,NEWLAN1D,1,WORK(KLAN1E),1)
C
             CALL DAXPY(LNEWLAN,DOTP,WORK(KLAN1D),1,NEWLAN1E,1)
             CALL DAXPY(LNEWLAN,DOTP,WORK(KLAN1E),1,NEWLAN1D,1)
C
          END DO
C
C---------------------------------------------------------------
C        Check the norm of the new Lanczos vector after
C        bi-orthogonalization. Determine if needs to be bi-orth.
C        again.
C---------------------------------------------------------------
C
          SQNORM = DDOT(LNEWLAN,NEWLAN1E,1,NEWLAN1E,1)
     &           - DDOT(LNEWLAN,NEWLAN1D,1,NEWLAN1D,1)
C          
          NORM_RATIO = SQNORM/OLD_NORM
          IF (DABS(NORM_RATIO) .GT. THROUND) EXIT
C
          OLD_NORM = SQNORM
C
      END DO
C
C---------------------------------------
C     Normalize the new Lanczos vector.
C     Print the new Lanczos vector.
C---------------------------------------
C
      NORM = DSQRT(DABS(SQNORM))
      INVNORM = ONE/NORM
C
      CALL DSCAL(LNEWLAN,INVNORM,NEWLAN1E,1)
      CALL DSCAL(LNEWLAN,INVNORM,NEWLAN1D,1)
C
C
      IF ( IPRSOP .GE. 5 ) THEN
C
         CALL AROUND('Orthonormalized Lanczos vector in 
     &                               RP_LANCZOS_BIORTH')
C
         WRITE(LUPRI,'(I8,1X,F16.10,5X,F16.10)')
     &        (I,NEWLAN1E(I),NEWLAN1D(I),I=1,LNEWLAN)
         WRITE(LUPRI,'(A,D24.10)') 'Norm of the new Lanczos vector after
     &        RP_LANCZOS_BIORTH', SQNORM
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_BIORTH')
C
      RETURN
C
C
      END
