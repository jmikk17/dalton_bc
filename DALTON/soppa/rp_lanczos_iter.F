C
C  /* Deck rp_lanczos_iter */
      SUBROUTINE RP_LANCZOS_ITER(NOLDTR,ISYMTR,
     &                            e_diag, d_diag, a_offdiag, b_offdiag,
     &                            LANC_CHAIN_MAX,BREAKDOWN,     
     &                            WORK,LWORK)
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Nov 2019
C
C     PURPOSE: Execute one Lanczos iteration for the RPA eigenvalue problem. 
C              Read the (ERES transformed) result vectors, then two previous
C              Lanczos vectors and execute the Lanczos recursion.
C              Normalize the new Lanczos vecor and write to file.
C              Save the e, d, a, and b elements of the Lanczos block-tridiagonal
C              matrix.
C
C     NOLDTR            # of previous Lanczos vectors
C     LANC_CHAIN_MAX    # of Lanczos iterations to be executed (k_max)
C     BREAKDOWN         Set to true if a serious break-down is encountered
C     ISYMTR            Symmetry of the property
C     e_diag            diagonal elements in the diagonal blocks of the
C                       block-tridiag. matrix T (diag. for A' in reordered T)
C     d_diag            off-diagonal elements in the diagonal blocks of the
C                       block-tridiag. matrix T (diag. for B' in reordered T)
C     a_offdiag         diagonal elements in the off-diagonal blocks of the 
C                       block-tridiag. matrix T (sub- and super-diag. for A' in
C                       reordered T) 
C     b_offdiag         off-diagonal elements in the off-diagonal blocks of 
C                       the block-tridiag. matrix T (sub- and super-diag. for B'
C                       in reordered T)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
          use so_info, only: sop_dp, sop_lanc_chain_len
C          
      IMPLICIT NONE
C      
C#include "implicit.h"
#include "priunit.h"
#include "soppinf.h"
#include "ccsdsym.h"
C
C----------------
C     Parameters
C----------------
C      
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D+00, ONE = 1.0D0
      REAL(SOP_DP), PARAMETER :: THRLAN = 1.0D-12
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      REAL(SOP_DP), INTENT(INOUT) :: WORK(LWORK)
      INTEGER, INTENT(IN) :: NOLDTR, ISYMTR, LWORK
      INTEGER, INTENT(INOUT) :: LANC_CHAIN_MAX
      REAL(SOP_DP), INTENT(INOUT), DIMENSION(LANC_CHAIN_MAX) :: 
     &                                               e_diag, d_diag 
      REAL(SOP_DP), INTENT(INOUT), DIMENSION(LANC_CHAIN_MAX-1) :: 
     &                                         a_offdiag, b_offdiag 
      LOGICAL, INTENT(INOUT) :: BREAKDOWN
C      
C---------------------
C     Local variables
C---------------------
C
      REAL(SOP_DP) :: SqNorm, abs_Norm, Norm, InvNorm 
      INTEGER :: LLAN1E
      INTEGER :: KIS, KRES1E, KRES1D, KLAN1E, KLAN1D, KEND1, LWORK1
C
C--------------------
C     BLAS routines
C--------------------
C
      REAL(SOP_DP) :: DDOT
C
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LANCZOS_ITER')
C
C
C--------------------------------------------------------------------
C     Work space allocation.
C     LLAN1E : length of the excitation part of Lanczos vector, same 
C     for the deexcitation part.
C--------------------------------------------------------------------
C
      LLAN1E    = NT1AM(ISYMTR)
C
      KIS     = 1
      KRES1E  = KIS + LLAN1E
      KRES1D  = KRES1E + LLAN1E
      KLAN1E  = KRES1D + LLAN1E
      KLAN1D  = KLAN1E + LLAN1E
      KEND1   = KLAN1D + LLAN1E
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_LANCZOS_ITER.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_LANCZOS_ITER.1',' ',
     &                                KEND1,LWORK)
C
C
C===================================================================
C    Execute one Lanczos iteration to create k+1^th Lanczos vector. 
C===================================================================
C
C------------------------------------------------------------------------
C     Open the files with result (E'-transformed) and Lanczos 
C     vectors.
C     Read the (SO_ERES) transformed result vectors.
C     Read the Lanczos vectors from the previous iteration (NOLDTR+1=k).
C------------------------------------------------------------------------
C
      CALL SO_OPEN(LURS1E,FNRS1E,LLAN1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LLAN1E)
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LLAN1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LLAN1E)
C
      CALL SO_READ(WORK(KRES1E), LLAN1E, LURS1E, FNRS1E, NOLDTR+1)
      CALL SO_READ(WORK(KRES1D), LLAN1E, LURS1D, FNRS1D, NOLDTR+1)
C
      CALL SO_READ(WORK(KLAN1E), LLAN1E, LUTR1E, FNTR1E, NOLDTR+1)
      CALL SO_READ(WORK(KLAN1D), LLAN1E, LUTR1D, FNTR1D, NOLDTR+1)
C
C---------------------------------------------------------
C     Scale the deexcitation part of the result vector. 
C---------------------------------------------------------
C
      CALL dscal(LLAN1E, -ONE, WORK(KRES1D), 1)
C
C--------------------------------------------------------------
C     Compute the new diagonal elements, e_k and d_k, of the 
C     block-tridiagonal matrix blocks.
C     Quit if the k_max (desired nr of iterations) is reached.      
C--------------------------------------------------------------
C
      e_diag(NOLDTR+1) = ddot(LLAN1E, WORK(KRES1E), 1, WORK(KLAN1E), 1)
     &                 - ddot(LLAN1E, WORK(KRES1D), 1, WORK(KLAN1D), 1)
      d_diag(NOLDTR+1) = ddot(LLAN1E, WORK(KRES1E), 1, WORK(KLAN1D), 1)
     &                 - ddot(LLAN1E, WORK(KRES1D), 1, WORK(KLAN1E), 1)
C
      IF (NOLDTR+1 .EQ. LANC_CHAIN_MAX) THEN
          CALL QEXIT('RP_LANCZOS_ITER')
          RETURN
      END IF
C
C----------------------------------------------------------------
C     Lanczos recursion. Subtract the components of the previous
C     Lanczos vector (X_k, Y_k).
C----------------------------------------------------------------
C
      CALL daxpy(LLAN1E,d_diag(NOLDTR+1),WORK(KLAN1D),1,
     &           WORK(KRES1E),1)
      CALL daxpy(LLAN1E,-e_diag(NOLDTR+1),WORK(KLAN1E),1,
     &           WORK(KRES1E),1)
C
      CALL daxpy(LLAN1E,d_diag(NOLDTR+1),WORK(KLAN1E),1,
     &           WORK(KRES1D),1)
      CALL daxpy(LLAN1E,-e_diag(NOLDTR+1),WORK(KLAN1D),1,
     &           WORK(KRES1D),1)
C
C------------------------------------------------------------------
C     Overwrite the Lanczos vector from the previous iteration (k)
C     with the Lanczos vectors from one iteration down (k-1).
C------------------------------------------------------------------
C
      IF (NOLDTR .NE. 0) THEN
C          
         CALL SO_READ(WORK(KLAN1E),LLAN1E,LUTR1E,FNTR1E,NOLDTR)
         CALL SO_READ(WORK(KLAN1D),LLAN1E,LUTR1D,FNTR1D,NOLDTR)
C
C-----------------------------------------------------------------
C     Lanczos recursion continued. Subtract the components of the 
C     previous Lanczos vector (X_{k-1}, Y_{k-1}).
C-----------------------------------------------------------------
C
          CALL daxpy(LLAN1E,b_offdiag(NOLDTR),WORK(KLAN1D),1,
     &                                         WORK(KRES1E),1)
          CALL daxpy(LLAN1E,-a_offdiag(NOLDTR),WORK(KLAN1E),1,
     &                                         WORK(KRES1E),1)
C
          CALL daxpy(LLAN1E,b_offdiag(NOLDTR),WORK(KLAN1E),1,
     &                                         WORK(KRES1D),1)
          CALL daxpy(LLAN1E,-a_offdiag(NOLDTR),WORK(KLAN1D),1,
     &                                         WORK(KRES1D),1)
C          
      END IF
C
      IF ( IPRSOP .GE. 9 ) THEN
C
         CALL AROUND('New Lanczos vector after daxpy in lanczos_iter')
         WRITE(LUPRI,'(I8,1X,F16.10,5X,F16.10)')
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LLAN1E)
C
      END IF
C
C-------------------------------------------------------------------
C     Bi-orthogonalize the new Lanczos vector against all previous 
C     Lanczos vectors (also normalize). 
C-------------------------------------------------------------------
C
      CALL RP_LANCZOS_BIORTH(NOLDTR,LLAN1E,WORK(KRES1E),WORK(KRES1D),
     &                                SqNorm,Norm,WORK(KEND1),LWORK1)
C
C-----------------------------------------------------------------------------
C     Compute the norm of the new Lanczos vector and check for a serious
C     break-down. 
C     Check if the excitation or de-exc. symmetry vector is found.
C     Save the off-diagonal elemnts a_k, b_k of the block-tridiagonal matrix.
C     Write the new Lanczos vector to file.
C-----------------------------------------------------------------------------
C      
      abs_Norm = dabs(SqNorm)
C      
      IF (abs_Norm .LE. THRLAN) THEN
          WRITE(LUPRI,*) 'WARNING: After RP_LANCZOS_BIORTH, at it.= ',
     &                             NOLDTR+1, ' abs_Norm = ', abs_Norm,
     &    'Invariant subspace or serious break-down reached, EXITING'
          BREAKDOWN = .TRUE.
C
          CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
          CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
C
          CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
          CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
          CALL FLSHFO(LUPRI)
          CALL QEXIT('RP_LANCZOS_ITER')
          RETURN
      END IF
C
      IF (SqNorm .GE. ZERO) THEN
          a_offdiag(NOLDTR+1) = Norm
          b_offdiag(NOLDTR+1) = ZERO
C
          CALL SO_WRITE(WORK(KRES1E),LLAN1E,LUTR1E,FNTR1E,NOLDTR+2)
          CALL SO_WRITE(WORK(KRES1D),LLAN1E,LUTR1D,FNTR1D,NOLDTR+2)
C          
      ELSE
          b_offdiag(NOLDTR+1) = Norm
          a_offdiag(NOLDTR+1) = ZERO
C
          CALL dscal(LLAN1E, -ONE, WORK(KRES1E), 1)
          CALL dscal(LLAN1E, -ONE, WORK(KRES1D), 1)
C          
          CALL SO_WRITE(WORK(KRES1D),LLAN1E,LUTR1E,FNTR1E,NOLDTR+2)
          CALL SO_WRITE(WORK(KRES1E),LLAN1E,LUTR1D,FNTR1D,NOLDTR+2)
C          
      END IF
C      
      IF ( IPRSOP .GE. 7 ) THEN
C
         CALL AROUND('New normalized Lanczos vector')
         WRITE(LUPRI,'(I8,1X,F16.10,5X,F16.10)')
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LLAN1E)
      END IF
C
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C---------------------------------
C     Remove from trace and exit.
C---------------------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_LANCZOS_ITER')
C
      RETURN
      END
