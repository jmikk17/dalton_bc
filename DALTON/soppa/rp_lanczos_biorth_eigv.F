C
C  /* Deck rp_lanczos_biorth_eigv */
      SUBROUTINE RP_LANCZOS_BIORTH_EIGV(EIGVECR,EIGVECL,LANC_CHAIN_MAX,
     &                                                  EIGVAL,EIGVALC)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Feb 2020
C
C     PURPOSE: Given left and right Lanczos eigenvectors from a LAPACK 
C              DGEEV output, biorthogonalize them to each other so that
C              <R_i|R_i>=1 and <L_i|R_i>=1, by scaling both of them with 
C              1/sqrt(<L_i|R_i>).
C              TODO: adapt to the case of complex eigenpairs.
C          
C     EIGVECR           Right Lanczos eigenvectors, dim=(2*LANC_CHAIN_MAX)^2
C     EIGVECL           Left Lanczos eigenvectors, dim=(2*LANC_CHAIN_MAX)^2
C     EIGVAL            (Re) eigenvalues, dim=2*LANC_CHAIN_MAX 
C     EIGVALC           (Im) eigenvalues, dim=2*LANC_CHAIN_MAX 
C     LANC_CHAIN_MAX    Length of Lanczos chain       
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          
         use so_info, only: sop_dp, sop_lanc_chain_len
C          
      IMPLICIT NONE
C            
#include "priunit.h"
C
#include "ccsdsym.h"
#include "ccorb.h"
#include "soppinf.h"
C
C----------------
C     Parameters
C----------------
C      
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D+00, ONE = 1.0D+00
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      INTEGER,INTENT(IN) ::      LANC_CHAIN_MAX
      REAL(SOP_DP),INTENT(IN),DIMENSION
     & (2*LANC_CHAIN_MAX,2*LANC_CHAIN_MAX) ::
     &                           EIGVECL, EIGVECR
      REAL(SOP_DP),INTENT(IN),DIMENSION(2*LANC_CHAIN_MAX) ::
     &                           EIGVAL, EIGVALC
C
C---------------------
C     Local variables
C---------------------
C
C      comment out JVEC and DOTP if bi-orthogonalization is to be used
C      INTEGER :: JVEC
C      REAL(SOP_DP) :: DOTP
      REAL(SOP_DP) :: DOTLR, SQRDOTLR, INV
C
C------------------------
C     BLAS functions used
C------------------------
C      
      REAL(SOP_DP) :: DDOT
C      
C-----------------------------------------------------
C     Add to trace.
C-----------------------------------------------------
C
      CALL QENTER('RP_LANCZOS_BIORTH_EIGV')
C
C-------------------------------------------------------------------
C     Loop over Lanczos left and right eigenvectors.
C     Compute dot products <L_i|R_j>, bi-orthogonalize.
C     !!! This is left here in case the diagonalization routine is 
C     changed and one needs to explicitly bi-orthogonalize.
C-------------------------------------------------------------------
C
C      DO J = 1, 2*LANC_CHAIN_MAX-1
C        DO I = J+1, 2*LANC_CHAIN_MAX
C          DOTP = DDOT(2*LANC_CHAIN_MAX,EIGVECL(1,I),1,EIGVECR(1,J),1)
C          CALL DAXPY(2*LANC_CHAIN_MAX,-DOTP,EIGVECR(1,I),1,
C     &               EIGVECR(1,J),1)
C
C        END DO
C      END DO
C
C----------------------------------------------------------------
C     Loop over Lanczos left and right eigenvectors.
C     Compute the dot product <L_i|R_i> and it's square root, 
C     use it to scale both eigenvectors. In case the dot product
C     is negative, the sign is assigned to the left eigenvector.     
C----------------------------------------------------------------
C
      DO J = 1,2*LANC_CHAIN_MAX
        DOTLR = DDOT(2*LANC_CHAIN_MAX,EIGVECL(1,J),1,EIGVECR(1,J),1)
        SQRDOTLR = DSQRT(DABS(DOTLR))
        INV = ONE/SQRDOTLR
C        
C
        IF (DOTLR .GT. ZERO) THEN
            CALL DSCAL(2*LANC_CHAIN_MAX,INV,EIGVECL(1,J),1 )
            CALL DSCAL(2*LANC_CHAIN_MAX,INV,EIGVECR(1,J),1 )
        ELSE
            CALL DSCAL(2*LANC_CHAIN_MAX,-INV,EIGVECL(1,J),1 )
            CALL DSCAL(2*LANC_CHAIN_MAX,INV,EIGVECR(1,J),1 )
        END IF
      END DO
C
C
C---------------------------------------------------------------
C     Write the newly bi-othonormalized eigenvectors to output.
C---------------------------------------------------------------
C      
      IF ( IPRSOP .GE. 5 ) THEN
        CALL AROUND('Bi-orthogonalized Lanczos eigenvectors')
        DO I = 1,2*LANC_CHAIN_MAX
            WRITE (LUPRI,'(8X,A,1X,I8)') 'EIGENVECTOR No.',I
C            WRITE (LUPRI,'(8X,A,1X,I8,30X,A)')'|  EIGENVECTOR No.',I,'|'
            WRITE (LUPRI,'(8X,A,A)') '-----------------------------',
     &      '--------------------'
            WRITE (LUPRI,'(19X,A,16X,A)')'Right','Left'
            WRITE (LUPRI,'(8X,A,A)') '-----------------------------',
     &      '--------------------'
            WRITE (LUPRI,'(6X,I8,1X,F16.10,5X,F16.10)')
     &      (K,EIGVECR(K,I),EIGVECL(K,I)
     &      ,K=1,2*LANC_CHAIN_MAX)
            WRITE (LUPRI,'(8X,A,A)') '-----------------------------',
     &      '--------------------'
            WRITE (LUPRI,*)
            WRITE (LUPRI,*)
        END DO
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_BIORTH_EIGV')
C
      RETURN
C
C
      END

