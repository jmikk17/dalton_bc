C
C  /* Deck rp_lanczos_trs_str */
      SUBROUTINE RP_LANCZOS_TRS_STR(EIGVECR,EIGVECL,LANC_CHAIN_MAX,
     &                             TRS_STR_LAN,GPNORM)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Feb 2020
C
C     PURPOSE: Given the left and right Lanczos eigevectors and the norm of
C              the gradient prop. vector that started the chain, compute
C              the transition strengths. The eigenvectors are
C              from LAPACK dgeev routine, stored as columns, sorted in
C              order of ascending eigenvalues (deexcit.-> excit.).
C              Don't use if complex eigenvalue pairs present! todo fix this
C          
C     EIGVECR           Right Lanczos eigenvectors, dim=(2*LANC_CHAIN_MAX)^2
C     EIGVECL           Left Lanczos eigenvectors, dim=(2*LANC_CHAIN_MAX)^2
C     LANC_CHAIN_MAX    Length of Lanczos chain       
C     TRS_STR_LAN       Array of transition moments, length=LANC_CHAIN_MAX
C     GPNORM            Norm of the gradient prop. vector (start vector)
C          
C The eq. used for i^th trs. strength (for a given component): 
C ---------------------------------------------------------------------
C |  ||mu_i||^2 = ( R_1i - R_(K+1)i ) * ( L1i + L(K+1)i  ) / GPNORM^2  |
C ---------------------------------------------------------------------
C R and L are matrices of the right and left eigenvectors (column sorted)          
C K = LANC_CHAIN_MAX          
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
     &                           (4*LANC_CHAIN_MAX*LANC_CHAIN_MAX) :: 
     &                           EIGVECL, EIGVECR
      REAL(SOP_DP),INTENT(INOUT),DIMENSION(LANC_CHAIN_MAX) :: 
     &                           TRS_STR_LAN
      REAL(SOP_DP),INTENT(IN) :: GPNORM
C
C---------------------
C     Local variables
C---------------------
C
      INTEGER :: OFFSET
      REAL(SOP_DP) :: SQGPNORM
C
C      
C-----------------------------------------------------
C     Add to trace. Assign values to scalar variables.
C-----------------------------------------------------
C
      CALL QENTER('RP_LANCZOS_TRS_STR')
C
      SQGPNORM = GPNORM * GPNORM
C
C--------------------------------------------------
C     Initialize the trs. moments array to zeros.
C--------------------------------------------------
C
      CALL DZERO(TRS_STR_LAN,LANC_CHAIN_MAX)
C
C--------------------------------------------------------------------
C     Loop over Lanczos excitation eigenvectors.
C     Use the 1st and the (K+1)^th (full dim. is (2K)^2 )
C     elements of an eigenvectors to compute trans. strengths. 
C     (EIGVECR and EIGVECL are offset by 2*LANC_CHAIN_MAX*LANC_CHAIN_MAX
C     to begin with since we only need excitation eigenvectors)      
C     Scale the transition moments vector.
C--------------------------------------------------------------------
C
      OFFSET = 2*LANC_CHAIN_MAX*LANC_CHAIN_MAX
      DO J = 1,LANC_CHAIN_MAX
        TRS_STR_LAN(J) = 
     &  ( EIGVECR(OFFSET+1) - EIGVECR(OFFSET+LANC_CHAIN_MAX+1) ) *
     &  ( EIGVECL(OFFSET+1) + EIGVECL(OFFSET+LANC_CHAIN_MAX+1) )
        OFFSET = OFFSET + 2*LANC_CHAIN_MAX
C        
      END DO
C
      CALL DSCAL(LANC_CHAIN_MAX,SQGPNORM,TRS_STR_LAN,1)
C      
C--------------------------------------------------------------
C     Print the transition strengths.
C--------------------------------------------------------------
C      
      IF ( IPRSOP .GE. 2 ) THEN
C
         CALL AROUND('Transition strengths in RP_LANCZOS_TRS_STR')
C
         WRITE(LUPRI,'(I8,1X,F16.10)')
     &        (I,TRS_STR_LAN(I),I=1,LANC_CHAIN_MAX)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_TRS_STR')
C
      RETURN
C
C
      END
