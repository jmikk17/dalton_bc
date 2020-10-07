C
C  /* Deck rp_lanczos_osc_str */
      SUBROUTINE RP_LANCZOS_OSC_STR(EIGVAL,EIGVALC,LANC_CHAIN_MAX,
     &                             TRS_STR_LAN,OSC_STR_LAN)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Feb 2020
C
C     PURPOSE: Given Lanczos transition strengths and eigenvalues 
C              (x, y, or z-component), compute the oscillator strengths.
C              Both arrays are sorted in ascending order of eigenvalues.
C              Eigenvalues computed with LAPACK dgeev.
C              Don't use if complex eigenvalue pairs present! todo fix this
C          
C     EIGVAL            (Re) Lanczos eigenvalues, length=(2*LANC_CHAIN_MAX)
C     EIGVALC           (Im) Lanczos eigenvalues, length=(2*LANC_CHAIN_MAX)
C     LANC_CHAIN_MAX    Length of Lanczos chain       
C     TRS_STR_LAN       Array of transition strengths, length=LANC_CHAIN_MAX
C     OSC_STR_LAN       Array of oscillator strengths, length=LANC_CHAIN_MAX
C
C The eq. used for i^th osc. strength (here for a given component):
C   ---------------------------------
C  |  f_i0 = 2/3 * w_i * ||mu_i||^2  |        
C   ---------------------------------         
C ||mu_i||^2 is the transition strength
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
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D+00, TWO = 2.0D+00, 
     &              THREE = 3.0D+00
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C     
      INTEGER,INTENT(IN) :: LANC_CHAIN_MAX
      REAL(SOP_DP),INTENT(IN),DIMENSION(2*LANC_CHAIN_MAX) :: EIGVAL,
     &                                                       EIGVALC
      REAL(SOP_DP),INTENT(IN),DIMENSION(LANC_CHAIN_MAX)  :: TRS_STR_LAN
      REAL(SOP_DP),INTENT(OUT),DIMENSION(LANC_CHAIN_MAX) :: OSC_STR_LAN
C
C---------------------
C     Local variables
C---------------------
C
      INTEGER :: IOSC
      REAL(SOP_DP) :: OSC_FACTOR
C
C      
C-----------------------------------------------------
C     Add to trace. Assign values to scalar variables.
C-----------------------------------------------------
C
      CALL QENTER('RP_LANCZOS_OSC_STR')
C
      IOSC = 0
      OSC_FACTOR = TWO/THREE
C
C--------------------------------------------------
C     Initialize the osc. strengths array to zeros.
C--------------------------------------------------
C
      CALL DZERO(OSC_STR_LAN,LANC_CHAIN_MAX)
C
C------------------------------------------------------------
C     Loop over Lanczos trans. moments and eigenvalues.
C     (EIGVAL offset by LANC_CHAIN_MAX since we only need 
C     excitation eigenvalues) 
C     Compute oscillator strengths.
C------------------------------------------------------------
C
      DO J = 1,LANC_CHAIN_MAX
        OSC_STR_LAN(J) = 
     &                 EIGVAL(LANC_CHAIN_MAX+J) * TRS_STR_LAN(J)

      END DO
C
      CALL DSCAL(LANC_CHAIN_MAX,OSC_FACTOR,OSC_STR_LAN,1)
C      
C-------------------------------------
C     Print the oscillator strengths.
C-------------------------------------
C      
      IF ( IPRSOP .GE. 0 ) THEN 
C
         CALL AROUND('Oscillator strengths in RP_LANCZOS_OSC_STR')
C
         WRITE(LUPRI,'(I8,1X,F13.6)')
     &        (I,OSC_STR_LAN(I),I=1,LANC_CHAIN_MAX)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_OSC_STR')
C
      RETURN
C
C
      END
