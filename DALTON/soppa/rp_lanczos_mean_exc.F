C
C  /* Deck rp_lanczos_mean_exc */
      SUBROUTINE RP_LANCZOS_MEAN_EXC(EIGVAL,EIGVALC,LANC_CHAIN_MAX,
     &                              S_SUM,L_SUM,I_SUM,OSC_STR_LAN)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Feb 2020
C
C     PURPOSE: Given Lanczos oscillator strengths  and eigenvalues 
C              (x, y or z-component), compute the S_0 and L_0 sums (1 comp.).
C              Both arrays are sorted in ascending order of eigenvalues.
C              Eigenvalues computed with LAPACK dgeev.
C              Don't use if complex eigenvalue pairs present! todo fix this
C          
C     EIGVAL            (Re) Lanczos eigenvalues, length=(2*LANC_CHAIN_MAX)
C     EIGVALC           (Im) Lanczos eigenvalues, length=(2*LANC_CHAIN_MAX)
C     LANC_CHAIN_MAX    Length of Lanczos chain       
C     OSC_STR_LAN       Array of oscillator strengths, length=LANC_CHAIN_MAX
C     S_SUM             Sum of oscillator strengths (scalar).
C     L_SUM             Sum of ln(w) weighted osc. strengths (scalar).
C     I_SUM             (1 component of) Mean exc. energy, I(0) (scalar).
C
C See Eq.(7.95) in Molecular Electromagnetism by S.P.A. Sauer
C Here calculated for 1 component, fx x:
C   ------------------------------------------------------------------
C  |  S^xx(0) = sum_n(f^xx_n0)  L^xx(0) = sum_n(ln(w^xx_n) * f^xx_n0) |
C   ------------------------------------------------------------------
C  |  ln(I^xx(0)) = L^xx(0) / S^xx(0)                                 |                   
C   ------------------------------------------------------------------
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
C codata.h has the constant, XTEV, to convert from Hartree to eV      
#include "codata.h"
C
C----------------
C     Parameters
C----------------
C      
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D+00, THREE = 3.0D+00  
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      INTEGER,INTENT(IN) ::  LANC_CHAIN_MAX
      REAL(SOP_DP),INTENT(IN),DIMENSION(2*LANC_CHAIN_MAX) :: EIGVAL,
     &                                                       EIGVALC
      REAL(SOP_DP),INTENT(IN),DIMENSION(LANC_CHAIN_MAX) :: OSC_STR_LAN
      REAL(SOP_DP) :: S_SUM, L_SUM, I_SUM
C
C---------------------
C     Local variables
C---------------------
C
      REAL(SOP_DP) :: LNW
C
C-----------------------------------------------------
C     Add to trace. Assign values to scalar variables.
C-----------------------------------------------------
C
      CALL QENTER('RP_LANCZOS_MEAN_EXC')
C
      S_SUM = ZERO
      L_SUM = ZERO
      I_SUM = ZERO
C
C------------------------------------------------------------
C     Compute the sum of Lanczos osc. strengths.
C------------------------------------------------------------
C
      DO I = 1,LANC_CHAIN_MAX
        S_SUM = S_SUM + OSC_STR_LAN(I)
      END DO
C
C----------------------------------------------------------------------
C     Loop over Lanczos osc. strengths and eigenvalues.
C     Compute the sum of ln(w) weighted osc. strengths.
C     Scaling by 3 is due to the 1/3 factor present in osc. strengths.      
C----------------------------------------------------------------------
C
      DO J = 1,LANC_CHAIN_MAX
        LNW = DLOG(EIGVAL(LANC_CHAIN_MAX+J))
        L_SUM = L_SUM + ( LNW * OSC_STR_LAN(J) )
      END DO
C
      L_SUM = L_SUM*THREE
      S_SUM = S_SUM*THREE
      I_SUM = DEXP(L_SUM/S_SUM) * XTEV
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_MEAN_EXC')
C
      RETURN
C
C
      END
