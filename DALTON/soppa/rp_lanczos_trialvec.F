C
C  /* Deck rp_lanczos_trialvec */
      SUBROUTINE RP_LANCZOS_TRIALVEC(GPNORM,
     &                    LABEL,ISYMTR,IMAGPROP,MODEL,
     &                    T2MP,LT2MP,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    DENSAI,LDENSAI,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, 2019-20

C     PURPOSE: Make the initial trial vector for solving the linear
C              response eigenvalue problem using Lanczos RPA solver. The
C              iterations are initiated by a vector consisting of the upper part
C              of the gradient propery vector, the lower part is set to zeros.
C              The norm is saved as it is later used to compute the transition
C              strengths in Lanczos basis.
C     
C     GPNORM        The norm of the trial vector
C     LABEL         Label of the property to be computed
C     ISYMTR        Symmetry of the trial vector
C     IMAGPROP      True if the property is imaginary                  
C     MODEL         Model, in this case RPA only
C     T2MP          Array of the T2 amplitudes
C     LT2MP         Length of the T2 amplitudes array              
C     DENSIJ        Occ.-occ. part of MP2 density
C     LDENSIJ       Length of the occ.-occ. part of MP2 density
C     DENSAB        Virt.-virt. part of MP2 density          
C     LDENSAB       Length of virt.-virt. part of MP2 density                     
C     DENSAI        Virt.-occ. part of MP2 density
C     LDENSAI       Lenght of virt.-occ. part of MP2 density 
C                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C                   
          use so_info, only: sop_dp
C
      IMPLICIT NONE
#include "priunit.h"
#include "soppinf.h"
#include "ccsdsym.h"
C
C-----------------
C    Parameters
C-----------------
C          
      REAL(SOP_DP), PARAMETER :: ZERO = 0.0D+00, ONE = 1.0D+00, EPSI = 1.0D-8
C
C-----------------------------
C    Formal (dummy) arguments
C-----------------------------
C          
      LOGICAL, INTENT(OUT) :: IMAGPROP
C
      CHARACTER(LEN=5), INTENT(IN) :: MODEL
      CHARACTER(LEN=8), INTENT(IN) :: LABEL
C
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LDENSIJ) :: DENSIJ
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LDENSAB) :: DENSAB
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LDENSAI) :: DENSAI
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LT2MP) :: T2MP
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LWORK) :: WORK
C
      INTEGER, INTENT(IN) :: ISYMTR, LDENSIJ, LDENSAB, LDENSAI
      INTEGER, INTENT(IN) :: LT2MP, LWORK
C
      REAL(SOP_DP), INTENT(OUT) :: GPNORM
C
C---------------------
C    Local variables
C---------------------
C          
      INTEGER :: LGPVEC, LGPVEC2, KGPVEC, KGPVEC2, KEND1, LWORK1
C
C---------------------
C    BLAS functions
C---------------------
C          
      REAL(SOP_DP) :: dnrm2
C          
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LANCZOS_TRIALVEC')
C
C--------------------------------
C     Allocation of work space. 
C--------------------------------
C
      LGPVEC = NT1AM(ISYMTR)
C     change this length accordingly if needed
      LGPVEC2 = 0
C
      KGPVEC  = 1
      KGPVEC2  = KGPVEC + LGPVEC
      KEND1   = KGPVEC2 + LGPVEC2
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_LANCZOS_TRIALVEC',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('RP_LANCZOS_TRIALVEC.1',' '
     &                              ,KEND1,LWORK)
C
C----------------------------------
C     Open files for trial vector.
C----------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LGPVEC)
      CALL SO_OPEN(LUTR1D,FNTR1D,LGPVEC)
C
C------------------------------------------
C     Compute the gradient property vector.
C------------------------------------------
C
      CALL SO_GETGP(WORK(KGPVEC),LGPVEC,WORK(KGPVEC2),LGPVEC2,
     &                    LABEL,ISYMTR,IMAGPROP,MODEL,
     &                    T2MP,LT2MP,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    DENSAI,LDENSAI,WORK(KEND1),LWORK1)

C
C-----------------------------------------------------------------
C     Normalize the upper part of the GP vector and save the norm.
C-----------------------------------------------------------------
C
      GPNORM = DNRM2(LGPVEC, WORK(KGPVEC), 1)
      CALL DSCAL(LGPVEC, ONE/GPNORM, WORK(KGPVEC), 1)
C      WRITE(LUPRI,'(A,F16.10)')'Norm of the gradient property vector = '
C     &                         ,GPNORM
C
C-----------------------------------------------------------------------
C    Write the upper part of the GP vector as the upper/excitation part 
C    of the trial vector.
C-----------------------------------------------------------------------
C     
      CALL SO_WRITE(WORK(KGPVEC),LGPVEC,LUTR1E,FNTR1E,1)
C
C-------------------------------------------------------------------
C    Write zeros as the lower/deexcitation part of the trial vector. 
C-------------------------------------------------------------------
C
      CALL DZERO(WORK(KGPVEC),LGPVEC)
C
      CALL SO_WRITE(WORK(KGPVEC),LGPVEC,LUTR1D,FNTR1D,1)
C
C-----------------------------------------------------
C    Write the othonormalized trial vector to output.
C-----------------------------------------------------
C
      IF ( IPRSOP .GE. 5 ) THEN
C
         WRITE(LUPRI,9001)
         WRITE(LUPRI,'(7X,A)') 'The orthonormalized trial vector'
         WRITE(LUPRI,9001)
C
         CALL SO_READ(WORK(KGPVEC),LGPVEC,LUTR1E,FNTR1E,1)
         WRITE (LUPRI,'(I8,1X,F16.10)')
     &          (I,WORK(KGPVEC+I-1),I=1,LGPVEC)
         CALL SO_READ(WORK(KGPVEC),LGPVEC,LUTR1D,FNTR1D,1)
         WRITE (LUPRI,'(I8,1X,F16.10)')
     &          (I,WORK(KGPVEC+I-1),I=1,LGPVEC)
C
      END IF
C
C--------------------------------------------------
C     Close the files containing the trial vector.
C--------------------------------------------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_LANCZOS_TRIALVEC')
C
      RETURN
C
 9001 FORMAT(1X,'---------------------------------------------',
     &       '----------------')
C
      END
