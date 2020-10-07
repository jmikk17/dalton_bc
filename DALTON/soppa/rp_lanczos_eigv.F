C
C  /* Deck rp_lanczos_eigv */
      SUBROUTINE RP_LANCZOS_EIGV(E_DIAG,D_DIAG,A_OFFDIAG,B_OFFDIAG,
     &                       LLANCHAIN,EIGVAL,EIGVALC,EIGVECR,EIGVECL,
     &                                              T_MATR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Feb 2020
C
C     PURPOSE: Build Johnson et al. matrix T from coefficients computed
C              by Lanczos iterations. Diagonalize it, reorder the left and
C              right eigenvectors. 
C              TODO: introduce handling of complex eigenvalues. In case complex
C              eigenvalues are found, DGEEV returns the complex
C              eigenvalues and eigenvectors as pairs (see LAPACK documentation)
C              and this routine will not handle them correctly.
C          
C     E_DIAG       Diagonal elements in the diagonal blocks A' in T matrix. 
C     D_DIAG       Offdiagonal elements in the diagonal blocks A' in T matrix
C     A_OFFDIAG    Diagonal elements in the offdiagonal blocks B' in T matrix
C     B_OFFDIAG    Offdiagonal elements in the offdiagonal blocks B' in T matrix 
C     LLANCHAIN    The length of Lanczos chain (equal to 1/2 dimesnion of T)      
C     EIGVAL       The real part of eigenvalues      
C     EIGVALC      The imaginary part of eigenvalues
C     EIGVECR      The right eigenvectors
C     EIGVECL      The left eigenvectors
C     T_MATR       The Johnson T matrix
C          
C   ----------------------------------------       
C  |                                       |
C  |   |  A   B |   Lanczos  |  A'  B' |   |
C  |   | -B  -A |   ------>  | -B' -A' |   |
C  |                                       |
C   ----------------------------------------       
C  See more in: https://doi.org/10.1016/S0010-4655(99)00248-9          
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          
         use so_info, only: sop_dp, sop_lanc_chain_len,
     &                      sop_lanc_conv_check         
C          
      IMPLICIT NONE
C            
#include "priunit.h"
C
C ccsdsym has integer i declared
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
C---------------------------------
C     Dimensions of the arguments
C---------------------------------
C
      INTEGER, INTENT(IN) :: LLANCHAIN,LWORK
      REAL(SOP_DP), INTENT(INOUT), DIMENSION(LWORK) :: WORK
      REAL(SOP_DP), INTENT(IN), DIMENSION(LLANCHAIN) :: E_DIAG,D_DIAG
      REAL(SOP_DP), INTENT(IN), DIMENSION(LLANCHAIN-1) :: 
     &                                              A_OFFDIAG,B_OFFDIAG 
      REAL(SOP_DP), INTENT(OUT), DIMENSION(2*LLANCHAIN) :: 
     &                                                    EIGVAL,EIGVALC
      REAL(SOP_DP), INTENT(OUT), DIMENSION(2*LLANCHAIN,2*LLANCHAIN) ::
     &                                            T_MATR
      REAL(SOP_DP), INTENT(OUT), DIMENSION(4*LLANCHAIN*LLANCHAIN) ::
     &                                            EIGVECR,EIGVECL
C
C---------------------
C     Local variables
C---------------------
C
      INTEGER :: DIM_T,DIM_EIGVAL
      INTEGER :: ICOMPLX,IOFFSET,IERR
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LANCZOS_EIGV')
C
C----------------------------------------
C     Set the dimensions for eigenpairs.
C----------------------------------------
C
      DIM_EIGVAL = 2*LLANCHAIN
      DIM_T =      4*LLANCHAIN*LLANCHAIN
C
C--------------------------------------------------
C     Initialize the matrix and vectors to zeros.
C--------------------------------------------------
C
      CALL DZERO(T_MATR,DIM_T)
      CALL DZERO(EIGVAL,DIM_EIGVAL)
      CALL DZERO(EIGVALC,DIM_EIGVAL)
      CALL DZERO(EIGVECR,DIM_T)
      CALL DZERO(EIGVECL,DIM_T)
C
C
C--------------------------------------------
C      Build the block-tridiagonal matrix T.
C--------------------------------------------
C
      DO i = 1, LLANCHAIN
         T_MATR(i,i) = E_DIAG(i)      
         T_MATR(LLANCHAIN+i,LLANCHAIN+i) = -E_DIAG(i)      
         T_MATR(i,LLANCHAIN+i) = D_DIAG(i)      
         T_MATR(LLANCHAIN+i,i) = -D_DIAG(i)      

         IF (I .LT. LLANCHAIN) THEN
             T_MATR(i,i+1) = A_OFFDIAG(i)
             T_MATR(i+1,i) = A_OFFDIAG(i)
             T_MATR(LLANCHAIN+i,LLANCHAIN+1+i) = -A_OFFDIAG(i)
             T_MATR(LLANCHAIN+1+i,LLANCHAIN+i) = -A_OFFDIAG(i)

             T_MATR(i,LLANCHAIN+1+i) = B_OFFDIAG(i)
             T_MATR(i+1,LLANCHAIN+i) = B_OFFDIAG(i)
             T_MATR(LLANCHAIN+i,i+1) = -B_OFFDIAG(i)
             T_MATR(LLANCHAIN+1+i,i) = -B_OFFDIAG(i)
         END IF
      END DO
C
C--------------------------------------------
C     Print the block-tridiagonal matrix T.
C--------------------------------------------
C
      IF ( IPRSOP .GE. 5 ) THEN
         CALL AROUND('--- The reordered Lanczos matrix --- ')
         CALL OUTPUT(T_MATR,1,DIM_EIGVAL,1,DIM_EIGVAL,
     &                     DIM_EIGVAL,DIM_EIGVAL,1,LUPRI)
      END IF
C      
C-------------------------------------------------
C     Diagonalize the block-tridiagonal matrix T.
C     Check for imaginary part in eigenvalues.
C-------------------------------------------------
C 
      IERR = 0
      CALL dgeev('V','V',DIM_EIGVAL,T_MATR,
     &           DIM_EIGVAL,EIGVAL,EIGVALC,
     &           EIGVECL,DIM_EIGVAL,
     &           EIGVECR,DIM_EIGVAL,
     &           WORK, LWORK, IERR)
C
      ICOMPLX = 0
      DO I = 1,DIM_EIGVAL
         IF (EIGVALC(I) .NE. ZERO) THEN
            ICOMPLX = ICOMPLX + 1
            WRITE(LUPRI,'(I10,1P,2D15.8,A/)') I,EIGVAL(I),
     &            EIGVALC(I),
     &            ' *** LANCZOS WARNING **** COMPLEX VALUE.'
         END IF
      END DO
C
C-------------------------------------------------------------------------
C     Reorder eigvalues and (left and right) eigvectors (ascending order).
C-------------------------------------------------------------------------
C
      CALL RGORD2(DIM_EIGVAL,DIM_EIGVAL,EIGVAL,EIGVALC,EIGVECR,EIGVECL,
     &            .FALSE.)
C
C--------------------------------------------------------------
C     Write eigvalues and (left and right) eigvectors to output.
C---------------------------------------------------------------
C 
      IF ( IPRSOP .GE. 2 ) THEN
        WRITE (LUPRI,*) ' '
        WRITE (LUPRI,'(A,I8,A)') ' After Lanczos iteration #',LLANCHAIN,
     &      ':  Eigenvalues in Lanczos basis '
        WRITE (LUPRI,*) '---------------------------------------------',
     &       '-----------------------------'
C       
        IOFFSET = 0
        DO I= 1,DIM_EIGVAL
           WRITE (LUPRI,'(A,I4,A,F12.6,A,F12.6)') 
     &   ' - EIGENVALUE NO.',I,
     &   ' - EIGENVALUE (au) ',EIGVAL(I),
     &   ' - Im part ', EIGVALC(I)
C           IF ( IPRSOP .GE. 5 ) THEN          
C               WRITE (LUPRI,*) '      Right eigenvector    Left eigenvector'
C               WRITE (LUPRI,'(I8,1X,F14.8,1X,F14.8)')
C     &         (J,EIGVECR(IOFFSET+J),EIGVECL(IOFFSET+J)
C     &         ,J=1,DIM_EIGVAL)
C           END IF
           IOFFSET = IOFFSET + DIM_EIGVAL
        END DO
        WRITE (LUPRI,*) ' '
      END IF
C
C----------------------------------------------------------------------
C    Bi-orthogonalize the left and right eigenvectors.
C----------------------------------------------------------------------
C     
      CALL RP_LANCZOS_BIORTH_EIGV(EIGVECR,EIGVECL,LLANCHAIN,
     &                                       EIGVAL,EIGVALC)
C     

C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_EIGV')
C
      RETURN
C
C
      END
