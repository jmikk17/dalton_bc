C
C  /* Deck rp_lanczos_eigv_block.f */
      SUBROUTINE RP_LANCZOS_EIGV_BLOCK(E_DIAG,D_DIAG,A_OFFDIAG,B_OFFDIAG
     &                       ,LLANCHAIN,EIGVAL,EIGVALC,EIGVECR,EIGVECL,
     &                                              T_MATR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Dec 2019
C
C     PURPOSE: Build the block-tridiagonal matrix T from the coefficients 
C              computed in Lanczos iterations. Diagonalize it, reorder
C              the left and right eigenvectors and bi-orthonormalize them.
C              ATTENTION! Using this routine instead of the
C              rp_lanczos_eigv.f requires adjusting the
C              rp_lanczos_trs_str routine accordingly - the elements
C              picked from the eigenvectors will have different indices.            
C              TODO: introduce handling of complex eigenvalues. In case complex
C              eigenvalues are found, DGEEV returns the complex
C              eigenvalues and eigenvectors as pairs (see LAPACK documentation)
C              and this routine will not handle them correctly. Also:
C              could benefit from a diagonalization routine that
C              exploits the sparsity (nonsymmetric banded matrix
C              diagonalization routine).
C          
C     E_DIAG       Diagonal elements in the diagonal blocks in T matrix. 
C     D_DIAG       Offdiagonal elements in the diagonal blocks in T matrix. 
C     A_OFFDIAG    Diagonal elements in the offdiagonal blocks in T matrix. 
C     B_OFFDIAG    Offdiagonal elements in the offdiagonal blocks in T matrix. 
C     LLANCHAIN    The length of Lanczos chain (equal to 1/2 dimesnion of T).         
C     EIGVAL       The real part of eigenvalues      
C     EIGVALC      The imaginary part of eigenvalues
C     EIGVECR      The right eigenvectors
C     EIGVECL      The left eigenvectors
C     T_MATR       The block-tridiagonal T matrix
C
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
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      INTEGER, INTENT(IN) :: LLANCHAIN,LWORK
      REAL(SOP_DP), INTENT(INOUT), DIMENSION(LWORK) :: WORK
      REAL(SOP_DP), INTENT(IN), DIMENSION(LLANCHAIN) :: E_DIAG,D_DIAG
      REAL(SOP_DP), INTENT(IN), DIMENSION(LLANCHAIN-1) :: 
     &                                              A_OFFDIAG,B_OFFDIAG 
      REAL(SOP_DP), INTENT(OUT), DIMENSION(2*LLANCHAIN) :: 
     &                                                    EIGVAL,EIGVALC
      REAL(SOP_DP), INTENT(OUT), DIMENSION(4*LLANCHAIN*LLANCHAIN) ::
     &                                            T_MATR,EIGVECR,EIGVECL
C
C---------------------
C     Local variables
C---------------------
C
      INTEGER :: DIM_T,DIM_EIGVAL
      INTEGER :: ICOMPLX,IOFFSET,IERR
C
C------------------------
C     BLAS functions used
C------------------------
C
C      
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LANCZOS_EIGV_BLOCK')
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
         T_MATR((2*i-1) + DIM_EIGVAL*(2*i-2)) = E_DIAG(I)      
         T_MATR((2*i) + DIM_EIGVAL*(2*i-1)) = -E_DIAG(I)      
         T_MATR((2*i-1) + DIM_EIGVAL*(2*i-1)) = D_DIAG(I)      
         T_MATR((2*i) + DIM_EIGVAL*(2*i-2)) = -D_DIAG(I)      

         IF (I .LT. LLANCHAIN) THEN
             T_MATR((2*i-1) + DIM_EIGVAL*(2*i)) = A_OFFDIAG(I)
             T_MATR((2*i) + DIM_EIGVAL*(2*i+1)) = -A_OFFDIAG(I)
             T_MATR((2*i-1) + DIM_EIGVAL*(2*i+1)) = B_OFFDIAG(I)
             T_MATR((2*i) + DIM_EIGVAL*(2*i)) = -B_OFFDIAG(I)

             T_MATR((2*i+1) + DIM_EIGVAL*(2*i-2)) = A_OFFDIAG(I)
             T_MATR((2*i+2) + DIM_EIGVAL*(2*i-1)) = -A_OFFDIAG(I)
             T_MATR((2*i+1) + DIM_EIGVAL*(2*i-1)) = B_OFFDIAG(I)
             T_MATR((2*i+2) + DIM_EIGVAL*(2*i-2)) = -B_OFFDIAG(I)
         END IF
      END DO
C
C--------------------------------------------
C     Print the block-tridiagonal matrix T.
C--------------------------------------------
C
      IF ( IPRSOP .GE. 2 ) THEN
         CALL AROUND('--- The block-tridiagonal Lanczos matrix --- ')
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
C------------------------------------------------------------------------
C     Reorder eigvalues and (left and right) eigvectors (ascending order).
C------------------------------------------------------------------------
C
      CALL RGORD2(DIM_EIGVAL,DIM_EIGVAL,EIGVAL,EIGVALC,EIGVECR,EIGVECL,
     &            .FALSE.)
C
C--------------------------------------------------------------
C     Write eigvalues and (left and right) eigvectors to output.
C---------------------------------------------------------------
C 
      IF ( IPRSOP .GE. 2 ) THEN
        WRITE (LUPRI,'(A,I6,A)') 'After Lanczos iteration #',LLANCHAIN,
     &      ':  Eigenvalues in Lanczos basis '
C       
        IOFFSET = 0
        DO I= 1,DIM_EIGVAL
           WRITE (LUPRI,'(A,I4,A,F12.6,A,F12.6)') !F13.7>22.16
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
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_EIGV_BLOCK')
C
      RETURN
C
C
      END
