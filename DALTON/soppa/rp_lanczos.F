C
C  /* Deck rp_lanczos */
      SUBROUTINE RP_LANCZOS(MODEL,LABEL,ISYMTR,
     &                     DENSIJ,LDENSIJ,
     &                     DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,LT2MP,
     &                     FOCKD,LFOCKD,
     &                     LCONV,S_0,L_0,I_0,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Feb 2020
C
C     PURPOSE: Solve the RPA eigenvalue problem and compute mean
C     exitation energy, I(0), in Lanczos basis. Lanczos iterations are
C     run here for one component (x, y or z).
C
C     MODEL         Model, in this case RPA only 
C     LABEL         Label of the property to be computed
C     ISYMTR        Symmetry
C     DENSIJ        Occ.-occ. part of MP2 density
C     LDENSIJ       Length of the occ.-occ. part of MP2 density
C     DENSAB        Virt.-virt. part of MP2 density          
C     LDENSAB       Lenght of virt.-virt. part of MP2 density                     
C     DENSAI        Virt.-occ. part of MP2 density
C     LDENSAI       Lenght of virt.-occ. part of MP2 density 
C     T2MP          Array of the T2 amplitudes
C     LT2MP         Length of the T2 amplitudes array           
C     FOCKD         Fock diagonal 
C     LFOCKD        Length of Fock diagonal  
C     LCONV         # of diagonalizations performed
C     S_0           S(0) value array (length=1 if no intermediate
C     L_0           L(0) value array (length=1 if no intermediate
C     I_0           I(0) value array (length=1 if no intermediate
C     (arrays length is 1 if no intermediate diagonalization)          
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
C
#ifdef VAR_MPI
      use so_parutils, only: parsoppa_do_eres, my_mpi_integer,
     &                       soppa_nint
#endif
      use so_info, only: sop_linres, so_full_name, so_has_doubles,
     &                   so_needs_densai, sop_dp, sop_stat_trh,
     &                   so_model_number,
     &                   sop_mp2ai_done, sop_lanczos,
     &                   sop_lanc_chain_len, sop_lanc_conv_check,
     &                   sop_lanc_conv_nr 
C
      IMPLICIT NONE
C#include "implicit.h"
#ifdef VAR_MPI
#include "mpif.h"
C  IRAT in order to assign space for load-balancing
#include "iratdef.h"
#endif
#include "priunit.h"
C
#include "soppinf.h"
#include "codata.h"
#include "ccsdsym.h"
C
C-----------------
C     Parameters
C-----------------
C
      REAL(sop_dp), PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0
C
C----------------------
C     Formal arguments
C----------------------
C
      CHARACTER(LEN=5), INTENT(IN) :: MODEL 
      CHARACTER(LEN=8), INTENT(IN) :: LABEL
C
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LDENSIJ) :: DENSIJ
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LDENSAB) :: DENSAB 
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LDENSAI) :: DENSAI 
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LT2MP) :: T2MP 
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LFOCKD) :: FOCKD 
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LWORK) :: WORK 
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LCONV) :: S_0, L_0, I_0
C
      INTEGER, INTENT(IN) :: ISYMTR, LDENSIJ, LDENSAB, LDENSAI
      INTEGER, INTENT(IN) :: LT2MP, LFOCKD, LWORK 
      INTEGER, INTENT(IN) :: LCONV 
C
C      
#ifdef VAR_MPI
      INTEGER ::  CP_ISYMTR
      INTEGER :: MAXNUMJOBS
C     This array is only there to ensure that the four above variables
C     are allocated consecutively, so that it can be send together. Only
C     use it for this purpose.
C     The definition must match that in soppa_nodedriver
      INTEGER ::  INFO_ARRAY(6)
      EQUIVALENCE (info_array(1), cp_isymtr), (info_array(2),nit),
     &            (info_array(3), nnewtr),    (info_array(4),noldtr),
     &            (info_array(5), idtype), (info_array(6),imodel)
      INTEGER(MPI_INTEGER_KIND) :: IERR, numprocs
      INTEGER :: lAssignedIndices, kAssignedIndices
#endif
C
C---------------------
C     Local variables
C---------------------
C
      LOGICAL :: IMAGPROP
C
      LOGICAL :: DOUBLES, BREAKDOWN
      CHARACTER(LEN=11) :: FULL_NAME
      CHARACTER(LEN=1)  :: CHAR_LABEL
      INTEGER :: LANC_CHAIN_MAX
      INTEGER :: LLAN_CHAIN
      INTEGER :: NEXCI, NIT, NOLDTR, NNEWTR, IDTYPE, IMODEL
      INTEGER :: KEDIAG, KDDIAG, KAOFFDIAG, KBOFFDIAG, KEND1, LWORK1
      INTEGER :: KBTRID, KWI1, KBTRIDEIVAL, KBTRIDEIVECR, KBTRIDEIVECL,
     &           KEND2, LWORK2, KTRSSTR, KOSCSTR, KEND3, LWORK3
      INTEGER :: IERR1, ICOMPLX, IOFFSET, ICONV, ITURN
      REAL(SOP_DP) :: DTIME, DUMMY, WTIMES, WTIMEE, SECOND
      REAL(SOP_DP) :: GPNORM
C
      NNEWTR = 1
      NEXCI = 1
C      
C-------------------------------------------------------------------
C     Add to trace. Set model info.
C     Determine if Lanczos chain length is smaller than dim of the
C     irrep, adjust accordingly.
C-------------------------------------------------------------------
C
      CALL QENTER('RP_LANCZOS')
C
C     Basic info on what we do here
      DOUBLES = so_has_doubles(model)
      FULL_NAME = SO_FULL_NAME(MODEL)
      IMODEL = SO_MODEL_NUMBER(MODEL)
C
C      Left here for possible future use      
CC     Set up to update the GP vector if needed
C      UPDATE_GP = (.NOT. sop_mp2ai_done) .AND.
C     &            SO_NEEDS_DENSAI(MODEL)
C
      LANC_CHAIN_MAX = MIN(NT1AM(ISYMTR), sop_lanc_chain_len)
C      
C----------------------------------------------------------------------------
C     Allocate memory for the elements (e, d, a, b) of the lanczos T matrix.
C----------------------------------------------------------------------------
C
      KEDIAG    = 1
      KDDIAG    = KEDIAG + LANC_CHAIN_MAX
      KAOFFDIAG = KDDIAG + LANC_CHAIN_MAX
      KBOFFDIAG = KAOFFDIAG + LANC_CHAIN_MAX-1
      KEND1     = KBOFFDIAG + LANC_CHAIN_MAX-1
C
#ifdef VAR_MPI
C--------------------------------------------------------------------
C     For MPI, we need some space in which to store the indices each
C     process is to work with in SO_ERES.
C--------------------------------------------------------------------
C
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      maxnumjobs = soppa_nint - min(soppa_nint, numprocs) + 1
      if ( numprocs .eq. 1 ) then
C Not a real parallel job, don't bother
         lAssignedIndices = 1
         kAssignedIndices = 0
      else
         lAssignedIndices = (maxnumjobs + 1) /IRAT
         kAssignedIndices = KEND1
         KEND1 = kAssignedIndices + lAssignedIndices
      endif
#endif
      LWORK1 = LWORK - KEND1
C
      CALL SO_MEMMAX('RP_LANCZOS.1',LWORK1)
      IF(LWORK1.LT.0) CALL STOPIT('RP_LANCZOS.1','  ',
     &                             KEND1,LWORK)
C---------------------------------------
C     Write zeros to e, d, a, b arrays.
C---------------------------------------
C
      CALL DZERO(WORK(KEDIAG),LANC_CHAIN_MAX)
      CALL DZERO(WORK(KDDIAG),LANC_CHAIN_MAX)
      CALL DZERO(WORK(KAOFFDIAG),LANC_CHAIN_MAX-1)
      CALL DZERO(WORK(KBOFFDIAG),LANC_CHAIN_MAX-1)
C
C---------------------
C     Write a banner.
C---------------------
C
      WRITE(LUPRI,'(A)') ''
      WRITE(LUPRI,9001)
      WRITE(LUPRI,'(7X,3A)') ADJUSTR('Lanczos'),
     &        ' iterations, Property ',LABEL
      WRITE(LUPRI,9001)
      WRITE(LUPRI,'(A)') ''
C
C----------------------------------------
C     Create an initital Lanczos vector. 
C----------------------------------------
C      
      DTIME      = SECOND()
      IF ((MODEL .EQ. 'AORPA') .AND. SOP_LANCZOS ) THEN
         CALL  RP_LANCZOS_TRIALVEC(GPNORM,
     &                    LABEL,ISYMTR,IMAGPROP,MODEL,
     &                    T2MP,LT2MP,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    DENSAI,LDENSAI,WORK(KEND1),LWORK1)
      END IF
C
      DTIME      = SECOND()   - DTIME
      SOTIME(32) = SOTIME(32) + DTIME
C
C--------------------------------------------------------------------
C     Initialize iteration counter, number of old Lanczos vectors,
C     idtype, logical for a serious break-down and how many Lanczos
C     to run (k_max).
C     When creating a new Lanczos vectors with a help of ERES: the
C     de-excitation part needs to be done explicitly, ergo idtype = 0
C--------------------------------------------------------------------
C
      NIT    = 0
C
      NOLDTR = 0
C
      BREAKDOWN = .FALSE.
C
      IDTYPE = 0
C
      LLAN_CHAIN = LANC_CHAIN_MAX
C
C----------------------------------
C     Loop for Lanczos iterations.
C----------------------------------
C
      DO NIT = 1,LANC_CHAIN_MAX
C
C--------------------------------------------------------------
C        Count number of iterations and write header to output.
C--------------------------------------------------------------
C
         IF ( IPRSOP .GE. 5 ) THEN
C
            WRITE(LUPRI,'(/,2X,A)') '================================'//
     &                              '=================================='
            WRITE(LUPRI,'(11X,I3,3A,I1,A)') NIT,'. ',
     &              TRIM(FULL_NAME),
     &              '  iteration, Symmetry ', ISYMTR
            WRITE(LUPRI,'(2X,A,/)') '================================'//
     &                              '=================================='
C
         END IF
C
C-----------------------------------------------------------------
C        Make E linear transformation of a Lanczos vector giving
C        a resultvector.
C----------------------------------------------------------------
#ifdef VAR_MPI
C In parallel, send slaves to so_eres
C
         call mpi_bcast( parsoppa_do_eres, 1, my_mpi_integer, 0,
     &                   mpi_comm_world, ierr )
C ISYMTR is a non-local parameter, we need to copy it to the info-array
         CP_ISYMTR = ISYMTR
         CALL MPI_BCAST( INFO_ARRAY, 6, MY_MPI_INTEGER, 0,
     &                   MPI_COMM_WORLD, IERR)
#endif
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_ERES(MODEL,NOLDTR,NNEWTR,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                T2MP,LT2MP,FOCKD,LFOCKD,DENSAI,LDENSAI,NIT,ISYMTR,
     &                IDTYPE,
#ifdef VAR_MPI
     &                WORK(kAssignedIndices),maxnumjobs,
#endif
     &                WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(1)  = SOWTIM(1)  + WTIMEE - WTIMES
C
C-----------------------------------
C        Perform Lanczos recursion.
C-----------------------------------
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
C
         CALL RP_LANCZOS_ITER(NOLDTR,ISYMTR,
     &                         WORK(KEDIAG),WORK(KDDIAG),
     &                         WORK(KAOFFDIAG),WORK(KBOFFDIAG),
     &                         LANC_CHAIN_MAX,BREAKDOWN,
     &                         WORK(KEND1),LWORK1)
C
         DTIME      = SECOND()   - DTIME
         SOTIME(33) = SOTIME(33) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(3)  = SOWTIM(3)  + WTIMEE - WTIMES
C
C--------------------------------------------------------------------------------
C        Increase the number of created Lanczos vectors and check if a
C        serious break-down was encountered.
C        Diagonalize either intermediate or final Lanczos T matrix (the
C        left and right eigenvectors need to be bi-orthonormalized),
C        compute transition strengths, then oscillator strengths and
C        calculate the mean excitation energy for the current component (x,y or
C        z). Exit the loop if a serious break-dow was encountered. 
C--------------------------------------------------------------------------------
C         
         NOLDTR = NOLDTR + 1
C         
         IF (BREAKDOWN) THEN
             ! since the noldtr already incremented
             LLAN_CHAIN = NOLDTR  
         END IF 
C         
         IF (sop_lanc_conv_check) THEN
             ICONV = MOD(NOLDTR, sop_lanc_conv_nr)
             ITURN = INT(NOLDTR/ sop_lanc_conv_nr)
C         
             IF ((NOLDTR .EQ. LLAN_CHAIN) .AND. (ICONV .NE. 0) ) THEN
                 ITURN = ITURN + 1
             END IF 
C         
         ELSE
             ICONV = 1
             ITURN = LCONV
         END IF
C         
         IF  ( (ICONV .EQ. 0) .OR. (NOLDTR .EQ. LLAN_CHAIN) ) THEN
            KBTRID = KEND1
            KWI1 = KBTRID + 4*NOLDTR*NOLDTR
            KBTRIDEIVAL = KWI1 + 2*NOLDTR
            KBTRIDEIVECR = KBTRIDEIVAL+ 2*NOLDTR
            KBTRIDEIVECL = KBTRIDEIVECR+ 
     &               4*NOLDTR*NOLDTR
C
            KEND2 = KBTRIDEIVECL+ 4*NOLDTR*NOLDTR 
C     
            LWORK2 = LWORK - KEND2
C
            CALL SO_MEMMAX('RP_LANCZOS.1',LWORK2)
            IF(LWORK2.LT.0) CALL STOPIT('RP_LANCZOS.1','  ',
     &                             KEND2,LWORK)
C            
            CALL RP_LANCZOS_EIGV(WORK(KEDIAG),WORK(KDDIAG),
     &                           WORK(KAOFFDIAG),WORK(KBOFFDIAG),NOLDTR,
     &                                    WORK(KBTRIDEIVAL),WORK(KWI1),
     &                          WORK(KBTRIDEIVECR),WORK(KBTRIDEIVECL),
     &                                  WORK(KBTRID),WORK(KEND2),LWORK2)
C            
            CALL RP_LANCZOS_BIORTH_EIGV(WORK(KBTRIDEIVECR),
     &                                  WORK(KBTRIDEIVECL),
     &                     NOLDTR,WORK(KBTRIDEIVAL),WORK(KWI1))
C            
            KTRSSTR = KEND2
            KOSCSTR = KTRSSTR + NOLDTR
            KEND3   = KOSCSTR + NOLDTR
            LWORK3  = LWORK - KEND3
C
            CALL SO_MEMMAX('RP_LANCZOS.1',LWORK3)
            IF(LWORK3.LT.0) CALL STOPIT('RP_LANCZOS.1','  ',
     &                             KEND3,LWORK)
C
            CALL RP_LANCZOS_TRS_STR(WORK(KBTRIDEIVECR),
     &                   WORK(KBTRIDEIVECL),NOLDTR,WORK(KTRSSTR),GPNORM)
C
            CALL RP_LANCZOS_OSC_STR(WORK(KBTRIDEIVAL),WORK(KWI1),
     &                        NOLDTR,WORK(KTRSSTR),
     &                        WORK(KOSCSTR))
C
C
            CALL RP_LANCZOS_MEAN_EXC(WORK(KBTRIDEIVAL),WORK(KWI1),
     &         NOLDTR,S_0(ITURN),L_0(ITURN),I_0(ITURN),WORK(KOSCSTR))
C            
            IF (BREAKDOWN) THEN
                WRITE(LUPRI,*)'EXITING after a serious break-down'
                EXIT
            END IF
C            
         END IF         
C            
C
C     Lanczos loop         
      END DO
C   
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_LANCZOS')
C
      RETURN
C
 9001 FORMAT(1X,'=============================================',
     &       '================')
C
      END
