C
C  /* Deck rp_lanczos_drv */
      SUBROUTINE RP_LANCZOS_DRV(WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Luna Zamokaite, Sep 2020
C
C     PURPOSE: Main driver routine for solving the RPA eigenvalue
C              problem and computing the mean excitation energy in
C              Lanczos basis.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      use so_info, only: sop_num_models,
     &                   sop_models, sop_mod_fullname,
     &                   so_num_active_models,
     &                   so_get_active_models,
     &                   so_has_doubles,
     &                   so_singles_second,
     &                   sop_dens_label, sop_model_rpad,
     &                   sop_model_rpa,
     &                   sop_conv_thresh,
     &                   so_needs_densai,
     &                   sop_mp2ai_done,
     &                   sop_lanczos,sop_lanc_chain_len,
     &                   sop_lanc_conv_check,
     &                   sop_lanc_conv_nr, sop_dp 
C
#ifdef VAR_MPI
      use so_parutils, only: soppa_initialize_slaves,
     &                       soppa_release_slaves
#endif
C
      IMPLICIT NONE
C      
#include "priunit.h"
C
C XTEV for converting a.u. to eV
#include "codata.h"
C
#include "soppinf.h"
#include "cbilnr.h"
#include "cbiexc.h"
#include "absorp.h"
#include "ccsdsym.h"
#include "inforb.h"
C for irat
#include "iratdef.h"
#include "maxaqn.h"
C for mxshel, mxprim
#include "maxorb.h"
C for mxcont
#include "aovec.h"
C NEED DOCCSD parameter from gnring.h
#include "gnrinf.h"
C
C----------------
C     Parameters
C----------------
C      
      REAL(SOP_DP), PARAMETER :: THREE = 3.0D+00
C
C---------------------------------
C     Dimensions of the arguments
C---------------------------------
C
      INTEGER, INTENT(IN) :: LWORK
      REAL(sop_dp), INTENT(INOUT), DIMENSION(LWORK) :: WORK
C
C---------------------
C     Local variables
C---------------------
C
      CHARACTER(LEN=5) :: MODEL
      CHARACTER(LEN=8) :: LABEL1
      CHARACTER(LEN=11) :: FULLNAME
C
      INTEGER :: LFOCKD, LT2AM, LDENSIJ, LDENSAB, LDENSAI, LWORK1
      INTEGER :: KFOCKD, KT2AM, KDENSIJ, KDENSAB, KDENSAI, KEND1
      INTEGER :: NVARPT
C      
      INTEGER :: LCONV1(3),KSSUM(3),KLSUM(3),KISUM(3)
      INTEGER :: ISYM, ISYMTR, IOPER, IDIP, IMODEL
C
      REAL(SOP_DP) :: STOTAL, LTOTAL, ITOTAL
C
      REAL(SOP_DP) :: DTIME, SECOND, DUMMY, WTIME, TIMWTO, TIMEIN
      REAL(SOP_DP) :: TIMOUT, TIMTOT
C      
      INTEGER :: LPARRA, KPARRA, KEND2, LWORK2
C      
      integer :: num_active
      logical :: active_models(sop_num_models)
      logical :: doubles, need_t2, rpa_only
      character(len=4) :: old_denslab, denslab
      logical :: new_densai, imagprop
#ifdef VAR_MPI
!
! This variable ensures that common blocks are only sent to the slaves
! once.
      LOGICAL update_common_blocks

      update_common_blocks = .true.
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LANCZOS_DRV')
C
C--------------------------------------------------------------------
C     Get the active models and their numbers, set no RPA(D). Quit if 
C     more than RPA requested.
C     Set some logicals (relevant if the routine will be extended to
C     SOPPA model in the future).      
C--------------------------------------------------------------------
C
      CALL SO_GET_ACTIVE_MODELS (ACTIVE_MODELS)
      NUM_ACTIVE = COUNT ( ACTIVE_MODELS )
C
      ACTIVE_MODELS(SOP_MODEL_RPAD) = .FALSE.
C
      RPA_ONLY = ACTIVE_MODELS(SOP_MODEL_RPA) .AND. 
     &           (NUM_ACTIVE.EQ.1) 
C
      IF (.NOT. RPA_ONLY) 
     &   CALL QUIT('LANCZOS solver only implemented for RPA')
C
      TRIPLET = .FALSE.
C
      OLD_DENSLAB = 'NONE'
C      
C----------------------------------
C     Start timing the calculation.
C----------------------------------
C
      CALL TIMER('START ',TIMEIN,TIMOUT)
C
      DTIME   = SECOND()
      TIMTOT  = DTIME
      CALL GETTIM (DUMMY,WTIME)
      TIMWTO  = WTIME
C
C----------------------------------------------------
C     Initialize a few pointer arrays and dimensions.
C----------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_INIT(WORK,LWORK)
      DTIME     = SECOND()  - DTIME
      SOTIME(2) = SOTIME(2) + DTIME
C
C---------------------
C     Set print level.
C---------------------
C
      IPRSOP = IPRLNR
C
      IF (IPRLNR .GE. 0)
     &     CALL TITLER('BLOCK LANCZOS RPA SOLVER','#',118)
C
C-------------------------------------------------------------
C     Allocate work space for Fock diagonal, MP2 densities, T2
C     amplitudes.
C-------------------------------------------------------------
C
      LFOCKD  = NORBT
C
      KFOCKD  = 1
      KEND1   = KFOCKD  + LFOCKD
      LWORK1  = LWORK   - KEND1
C
C
      IF (.NOT.RPA_ONLY) THEN
C
         LT2AM   = NT2AMX
         LDENSIJ = NIJDEN(1)
         LDENSAB = NABDEN(1)
         LDENSAI = NAIDEN(1)
C
         KT2AM   = KFOCKD  + LFOCKD
         KDENSIJ = KT2AM   + LT2AM
         KDENSAB = KDENSIJ + LDENSIJ
         KDENSAI = KDENSAB + LDENSAB
         KEND1   = KDENSAI + LDENSAI
         LWORK1  = LWORK   - KEND1
      ELSE
         LT2AM   = 0
         LDENSIJ = 0
         LDENSAB = 0
         LDENSAI = 0
C        Force a crash if this is used
         KT2AM   = HUGE(LWORK)
         KDENSIJ = HUGE(LWORK)
         KDENSAB = HUGE(LWORK)
         KDENSAI = HUGE(LWORK)
      END IF
C
C--------------------------------------------------------------------
C     If intermediate diagonalization is requested: calculate the
C     length of S(0), L(0) and I(0) arrays for x, y and z components. 
C     Otherwise these arrays have length 1. 
C     Allocate work space for these arrays accordingly and initiate
C     them to zeros.      
C--------------------------------------------------------------------
C
      IF (SOP_LANC_CONV_CHECK) THEN
          LCONV1 = INT(SOP_LANC_CHAIN_LEN / sop_lanc_conv_nr) 
          IF ( MOD(SOP_LANC_CHAIN_LEN,sop_lanc_conv_nr).NE. 0 ) THEN
            LCONV1 = LCONV1 + 1
          END IF
          DO ISYM = 1, NSYM
            DO IOPER = 1, NOPER(ISYM)
                LABEL1 = LABOP(IOPER,ISYM)
                IF (LABEL1(2:7).NE.'DIPLEN') CYCLE
                IDIP = ICHAR(LABEL1(1:1)) - ICHAR('X') + 1
                ! check if irrep shorter than k_max
                IF (NT1AM(ISYM) .LT. sop_lanc_chain_len) THEN
                    LCONV1(IDIP) = INT(NT1AM(ISYM) / sop_lanc_conv_nr)
                    IF ( MOD(NT1AM(ISYM),sop_lanc_conv_nr) .NE. 0 ) THEN
                        LCONV1(IDIP) = LCONV1(IDIP) + 1
                    END IF
                END IF
            END DO
          END DO
      ELSE
          LCONV1 = 1
C          
      END IF 
C      
      DO IDIP = 1,3
        KSSUM(IDIP) = KEND1
        KEND1 = KEND1 + LCONV1(IDIP)
        KLSUM(IDIP) = KEND1
        KEND1 = KEND1 + LCONV1(IDIP)
        KISUM(IDIP) = KEND1
        KEND1 = KEND1 + LCONV1(IDIP)
      END DO
C
      LWORK1  = LWORK   - KEND1
C

      CALL SO_MEMMAX ('RP_LANCZOS_DRV.1',LWORK1)
      IF (LWORK1 .LT .0) CALL STOPIT('RP_LANCZOS_DRV.1',' ',KEND1,LWORK)
C
      DO IDIP = 1,3
          CALL DZERO(WORK(KSSUM(IDIP)),LCONV1(IDIP))
          CALL DZERO(WORK(KLSUM(IDIP)),LCONV1(IDIP))
          CALL DZERO(WORK(KISUM(IDIP)),LCONV1(IDIP))
      END DO
C
C------------------------------------------------
C     Get MO-energies (the Fock matrix diagonal).
C------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_MOENERGY(WORK(KFOCKD),WORK(KEND1),LWORK1)
      DTIME     = SECOND()  - DTIME
      SOTIME(3) = SOTIME(3) + DTIME
C
C      
#ifdef VAR_MPI
C
C---------------------------------------------------------
C              Ready the slaves for parallel calculations.
C---------------------------------------------------------
C
      call soppa_initialize_slaves ( update_common_blocks, rpa_only)
      update_common_blocks = .false.
#endif
C
C======================================
C     Initiate the Lanczos calculation.
C======================================
C
C--------------------------------------------------
C     Loop over symmetry of the property operators.
C--------------------------------------------------
C
      DO ISYM = 1, NSYM
C
         DO IOPER = 1, NOPER(ISYM)
C
            LABEL1 = LABOP(IOPER,ISYM)
            IF (LABEL1(2:7).NE.'DIPLEN') CYCLE

            IDIP = ICHAR(LABEL1(1:1)) - ICHAR('X') + 1

            ISYMTR = ISYM
C            
C            
C--------------------------------------------
C           Loop over the posible methods.
C           (left for future implementations)
C--------------------------------------------
C
            DO IMODEL = 1, SOP_NUM_MODELS
C
C              Skip methods not treated in this calcultion
               IF (.NOT.ACTIVE_MODELS(IMODEL) ) CYCLE
C
C              Look up info on this model
C              (the model number -> model conversion is defined in so_info.F90)
C
               MODEL = SOP_MODELS(IMODEL)
               FULLNAME = SOP_MOD_FULLNAME(IMODEL)
               DOUBLES = SO_HAS_DOUBLES(IMODEL)
               NEED_T2 = DOUBLES.OR.SO_SINGLES_SECOND(MODEL)
               DENSLAB = SOP_DENS_LABEL(IMODEL)
C
               IF (IPRSOP .GE. 0) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A,2(/A,I8))') TRIM(FULLNAME)//':',
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h + h-p variables       :',2*NT1AM(ISYM)
                  IF (DOUBLES) THEN
                     NVARPT = 2*(NT1AM(ISYM)+NT2AM(ISYM))
                     WRITE(LUPRI,'(A,I8,/A,I8)')
     &               ' 2p-2h + 2h-2p variables   :',2*NT2AM(ISYM),
     &               ' Total number of variables :',NVARPT
                  END IF
                  WRITE(LUPRI,9001)
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C              (left for future implementations)
C-----------------------------------------------------
C
               IF ((NEED_T2).AND.
     &             (DENSLAB.NE.OLD_DENSLAB)) THEN
                  WRITE(LUPRI,*) 'WARNING: NEED_T2 should not be True!'
                  CALL SO_GETT2(DENSLAB,
     &                          WORK(KFOCKD),LFOCKD,WORK(KT2AM),LT2AM,
     &                          WORK(KDENSAI),
     &                          LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                          WORK(KDENSAB),LDENSAB,
     &                          WORK(KEND1),LWORK1)
                  OLD_DENSLAB = DENSLAB
               ENDIF
C
C-----------------------------------------------
C              Initialize DENSAI if needed.
C              (left for future implementations)
C-----------------------------------------------
C
               IF ( SO_NEEDS_DENSAI(MODEL) .AND.
     &              (.NOT.SOP_MP2AI_DONE) ) THEN
                  CALL DZERO(WORK(KDENSAI),LDENSAI)
                  NEW_DENSAI = .TRUE.
                  WRITE(LUPRI,*) 'WARNING: NEW_DENSAI true!'
               ELSE
                  NEW_DENSAI = .FALSE.
               ENDIF
C
C-----------------------------------------------------------------------
C              Run the requested # of Lanczos iterations for the current
C              component/symmetry. Compute S(0), L(0), I(0).               
C-----------------------------------------------------------------------
C
               CALL RP_LANCZOS(MODEL,LABEL1,ISYMTR,
     &             WORK(KDENSIJ), LDENSIJ, WORK(KDENSAB), LDENSAB,
     &             WORK(KDENSAI), LDENSAI, WORK(KT2AM), LT2AM,
     &             WORK(KFOCKD), LFOCKD, LCONV1(IDIP),
     &             WORK(KSSUM(IDIP)),          
     &             WORK(KLSUM(IDIP)),          
     &             WORK(KISUM(IDIP)),          
     &             WORK(KEND1),LWORK1)
C
C------------------------------------------------------------------
C              Save AI density matrix, if it has been recalculated.
C              (left for future implementations)
C------------------------------------------------------------------
C
               IF (NEW_DENSAI) THEN
                  CALL SO_DUMP_DENSAI(DENSLAB,WORK(KDENSAI),LDENSAI)
               END IF
C
C           Loop over models    
            END DO
C
C        Loop over operators
         END DO
C
C     Loop over symmetry groups
      END DO
C 
C
C--------------------------------------------------
C     Print a table of S(0), L(0) and I(0)i values.
C--------------------------------------------------
C
      CALL AROUND('Oscillator Strength Sum Rules in Lanczos basis')
C
      WRITE (LUPRI,'(//,18X,A,/,5X,A,10X,A,10X,A,10X,A,/)')
     &   'S(0) Sum Rule : Dipole Length Approximation in a.u.',
     &   'xx - component','yy - component','zz - component','total'
      DO J = 1, MAXVAL(LCONV1)
        STOTAL = WORK(KSSUM(1)+MIN(J,LCONV1(1))-1) + 
     &           WORK(KSSUM(2)+MIN(J,LCONV1(2))-1) +
     &           WORK(KSSUM(3)+MIN(J,LCONV1(3))-1) 
        STOTAL = STOTAL/THREE
        WRITE (LUPRI,'(4(F13.6,11X))') 
     $        (WORK(KSSUM(K)+MIN(J,LCONV1(K))-1),K=1,3), STOTAL 
      END DO
C
      WRITE (LUPRI,'(//,18X,A,/,5X,A,10X,A,10X,A,10X,A,/)')
     &   'L(0) Sum Rule : Dipole Length Approximation in a.u.',
     &   'xx - component','yy - component','zz - component','total'
      DO J = 1, MAXVAL(LCONV1)
        LTOTAL = WORK(KLSUM(1)+MIN(J,LCONV1(1))-1) + 
     &           WORK(KLSUM(2)+MIN(J,LCONV1(2))-1) +
     &           WORK(KLSUM(3)+MIN(J,LCONV1(3))-1) 
        LTOTAL = LTOTAL/THREE
        WRITE (LUPRI,'(4(F13.6,11X))') 
     $        (WORK(KLSUM(K)+MIN(J,LCONV1(K))-1),K=1,3), LTOTAL 
      END DO
C
      WRITE (LUPRI,'(//,18X,A,/,4X,A,10X,A,10X,A,10X,A,/)')
     &   'I(0) Sum Rule : Dipole Length Approximation in eV',
     &   'xx - component','yy - component','zz - component','total'
      DO J = 1, MAXVAL(LCONV1)
        LTOTAL = WORK(KLSUM(1)+MIN(J,LCONV1(1))-1) + 
     &           WORK(KLSUM(2)+MIN(J,LCONV1(2))-1) +
     &           WORK(KLSUM(3)+MIN(J,LCONV1(3))-1) 
        STOTAL = WORK(KSSUM(1)+MIN(J,LCONV1(1))-1) + 
     &           WORK(KSSUM(2)+MIN(J,LCONV1(2))-1) +
     &           WORK(KSSUM(3)+MIN(J,LCONV1(3))-1) 
        ITOTAL = DEXP(LTOTAL/STOTAL) * XTEV
        WRITE (LUPRI,'(4(F13.6,11X))') 
     $        (WORK(KISUM(K)+MIN(J,LCONV1(K))-1),K=1,3), ITOTAL 
      END DO
C
C
#ifdef VAR_MPI
C-------------------------------------------------------
C              Release slaves to the global node-driver.
C-------------------------------------------------------
      call soppa_release_slaves()
#endif
C
C---------------------------------------------------
C     Allocate work space for printing timing stats?
C---------------------------------------------------
C
      LPARRA = LSOTIM
C
      KPARRA = KEND1
      KEND2  = KPARRA + LPARRA
      LWORK2 = LWORK  - KEND2
C
      CALL SO_MEMMAX ('RP_LANCZOS_DRV.2     ',LWORK2)
      IF (LWORK2 .LT.0) CALL STOPIT('RP_LANCZOS_DRV.2',' ',KEND2,LWORK)
C
C---------------------------------------------------
C     Print memory statistics for SOPPA subroutines.
C---------------------------------------------------
C
      CALL SO_MEMMAX('STATISTICS      ',0)
C
C-----------------------------------------
C     Print timings for SOPPA subroutines.
C-----------------------------------------
C
      CALL SO_TIME(TIMTOT,TIMWTO,WORK(KPARRA),LPARRA)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LANCZOS_DRV')
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END
