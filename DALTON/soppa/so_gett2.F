C
C  /* Deck so_gett2 */
      SUBROUTINE SO_GETT2(KEY,FOCKD,LFOCKD,T2AM,LT2AM,T2SAM,LT2SAM,
     &                    DENSAI,LDENSAI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    DENS3IJ,DENS3AB,WORK,LWORK)
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Eiter read from disk or calculate the appropriate T2 and possibly T1 amplitudes
C              for RPA(D), SOPPA, SOPPA(CCSD) and TOPPA and calculate
C              the corresponding densities.
#ifdef VAR_MPI
      use so_parutils, only: soppa_update_amplitudes
#endif
C
      use so_info, only : sop_mp2ai_done, sop_dp
C#include "implicit.h"
      implicit none

#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      REAL(SOP_DP), PARAMETER :: HALF = 0.5D0
C
C     Arguments
      character(len=4), intent(in) :: KEY
      real(sop_dp), intent(in) :: FOCKD(LFOCKD)
      real(sop_dp), intent(out) :: T2AM(LT2AM), T2SAM(LT2SAM),
     &                             DENSAI(LDENSAI),
     &                             DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                             DENS3IJ(LDENSIJ), DENS3AB(LDENSAB)
      real(sop_dp), intent(inout) :: WORK(LWORK)
      integer, intent(in) :: LFOCKD, LT2AM, LT2SAM, LDENSAI, LDENSIJ, 
     &                       LDENSAB, LWORK

C     Local Variables
      LOGICAL DOCCSD, LEXIST
      CHARACTER(LEN=8) :: FNAME, FNAME2
      LOGICAL :: FOUND_DENS
      INTEGER :: IOSTATR, IDUMMY, LUTAM
      real(sop_dp) :: dtime
C     Time-function.
      real(sop_dp) :: second
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_GETT2')
C
      IF ( KEY(1:2) .EQ. 'MP' ) THEN
         FNAME = 'MP2__TAM'
         IF ( KEY .EQ. 'MP3 ' ) THEN
            FNAME2 = 'MP3__TAM'
         ENDIF
      ELSE IF ( KEY .EQ. 'CC2 ' ) THEN
         FNAME = 'CC2__TAM'
      ELSE IF ( KEY .EQ. 'CCSD' ) THEN
         FNAME = 'CCSD_TAM'
      ELSE
         CALL QUIT('Invalid KEY in SO_GETT2: '//KEY)
      END IF
C
C----------------------------------------------------------
C        Check if T2-amplitudes has been calculated in CC:
C----------------------------------------------------------
C
      INQUIRE (FILE=FNAME,EXIST=LEXIST)
      IF ((.NOT.LEXIST) .AND. (KEY .EQ. 'MP2 ')) THEN
         WRITE(LUPRI, '(2A)') 'The (ai|bj) integrals for the T2 ',
     &                        'amplitudes are calculated in SOPPA'
C
C----------------------------------------
C        Calculate (ia|jb) integrals .
C----------------------------------------
C
         DTIME     = SECOND()
         CALL SO_IAJB(T2AM,WORK,LWORK)
         DTIME     = SECOND()  - DTIME
         SOTIME(5) = SOTIME(5) + DTIME
C
C-----------------------------------------
C        Calculate MP2 T2 amplitudes
C-----------------------------------------
C
         DTIME     = SECOND()
         CALL SO_MP2T2(FOCKD,LFOCKD,T2AM,LT2AM,WORK,LWORK)
         DTIME     = SECOND()  - DTIME
         SOTIME(5) = SOTIME(5) + DTIME
         FOUND_DENS = .FALSE.
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,FNAME,'NEW',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUTAM
         WRITE (LUTAM) T2AM
C
      ELSE IF (LEXIST) THEN
         WRITE(LUPRI, '(3A)') 'The T2 ampliltudes are read ',
     &                        'from the file: ',FNAME
C-----------------------------------------------
C        Read T1 and T2-amplitudes from disk.
C-----------------------------------------------
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,FNAME,'OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUTAM
C
C        T1 amplitudes only if from CC
C
         IF ( KEY(1:2) .EQ. 'CC' ) THEN
            READ(LUTAM,END=901) (DENSAI(I), I = 1,LDENSAI)
            CALL DSCAL(LDENSAI,HALF,DENSAI,1)
         END IF
C
         READ(LUTAM,END=902) (T2AM(I), I = 1,LT2AM)
C
C-------------------------------------------
C        Write T2-amplitudes to output file.
C-------------------------------------------
C
         IF ( IPRSOP .GT. 10 ) THEN
            IF ( KEY(1:2) .EQ. 'MP' ) THEN
               CALL HEADER('MP2 Amplitudes',-1)
            ELSE IF ( KEY .EQ. 'CC2 ' ) THEN
               CALL HEADER('CC2 Amplitudes',-1)
            ELSE IF ( KEY .EQ. 'CCSD' ) THEN
               CALL HEADER('CCSD Amplitudes',-1)
            END IF
            WRITE(LUPRI, '(1X,D20.13)') (T2AM(I),I=1, LT2AM)
         END IF
C
C------------------------------------------------------------
C        Attempt to read the D(ij) and D(ab) density matrices
C------------------------------------------------------------
C
         READ(LUTAM,IOSTAT=IOSTATR) DENSIJ
         IF (IOSTATR.GT.0) THEN
            CALL QUIT('Error reading file '//FNAME)
         ELSE IF (IOSTATR.LT.0) THEN
            FOUND_DENS = .FALSE.
            BACKSPACE (LUTAM)
         ELSE
            READ(LUTAM) DENSAB
C
C        And also for MP2 get D(ai) block
            IF ( (KEY(1:2).NE.'CC').AND.sop_mp2ai_done) THEN
               READ(LUTAM,END=901) (DENSAI(I), I = 1,LDENSAI)
            END IF
            FOUND_DENS = .TRUE.
         END IF
C
      ELSE
         CALL QUIT('SOPPA(CC) requires running the CC module')
      END IF
C
C----------------------------------------------------
C     Calculate the density matrices D(ij) and D(ab).
C----------------------------------------------------
C
      IF (.NOT.FOUND_DENS) THEN
C
         DTIME     = SECOND()
         CALL SO_DENS(DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                T2AM,LT2AM,WORK,LWORK)
         DTIME     = SECOND()  - DTIME
         SOTIME(4) = SOTIME(4) + DTIME
C
         WRITE (LUTAM) DENSIJ
         WRITE (LUTAM) DENSAB
C
      END IF
C
      CALL GPCLOSE(LUTAM,'KEEP')
C
C Juliane 24/05/22
C----------------------------------------------------
C     The second order T2-amplitudes only if TOPPA
C----------------------------------------------------

      IF (KEY.EQ.'MP3 ') THEN
         INQUIRE (FILE=FNAME2,EXIST=LEXIST)
         IF (LEXIST) THEN
            WRITE(LUPRI, '(3A)') 'The T2 MP3 ampliltudes are read ',
     &                        'from the file: ',FNAME2
C
C-----------------------------------------------
C        Read second order T2-amplitudes from disk.
C-----------------------------------------------
C
            LUTAM = -1
            CALL GPOPEN(LUTAM,FNAME2,'OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
            REWIND LUTAM
C
            READ(LUTAM,END=902) (T2SAM(I), I = 1,LT2SAM)
C
C-------------------------------------------
C        Write second order T2-amplitudes to output file.
C-------------------------------------------
C
            IF ( IPRSOP .GT. 10 ) THEN
               CALL HEADER('MP3 Amplitudes',-1)
               WRITE(LUPRI, '(1X,D20.13)') (T2SAM(I),I=1, LT2SAM)
            END IF
C
C------------------------------------------------------------
C        Attempt to read the D3(ij) and D3(ab) density matrices
C------------------------------------------------------------
C
            READ(LUTAM,IOSTAT=IOSTATR) DENS3IJ
            IF (IOSTATR.GT.0) THEN
               CALL QUIT('Error reading file '//FNAME2)
            ELSE IF (IOSTATR.LT.0) THEN
               FOUND_DENS = .FALSE.
               BACKSPACE (LUTAM)
            ELSE
               READ(LUTAM) DENS3AB
               FOUND_DENS = .TRUE.
            END IF
C
            IF (.NOT.FOUND_DENS) THEN
               DTIME     = SECOND()
               CALL SO_DENS3(DENS3IJ,LDENSIJ,DENS3AB,LDENSAB,
     &                       T2AM, LT2AM, T2SAM, LT2SAM,
     &                       WORK,LWORK)
               DTIME     = SECOND()  - DTIME
               SOTIME(42) = SOTIME(42) + DTIME
C
               WRITE (LUTAM) DENS3IJ
               WRITE (LUTAM) DENS3AB
C
            ENDIF
C
            CALL GPCLOSE(LUTAM,'KEEP')
C
         ELSE
            CALL QUIT('TOPPA requires running the CC module')
         ENDIF
      ENDIF
C
C------------------------------------------------------
C     Replace T2-amplitudes in T2AM with
C     (2 * coulomb - exchange of T2-amplitudes).
C------------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_T2MP(T2AM,LT2AM)
      DTIME     = SECOND()  - DTIME
      SOTIME(5) = SOTIME(5) + DTIME
C
      IF (KEY.EQ.'MP3 ') THEN
         DTIME     = SECOND()
         CALL SO_T2MP(T2SAM,LT2SAM)
         DTIME     = SECOND()  - DTIME
         SOTIME(5) = SOTIME(5) + DTIME
      END IF
C
C-----------------------------------------------------
C     Send the amplitudes to the slaves
C-----------------------------------------------------
#ifdef VAR_MPI
      call soppa_update_amplitudes(t2am, lt2am, t2sam, lt2sam, key)
#endif
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_GETT2')
C
      RETURN

901   CALL QUIT('End of file when reading singles amplitudes')
902   CALL QUIT('End of file when reading doubles amplitudes')


      END

      SUBROUTINE SO_DUMP_DENSAI(KEY,DENSAI,LDENSAI)
C
C  Dumps the AI density to file
C
      use so_info, only : sop_dp
      implicit none
      character(len=4), intent(in) :: KEY
      real(sop_dp), intent(in) :: DENSAI(LDENSAI)
      integer, intent(in) :: LDENSAI
C
      character(len=8) :: FNAME
      INTEGER :: LUTAM, I, IDUMMY
C
C     Set filename from key
      IF (KEY .EQ. 'MP3 ') THEN
         FNAME = 'MP2__TAM'
      ELSE
         FNAME = KEY//'_TAM'
         IF (FNAME(4:4).EQ.' ') THEN
            FNAME(4:4) = '_'
         END IF
      END IF
      LUTAM = -1
C  GPOPEN SEEMS TO NOT ALLOW US TO SPECIFY OPEN FOR APPENDING,
      CALL GPOPEN(LUTAM,FNAME,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUTAM)
      READ(LUTAM) ! T2 amplitudes
      READ(LUTAM) ! D(ij)
      READ(LUTAM) ! D(ab)
C
      WRITE(LUTAM) (DENSAI(I), I = 1,LDENSAI)
      CALL GPCLOSE(LUTAM,'KEEP')
C
      END

      SUBROUTINE SO_READ_DENS(KEY,DENSIJ,LDENSIJ,
     &                        DENSAB,LDENSAB,DENSAI,LDENSAI)
C
      use so_info, only: sop_dp
      implicit none
      character(len=4), intent(in) :: KEY
      real(sop_dp), intent(out) :: DENSIJ(LDENSIJ),
     &                            DENSAB(LDENSAB), DENSAI(LDENSAI)
      integer, intent(in) :: LDENSIJ, LDENSAB, LDENSAI
C
      real(sop_dp), parameter :: HALF = 0.5D0
      character(len=8) :: FNAME
      INTEGER :: LUTAM, I, IDUMMY
C
C     Set filename from key
      IF (KEY .EQ. 'MP3 ') THEN
         FNAME = 'MP2__TAM'
      ELSE
         FNAME = KEY//'_TAM'
         IF (FNAME(4:4).EQ.' ') THEN
            FNAME(4:4) = '_'
         END IF
      END IF
      LUTAM = -1
C
      CALL GPOPEN(LUTAM,FNAME,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUTAM)

      IF (KEY(1:2).EQ.'CC') THEN
C    CC method - T1 amplitudes are first
         READ(LUTAM) DENSAI
         CALL DSCAL(LDENSAI,HALF,DENSAI,1)
      END IF
C     Skip T2 amplitudes
      READ(LUTAM)
C     Read densij and densab
      READ(LUTAM) DENSIJ
      READ(LUTAM) DENSAB
C
C     In case of MP2 and MP3: densai is last
      IF (KEY(1:2).EQ.'MP') THEN
         READ(LUTAM) DENSAI
      END IF
      CALL GPCLOSE(LUTAM,'KEEP')

      RETURN

      END
