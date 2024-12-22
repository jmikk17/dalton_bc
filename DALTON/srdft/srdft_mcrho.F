
#ifdef JKP_DEBUG
C*****************************************************************************
      SUBROUTINE MCRHO(CMO,INDXCI,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Jan. 2003
C
C     Purpose : Generate active density from multiconfigurationel
C               wavefunction on a grid.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "inforb.h"
#include "infvar.h"
#include "dftinf.h"
#include "dftcom.h"

      LOGICAL DOGGA,DOBCK,DOLYP,DOLND,DO_MOLGRAD,DODCAO
      DIMENSION CMO(*),WORK(LWORK)
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C     =================================
C     Initializations :
C     Number of AOs and their addresses
C     =================================
      DOLND = .FALSE.
      DOBCK = .FALSE.
      DOLYP = .FALSE.
      DO_MOLGRAD = .FALSE.
      DOGGA = DOBCK .OR. DOLYP
      
      NDER  = 0
      IF (DO_MOLGRAD) NDER = 1
      IF (DOGGA) THEN
         NDER = NDER + 1
         IF (DFTPOT) NDER = NDER + 1
      END IF

      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1= (NSOB1-1)*NBAST + 1
C     =============================
C     Allocations
C     =============================
      DODCAO = (NISHT.GT.0)
      NBUF   = 1000000
      KCREF  = 1
      KX     = KCREF + NCONF
      KY     = KX    + NBUF
      KZ     = KY    + NBUF
      KW     = KZ    + NBUF
      KGSO   = KW    + NBUF
      KCNT   = KGSO  + NTYPSO*NBAST
      KDGA   = KCNT  + NBAST
      KDV    = KDGA  + NBAST
      KDVAO  = KDV   + NNASHX
      KDCAO  = KDVAO + N2BASX
      IF (DODCAO) THEN
          KWRK  = KDCAO  + N2BASX
      ELSE
          KWRK  = KDCAO
      ENDIF
      LWRK  = LWORK + 1 - KWRK
      CALL MCRHO1(CMO,WORK(KDV),DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGA),DOGGA,WORK(KCREF),INDXCI,DODCAO,
     &            WORK(KDCAO),WORK(KDVAO),WORK(KWRK),LWRK,IPRINT)

      CALL TIMER('MCRHO',TIMSTR,TIMEND)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE MCRHO1(CMO,DV,DOLND,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,
     &                  NCNT,DMATGAO,DOGGA,CREF,INDXCI,DODCAO,DCAO,DVAO,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Jan. 2003
C
C     Purpose : Generate density from multiconfigurationel
C               wavefunction on a grid.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "infpar.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
      LOGICAL   DOLND,DOGGA,DODCAO
      DIMENSION CMO(*),DV(NNASHX),CORX(NBUF),CORY(NBUF),
     &          CORZ(NBUF),WEIGHT(NBUF),GSO(NBAST*NTYPSO),
     &          NCNT(NBAST),DMATGAO(NBAST),CREF(NCONF),
     &          DCAO(NBAST,NBAST),DVAO(NBAST,NBAST),WORK(LWORK), COR(3)
      CHARACTER*17 QUADNAME
#include "chrnos.h"

C     =============================
C     Open file for dumping density
C     =============================
      CALL GPOPEN(LUDDUMP,'DENSITY.GRID','UNKNOWN',
     &                ' ','FORMATTED',IDUMMY,.FALSE.)

C     =============================
C     Retrieve CREF
C     =============================
C
      REWIND LUIT1
      CALL MOLLAB('STARTVEC ',LUIT1,LUERR)
      DO 100 I = 1, (ISTACI-1)
         READ (LUIT1)
  100 CONTINUE
      CALL READT(LUIT1,NCONF,CREF)

C     =============================
C     Calculate DV
C     =============================
C
      CALL MAKDV(CREF,DV,INDXCI,WORK,LWORK)

C     =============================
C     Calculate inactive and active
C     density matrises in AO-basis
C     =============================
C
      CALL FCKDEN(DODCAO,.TRUE.,DCAO,DVAO,CMO,DV,WORK,LWORK)
C
C     =============================
C     Calculate grid
C     =============================

      IF (.NOT.srDFTGRID_done) THEN
         CALL MAKE_srDFTGRID(WORK,LWORK,NTOT_DFTGRID,1,.FALSE.)
         CALL CONDFT
         srDFTGRID_done = .TRUE.
      END IF

      CELCTRN = D0
      VELCTRN = D0

      LUQUAD = -1

C     Make quadname: Can take 10000 procs

      IF (MYNUM .GT. 0) THEN
         QUADNAME = 'DALTON.QUAD.n'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
C
         CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)
      ELSE
         CALL GPOPEN(LUQUAD,'DALTON.QUAD','OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      END IF

      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL REAQUA_srdft(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
#if defined (VAR_SRDFT_MPI)
         DO 300 IPNT = 1+MYNUM, NPOINT, NODTOT+1
#else
         DO 300 IPNT = 1, NPOINT
#endif
            IF (IPRINT .GT. 100) THEN
               WRITE (LUPRI,'(2X,I6,4F12.6)')
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF

            WGHT = WEIGHT(IPNT)

C           AOs
C           ===

            THRINT = DFTHRI/WGHT
            COR(1) = CORX(IPNT)
            COR(2) = CORY(IPNT)
            COR(3) = CORZ(IPNT)
            CALL GETSOS(GSO,COR,WORK,LWORK,
     &                  NBAST,DOLND,DOGGA,THRINT,IPRINT)

C           Density & Number of electrons
C           =============================

C           ... inactive part
            IF (DODCAO) THEN
               CALL GETRHO_srdft(DCAO,GSO,RHO,DMATGAO,THRINT,IPRINT)
               CELCTRN = CELCTRN + WGHT*RHO
            ENDIF
C           ... active part
            CALL GETRHO_srdft(DVAO,GSO,RHO,DMATGAO,THRINT,IPRINT)
            VELCTRN = VELCTRN + WGHT*RHO

C           Dump active density to file
C           =============================

            WRITE(LUDDUMP,'(4F25.13)') CORX(IPNT),CORY(IPNT),
     &           CORZ(IPNT),RHO

  300    CONTINUE

         GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF

      CALL GPCLOSE(LUQUAD,'KEEP')

C     Test on the number of electrons

      ELCTRX = FLOAT(2*NRHFT)
      ELCTRN = CELCTRN + VELCTRN
      ERROR  = ELCTRN - ELCTRX

C     Print section

      WRITE (LUPRI,'(3(/2X,A,F14.6))')
     &' Nr. of inactive electrons from numerical integration :',CELCTRN,
     &' Nr. of active electrons from numerical integration   :',VELCTRN,
     &' Total nr. of electrons from numerical integration    :',ELCTRN
      WRITE (LUPRI,'(2(/2X,A,F14.6))')
     &' Number of electrons from orbial occupations          :',ELCTRX,
     &' Error in the number of electrons                     :',ERROR
      IF (ABS(ERROR) .GT. DFTELS) THEN
         WRITE (LUPRI,'(/2X,A,F14.6,/2X,A)')
     &' Error larger than DFTELS (set input)                 :',DFTELS,
     &   ' Calculation aborted.'
         CALL QUIT
     &    ('Wrong number of electrons in DFTDRV. Calculation aborted.')
      END IF

      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(/2X,A,F14.7,6X,D8.2,I14)')
     &      ' Number of electrons/abscissas:  ',
     &        ELCTRN,ELCTRN-ELCTRX,NPNTS
      END IF
      RETURN
      END
#endif  /* JKP_DEBUG */
