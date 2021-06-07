!
!     File: lrescinf.h
!     Purpose: Control of what to do in LRESC module
!
!     NOTE:
!
!cx         jim-gesc -1: RNLRSC login included in abainf.h, jimprt for debugging prints
!cx         jja  : edited on April 2021 to add variables for EFG
      REAL*8 calfa, CFCZK, CSDZK, CFCBS, CSDBS, CPSOK, CLKIN, CDIAM,
     &       CDIAD, CDIAK, CANGP,CFCAV, CEFGDW,
     &       CEFGMV, CEFGpqp, CEFGkin, C4EFGSO, C4EFGMV, C4EFGDW

      PARAMETER (calfa=1.0/137.036, CFCZK=1.0/3.0, CSDZK=-0.25D0,
     &            CFCBS=-0.25D0, CSDBS=-0.25D0, CPSOK=0.50D0,
     &            CLKIN=1.0D0, CDIAM=-1.0D0, CDIAD=-1.0D0,
     &            CDIAK=10.0/137.036, CANGP=-0.5D0,
     &            CFCAV=-1.0*7.0/16.0, CEFGDW=-1.0D0, 
     &            CEFGMV=-1.0D0, CEFGpqp=0.25D0, CEFGkin=-0.5D0,
     &		  CEFGlap=1.0/8.0)

      INTEGER JIMPRT, LRATOM, JJAPRT

    !LRESC corrections to shielding
      LOGICAL SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &         SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG,
     &         PRTALL1

    !LRESC corrections to EFG and common logicals
      LOGICAL ORBCON,PQPKINLRESC,PRTALL2

      DOUBLE PRECISION LRFCAV(3,3), LRDIAK(3,3), LRANGP(3,3),
     &                  LRDIAM(3,3), LRDIAD(3,3), LRLKIN(3,3),
     &                  LRPSOK(3,3), LRPSKI(3,3), LRFCZK(3,3),
     &                  LRSDZK(3,3), LRFCBS(3,3), LRSDBS(3,3),
     &                  LRGAUG(3),EFGC0(100),EFGC2(100,5)

      COMMON /LRESCINF/ SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &    SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG, JIMPRT, LRATOM,
     &    LRFCAV, LRDIAK, LRANGP, LRDIAM, LRDIAD, LRLKIN, LRPSOK,
     &    LRPSKI, LRFCZK, LRSDZK, LRFCBS, LRSDBS, LRGAUG, PRTALL1,
     &    EFGC0,EFGC2, ORBCON,PQPKINLRESC, PRTALL2, JJAPRT
! -- end of lrescinf.h --
