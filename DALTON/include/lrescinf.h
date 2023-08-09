!
!     File: lrescinf.h
!     Purpose: Control of what to do in LRESC module
!
!     NOTE:
!
!cx         jim-gesc -1: RNLRSC login included in abainf.h, jimprt for debugging prints
!cx         jja  : edited on April 2021 to add variables for EFG
!cx         jja  : edited on February 2023 to add and modify variables for shielding
      REAL*8 calfa, CPNR, CDNR, CFCZK, CFCDE, CSDK, CSDD, CSDBS, CPSOK,
     &        CFCBS, COZK, CDIAM, CDIAD, CDIAK, CSDAV, CANGP,CFCAV, CPSOOZMV, 
     &        CPSOOZDW, CEFGDW, CEFGMV, CEFGpqp, CEFGkin, C4EFGSO, C4EFGMvpqp,
     &        CEFGlap, C4EFGDwpqp

      PARAMETER (calfa=1.0/137.036, CFCZK=-3.0/8.0, CFCDE=0.25D0,
     &            CSDK=-0.375D0, CSDD=1/4.D0, CPNR=-0.5D0, CDNR=1.0D0,
     &            CFCBS=0.25D0, CSDBS=0.25D0, 
     &            CPSOK=-0.125D0, COZK=-1.0D0, CDIAM=1.0D0, 
     &            CDIAD=1.0D0, CDIAK=1/6.0, CANGP=-0.5D0,
     &            CFCAV=-7.0/16.0, CSDAV=-1.0/4.0, CPSOOZMV=-0.5D0,
     &            CPSOOZDW=-0.5D0, COZFCSO=-0.25D0, COZSDSO=-0.25D0, 
     &            CEFGDW=-1.0D0, CEFGMV=-1.0D0, CEFGpqp=0.25D0, 
     &            CEFGkin=-0.5D0, CEFGlap=1.0/8.0, C4EFGSO=0.25D0,
     &            C4EFGMvpqp=0.25D0, C4EFGDwpqp=0.25D0)

      INTEGER JIMPRT, LRATOM, JJAPRT

    !LRESC corrections to shielding
      LOGICAL SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &         SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG,
     &         PRTALL1, LRATOM_changed, LRANISO, SIGMAPNR, SIGMADNR

    !LRESC corrections to EFG and common logicals
      LOGICAL ORBCON,PQPKINLRESC,PRTALL2,,EFGC2FLAG,EFGC4FLAG

      DOUBLE PRECISION LRFCAV(3,3), LRDIAK(3,3), LRANGP(3,3), 
     &                  LRSDAV(3,3), LRDIAM(3,3), LRDIAD(3,3),
     &                  LROZK(3,3), LRPSOK(3,3), LRPSKI(3,3),
     &                  LRFCZK(3,3),
     &                  LRFCDE(3,3), LRSDK(3,3),
     &                  LRSDD(3,3), LRFCBS(3,3), 
     &                  LRSDBS(3,3), LRPNR(3,3),LRDNR(3,3),SGDNR(3,3),
     &			SLRESC(3,3), SGPNR(3,3), SGD0S(3,3), SGD1S(3,3),
     &                  SGP1S(3,3), SGP1T(3,3), SGP3S(3,3),SGP3T(3,3),
     &                  QRPSOOZMV(3,3), QRPSOOZDW(3,3),
     &                  QROZFCSO(3,3), QROZSDSO(3,3),
     &                  LRGAUG(3),EFGC0(100),EFGC2(100,5),LRNUCEXP,
     &                  EFGC4(100,6)

      COMMON /LRESCINF/ SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &    SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG, JIMPRT, LRATOM,
     &    LRFCAV, LRDIAK, LRANGP, LRSDAV, LRDIAM, LRDIAD, LROZK,
     &    LRPSOK, LRPSKI, LRFCZK, LRFCDE, LRSDK, LRSDD, SGDNR,
     &    LRFCBS, LRSDBS, LRPNR, LRDNR, SLRESC, SGPNR, SGD0S, 
     &    SGD1S, SGP1S, SGP1T, SGP3S, SGP3T, QRPSOOZMV, QRPSOOZDW,
     &    QROZFCSO, QROZSDSO, LRGAUG, PRTALL1, LRANISO,
     &    EFGC0,EFGC2, ORBCON,PQPKINLRESC, PRTALL2, JJAPRT,
     &    LRNUCEXP, LRATOM_changed, SIGMAPNR, SIGMADNR,EFGC2FLAG,
     &    EFGC4FLAG, EFGC4
! -- end of lrescinf.h --
