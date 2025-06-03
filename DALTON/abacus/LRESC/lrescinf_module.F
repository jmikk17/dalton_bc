!
!  Dalton, a molecular electronic structure program
!  Copyright (C) The Dalton Authors (see AUTHORS file for details).
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
! Module: lrescinf_module
      MODULE lrescinf_module
      IMPLICIT NONE

!--------------------------------------------------------
!    1. Default access is PRIVATE
!--------------------------------------------------------
      PRIVATE

!--------------------------------------------------------
! 2. Public constants (read-only)
!--------------------------------------------------------
      PUBLIC :: calfa, CPNR, CDNR, CFCZK, CFCDE, CSDK, CSDD
      PUBLIC :: CFCBS, CSDBS, CPSOK, CPSOKaba, COZK, CSDAV
      PUBLIC :: CDIAM, CDIAD, CDIAK, CANGP
      PUBLIC :: CFCAV, CFCAV0, CFCAV1, CFCAV2, CFCAV3
      PUBLIC :: CPSOOZMV, CPSOOZDW, COZFCSO, COZSDSO
      PUBLIC :: CEFGDW, CEFGMV, CEFGpqp, CEFGkin, CEFGlap
      PUBLIC :: C4EFGSO, C4EFGSOL, C4EFGMvpqp, C4EFGDwpqp

      REAL(8), PARAMETER :: calfa = 1.0 / 137.036

      REAL(8), PARAMETER :: CPNR = -0.5D0
      REAL(8), PARAMETER :: CDNR = 1.0D0
      REAL(8), PARAMETER :: CFCZK = -3.0/8.0
      REAL(8), PARAMETER :: CFCDE = 0.25D0
      REAL(8), PARAMETER :: CSDK = -0.375D0
      REAL(8), PARAMETER :: CSDD = 1.0 / 4.0
      REAL(8), PARAMETER :: CFCBS = 0.25D0
      REAL(8), PARAMETER :: CSDBS = 0.25D0
      REAL(8), PARAMETER :: CPSOK = -1.0 / 8.0
      REAL(8), PARAMETER :: CPSOKaba = -0.5D0
      REAL(8), PARAMETER :: COZK = -1.0D0
      REAL(8), PARAMETER :: CSDAV = -3.0 / 8.0
      REAL(8), PARAMETER :: CDIAM = 1.0D0
      REAL(8), PARAMETER :: CDIAD = 1.0D0
      REAL(8), PARAMETER :: CDIAK = 1.0 / 6.0
      REAL(8), PARAMETER :: CANGP = -0.5D0
      REAL(8), PARAMETER :: CFCAV = -7.0 / 16.0
      REAL(8), PARAMETER :: CFCAV0 = -3.0 / 8.0
      REAL(8), PARAMETER :: CFCAV1 = 1.0 / 8.0
      REAL(8), PARAMETER :: CFCAV2 = 1.0 / 4.0
      REAL(8), PARAMETER :: CFCAV3 = -3.0 / 8.0
      REAL(8), PARAMETER :: CPSOOZMV = -0.5D0
      REAL(8), PARAMETER :: CPSOOZDW = -0.5D0
      REAL(8), PARAMETER :: COZFCSO = -0.25D0
      REAL(8), PARAMETER :: COZSDSO = -0.25D0

      REAL(8), PARAMETER :: CEFGDW = 1.0D0
      REAL(8), PARAMETER :: CEFGMV = 1.0D0
      REAL(8), PARAMETER :: CEFGpqp = 0.25D0
      REAL(8), PARAMETER :: CEFGkin = -0.5D0
      REAL(8), PARAMETER :: CEFGlap = 1.0 / 8.0
      REAL(8), PARAMETER :: C4EFGSO = 0.25D0
      REAL(8), PARAMETER :: C4EFGSOL = 1.0D0
      REAL(8), PARAMETER :: C4EFGMvpqp = 0.25D0
      REAL(8), PARAMETER :: C4EFGDwpqp = 0.25D0

!--------------------------------------------------------
! 3. Public variables (read-write access)
!--------------------------------------------------------      
      PUBLIC :: JIMPRT, LRATOM, JJAPRT, NLRATOMS, LRATOMSS
      PUBLIC :: SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S
      PUBLIC :: SIGMAP2S, SIGMAP2T, LRESCALL, GAUCHANG, PRTALL1
      PUBLIC :: LRATOM_changed, LRANISO, SIGMAPNR, SIGMADNR
      PUBLIC :: PYINTS, EXTRAINFO, WARNLRS, WARNLRT, WARNQR, LRINTS
      PUBLIC :: ORBCON, PQPKINLRESC, PRTALL2, EFGC2FLAG, EFGC4FLAG
      PUBLIC :: LRFCAV, LRFCAVpy, LRFCAV0, LRFCAV1, LRFCAV2, LRFCAV3
      PUBLIC :: LRPNR, LRDNR, SGDNR, SGPNR
      PUBLIC :: SGD0S, SGD1S, SGP1S, SGP1T, SGP3S, SGP3T
      PUBLIC :: SLRESC, LRSDAV
      PUBLIC :: LRDIAK, LRANGP, LRANGPavg
      PUBLIC :: LRDIAM, LRDIAD, LROZK
      PUBLIC :: LRPSOK, LRPSKI, LRFCZK
      PUBLIC :: LRPSOKaba, LRFCDE, LRSDK, LRSDD
      PUBLIC :: LRSDZK, LRFCBS, LRSDBS
      PUBLIC :: QRPSOOZMV, QRPSOOZDW, QROZFCSO, QROZSDSO
      PUBLIC :: LRGAUG, EFGC0, EFGC2, EFGC4, EFGTIME, LRNUCEXP


      INTEGER :: JIMPRT, LRATOM, JJAPRT, NLRATOMS, LRATOMSS(100)
      LOGICAL :: SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S
      LOGICAL :: SIGMAP2S, SIGMAP2T, LRESCALL, GAUCHANG, PRTALL1
      LOGICAL :: LRATOM_changed, LRANISO, SIGMAPNR, SIGMADNR
      LOGICAL :: PYINTS, EXTRAINFO, WARNLRS, WARNLRT, WARNQR, LRINTS
      LOGICAL :: ORBCON, PQPKINLRESC, PRTALL2, EFGC2FLAG, EFGC4FLAG


      DOUBLE PRECISION :: LRFCAV(3, 3), LRFCAVpy(3,3)
      DOUBLE PRECISION :: LRFCAV0(3,3), LRFCAV1(3,3), LRFCAV2(3,3)
      DOUBLE PRECISION :: LRFCAV3(3,3)
      DOUBLE PRECISION :: LRPNR(3,3),LRDNR(3,3),SGDNR(3,3), SGPNR(3,3)
      DOUBLE PRECISION :: SGD0S(3,3), SGD1S(3,3), SGP1S(3,3), SGP1T(3,3)
      DOUBLE PRECISION :: SGP3S(3,3), SGP3T(3,3)
      DOUBLE PRECISION :: SLRESC(3,3)
      DOUBLE PRECISION :: LRSDAV(3,3)
      DOUBLE PRECISION :: LRDIAK(3, 3), LRANGP(3, 3), LRANGPavg(3,3)
      DOUBLE PRECISION :: LRDIAM(3, 3), LRDIAD(3, 3), LROZK(3, 3)
      DOUBLE PRECISION :: LRPSOK(3, 3), LRPSKI(3, 3), LRFCZK(3, 3)
      DOUBLE PRECISION :: LRPSOKaba(3,3)
      DOUBLE PRECISION :: LRFCDE(3,3), LRSDK(3,3), LRSDD(3,3)
      DOUBLE PRECISION :: LRSDZK(3, 3), LRFCBS(3, 3), LRSDBS(3, 3)
      DOUBLE PRECISION :: QRPSOOZMV(3,3), QRPSOOZDW(3,3)
      DOUBLE PRECISION :: QROZFCSO(3,3), QROZSDSO(3,3)

      DOUBLE PRECISION :: LRGAUG(3), EFGC0(100), EFGC2(100, 6)
      DOUBLE PRECISION :: EFGC4(100,7), EFGTIME(5,2), LRNUCEXP

      !CONTAINS

      ! Functions and/or subroutines may be added here if neccesary.

      END MODULE lrescinf_module

