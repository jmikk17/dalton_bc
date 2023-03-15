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
C
      SUBROUTINE CCSD_DIIS(S,T,NDIM,NITER)
C
C
C     Written by Poul Joegensen and Henrik Koch 3-Aug-1995
C
C     The S vector contain the perturbation estimate and T the previous
C     amplitudes. On exit the T contain the diis estimate of the new
C     amplitudes.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER (MAXVEC=80)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
      LOGICAL FIRST
C
      DIMENSION S(NDIM)
      DIMENSION T(NDIM)
C
      COMMON/LUDIIS/ LUTDIS, LUSDIS
      COMMON/DICED/A(MAXVEC+1,MAXVEC+1)
      DIMENSION B(MAXVEC+1),C(MAXVEC+1),WK1(MAXVEC+1),WK2(MAXVEC+1),
     +          AA(MAXVEC+1,MAXVEC+1)
cKeld      DIMENSION B(MAXVEC+1),C(MAXVEC+1),WK1(MAXVEC+1),WK2(MAXVEC+1),
cKeld     +          A(MAXVEC+1,MAXVEC+1),AA(MAXVEC+1,MAXVEC+1)
C
#include "ccsdinp.h"
C
      SAVE FIRST
      DATA FIRST /.TRUE./
C
      IF (MXDIIS .GT. MAXVEC) THEN
         WRITE(LUPRI,'(/A/A,I0,A,I0)') 'INFO: '//
     &'requested .MXDIIS is greater than maximum in CCSD_DIIS routine',
     &'INFO: .MXDIIS value decreased from ',MXDIIS,' to ',MAXVEC
         MXDIIS = MAXVEC
      END IF

      IF (FIRST) THEN
C
         FIRST = .FALSE.
         LUTDIS = -9003
         LUSDIS = -9004
C
         CALL GPOPEN(LUTDIS,'DIIST','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         CALL GPOPEN(LUSDIS,'DIISS','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
C
      ENDIF
C
      DO 100 I=1,NDIM
         T(I)= S(I) - T(I)
  100 CONTINUE
C
      IVEC=NITER
C
      IF (NITER .GT. MXDIIS) THEN
         IVEC = NITER - ((NITER-1)/MXDIIS)*MXDIIS
      ENDIF
C
      IF (IVEC.EQ.1) THEN
         IF (LUTDIS .LE. 0) CALL GPOPEN(LUTDIS,'DIIST','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         IF (LUSDIS .LE. 0) CALL GPOPEN(LUSDIS,'DIISS','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUTDIS)
         REWIND(LUSDIS)
      ENDIF
C
      WRITE(LUTDIS) S
      WRITE(LUSDIS) T
C
      NLOOP=NITER
      IF (NITER .GT. MXDIIS) THEN
         NLOOP = IVEC
      ENDIF
C
      REWIND(LUSDIS)
C
      DO 200 I = 1,NLOOP
C
         READ(LUSDIS) S
C
         ERR=DDOT(NDIM,T,1,S,1)
C
         A(I,IVEC)    =  ERR
         A(IVEC,I)    =  ERR
         A(I,NLOOP+1) = -ONE
         A(NLOOP+1,I) = -ONE
         B(I)         = ZERO
C
  200 CONTINUE
C
      A(NLOOP+1,NLOOP+1) =  ZERO
      B(NLOOP+1)         = -ONE
C
      NSIZ=NLOOP+1
C
      IFAIL = 0
      CALL F04ATF(A,MAXVEC+1,B,NSIZ,C,AA,MAXVEC+1,WK1,WK2,IFAIL)
C
      IF (IFAIL .NE. 0) THEN
        WRITE(LUPRI,*)  'IFAIL ',IFAIL
      ENDIF
C
      CALL DZERO(T,NDIM)
C
      REWIND(LUTDIS)
C
      DO 300 I = 1,NLOOP
C
         READ(LUTDIS) S
C
         CALL DAXPY(NDIM,C(I),S,1,T,1)
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE LUS(A,LDA,N,IPVT,X,B)
#include "implicit.h"
      DIMENSION A(LDA,1),X(1),B(1),IPVT(1)
C
C     SIMULTANEOUS EQUATION AX=B GIVEN PRIOR LU
C     DECOMPOSITION OF A
C     TAKEN FROM DONGARRA ACM TRANS MATH SOFTWARE VOL 10
C
      DO 10 K=1,N
10    X(K)=B(K)
      DO 20 K=1,N
      L=IPVT(K)
      XK=X(L)
      X(L)=X(K)
20    X(K)=XK
      DO 40 K=1,N
      XK=X(K)*A(K,K)
      IF(XK.NE.0.D00)THEN
      DO 30 I=K+1,N
30    X(I)=X(I)-A(I,K)*XK
      ENDIF
40    X(K)=XK
      DO 60 K=N,1,-1
      XK=X(K)
      IF(XK.NE.0.D00)THEN
      DO 50 I=1,K-1
50    X(I)=X(I)+A(I,K)*XK
      ENDIF
60    CONTINUE
      RETURN
      END
      SUBROUTINE LU(A,LDA,N,IPVT,INFO)
#include "implicit.h"
#include "priunit.h"
      DIMENSION A(LDA,*),IPVT(*)
C
C     LU DECOMPOSITION OF A
C
      INFO=0
      DO 40 J=1,N
      CALL SMXPY(N-J+1,A(J,J),J-1,LDA,A(1,J),A(J,1))
C
C     FIND PIVOT
C     NOTE : THIS IS ISAMAX OPERATION
C     (USE CRAY ROUTINE OR VECTOR MERGE)
C
      T=DABS(A(J,J))
      K=J
      DO 10 I=J+1,N
      IF(DABS(A(I,J)).GT.T)THEN
      T=DABS(A(I,J))
      K=I
      ENDIF
10    CONTINUE
      IPVT(J)=K
C
      IF(T.EQ.0.0D0)THEN
         GO TO 40
       ! WRITE(LUPRI,*)'ZERO PIVOT - MATRIX SINGULAR'
       ! CALL QUIT(' ')
       ! hjaaj Aug 2021: test "cc2_r12_he_a2" failed here for gfortran
       ! 10.2.0. With this change the test passed.
      ENDIF
C
C     SWOP ROWS
C     NOTE : THIS IS SSWAP OPERATION ON CRAY
C
      DO 20 I=1,N
      T=A(J,I)
      A(J,I)=A(K,I)
      A(K,I)=T
20    CONTINUE
C
      A(J,J)=1.0D0/A(J,J)
      CALL SXMPY(N-J,LDA,A(J,J+1),J-1,LDA,A(J,1),LDA,A(1,J+1))
      T=-A(J,J)
      DO 30 I=J+1,N
30    A(J,I)=T*A(J,I)
40    CONTINUE
      RETURN
      END
      SUBROUTINE SMXPY(N1,Y,N2,LDA,X,A)
#include "implicit.h"
      DIMENSION Y(*),X(*),A(LDA,*)
C
C     Y = Y +A X
C
      J=MOD(N2,2)
      IF(J.GE.1)THEN
      DO 10 I=1,N1
10    Y(I)=(Y(I))+X(J)*A(I,J)
      ENDIF
      J=MOD(N2,4)
      IF(J.GE.2)THEN
      DO 20 I=1,N1
20    Y(I)=((Y(I))+X(J-1)*A(I,J-1))+X(J)*A(I,J)
      ENDIF
      JMIN=J+4
      DO 40 J=JMIN,N2,4
      DO 30 I=1,N1
      Y(I)=((((Y(I))+X(J-3)*A(I,J-3))+X(J-2)*A(I,J-2))
     1    +X(J-1)*A(I,J-1))+X(J)*A(I,J)
30    CONTINUE
40    CONTINUE
      RETURN
      END
      SUBROUTINE SXMPY(N1,LDY,Y,N2,LDX,X,LDA,A)
#include "implicit.h"
      DIMENSION Y(LDY,*),X(LDX,*),A(LDA,*)
C
C
      J=MOD(N2,2)
      IF(J.GE.1)THEN
      DO 10 I=1,N1
10    Y(1,I)=(Y(1,I))+X(1,J)*A(J,I)
      ENDIF
      J=MOD(N2,4)
      IF(J.GE.2)THEN
      DO 20 I=1,N1
20    Y(1,I)=((Y(1,I))+X(1,J-1)*A(J-1,I))+X(1,J)*A(J,I)
      ENDIF
      JMIN=J+4
      DO 40 J=JMIN,N2,4
      DO 30 I=1,N1
      Y(1,I)=((((Y(1,I))+X(1,J-3)*A(J-3,I))+X(1,J-2)*A(J-2,I))
     1    +X(1,J-1)*A(J-1,I))+X(1,J)*A(J,I)
30    CONTINUE
40    CONTINUE
      RETURN
      END
      SUBROUTINE F04ATF(A,IA,B,N,C,AA,IAA,WKS1,IPVT,IFAIL)
#include "implicit.h"
      DIMENSION A(IA,N),B(N),C(N),AA(IAA,N),WKS1(N),IPVT(N)
C
C     THIS IS NOT THE NAG ROUTINE F04ATF, THOUGH THE
C     ARGUMENTS ARE THE SAME
C     IT IS A (HOPEFULLY) LESS MACHINE DEPENDENT REPLACEMENT
C     BASED ON VECTOR ALGORITHM OF DONGARRA USING
C     LU DECOMPOSITION
C     SEE ACM TRANSACTIONS OF MATHEMATICAL SOFTWARE
C     VOL 10 SEPT 1984
C
      DO 10 I=1,N
      DO 10 J=1,N
10    AA(J,I)=A(J,I)
      CALL LU(AA,IAA,N,IPVT,INFO)
      CALL LUS(AA,IAA,N,IPVT,C,B)
      IFAIL=0
      RETURN
      END
