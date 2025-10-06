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
!  code, you can obtain one at
!  https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
C

!     As the subroutine duplicates the matrix, its use is recommended
!     just for small matrices.
!      subroutine diag_lresc(matrix,nmatrix,eigenvalues,autoestados)


 !     INTERFACE
 !     SUBROUTINE tred2LRESC(a,n,np,d,e)
 !     INTEGER n,np
 !     DOUBLE PRECISION, allocatable :: a(:,:),d(:),e(:)
 !     INTEGER i,j,k,l
 !     DOUBLE PRECISION f,g,h,hh,scale
 !     END SUBROUTINE
 !     END INTERFACE

 !     INTERFACE
 !      SUBROUTINE tqliLRESC(d,e,n,np,z)
 !     INTEGER n,np
 !     DOUBLE PRECISION, allocatable :: d(:),e(:),z(:,:)
 !     INTEGER i,iter,k,l,m
 !!     DOUBLE PRECISION b,c,dd,f,g,p,r,s,pythagLRESC
 !     END SUBROUTINE
 !     END INTERFACE

 !     INTERFACE
 !     SUBROUTINE eigsrtLRESC(d,v,n,np)
 !     INTEGER n,np
 !     DOUBLE PRECISION, allocatable :: d(:),v(:,:)
 !     INTEGER i,j,k
 !     REAL p
 !     END SUBROUTINE
 !     END INTERFACE

 !      DOUBLE PRECISION, allocatable :: a(:,:),d(:),e(:)
 !      DOUBLE PRECISION, intent(in) :: matrix(nmatrix,nmatrix)
 !      DOUBLE PRECISION, intent(out) :: eigenvalues(nmatrix)
 !      DOUBLE PRECISION, intent(out) :: autoestados(nmatrix,nmatrix)
 !      INTEGER, intent(in) :: nmatrix
 !      INTEGER n, np
 !      n=nmatrix

 !       allocate(a(nmatrix,nmatrix))
 !       allocate(e(nmatrix))
 !       allocate(d(nmatrix))

!        a=matrix
!        d=0.d0
!        e=0.d0
c       A is the matrix to tridiagonalize. D the resultant diagonal and e
c       a secondary array.
!        call tred2LRESC(a,nmatrix,nmatrix,d,e)
c       A is "Q", as its the tred2 output. d-out the eigenvalues and e-out is useless.
!        call tqliLRESC(d,e,nmatrix,nmatrix,a)
!        call eigsrtLRESC(d,a,nmatrix,nmatrix)

c       eigenvalues are ordered from smallest to largest
!        eigenvalues=d
c       A is now the matrix whose columns are the eigenstates
!        autoestados=a
!        deallocate (a)
!        deallocate (d)
!        deallocate (e)

!        end




      subroutine shielddiag(matrix,n,parameters,principal_axes,eigenv)
C     Constructed according to DIRAC23: prp/pamprp.F
C     Parameters calculated in accordance with J.Mason, Solid State
C     Nuc.Magn.Resonance, 2 (1993), 285
C     Tensor treatment according to F. A. Anet et. alter, Concepts in
C     Magnetic Resonance, 1991, vol. 3, no 4.
#include "implicit.h"
#include "dummy.h"
#include "wrkrsp.h"
      DOUBLE PRECISION, intent(in) :: matrix(n,n)
      DOUBLE PRECISION, intent(out) :: parameters(4)
      DOUBLE PRECISION, intent(out) :: principal_axes(3,3)
      DOUBLE PRECISION, intent(out) :: eigenv(3)
      DOUBLE PRECISION :: ISO
      INTEGER, intent(in) :: n

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0)
      DIMENSION PVAL(6), WRK(3), IWRK(3), AXES(3,3)
      DIMENSION SYM(3,3), ANTI(3,3)
      PARAMETER (thresh = 1D-06)

!     Calculate ANTI and SYM matrices
      DO I = 1, 3
      DO J = 1, 3
         ANTI(I, J) = (matrix(I, J) - matrix(J, I)) / D2
         SYM(I, J) = (matrix(I, J) + matrix(J, I)) / D2
      ENDDO
      ENDDO

      CALL DUNIT(AXES, 3)

!     Fill PVAL vector with elements from the symmetric matrix SYM
      IJ = 1
      DO I = 1, 3
      DO J = 1, I
         PVAL(IJ) = SYM(I, J)
         IJ = IJ + 1
      ENDDO
      ENDDO

!     Diagonalize SYM matrix
      CALL JACO(PVAL, AXES, 3, 3, 3, WRK, IWRK)

!     Rearrange PVAL elements
      PVAL(1) = PVAL(1)
      PVAL(2) = PVAL(3)
      PVAL(3) = PVAL(6)

!     Find maximum, medium and minimum indices in PVAL
      IZ = IDMAX(3, PVAL, 1)
      IX = IDMIN(3, PVAL, 1)
      DO I = 1, 3
         IF (I /= IZ .AND. I /= IX) IY = I
      ENDDO
      principal_axes(:,1) = AXES(:,IX)  ! min eigenvalue axis
      principal_axes(:,2) = AXES(:,IY)  ! medium
      principal_axes(:,3) = AXES(:,IZ)  ! max eigenvalue axis

!     Calculate parameters
      ISO  = ( PVAL(IX) + PVAL(IY) + PVAL(IZ) ) / 3
      SPAN = PVAL(IZ) - PVAL(IX)
      SKEW = D3 * (ISO - PVAL(IY)) / SPAN
      ANIS = PVAL(IZ) - (PVAL(IX) + PVAL(IY)) / D2
      ASYM = (PVAL(IY) - PVAL(IX)) / (PVAL(IZ) - ISO)

!     Adjust values to avoid numerical instabilities
      IF (ABS(SPAN) < thresh) THEN
         SPAN = D0
         SKEW = D0
         ANIS = D0
         ASYM = D1
      ELSE
         IF (ABS(PVAL(IY) - PVAL(IX)) .LT. thresh) SKEW = D1
         IF (ABS(PVAL(IZ) - PVAL(IY)) .LT. thresh) SKEW = -D1
      ENDIF
      IF (ABS((PVAL(IZ) - iso)) .LT. thresh) ASYM = D1

c     eigenvalues are ordered from smallest to largest (in case needed)
      eigenv = [PVAL(IX), PVAL(IY), PVAL(IZ)]
      parameters = [SPAN, SKEW, ANIS, ASYM]

      end 
