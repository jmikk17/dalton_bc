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
      subroutine diagonalizar(matrix,nmatrix,autovalores,autoestados)


      INTERFACE
      SUBROUTINE tred2LRESC(a,n,np,d,e)
      INTEGER n,np
      DOUBLE PRECISION, allocatable :: a(:,:),d(:),e(:)  
      INTEGER i,j,k,l
      DOUBLE PRECISION f,g,h,hh,scale
      END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE tqliLRESC(d,e,n,np,z)
      INTEGER n,np
      DOUBLE PRECISION, allocatable :: d(:),e(:),z(:,:)
      INTEGER i,iter,k,l,m
      DOUBLE PRECISION b,c,dd,f,g,p,r,s,pythagLRESC
      END SUBROUTINE
      END INTERFACE

      INTERFACE
      SUBROUTINE eigsrtLRESC(d,v,n,np)
      INTEGER n,np
      DOUBLE PRECISION, allocatable :: d(:),v(:,:)
      INTEGER i,j,k
      REAL p
      END SUBROUTINE
      END INTERFACE

       DOUBLE PRECISION, allocatable :: a(:,:),d(:),e(:)
       DOUBLE PRECISION, intent(in) :: matrix(nmatrix,nmatrix)
       DOUBLE PRECISION, intent(out) :: autovalores(nmatrix)
       DOUBLE PRECISION, intent(out) :: autoestados(nmatrix,nmatrix)
       INTEGER, intent(in) :: nmatrix
       INTEGER n, np
       n=nmatrix

        allocate(a(nmatrix,nmatrix))
        allocate(e(nmatrix))
        allocate(d(nmatrix))

        a=matrix
        d=0.d0
        e=0.d0
c       A is the matrix to tridiagonalize. D the resultant diagonal and e
c       a secondary array.
        call tred2LRESC(a,nmatrix,nmatrix,d,e)
c       A is "Q", as its the tred2 output. d-out the eigenvalues and e-out is useless.
        call tqliLRESC(d,e,nmatrix,nmatrix,a)
        call eigsrtLRESC(d,a,nmatrix,nmatrix)

c       eigenvalues are ordered from smallest to largest
        autovalores=d
c       A is now the matrix whose columns are the eigenstates
        autoestados=a
        deallocate (a)
        deallocate (d)
        deallocate (e)

        end 





