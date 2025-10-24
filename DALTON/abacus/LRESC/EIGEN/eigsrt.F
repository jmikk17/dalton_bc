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

      SUBROUTINE eigsrtLRESC(d,v,n,np)
      INTEGER n,np
      DOUBLE PRECISION, allocatable :: d(:),v(:,:)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
