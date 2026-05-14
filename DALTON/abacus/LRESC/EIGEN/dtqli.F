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

      FUNCTION pythagLRESC(a,b)
      DOUBLE PRECISION a,b,pythagLRESC
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythagLRESC=absa*sqrt(1.d0+(absb/absa)**2)
      else
        if(absb.eq.0.d0)then
          pythagLRESC=0.d0
        else
          pythagLRESC=absb*sqrt(1.d0+(absa/absb)**2)
        endif
      endif
      return
      END


C Tridiagonalization subroutine
      SUBROUTINE tqliLRESC(d,e,n,np,z)
      INTEGER n,np
      DOUBLE PRECISION, allocatable :: d(:),e(:),z(:,:)
      INTEGER i,iter,k,l,m
      DOUBLE PRECISION b,c,dd,f,g,p,r,s,pythagLRESC
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.d0
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.30)stop 'too many iterations in tqli'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.d0*e(l))
          r=pythagLRESC(g,1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.d0
          c=1.d0
          p=0.d0
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythagLRESC(f,g)
            e(i+1)=r
            if(r.eq.0.d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.d0*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
!     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
!     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.d0
          goto 1
        endif
15    continue
      return
      END
