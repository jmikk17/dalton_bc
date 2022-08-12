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

      SUBROUTINE EFGdrv(xmos,PRT,cmos,OVERLP)
C
C  Based on AVG program              
C  Melo, azua, giribet, alejandrito, 
C  G. Aucar, Danian, romero, homero,bart y lisa

C  Revisited and modified: J. Aucar 2020/2021
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inforb.h"
#include "nuclei.h"
#include "inftap.h"
#include "lrescinf.h"
#include "chrnos.h"


      CHARACTER*8 B0
      integer PRT
      CHARACTER*8  STARS
      CHARACTER*8 RTNLBL(2)
      PARAMETER (STARS = '********')
      double precision, dimension(NBAST*NORBT), intent(in) :: xmos
      double precision, dimension(NBAST,NBAST), intent(inout) :: OVERLP
      double precision cmos(NORBT,NBAST)
      CHARACTER*8 :: ZZEFG(NUCIND),GZZEFG(3*NUCIND)
      double precision :: ORBC(NOCCT),Av,Aux
      
      double precision, dimension(NBAST,NBAST) :: aux1,aux2,aux3
      double precision, dimension(NBAST,NBAST) :: aux0
      double precision, dimension(NBAST,NBAST) :: XDIPVEL
      double precision, dimension(NBAST,NBAST) :: YDIPVEL,ZDIPVEL
      double precision, dimension(NBAST,NBAST) :: GZZEFGx,GZZEFGy
      double precision, dimension(NBAST,NBAST) :: GZZEFGz
      double precision, dimension(NBAST,NBAST) :: X1SPNORB
      double precision, dimension(NBAST,NBAST) :: Y1SPNORB,Z1SPNORB
      double precision, dimension(NBAST,NBAST) :: ZZEFGm


C ------------------------------------------------------
C       DALTON.CM variables
C ------------------------------------------------------
       naos = NBAST ! NBAST: Total number of (atomic) basis functions
       nmos= NORBT !NORBT: Total number of molecular orbitals
       !Some other variables (por possible future use):
       !NRHFT:  Occupied SCF orbitals (of closed shells¿?)
       !NOCCT: Total number of occupied orbitals
       !NASHT: Open shell SCF orbitals. See on sirinp.F.
C      cmos(i,j)   :  Molecular orbitals coeficients matrix elements (from SIRIUS.RST)

       k = 1

       cmos=0.0
       do i=1,nmos                  !i corresponds to a MO
          do j=1,naos				   !j corresponds to j coef. of the ith MO
             cmos(i,j) = xmos(k)		! Transforms vector into matrix
             k = k + 1
          end do
       end do

       
        !write(lupri,*) 'CMOS (avgdrv):' 
        !call OUTPUT(CMOS,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)  ! prints CMOS

       IF (ORBCON) B0='ORBCON  '
       call READAOP(OVERLP, naos, 'OVERLAP ') !gets OVERLAP matrix
!      write(lupri,*) 'OVERLAP (avgdrv):'
!      call OUTPUT(OVERLP,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints OVERLP por debugging

C     Get the operators (its representation in the AO basis)
      call READAOP(XDIPVEL, naos, 'XDIPVEL ') !gets XDIPVEL matrix
      call READAOP(YDIPVEL, naos, 'YDIPVEL ') !gets YDIPVEL matrix
      call READAOP(ZDIPVEL, naos, 'ZDIPVEL ') !gets ZDIPVEL matrix
      call READAOP(X1SPNORB, naos,'X1SPNORB') !gets X1SPNORB matrix
      call READAOP(Y1SPNORB, naos,'Y1SPNORB') !gets Y1SPNORB matrix
      call READAOP(Z1SPNORB, naos,'Z1SPNORB') !gets Z1SPNORB matrix

      

      DO I=1,NUCIND
         ZZEFG(I)= 'ZZEFG'//CHRNOS(I/10)//CHRNOS(MOD(I,10))//CHRNOS(1)
         DO J=1,3
            INDEX=3*(I-1)+J
            GZZEFG(INDEX)= 'GEFG '//CHRNOS(0)//CHRNOS(0)//CHRNOS(INDEX) !
         ENDDO
      
      call READAOP(GZZEFGx, naos,GZZEFG(3*(I-1)+1)) !gets GZZEFG_X(I) matrix
      call READAOP(GZZEFGy, naos,GZZEFG(3*(I-1)+2)) !gets GZZEFG_Y(I) matrix
      call READAOP(GZZEFGz, naos,GZZEFG(3*(I-1)+3)) !gets GZZEFG_Z(I) matrix
      call READAOP(ZZEFGm,naos,ZZEFG(I))            !gets ZZEFG0(I)1 matrix
       


C     Correccion pqp
c     pqp - no chain rule applied
       call AVGABC(cmos,'XDIPVEL ',ZZEFG(I),'XDIPVEL ',B0,PRT,ORBC,Av)
       Aux=Av
       call AVGABC(cmos,'YDIPVEL ',ZZEFG(I),'YDIPVEL ',B0,PRT,ORBC,Av)
       Aux=Aux+Av
       call AVGABC(cmos,'ZDIPVEL ',ZZEFG(I),'ZDIPVEL ',B0,PRT,ORBC,Av)
       Aux=Aux+Av
       EFGC2(I,3)=-CEFGpqp*calfa*calfa*Aux


c     -(nabla qzz)p +qp^2 - chain rule applied - for comparisson
c      call AVGAB(cmos,GZZEFG(3*(I-1)+1),'XDIPVEL ',B0,PRT,ORBC,Av)
c      Aux=-Av
c      call AVGAB(cmos,GZZEFG(3*(I-1)+2),'YDIPVEL ',B0,PRT,ORBC,Av)
c      Aux=Aux-Av
c      call AVGAB(cmos,GZZEFG(3*(I-1)+3),'ZDIPVEL ',B0,PRT,ORBC,Av)
c      Aux=Aux-Av
c      call AVGAB(cmos,ZZEFG(I),'KINENERG',B0,PRT,ORBC,Av)
c      Aux=Aux+2*Av
c      EFGC2(I,3)=CEFGpqp*calfa*calfa*Aux

C     Correccion kin      
       call AVGAB(cmos,'KINENERG',ZZEFG(I),B0,PRT,ORBC,Av)
       EFGC2(I,4)=CEFGkin*calfa*calfa*Av

       
C     Correcciones orden c^(-4)
      


      aux0=matmul(transpose(cmos),cmos)

c     Hago (pq).I x p     
      aux1=matmul(matmul(GZZEFGy,aux0),ZDIPVEL)-
     &     matmul(matmul(GZZEFGz,aux0),YDIPVEL)
      aux2=matmul(matmul(GZZEFGz,aux0),XDIPVEL)-
     &     matmul(matmul(GZZEFGx,aux0),ZDIPVEL)
      aux3=matmul(matmul(GZZEFGx,aux0),YDIPVEL)-
     &     matmul(matmul(GZZEFGy,aux0),XDIPVEL)


      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(1)='1 A  21 ' 
C     Replace time information with symmetry information
      RTNLBL(2)='SQUARE  '


      call WRITEAOPJUAN(aux1,naos,
     &      'E1GLRSx'//CHRNOS(I),RTNLBL) !Write Matrix in AOPROPER
      call WRITEAOPJUAN(aux2,naos,
     &      'E1GLRSy'//CHRNOS(I),RTNLBL) !Write Matrix in AOPROPER
      call WRITEAOPJUAN(aux3,naos,
     &      'E1GLRSz'//CHRNOS(I),RTNLBL) !Write Matrix in AOPROPER

      call WRITEAOPJUAN(aux1+aux2+aux3,naos,
     &      'E1GLRS '//CHRNOS(I),RTNLBL) !Write Matrix in AOPROPER


      ENDDO

      RTNLBL(2)='ANTISYMM'
      call WRITEAOPJUAN(X1SPNORB+Y1SPNORB+Z1SPNORB,naos
     &                  ,'E1G  SO1',RTNLBL) !Write Matrix in AOPROPER

      

      return
      END SUBROUTINE
      
      
      SUBROUTINE AVGAB(cmos,char1,char2,B0,PRT,ORBC,Av)
C     Calculates  the expectation value <A1.A2> as Av, where A1 and A2
C     are defined by char1 and char2. It also calculates the molecular
C     orbital contributions ORBC.
C     J. Aucar 2020
C     ******  ******
#include "inforb.h"
#include "priunit.h"

      double precision cmos(NORBT,NBAST),ORBC(NOCCT)
      CHARACTER*8 B0,char1,char2
      integer PRT, naos
      double precision,intent(out) :: Av
      double precision A1(NBAST,NBAST), A2(NBAST,NBAST)

      naos = NBAST
      nmos=NORBT
      if (B0.Eq.'ORBCON  ') then
       write(lupri,*) ' ******ORBITAL CONTRIBUTIONS DETECTED!*******'
      endif 

      call READAOP(A1,naos,char1)
      call READAOP(A2,naos,char2)
      call cuentaAB(ORBC,A1,A2,cmos,NOCCT*2,naos,NORBT)
!      write(lupri,*) 'CMOS (avgdrv.F):'
!     call OUTPUT(CMOS,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)  ! prints CMOS
!      call OUTPUT(A1,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints A1 for debugging
!      call OUTPUT(A2,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints A2 for debugging
!      call OUTPUT(ORBC,1,1,1,1,NOCCT,NOCCT,1,LUPRI)! prints ORBC for debugging

      if (B0.Eq.'ORBCON  ') then
         write(lupri,*) '< ',char1,' >< ',char2,'>' 
      endif 
      S1 = 0.0

      do i1 = 1,NOCCT
         if (B0.Eq.'ORBCON  ') then
            write(lupri,*) 'Orbital : ',i1, '  Valor', 2*ORBC(i1)
         endif 
         S1 = S1 + ORBC(i1)
      enddo
      Av = 2.0*S1
      if (PRT.GT.1) then
         write(lupri,*) '< ',char1,' >< ',char2,' >  = ', Av
         if (B0.Eq.'ORBCON  ') then
            write(lupri,*) '======================================= '
         endif
      endif

      return
      END
      
      
      SUBROUTINE AVGABC(cmos,char1,char2,char3,B0,PRT,ORBC,Av)
C     Calculates  the expectation value <A1.A2.A3> as Av, where A1,A2, A3
C     are defined by char1,char2,char3. It also calculates the molecular
C     orbital contributions ORBC.
C     J. Aucar 2020
C     ******  ******
#include "inforb.h"
#include "priunit.h"

      double precision cmos(NORBT,NBAST),ORBC(NOCCT)
      CHARACTER*8 B0,char1,char2,char3
      integer PRT, naos,mu,nu
      double precision,intent(out) :: Av
      double precision A1(NBAST,NBAST), A2(NBAST,NBAST),A3(NBAST,NBAST)

      naos = NBAST
      nmos=NORBT
      if (B0.Eq.'ORBCON  ') then
       write(lupri,*) ' ******ORBITAL CONTRIBUTUIONS DETECTED!*******'
      endif 

      call READAOP(A1,naos,char1)
      call READAOP(A2,naos,char2)
      call READAOP(A3,naos,char3)
      call cuentaABC(ORBC,A1,A2,A3,cmos,NOCCT*2,naos,NORBT)
!      call OUTPUT(A1,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints A1 for debugging
!      call OUTPUT(A2,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints A2 for debugging
!      call OUTPUT(A3,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints A3 for debugging
!      call OUTPUT(ORBC,1,1,1,1,NOCCT,NOCCT,1,LUPRI)! prints ORBC for debugging

      if (B0.Eq.'ORBCON  ') then
         write(lupri,*) '< ',char1,' >< ',char2,'>','< ',char3,'>'
      endif 
      S1 = 0.0
      
      do i1 = 1,NOCCT
         if (B0.Eq.'ORBCON  ') then
            write(lupri,*) 'Orbital : ',i1, '  Valor', 2*ORBC(i1)
         endif 
         S1 = S1 + ORBC(i1)
      enddo
      Av = 2.0*S1
      if (PRT.GT.1) then
       write(lupri,*) '< ',char1,' >< ',char2,' >< ',char3,' >  = ',Av
       if (B0.Eq.'ORBCON  ') then
        write(lupri,*) '======================================= '
       endif
      endif

      return
      END
