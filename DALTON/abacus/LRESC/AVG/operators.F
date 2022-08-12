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
        SUBROUTINE READAOP(A, naos, LABEL)
#include "priunit.h"
C
C Linux version, Melo, Azua, Giribet, Alejandrito, Aucar(x2), Danian, Romero, Homero,
C  Bart y Lisa
C
C       A : Coefficient Matrix in Atomic Basis
C           of LABEL operator
C       n : Dim(A)   !<- n is the number of atomic orbitals (naos)
C   LABEL : Searched Operator Label
C
C          03/2020. J. Aucar
C          Dynamic version
C
C -------------------------------------------------
C       LEO PARAMETROS DE AOPROPER
C -------------------------------------------------
      INTEGER*4 nelem, mu, nu, naos, IREC, IERR
      DOUBLE PRECISION,  DIMENSION(naos,naos), INTENT(OUT) :: A
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: integrals
      CHARACTER*8  STARS, LABEL, B(4)
      PARAMETER (STARS = '********')
      PARAMETER (LUINP222 = 222)
C ===================================================
C ------------------------------------------------
C       Reads AOPROPER for building xx vector
C       whose elements are the LABEL operator
C       integrals and put them into A matrix
C ------------------------------------------------

       IF (LABEL.Eq.'IDENTITY') then
           do nu = 1,naos
              do mu = 1,naos
                 A(mu,nu) = 0.0
                 IF (mu.Eq.nu) then
                 A(mu,nu) = 1.0
                 ENDIF
              enddo
           enddo
       ENDIF
       open(LUINP222,file='AOPROPER',status='OLD',FORM='UNFORMATTED')
c      write(*,*) 'Opening AOPROPER file to look for :', LABEL     ! For debugging
       REWIND LUINP222
       IREC = 0
       IERR = 0
  1    READ (LUINP222,END=3,ERR=2) B
       IREC = IREC + 1
       IF (B(1) .NE. STARS) GO TO 1
!       write(*,*) B(1), B(2), B(3), B(4) !FOR DEBUGGING
C      WRITE (*, '(5X,I5,3X,4(2X,A8))')  IREC, B       ! For debugging
       IF (B(4).EQ.LABEL) then
        IF (B(3).EQ.'SYMMETRI') then
         nelem = naos*(naos+1)/2          ! Independent elements
         ALLOCATE(integrals(nelem))
         read (LUINP222,end=3,err=2) (integrals(i),i=1,nelem)
         k = 1
         do nu=1,naos
         do mu=1,nu
            A(mu,nu) = integrals(k)
           A(nu,mu) = A(mu,nu)
           k = k + 1
          enddo
         enddo

      ELSEIF(B(3).EQ.'SQUARE  ') then
            nelem = naos*naos          ! Independent elements
            ALLOCATE(integrals(nelem))
            read (LUINP222,end=3,err=2) (integrals(i),i=1,nelem)
            k = 1
            do nu=1,naos
            do mu=1,naos
                  A(mu,nu) = integrals(k)
              k = k + 1
             enddo
            enddo

        ELSE
         nelem = naos*(naos+1)/2
         ALLOCATE(integrals(nelem))
         read (LUINP222,end=3,err=2) (integrals(i),i=1,nelem)
         k = 1
         do nu = 1,naos
          do mu = 1,nu
            A(mu,nu) = integrals(k)
           A(nu,mu) = -1.0*A(mu,nu)
           k = k + 1
          enddo
         enddo
        ENDIF
        DEALLOCATE(integrals)
       ENDIF
       GO TO 1
C
   2  CONTINUE
      !IREC = IREC + 1
      IERR = IERR + 1
      WRITE (*, '(/A,I5/)') ' ERROR READING RECORD NO.',IERR
      REWIND LUINP222
      DO 102 I = 1,IREC
         READ (LUINP222) J
  102 CONTINUE
      IF (IERR .LE. 2) GO TO 1
  202 CONTINUE
         READ (LUINP222,END=3) J
         IREC = IREC + 1
      GO TO 202
C
   3  CONTINUE
C     WRITE (*,'(/I10,A)') IREC,
C    *   ' records read before EOF on file. AOPROPER'   ! For debugging
      close(LUINP222)
C      write(*,*) 'Closing File AOPROPER '              ! For debugging
C      write(*,*) '*****************************************
      return
      END
C
C ***********************************************************************
C   Subroutines for <A.B> and <A.B.C> calculations
C   Version with built-in Fortran functions - August 2020
C ==================================================================
C       PARAMETROS  
C   ORBC(a) : Orbital Contribution to <A.B> or <A.B.C>
C A1,A2,A3  : Matrices of integrals in atomic basis
C      CMOS : Molecular Coefficients Matrix
C     nelec : Electrons Number
C      naos : Number of orbitals in atomic basis
C      nmos : Number of molecular orbitals
C       a   : Number of occupied molecular orbitals = nelec/2
C -----------------------------------------------------------------
      SUBROUTINE cuentaAB(ORBC,A1,A2,CMOS,nelec,naos,nmos)


       INTEGER*4  naos, nmos, nelec, a
       DOUBLE PRECISION,  DIMENSION(naos,naos), INTENT(IN) :: A1,A2
       DOUBLE PRECISION,  DIMENSION(nmos,naos), INTENT(IN) :: CMOS
       DOUBLE PRECISION,  DIMENSION(nmos,nmos) :: gama, delta
       DOUBLE PRECISION,  DIMENSION(nelec/2), INTENT(OUT) :: ORBC
C -----------------------------------------------------------------
      gama=matmul(CMOS,matmul(A1,transpose(CMOS)))
      delta=matmul(CMOS,matmul(A2,transpose(CMOS)))
      ORBC=0.0

!      call OUTPUT(A1,1,naos,1,naos,naos,naos,1,LUPRI)! prints A1 for debugging
!      call OUTPUT(A2,1,naos,1,naos,naos,naos,1,LUPRI)! prints A2 for debugging
c       call OUTPUT(matmul(transpose(cmos),cmos)
c     &      ,1,naos,1,naos,naos,naos,1,LUPRI)! prints CMOS for debugging
!      call OUTPUT(delta,1,naos,1,naos,naos,naos,1,LUPRI)! prints CMOS for debugging
      do a = 1,nelec/2
       do i = 1,nmos
        ORBC(a) = ORBC(a) + gama(a,i)*delta(i,a)
       enddo
      enddo
      return
      end subroutine
C
C *****************************************************************
       SUBROUTINE cuentaABC(ORBC,A1,A2,A3,CMOS,nelec,naos,nmos)
       INTEGER*4  naos, nmos, nelec, a
       DOUBLE PRECISION,  DIMENSION(naos,naos), INTENT(IN) :: A1,A2,A3
       DOUBLE PRECISION,  DIMENSION(nmos,naos), INTENT(IN) :: CMOS
       DOUBLE PRECISION,  DIMENSION(nmos,nmos) :: alpha,beta,gama
       DOUBLE PRECISION,  DIMENSION(nelec/2), INTENT(OUT) :: ORBC
C -----------------------------------------------------------------
      alpha=matmul(matmul(CMOS,A1),transpose(CMOS))
      beta=matmul(matmul(CMOS,A2),transpose(CMOS))
      gama=matmul(matmul(CMOS,A3),transpose(CMOS))
      ORBC=0.0

      do a = 1,nelec/2
       do i = 1,nmos
        do j=1,nmos
         ORBC(a) = ORBC(a) + alpha(a,i)*beta(i,j)*gama(j,a)
        enddo
       enddo
      enddo

      return
      end

C -----------------------------------------------------------------
      SUBROUTINE WRITEAOPJUAN(A, naos, LABEL,RTNLBL)
#include "priunit.h"
#include "iratdef.h"
#include "inforb.h"
C
C J. Aucar - 2021
C
C       A : Coefficient Matrix in Atomic Basis
C           of LABEL operator
C       naos : Dim(A)
C   LABEL : Searched Operator Label
C
C -------------------------------------------------
C       READ AOPROPER PARAMETERS
C -------------------------------------------------
      INTEGER*4 nelem, mu, nu, naos, IREC, IERR
      DOUBLE PRECISION,  DIMENSION(naos,naos), INTENT(IN) :: A
      CHARACTER*8  STARS, LABEL, RTNLBL(2)
      PARAMETER (LUINP222 = 222)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: integrals
      PARAMETER (STARS = '********')

C ===================================================
C ------------------------------------------------
C       Reads AOPROPER for building xx vector
C       whose elements are the LABEL operator
C       integrals and put them into A matrix

C ------------------------------------------------
      OPEN (LUINP222,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER',position="append")


      BACKSPACE LUINP222
      
c      write(*,*) 'Opening AOPROPER file to write ', LABEL,' on it'     ! For debugging
      
      IREC = 0
      IERR = 0


      CALL NEWLB2(LABEL,RTNLBL,LUINP222,LUPRI)
c      write(*,*) STARS,RTNLBL,LABEL ! for debugging
      
      IF (RTNLBL(2).EQ.'SYMMETRI') then
            nelem = naos*(naos+1)/2      !Independent Elements
            ALLOCATE(integrals(nelem))
            k = 1
            do nu=1,naos
            do mu=1,nu
                  integrals(k) = A(mu,nu)
                  k = k + 1
            enddo
            enddo
            
      ELSEIF(RTNLBL(2).EQ.'ANTISYMM') THEN
            nelem = naos*(naos+1)/2        !Independent Elements
            ALLOCATE(integrals(nelem))
            k = 1
            do nu = 1,naos
            do mu = 1,nu
                  integrals(k) = A(mu,nu)
                  k = k + 1
            enddo
            enddo
      ELSE
            nelem = naos*naos       !Independent Elements
            ALLOCATE(integrals(nelem))
            k = 1
            do nu = 1,naos
            do mu = 1,naos
                  integrals(k) = A(mu,nu)
                  k = k + 1
            enddo
            enddo
      ENDIF


      LEN = MAX(4,nelem) !Same as in gp/gphjj.F

C     Write integrals in AOPROPER. IRAT=2 for 64bits integers and 1
C     otherwise
      CALL WRITI(LUINP222,IRAT*LEN,integrals)
      
C     We add an extra label to signify EOF
      CALL NEWLB2('EOFLABEL',RTNLBL,LUINP222,LUPRI)
      close(LUINP222)


c      write(*,*) 'Closing File AOPROPER '              ! For debugging
      DEALLOCATE(integrals)



      
      RETURN
      END
