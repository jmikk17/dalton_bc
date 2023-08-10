
!  Dalton, a molecular electronic structure program
!  Copyright (C) by the authors of Dalton.
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

      SUBROUTINE import_matrix(noperators,filename)
C
C  Created by: J. Aucar 2022
C
C Check include files              
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inforb.h"
#include "nuclei.h"
#include "inftap.h"
#include "lrescinf.h"
#include "chrnos.h"

      CHARACTER*8 RTNLBL(2)
      integer :: noperators
      double precision,allocatable :: matrices(:,:,:)
      character(len=8), dimension(noperators) :: newlabels
      character (len=8) :: label
      character (len=3) :: filename
      double precision,allocatable :: matrix(:,:)
      integer :: operadores

C ------------------------------------------------------
C       DALTON.CM variables
C ------------------------------------------------------
      naos = NBAST ! NBAST: Total number of (atomic) basis functions
      nmos = NORBT !NORBT: Total number of molecular orbitals
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(1)='1 A  21 '

c     Import new matrices from PySCF
      allocate (matrices(naos,naos,noperators))
      call import_numpytensor(matrices,newlabels,naos,noperators
     &                        ,operadores,filename)
      do ii = 1, operadores
         allocate (matrix(naos,naos))
         matrix = matrices(:,:,ii)
         label = newlabels(ii)
         RTNLBL(2)='SQUARE  '
         call WRITEAOPJUAN(matrix,naos,label,RTNLBL) !Write Matrix in AOPROPER
c        write(LUPRI,*) label
c        call OUTPUT(matrix,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)! prints matrix for debugging
         deallocate (matrix)
      end do
      deallocate(matrices)


      return
      END SUBROUTINE



      subroutine import_numpytensor(matrices,labels,naos,noperators
     &                              ,operadores,filename)

      use iso_fortran_env
      double precision,intent(out) :: matrices(naos,naos,noperators)
      double precision :: matrix(naos,naos)
      character(len=8), dimension(noperators),intent(out) :: labels
      character(len=3), intent(in) :: filename
      integer, intent(out) :: operadores
      integer, intent(in) :: naos,noperators
      integer :: naospyscf
      character (len=8) :: label
      open(newunit=iu,file=filename//".bin"
     &   ,access="stream",status="old",action="read")
      REWIND iu
      read(iu) operadores
      read(iu) naospyscf
      if (naos.ne.naospyscf) then
      write(*,*) "naos(dal,pyscf):",naos,naospyscf
      write(*,*) "Operators to import from PySCF:",operadores
         CALL QUIT('naos not equal to naospyscf')
      end if
      do i=1, operadores
         read(iu) label
c        write(*,*) "reading:",label
         read(iu) matrix
         matrices(:,:,i)=matrix
         labels(i)=label
      end do
      close(iu)


      end subroutine

