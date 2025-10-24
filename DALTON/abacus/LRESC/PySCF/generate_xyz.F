
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

      SUBROUTINE generate_xyz()
C
C  Created by: J. Aucar 2023
C
C Check include files              
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "nuclei.h"

       PARAMETER (TOANGS=0.52917724924D+00)

C ------------------------------------------------------
      write(LUPRI,*) NATOMS
      write(LUPRI,*) "XYZ file generated from DALTON"
      do ii= 1, NATOMS
         write(LUPRI,*) NAMN(ii), TOANGS*CORD(1:3,ii)
      end do

      open(12, file="MOLECULE.xyz", status="new", action="write")
      write(12,'(I0)') NATOMS
      write(12,'(A)') "XYZ file generated from DALTON"
      do ii= 1, NATOMS
        write(12,'(I0,3F13.9)') int(CHARGE(ii)),TOANGS*CORD(1:3,ii)
      end do
      close(12)



      return
      END SUBROUTINE
