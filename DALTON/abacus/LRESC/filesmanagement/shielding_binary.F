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


      subroutine create_array_oneatom(lresc)
      USE lrescinf_module
      real(8), intent(out) :: lresc(3,3,25)
      lresc=0.d0


C     Acronyms according to https://doi.org/10.1063/5.0264596
C      
      lresc(:,:,1)=SGPNR                         ! nr-para
      lresc(:,:,2)=SGDNR                         ! nr-dia
      lresc(:,:,3)=SGPNR+SGDNR                   ! nr-total
      lresc(:,:,4)=SGP1S+SGP1T+SGP3S+SGP3T       ! lr-para
      lresc(:,:,5)=SGD0S+SGD1S                   ! lr-dia
      lresc(:,:,6)=SLRESC                        ! lr-total
      lresc(:,:,7)=SGPNR+SGP1S+SGP1T+SGP3S+SGP3T ! rel-para
      lresc(:,:,8)=SGDNR+SGD0S+SGD1S             ! rel-dia
      lresc(:,:,9)=SGPNR+SGDNR+SLRESC            ! rel-total
      lresc(:,:,10)=LROZK                        ! OZK
      lresc(:,:,11)=LRPSOK                       ! PSOK
      lresc(:,:,12)=LRFCZK+LRFCDE                ! SZK(Fc)
      lresc(:,:,13)=LRSDK+LRSDD                  ! SZK(SD)
      lresc(:,:,14)=LRFCBS                       ! Bso(FC)
      lresc(:,:,15)=LRSDBS                       ! BSO(SD)
      lresc(:,:,16)=QRPSOOZMV                    ! Mv/p
      lresc(:,:,17)=QRPSOOZDW                    ! Dw/p
      lresc(:,:,18)=QROZFCSO                     ! SO(FC)
      lresc(:,:,19)=QROZSDSO                     ! SO(SD)
      IF (PYINTS) lresc(:,:,20)=LRFCAVpy         ! FC/d
      IF(.NOT.PYTINTS) lresc(:,:,20)=LRFCAV      ! FC/d
      lresc(:,:,21)=LRANGPavg                    ! PSO-OZ
      lresc(:,:,22)=LRDIAK                       ! DIAK
      lresc(:,:,23)=LRSDAV                       ! SD/d
      lresc(:,:,24)=LRDIAM                       ! Mv/d
      lresc(:,:,25)=LRDIAD                       ! Dw/d
      end


      subroutine binary_lresc(array)
      USE lrescinf_module
#include "mxcent.h"
#include "nuclei.h"
      real(8), intent(in) :: array(NLRATOMS,3,3,25)
      character(len=*), parameter :: archivo = "lresc-shield-array.npy"
      character(len=10), dimension(25) :: labels
      data labels / 'nr_para  ', 'nr_dia   ', 'nr_total ', 'lr_para  ',
     &         'lr_dia   ', 'lr_total ', 'rel_para ', 'rel_dia  ',
     &         'rel_total', 'OZK      ', 'PSOK     ', 'SZK-FC   ',
     &         'SZK-SD   ', 'BSO-FC   ', 'BSO-SD   ', 'Mv/p     ',
     &         'Dw/p     ', 'SO-FC    ', 'SO-SD    ', 'FC/d     ',
     &         'PSO-OZ   ', 'DIAK     ', 'SD/d     ', 'Mv/d     ',
     &         'Dw/d     ' /

        
      ! Open binary file in stream mode
      open(newunit=iu, file=archivo, action="write", 
     &     access="stream", form="unformatted")


      ! Write info about nuclei names
      ! write(iu) NAMN(LRATOM)


      ! Write info about dimensions
      write(iu) NLRATOMS,3,3,25

      ! Write data
      write(iu) array

      ! Write labels
      write(iu) labels


      close(iu)

      end subroutine
