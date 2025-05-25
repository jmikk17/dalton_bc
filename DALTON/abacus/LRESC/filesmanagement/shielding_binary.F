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


      
      lresc(:,:,1)=SGPNR                         ! nr-para
      lresc(:,:,2)=SGDNR                         ! nr-dia
      lresc(:,:,3)=SGPNR+SGDNR                   ! nr-total
      lresc(:,:,4)=SGP1S+SGP1T+SGP3S+SGP3T       ! lr-para
      lresc(:,:,5)=SGD0S+SGD1S                   ! lr-dia
      lresc(:,:,6)=SLRESC                        ! lr-total
      lresc(:,:,7)=SGPNR+SGP1S+SGP1T+SGP3S+SGP3T ! rel-para
      lresc(:,:,8)=SGDNR+SGD0S+SGD1S             ! rel-dia
      lresc(:,:,9)=SGPNR+SGDNR+SLRESC            ! rel-total
      lresc(:,:,10)=LROZK                        ! L-kin
      lresc(:,:,11)=LRPSOK                       ! Pso-kin
      lresc(:,:,12)=LRFCZK+LRFCDE                ! FC
      lresc(:,:,13)=LRSDK+LRSDD                  ! SD
      lresc(:,:,14)=LRFCBS                       ! FC-Bso
      lresc(:,:,15)=LRSDBS                       ! SD-Bso
      lresc(:,:,16)=QRPSOOZMV                    ! PSO-L-Mv
      lresc(:,:,17)=QRPSOOZDW                    ! PSO-L-Dw
      lresc(:,:,18)=QROZFCSO                     ! FC-L-So
      lresc(:,:,19)=QROZSDSO                     ! SD-L-So
      lresc(:,:,20)=LRFCAV                       ! FC
      lresc(:,:,21)=LRANGP                       ! AngPso
      lresc(:,:,22)=LRDIAK                       ! Diakin
      lresc(:,:,23)=LRSDAV                       ! SD
      lresc(:,:,24)=LRDIAM                       ! DiaMv
      lresc(:,:,25)=LRDIAD                       ! DiaDw
      end


      subroutine binary_lresc(array)
      USE lrescinf_module
#include "mxcent.h"
#include "nuclei.h"
      real(8), intent(in) :: array(NLRATOMS,3,3,25)
      character(len=*), parameter :: archivo = "lresc-array.npy"
        
      ! Open binary file in stream mode
      open(newunit=iu, file=archivo, action="write", 
     &     access="stream", form="unformatted")


      ! Write info about nuclei names
      ! write(iu) NAMN(LRATOM)


      ! Write info about dimensions
      write(iu) NLRATOMS,3,3,25

      ! Write data
      write(iu) array

      close(iu)

      end subroutine
