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

C  /* Deck searchword */
      subroutine SearchWord(file_unit, word_to_search, found)
C     Search certain word in a specific file.
C     The file must be open before calling to this subroutine.
C     J. Aucar 2023
C
#include "priunit.h"
      integer, intent(in) :: file_unit
      character(len=*), intent(in) :: word_to_search
      logical, intent(out) :: found
      character(len=100) :: line
      integer :: i

!     Initialize the control variable
      found = .false.

!     Read line by line from the file
      do
         read(file_unit, '(A)', iostat=i) line
         if (i /= 0) exit ! End of file
!        Search for the word in the line
         if (index(line, word_to_search) /= 0) then
            found = .true.
            exit
         end if
      end do

      end subroutine SearchWord


C  /* Deck find_warning_lr */
      subroutine find_warning_lr(filename, operator, has_warning)
C     Search for warnings in LR related to a specific operator, in filename.
C     The file must be open before calling to this subroutine.
C     J. Aucar 2025
C
      implicit none
      character(len=*), intent(in)  :: filename
      character(len=*), intent(in)  :: operator
      logical, intent(out)          :: has_warning
      integer :: pos
      character(len=512) :: extracted


      character(len=512) :: line
      logical :: in_block
      integer :: ios, unit

      has_warning = .false.
      in_block = .false.

      open(newunit=unit, file=filename, status='old',
     &     action='read', iostat=ios)
      if (ios /= 0) then
          print *, "Error while trying to open file ", filename
          return
      end if

      do
          read(unit,'(A)',iostat=ios) line
          if (ios /= 0) exit

          ! Detect line
          pos = index(line, "RSPLR -- operator label :")
          if (pos > 0) then
          ! Extract text after ':'
              extracted = adjustl(line(pos+len
     &                    ("RSPLR -- operator label :"):))
              if (index(extracted, trim(operator)) > 0) then
              in_block = .true.
              end if
          end if

          if (in_block) then
          ! Find WARNING
          if (index(line, "WARNING") > 0) then
              ! Deal with PSO vs PSOKJ as an exception
              if (operator.eq."PSO") then
              if (index(line, "PSOKJ") > 0) then
                  has_warning = .false.
              end if
              else
              has_warning = .true.
              end if
          end if

          ! Detect end of block
          if (index(line,
     &    "Convergence of RSP solution vectors") > 0) then
              in_block = .false.
          end if
          end if

      end do

      close(unit)

      end subroutine find_warning_lr


C  /* Deck find_warning_qr */
      subroutine find_warning_qr(filename, operator, has_warning)
C     Search for warnings in QR related to a specific operator, in filename.
C     The file must be open before calling to this subroutine.
C     J. Aucar 2025
C
      implicit none
      character(len=*), intent(in)  :: filename
      character(len=*), intent(in)  :: operator
      logical, intent(out)          :: has_warning
      integer :: pos, endp, start
      character(len=512) :: extracted, line
      integer :: ios, unit

      has_warning = .false.

      open(newunit=unit, file=filename, status='old',
     &     action='read', iostat=ios)
      if (ios /= 0) then
        print *, "Error while trying to open file ", filename
        return
      end if

      do
        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) exit   ! fin del archivo

        ! Find WARNING
        pos = index(line, "@ WARNING---- Response label ")
        if (pos > 0) then
          start = pos + len("@ WARNING---- Response label ")

          ! Deal with spaces
          do while (start <= len_trim(line) .and.
     &              line(start:start) == ' ')
            start = start + 1
          end do

          ! Extract word until space or end of line
          endp = start
          do while (endp <= len_trim(line) .and. line(endp:endp) /= ' ')
            endp = endp + 1
          end do
          extracted = line(start:endp-1)

          ! Comparisson
          if (trim(extracted) == trim(operator)) then
            has_warning = .true.
            exit
          end if
        end if

      end do

      close(unit)

      end subroutine find_warning_qr



C  /* Deck check_lr_s */
      subroutine check_lr_s
C     Check which operators are associated to warnings in
C     singlet linear responses. Print which LRESC corrections
C     are associated to those operators.
C     J. Aucar 2025
C
      USE lrescinf_module
#include "implicit.h"
#include "priunit.h"


      logical :: warn_darwin, warn_massvelo
      logical :: warn_psokj, warn_nsco
      logical :: warn_pso, warn_ozke, warn_angmom


      call find_warning_lr("LRESC-LR-S", "DARWIN", warn_darwin)
      call find_warning_lr("LRESC-LR-S", "NSCO", warn_nsco)
      call find_warning_lr("LRESC-LR-S", "MASSVELO", warn_massvelo)

      call find_warning_lr("LRESC-LR-S", "PSOKJ", warn_psokj)
      call find_warning_lr("LRESC-LR-S", "PSO", warn_pso)
      call find_warning_lr("LRESC-LR-S", "OZKE", warn_ozke)
      call find_warning_lr("LRESC-LR-S", "ANGMOM", warn_angmom)



      IF(SIGMAD1S) THEN
      if ((warn_darwin).or.(warn_nsco)) then
        write(LUPRI,'(7X,A)')
     &  "In Dw/d correction. Operator(s) involved:"
        if (warn_darwin) write(LUPRI,'(7X,A)') "DARWIN"
        if (warn_nsco) write(LUPRI,'(7X,A)') "NSCO"
      end if

      if ((warn_massvelo).or.(warn_nsco)) then
        write(LUPRI,'(7X,A)')
     &  "In Mv/d correction. Operator(s) involved:"
        if (warn_massvelo) write(LUPRI,'(7X,A)') "MASSVELO"
        if (warn_nsco) write(LUPRI,'(7X,A)') "NSCO"
      end if
      END IF

      IF(SIGMAP1S) THEN
      if ((warn_pso).or.(warn_ozke)) then
        write(LUPRI,'(7X,A)')
     &  "In OZK correction. Operator(s) involved:"
        if (warn_pso) write(LUPRI,'(7X,A)') "PSO"
        if (warn_ozke) write(LUPRI,'(7X,A)') "OZKE"
      end if

      if ((warn_psokj).or.(warn_angmom)) then
        write(LUPRI,'(7X,A)')
     &  "In PSOK correction. Operator(s) involved:"
        if (warn_psokj) write(LUPRI,'(7X,A)') "PSOKJ"
        if (warn_angmom) write(LUPRI,'(7X,A)') "ANGMOM"
      end if
      END IF



      end subroutine check_lr_s


C  /* Deck check_lr_t */
      subroutine check_lr_t
C     Check which operators are associated to warnings in
C     triplet linear responses. Print which LRESC corrections
C     are associated to those operators.
C     J. Aucar 2025
C
      USE lrescinf_module
#include "implicit.h"
#include "priunit.h"

      logical :: warn_fc, warn_sd, warn_kinenerg, warn_ddd
      logical :: warn_somf



      call find_warning_lr("LRESC-LR-T", "FC", warn_fc)
      call find_warning_lr("LRESC-LR-T", "SD", warn_sd)
      call find_warning_lr("LRESC-LR-T", "KINENERG", warn_kinenerg)
      call find_warning_lr("LRESC-LR-T", "dd/d", warn_ddd)
      call find_warning_lr("LRESC-LR-T", "SOMF", warn_somf)


      IF(SIGMAP1T) THEN
      if ((warn_fc).or.(warn_kinenerg).or.(warn_ddd)) then
        write(LUPRI,'(7X,A)')
     &  "In SZK(Fc) correction. Operator(s) involved:"
        if (warn_fc) write(LUPRI,'(7X,A)') "FC"
        if (warn_kinenerg) write(LUPRI,'(7X,A)') "KINENERG"
        if (warn_ddd) write(LUPRI,'(7X,A)') "dd/d"
      end if

      if ((warn_sd).or.(warn_kinenerg).or.(warn_ddd)) then
        write(LUPRI,'(7X,A)')
     &  "In SZK(Sd) correction. Operator(s) involved:"
        if (warn_sd) write(LUPRI,'(7X,A)') "SD"
        if (warn_kinenerg) write(LUPRI,'(7X,A)') "KINENERG"
        if (warn_ddd) write(LUPRI,'(7X,A)') "dd/d"
      end if

      if ((warn_fc).or.(warn_somf)) then
        write(LUPRI,'(7X,A)')
     &  "In BSO(Fc) correction. Operator(s) involved:"
        if (warn_fc) write(LUPRI,'(7X,A)') "FC"
        if (warn_somf) write(LUPRI,'(7X,A)') "SOMF"
      end if

      if ((warn_sd).or.(warn_somf)) then
        write(LUPRI,'(7X,A)')
     &  "In BSO(Sd) correction. Operator(s) involved:"
        if (warn_sd) write(LUPRI,'(7X,A)') "SD"
        if (warn_somf) write(LUPRI,'(7X,A)') "SOMF"
      end if
      END IF



      end subroutine check_lr_t


C  /* Deck check_qr */
      subroutine check_qr
C     Check which operators are associated to warnings in
C     quadratic responses. Print which LRESC corrections
C     are associated to those operators.
C     J. Aucar 2025
      USE lrescinf_module
#include "implicit.h"
#include "priunit.h"

      logical :: warn_xso, warn_yso, warn_zso
      logical :: warn_darwin_qr, warn_pso_qr, warn_oz_qr, warn_mv_qr
      logical :: warn_so_qr, warn_fc_qr, warn_sd_qr

      call find_warning_qr("LRESC-QR", "DARWIN", warn_darwin_qr)
      call find_warning_qr("LRESC-QR", "PSO", warn_pso_qr)
      call find_warning_qr("LRESC-QR", "OZ", warn_oz_qr)
      call find_warning_qr("LRESC-QR", "MASSVELO", warn_mv_qr)

      warn_so_qr = .false.
      call find_warning_qr("LRESC-QR", "X1SPNORB", warn_xso)
      call find_warning_qr("LRESC-QR", "Y1SPNORB", warn_yso)
      call find_warning_qr("LRESC-QR", "Z1SPNORB", warn_zso)
      if ((warn_xso).or.(warn_yso).or.(warn_zso)) then
        warn_so_qr = .true.
      end if

      call find_warning_qr("LRESC-QR", "FC", warn_fc_qr)
      call find_warning_qr("LRESC-QR", "SD", warn_sd_qr)

      IF(SIGMAP2S) THEN
      if ((warn_pso_qr).or.(warn_oz_qr).or.(warn_mv_qr)) then
        write(LUPRI,'(7X,A)')
     &  "In Mv/p correction. Operator(s) involved:"
        if (warn_pso_qr) write(LUPRI,'(7X,A)') "PSO"
        if (warn_oz_qr) write(LUPRI,'(7X,A)') "OZ"
        if (warn_mv_qr) write(LUPRI,'(7X,A)') "MASSVELO"
      end if

      if ((warn_pso_qr).or.(warn_oz_qr).or.(warn_darwin_qr)) then
        write(LUPRI,'(7X,A)')
     &  "In Dw/p correction. Operator(s) involved:"
        if (warn_pso_qr) write(LUPRI,'(7X,A)') "PSO"
        if (warn_oz_qr) write(LUPRI,'(7X,A)') "OZ"
        if (warn_darwin_qr) write(LUPRI,'(7X,A)') "DARWIN"
      end if
      END IF

      IF(SIGMAP2T) THEN
      if ((warn_so_qr).or.(warn_oz_qr).or.(warn_fc_qr)) then
        write(LUPRI,'(7X,A)')
     &  "In SO(Fc) correction. Operator(s) involved:"
        if (warn_so_qr) write(LUPRI,'(7X,A)') "SO"
        if (warn_oz_qr) write(LUPRI,'(7X,A)') "OZ"
        if (warn_fc_qr) write(LUPRI,'(7X,A)') "FC"
      end if

      if ((warn_so_qr).or.(warn_oz_qr).or.(warn_sd_qr)) then
        write(LUPRI,'(7X,A)')
     &  "In SO(Sd) correction. Operator(s) involved:"
        if (warn_so_qr) write(LUPRI,'(7X,A)') "SO"
        if (warn_oz_qr) write(LUPRI,'(7X,A)') "OZ"
        if (warn_sd_qr) write(LUPRI,'(7X,A)') "SD"
      end if
      END IF

      end subroutine check_qr
