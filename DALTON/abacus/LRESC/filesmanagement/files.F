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

      subroutine copy_file(oldID,newID,old_file, new_file)
      character(*)   :: old_file ,new_file
      logical        :: file_exists, file_open
      character(100) :: line
      INTEGER        :: oldID,newID
      
#include "priunit.h"
      
c     new file management
      inquire(file = new_file, exist = file_exists,
     &            opened = file_open)
      if (file_exists .AND. file_open) then
         CALL GPCLOSE(newID ,'DELETE')
      else
         if (.not. file_open) then
            CALL GPOPEN(newID,new_file,'NEW',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         end if
      end if
      
c     old file management
      inquire(file = old_file, exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            CALL GPOPEN(oldID,old_file,'OLD',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         end if
         rewind(oldID)
      else
         call quit('file not found')
      end if
      
      
      do while (.true.)
         read(oldID,'(a)', end=31) line
         write(newID,'(a)') line
      end do
 31   continue
      
      if (.not. file_open) CALL GPCLOSE(oldID,'KEEP')
      if (.not. file_open) CALL GPCLOSE(newID,'KEEP')
      
      RETURN
      END
      
      subroutine move_file(oldID,newID,old_file, new_file)
      character(*) :: old_file ,new_file
      logical      :: file_exists, file_open
      character(100) :: line
      INTEGER :: oldID,newID
      
#include "priunit.h"
      
c     new file management
      inquire(file = new_file, exist = file_exists,
     &            opened = file_open)
      if (file_exists .AND. file_open) then
         CALL GPCLOSE(newID ,'DELETE')
      else
         if (.not. file_open) then
            CALL GPOPEN(newID,new_file,'NEW',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         end if
      end if
      
c     old file management
      inquire(file = old_file, exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            CALL GPOPEN(oldID,old_file,'OLD',' ','FORMATTED',
     &                IDUMMY,.FALSE.)
         end if
         rewind(oldID)
      else
         call quit('file not found')
      end if
      
      
      do while (.true.)
         read(oldID,'(a)', end=31) line
         write(newID,'(a)') line
      end do
 31   continue
      
      if (.not. file_open) CALL GPCLOSE(oldID,'REMOVE')
      if (.not. file_open) CALL GPCLOSE(newID,'KEEP')
      
      RETURN
      END
