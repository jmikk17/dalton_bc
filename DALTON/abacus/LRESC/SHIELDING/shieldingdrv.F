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
C
C  /* Deck lrscinp_shi */      
C This module read nesserary input for LRESC calculations.
C jim-gesc : May-2012
C
C      
      SUBROUTINE LRSCINP_SHI(WORD)
C
C Purpose:
C     Options initialization of LRESC.
C     This routine is called by ABAINP : abacus/abadrv.F
C Author:
C     J. I. Melo
C  Edited:
C     J. J. Aucar 2021 (LRESC Structure)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "lrescinf.h"
#include "nuclei.h"
#include "abainf.h"
      PARAMETER (NTABLE = 12, D0 = 0.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
C#include "dorps.h"
#include "spnout.h"
c#include "lrescinf.h"
C                      1         2         3	    4	      5
      DATA TABLE /'.GAUGEO','.PRINT ','.SELECT','.PRTALL','.xXXXXx',
C                      6         7         8	    9	      10
     &            '.PARA1S','.PARA1T','.PARA3S','.PARA3T','.DIAM0S',
C                      11       12
     &            '.DIAM1S','.ANISO '/
C
c      write(lupri,*) ' LRSCOPTS :  called with : ' , WORD
      NEWDEF = (WORD .EQ. '*LROPTS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
c           write(lupri,*) ' @LRSCINP : reading : ' , WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12), I
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//
     &             ' input keywords',LUPRI)
                  GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &               '" not recognized in LRESC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//
     &         ' input keywords',LUPRI)
               CALL QUIT('Input keywords in LRESC Input, LRSCINP')
C   READ GAUGEO : DEFINIR ALGUNA VAR PARA GAUGEO
 1             CONTINUE
cv                write(lupri,*)' ..aca leo gaugeo ', WORD
                  READ (LUCMD,*) (LRGAUG(IS), IS = 1, 3)
                  DO ICENT = 1, NUCIND
c                 NAME =  NAMEX(3*ICENT)(1:4)
                     WRITE (LUPRI,'(2X,A,3X," : ",3(A1,2X,A,F15.10))')
     &                  NAMEX(3*ICENT)(1:4), '1' , 'x' , CORD(1,ICENT),
     &                  '2' , 'y' , CORD(2,ICENT),
     &                  '3' , 'z' , CORD(3,ICENT)
                  ENDDO
                  GAUCHANG =.TRUE.
                  ICHANG = ICHANG + 1
                  GOTO 100
C   SET PRINT LEVEL
 2             CONTINUE
                  READ (LUCMD,*) JIMPRT
                  ICHANG = ICHANG + 1
                  GOTO 100
C   SET WHICH NUC TO DO LRESC
 3             CONTINUE

cx                write(lupri,*)'@LRINP  antes de leer'
                  READ (LUCMD,*) LRATOM
                  ICHANG = ICHANG + 1
cx                write(lupri,*)'LRINP  despues de leer', LRATOM
cx                write(lupri,*)'Your selection is atom #:',LRATOM,
cx     $          'named : ', NAMN(LRATOM), 'on Molecule.mol file'
                  IF (LRATOM.GT.NATOMS)THEN
                     write(lupri,*) 'Input Error on LRSCINP'
                     write(lupri,*) 'Your selection exceeds total ',
     &                    'number of atoms on input file'
                     CALL QUIT('LRATOM greater than NATOMS, LRSCINP')
                  END IF

                  GOTO 100
 4             CONTINUE
                  PRTALL1 = .TRUE.
                  GOTO 100
 5             CONTINUE
c                  do algo 
                  GOTO 100
 6             CONTINUE
                  SIGMAP1S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 7             CONTINUE
                  SIGMAP1T = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 8             CONTINUE
                  SIGMAP3S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 9             CONTINUE
                  SIGMAP3T = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 10            CONTINUE
                  SIGMAD0S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 11            CONTINUE
                  SIGMAD1S=.TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 12            CONTINUE
                  LRANISO = .TRUE.
                  GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &                             '" not recognized in LRESC'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in LRESC')
            END IF
         END IF
 300  CONTINUE
c jim-dbg ICHANG > 0 cambios en default variables
c jim-dbg ICHANG < 0 no cambios


      RETURN
      END
C -----------------------------------------------------------------------
C
C  /* Deck lrscini_shi */
      SUBROUTINE LRSCINI_SHI
C
C Purpose:
C     Initialize the variables common LRESCINI: include/lrescinf.h.
C     This routine is call by ABAINIALL : abacus/abadrv.F
C Author:
C     J. I. Melo
C Edited:
C     J. J. Aucar. 2021
C
#include "implicit.h"
#include "lrescinf.h"
#include "mxcent.h"
#include "cbiqr.h"
#include "abainf.h"
C
      PARAMETER(D0 = 0.0D+0)

      JIMPRT    = 0

      IF (.not.LRATOM_changed) THEN
         LRATOM    = 1
      END IF

      LRGAUG(1) = D0
      LRGAUG(2) = D0
      LRGAUG(3) = D0

      GAUCHANG = .FALSE.

      SIGMAP1S = .FALSE.
      SIGMAP1T = .FALSE.
      SIGMAD1S = .FALSE.
      SIGMAD0S = .FALSE.
      SIGMAP3S = .FALSE.
      SIGMAP3T = .FALSE.

      LRESCALL = .TRUE.
C
C     Initialize /CUADRA/
C
      IPRINT = IPRDEF
      IPRQR  = IPRINT
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      THRESH = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
      NABAPP = 0
      LBFREQ = 1
      LCFREQ = 1
      CALL DZERO (QBFREQ,NFMAX)
      CALL DZERO (QCFREQ,NFMAX)
C
C  Init Result Matrices
C jim 4.0 : anddy aca ver que matrices encesita gesc.

!LRESC corrections to dia and para terms of shielding
      LRFCAV = D0
      LRDIAK = D0
      LRANGP = D0
      LRSDAV = D0

      LRDIAM = D0
      LRDIAD = D0

      LROZK  = D0
      LRPSOK = D0
      LRPSKI = D0  ! this is just for debugging
      LRFCZK = D0
      LRFCDE = D0
      LRSDK = D0
c      LRSDK2 = D0
      LRSDD = D0
c      LRSDD2 = D0
      LRFCBS = D0
      LRSDBS = D0

      RETURN
      END
c    ---------------------------------------
C  /* Deck lrscdrv_shi */
      SUBROUTINE LRSCDRV_SHI(WORK,LWORK)
C
C Purpose:
C     Main routine of LRESC module.
C     This routine is call by ABACTL : abacus/abadrv.F
C Author:
C     J. I. Melo
C Edited:
C     J. J. Aucar 2021
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
#include "nuclei.h"
#include "chrnos.h"              
#include "gnrinf.h"
C
C
      LOGICAL PING, FOUND
      DIMENSION WORK(LWORK)
      logical :: file_exists, file_open
      INTEGER, allocatable :: LRATOMS(:)
      PING=.FALSE.
      CALL QENTER('LRSCDRV_SHI')

C ----------------------------
      IF (LRESCALL) THEN
         SIGMAD0S = .TRUE.
         SIGMAD1S = .TRUE.
         SIGMAP1S = .TRUE.
         SIGMAP1T = .TRUE.
         SIGMAP3S = .TRUE.
         SIGMAP3T = .TRUE.
      ENDIF

C     Print information about LRESC metodology
      CALL HEADER_LRSHIELD
C
C     Diamagnetic corrections to Shielding:
C     Zeroth Order, Singlet of course
C
c     Prepared for selecting a nuclei sublist
      IF (LRATOM_changed) THEN
         allocate(LRATOMS(1))
         LRATOMS(1) = LRATOM
      ELSE
         allocate(LRATOMS(NUCIND))
         DO I = 1, NUCIND
            LRATOMS(I)=I
         END DO
      END IF
C --------------------------------------------------
C Some files have to be closed and/or deleted
c when the RESPONSE module is called before LRESC
      inquire(file = 'RSPVEC', exist = file_exists,
     &            opened = file_open)
      IF (file_exists .AND. file_open) 
     &      CALL GPCLOSE(LURSP ,'DELETE') !RSPVEC file


      inquire(file = 'RSP_BVEC', exist = file_exists,
     &            opened = file_open)
      IF (file_exists .AND. file_open) 
     &      CALL GPCLOSE(LURSP3 ,'DELETE')


      inquire(file = 'RSP_DIAG', exist = file_exists,
     &            opened = file_open)
      IF (file_exists .AND. file_open) 
     &      CALL GPCLOSE(LURSP4 ,'DELETE')


      inquire(file = 'RSP_SVEC', exist = file_exists,
     &            opened = file_open)
      IF (file_exists .AND. file_open) 
     &      CALL GPCLOSE(LURSP5 ,'DELETE')
C --------------------------------------------------
      call copy_file(LUESC,LUESC2,'DALTON.INP', 'DALTON_BACKUP.INP')
C --------------------------------------------------
c     Construct LR singlet inputs and run responses (all at once)
      IF (SIGMAP1S.OR.SIGMAD1S) THEN
         call LRINPUT_ALL(LRATOMS,SIZE(LRATOMS),.FALSE.)
c        Run Linear Responses
         LUPRI = 6232
         OPEN(UNIT=LUPRI, FILE='LRESC-LR-S',status='unknown',
     &        position='append')
         call RSPDRV(WORK,LWORK)
         CLOSE(LUPRI)
         LUPRI  = 2
      END IF
C --------------------------------------------------
C Some files have to be deleted
c when calling for TRIPLET responses after singlet responses
      inquire(file = 'RSPVEC', exist = file_exists)
      IF (file_exists) THEN
         CALL GPOPEN(LURSP,'RSPVEC','OLD',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         CALL GPCLOSE(LURSP ,'DELETE')
      END IF

      LURSP1 = -1
      inquire(file = 'RSPRST.E2C', exist = file_exists)
      IF (file_exists) THEN
         CALL GPOPEN(LURSP1,'RSPRST.E2C','OLD',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         CALL GPCLOSE(LURSP1 ,'DELETE')
      END IF
C --------------------------------------------------
c     Construct LR triplet inputs and run responses (all at once)
      IF (SIGMAP1T) THEN
         call LRINPUT_ALL(LRATOMS,SIZE(LRATOMS),.TRUE.)
c        Run Linear Responses
         LUPRI = 6232
         OPEN(UNIT=LUPRI, FILE='LRESC-LR-T',status='unknown',
     &        position='append')
         call RSPDRV(WORK,LWORK)
         CLOSE(LUPRI)
         LUPRI  = 2
      END IF
C --------------------------------------------------
C NUEVOOOOO 1/r3
c      KFREE  = 1
c      LFREE  = LWORK

c      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK ,KFREE ,LFREE)
c      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO))

c      IF (.NOT.FOUND)
c     & CALL QUIT('***Error AVELRSC: CMO is not in SIRIFC***')


c         WRITE(lupri,*)' CMOS : '
c         CALL OUTPUT(WORK(KCMO),1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
c         WRITE(lupri,*) '   '
c      call construct_r3(WORK(KCMO),PRT)
C --------------------------------------------------
      DO I = 1, size(LRATOMS)
      LRATOM = LRATOMS(I)

C     Next lines will be used for setting GO at each nucleus. Should be
C     done also for LR and QR.
c      call REGENERATE_INTEGRALS(CORD(:,LRATOM)) !Placing GO at each nucleus
c      NEWPRP = .TRUE.  !So that AOPROPER is re-written
c      HRINPC = .FALSE. !So that HERMIT INPUTS IS READ AGAIN FROM DALTON.INP
c      CALL HERCTL(WORK,LWORK)

      IF (SIGMAD0S) THEN
         TRPLET = .FALSE.
         TRPFLG = .FALSE.
         CALL AVELRSC('ANGP',WORK,LWORK) ! PSO-OZ
         CALL AVELRSC('FCAV',WORK,LWORK) ! FC
         CALL AVELRSC('DIAK',WORK,LWORK) ! NSKE
         CALL AVELRSC('SDAV',WORK,LWORK) ! SD_ij
      END IF
C --------------------------------------------------
C     Paramagnetic corrections to Shielding:
C     First Order Singlet
      IF (SIGMAP1S)  THEN
         CALL GET_LR('PSOK',.FALSE.) ! ANGMOM | PSO-K
         CALL GET_LR('OZK ',.FALSE.) !   PSO  | OZ-KE
      END IF

C     Diamagnetic corrections to Shielding:
C     First Order Singlet
C
      IF (SIGMAD1S) THEN
         CALL GET_LR('DIAM',.FALSE.) ! NSCO | Mv
         CALL GET_LR('DIAD',.FALSE.) ! NSCO | Dw
      END IF
C
C
C     Paramagnetic corrections to Shielding:
C     First Order TRIPLET
      IF (SIGMAP1T) THEN
         CALL GET_LR('FCZK',.TRUE.) !   FC   | KINENERG
         CALL GET_LR('FCDE',.TRUE.) !   FC   | dd/didj
         CALL GET_LR('SDK ',.TRUE.) !  SD_ij | KINENERG
         CALL GET_LR('SDD ',.TRUE.) !   SD   | dd/d d
         CALL GET_LR('FCBS',.TRUE.) !   FC   | SOMF_ij
         CALL GET_LR('SDBS',.TRUE.) !   SD   | SOMF
c         CALL GET_LR('SDK2',.TRUE.)
c         CALL GET_LR('SDD2',.TRUE.)
      END IF

c     Paramagnetic corrections to Shielding:
c     Third Order

c     Singlets
      IF (SIGMAP3S) THEN
         IF (JIMPRT.GT.1) THEN
         WRITE (LUPRI,'(A)') ' Calling to SIGMAP3S for : Paramagnetic'//
     &                     ' Singlet Third Order Shielding corrections'
         END IF
         CALL QUADRATICLR('OZ-PSO-MV',WORK,LWORK) !OK
         CALL QUADRATICLR('OZ-PSO-Dw',WORK,LWORK) !OK
      END IF

c     Triplets
      IF (SIGMAP3T) THEN
         IF (JIMPRT.GT.1) THEN
         WRITE (LUPRI,'(A)')' Calling to SIGMAP3T for : Paramagnetic'//
     &                     ' Triplet Third Order Shielding corrections'
         END IF
         CALL QUADRATICLR('OZ-FC-SO ',WORK,LWORK) !OK
         CALL QUADRATICLR('OZ-SD-SO ',WORK,LWORK) !OK
      END IF

c      Null (as expected)!
c       CALL QUADRATICLR('SZ-PSO-MV',WORK,LWORK)
c       CALL QUADRATICLR('SZ-PSO-Dw',WORK,LWORK)

c      Null (shouldnt be)
c       CALL QUADRATICLR('SZ-FC-SO ',WORK,LWORK) ! Según paper 2011, esto da cero
c       CALL QUADRATICLR('SZ-SD-SO ',WORK,LWORK) !ONIT



      SLRESC  = 0.0
      SGD0S   = 0.0  ! sigma 0order singlet DIAM  : Fc, AngPso, DiaKin
      SGD1S   = 0.0  ! sigma 1order singlet DIAM  : DiaMv, DiaDw
      SGP1S   = 0.0  ! sigma 1order singlet PARAM : Lkin, Psokin
      SGP1T   = 0.0  ! sigma 1order triplet PARAM : FcKin, SdKin, SdBso, FcBso
      SGP3S   = 0.0  ! sigma 3order singlet PARAM : PSOOZMv, PSOOZDw
      SGP3T   = 0.0  ! sigma 3order triplet PARAM : OZFCSO, OZSDSO
C     Add constants and print all in ppm's
!     Diamagnetic corrections
         LRFCAV = calfa*calfa*calfa*calfa*1.0D+6* CFCAV*LRFCAV
         LRDIAK = calfa*calfa*calfa*calfa*1.0D+6* CDIAK*LRDIAK
         LRSDAV = calfa*calfa*calfa*calfa*1.0D+6* CSDAV*LRSDAV
         LRANGP = calfa*calfa*calfa*calfa*1.0D+6* CANGP*LRANGP
         LRDIAD = calfa*calfa*1.0D+6* CDIAD* LRDIAD
         LRDIAM = calfa*calfa*1.0D+6* CDIAM* LRDIAM
!     Paramagnetic corrections
         LRFCZK = calfa*calfa*calfa*calfa*1.0D+6* CFCZK*LRFCZK
         LRFCDE = calfa*calfa*calfa*calfa*1.0D+6* CFCDE*LRFCDE
         LRSDK = calfa*calfa*calfa*calfa*1.0D+6* CSDK*LRSDK
c         LRSDK2 = calfa*calfa*calfa*calfa*1.0D+6* CSDK2*LRSDK2
         LRSDD = calfa*calfa*calfa*calfa*1.0D+6* CSDD*LRSDD
c         LRSDD2 = calfa*calfa*calfa*calfa*1.0D+6* CSDD2*LRSDD2
         LRFCBS = calfa*calfa*calfa*calfa*1.0D+6* CFCBS*LRFCBS
         LRSDBS = calfa*calfa*calfa*calfa*1.0D+6* CSDBS*LRSDBS
         LROZK  = calfa*calfa*calfa*calfa*1.0D+6* COZK*LROZK
         LRPSOK = calfa*calfa*calfa*calfa*1.0D+6* CPSOK*LRPSOK
c
         QRPSOOZMV = calfa*calfa*1.0D+6* CPSOOZMV*QRPSOOZMV
         QRPSOOZDW = calfa*calfa*1.0D+6* CPSOOZDW*QRPSOOZDW
         QROZFCSO =  calfa*calfa*1.0D+6* COZFCSO*QROZFCSO
         QROZSDSO =  calfa*calfa*1.0D+6* COZSDSO*QROZSDSO
C

!     Diamagnetic corrections
         SGD0S = SGD0S + LRFCAV + LRDIAK + LRANGP + LRSDAV
         SGD1S = SGD1S + LRDIAD + LRDIAM
!     Paramagnetic corrections
         SGP1S = SGP1S + LROZK + LRPSOK
         SGP1T = SGP1T + LRFCZK + LRFCDE + LRSDK
     &         + LRSDD  + LRFCBS + LRSDBS
         SGP3S = SGP3S + QRPSOOZMV + QRPSOOZDW
         SGP3T = SGP3T + QROZFCSO + QROZSDSO
!     Total
      SLRESC  =  SGD0S + SGD1S + SGP1S + SGP1T + SGP3S + SGP3T
C
C     LRESC-SHIELDING results section
      call PRINT_LRSHIELD
      END DO

c     SIRIFC has to remain open after all
      IF (SIGMAD1S.OR.SIGMAP1S.OR.SIGMAP1T) THEN
      IF ((.NOT.SIGMAP3S).OR.(.NOT. SIGMAP3T)) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',' ',IDUMMY, !BORRAR
     &                     .FALSE.)
      END IF
      END IF

      call move_file(LUESC2,LUESC,'DALTON_BACKUP.INP', 'DALTON.INP')
      deallocate(LRATOMS)
      CALL QEXIT('LRSCDRV_SHI')
      RETURN
      END
c    ---------------------------------------
c    ---------------------------------------
