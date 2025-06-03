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
      USE lrescinf_module
c#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "abainf.h"
      PARAMETER (NTABLE = 16, D0 = 0.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
C#include "dorps.h"
#include "spnout.h"
C                      1         2         3	    4	      5
      DATA TABLE /'.GAUGEO','.PRINT ','.SELECT','.PRTALL','.xXXXXx',
C                      6         7         8	    9	      10
     &            '.PARA1S','.PARA1T','.PARA3S','.PARA3T','.DIAM0S',
C                      11       12        13       14         15
     &            '.DIAM1S','.ANISO ','.PARANR','.DIANR ','.PYINTS',
C                      16
     &            '.EXTRAI'/
C
      NEWDEF = (WORD .EQ. '*LROPTS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
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
C              Redefining Gauge Origin
 1             CONTINUE
               CALL QUIT('(.GAUGEO) was deprecated.')
C                  READ (LUCMD,*) (LRGAUG(IS), IS = 1, 3)
c                  GAUCHANG =.TRUE.
                  ICHANG = ICHANG + 1
                  GOTO 100
C              SET PRINT LEVEL for LRESC-SHIELDING
 2             CONTINUE
                  READ (LUCMD,*) JIMPRT
                  ICHANG = ICHANG + 1
                  GOTO 100
C              TO BE DEVELOPED. It's defined under *LR-SHI
 3             CONTINUE
                  READ (LUCMD,*) LRATOM
                  ICHANG = ICHANG + 1
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
c                 Free option
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
                  SIGMAP2S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 9             CONTINUE
                  SIGMAP2T = .TRUE.
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
 13            CONTINUE
                  SIGMAPNR = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 14            CONTINUE
                  SIGMADNR = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 15            CONTINUE
                  PYINTS = .TRUE.
                  ICHANG = ICHANG + 1
                  GOTO 100
 16            CONTINUE
                  EXTRAINFO = .TRUE.
                  ICHANG = ICHANG + 1
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


      RETURN
      END
C -----------------------------------------------------------------------
C
C  /* Deck lrscini_shi */
      SUBROUTINE LRSCINI_SHI
C
C Purpose:
C     Initialize the variables: lrescinf_module.
C     This routine is call by ABAINIALL : abacus/abadrv.F
C Author:
C     J. I. Melo
C Edited:
C     J. J. Aucar. 2021
C
      USE lrescinf_module
#include "implicit.h"
#include "mxcent.h"
#include "cbiqr.h"
#include "abainf.h"
C
      PARAMETER(D0 = 0.0D+0)

      JIMPRT    = 0


      LRGAUG(1) = D0
      LRGAUG(2) = D0
      LRGAUG(3) = D0

C      GAUCHANG = .FALSE.

      SIGMAPNR = .FALSE.
      SIGMADNR = .FALSE.
      SIGMAP1S = .FALSE.
      SIGMAP1T = .FALSE.
      SIGMAD1S = .FALSE.
      SIGMAD0S = .FALSE.
      SIGMAP2S = .FALSE.
      SIGMAP2T = .FALSE.

      LRESCALL = .TRUE.

      PYINTS    = .FALSE.
      EXTRAINFO = .FALSE.

      WARNLRS = .FALSE.
      WARNLRT = .FALSE.
      WARNQR  = .FALSE.
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
!     LRESC corrections
      LRDNR  = D0
      LRFCAV = D0
!
      LRFCAV0 = D0
      LRFCAV1 = D0
      LRFCAV2 = D0
      LRFCAV3 = D0
!
      LRDIAK = D0
      LRANGP = D0
      LRANGPavg = D0
      LRSDAV = D0

      LRDIAM = D0
      LRDIAD = D0

      LRPNR  = D0
      LROZK  = D0
      LRPSOK = D0
      LRPSOKaba = D0
      LRFCZK = D0
      LRFCDE = D0
      LRSDK = D0
      LRSDD = D0
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
      USE lrescinf_module
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "nuclei.h"
#include "chrnos.h"              
#include "gnrinf.h"
#include "cbiher.h"
C
C
      LOGICAL PING, FOUND
      DIMENSION WORK(LWORK)
      logical :: file_exists, file_open
      logical :: warning_found
      INTEGER, allocatable :: LRATOMS(:)
      real(8) ::  array_oneatom(3,3,25)
      real(8), allocatable :: lresc_array(:,:,:,:)
      PING=.FALSE.
      CALL QENTER('LRSCDRV_SHI')

C ----------------------------
      IF (LRESCALL) THEN
         SIGMADNR = .TRUE.
         SIGMAPNR = .TRUE.
         SIGMAD0S = .TRUE.
         SIGMAD1S = .TRUE.
         SIGMAP1S = .TRUE.
         SIGMAP1T = .TRUE.
         SIGMAP2S = .TRUE.
         SIGMAP2T = .TRUE.
      ENDIF

C     Print information about LRESC metodology
      CALL HEADER_LRSHIELD
C
C     Diamagnetic corrections to Shielding:
C     Zeroth Order, Singlet of course
C
c     Nuclei list
      IF ( (LRATOM_changed).and.(NLRATOMS.EQ.1) ) THEN
         allocate(LRATOMS(1))
         LRATOMS(1) = LRATOM
      ELSEIF (.not.LRATOM_changed ) THEN
         allocate(LRATOMS(NUCIND))
         DO I = 1, NUCIND
            LRATOMS(I)=I
         END DO
      ELSE
         allocate(LRATOMS(NLRATOMS))
         DO I = 1, NLRATOMS
            LRATOMS(I)=LRATOMSS(I)
         END DO
      END IF


      allocate(lresc_array(NLRATOMS,3,3,25))


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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     PySCF(libcint) integrals
      if (PYINTS) THEN
         call import_matrix(19*NUCIND,"shi")
      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     PSOK (AVG) integrals
      IF (SIGMAP1S)  THEN
         CALL construct_operator('PSOK') ! PSOK with AVG
      END IF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C --------------------------------------------------
c     Construct LR singlet inputs and run responses (all at once)
      IF (SIGMAPNR.OR.SIGMAP1S.OR.SIGMAD1S) THEN
         call LRINPUT_ALL(LRATOMS,SIZE(LRATOMS),.FALSE.)
         CALL GETTIM(TIMLRS,WALLRS)
         IF (JIMPRT.GE.1)
     &   WRITE (LUPRI,'(5X,A)') 'Singlet Linear Responses'//
     &                       ' are under calculation'
c        Run Linear Responses
         LUPRI = 6232
         OPEN(UNIT=LUPRI, FILE='LRESC-LR-S',status='unknown',
     &        position='append')
         call RSPDRV(WORK,LWORK)
         CLOSE(LUPRI)
         LUPRI  = 2

C        Check for Warnings
         open(unit=6232, file=trim('LRESC-LR-S'), status='old')
         call SearchWord(6232, 'WARNING', warning_found)
         if (warning_found) WARNLRS = .true.
         close(6232)

         CALL GETTIM(TEND,WEND)
         TIMLRS = TEND - TIMLRS
         WALLRS = WEND - WALLRS
         IF (JIMPRT .GE. 2) THEN
            CALL TIMTXT('    Total CPU  time used in LRESC-LR-S:',
     &                       TIMLRS,LUPRI)
            CALL TIMTXT('    Total wall time used in LRESC-LR-S:',
     &                       WALLRS,LUPRI)
         END IF

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
         CALL GETTIM(TIMLRT,WALLRT)
         IF (JIMPRT.GE.1)
     &   WRITE (LUPRI,'(5X,A)') 'Triplet Linear Responses'//
     &                       ' are under calculation'
c        Run Linear Responses
         LUPRI = 6232
         OPEN(UNIT=LUPRI, FILE='LRESC-LR-T',status='unknown',
     &        position='append')
         call RSPDRV(WORK,LWORK)
         CLOSE(LUPRI)
         LUPRI  = 2

C        Check for Warnings
         open(unit=6232, file=trim('LRESC-LR-T'), status='old')
         call SearchWord(6232, 'WARNING', warning_found)
         if (warning_found) WARNLRT = .true.
         close(6232)

         CALL GETTIM(TEND,WEND)
         TIMLRT = TEND - TIMLRT
         WALLRT = WEND - WALLRT
         IF (JIMPRT .GE. 2) THEN
            CALL TIMTXT('    Total CPU  time used in LRESC-LR-T:',
     &                       TIMLRT,LUPRI)
            CALL TIMTXT('    Total wall time used in LRESC-LR-T:',
     &                       WALLRT,LUPRI)
         END IF

      END IF
C --------------------------------------------------
      DO I = 1, size(LRATOMS)
      LRATOM = LRATOMS(I)

C     Next lines will be used for setting GO at each nucleus. Should be
C     done also for LR and QR.
c      call REGENERATE_INTEGRALS(CORD(:,LRATOM)) !Placing GO at each nucleus
c      NEWPRP = .TRUE.  !So that AOPROPER is re-written
c      HRINPC = .FALSE. !So that HERMIT INPUTS IS READ AGAIN FROM DALTON.INP
c      CALL HERCTL(WORK,LWORK)
c      HRINPC = .TRUE.

      IF (SIGMADNR) THEN
         TRPLET = .FALSE.
         TRPFLG = .FALSE.
         CALL AVELRSC('DINR',WORK,LWORK) ! NSCO
      END IF

      IF (SIGMAD0S) THEN
         TRPLET = .FALSE.
         TRPFLG = .FALSE.
         IF (EXTRAINFO) CALL AVELRSC('ANGP',WORK,LWORK) ! PSO-OZ
         CALL AVELRSC('NEWJ',WORK,LWORK) ! PSO-OZ with AVG
         CALL AVELRSC('FCAV',WORK,LWORK) ! FC
         CALL AVELRSC('DIAK',WORK,LWORK) ! NSKE
         CALL AVELRSC('SDAV',WORK,LWORK) ! SD_ij
!      END IF
         IF (PYINTS) THEN
            CALL AVELRSC('FCt0',WORK,LWORK) ! FC type (FC integrals)
            CALL AVELRSC('FCt1',WORK,LWORK) ! FC type # (r . nabla) / (r^3)
            CALL AVELRSC('FCt2',WORK,LWORK) ! FC type # (r_i . nabla_j) /(r^3)
            CALL AVELRSC('FCt3',WORK,LWORK) ! FC type # (r_i r_j)/(r^5) r.nabla
         END IF

      END IF
C --------------------------------------------------
C     Paramagnetic corrections to Shielding:
C     Non-relativistic singlet
      IF (SIGMAPNR)  THEN
         CALL GET_LR('NRPA',.FALSE.) ! PSO | ANGMOM
      END IF
C     First Order Singlet
      IF (SIGMAP1S)  THEN
         CALL GET_LR('PSOK',.FALSE.) ! ANGMOM | PSO-K with AVG
         CALL GET_LR('OZK ',.FALSE.) !   PSO  | OZ-KE
         IF (EXTRAINFO) THEN
            CALL GET_LR('PSKK',.FALSE.) ! ANGMOM | PSO-K with abacus integrals
         END IF
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
      END IF

c     Paramagnetic corrections to Shielding:
c     Third Order

c     Singlets
      IF (SIGMAP2S) THEN
         IF (JIMPRT.GE.1) THEN
         WRITE (LUPRI,'(5X,A)') 'Singlet quadratic responses'//
     &                     ' are under calculation'
         END IF
         CALL QUADRATICLR('OZ-PSO-MV',WORK,LWORK) !OK
         CALL QUADRATICLR('OZ-PSO-Dw',WORK,LWORK) !OK
      END IF

c     Triplets
      IF (SIGMAP2T) THEN
         IF (JIMPRT.GE.1) THEN
         WRITE (LUPRI,'(5X,A)') 'Singlet quadratic responses'//
     &                     ' are under calculation'
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

C     Check for Warnings
      if ( (SIGMAP2S).or.(SIGMAP2T) ) then
         open(unit=6232, file=trim('LRESC-QR'), status='old')
         call SearchWord(6232, 'WARNING', warning_found)
         if (warning_found) WARNQR = .true.
         close(6232)
      end if


C --Warnings
      IF ( (WARNLRS).OR.(WARNLRT).OR.(WARNQR) )
     & CALL HEADER('Warnings :',3)
      IF (WARNLRS) WRITE (LUPRI,'(4X,A)')
     &      '-  Check LRESC-LR-S file.'
      IF (WARNLRT) WRITE (LUPRI,'(4X,A)')
     &      '-  Check LRESC-LR-T file.'
      IF (WARNQR) WRITE (LUPRI,'(4X,A)')
     &      '-  Check LRESC-QR file.'




      SLRESC  = 0.0
      SGPNR   = 0.0  ! non-relativistic paramagnetic
      SGDNR   = 0.0  ! non-relativistic diamagnetic
      SGD0S   = 0.0  ! sigma 0order singlet DIAM  : Fc, AngPso, DiaKin
      SGD1S   = 0.0  ! sigma 1order singlet DIAM  : DiaMv, DiaDw
      SGP1S   = 0.0  ! sigma 1order singlet PARAM : Lkin, Psokin
      SGP1T   = 0.0  ! sigma 1order triplet PARAM : FcKin, SdKin, SdBso, FcBso
      SGP3S   = 0.0  ! sigma 3order singlet PARAM : PSOOZMv, PSOOZDw
      SGP3T   = 0.0  ! sigma 3order triplet PARAM : OZFCSO, OZSDSO
C     Add constants and print all in ppm's
!     Diamagnetic non-relativistic value
         LRDNR = calfa*calfa*1.0D+6*CDNR*LRDNR
!     Diamagnetic corrections
         LRFCAV = calfa*calfa*calfa*calfa*1.0D+6* CFCAV*LRFCAV
         LRDIAK = calfa*calfa*calfa*calfa*1.0D+6* CDIAK*LRDIAK
         LRSDAV = calfa*calfa*calfa*calfa*1.0D+6* CSDAV*LRSDAV
         LRANGP = calfa*calfa*calfa*calfa*1.0D+6* CANGP*LRANGP
         LRANGPavg = calfa*calfa*calfa*calfa*1.0D+6* CANGP*LRANGPavg
         LRDIAD = calfa*calfa*1.0D+6* CDIAD* LRDIAD
         LRDIAM = calfa*calfa*1.0D+6* CDIAM* LRDIAM
!
         LRFCAV0 = calfa*calfa*calfa*calfa*1.0D+6* CFCAV0*LRFCAV0
         LRFCAV1 = calfa*calfa*calfa*calfa*1.0D+6* CFCAV1*LRFCAV1
         LRFCAV2 = calfa*calfa*calfa*calfa*1.0D+6* CFCAV2*LRFCAV2
         LRFCAV3 = calfa*calfa*calfa*calfa*1.0D+6* CFCAV3*LRFCAV3
         LRFCAVpy = LRFCAV0+LRFCAV1+LRFCAV2+LRFCAV3 !Imported from pyscf
!
!     Paramagnetic non-relativistic value
         LRPNR = calfa*calfa*1.0D+6*CPNR*LRPNR
!     Paramagnetic corrections
         LRFCZK = calfa*calfa*calfa*calfa*1.0D+6* CFCZK*LRFCZK
         LRFCDE = calfa*calfa*calfa*calfa*1.0D+6* CFCDE*LRFCDE
         LRSDK = calfa*calfa*calfa*calfa*1.0D+6* CSDK*LRSDK
         LRSDD = calfa*calfa*calfa*calfa*1.0D+6* CSDD*LRSDD
         LRFCBS = calfa*calfa*calfa*calfa*1.0D+6* CFCBS*LRFCBS
         LRSDBS = calfa*calfa*calfa*calfa*1.0D+6* CSDBS*LRSDBS
         LROZK  = calfa*calfa*calfa*calfa*1.0D+6* COZK*LROZK
         LRPSOK = calfa*calfa*calfa*calfa*1.0D+6* CPSOK*LRPSOK
         LRPSOKaba = calfa*calfa*calfa*calfa*1.0D+6* CPSOKaba*LRPSOKaba
c
         QRPSOOZMV = calfa*calfa*1.0D+6* CPSOOZMV*QRPSOOZMV
         QRPSOOZDW = calfa*calfa*1.0D+6* CPSOOZDW*QRPSOOZDW
         QROZFCSO =  calfa*calfa*1.0D+6* COZFCSO*QROZFCSO
         QROZSDSO =  calfa*calfa*1.0D+6* COZSDSO*QROZSDSO
C

!     Diamagnetic corrections
         IF (PYINTS) THEN
            SGD0S = SGD0S + LRFCAVpy + LRDIAK + LRANGPavg + LRSDAV
         ELSE
            SGD0S = SGD0S + LRFCAV + LRDIAK + LRANGPavg + LRSDAV
         END IF
         SGD1S = SGD1S + LRDIAD + LRDIAM
!     Paramagnetic corrections
         SGPNR = LRPNR
         SGDNR = LRDNR
         SGP1S = SGP1S + LROZK + LRPSOK
         SGP1T = SGP1T + LRFCZK + LRFCDE + LRSDK
     &         + LRSDD  + LRFCBS + LRSDBS
         SGP3S = SGP3S + QRPSOOZMV + QRPSOOZDW
         SGP3T = SGP3T + QROZFCSO + QROZSDSO
!     Total (doesnt include non-rel contributions)
      SLRESC  = SGD0S + SGD1S + SGP1S + SGP1T + SGP3S + SGP3T
C
C     LRESC-SHIELDING results section
      call PRINT_LRSHIELD


C     Save results to array
      call create_array_oneatom(array_oneatom)
      IF (NLRATOMS.EQ.1) lresc_array(1, :, :,:) = array_oneatom
      IF (NLRATOMS.NE.1) lresc_array(LRATOM, :, :,:) = array_oneatom

      END DO

C     Export results to binary file
      call binary_lresc(lresc_array)
      deallocate(lresc_array)


c     SIRIFC has to remain open after all
      IF (SIGMAD1S.OR.SIGMAPNR.OR.SIGMAP1S.OR.SIGMAP1T) THEN
      IF ((.NOT.SIGMAP2S).OR.(.NOT. SIGMAP2T)) THEN
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
