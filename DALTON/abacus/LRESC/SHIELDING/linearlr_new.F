      SUBROUTINE LINEARLR_NEW(KEY,WORK,LWORK)

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
#include "chrxyz.h"

      DIMENSION WORK(LWORK)
      REAL*8 LR,TOTALLR
      REAL*8 LRtensor(3,3)
      CHARACTER*4 :: KEY
      CHARACTER*8, allocatable :: operadores_A(:)
      CHARACTER*8, allocatable :: operadores_B(:)
      CHARACTER*8 :: operadores_duo(2)
      INTEGER :: SUMA_B
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      LOGICAL :: TRPLET_FLAG

      CALL QENTER('LRSLINEAR')
      CALL TIMER('START ',TIMEIN,TIMOUT)


C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'DIAM'  DIAMAGNETIC <<DIA,MV>>
C         'DIAD'  DIAMAGNETIC <<DIA,DW>>
C         'PSOK'  PARAMAGNETIC with built-in integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'LKIN'  PARAMAGNETIC
C         'FCBS' 'SDBS'
C         'FCZK' 'FCZK'
      IF (JIMPRT .GE. 1) THEN
C        WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
C        WRITE(LUPRI,*)
         IF (KEY.EQ.'DIAM') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic Second Order Singlet (DIAM)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'DIAD') WRITE(LUPRI,'(5X,A,A)')
     &   'Diamagnetic Second Order Singlet (DIAD)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'PSOK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (PSOK)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'PSOKI') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (PSOKI)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'LKIN') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Singlet (LKIN)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'FCZK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCZK)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'SDZK') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDZK)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'FCBS') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (FCBS)',
     &   ' correction is under calculation...'
         IF (KEY.EQ.'SDBS') WRITE(LUPRI,'(5X,A,A)')
     &   'Paramagnetic Second Order Triplet (SDBS)',
     &   ' correction is under calculation...'
C         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
      END IF
      TOTALLR = 0.d0
      LR = 0.d0

      do i = 1, 3
         ICOORD = i !(x,y,z)
         ISCOOR = LRATOM*3 + (ICOORD - 3)
         IFIRST=ISCOOR/100
         ISECND = MOD(ISCOOR,100)/10
         ITHIRD = MOD(MOD(ISCOOR,100),10)
         nuc_coord(i)=CHRNOS(IFIRST)//CHRNOS(ISECND)//CHRNOS(ITHIRD)
      end do

      nuc_number=CHRNOS(LRATOM/100)//CHRNOS(MOD(LRATOM,100)/10)//
     &              CHRNOS(MOD(LRATOM,10))

      IF((KEY.EQ."DIAM").OR.(KEY.EQ."DIAD")) THEN
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(operadores_A(9))
c                  3*i+j -3                                 ij
         operadores_A(1) = nuc_coord(1)//'NSCO'//CHRXYZ(1) !11
         operadores_A(2) = nuc_coord(1)//'NSCO'//CHRXYZ(2) !12
         operadores_A(3) = nuc_coord(1)//'NSCO'//CHRXYZ(3) !13
         operadores_A(4) = nuc_coord(2)//'NSCO'//CHRXYZ(1) !21
         operadores_A(5) = nuc_coord(2)//'NSCO'//CHRXYZ(2) !22
         operadores_A(6) = nuc_coord(2)//'NSCO'//CHRXYZ(3) !23
         operadores_A(7) = nuc_coord(3)//'NSCO'//CHRXYZ(1) !31
         operadores_A(8) = nuc_coord(3)//'NSCO'//CHRXYZ(2) !32
         operadores_A(9) = nuc_coord(3)//'NSCO'//CHRXYZ(3) !33
C        -----------------------------------
         allocate(operadores_B(1))
         IF (KEY.EQ."DIAD") operadores_B(1) ='DARWIN  '
         IF (KEY.EQ."DIAM") operadores_B(1) ='MASSVELO'
      END IF



c      IF (KEY.EQ.'SZ-PSO-Dw') THEN
c         allocate(operadores(3,3))
c         allocate(spin_rank(3,3))
c         operadores(1,1) = 'X1SPNORB'
c         operadores(2,1) = 'Y1SPNORB'
c         operadores(3,1) = 'Z1SPNORB'
c         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
c         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
c         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
c         operadores(:,3) = 'DARWIN  '
c         spin_rank(:,1) = '1'
c         spin_rank(:,2) = '0'
c         spin_rank(:,3) = '0'
c      END IF

c      IF (KEY.EQ.'OZ-FC-SO ') THEN
c         allocate(operadores(3,3))
c         allocate(spin_rank(3,3))
c         operadores(1,1) = 'XANGMOM '
c         operadores(2,1) = 'YANGMOM '
c         operadores(3,1) = 'ZANGMOM '
c         spin_rank(:,1) = '0'
c         operadores(1,2) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
c         operadores(2,2) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
c         operadores(3,2) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
c         spin_rank(:,2) = '1'
c         operadores(1,3) = 'X1SPNORB'
c         operadores(2,3) = 'Y1SPNORB'
c         operadores(3,3) = 'Z1SPNORB'
c         spin_rank(:,3) = '1'
c      END IF

c      IF (KEY.EQ.'OZ-SD-SO ') THEN
c         allocate(operadores(9,3))
c         allocate(spin_rank(9,3))
c         operadores(1,1) = 'XANGMOM '
c         operadores(2,1) = 'YANGMOM '
c         operadores(3,1) = 'ZANGMOM '
c         spin_rank(:,1) = '0'
c         operadores(1,2) = 'SD '//nuc_coord(1)//' x'
c         operadores(2,2) = 'SD '//nuc_coord(1)//' y'
c         operadores(3,2) = 'SD '//nuc_coord(1)//' z'
c         operadores(4,2) = 'SD '//nuc_coord(2)//' x'
c         operadores(5,2) = 'SD '//nuc_coord(2)//' y'
c         operadores(6,2) = 'SD '//nuc_coord(2)//' z'
c         operadores(7,2) = 'SD '//nuc_coord(3)//' x'
c         operadores(8,2) = 'SD '//nuc_coord(3)//' y'
c         operadores(9,2) = 'SD '//nuc_coord(3)//' z'
c         spin_rank(:,2) = '1'
c         operadores(1,3) = 'X1SPNORB'
c         operadores(2,3) = 'Y1SPNORB'
c         operadores(3,3) = 'Z1SPNORB'
c         spin_rank(:,3) = '1'
c      END IF



c     Elemento ij del tensor de shielding
c      SUMA_B = SIZE(operadores,DIM=1)/3
      do i = 1, 3 
      do j = 1, 3
         LR = 0.d0
         TOTALLR = 0.d0

!DIA         
         operadores_duo(1) = operadores_A(3*i+j -3)
         operadores_duo(2) = operadores_B(1)
c         IF (operadores(1,3)(3:8).EQ.'SPNORB') THEN
c            k = 3
c         else
c            k = 1
c         END IF
      k = 1

c        IF (KEY(4:5).EQ.'SD') THEN
c           l = 3
c        else
c           l = 1
c        END IF
      l = 1

cc         do m = 1, k
cc            operadores_d(3)=operadores(m,3)
cc            spin_trio(3)=spin_rank(m,3)
c         do n = 1, l
c            operadores_duo(2) = operadores(SUMA_B*(j-1)+n,2)
c            operadores_duo(2) = operadores(j,2)
            call create_LRinput(operadores_duo,TRPLET_FLAG)
            LUPRI = 6232
            OPEN(UNIT=LUPRI, FILE='LRESC-LR',status='unknown')
            call RSPDRV(WORK,LWORK)
            CLOSE(LUPRI)
            LUPRI  = 2
         
c            CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
c            CALL GPCLOSE(LURSP3,'DELETE')
c            CALL GPCLOSE(LURSP4,'DELETE')
c            CALL GPCLOSE(LURSP5,'DELETE')
            call get_LRvalue(LR,operadores_duo)
c           las constantes tienen aca un (-) extra xq saco la respuesta
c           de lo que imprime como salida
c            WRITE(LUPRI,*) "AUCAR-LR",operadores_duo,i,j,calfa*calfa*1.0D6*LR
            WRITE(LUPRI,*) "AUCAR-LR",calfa*calfa*1.0D6*LR
c            TOTALLR = TOTALLR + LR
            IF (JIMPRT .GT. 2) write(lupri,*) operadores_duo(1),
     &         ' ',operadores_duo(2)," LR: ",LR
c         end do
cc         end do

c         LRtensor(i,j) = TOTALLR
         IF (JIMPRT .GT. 2) THEN
            write(lupri,*) KEY,"(",i,",",j,"):",LRtensor(i,j)
            write(lupri,*)
         END IF

      end do
      end do

c      IF (KEY.EQ.'OZ-PSO-MV') QRPSOOZMV = QRtensor
c      IF (KEY.EQ.'OZ-PSO-Dw') QRPSOOZDW = QRtensor
c      IF (KEY.EQ.'OZ-FC-SO ') QROZFCSO = QRtensor
c      IF (KEY.EQ.'OZ-SD-SO ') QROZSDSO = QRtensor

      deallocate(operadores_A)
      deallocate(operadores_B)


      IF (JIMPRT .GT. 1) CALL TIMER ('LRSLINEAR',TIMEIN,TIMOUT)
      CALL QEXIT('LRSLINEAR')

      RETURN

      END 


      subroutine create_LRinput(operadores,TRPLET_FLAG)
      character*10 :: namefile
      logical            :: file_exists, file_open
      logical, intent(in)     :: TRPLET_FLAG
      character(100) :: line
      INTEGER :: fileID
      CHARACTER*8, intent(in) :: operadores(2)

#include "priunit.h"
      namefile='DALTON.INP'
      fileID = -1
      CALL GPOPEN(fileID,namefile,'NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      write(fileID,'(a)') '**DALTON INPUT'
      write(fileID,'(a)') '**RESPONSE'
      write(fileID,'(a)') '*LINEAR'
      write(fileID,'(a)') '.PRINT'
      write(fileID,'(a)') ' 0'
      write(fileID,'(a)') '.PROPRT'
      write(fileID,'(a)') operadores(1)
      write(fileID,'(a)') '.PROPRT'
      write(fileID,'(a)') operadores(2)
      IF (TRPLET_FLAG) write(fileID,'(a)') '.TRIPLET'
      write(fileID,'(a)') '**END OF DALTON INPUT'

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

      RETURN
      END


      subroutine get_LRvalue(LRVAL,labels_in)
      logical :: file_exists, file_open
      INTEGER :: fileID
      CHARACTER*164 STRING
      CHARACTER*8, intent(in) :: labels_in(2)
      CHARACTER*8:: labels_read(2)
      REAL*8, intent(out) :: LRVAL
c

#include "priunit.h"

C Based on get_QRvalue

      inquire(file = 'LRESC-LR', exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            fileID = -1
            CALL GPOPEN(fileID,'LRESC-LR','OLD',' ','FORMATTED',
     &           IDUMMY,.FALSE.)
         end if
         rewind(fileID)
      else
         call quit('LRESC-LR file not found')
      end if

 346  READ (fileID,'(A164)',END=347) STRING
      IF (STRING(2:25) .EQ. 'Linear response function') THEN
          READ (fileID,'(36X,A8)') labels_read(1)
          READ (fileID,'(36X,A8)') labels_read(2)
      IF (labels_read(1) .EQ. labels_in(1) .AND.
     &    labels_read(2) .EQ. labels_in(2) ) then
         READ (fileID,'(/,49X,F15.8)') LRVAL
         WRITE(LUPRI,*) "AUCAR ",labels_read,LRVAL
      ELSE
         GOTO 346
      END IF

      ELSE
         GOTO 346
      END IF

 347         CONTINUE

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

      RETURN
      END
