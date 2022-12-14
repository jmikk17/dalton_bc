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
      REAL*8 LR(3,3),TOTALLR(3,3)
      REAL*8 LRtensor(3,3)
      CHARACTER*4 :: KEY
      CHARACTER*8, allocatable :: operadores_A(:,:)
      CHARACTER*8, allocatable :: operadores_B(:,:)
      INTEGER                  :: NopA, NopB
      CHARACTER*8, dimension(2) :: operadores_duo
      INTEGER :: SUMA_B
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number
      LOGICAL :: TRPLET_FLAG

      CALL QENTER('LRSLINEAR')
c      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL GETTIM(TIMHER,WALHER)

c     NopA: number of A operators for each shield tensor element
c     NopB: number of B operators for each shield tensor element
c      LUDV   = N2ASHX
c      LPVX   = 0

c      LUSOVE = 456
c      LUGDVE = 457
c      LUREVE = 458
c      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
c      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
c      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)

C  LINEAR LRESC SINGLET ROUTINE
C   KEY = 'DIAM'  DIAMAGNETIC <<DIA,MV>>
C         'DIAD'  DIAMAGNETIC <<DIA,DW>>
C         'PSOK'  PARAMAGNETIC with built-in integrals
C         'PSKI'  PARAMAGNETIC with angkin routine <<PSO.1.Kinenerg;PSO>>
C         'LKIN'  PARAMAGNETIC
C         'FCBS' 'SDBS'
C         'FCZK' 'FCZK'
      IF (JIMPRT .GE. 1) THEN
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
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                    ij
         operadores_A(1,1) = nuc_coord(1)//'NSCO'//CHRXYZ(1) !11
         operadores_A(2,1) = nuc_coord(1)//'NSCO'//CHRXYZ(2) !12
         operadores_A(3,1) = nuc_coord(1)//'NSCO'//CHRXYZ(3) !13
         operadores_A(4,1) = nuc_coord(2)//'NSCO'//CHRXYZ(1) !21
         operadores_A(5,1) = nuc_coord(2)//'NSCO'//CHRXYZ(2) !22
         operadores_A(6,1) = nuc_coord(2)//'NSCO'//CHRXYZ(3) !23
         operadores_A(7,1) = nuc_coord(3)//'NSCO'//CHRXYZ(1) !31
         operadores_A(8,1) = nuc_coord(3)//'NSCO'//CHRXYZ(2) !32
         operadores_A(9,1) = nuc_coord(3)//'NSCO'//CHRXYZ(3) !33
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         IF (KEY.EQ."DIAD") operadores_B(:,1) ='DARWIN  '
         IF (KEY.EQ."DIAM") operadores_B(:,1) ='MASSVELO'
      END IF


      IF(KEY.EQ."PSOK") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                    ij
         operadores_A(1,1) = 'PSOKE'//nuc_coord(1)           !11 
         operadores_A(2,1) = 'PSOKE'//nuc_coord(1)           !12 editar
         operadores_A(3,1) = 'PSOKE'//nuc_coord(1)           !13 editar
         operadores_A(4,1) = 'PSOKE'//nuc_coord(1)           !21 editar
         operadores_A(5,1) = 'PSOKE'//nuc_coord(2)           !22
         operadores_A(6,1) = 'PSOKE'//nuc_coord(1)           !23 editar
         operadores_A(7,1) = 'PSOKE'//nuc_coord(1)           !31 editar
         operadores_A(8,1) = 'PSOKE'//nuc_coord(1)           !32 editar
         operadores_A(9,1) = 'PSOKE'//nuc_coord(3)           !33
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         operadores_B(1,1) = 'XANGMOM '                      !11 
         operadores_B(2,1) = 'XANGMOM '                      !12 
         operadores_B(3,1) = 'XANGMOM '                      !13 
         operadores_B(4,1) = 'XANGMOM '                      !21 
         operadores_B(5,1) = 'YANGMOM '                      !22 
         operadores_B(6,1) = 'XANGMOM '                      !23 
         operadores_B(7,1) = 'XANGMOM '                      !31 
         operadores_B(8,1) = 'XANGMOM '                      !32 
         operadores_B(9,1) = 'ZANGMOM '                      !33 
      END IF

      IF(KEY.EQ."LKIN") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .FALSE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                    ij
         operadores_A(1,1) = 'PSO '//nuc_coord(1)//' '       !11 
         operadores_A(2,1) = 'PSO '//nuc_coord(1)//' '       !12
         operadores_A(3,1) = 'PSO '//nuc_coord(1)//' '       !13
         operadores_A(4,1) = 'PSO '//nuc_coord(2)//' '       !21
         operadores_A(5,1) = 'PSO '//nuc_coord(2)//' '       !22
         operadores_A(6,1) = 'PSO '//nuc_coord(2)//' '       !23
         operadores_A(7,1) = 'PSO '//nuc_coord(3)//' '       !31
         operadores_A(8,1) = 'PSO '//nuc_coord(3)//' '       !32
         operadores_A(9,1) = 'PSO '//nuc_coord(3)//' '       !33
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         operadores_B(1,1) = CHRXYZ(1)//'OZKE   '            !11 
         operadores_B(2,1) = CHRXYZ(2)//'OZKE   '            !12 
         operadores_B(3,1) = CHRXYZ(3)//'OZKE   '            !13 
         operadores_B(4,1) = CHRXYZ(1)//'OZKE   '            !21 
         operadores_B(5,1) = CHRXYZ(2)//'OZKE   '            !22 
         operadores_B(6,1) = CHRXYZ(3)//'OZKE   '            !23 
         operadores_B(7,1) = CHRXYZ(1)//'OZKE   '            !31 
         operadores_B(8,1) = CHRXYZ(2)//'OZKE   '            !32 
         operadores_B(9,1) = CHRXYZ(3)//'OZKE   '            !33 
      END IF

      IF(KEY.EQ."FCZK") THEN
         NopA = 1
         NopB = 1
         TRPLET_FLAG = .TRUE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                          ij
         operadores_A(:,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number  !all
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         operadores_B(:,1) = 'KINENERG'                            !all
      END IF

c     Under construction
      IF(KEY.EQ."SDZK") THEN
         NopA = 2
         NopB = 1
         TRPLET_FLAG = .TRUE.
C        -----------------------------------
         allocate(operadores_A(9,NopA))
c                  3*i+j-3                                          ij
         operadores_A(1,1) = 'SD '//nuc_coord(1)//' x'  !11
C        -----------------------------------
         allocate(operadores_B(9,NopB))
         operadores_B(:,1) = 'KINENERG'                            !all
      END IF

c     Elemento ij del tensor de shielding
      do i = 1, 3 
      do j = 1, 3
c         LR = 0.d0
c         TOTALLR = 0.d0
         do k = 1, NopA
         do l = 1, NopB
            operadores_duo(1) = operadores_A(3*i+j-3,k)
            operadores_duo(2) = operadores_B(3*i+j-3,l)

            call create_LRinput(operadores_duo,TRPLET_FLAG)
            LUPRI = 6232
            OPEN(UNIT=LUPRI, FILE='LRESC-LR',status='unknown')
c            OPEN(UNIT=LUPRI, FILE='LRESC-LR '//KEY//
c     &           ' '//NAMN(LRATOM),status='unknown')
            call RSPDRV(WORK,LWORK)
            CLOSE(LUPRI)
            LUPRI  = 2
         
c            CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
c            CALL GPCLOSE(LURSP3,'DELETE')
c            CALL GPCLOSE(LURSP4,'DELETE')
c            CALL GPCLOSE(LURSP5,'DELETE')
            call get_LRvalue(LR(i,j),operadores_duo) !LR from here has an extra minus sign
c            WRITE(LUPRI,*) "AUCAR-LR",operadores_duo,i,j,calfa*calfa*1.0D6*LR
c            WRITE(LUPRI,*) "AUCAR-LR",calfa*calfa*1.0D6*LR
            TOTALLR(i,j) = TOTALLR(i,j) + LR(i,j)
            IF (JIMPRT .GT. 2) write(lupri,*) operadores_duo(1),
     &         ' ',operadores_duo(2)," LR: ",LR(i,j)
         end do
         end do

         LRtensor(i,j) = TOTALLR(i,j)
         IF (JIMPRT .GT. 2) THEN
            write(lupri,*) KEY,"(",i,",",j,"):",LRtensor(i,j)
            write(lupri,*)
         END IF

      end do
      end do

      IF (KEY.EQ.'DIAM') LRDIAM = LRtensor
      IF (KEY.EQ.'DIAD') LRDIAD = LRtensor
      IF (KEY.EQ.'PSOK') LRPSOK = LRtensor
      IF (KEY.EQ.'LKIN') LRLKIN = LRtensor
      IF (KEY.EQ.'FCZK') LRFCZK = LRtensor

      deallocate(operadores_A)
      deallocate(operadores_B)

c      CALL GPCLOSE(LUSOVE,'DELETE')
c      CALL GPCLOSE(LUGDVE,'DELETE')
c      CALL GPCLOSE(LUREVE,'DELETE')

c      CALL TIMER ('LRSLINEAR '//KEY//' '//NAMN(LRATOM),TIMSTR,TIMEND)
      CALL GETTIM(TEND,WEND)
      TIMHER = TEND - TIMHER
      WALHER = WEND - WALHER
      CALL TIMTXT(' Total CPU  time used in '//KEY//' '//NAMN(LRATOM)//
     &            ' LR:',TIMHER,LUPRI)
      CALL TIMTXT(' Total wall time used in '//KEY//' '//NAMN(LRATOM)//
     &            ' LR:',WALHER,LUPRI)
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
c      write(fileID,'(a)') '.PRINT'
c      write(fileID,'(a)') ' 0'
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
      CHARACTER*164 STRING, TMPSTRING
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

 1    READ (fileID,'(A164)',END=3) STRING
c     Getting LR when its called for the first time
      IF (STRING(1:30) .EQ. '@ FREQUENCY INDEPENDENT SECOND') THEN
          READ (fileID,'(A164)') TMPSTRING !Empty line
          READ (fileID,'(A164)') TMPSTRING !LR of A,A
          READ (fileID,'(31X,F19.12)') LRVAL
      ELSE
         GOTO 2  
      END IF

c     Getting LR when its not called for the first time
 2    IF (STRING(2:25) .EQ. 'Linear response function') THEN
          READ (fileID,'(36X,A8)') labels_read(1)
          READ (fileID,'(36X,A8)') labels_read(2)
      IF (labels_read(1) .EQ. labels_in(1) .AND.
     &    labels_read(2) .EQ. labels_in(2) ) then
         READ (fileID,'(/,48X,F15.12)') LRVAL
      ELSE
         GOTO 1
      END IF

      ELSE
         GOTO 1
      END IF

 3           CONTINUE

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

      RETURN
      END
