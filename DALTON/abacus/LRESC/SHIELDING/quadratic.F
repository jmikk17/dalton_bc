      SUBROUTINE QUADRATICLR(KEY,WORK,LWORK)

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

      DIMENSION WORK(LWORK)
      logical :: file_exists, file_open
      character(100) :: line
      REAL*8 QR,TOTALQR,constante
      REAL*8 QRtensor(3,3)
      CHARACTER*9 :: KEY
      CHARACTER*8, allocatable :: operadores(:,:)
      CHARACTER*8 :: operadores_trio(3)
      CHARACTER*1 :: spin_trio(3)
      CHARACTER*1, allocatable :: spin_rank(:,:)
      INTEGER :: SUMA_B

      double precision, dimension(NBAST,NBAST) :: Identity,OVERLP,aux1
      CHARACTER*8 RTNLBL(2)

      double precision, dimension(NBAST*NORBT) :: xmos
      double precision   ::   cmos(NORBT,NBAST)
      integer :: naos,nmos
      LOGICAL FOUND
      character*3, dimension(3) :: nuc_coord
      character*3               :: nuc_number

c     operadores contiene en el primer indice la cantidad
c     de trio de operadores para las respuestas cuadráticas
      
      CALL QENTER('LRSCQUAD')
      CALL TIMER('START ',TIMEIN,TIMOUT)

      naos = NBAST ! NBAST: Total number of (atomic) basis functions
      nmos= NORBT !NORBT: Total number of molecular orbitals

      TOTALQR = 0.d0
      QR = 0.d0

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

      IF (KEY.EQ.'OZ-PSO-MV') THEN
c        constante = -1000000 * calfa**2 / 2
         allocate(operadores(3,3))
         allocate(spin_rank(3,3))
         operadores(1,1) = 'XANGMOM '
         operadores(2,1) = 'YANGMOM '
         operadores(3,1) = 'ZANGMOM '
         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
         operadores(1,3) = 'MASSVELO'
         operadores(2,3) = 'MASSVELO'
         operadores(3,3) = 'MASSVELO'
         spin_rank(:,:) = '0'
      END IF

      IF (KEY.EQ.'OZ-PSO-Dw') THEN
c        constante = 1000000/(2*(137.036**2))
         allocate(operadores(3,3))
         allocate(spin_rank(3,3))
         operadores(1,1) = 'XANGMOM '
         operadores(2,1) = 'YANGMOM '
         operadores(3,1) = 'ZANGMOM '
         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
         operadores(:,3) = 'DARWIN  '
         spin_rank(:,:) = '0'
      END IF


      IF (KEY.EQ.'SZ-PSO-MV') THEN
c        constante = 1000000/(137.036**2)
         allocate(operadores(3,3))
         allocate(spin_rank(3,3))
         operadores(1,1) = 'X1SPNORB'
         operadores(2,1) = 'Y1SPNORB'
         operadores(3,1) = 'Z1SPNORB'
         spin_rank(:,1) = '1'
         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
         spin_rank(:,2) = '0'
         operadores(:,3) = 'MASSVELO'
         spin_rank(:,3) = '0'
      END IF

      IF (KEY.EQ.'SZ-PSO-Dw') THEN
c        constante = 1000000/(137.036**2)
         allocate(operadores(3,3))
         allocate(spin_rank(3,3))
         operadores(1,1) = 'X1SPNORB'
         operadores(2,1) = 'Y1SPNORB'
         operadores(3,1) = 'Z1SPNORB'
         operadores(1,2) = 'PSO '//nuc_coord(1)//' '
         operadores(2,2) = 'PSO '//nuc_coord(2)//' '
         operadores(3,2) = 'PSO '//nuc_coord(3)//' '
         operadores(:,3) = 'DARWIN  '
         spin_rank(:,1) = '1'
         spin_rank(:,2) = '0'
         spin_rank(:,3) = '0'
      END IF

      IF (KEY.EQ.'OZ-FC-SO ') THEN
c        constante = 1000000/(4*(137.036**2))
         allocate(operadores(3,3))
         allocate(spin_rank(3,3))
         operadores(1,1) = 'XANGMOM '
         operadores(2,1) = 'YANGMOM '
         operadores(3,1) = 'ZANGMOM '
         spin_rank(:,1) = '0'
         operadores(1,2) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         operadores(2,2) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         operadores(3,2) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         spin_rank(:,2) = '1'
         operadores(1,3) = 'X1SPNORB'
         operadores(2,3) = 'Y1SPNORB'
         operadores(3,3) = 'Z1SPNORB'
         spin_rank(:,3) = '1'
      END IF

      IF (KEY.EQ.'OZ-SD-SO ') THEN
c        constante = 1000000/(4*(137.036**2))
         allocate(operadores(9,3))
         allocate(spin_rank(9,3))
         operadores(1,1) = 'XANGMOM '
         operadores(2,1) = 'YANGMOM '
         operadores(3,1) = 'ZANGMOM '
         spin_rank(:,1) = '0'
         operadores(1,2) = 'SD '//nuc_coord(1)//' x'
         operadores(2,2) = 'SD '//nuc_coord(1)//' y'
         operadores(3,2) = 'SD '//nuc_coord(1)//' z'
         operadores(4,2) = 'SD '//nuc_coord(2)//' x'
         operadores(5,2) = 'SD '//nuc_coord(2)//' y'
         operadores(6,2) = 'SD '//nuc_coord(2)//' z'
         operadores(7,2) = 'SD '//nuc_coord(3)//' x'
         operadores(8,2) = 'SD '//nuc_coord(3)//' y'
         operadores(9,2) = 'SD '//nuc_coord(3)//' z'
         spin_rank(:,2) = '1'
         operadores(1,3) = 'X1SPNORB'
         operadores(2,3) = 'Y1SPNORB'
         operadores(3,3) = 'Z1SPNORB'
         spin_rank(:,3) = '1'
      END IF

c     Está entrando dos veces a la condición, no?
      IF (KEY.EQ.'SZ-SD-SO '.OR. KEY.EQ.'SZ-FC-SO ') THEN
        CALL RD_SIRIFC('CMO',FOUND,xmos)
        cmos=0.0
        k = 1
        do i=1,nmos                  !i corresponds to a MO
          do j=1,naos                              !j corresponds to j coef. of the ith MO
            cmos(i,j) = xmos(k)                ! Transforms vector into matrix
c           write(lupri,*) "JUAN-cmos: ",i,j,cmos(i,j)
c           write(lupri,*)
            k = k + 1
          end do
        end do
c
        CALL GETDAT(RTNLBL(1),RTNLBL(2))
        RTNLBL(1)='1 A  21 '
C       Replace time information with symmetry information
        RTNLBL(2)='SQUARE  '
c        RTNLBL(2)='ANTISYMM'
        call WRITEAOPJUAN(matmul(transpose(cmos),cmos),NBAST,
     &        'IDENTIDA',RTNLBL) !Write Matrix in AOPROPER
c        write(LUPRI,*) "MATRIZ IDENTIDAD:"
c        call OUTPUT(matmul(transpose(cmos),cmos),1,NORBT,1,NORBT,
c     &         NORBT,NORBT,1,LUPRI) ! print matrix

        END IF
c
      IF (KEY.EQ.'SZ-FC-SO ') THEN
c        Constante: (137,036^-2/2)*1000000
c         constante = 1000000/(2*(137.036**2))
         allocate(operadores(3,3))
         allocate(spin_rank(3,3))
c         operadores(1,2) = 'X1SPNORB'
c         operadores(2,2) = 'Y1SPNORB'
c         operadores(3,2) = 'Z1SPNORB'
c         spin_rank(1,2) = '1'
c         spin_rank(2,2) = '1'
c         spin_rank(3,2) = '1'
         operadores(:,2) = 'IDENTIDA'
         spin_rank(:,2) = '1'

         operadores(1,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         operadores(2,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         operadores(3,1) = 'FC '//NAMN(LRATOM)(1:2)//nuc_number
         spin_rank(:,1) = '1'

         operadores(1,3) = 'X1SPNORB'
         operadores(2,3) = 'Y1SPNORB'
         operadores(3,3) = 'Z1SPNORB'
         spin_rank(:,3) = '1'

      END IF

      IF (KEY.EQ.'SZ-SD-SO ') THEN
c         CALL GETDAT(RTNLBL(1),RTNLBL(2))
c         RTNLBL(1)='1 A  21 '
C        Replace time information with symmetry information
c         RTNLBL(2)='SYMMETRI'
c         Identity = 0                          ! Initialize the array.
c         forall(j = 1:NBAST) Identity(j,j) = 1 ! Set the diagonal.
c         call WRITEAOPJUAN(Identity,NBAST,'IDENTITY',RTNLBL) !Write Matrix in AOPROPER

c         WRITE(lupri,*)' IDENTITY : '
c         CALL OUTPUT(Identity,1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
c         do i=1,NBAST
c         do j=1,NBAST
c            write(lupri,*) i,j,OVERLP(i,j)
c         enddo
c         enddo
c         WRITE(lupri,*) '   '

c         WRITE(lupri,*)' OVERLAP : '
c         call READAOP(OVERLP, NBAST, 'OVERLAP ') !gets OVERLAP matrix
c         CALL OUTPUT(OVERLP,1,NORBT,1,5,NORBT,NORBT,1,LUPRI)
c         WRITE(lupri,*)
c         do i=1,NBAST
c         do j=1,NBAST
c            write(lupri,*) i,OVERLP(i,i)
c         enddo
c         enddo


c         CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO))
cc         CALL RD_SIRIFC('CMO',FOUND,xmos)
cc                cmos=0.0
cc       k = 1
cc       do i=1,nmos                  !i corresponds to a MO
cc          do j=1,naos                              !j corresponds to j coef. of the ith MO
cc             cmos(i,j) = xmos(k)                ! Transforms vector into matrix
c             write(lupri,*) "JUAN-cmos: ",i,j,cmos(i,j)
c             write(lupri,*)
cc            k = k + 1
cc         end do
cc       end do

c         do i=1,NBAST
c         do j=1,NBAST
c            write(lupri,*) "JUAN-Ident: ",
c     &                                i,j,matmul(transpose(cmos),cmos)!(i,j)
c             write(lupri,*)
c         enddo
c         enddo


cc         CALL GETDAT(RTNLBL(1),RTNLBL(2))
cc         RTNLBL(1)='1 A  21 '
C        Replace time information with symmetry information
cc         RTNLBL(2)='SQUARE  '
c        RTNLBL(2)='ANTISYMM'
cc      call WRITEAOPJUAN(matmul(transpose(cmos),cmos),NBAST,
cc    &        'IDENTIDA',RTNLBL) !Write Matrix in AOPROPER
c        Constante: (137,036^-2/2)*1000000
c         constante = 1000000/(2*(137.036**2))
         allocate(operadores(9,3))
         allocate(spin_rank(9,3))
c         operadores(1,1) = 'OVERLAP '
c         operadores(2,1) = 'OVERLAP '
c         operadores(3,1) = 'OVERLAP '

         operadores(:,1) = 'IDENTIDA'
         spin_rank(:,1) = '1'
c                              nuc coord
         operadores(1,2) = 'SD '//nuc_coord(1)//' x'
         operadores(2,2) = 'SD '//nuc_coord(1)//' y'
         operadores(3,2) = 'SD '//nuc_coord(1)//' z'
         operadores(4,2) = 'SD '//nuc_coord(2)//' x'
         operadores(5,2) = 'SD '//nuc_coord(2)//' y'
         operadores(6,2) = 'SD '//nuc_coord(2)//' z'
         operadores(7,2) = 'SD '//nuc_coord(3)//' x'
         operadores(8,2) = 'SD '//nuc_coord(3)//' y'
         operadores(9,2) = 'SD '//nuc_coord(3)//' z'
         spin_rank(:,2) = '1'

         operadores(1,3) = 'X1SPNORB'
         operadores(2,3) = 'Y1SPNORB'
         operadores(3,3) = 'Z1SPNORB'
         spin_rank(:,3) = '1'

c         call READAOP(aux1, NBAST, 'X1SPNORB') !TEST
c         call READAOP(aux1, NBAST, 'FC X1001') !TEST
c         call READAOP(aux1, NBAST, 'OVERLAP ') !TEST
        END IF


c     Elemento ij del tensor de shielding
      SUMA_B = SIZE(operadores,DIM=1)/3
      do i = 1, 3 !SIZE(operadores,DIM=1)
      do j = 1, 3 !SIZE(operadores,DIM=1)
         QR = 0.d0
         TOTALQR = 0.d0
         operadores_trio(1) = operadores(i,1)
c         operadores_trio (2) = operadores(j,2)
         spin_trio(1) = spin_rank(i,1)
c         spin_trio(2) = spin_rank(j,2)

         IF (operadores(1,3)(3:8).EQ.'SPNORB') THEN
            k = 3
         else
            k = 1
         END IF

        IF (KEY(4:5).EQ.'SD') THEN
           l = 3
        else
           l = 1
        END IF

         do m = 1, k
            operadores_trio(3)=operadores(m,3)
            spin_trio(3)=spin_rank(m,3)
         do n = 1, l
c            operadores_trio(2) = operadores(j+n-1,2)
c            write(lupri,*) "parametro",SUMA_B,"indice: ",SUMA_B*(j-1)+n
            operadores_trio(2) = operadores(SUMA_B*(j-1)+n,2)
            spin_trio(2)=spin_rank(SUMA_B*(j-1)+n,2)
c            write(lupri,*) 'Juan-1-indices ',i,' ',j,' ',operadores(1,1)
c            write(lupri,*) 'JUAN ',operadores_trio(1),' ',
c     &                     operadores_trio(2),' ',operadores_trio(3)
            call create_QRinput(operadores_trio,spin_trio)
            LUPRI = 6232
            OPEN(UNIT=LUPRI, FILE='LRESC-QR',status='unknown')
            call RSPDRV(WORK,LWORK)
            CLOSE(LUPRI)
            LUPRI  = 2
         
            CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
            CALL GPCLOSE(LURSP3,'DELETE')
            CALL GPCLOSE(LURSP4,'DELETE')
            CALL GPCLOSE(LURSP5,'DELETE')
            call get_QRvalue(QR,operadores_trio)
            TOTALQR = TOTALQR + QR
c        Change print level...
            write(lupri,*) operadores_trio(1),' ',operadores_trio(2),!debug
     &         ' ',operadores_trio(3),' ',spin_trio(1),' ',spin_trio(2),
     &         ' ',spin_trio(3)," QR: ",QR
         end do
         end do

         QRtensor(i,j) = TOTALQR
c         write(lupri,*) KEY,"QR-tensor:",constante*QRtensor(i,j)
         write(lupri,*) KEY,"(",i,",",j,"):",QRtensor(i,j)
         write(lupri,*) 
      end do
      end do

      IF (KEY.EQ.'OZ-PSO-MV') QRPSOOZMV = QRtensor
      IF (KEY.EQ.'OZ-PSO-Dw') QRPSOOZDW = QRtensor
      IF (KEY.EQ.'OZ-FC-SO ') QROZFCSO = QRtensor
      IF (KEY.EQ.'OZ-SD-SO ') QROZSDSO = QRtensor


c      do i = 1, SIZE(operadores,DIM=1)
c         call create_QRinput(operadores(i,:),spin_rank(i,:))
c        Imprimo la parte de respuestas fuera del out (sino se llena mucho)
c         LUPRI = 6232
c         OPEN(UNIT=LUPRI, FILE='LRESC-QR',status='unknown')
c         call RSPDRV(WORK,LWORK)
c         CLOSE(LUPRI)
c         LUPRI  = 2
         
c         CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
c         CALL GPCLOSE(LURSP3,'DELETE')
c         CALL GPCLOSE(LURSP4,'DELETE')
c         CALL GPCLOSE(LURSP5,'DELETE')
c      end do

c      do i = 1, SIZE(operadores,DIM=1)
c         call get_QRvalue(QR,operadores(i,:))
c         TOTALQR = TOTALQR + QR
c      end do

c      write(lupri,*) KEY," TOTAL-QR: ",TOTALQR*constante

      deallocate(operadores)
      deallocate(spin_rank)




C ACA ARRANCA PARTE BIEN VIEJA
cc      LUESC = -1
c      call copy_file(LUESC,LUESC2,'d22_TEST_QR1.dal', 'DALTON.INP')
c      LUPRI = 6232
c      OPEN(UNIT=LUPRI, FILE='trash.txt',status='unknown')
c      call RSPDRV(WORK,LWORK)
c      CLOSE(LUPRI)
c      operadores(1) = 'XANGMOM '
c      operadores(2) = 'FC X1001'
c      operadores(3) = 'X1SPNORB'
c      LUPRI  = 2
c      call get_QRvalue(QUAVAL,operadores)

c      write(lupri,*) "AUCAR-LRESC QR VALUE:",QUAVAL

c      CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
c      CALL GPCLOSE(LURSP3,'DELETE')
c      CALL GPCLOSE(LURSP4,'DELETE')
c      CALL GPCLOSE(LURSP5,'DELETE')

c      call copy_file(LUESC,LUESC2,'d22_TEST_QR2.dal', 'DALTON.INP')
c      LUPRI = 6232
c      OPEN(UNIT=LUPRI, FILE='trash.txt',status='unknown')
c      call RSPDRV(WORK,LWORK)
c      CLOSE(6232)
c      operadores(1) = 'XANGMOM '
c      operadores(2) = 'DARWIN  '
c      operadores(3) = 'MASSVELO'
c      LUPRI  = 2
c      call get_QRvalue(QUAVAL,operadores)
c      write(lupri,*) "AUCAR-LRESC QR VALUE:",QUAVAL
c      CALL GPCLOSE(LURSP ,'KEEP') !RSPVEC file
c      CALL GPCLOSE(LURSP3,'DELETE')
c      CALL GPCLOSE(LURSP4,'DELETE')
c      CALL GPCLOSE(LURSP5,'DELETE')


      CALL QEXIT('LRSCQUAD')
      CALL TIMER ('LRSCQUAD',TIMEIN,TIMOUT) !Debiera agregar un if segun print level
      RETURN

      END 


      subroutine create_QRinput(operadores,spin_rank)
      character*10 :: namefile
      logical            :: file_exists, file_open
      character(100) :: line
      INTEGER :: fileID
      CHARACTER*8, intent(in) :: operadores(3)
      CHARACTER*1, intent(in) :: spin_rank(3)

#include "priunit.h"
      namefile='DALTON.INP'
      fileID = -1
      CALL GPOPEN(fileID,namefile,'NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

c      line = 'whatever I want to write'
c      write(fileID,'(a)') line
c      write(fileID,'(a)') 'easier way'
      write(fileID,'(a)') '**DALTON INPUT'
      write(fileID,'(a)') '**RESPONSE'
      write(fileID,'(a)') '*LINEAR'
      write(fileID,'(a)') '.PRINT'
      write(fileID,'(a)') ' 0'
      write(fileID,'(a)') '*QUADRA'
      write(fileID,'(a)') '.PRINT'
      write(fileID,'(a)') ' 0'
      write(fileID,'(a)') '.ISPABC'
      write(fileID,'(a)') spin_rank(1)//'    '//
     &                    spin_rank(2)//'    '//spin_rank(3)
      write(fileID,'(a)') '.APROP'
      write(fileID,'(a)') operadores(1)//' '//spin_rank(1)
      write(fileID,'(a)') '.BPROP'
      write(fileID,'(a)') operadores(2)//' '//spin_rank(2)
      write(fileID,'(a)') '.CPROP'
      write(fileID,'(a)') operadores(3)//' '//spin_rank(3)
      write(fileID,'(a)') '**END OF DALTON INPUT'

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP')

c      CLOSE(fileID)


      RETURN
      END


      subroutine copy_file(oldID,newID,old_file, new_file)
      character(*) :: old_file ,new_file
      logical            :: file_exists, file_open
      character(100) :: line
      INTEGER :: oldID,newID

#include "priunit.h"

c      OPEN(UNIT=oldID, FILE=old_file, STATUS='OLD'
c     &     ,ACTION='READ', IOSTAT=IERR, RECL=1)
c      REWIND(oldID)
c      OPEN(UNIT=newID, FILE=new_file,STATUS='REPLACE',
c     &     ACTION='WRITE', IOSTAT=IERR, RECL=1, FORM='FORMATTED')

c      newID=8756323
c      OPEN(UNIT=newID, FILE=new_file,status='unknown')
       CALL GPOPEN(newID,new_file,'OLD',' ','FORMATTED',
     &             IDUMMY,.FALSE.)

      inquire(file = old_file, exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            oldID = -1
            CALL GPOPEN(oldID,old_file,'OLD',' ','FORMATTED',
     &           IDUMMY,.FALSE.)
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
c      if (file_open) CALL GPCLOSE(oldID,'KEEP')
c      if (file_open) CALL GPCLOSE(newID,'KEEP')


c      CLOSE(newID)

      RETURN
      END

      subroutine get_QRvalue(QUAVAL,labels_in)
      logical :: file_exists, file_open
      INTEGER :: fileID
      CHARACTER STRING*75
      CHARACTER*8, intent(in) :: labels_in(3)
      CHARACTER*8:: labels_read(3)
c      CHARACTER*8, intent(in) :: PRUEBA
      REAL*8 FRQB, FRQC
      REAL*8, intent(out) :: QUAVAL
c

#include "priunit.h"

C Based on some lines of SUBROUTINE BCCHK

      inquire(file = 'RESULTS.RSP', exist = file_exists,
     &            opened = file_open)
      if (file_exists) then
         if (.not. file_open) then
            fileID = -1
            CALL GPOPEN(fileID,'RESULTS.RSP','OLD',' ','FORMATTED',
     &           IDUMMY,.FALSE.)
         end if
         rewind(fileID)
      else
         call quit('RESULTS.RST file not found')
      end if


 346  READ (fileID,'(A75)',END=347) STRING
      IF (STRING(1:8) .EQ. '@ Quadra') THEN
c De referencia         
c               READ (LURSPRES,'(30X,A8,2I5)') LABEL(1),KSYMA,KSPINA
c               READ (LURSPRES,'(30X,A8,2I5)') LABEL(2),KSYMB,KSPINB
c               READ (LURSPRES,'(30X,A8,2I5)') LABEL(3),KSYMC,KSPINC

c Primer intento (hay que ir borrando el archivo de lectura)
c         READ (fileID,'(A75)') STRING
c         READ (fileID,'(A75)') STRING
c         READ (fileID,'(A75)') STRING

          READ (fileID,'(30X,A8)') labels_read(1)
          READ (fileID,'(30X,A8)') labels_read(2)
          READ (fileID,'(30X,A8)') labels_read(3)

      IF (labels_read(1) .EQ. labels_in(1) .AND.
     &    labels_read(2) .EQ. labels_in(2) .AND.
     &    labels_read(3) .EQ. labels_in(3) ) then
c          write(2,*) "TEST-AUCAR-read",labels_read
c          write(2,*) "TEST-AUCAR-in",labels_in
         READ (fileID,'(/,30X,3F15.8)') FRQB, FRQC, QUAVAL
      ELSE
         GOTO 346
      END IF

      ELSE
         GOTO 346
      END IF

 347         CONTINUE
C         QUAVAL = -QUAVAL !ESTO NO VA... Lo había agregado pero está  mal


cdo while (.true.)
c         read(oldID,'(a)', end=31) line
c         write(newID,'(a)') line
c      end do
c 31   continue

      if (.not. file_open) CALL GPCLOSE(fileID,'KEEP') !leer trayendo labels
c      if (.not. file_open) CALL GPCLOSE(fileID,'REMOVE') !elimino RESULTS.RSP luego de leer


      RETURN
      END
