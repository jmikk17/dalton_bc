      subroutine numint_driver(ao_density,work,lwork) !,iorder,propri,iprint)
!*****************************************************************************
!     JESM 16-10-2024
!     Driver routine for numerical integration of atomic moments
!     
!     Purpose:
!     Calculate AO integrals for all atomic moments specificed in input and
!     write to AOPROPER
!*****************************************************************************
#include "implicit.h"
#include "dftinf.h"
#include "dftcom.h"
#include "inforb.h"
#include "inftap.h"
#include "dummy.h"
#include "rspprp.h"
!
! dftinf.h:       nder, ntypso 
! dftcom.h:       dftgrid_done_old, dfthri
! inforb.h:       nbast, nnashx, n2basx, nisht, nasht, ncmot (TODO: ccorb.h is a seperate file, is that one neccesary?)
! inftap.h:       lusifc, lbsifc
! dummy.h:        idummy
!
      integer luddump, luquad
      dimension work(lwork), ao_density(n2basx)

      character(len=8) label

      call qenter('numint')
      call newtimer('START')

      kfrsav = 1
      kfree  = kfrsav
      lfree  = lwork 

! open file for dumping density and basis functions in each grid point
      luddump = -1
      call gpopen(luddump,'DENSITY.GRID','UNKNOWN',
     &            ' ','UNFORMATTED',idummy,.false.)
! file for dumping population weights
      luwdump = -1
      call gpopen(luwdump,'WEIGHT.GRID','UNKNOWN',
     &            ' ','UNFORMATTED',idummy,.false.)
! file for ao property integrals
      luprop = -1
      call gpopen(luprop,'AOPROPER','UNKNOWN',
     &            ' ','UNFORMATTED',idummy,.false.)

! set up dft orbital common block so we don't evaluate uneccesary terms 
! TODO: maybe these should be saved and set back to original val?
! TODO: look into calling DFTINI or SETUPSOS instead of this
      DFTHR0 = 1.0D-8
      DFTHRL = 1.0D-10
!     DFTHRI = 1.0D-13
      DFTHRI = 2.0D-12
      nder = 0
      ntypso = 1
      nso0 = 1
      kso0 = 1

! make dftgrid BEFORE allocations, this requires a lot of mem but releases it again
      if (.not. dftgrid_done_old) then
         call make_dftgrid(work,lwork,ntot_dftgrid,1,.false.)
         dftgrid_done_old = .true.
      end if

! allocate static memory for 1 000 000 grid points
! this seems like its nbuf pr atom, so can probably be lowered?
! but is hardcoded to 1 000 000 in grid routine
      nbuf = 1 000 000

      write(*,*) 'lfree before allocation', lfree

! grid coords
      call memget('REAL',kx,nbuf,work,kfree,lfree)
      call memget('REAL',ky,nbuf,work,kfree,lfree)
      call memget('REAL',kz,nbuf,work,kfree,lfree)
! integrations weights
      call memget('REAL',kw,nbuf,work,kfree,lfree)
! basis functions evaluated in a point 
      call memget('REAL',kgso,nbast,work,kfree,lfree)
! for ncnt, goes into getsos, propbalby related to which center the basis func. is on 
      call memget('INT4',kcnt,nbast,work,kfree,lfree) 
! DMATGAO, unneccesary output from getrho_srdft
      call memget('REAL',kdga,nbast,work,kfree,lfree) 

      write(*,*) 'lfree after allocation', lfree
      
! print evaluated basis functions and density in each point to file
      call evaluate_grid(work(kx),work(ky),work(kz),work(kw),
     &   work(kgso),work(kcnt),work(kdga),ao_density,
     &   luddump,nbuf,work,lwork)

      call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)

! allocate for mbis routine 
      max_pnt = 1 000 000 
      max_atm = 50 !max atoms
         
      write(*,*) 'lfree before mbis', lfree

      call memget('REAL',kx,max_pnt,work,kfree,lfree)
      call memget('REAL',ky,max_pnt,work,kfree,lfree)
      call memget('REAL',kz,max_pnt,work,kfree,lfree)
      call memget('REAL',kw,max_pnt,work,kfree,lfree)

      write(*,*) 'lfree mid mbis', lfree

      call memget('REAL',krho,max_pnt,work,kfree,lfree)
      call memget('REAL',kw_mbis,max_atm,work,kfree,lfree)

      write(*,*) 'lfree after mbis', lfree

! calculate mbis population weights in each point and print to file
      call mbis_pop(work(kx),work(ky),work(kz),work(kw),
     &   work(krho),work(kw_mbis),
     &   max_pnt,max_atm,luddump,luwdump)

      call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)

      
! allocate for ao integrals

      call memget('REAL',kint,nnbasx,work,kfree,lfree)
      call memget('REAL',kgso,nbast,work,kfree,lfree)

! calculate ao-integrals and print to AOPROPER
! TODO: move this label loop inside the numerical integration
      do i = 1,nprlbl
         label = prplbl(i)
         if (label(1:2) .eq. 'AM') then
            call num_aoint(label,work(kint),work(kgso),
     &         luddump,luwdump)
         end if
      enddo

      call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)

      call gpclose(luprop,'KEEP')
      call gpclose(luddump,'KEEP')
      call gpclose(luwdump,'KEEP')

      call newtimer('numint')
      call qexit('numint')

      return
      end


      subroutine evaluate_grid(corx,cory,corz,int_weight,
     &   gso,ncnt,dmatgao,ao_density,luddump,nbuf,
     &   work,lwork)
#include "implicit.h"
#include "inforb.h"
#include "dummy.h"
#include "rspprp.h"
#include "inflr.h"
#include "dftcom.h"
#include "priunit.h"

      parameter (d0 = 0.0d0)

! ingoing
      real*8    corx(nbuf),cory(nbuf),corz(nbuf),int_weight(nbuf),
     &          gso(nbast),ncnt(nbast),dmatgao(nbast),
     &          ao_density(nbast,nbast),work(lwork)
      integer   luddump,nbuf

! defined within routine
      real*8    cor(3)
      integer   npnts,npoint,ipnt
      real*8    nelctrn,weight,rho

! open file where dft grid is stored 
      luquad = -1
      call gpopen(luquad,'DALTON.QUAD','OLD','SEQUENTIAL',
     &   'UNFORMATTED',idummy,.false.)

      write(*,*) 'ao dens in eval grid'
      write(*,*) ao_density

      nelctrn = d0

      rewind luddump

      npnts = 0
  100 continue
      read(luquad) npoint
      if (npoint .gt. 0) then
         write(luddump) npoint
         npnts = npnts + npoint
         call reaqua_srdft(corx,cory,corz,int_weight,luquad,npoint)
         do 200 ipnt = 1, npoint

            if (iprlr .gt. 150) then
               write (lupri,'(I8,4F12.6)')
     &         ipnt,corx(ipnt),cory(ipnt),corz(ipnt),int_weight(ipnt)
            END IF

            weight = int_weight(ipnt)
            cor(1) = corx(ipnt)
            cor(2) = cory(ipnt)
            cor(3) = corz(ipnt)
            ! dfthri hardcoded to 10^-13 in driver
            ! TODO: add a check on weight being 0
            thrint = dfthri/weight 

! evaluate orbitals in the grid point
            call getsos(gso,ncnt,cor,work,lwork,
     &                  nbast,.false.,.false.,thrint,iprlr)
!           CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,
!    &                  NBAST,DO_LONDON,DOGGA,THRINT,IPRINT)

! get density in point
            rho = d0
            call getrho_srdft(ao_density,gso,rho,dmatgao,
     &         thrint,iprlr)
            nelctrn = nelctrn + weight*rho

            write(luddump) corx(ipnt),cory(ipnt),corz(ipnt),
     &           int_weight(ipnt),rho,(gso(j),j=1,nbast)

  200    continue
         goto 100
      else if (npoint .eq. 0) then
         goto 100
      end if
      write(luddump) npoint ! writing the negative npoint val as a footer?

      call gpclose(luquad,'KEEP')

! testing we integrate correct num. of electrons
      elctrx = float(2*nrhft)
      error  = nelctrn - elctrx

! print section for num int
      write (lupri,'(1(/2X,A,F14.6))')
     &' Nr. of electrons from numerical integration          :',nelctrn
      write (lupri,'(2(/2X,A,F14.6))')
     &' Number of electrons from orbial occupations          :',elctrx,
     &' Error in the number of electrons                     :',error
      if (abs(error) .gt. dftels) then
         write (lupri,'(/2X,A,F14.6,/2X,A)')
     &' Error larger than DFTELS (set input)                 :',dftels,
     &   ' Calculation aborted.'
         call quit
     &    ('Wrong number of electrons in numint_driver.')
      end if

      return
      end

! Based on the MBIS article by Verstraelen et.al - doi.org/10.1021/acs.jctc.6b00456
      subroutine mbis_pop(corx,cory,corz,int_weight,rho,
     &                 mbis_weight,max_pnt,max_atm,luddump,luwdump)
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "priunit.h"
#include "infrsp.h"
#include "inforb.h"

      parameter (d0 = 0.0d0, pi=acos(-1.0d0), epi=8.0d0*pi,
     &           qconv=1.0d-9, eps=1.0d-14, dencut=1.0d-10)
! incoming
      real*8    corx(max_pnt),cory(max_pnt),corz(max_pnt),
     &          int_weight(max_pnt),rho(max_pnt),
     &          mbis_weight(max_atm)
      integer   max_pnt,max_atm,luddump
! defined in here
! TODO: some should be allocated with memget instead
! 10 used as max l, but currently only imp. for 3 in zeffini
      integer   npoint,lmax(max_atm),ipt
      real*8    zan(max_atm),qatm(max_atm),qtmp(max_atm),
     &          zeff(10,max_atm),zstmp(10,max_atm),
     &          qshell(10,max_atm),qstmp(10,max_atm),
     &          wtatm(10,max_atm),
     &          suma,tote,qtot,dx,dy,dz,rhom,fwt,tmp,ztmp,xtmp,rmsq,dis,
     &          coor(3,max_atm), wtot

      write(*,*) 'check from mbis_pop', qconv

      kmax = 500
      kpar = max_atm
      kpt = max_pnt

      nat=nucdep
      if (nat.gt.kpar) then
        write(lupri,*)'MBIS: too many atoms',nat,kpar
        call quit('MBIS: too many atoms')
      endif

! saving variables in other names to fit with old code
! TODO: change zan to atmchg from infrsp.h. Atmchg includes degenerate atoms
      do i=1,nat
        zan(i)=charge(i) !nuclear charge of atom, nuclei.h
        qatm(i)=d0
        do j=1,3
          coor(j,i)=cord(j,i) !cord(3,mxcent), nuclei.h
        end do
      enddo

      rewind(luddump)
      rewind(luwdump)

      suma = d0
      ipt=1
      start=1
10    continue
      read(luddump) npoint
      write(*,*) 'npoint', npoint
      if (npoint.gt.0) then
         do 200 i=start,start+npoint-1
            read(luddump) corx(ipt),cory(ipt),corz(ipt),int_weight(ipt),
     &                    rho(ipt)
            suma=suma+int_weight(ipt)*rho(ipt)
            ipt=ipt+1
            if (ipt.eq.kpt) then
               write(lupri,*) 'MBIS: too many points',ipt,kpt
               call quit('MBIS: too many points')
            endif
 200     continue
         write(*,*) 'suma', suma
         start=start+npoint
         goto 10
      else if (npoint.eq.0) then
         goto 10
      end if

! tot grid points = ipt-1 due to the way the reading is constructed
      ipt=ipt-1
      write(lupri,*) ' grid points, tot elec.'
      write(lupri,14) ipt,suma
      write(lupri,*)' '
14    format(i10,2f20.10)
30    format(2i5,5f12.4)
      
! guess initial parameters for pro-atomic density 
      call zeffini(kpar,nat,zan,zeff,qshell,lmax)
      write(*,*)'Zeff initial' ! 1/sigma
      do i=1,nat
        write(*,30)i,lmax(i),(zeff(j,i),j=1,3)
      enddo
      write(*,*)'Qshell initial' ! N_Ai 
      do i=1,nat
        write(*,30)i,lmax(i),(qshell(j,i),j=1,3)
      enddo

! iterate until MBIS charges are converged
      do kiter=1,kmax
         write(*,*)' '
         write(*,*)' Iteration ',kiter

! in the first part of the itteration we calculate the population of the pro-atomic densities
         do n=1,nat
            do l=1,lmax(n)
               qstmp(l,n)= d0
            enddo
         enddo

c       loop over points
         do i=1,ipt
            rhom = rho(i)
            fwt = int_weight(i) 

            do n=1,nat
               do l=1,lmax(n)
                  wtatm(l,n)=0.0d0
               enddo
            enddo
c         construct atomic weight function in the given point
            do n=1,nat
               dx = corx(i) - coor(1,n)
               dy = cory(i) - coor(2,n)
               dz = corz(i) - coor(3,n) 
               dis = dx*dx + dy*dy + dz*dz
               dis = sqrt(dis)
               do l=1,lmax(n)
                  ztmp = zeff(l,n)
                  znorm = (ztmp**3)/epi
                  xtmp = znorm*exp(-ztmp*dis)
                  rhoa = qshell(l,n)
                  tmp = xtmp*rhoa*fwt
                  if (abs(tmp).lt.dencut) tmp=d0
                  wtatm(l,n) = wtatm(l,n) + tmp
               enddo
            enddo

!         assign ISA fractions and partition electrons
            wtot = d0
            do n=1,nat
               do l=1,lmax(n)
                  wtot = wtot + wtatm(l,n) 
               enddo
            enddo
            if (abs(wtot).gt.d0 .and. (abs(rhom*fwt).gt.eps)) then
               do n=1,nat
                  do l=1,lmax(n)
                     qstmp(l,n) = qstmp(l,n) + wtatm(l,n)*fwt*rhom/wtot !N(Ai) population, eq 18 
                  enddo
               enddo
            endif
         enddo
!        end ipt loop

         write(*,*)' Qshell'
         tote=d0
         do n=1,nat
            do l=1,lmax(n)
               tote = tote + qstmp(l,n)
               tmp = qstmp(l,n)-qshell(l,n)
               write(*,21)n,l,qshell(l,n),qstmp(l,n),tmp 
            enddo
         enddo
         write(*,*)' Electrons from integration:', suma
         write(*,*)' Electrons from MBIS       :', tote
         rmsq=d0
         qtot=d0
         qele=d0
         write(*,*)' Qatom'
         do n=1,nat
            qtmp(n)=d0
            do l=1,lmax(n)
               qtmp(n)=qtmp(n)+qstmp(l,n)
            enddo
            qele = qele + qtmp(n) 
            qtmp(n) = zan(n) - qtmp(n)
            tmp = qtmp(n)-qatm(n) 
            rmsq = rmsq + tmp*tmp 
            qtot = qtot + qtmp(n)
            write(*,20)n,qatm(n),qtmp(n),tmp
         enddo
         qxs=nint(qele)
         qxs=qele+qtot-qxs
         rmsq=sqrt(rmsq)
         write(*,22)rmsq,qele,qtot,qxs
         if (rmsq.lt.qconv) goto 900
         do n=1,nat
            qatm(n) = qtmp(n) !N(A) saved for next iteration (for eq 20)
            do l=1,lmax(n)
               qshell(l,n)=qstmp(l,n) !New N(Ai) saved for next iteration to make p0(Ai)
            enddo
         enddo

c        update zeff
         do n=1,nat
            do l=1,lmax(n)
               zstmp(l,n)=d0
            enddo
         enddo

c        loop over points
         do i=1,ipt
            rhom = rho(i)
            fwt = int_weight(i)

            do n=1,nat
               do l=1,lmax(n)
                  wtatm(l,n)=d0
               enddo
            enddo

c         reconstruct atomic weights functions - same as before
            do n=1,nat
               dx = corx(i) - coor(1,n)
               dy = cory(i) - coor(2,n)
               dz = corz(i) - coor(3,n) 
               dis = dx*dx + dy*dy + dz*dz
               dis = sqrt(dis)
               do l=1,lmax(n)
                  ztmp = zeff(l,n)
                  znorm = (ztmp**3)/epi
                  xtmp = znorm*exp(-ztmp*dis)
                  rhoa = qshell(l,n)
                  tmp = xtmp*rhoa*fwt
                  if ((abs(tmp)).lt.dencut) tmp=d0
                  wtatm(l,n) = wtatm(l,n) + tmp
               enddo
            enddo
c         assign ISA fractions and partition electrons 
            wtot = d0
            do n=1,nat
               do l=1,lmax(n)
                  wtot = wtot + wtatm(l,n) !does the abs change these?
               enddo
            enddo
            if (abs(wtot).gt.d0 .and. abs(rhom*fwt).gt.eps) then
               do n=1,nat
                  dx = corx(i) - coor(1,n)
                  dy = cory(i) - coor(2,n)
                  dz = corz(i) - coor(3,n) 
                  dis = dx*dx + dy*dy + dz*dz
                  dis = sqrt(dis)
                  vtmp = zan(n)/dis
                  do l=1,lmax(n)
                     zstmp(l,n) = zstmp(l,n) + 
     &                  wtatm(l,n)*rhom*dis*fwt/wtot 
                  enddo
               enddo
            endif
         enddo
c        end ipt loop

         write(*,*)' Zshell'
!        write(lupri,*)' Zshell'
         do n=1,nat
            do l=1,lmax(n)
               zstmp(l,n) = (3.0d0*qshell(l,n))/zstmp(l,n)
               tmp = zstmp(l,n)-zeff(l,n)
               write(*,21)n,l,zeff(l,n),zstmp(l,n),tmp
!              write(lupri,21)n,l,zeff(l,n),zstmp(l,n),tmp
               zeff(l,n) = zstmp(l,n)
            enddo
         enddo

        write(*,*)' '
      enddo 
c     end kitter loop

20    format(i5,3f15.8)
21    format(2i5,3f15.8)
22    format(' RMS-dQ , Qele, Qtot, Qxs',4f15.8)

900   continue

      write(*,*) ' MBIS converged!'

      write(*,*)' '
      write(*,*)' Final converged results'
      write(*,*)' Qatom'
      do n=1,nat
        qtmp(n)=d0
        do l=1,lmax(n)
          qtmp(n)=qtmp(n)+qstmp(l,n)
        enddo
        qtmp(n) = nint(zan(n)) - qtmp(n)
        write(*,20)n,qtmp(n)
      enddo

! calculate final mbis weights for print in each point
! using qstmp since qshell is not saved in final iteration

      write(luwdump) ipt

      do i=1,ipt
         rhom = rho(i)
         fwt  = int_weight(i) 
         do n=1,nat
            mbis_weight(n) = d0
            do l=1,lmax(n)
               wtatm(l,n) = d0
            enddo
         enddo
         do n=1,nat
            dx = corx(i) - coor(1,n)
            dy = cory(i) - coor(2,n)
            dz = corz(i) - coor(3,n) 
            dis = dx*dx + dy*dy + dz*dz
            dis = sqrt(dis)
            do l=1,lmax(n)
               ztmp = zeff(l,n)
               znorm = (ztmp**3)/epi
               xtmp = znorm*exp(-ztmp*dis)
               rhoa = qstmp(l,n)
               tmp = xtmp*rhoa*fwt
               if (abs(tmp).lt.dencut) then
                  tmp=d0
               end if
               wtatm(l,n) = wtatm(l,n) + tmp
            enddo
         enddo

         wtot = d0
         do n=1,nat
            do l=1,lmax(n)
               wtot = wtot + wtatm(l,n)
            enddo
         enddo
         if (abs(wtot).gt.d0 .and. (abs(rhom*fwt).gt.eps)) then
            do n=1,nat
               do l=1,lmax(n)
                  mbis_weight(n) = mbis_weight(n) + wtatm(l,n)/wtot
               enddo
            enddo
         endif
         write(luwdump) (mbis_weight(j),j=1,nat)
      enddo !ipt

      write(luwdump) -1

      return
      end

      subroutine zeffini(kpar,nat,zan,zeff,qshell,lmax)
!     Purpose : Obtain the initial values of zeff and qshell 
#include "implicit.h"
#include "mxcent.h"
      parameter(d0=0.0d0)
#include "nuclei.h"
      real*8    zan(kpar),zeff(10,kpar),qshell(10,kpar),znuc
      integer   lmax(kpar)

      zeff = 0.0d0
      qshell = 0.0d0
      lmax = 0

      do i=1,nat
         znuc = zan(i)
         zeff(1,i) = 2.0d0*znuc
         if (znuc.le.2.0d0) then
            lmax(i) = 1
            qshell(1,i)=znuc
            qshell(2,i)=d0
         endif
         if (znuc.gt.2.0d0 .and. znuc.le.10.0d0) then
            lmax(i) = 2
            zeff(2,i) = 2.0d0
            qshell(1,i)=2.0d0
            qshell(2,i)=znuc-2.0d0 
         endif
         if (znuc.gt.10.0d0 .and. znuc.le.18.0d0) then
            lmax(i) = 3
            zeff(2,i) = 2.0d0*sqrt(znuc)
            zeff(3,i) = 2.0d0
            qshell(1,i)=2.0d0
            qshell(2,i)=8.0d0 
            qshell(3,i)=znuc-10.0d0 !Z-core=valens, neccesary since we dont know if valence is full
         endif
         if (znuc.gt.18) call quit('Z > 18 not implemented in zeffini')
      enddo
 
      return
      end


      subroutine num_aoint(integral_label,ao_int,gso,luddump,luwdump)
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "priunit.h"
#include "rspprp.h"
#include "inflr.h"
      parameter (d0 = 0.0d0)
!incoming
      character*8 integral_label
      integer     luddump,luwdump
      real*8      ao_int(nnbasx),gso(nbast)
!in here
      character*8 integral_type, rtnlbl(2)
      integer     ipt, ipt_check, atm_idx, atm_chg, cor_idx
      real*8      tmp, int_weight, mbis_weight, cor(3),
     &            ao_test(nbast,nbast)

! we could in principle calc. all integrals the first time we are here, but
! with the 'select' feature we will not know which integrals we need

! label: AM + two char for index + two char for charge + two char for direc. and moment

      read(integral_label(3:4),*) atm_idx
      read(integral_label(5:6),*) atm_chg

      cor_idx = 0

      integral_type = integral_label(7:8)
      select case (integral_type)
! this can be extended to two cor_idx for quadrupole integrals
         case ('00')
            cor_idx = 0
         case ('0x')
            cor_idx = 1
         case ('0y')
            cor_idx = 2
         case ('0z')
            cor_idx = 3
         case default
            write(*,*) 'Label', integral_type, 'not reognized'
     &         // 'or implemented yet'
            call quit("Label not recognized in numint_write")
         end select

      rewind(luddump)
      rewind(luwdump)

      read(luwdump) ipt_check

      ao_int = d0

      start=1
      ipt=0

      if (cor_idx.eq.0) then
      
!luddump: corx(r),cory(r),corz(r),int_weight(r),rho(r),(chi_a(r),a=1,nbast)
!luwdump: (mbis_weight_a(r),a=1,nucdep)

10       continue
         read(luddump) npoint
         if (npoint.gt.0) then
            do i=start,start+npoint-1
               ipt=ipt+1
               read(luddump) tmp,tmp,tmp,int_weight,tmp,
     &                        (gso(j),j=1,nbast)
               read(luwdump) (tmp,j=1,atm_idx-1), mbis_weight
               do a=1,nbast
                  ! index of the last element in previous row, i.e. the diagonal
                  idx=(a*a-a)/2 
                  do b=1,a
                     ao_int(idx+b) = ao_int(idx+b) + 
     &                  gso(a)*gso(b)*mbis_weight*int_weight
                  enddo
               enddo
            enddo
            start=start+npoint
            goto 10
         else if (npoint.eq.0) then
            goto 10
         end if

      else if (cor_idx.gt.0) then

20       continue
         read(luddump) npoint
         if (npoint.gt.0) then
            do i=start,start+npoint-1
               ipt=ipt+1
               read(luddump) cor(1),cor(2),cor(3),int_weight,tmp,
     &                     (gso(j),j=1,nbast)
               read(luwdump) (tmp,j=1,atm_idx-1), mbis_weight
               dist = cor(cor_idx) - cord(cor_idx,atm_idx)
               do a=1,nbast
                  idx=(a*a-a)/2
                  do b=1,a
                     ao_int(idx+b) = ao_int(idx+b) + 
     &                  gso(a)*gso(b)*mbis_weight*int_weight*dist
                  enddo
               enddo
            enddo
            start=start+npoint
            goto 20
         else if (npoint.eq.0) then
            goto 20
         end if

      end if 

      if (ipt .ne. ipt_check) 
     &   call quit("Error in reading of files in numint")

      read(luwdump) ipt_check
      if (ipt_check .ne. -1) 
     &   call quit("Error in reading of files in numint")
      
!     CALL DSPTSI(NBAST,ao_int,ao_test)
!     write(lupri,*) 'ao_int test'
!     CALL OUTPUT(ao_test,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)

      call getdat(rtnlbl(1),rtnlbl(2))

      rtnlbl(2) = 'SYMMETRI'
      !rtnlbl(1) holds opperator sym, we ignore that for now?
         !this is not read in in getgpv, seems like its inferred from integral vals.?
      
      call wrtpro(ao_int,nnbasx,integral_label,rtnlbl,iprlr)

      return
      end
