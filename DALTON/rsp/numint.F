      subroutine numint_driver(ao_density,work,lwork,cc_flag)
!*****************************************************************************
!     JESM 16-10-2024
!     Driver routine for numerical integration of atomic moments
!     
!     Purpose:
!     Calculate AO integrals for all atomic moments specificed in input and
!     write to AOPROPER
!*****************************************************************************
#include "implicit.h"
#include "dftcom.h"
#include "inforb.h"
#include "dummy.h"
#include "mxcent.h"
#include "nuclei.h"
#include "rspprp.h"
#include "ccrspprp.h"
#include "inftap.h"
#include "inflr.h"
#include "priunit.h"
! Used from include:
! implicit.h:     bad coding practice
! dftcom.h:       dftgrid_done_old
! inforb.h:       nbast, nnbasx, n2basx, nisht, nasht (ccorb.h for cc?)
! dummy.h:        idummy
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       nucdep
! rspprp.h:       prplbl/prplbl_cc
! inftap.h:       luprop
! inflr.h:        iprlr
! priunit.h:      lupri
!
! Files used for:
! luquad:         grid points of the quadrature
! luddump:        dumping density and orbital values in grid points
! luwdump:        dumping mbis weights in grid points
! luprop:         property integrals
! lupri:          output file
!

! incoming:
      integer lwork
      logical cc_flag
      dimension work(lwork), ao_density(n2basx)
! defined here:
      integer luddump, luwdump, kfrsav, kfree, lfree, nbuf
      logical save_mem
      character(len=8) label

! option for old or new num_aoint
      save_mem = .false.

! timers and memory pointers
      call qenter('numint_driver')
      call newtimer('START')

      kfrsav = 1
      kfree  = kfrsav
      lfree  = lwork 

! open files
      luddump = -1
      call gpopen(luddump,'DENSITY.GRID','UNKNOWN',
     &            ' ','UNFORMATTED',idummy,.false.)
      luwdump = -1
      call gpopen(luwdump,'WEIGHT.GRID','UNKNOWN',
     &            ' ','UNFORMATTED',idummy,.false.)
      luprop = -1
      call gpopen(luprop,'AOPROPER','UNKNOWN',
     &            ' ','UNFORMATTED',idummy,.false.)

! set up dft orbital common block 
      call setupsos(0,.false.,idum1,idum2)

! make dftgrid if first run
      if (.not. dftgrid_done_old) then
         call make_dftgrid(work,lwork,ntot_dftgrid,iprlr,.false.)
         dftgrid_done_old = .true.
      end if

! allocate for proccesing of grid data
! make_dftgrid hardcodes limit of points per atom to 1.0d6
! however, can never be more than total grid
      nbuf = min(1 000 000,ntot_dftgrid)

      call memget('REAL',kx,nbuf,work,kfree,lfree)
      call memget('REAL',ky,nbuf,work,kfree,lfree)
      call memget('REAL',kz,nbuf,work,kfree,lfree)
      call memget('REAL',kw,nbuf,work,kfree,lfree)

      call memget('REAL',kgso,nbast,work,kfree,lfree)
      call memget('REAL',kgso_prod,nnbasx,work,kfree,lfree)
      call memget('INT4',kcnt,nbast,work,kfree,lfree) 
      call memget('REAL',kdga,nbast,work,kfree,lfree) 

! evaluate orbitals and density in grid points 
      call evaluate_grid(work(kx),work(ky),work(kz),work(kw),
     &   work(kgso),work(kgso_prod),work(kcnt),work(kdga),ao_density,
     &   luddump,nbuf,work,lwork,cc_flag,save_mem,ntot_screened)

      call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)

! allocate for MBIS populations routine 
      ! maybe move these to a common block? or remove all together?
      max_pnt = min(1 000 000,ntot_screened)
      max_atm = min(50,nucdep)
         
      write(*,*) 'lfree before mbis', lfree

      call memget('REAL',kx,max_pnt,work,kfree,lfree)
      call memget('REAL',ky,max_pnt,work,kfree,lfree)
      call memget('REAL',kz,max_pnt,work,kfree,lfree)
      call memget('REAL',kw,max_pnt,work,kfree,lfree)

      write(*,*) 'lfree mid mbis', lfree

      call memget('REAL',krho,max_pnt,work,kfree,lfree)
      call memget('REAL',kw_mbis,max_atm,work,kfree,lfree)

      write(*,*) 'lfree after mbis', lfree

! evaluate MBIS pop. weights in each point
      call mbis_pop(work(kx),work(ky),work(kz),work(kw),
     &   work(krho),work(kw_mbis),
     &   max_pnt,max_atm,luddump,luwdump)

      call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)
! allocate for and calculate ao integrals, memory saving version
      if (save_mem) then

         call memget('REAL',kint,nnbasx,work,kfree,lfree)
         call memget('REAL',kgso,nbast,work,kfree,lfree)

         if (cc_flag) then
            do i = 1,nprlbl_cc
               label = prplbl_cc(i)
               write(*,*) 'lbl check', label
               if (label(1:2) .eq. 'AM') then
                  call num_aoint(label,work(kint),work(kgso),
     &               luddump,luwdump)
               end if
            enddo
         else
            do i = 1,nprlbl
               label = prplbl(i)
               write(*,*) 'lbl check', label
               if (label(1:2) .eq. 'AM') then
                  call num_aoint(label,work(kint),work(kgso),
     &               luddump,luwdump)
               end if
            enddo
         end if

! same as above, but hopefully faster 
      else
         nbuf2 = 30 000
         natmlbl = 0
         if (cc_flag) then
            do i = 1,nprlbl_cc
               label = prplbl_cc(i)
               if (label(1:2) .eq. 'AM') then
                  natmlbl = natmlbl + 1
               end if
            enddo
         else
            do i = 1,nprlbl
               label = prplbl(i)
               if (label(1:2) .eq. 'AM') then
                  natmlbl = natmlbl + 1
               end if
            enddo
         end if
         call memget('REAL',kint,nnbasx*natmlbl,work,kfree,lfree)
         call memget('REAL',kgso_prod,nnbasx*nbuf2,work,kfree,lfree)
         call memget('REAL',kw,nbuf2,work,kfree,lfree)
         call memget('REAL',kw_mbis,natmlbl*nbuf2,work,kfree,lfree)

         call num_aoint_opt(work(kint),work(kgso_prod),
     &                  work(kw),work(kw_mbis),
     &                  luddump,luwdump,cc_flag,nbuf2,natmlbl)

      end if
      call memrel('inside numint',work,kfrsav,kfrsav,kfree,lfree)

      call gpclose(luprop,'KEEP')
      call gpclose(luddump,'KEEP')
      call gpclose(luwdump,'KEEP')

      call newtimer('numint_driver')
      call qexit('numint_driver')

      return
      end


      subroutine evaluate_grid(corx,cory,corz,int_weight,
     &   gso,gso_prod,ncnt,dmatgao,ao_density,luddump,nbuf,
     &   work,lwork,cc_flag,save_mem,ntot_screened)
#include "implicit.h"
#include "dftcom.h"
#include "inforb.h"
#include "dummy.h"
#include "rspprp.h"
#include "inflr.h"
#include "priunit.h"
! Used from include:
! implicit.h:     
! dftcom.h        dftels, dfthri, dfthr0
! inforb.h:       nrhft, nnbasx (ccorb.h for cc?)
! dummy.h:        idummy
! rspprp.h:       MXLLR (needed for inflr.h)
! inflr.h:        iprlr
! priunit.h:      lupri
!

      parameter (d0 = 0.0d0)
! ingoing
      real*8    corx(nbuf),cory(nbuf),corz(nbuf),int_weight(nbuf),
     &          gso(nbast),ncnt(nbast),dmatgao(nbast),
     &          ao_density(nbast,nbast),work(lwork), gso_prod(nnbasx)
      integer   luddump,nbuf
      logical   cc_flag, save_mem

! defined within routine
      real*8    cor(3)
      integer   ntot_screened,atm_points,ipnt,idx,i,j
      real*8    nelctrn,weight,rho

! open file where dft grid is stored 
      luquad = -1
      call gpopen(luquad,'DALTON.QUAD','OLD','SEQUENTIAL',
     &   'UNFORMATTED',idummy,.false.)

!     write(*,*) 'ao dens in eval grid'
!     do i = 1,nbast
!        write(*,*) (ao_density(i,j),j=1,nbast)
!     enddo

      if (cc_flag) then
         do i = 1,nbast
            do j = 1,nbast
               ao_density(i,j) = 0.5d0 * 
     &            (ao_density(i,j) + ao_density(j,i))
               ao_density(j,i) = ao_density(i,j)
            end do
         end do
!        write(*,*) 'CC ao dens in evaluate_grid after sym'
!        do i = 1,nbast
!           write(*,*) (ao_density(i,j),j=1,nbast)
!        end do
      end if

      nelctrn = d0

      rewind luddump

!     write(*,*) 'i,rho,intweight print' !debug

      ntot_screened = 0
  100 continue
      read(luquad) atm_points
      if (atm_points .gt. 0) then
         call reaqua_srdft(corx,cory,corz,int_weight,luquad,atm_points)
         do 200 ipnt = 1,atm_points

            if (iprlr .gt. 150) then
               write (lupri,'(I8,4F12.6)')
     &         ipnt,corx(ipnt),cory(ipnt),corz(ipnt),int_weight(ipnt)
            end if

            weight = int_weight(ipnt)
            cor(1) = corx(ipnt)
            cor(2) = cory(ipnt)
            cor(3) = corz(ipnt)
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
            if (rho .lt. dfthr0) cycle 
            nelctrn = nelctrn + weight*rho

            if (save_mem) then
               write(luddump) corx(ipnt),cory(ipnt),corz(ipnt),
     &               int_weight(ipnt),rho,(gso(j),j=1,nbast)
            else
               do i = 1,nbast
                  idx = (i*i-i)/2
                  do j = 1,i
                     gso_prod(idx+j) = gso(i)*gso(j)*int_weight(ipnt)
                  enddo
               enddo
               write(luddump) corx(ipnt),cory(ipnt),corz(ipnt),
     &               int_weight(ipnt),rho,(gso_prod(j),j=1,nnbasx)
            endif

!           write(*,*) ipnt, rho, weight !debug

            ntot_screened = ntot_screened + 1

  200    continue
         goto 100
      else if (atm_points .eq. 0) then
         goto 100
      end if

      if (atm_points .ne. -1) then
         write(*,*) 'Error in reading, footer is', atm_points 
         call quit("Error reading file in evaluate_grid")
      end if

!     write(*,*) 'print done' !debug

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
       !zan(i)=atmchg(i) !changed to symmetry-proof charge array, infrsp.h ! doesnt work for cc, needs something general
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
      read(luddump,end=100) corx(ipt),cory(ipt),corz(ipt),
     &                      int_weight(ipt),rho(ipt)
      suma = suma + int_weight(ipt) * rho(ipt)
      ipt = ipt + 1
      if (ipt-1.gt.kpt) then
         write(lupri,*) 'MBIS: too many points',ipt,kpt
         call quit('MBIS: too many points')
      endif
      goto 10

100   continue

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
        write(*,30)i,lmax(i),(zeff(j,i),j=1,lmax(i))
      enddo
      write(*,*)'Qshell initial' ! N_Ai 
      do i=1,nat
        write(*,30)i,lmax(i),(qshell(j,i),j=1,lmax(i))
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
!     write(*,*) 'i,mbis print' !debug

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
!        write(luwdump) (mbis_weight(j)*int_weight(i),j=1,nat)
!        write(*,*) i,(mbis_weight(j),j=1,nat) !debug
      enddo !ipt

      write(luwdump) -1
!     write(*,*) 'print done' !debug

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
         zeff(1,i) = 2.0d0*znuc !used for everythin except ghost, where we overwrite
         if (znuc.le.1.0d-8) then 
         ! check for dummy atom - seems dangerous to check agains an exact zero
            lmax(i) = 1
            zeff(1,i) = 1.0d0
            qshell(1,i)= 1.0d-3
            qshell(2,i)=d0
         end if
         if (znuc.gt.1.0d-8 .and. znuc.le.2.0d0) then
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
         read(luddump,end=100) tmp,tmp,tmp,int_weight,tmp,
     &                 (gso(j),j=1,nbast)
         read(luwdump) (tmp,j=1,atm_idx-1), mbis_weight
         do a=1,nbast
            ! index of the last element in previous row, i.e. the diagonal
            idx=(a*a-a)/2 
            do b=1,a
               ao_int(idx+b) = ao_int(idx+b) + 
     &            gso(a)*gso(b)*mbis_weight*int_weight
            enddo
         enddo
         goto 10

      else if (cor_idx.gt.0) then

20       continue
         read(luddump,end=100) cor(1),cor(2),cor(3),int_weight,tmp,
     &                     (gso(j),j=1,nbast)
         read(luwdump) (tmp,j=1,atm_idx-1), mbis_weight
         dist = cor(cor_idx) - cord(cor_idx,atm_idx)
         do a=1,nbast
            idx=(a*a-a)/2
            do b=1,a
               ao_int(idx+b) = ao_int(idx+b) + 
     &            gso(a)*gso(b)*mbis_weight*int_weight*dist
            enddo
         enddo
         goto 20

      end if 

100   continue

!     if (ipt .ne. ipt_check) 
!    &   call quit("Error in reading of files in numint")

!     read(luwdump) ipt_check
!     if (ipt_check .ne. -1) 
!    &   call quit("Error in reading of files in numint")
      
!     CALL DSPTSI(NBAST,ao_int,ao_test)
!     write(lupri,*) 'ao_int test'
!     CALL OUTPUT(ao_test,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      write(*,*) 'ao int test'
      do i = 1,nnbasx
         write(*,*) ao_int(i)
      end do

      call getdat(rtnlbl(1),rtnlbl(2))

      rtnlbl(2) = 'SYMMETRI'
      !rtnlbl(1) holds opperator sym, we ignore that for now?
         !rtnlbl(1) is not read in in getgpv, seems like its inferred from integral vals.?
      
      call wrtpro(ao_int,nnbasx,integral_label,rtnlbl,iprlr)

      return
      end

      subroutine num_aoint_opt(ao_int, gso_prod,
     &           int_weight, mbis_weight,
     &           luddump, luwdump, cc_flag, nbuf2, natmlbl)
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "priunit.h"
#include "rspprp.h"
#include "inflr.h"
#include "ccrspprp.h"
      parameter (d0 = 0.0d0)
!incoming
      integer     luddump,luwdump,nbuf2,natmlbl
      real*8      ao_int(nnbasx,natmlbl),gso_prod(nnbasx,nbuf2),
     &            int_weight(nbuf2),mbis_weight(nbuf2,natmlbl)
      logical     cc_flag
!in here
      character*8 label, rtnlbl(2), prplbl_atm(natmlbl)
      integer     npoints, lbl_idx, block_size, points_read, footer, 
     &            atm_idx(natmlbl), atm_chg(natmlbl), cor_idx(natmlbl),
     &            i,j,k,n,idx,offset
      real*8      dist,scal,val, ao_test(nnbasx)

! parse all labels in here, so we can calculate all integrals with minimal file i/o
!     call newtimer('START')
      lbl_idx = 0
      if (cc_flag) then
         do i = 1,nprlbl_cc
            label = prplbl_cc(i)
            if (label(1:2) .eq. 'AM') then
               lbl_idx = lbl_idx + 1
               prplbl_atm(lbl_idx) = label
               read(label(3:4),*) atm_idx(lbl_idx)
               read(label(5:6),*) atm_chg(lbl_idx)
               if (label(7:8) .eq. '00') then
                  cor_idx(lbl_idx) = 0
               else if (label(7:8) .eq. '0x') then
                  cor_idx(lbl_idx) = 1
               else if (label(7:8) .eq. '0y') then
                  cor_idx(lbl_idx) = 2
               else if (label(7:8) .eq. '0z') then
                  cor_idx(lbl_idx) = 3
               else
                  write(*,*) 'Label ', label, ' not recognized'
                  call quit("Label not recognized in numint_write")
               end if
            endif
         enddo
      else
         do i = 1,nprlbl
            label = prplbl(i)
            if (label(1:2) .eq. 'AM') then
               lbl_idx = lbl_idx + 1
               prplbl_atm(lbl_idx) = label
               read(label(3:4),*) atm_idx(lbl_idx)
               read(label(5:6),*) atm_chg(lbl_idx)
               if (label(7:8) .eq. '00') then
                  cor_idx(lbl_idx) = 0
               else if (label(7:8) .eq. '0x') then
                  cor_idx(lbl_idx) = 1
               else if (label(7:8) .eq. '0y') then
                  cor_idx(lbl_idx) = 2
               else if (label(7:8) .eq. '0z') then
                  cor_idx(lbl_idx) = 3
               else
                  write(*,*) 'Label ', label, ' not recognized'
                  call quit("Label not recognized in numint_write")
               end if
            endif
         enddo
      end if

      ao_int = d0
      
! luwdump is file for MBIS weights
!  first line is npoints
!  then npoints lines with nucdep MBIS weights
!  final line footer is -1
      rewind(luwdump)
      read(luwdump) npoints
! luddump is file for grid 
!  npoints lines with coords, int. weight, density, orbitals
!  no header or footer, so reading is controlled through luwdump
      rewind(luddump)

      points_read = 0
      block_size = 0

      do while (points_read .lt. npoints)
         call read_block(luddump,luwdump,gso_prod,int_weight,
     &            mbis_weight,nbuf2,points_read,npoints,block_size,
     &            cor_idx,atm_idx,natmlbl)

!        call dgemm('N', 'N', nnbasx, natmlbl, block_size,
         call dgemm_test('N', 'N', nnbasx, natmlbl, nbuf2,
     &     1.0d0, gso_prod, nnbasx,
!    &     mbis_weight, block_size,
     &     mbis_weight, nbuf2,
     &     1.0d0, ao_int, nnbasx)
      enddo
!        do k=1,natmlbl
!        ! we first try for only cor_idx = 0
!        ! from changes to grid routine we now print gso_prod*int_weight(nnbasx,r)
!        ! so we can multiply with mbis_weight(r,nucdep) to get a (nnbasx,nucdep)
!           call dgemm('N', 'N', nnbasx, 1, block_size,
!    &         1.0d0, gso_prod, nnbasx,
!    &         mbis_weight(1,atm_idx(k)), nbuf2,
!    &         1.0d0, ao_int(1,k), nnbasx)
!           
!        enddo
!     enddo

!        do k=1,natmlbl
!           if (cor_idx(k) .eq. 0) then
!              do n=1,block_size
!                 offset = 1 + (n-1)*nbast
!                 scal = mbis_weight(atm_idx(k),n)*int_weight(n)
!                 call dcopy(nbast,gso(offset),1,gso_scaled(offset),1)
!                 call dscal(nbast,scal,gso_scaled(offset),1)
!              enddo
!           else
!              do n=1,block_size
!                 offset = 1 + (n-1)*nbast
!                 dist = cor(cor_idx(k),n) - cord(cor_idx(k),atm_idx(k))
!                 scal = mbis_weight(atm_idx(k),n)*int_weight(n)*dist
!                 call dcopy(nbast,gso(offset),1,gso_scaled(offset),1)
!                 call dscal(nbast,scal,gso_scaled(offset),1)
!              enddo
!           end if
!           do i=1,nbast
!              do j=1,i 
!                 idx = (i*(i-1))/2 + j
!                 val = ddot_opt(block_size,gso_scaled(i),nbast,
!    &                       gso(j),nbast) ! takes to long, need to redo mem. layout
!                 ao_int(idx,k) = ao_int(idx,k) + val
!              enddo
!           enddo
!        enddo
!     end do

      read(luwdump) footer
      if (footer .ne. -1) then
         write(*,*) 'Error in reading, footer is', footer
         call quit("Error reading file in num_aoint")
      end if

      call getdat(rtnlbl(1),rtnlbl(2))
      rtnlbl(2) = 'SYMMETRI'

      do k=1,natmlbl
         write(*,*) 'Writing lbl:', prplbl_atm(k)
         ao_test = ao_int(1:nnbasx,k)
         call wrtpro(ao_test,nnbasx,prplbl_atm(k),rtnlbl,iprlr)
!        call wrtpro(ao_int(1,k),nnbasx,prplbl_atm(k),rtnlbl,iprlr)
         write(*,*) 'ao int test'
         do i = 1,nnbasx
            write(*,*) ao_int(i,k)
         end do
      enddo

!     call newtimer('num_aoint_opt')

      return
      end

      subroutine read_block(luddump,luwdump,gso_prod,int_weight,
     &      mbis_weight,max_block,points_read,npoints,block_size,
     &      cor_idx,atm_idx,natmlbl)
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "priunit.h"
! incoming      
      integer luddump, luwdump, max_block, points, block_size,
     &         points_read, cor_idx(natmlbl), atm_idx(natmlbl),
     &         natmlbl
      real*8   int_weight(max_block),
     &         mbis_weight(max_block,natmlbl), 
     &         gso_prod(nnbasx,max_block)
! in here
      integer points_remaining, i, j
      real*8 tmp, tmp_weight(nucdep), cor(3)

      points_remaining = npoints - points_read
      max_read = min(points_remaining, max_block)

! 0 initilization important for the final read, which does not fill entire matrix
      mbis_weight = 0.0d0
      gso_prod = 0.0d0

      if (maxval(cor_idx) .eq. 0) then
         if (nucdep .eq. natmlbl) then
            do i=1,max_read
               read(luddump) cor(1), cor(2), cor(3),int_weight(i),
     &         tmp, (gso_prod(j,i), j=1,nnbasx) ! tmp for rho
               read(luwdump) (mbis_weight(i,j), j=1,nucdep)
            enddo
         else
            do i=1,max_read
               read(luddump) cor(1), cor(2), cor(3),int_weight(i),
     &         tmp, (gso_prod(j,i), j=1,nnbasx)
               read(luwdump) (tmp_weight(j), j=1,nucdep)
               do k=1,natmlbl
                  mbis_weight(i,k) = tmp_weight(atm_idx(k))
               enddo
            enddo
         endif
      else
         do i=1,max_read
            read(luddump) cor(1), cor(2), cor(3), int_weight(i),
     &      tmp, (gso_prod(j,i), j=1,nnbasx)
            read(luwdump) (tmp_weight(j), j=1,nucdep)
            do k=1,natmlbl
               if (cor_idx(k) .eq. 0) then
                  mbis_weight(i,k) = tmp_weight(atm_idx(k))
               else
                  dist = cor(cor_idx(k)) - cord(cor_idx(k),atm_idx(k))
                  mbis_weight(i,k) = tmp_weight(atm_idx(k)) * dist
               endif
            enddo
         enddo
      endif

      points_read = points_read + max_read
      ! block_size always = max_block except in the last read
      block_size = max_read
      write(*,*) 'block size:', block_size

      return
      end

      SUBROUTINE dgemm_test(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     +         BETA,C,LDC)
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA and NROWB  as the number of rows of  A
*     and  B  respectively.
*
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      IF (nota) THEN
          nrowa = m
      ELSE
          nrowa = k
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
*
*     Test the input parameters.
*
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND.
     +    (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND.
     +         (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And if  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (notb) THEN
          IF (nota) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (nota) THEN
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 170 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 130 i = 1,m
                          c(i,j) = zero
  130                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 140 i = 1,m
                          c(i,j) = beta*c(i,j)
  140                 CONTINUE
                  END IF
                  DO 160 l = 1,k
                      temp = alpha*b(j,l)
                      DO 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 200 j = 1,n
                  DO 190 i = 1,m
                      temp = zero
                      DO 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  180                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMM
*
      END
