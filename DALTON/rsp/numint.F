      subroutine numint_driver(integral_label,work,lwork) !,iorder,propri,iprint)
#include "implicit.h"
#include "dftinf.h"
#include "dftcom.h"
#include "inforb.h"
#include "inftap.h"
#include "dummy.h"
!
! dftinf.h:       nder, ntypso 
! dftcom.h:       dftgrid_done_old, dfthri
! inforb.h:       nbast, nnashx, n2basx, nisht, nasht, ncmot (TODO: ccorb.h is a seperate file, is that one neccesary?)
! inftap.h:       lusifc, lbsifc
! dummy.h:        idummy
!
      logical do_dcao, do_dvao
      integer luddump, luquad
      character integral_label*8
      dimension work(lwork)

      logical, save :: first_call = .true.

      call qenter('numint')

      if (first_call) then

! set up dft orbital common block so we don't evaluate uneccesary terms 
! TODO: maybe these should be saved and set back to original val?
! TODO: look into calling DFTINI or SETUPSOS instead of this
         DFTHR0 = 1.0D-8
         DFTHRL = 1.0D-10
!        DFTHRI = 1.0D-13
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

! check for active and inactive orbitals
         do_dcao = (nisht.gt.0)
         do_dvao = (nasht.gt.0)

! allocate static memory for 1 000 000 grid points
! this seems like its nbuf pr atom, so can probably be lowered?
! but is set to 1 000 000 in grid routine
         nbuf = 1 000 000
      !nbuf = 100 000 

! taken from prpget, not sure the purpose 
         kfrsav = 1
         kfree  = kfrsav
         lfree  = lwork 

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
! allocations for density mat. in ao basis
         if (do_dvao) then
            call memget('REAL',kdv,nnashx,work,kfree,lfree) 
            call memget('REAL',kdvao,n2basx,work,kfree,lfree) 
         else
            call memget('REAL',kdv,0,work,kfree,lfree) 
            call memget('REAL',kdvao,0,work,kfree,lfree) 
         end if
         if (do_dcao) then
            call memget('REAL',kdcao,n2basx,work,kfree,lfree) 
            call memget('REAL',kcmo,ncmot,work,kfree,lfree)
         else
            call memget('REAL',kdcao,0,work,kfree,lfree) 
         end if

         write(*,*) 'lfree after allocation', lfree
      
! read in mo coefficients and dv from sirius
! lusifc file structure explained in RD_SIRIFC
         if (lusifc .le. 0) call quit('SIRIFC not open in numint')
         rewind lusifc 
         call mollab(lbsifc,lusifc,lupri)
         read (lusifc)
         read (lusifc)
         if (do_dcao) then
! TODO: might need mo's for dvao calc. if we want uncoupled as well
            call readt (lusifc,ncmot,work(kcmo))
         else
            read (lusifc)
         end if
         if (do_dvao) then
            read (lusifc)
            call readt (lusifc,nnashx,work(kdv))
         end if

         write(*,*) 'cmo', (work(kcmo+i-1),i=1,ncmot)

! calculate dc from mo coeff. and backtransform dv and dc to ao basis
      ! if possible we can memrel kdv here
         call fckden(do_dcao,do_dvao,work(kdcao),work(kdvao),work(kcmo),
     &               work(kdv),work,lwork)

         write(*,*) 'dcao', (work(kdcao+i-1),i=1,n2basx)

! open file for dumping density and evaluated orbitals
         luddump = -1
         call gpopen(luddump,'DENSITY.GRID','UNKNOWN',
     &               ' ','UNFORMATTED',idummy,.false.)


         call numint_evaluate(work(kx),work(ky),work(kz),work(kw),
     &   work(kgso),work(kcnt),work(kdga),work(kdvao),work(kdcao),
     &   do_dvao,do_dcao,luddump,nbuf,work,lwork)

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
         call memget('REAL',kw_mbis,max_pnt*max_atm,work,kfree,lfree)

         write(*,*) 'lfree after mbis', lfree

         call mbis_pop(work(kx),work(ky),work(kz),work(kw),
     &                 work(krho),work(kw_mbis),
     &                 max_pnt,max_atm,luddump)

         call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)

         call gpclose(luddump,'KEEP')

         first_call = .false.


      end if ! first_call
      
!     call memget('REAL',kbas,max_pnt*nbast,work,kfree,lfree)
!     call numint_write

      call qexit('numint')

      call quit('end of numint, implementation not done')

      return
      end


      subroutine numint_evaluate(corx,cory,corz,int_weight,
     &   gso,ncnt,dmatgao,dvao,dcao,do_dvao,do_dcao,luddump,nbuf,
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
     &          dvao(nbast,nbast),dcao(nbast,nbast),
     &          work(lwork)
      logical   do_dvao,do_dcao
      integer   luddump,nbuf

! defined within routine
      real*8    cor(3)
      integer   npnts,npoint,ipnt
      real*8    celectrn,velectrn,weight,rho_inactive,rho_active,rho_tot

! open file where dft grid is stored 
      luquad = -1
      call gpopen(luquad,'DALTON.QUAD','OLD','SEQUENTIAL',
     &   'UNFORMATTED',idummy,.false.)

      celctrn = d0
      velctrn = d0

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
            write(*,*) 'getsos test'
            write(*,*) 'gso', gso
            write(*,*) 'cor', cor
!           CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,
!    &                  NBAST,DO_LONDON,DOGGA,THRINT,IPRINT)

! get density in point
            rho_inactive = d0
            rho_active = d0
            if (do_dcao) then
!           ... inactive part
               call getrho_srdft(dcao,gso,rho_inactive,dmatgao,
     &                           thrint,iprlr)
               celctrn = celctrn + weight*rho_inactive
            end if
!           ... active part
            if (do_dvao) then
               call getrho_srdft(dvao,gso,rho_active,dmatgao,
     &                           thrint,iprlr)
               velctrn = velctrn + weight*rho_active
            end if
            rho_tot = rho_inactive + rho_active

            write(luddump) corx(ipnt),cory(ipnt),corz(ipnt),
     &           int_weight(ipnt),rho_tot,(gso(j),j=1,nbast)

  200    continue
         goto 100
      else if (npoint .eq. 0) then
         goto 100
      end if
      write(luddump) npoint ! writing the negative npoint val as a footer?

      call gpclose(luquad,'KEEP')

! testing we integrate correct num. of electrons
      elctrx = float(2*nrhft)
      elctrn = celctrn + velctrn
      error  = elctrn - elctrx

! print section for num int
      write (lupri,'(3(/2X,A,F14.6))')
     &' Nr. of inactive electrons from numerical integration :',celctrn,
     &' Nr. of active electrons from numerical integration   :',velctrn,
     &' Total nr. of electrons from numerical integration    :',elctrn
      write (lupri,'(2(/2X,A,F14.6))')
     &' Number of electrons from orbial occupations          :',elctrx,
     &' Error in the number of electrons                     :',error
      if (abs(error) .gt. dftels) then
! dftels is 1.0d-3 as standard
         write (lupri,'(/2X,A,F14.6,/2X,A)')
     &' Error larger than DFTELS (set input)                 :',dftels,
     &   ' Calculation aborted.'
         call quit
     &    ('Wrong number of electrons in Numint. Calculation aborted.')
      end if

      return
      end


      subroutine mbis_pop(corx,cory,corz,int_weight,rho,
     &                 mbis_weight,max_pnt,max_atm,luddump)
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
     &          mbis_weight(max_pnt,max_atm)
      integer   max_pnt,max_atm,luddump
! defined in here, some should be allocated with memget instead
! 10 used as max l, but currently only imp. for 3
      integer   npoint,lmax(max_atm),ipt
      real*8    zan(max_atm),qatm(max_atm),qtmp(max_atm),
     &          zeff(10,max_atm),zstmp(10,max_atm),
     &          qshell(10,max_atm),qstmp(10,max_atm),
     &          vshell(10,max_atm),vatm(max_atm),wtatm(10,max_atm),
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

      write(*,*) 'rho,int_weight'
      write(*,*) (rho(i),i=1000,1005)
      write(*,*) (int_weight(i),i=1000,1005)

      write(*,*) 'luddump footer test', npoint

! tot grid points = ipt-1 due to the way the reading is constructed
      ipt=ipt-1
      write(lupri,*) ' grid points, tot elec.'
      write(lupri,14) ipt,suma
      write(lupri,*)' '
14    format(i10,2f20.10)
30    format(2i5,5f12.4)
      
c     estimate initial Zeff
      call zeffini(kpar,nat,zan,zeff,qshell,lmax)
      write(*,*)'Zeff initial' !relateret til sigma
      do i=1,nat
        write(*,30)i,lmax(i),(zeff(j,i),j=1,3)
      enddo
      write(*,*)'Qshell initial' !samme som population, N(Ai) - også estimeret i zeffini
      do i=1,nat
        write(*,30)i,lmax(i),(qshell(j,i),j=1,3)
      enddo

c     iterate until MBIS charges are converged
      do kiter=1,kmax
        mbis_weight = d0 
        write(*,*)' '
        write(*,*)' Iteration ',kiter
c       update qshell
        do n=1,nat
          do l=1,lmax(n)
            qstmp(l,n)= d0 !reset tmp vector for N(Ai) 
          enddo
        enddo
c       loop over points
      suma=d0
        do i=1,ipt
      suma=suma+int_weight(i)*rho(i)
          rhom = rho(i)
          fwt = int_weight(i) 
!         jcnt = icnt(i)
!         xj = coor(1,jcnt)
!         yj = coor(2,jcnt)
!         zj = coor(3,jcnt)

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
              znorm = (ztmp**3)/epi !epi??
              xtmp = znorm*exp(-ztmp*dis)
              rhoa = qshell(l,n)
              tmp = xtmp*rhoa*fwt
              if (abs(tmp).lt.dencut) then
                   tmp=d0 !eps = limit where we set 0 (10^-14) to not include noise
              end if
              wtatm(l,n) = wtatm(l,n) + tmp !making p0(Ai), eq 7
!        if(tmp.lt.zero) write(*,*) ' tmp, wtatm', tmp, wtatm(l,n)
!                   if(first.eqv..true.) then
!                       write (*,*) 'debugo',ztmp,znorm,rhoa
!                       first=.false.
!                   end if 
            enddo
          enddo

c         assign ISA fractions and partition electrons
          wtot = d0
          do n=1,nat
            do l=1,lmax(n)
              wtot = wtot + wtatm(l,n) !p0 made from p0(Ai), e.g. eq 5, maybe only p0i? or p0(r)
            enddo
          enddo
          if (abs(wtot).gt.d0 .and. (abs(rhom*fwt).gt.eps)) then
!           tmp = zero
            do n=1,nat
              do l=1,lmax(n)
                qstmp(l,n) = qstmp(l,n) + wtatm(l,n)*fwt*rhom/wtot !N(Ai) population, eq 18 
!               tmp = tmp + wtatm(l,n) 
!               wpnt(i,n) = wpnt(i,n) + wtatm(l,n)/wtot
                mbis_weight(i,n) = mbis_weight(i,n) + wtatm(l,n)/wtot

!               write(*,*) 'n,l,wpnt,wtatm,wtot,rhom,ipt', 
!    &             n, l, wpnt(i,n), wtatm(l,n), wtot,rhom, i

        if (mbis_weight(i,n) .lt. d0) write(*,*) 'wpnt lt zero'
              enddo
            enddo
!           tmp = tmp/wtot
!           if (abs(1.0d0-tmp).gt.1.0d-10) then
!             write(*,*)' weight sum not one'
!           endif
!         else
!           if (wtot.lt.zero) then
!             write(*,*)' wtot less than zero???'
!          !endif
!           if (rhom*fwt.lt.zero) then
!             write(*,*)' rhom*fwt less than zero???'
!             write(*,*)rhom,fwt
!           endif
          endif
        enddo
c       end ipt loop
!       write(*,*) 'mbis sum 1', wt1, ',mbis sum 2', wt2

      tote=suma!+sumb
      write(*,*)'tote', tote
        write(*,*)' Qshell'
        tote=d0
        do n=1,nat
          do l=1,lmax(n)
            tote = tote + qstmp(l,n) !summation over amount of elec.
            tmp = qstmp(l,n)-qshell(l,n) !change of N(Ai) in this iteration
            write(*,21)n,l,qshell(l,n),qstmp(l,n),tmp 
          enddo
        enddo
c       write(*,*)' electron count',tote
        rmsq=d0
        qtot=d0
        qele=d0
        write(*,*)' Qatom'
        do n=1,nat
          qtmp(n)=d0
          do l=1,lmax(n)
            qtmp(n)=qtmp(n)+qstmp(l,n) !N(Ai) to N(A)?, which is amount of electrons on atom A?
          enddo
          qele = qele + qtmp(n) !amount of electrons on A or n
          qtmp(n) = zan(n) - qtmp(n) !Q(A)=Z(A)-N(A) (eq. 18) because shape-functions are normalized (integral over f = 1, text after eq. 7)
          tmp = qtmp(n)-qatm(n) !Q from previouis iteration, 0 first time
          rmsq = rmsq + tmp*tmp !Eq. 20 but shifted with Z, so looking at the total charge instead of electron charge
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
22      format(' RMS-dQ , Qele, Qtot, Qxs',4f15.8)

c       update zeff
        do n=1,nat
          do l=1,lmax(n)
            zstmp(l,n)=d0
            vshell(l,n)=d0
          enddo
        enddo

c       loop over points
        do i=1,ipt
          rhom = rho(i)
          fwt = int_weight(i)
!         jcnt = icnt(i)
!         xj = coor(1,jcnt)
!         yj = coor(2,jcnt)
!         zj = coor(3,jcnt)

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
          if (abs(wtot).gt.d0 .and. abs(rhom*fwt).gt.eps) then !abs(rhom*fwt)
            do n=1,nat
               dx = corx(i) - coor(1,n)
               dy = cory(i) - coor(2,n)
               dz = corz(i) - coor(3,n) 
              dis = dx*dx + dy*dy + dz*dz
              dis = sqrt(dis)
              vtmp = zan(n)/dis
              do l=1,lmax(n)
                zstmp(l,n) = zstmp(l,n) + wtatm(l,n)*rhom*dis*fwt/wtot !looks similar to equation 18, just without 1/3N(Ai)
c               vshell(l,n) = vshell(l,n) + wtatm(l,n)*vtmp*fwt/wtot
             vshell(l,n) = vshell(l,n) + wtatm(l,n)*vtmp*rhom*fwt/wtot !what? something with valence maybe? or maybe just presented in another form, not actually used in itterations?
              enddo
            enddo
          endif
        enddo
c       end ipt loop

        write(*,*)' Zshell'
!       write(lupri,*)' Zshell'
        do n=1,nat
          do l=1,lmax(n)
            zstmp(l,n) = (3.0d0*qshell(l,n))/zstmp(l,n)
            tmp = zstmp(l,n)-zeff(l,n)
                write(*,*) 'zeff debug', l, n, zeff(l,n)
            write(*,21)n,l,zeff(l,n),zstmp(l,n),tmp
!           write(lupri,21)n,l,zeff(l,n),zstmp(l,n),tmp
            zeff(l,n) = zstmp(l,n) !zeff corresponds to 1/sigma in eq 19
c
            vshell(l,n) = vshell(l,n)/qshell(l,n) !vshell is leftover from old code 
c
          enddo
        enddo
        write(*,*)' Vshell'
        do n=1,nat
          vatm(n)=d0
          do l=1,lmax(n)
                write(*,*) 'vshell debug', l, n, vshell(l,n)
            vatm(n)=vatm(n)+vshell(l,n)
            write(*,21)n,l,vshell(l,n)
          enddo
        enddo
        write(*,*)' Vatom'
        do n=1,nat
          write(*,20)n,vatm(n)
        enddo

        write(*,*)' '
      enddo

20    format(i5,3f15.8)
21    format(2i5,3f15.8)

900   continue

      write(*,*) 'MBIS converged!'

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
      write(*,*)' Vatom'
      do n=1,nat
        write(*,20)n,vatm(n)
      enddo

c     ISA weighting
        do n=1,nat
          do l=1,lmax(n)
            vshell(l,n)=d0
          enddo
        enddo
        do i=1,ipt
          fwt = int_weight(i)
!         jcnt = icnt(i)
!         xj = coor(1,jcnt)
!         yj = coor(2,jcnt)
!         zj = coor(3,jcnt)

          do n=1,nat
            dx = corx(i) - coor(1,n)
            dy = cory(i) - coor(2,n)
            dz = corz(i) - coor(3,n)
            dis = dx*dx + dy*dy + dz*dz
            dis = sqrt(dis)
            vtmp = d0
            if(dis.gt.eps)vtmp = zan(n)/dis !eps or dencut?
            do l=1,lmax(n)
              ztmp = zeff(l,n)
              znorm = (ztmp**3)/epi
              xtmp = znorm*exp(-ztmp*dis)
              rhoa = qshell(l,n)
              tmp = xtmp*rhoa*fwt
              if (tmp.lt.dencut) tmp=d0 !eps or dencut?
              vshell(l,n) = vshell(l,n) + vtmp*tmp
            enddo
          enddo
        enddo
c       end ipt loop
        write(*,*)' Vshell-ISA   '
        do n=1,nat
          vatm(n)=d0
          do l=1,lmax(n)
            vshell(l,n) = vshell(l,n)/qshell(l,n)
            vatm(n)=vatm(n)+vshell(l,n)
            write(*,21)n,l,vshell(l,n)
          enddo
        enddo
        write(*,*)' Vatom-ISA   '
        do n=1,nat
          write(*,20)n,vatm(n)
        enddo

        do n=1,nat
          vatm(n)=d0
          tmp=d0
          do l=1,lmax(n)
            vatm(n)=vatm(n)+vshell(l,n)
            tmp=tmp+qshell(l,n)
          enddo
          vatm(n)=vatm(n)/tmp
        enddo
        write(*,*)' Vatom-ISA2   '
        do n=1,nat
          write(*,20)n,vatm(n)
        enddo

! CALCULATE GPV USING MBIS WEIGHT



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


