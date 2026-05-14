      subroutine numint_driver(ao_density,work,lwork,cc_flag)
!*****************************************************************************
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
#include "infpar.h"
! Used from include:
! dftcom.h:       srDFTGRID_done, mxbllen
! inforb.h:       nbast, nnbasx, n2basx, nisht, nasht
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
      dimension work(lwork), ao_density(nbast,nbast)
! defined here:
      integer luwdump, kfrsav, kfree, lfree, nbuf
      logical old_grid, dft_isgga, dogga
      character(len=8) label
      real*8 ele

      external evaluate_blocked_grid
      external dft_isgga

! option for old or new grid - currently only adjustable through this hardcoded variable
! old grid is less precise and less adjustable but can run mpi
      old_grid = .false.

! timers and memory pointers
      call qenter('numint_driver')

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

      if (old_grid) then

         ! set up dft orbital common block and grid
         call setupsos(0,.false.,idum1,idum2)
         if (.not. srDFTGRID_done) then
            call make_dftgrid(work,lwork,ntot_dftgrid,iprlr,.false.)
            srDFTGRID_done = .true.
         end if

         ! make_dftgrid routine limits points per atom to 1.0d6,
         ! however can never need more than total grid size
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
     &      work(kgso),work(kgso_prod),work(kcnt),work(kdga),ao_density,
     &      nbuf,work,lwork,cc_flag,ntot_screened)

         call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)
         ntot_dftgrid = ntot_screened

      else ! use blocked grid from regular DFT code 

         dogga = dft_isgga()
         ntot_dftgrid = 0

         if (cc_flag) then
            do i = 1,nbast
               do j = 1,nbast
                  ao_density(i,j) = 0.5d0 * 
     &               (ao_density(i,j) + ao_density(j,i))
                  ao_density(j,i) = ao_density(i,j)
               end do
            end do
         end if

         call memget('REAL',kprod,mxbllen*nnbasx,work,kfree,lfree)
         call newtimer('START')

         if (dogga) then
            call dftint(ao_density,1,-1,.false.,work(kfree),lfree,
     &               evaluate_blocked_grid,work(kprod),ele)
         else
            call dftint(ao_density,1,0,.false.,work(kfree),lfree,
     &               evaluate_blocked_grid,work(kprod),ele)
         endif

         write(luddump) -1, -1
         call newtimer('dftint')
         call memrel('dftint',work,kfrsav,kfrsav,kfree,lfree)

      end if

! allocate for MBIS populations routine 
      max_pnt = min(1 000 000,ntot_dftgrid)
      max_atm = min(50,nucdep)
         
      call memget('REAL',kx,max_pnt,work,kfree,lfree)
      call memget('REAL',ky,max_pnt,work,kfree,lfree)
      call memget('REAL',kz,max_pnt,work,kfree,lfree)
      call memget('REAL',kw,max_pnt,work,kfree,lfree)

      call memget('REAL',krho,max_pnt,work,kfree,lfree)
      call memget('REAL',kw_mbis,max_atm,work,kfree,lfree)

! calculate MBIS population weights in each point
      call newtimer('START')
      call mbis_pop(work(kx),work(ky),work(kz),work(kw),
     &   work(krho),work(kw_mbis),
     &   max_pnt,max_atm,luwdump,old_grid)
      call newtimer('MBIS')

      call memrel('numint',work,kfrsav,kfrsav,kfree,lfree)
! allocate for and calculate ao integrals, memory saving version
      if (old_grid) then

         call memget('REAL',kint,nnbasx,work,kfree,lfree)
         call memget('REAL',kgso,nbast,work,kfree,lfree)

         if (cc_flag) then
            do i = 1,nprlbl_cc
               label = prplbl_cc(i)
               if (label(1:2) .eq. 'AM') then
                  call num_aoint(label,work(kint),work(kgso),
     &               luddump,luwdump)
               end if
            enddo
         else
            do i = 1,nprlbl
               label = prplbl(i)
               if (label(1:2) .eq. 'AM') then
                  call num_aoint(label,work(kint),work(kgso),
     &               luddump,luwdump)
               end if
            enddo
         end if

! same as above, but with blocked grid 
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
         call memget('REAL',kint2,nnbasx*natmlbl,work,kfree,lfree)
         call memget('REAL',kgso_prod,nnbasx*mxbllen,work,kfree,lfree)
         call memget('REAL',kw_mbis,natmlbl*mxbllen,work,kfree,lfree)

         call newtimer('START')
         call num_aoint_opt(work(kint),work(kint2),work(kgso_prod),
     &                  work(kw_mbis),luwdump,cc_flag,nbuf2,natmlbl)
         call newtimer('num_aoint_opt')

      end if

      call memrel('inside numint',work,kfrsav,kfrsav,kfree,lfree)

      call gpclose(luprop,'KEEP')
      call gpclose(luddump,'KEEP')
      call gpclose(luwdump,'KEEP')

      call qexit('numint_driver')

      return
      end

      subroutine evaluate_blocked_grid(nblen,nblcnt,nblocks,ldaib,gao,
     &                                 rhoa,grada,dst,vfa,xcpot,
     &                                 coord,wght,gao_prod)
!*****************************************************************************
!     Callback for dftint, evaluating one block of DFT grid points
!
!     Purpose:
!     For each block of grid points passed by dftint, determine which AOs
!     are active using the symmetry block lists, compute weighted AO products
!     chi_a(r)*chi_b(r)*w(r) for all active upper-triangular (a,b) pairs,
!     and write coordinates, weights, density, and AO products to
!     DENSITY.GRID. Accumulates the total grid point count in ntot_dftgrid.
!*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "inforb.h"

#include "dftinf.h"
#include "energy.h"
#include "symmet.h"
#include "shells.h"
#include "dftcom.h"

#include "inftap.h"
! Used from include:
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       natoms (for dst array dimension)
! inforb.h:       nbast, nnbasx, nsym
! symmet.h:       maxrep, nsym
! dftcom.h:       mxbllen, ntot_dftgrid
! inftap.h:       luddump
!
! Files used for:
! luddump:        writing grid coordinates, weights, density, and AO products
!
      dimension gao(nblen,nbast,*), coord(3,nblen), wght(nblen),
     &          rhoa(nblen), grada(3,nblen),
     &          nblcnt(8), nblocks(2,ldaib,8),
     &          dst(natoms), vfa(nblen), xcpot(nblen)
      logical active(0:nbast)
      dimension gao_prod(mxbllen*nnbasx)
      integer full_idx(nnbasx),idx,joint_idx,idx_tot

      if (maxrep .ne. 0) 
     &   call quit("no symmetry in evaluate_blocked_grid")

      ntot_dftgrid = ntot_dftgrid + nblen

      do i = 0,nbast
         active(i) = .false.
      enddo
      do isym = 1,nsym
         do ibla = 1,nblcnt(isym)
            do i = nblocks(1,ibla,isym),nblocks(2,ibla,isym)
               active(i) = .true.
            enddo
         enddo
      enddo

      idx = 0
      do i=1,nbast
         if (active(i)) then
            do j=1,i
               if (active(j)) then
                  idx = idx + 1
                  full_idx(idx) = ((i*i-i)/2)+j
               endif
            enddo
         endif
      enddo

      idx_tot = idx

      write(luddump) nblen, idx_tot ! number of grid points, and number of ao products
      write(luddump) (full_idx(i),i=1,idx_tot)

      gao_prod = 0.0d0
      do k=1,nblen
         idx = 0
         do i=1,nbast
            if (active(i)) then
               do j=1,i
                  if (active(j)) then
                     idx = idx+1
                     !joint_idx = k + idx*(nblen-1)
                     joint_idx = idx + (k-1)*idx_tot
                     gao_prod(joint_idx) = gao(k,i,1)*gao(k,j,1)*wght(k)
                  endif
               enddo
            endif
         enddo
      enddo

      write(luddump) coord(1,:), coord(2,:), coord(3,:), 
     &               wght, rhoa, (gao_prod(i),i=1,nblen*idx)

      return
      end

      subroutine evaluate_grid(corx,cory,corz,int_weight,
     &   gso,gso_prod,ncnt,dmatgao,ao_density,nbuf,
     &   work,lwork,cc_flag,ntot_screened)
!*****************************************************************************
!     Evaluate orbitals and density on the old (non-blocked) DFT grid
!
!     Purpose:
!     Read grid points from DALTON.QUAD, evaluate AO basis functions and the
!     electron density at each point, discard points below the density
!     threshold (dfthr0), and write surviving points with coordinates,
!     weights, density, and AO values to DENSITY.GRID. Reports the integrated
!     electron count and aborts if it deviates from the reference by more
!     than DFTELS.
!*****************************************************************************
#include "implicit.h"
#include "dftcom.h"
#include "inforb.h"
#include "dummy.h"
#include "rspprp.h"
#include "inflr.h"
#include "priunit.h"
#include "inftap.h"
! Used from include:
! dftcom.h:       dftels, dfthri, dfthr0
! inforb.h:       nrhft, nnbasx
! dummy.h:        idummy
! rspprp.h:       MXLLR (needed for inflr.h)
! inflr.h:        iprlr
! priunit.h:      lupri
! inftap.h:       luddump
!
! Files used for:
! luquad:         reading grid points from DALTON.QUAD (opened here)
! luddump:        writing coordinates, weights, density, and AO values
! lupri:          printing electron count check
!

      parameter (d0 = 0.0d0)
! ingoing
      real*8    corx(nbuf),cory(nbuf),corz(nbuf),int_weight(nbuf),
     &          gso(nbast),ncnt(nbast),dmatgao(nbast),
     &          ao_density(nbast,nbast),work(lwork), gso_prod(nnbasx)
      integer   nbuf
      logical   cc_flag

! defined within routine
      real*8    cor(3)
      integer   ntot_screened,atm_points,ipnt,idx,i,j
      real*8    nelctrn,weight,rho

! open file where dft grid is stored 
      luquad = -1
      call gpopen(luquad,'DALTON.QUAD','OLD','SEQUENTIAL',
     &   'UNFORMATTED',idummy,.false.)

      if (cc_flag) then
         do i = 1,nbast
            do j = 1,nbast
               ao_density(i,j) = 0.5d0 * 
     &            (ao_density(i,j) + ao_density(j,i))
               ao_density(j,i) = ao_density(i,j)
            end do
         end do
      end if

      nelctrn = d0

      rewind luddump

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
            thrint = dfthri/weight 

! evaluate orbitals in the grid point
            call getsos(gso,ncnt,cor,work,lwork,
     &                  nbast,.false.,.false.,thrint,iprlr)
!           call getsos(gso,ncnt,cor,work,lwork,
!    &                  nbast,do_london,dogga,thrint,iprlr)

! get density in point
            rho = d0
            call getrho_srdft(ao_density,gso,rho,dmatgao,
     &         thrint,iprlr)
            if (rho .lt. dfthr0) cycle 
            nelctrn = nelctrn + weight*rho

            write(luddump) corx(ipnt),cory(ipnt),corz(ipnt),
     &            int_weight(ipnt),rho,(gso(j),j=1,nbast)

            ntot_screened = ntot_screened + 1

  200    continue
         goto 100
      else if (atm_points .eq. 0) then
         goto 100
      end if

      if (atm_points .ne. -1) then
         write(lupri,*) 'Error in reading, footer is', atm_points 
         call quit("Error reading file in evaluate_grid")
      end if

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
     &                 mbis_weight,max_pnt,max_atm,luwdump,old_grid)
!*****************************************************************************
!     Iterative MBIS population analysis
!
!     Purpose:
!     Read grid coordinates, weights, and electron density from DENSITY.GRID,
!     then iteratively optimize pro-atomic exponential parameters (zeff) and
!     shell populations (qshell) until the MBIS charges converge. On
!     convergence, compute per-point atom weights
!     w_A(r) = rho_A^0(r) / sum_B rho_B^0(r) and write them to WEIGHT.GRID
!     for use in the AO integral routines.
!*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "priunit.h"
#include "infrsp.h"
#include "inforb.h"
#include "inftap.h"
! Used from include:
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       nucdep, charge (nuclear charge), cord (atomic coordinates)
! priunit.h:      lupri
! inftap.h:       luddump (read from DENSITY.GRID)
!
! Files used for:
! luddump:        reading grid coordinates, weights, and electron density
! luwdump:        writing per-point MBIS atom weights (passed as argument)
! lupri:          printing convergence progress and final results
!
! TODO:
! * most defined in here should be allocated with memget instead 
! * 10 used as max l, but currently only imp. for 3 in zeffini
! * introduce common block for max vals
! * use Dalton variables instead of redefining charge etc.
! * cleanup of print -> conditional lupri print on print level
      parameter (d0 = 0.0d0, pi=acos(-1.0d0), epi=8.0d0*pi,
     &           qconv=1.0d-9, eps=1.0d-14, dencut=1.0d-10)
! incoming
      real*8    corx(max_pnt),cory(max_pnt),corz(max_pnt),
     &          int_weight(max_pnt),rho(max_pnt),
     &          mbis_weight(max_atm)
      integer   max_pnt,max_atm
      logical   old_grid
! defined in here
      integer   npoint,lmax(max_atm),ipt,tmp2
      real*8    zan(max_atm),qatm(max_atm),qtmp(max_atm),
     &          zeff(10,max_atm),zstmp(10,max_atm),
     &          qshell(10,max_atm),qstmp(10,max_atm),
     &          wtatm(10,max_atm),
     &          suma,tote,qtot,dx,dy,dz,rhom,fwt,tmp,ztmp,xtmp,rmsq,dis,
     &          coor(3,max_atm), wtot

      kmax = 1000
      kpar = max_atm
      kpt = max_pnt

      nat=nucdep
      if (nat.gt.kpar) then
        write(lupri,*)'MBIS: too many atoms',nat,kpar
        call quit('MBIS: too many atoms')
      endif

! saving variables in other names to fit with old code
      do i=1,nat
        zan(i)=charge(i) !nuclear charge of atom, nuclei.h
        qatm(i)=d0
        do j=1,3
          coor(j,i)=cord(j,i) !cord(3,mxcent), nuclei.h
        end do
      enddo

      rewind(luddump)
      rewind(luwdump)

      ipt=1
      start_pnt=1
      end_pnt=0
      suma = d0
      if (old_grid) then
10       continue
         read(luddump,end=100) corx(ipt),cory(ipt),corz(ipt),
     &                      int_weight(ipt),rho(ipt)
         suma = suma + int_weight(ipt) * rho(ipt)
         ipt = ipt + 1
         if (ipt-1.gt.kpt) then
            write(lupri,*) 'MBIS: too many points',ipt,kpt
            call quit('MBIS: too many points')
         endif
         goto 10
100      continue
      else
11       continue
         read(luddump) npoint, tmp2 ! points in block, number of ao products
         if (npoint.gt.0) then
            read(luddump) tmp2 ! gao index array
            end_pnt = end_pnt + npoint
            read(luddump) corx(start_pnt:end_pnt), 
     &                    cory(start_pnt:end_pnt),
     &                    corz(start_pnt:end_pnt),
     &                    int_weight(start_pnt:end_pnt),
     &                    rho(start_pnt:end_pnt)
            do i = 1,npoint
               suma = suma + int_weight(start_pnt+i-1) 
     &            * rho(start_pnt+i-1)
            enddo
            start_pnt = end_pnt + 1
            go to 11
         else if (npoint.eq.0) then
            go to 11
         end if
110      continue
         ipt = end_pnt + 1
      end if


! tot grid points = ipt-1 due to the way the initial reading was constructed
      ipt=ipt-1
      write(lupri,'(a, i10)') 'Total grid points:        ', ipt
      write(lupri,'(a, f10.5)') 'Total elec. from num int: ', suma
      write(lupri,*)' '
30    format(2i5,5f12.4)
      
! guess initial parameters for pro-atomic density 
      call zeffini(kpar,nat,zan,zeff,qshell,lmax)
      write(lupri,*)'Alpha initial val. (1/sigma)'
      do i=1,nat
        write(lupri,30)i,lmax(i),(zeff(j,i),j=1,lmax(i))
      enddo
      write(lupri,*)'Population initial val.' ! N_Ai 
      do i=1,nat
        write(lupri,30)i,lmax(i),(qshell(j,i),j=1,lmax(i))
      enddo

! iterate until MBIS charges are converged
      do kiter=1,kmax
         write(lupri,*)' '
         write(lupri,*)' Iteration ',kiter

! in the first part of the itteration we calculate the population of the pro-atomic densities
         do n=1,nat
            do l=1,lmax(n)
               qstmp(l,n)= d0
            enddo
         enddo

         do i=1,ipt
            rhom = rho(i)
            fwt = int_weight(i) 

            do n=1,nat
               do l=1,lmax(n)
                  wtatm(l,n)=0.0d0
               enddo
            enddo
!         construct atomic weight function in the given point
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

         write(lupri,*)' Qshell'
         tote=d0
         do n=1,nat
            do l=1,lmax(n)
               tote = tote + qstmp(l,n)
               tmp = qstmp(l,n)-qshell(l,n)
               write(lupri,21)n,l,qshell(l,n),qstmp(l,n),tmp 
            enddo
         enddo
         write(lupri,*)' Electrons from integration:', suma
         write(lupri,*)' Electrons from MBIS       :', tote
         rmsq=d0
         qtot=d0
         qele=d0
         write(lupri,*)' Qatom'
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
            write(lupri,20)n,qatm(n),qtmp(n),tmp
         enddo
         qxs=nint(qele)
         qxs=qele+qtot-qxs
         rmsq=sqrt(rmsq)
         write(lupri,22)rmsq,qele,qtot,qxs
         if (rmsq.lt.qconv) goto 900
         do n=1,nat
            qatm(n) = qtmp(n) !N(A) saved for next iteration (for eq 20)
            do l=1,lmax(n)
               qshell(l,n)=qstmp(l,n) !New N(Ai) saved for next iteration to make p0(Ai)
            enddo
         enddo

!        update zeff
         do n=1,nat
            do l=1,lmax(n)
               zstmp(l,n)=d0
            enddo
         enddo

         do i=1,ipt
            rhom = rho(i)
            fwt = int_weight(i)

            do n=1,nat
               do l=1,lmax(n)
                  wtatm(l,n)=d0
               enddo
            enddo

!         reconstruct atomic weights functions - same as before
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
!         assign ISA fractions and partition electrons 
            wtot = d0
            do n=1,nat
               do l=1,lmax(n)
                  wtot = wtot + wtatm(l,n) 
               enddo
            enddo
            !if (wtot.gt.d0 .and. abs(rhom*fwt).gt.eps) then
            if (abs(wtot).gt.d0 .and. abs(rhom*fwt).gt.eps) then
               do n=1,nat
                  dx = corx(i) - coor(1,n)
                  dy = cory(i) - coor(2,n)
                  dz = corz(i) - coor(3,n) 
                  dis = dx*dx + dy*dy + dz*dz
                  dis = sqrt(dis)
                  do l=1,lmax(n)
                     zstmp(l,n) = zstmp(l,n) + 
     &                  wtatm(l,n)*rhom*dis*fwt/wtot 
                  enddo
               enddo
            endif
         enddo

         write(lupri,*)' Zshell'
         do n=1,nat
            do l=1,lmax(n)
               if (abs(zstmp(l,n)).gt.eps) then ! 10^-14
                  zstmp(l,n) = (3.0d0*qshell(l,n))/zstmp(l,n)
               else 
                  zstmp(l,n) = (3.0d0*qshell(l,n))/sign(eps,zstmp(l,n))
               endif
               tmp = zstmp(l,n)-zeff(l,n)
               write(lupri,21)n,l,zeff(l,n),zstmp(l,n),tmp
               zeff(l,n) = zstmp(l,n)
            enddo
         enddo

        write(lupri,*)' '
      enddo 

20    format(i5,3f15.8)
21    format(2i5,3f15.8)
22    format(' RMS-dQ , Qele, Qtot, Qxs',4f15.8)

900   continue

      write(lupri,*) ' MBIS converged!'
      write(lupri,*)' '
      write(lupri,*)' Final converged results'
      write(lupri,*)' Qatom'

      do n=1,nat
        qtmp(n)=d0
        do l=1,lmax(n)
          qtmp(n)=qtmp(n)+qstmp(l,n)
        enddo
        qtmp(n) = nint(zan(n)) - qtmp(n)
        write(lupri,20)n,qtmp(n)
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
         !if (wtot.gt.d0 .and. (abs(rhom*fwt).gt.eps)) then
         if (abs(wtot).gt.d0 .and. (abs(rhom*fwt).gt.eps)) then
            do n=1,nat
               do l=1,lmax(n)
                  mbis_weight(n) = mbis_weight(n) + wtatm(l,n)/wtot
               enddo
            enddo
         endif
         write(luwdump) (mbis_weight(j),j=1,nat)
      enddo 

      write(luwdump) -1

      return
      end

      subroutine zeffini(kpar,nat,zan,zeff,qshell,lmax)
!*****************************************************************************
!     Initialize pro-atomic density parameters for MBIS
!
!     Purpose:
!     Set initial values of the exponential decay parameters (zeff) and
!     shell populations (qshell) for each atom based on nuclear charge.
!     Atoms are assigned one shell (H-He), two shells (Li-Ne), or three
!     shells (Na-Ar); Z > 18 is not implemented. These starting values feed
!     the iterative MBIS procedure in mbis_pop.
!*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
      parameter(d0=0.0d0)
#include "nuclei.h"
! Used from include:
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       (no variables used directly; included for implicit sizing)
!
      real*8    zan(kpar),zeff(10,kpar),qshell(10,kpar),znuc
      integer   lmax(kpar)

      zeff = 0.0d0
      qshell = 0.0d0
      lmax = 0

      do i=1,nat
         znuc = zan(i)
         zeff(1,i) = 2.0d0*znuc

         ! Check for dummy atom 
         if (znuc.le.1.0d-8) then 
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
            qshell(3,i)=znuc-10.0d0
         endif

         if (znuc.gt.18) call quit('Z > 18 not implemented in zeffini')
      enddo
 
      return
      end

      subroutine num_aoint(integral_label,ao_int,gso,luddump,luwdump)
!*****************************************************************************
!     Compute a single numerical AO integral from the old (non-blocked) grid
!
!     Purpose:
!     For one integral label (format: AM + 2-char atom index + 2-char charge
!     + 2-char type), read grid points and AO values from DENSITY.GRID and
!     MBIS weights from WEIGHT.GRID point by point, accumulate the integral
!     sum_r chi_a(r)*chi_b(r)*w_A(r)*weight(r)*[dist(r)] into ao_int, and
!     write the result to AOPROPER. Supports monopole (00) and first-moment
!     (0x/0y/0z) integrals. See num_aoint_opt for the blocked-grid version.
!*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "priunit.h"
#include "rspprp.h"
#include "inflr.h"
! Used from include:
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       cord (atomic coordinates for distance calculation)
! inforb.h:       nbast, nnbasx
! priunit.h:      lupri
! rspprp.h:       MXLLR (needed for inflr.h)
! inflr.h:        iprlr
!
! Files used for:
! luddump:        reading grid coordinates, weights, density, and AO values
! luwdump:        reading per-point MBIS atom weights
! luprop:         writing the AO integral matrix via wrtpro
!
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

! label: AM + two char for index + two char for charge + two char for direc. and moment

      read(integral_label(3:4),*) atm_idx
      read(integral_label(5:6),*) atm_chg

      cor_idx = 0

      integral_type = integral_label(7:8)
      select case (integral_type)
         case ('00')
            cor_idx = 0
         case ('0x')
            cor_idx = 1
         case ('0y')
            cor_idx = 2
         case ('0z')
            cor_idx = 3
         case default
            write(lupri,*) 'Label', integral_type, 'not reognized'
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


      call getdat(rtnlbl(1),rtnlbl(2))

      rtnlbl(2) = 'SYMMETRI'
      
      call wrtpro(ao_int,nnbasx,integral_label,rtnlbl,iprlr)

      return
      end

      subroutine num_aoint_opt(ao_int, ao_tmp, gso_prod,
     &           mbis_weight,luwdump, cc_flag, nbuf2, natmlbl)
!*****************************************************************************
!     Optimized computation of all numerical AO integrals (blocked grid)
!
!     Purpose:
!     Parse all AM-type integral labels at once, then loop over blocks from
!     DENSITY.GRID using read_block to collect AO products and MBIS weights.
!     Accumulate all integrals simultaneously with a single DGEMM call per
!     block, avoiding repeated file I/O by processing all labels in one pass.
!     Write each resulting AO integral matrix to AOPROPER via wrtpro.
!*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "priunit.h"
#include "rspprp.h"
#include "inflr.h"
#include "ccrspprp.h"
#include "inftap.h"
#include "dftcom.h"
! Used from include:
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       cord (passed through to read_block for distance calc.)
! inforb.h:       nbast, nnbasx
! priunit.h:      lupri
! rspprp.h:       prplbl, nprlbl (non-CC property labels)
! inflr.h:        iprlr
! ccrspprp.h:     prplbl_cc, nprlbl_cc (CC property labels)
! inftap.h:       luddump
! dftcom.h:       mxbllen, ntot_dftgrid
!
! Files used for:
! luddump:        reading blocked grid data (AO products and coordinates)
! luwdump:        reading per-point MBIS atom weights (passed as argument)
! luprop:         writing AO integral matrices via wrtpro
!
      parameter (d0 = 0.0d0)
!incoming
      integer     luwdump,nbuf2,natmlbl
      real*8      ao_int(nnbasx,natmlbl),gso_prod(mxbllen*nnbasx),
     &            ao_tmp(nnbasx,natmlbl),
     &            mbis_weight(mxbllen*natmlbl)
      logical     cc_flag
!in here
      character*8 label, rtnlbl(2), prplbl_atm(natmlbl)
      integer     npoints, lbl_idx, block_size, points_read, footer, 
     &            b_points, b_gsoprod,
     &            atm_idx(natmlbl), atm_chg(natmlbl), cor_idx(natmlbl),
     &            i,j,k,n,idx,offset, gso_prod_idx(nnbasx)
      real*8      dist,scal,val, ao_test(nnbasx)

! parse all labels in here, so we can calculate all integrals with minimal file i/o
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
                  write(lupri,*) 'Label ', label, ' not recognized'
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
                  write(lupri,*) 'Label ', label, ' not recognized'
                  call quit("Label not recognized in numint_write")
               end if
            endif
         enddo
      end if

      ao_int = d0
      
      rewind(luwdump)
      read(luwdump) npoints
      rewind(luddump)

      points_read = 0
      block_size = 0

      ao_test = d0

      do while (points_read .lt. npoints)
         read(luddump) b_points, b_gsoprod 
         read(luddump) (gso_prod_idx(k),k=1,b_gsoprod)
         if (b_points .eq. 0) cycle
         call read_block(luwdump,gso_prod,mbis_weight,
     &                    points_read,b_gsoprod,b_points,
     &                    cor_idx,atm_idx,natmlbl)

         ao_tmp = d0
         call dgemm('N', 'T', b_gsoprod, natmlbl, b_points,
     &     1.0d0, gso_prod, b_gsoprod,!b_points,!mxbllen,
     &     mbis_weight, natmlbl, !b_points,
     &     1.0d0, ao_tmp, nnbasx)

         points_read = points_read + b_points

         do i=1,b_gsoprod
            idx = gso_prod_idx(i)
            do j=1,natmlbl
               ao_int(idx,j) = ao_int(idx,j) + ao_tmp(i,j)
            enddo
         enddo         
      enddo

      read(luwdump) footer
      if (footer .ne. -1) then
         write(lupri,*) 'Error in reading, footer is', footer
         call quit("Error reading file in num_aoint")
      end if

      call getdat(rtnlbl(1),rtnlbl(2))
      rtnlbl(2) = 'SYMMETRI'

      do k=1,natmlbl
         ao_test = ao_int(1:nnbasx,k)
         call wrtpro(ao_test,nnbasx,prplbl_atm(k),rtnlbl,iprlr)
      enddo

      return
      end

      subroutine read_block(luwdump,gso_prod,mbis_weight,points_read,
     &      b_gsoprod,b_points,cor_idx,atm_idx,natmlbl)
!*****************************************************************************
!     Read one block of grid data for num_aoint_opt
!
!     Purpose:
!     Read one DFT block from DENSITY.GRID (AO products and, if needed,
!     grid coordinates) and the corresponding MBIS weights from WEIGHT.GRID.
!     For integrals that require a coordinate distance (cor_idx > 0), read
!     the coordinates and multiply each MBIS weight by the atom-centered
!     displacement, so the caller (num_aoint_opt) needs no per-point geometry.
!*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "priunit.h"
#include "inftap.h"
#include "dftcom.h"
! Used from include:
! mxcent.h:       MXCENT, needed for nuclei.h
! nuclei.h:       nucdep, cord (atomic coordinates for distance calculation)
! inforb.h:       nnbasx
! priunit.h:      lupri
! inftap.h:       luddump
! dftcom.h:       mxbllen
!
! Files used for:
! luddump:        reading AO products, coordinates, weights, and density
! luwdump:        reading per-point MBIS atom weights (passed as argument)
!
! incoming
      integer  luwdump, cor_idx(natmlbl), atm_idx(natmlbl),
     &         natmlbl, points_read, b_gsoprod, b_points
      real*8   mbis_weight(mxbllen*natmlbl), 
     &         gso_prod(mxbllen*nnbasx)
! in here
      integer  i, j, k, npoint, gso_prod_idx(nnbasx), joint_idx
      real*8   tmp(b_points), tmp_weight(nucdep), cor(3,b_points)


      gso_prod = 0.0d0
      mbis_weight = 0.0d0
      if (maxval(cor_idx) .eq. 0) then
         read(luddump) tmp, tmp, tmp, tmp, tmp, 
     &            (gso_prod(i),i=1,b_gsoprod*b_points)
         if (nucdep .eq. natmlbl) then
            do i=1,b_points
               read(luwdump) (mbis_weight(j+(i-1)*natmlbl), j=1,nucdep)
            enddo
         else
            do i=1,b_points
               read(luwdump) (tmp_weight(j), j=1,nucdep)
               do k=1,natmlbl
                  mbis_weight(k+(i-1)*natmlbl) = tmp_weight(atm_idx(k))
               enddo
            enddo
         endif
      else
         read(luddump) cor(1,:), cor(2,:), cor(3,:), tmp, tmp, 
     &            (gso_prod(i),i=1,b_gsoprod*b_points)
         do i=1,b_points
            read(luwdump) (tmp_weight(j), j=1,nucdep)
            do k=1,natmlbl
               if (cor_idx(k) .eq. 0) then
                  mbis_weight(k+(i-1)*natmlbl) = tmp_weight(atm_idx(k))
               else
                  dist = cor(cor_idx(k),i) 
     &                     - cord(cor_idx(k),atm_idx(k))
                  mbis_weight(k+(i-1)*natmlbl) = 
     &                     tmp_weight(atm_idx(k)) * dist
               endif
            enddo
         enddo
      endif

      return
      end

