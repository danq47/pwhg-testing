C ----------------------------------------------------------------- C
C - This is a parton level only analysis for ttbar production     - C
C - tops bs Ws etc are constructed using MC truth only - no jets  - C
C - The MC truth reconstruction has been tested (see sanity check - C
C - code in the analysis below).                                  - C
C - Since it is parton level you need to comment out from         - C
C - CALL HWDHOB down to CALL HWDHOB inclusive in main-HERWIG.f .  - C
C - Also at some point, the showering went into what looked like  - C
C - an infinite loop after 137K events - gdb said it was in       - C
C - HWHGUP. The same glitch did not occur with *** herwig6520.f *** C
C - I also eliminated the analysis as a possible cause (it occurs - C
C - with HWANAL removed). I did not see anything fishy with the   - C
C - Tevatron, semileptonic event that got caught.                 - C
C ----------------------------------------------------------------- C

c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer j,l,m,l1,l2,l3,lenocc
      character * 24 prefix1,prefix2,prefix3
     	external lenocc

      call inihists

      do j = 1,3
      	if(j.eq.1) then
      		prefix1 = '-incl'
      	elseif(j.eq.2) then
      		prefix1 = '-str'
      	elseif(j.eq.3) then
      		prefix1 = '-unstr'
C          elseif(j.eq.4) then
C             prefix1 = '-vstr'
C          elseif(j.eq.5) then
C             prefix1 = '-vunstr'
C          else
C       		prefix1 = '-qqb'
      	endif
      	
      	do l=1,2
      		if(l.eq.1) then
      			prefix2 = '-no-cuts'
            else
               prefix2 = '-gg'
C       		elseif(l.eq.2) then
C       			prefix2 = '-wa'
C       		elseif(l.eq.3) then
C       			prefix2 = '-coll'
      		endif

         	l1=lenocc(prefix1)
         	l2=lenocc(prefix2)

            if(j.ne.6) then
      		   call bookupeqbins('pT-j1-2GeV'//prefix1(1:l1)//prefix2(1:l2),2d0,0d0,500d0)
      		   call bookupeqbins('pT-j1-5GeV'//prefix1(1:l1)//prefix2(1:l2),5d0,0d0,1000d0)
      		   call bookupeqbins('pT-j1-10GeV'//prefix1(1:l1)//prefix2(1:l2),10d0,0d0,1500d0)
      		   call bookupeqbins('pT-j1-50GeV'//prefix1(1:l1)//prefix2(1:l2),50d0,0d0,2000d0)
      		   call bookupeqbins('pT-j1-200GeV'//prefix1(1:l1)//prefix2(1:l2),200d0,0d0,3000d0)
            elseif(l.eq.1.and.j.eq.6) then
               call bookupeqbins('pT-j1-2GeV'//prefix1(1:l1)//prefix2(1:l2),2d0,0d0,500d0)
               call bookupeqbins('pT-j1-5GeV'//prefix1(1:l1)//prefix2(1:l2),5d0,0d0,1000d0)
               call bookupeqbins('pT-j1-10GeV'//prefix1(1:l1)//prefix2(1:l2),10d0,0d0,1500d0)
               call bookupeqbins('pT-j1-50GeV'//prefix1(1:l1)//prefix2(1:l2),50d0,0d0,2000d0)
               call bookupeqbins('pT-j1-200GeV'//prefix1(1:l1)//prefix2(1:l2),200d0,0d0,3000d0)
            endif

      	enddo
      enddo

c      call bookupeqbins('pT-j1-2GeV',2d0,0d0,500d0)
c      call bookupeqbins('pT-j1-5GeV',5d0,0d0,1000d0)
c      call bookupeqbins('pT-j1-10GeV',10d0,0d0,1500d0)
c      call bookupeqbins('pT-j1-50GeV',50d0,0d0,2000d0)
c      call bookupeqbins('pT-j1-200GeV',200d0,0d0,3000d0)


      end

      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_rad.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      integer jpref
      character * 20 prefix1,prefix2,prefix3
      common/ccccprefix/jpref
      data whcprg/'NLO   '/
      real * 8  dsig0,dsig
      logical   ini
      data      ini/.true./
      save      ini
      integer   ihep                ! HEPEVT index.
      real * 8 p_top(4),p_tb(4),p_wp(4),p_wm(4),p_lwp(4),p_lwm(4),
     1         p_nuwp(4),p_nuwm(4),p_b(4),p_bb(4),y,eta,pt,mass,
     2         ptzmf(4),plzmf(4),p_jet(4)
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer mjets,jetvec(maxtracks)
      logical   isForClustering(maxtracks)
      real * 8 j_kt(maxjets),j_eta(maxjets),j_rap(maxjets),
     1     j_phi(maxjets),j_p(4,maxjets),y_t,y_tbar,deltay,y_ttb,y_j1
      real * 8 deltaphi_j_t
      integer j,id,i_top,i_atop,i_jet,i_bfromtop,i_abfromatop,
     1     i_wp,i_wm,i_lwp,i_lwm,i_nuwp,i_nuwm,i_bjet,i_abjet,jhep,
     1     i_part,njets20,njets30,njets40
      real * 8 mtop,mtb,mwp,mwm,mb,mbb,p_bmax,e_bmax,xb,
     1     p_bbmax,e_bbmax,xbb,ewp,pwp,ewm,pwm,xw,
     2     dy,deta,dphi,dr,cth1,cth2,ptj1,mttbar,ptt,pttb,yj_minus_yttb
      integer jcth1,i_jets,ixx,jzz,qxx,jet_index,kxx,counter1,counter2,lxx,mxx,jxx
      integer jet_position(maxjets),l1,l2,l3
      real * 8 w(4),pb(4),ptb,eta_t,eta_tb,eta_ttb,eta_j1
      real * 8 prodvec2,powheginput
      real * 8 azi,deltaphi,low_bin_edge,pseudorapidity
      logical incl,stretched,unstretched,qqb
      logical nocuts,wa,wam,wab,wabm
      logical pt10,pt25,pt50,condition1,condition2,condition3
      logical sonofid
      external sonofid,pseudorapidity
      integer in_jet,lenocc
      external in_jet,azi,deltaphi,lenocc
      integer ngenerations,inotfound,iprodrad
      common/cngenerations/ngenerations
      integer id1,id2
      real*8 random,rweight
      external random

      ngenerations = powheginput("#ngenerations")
      if(ngenerations.lt.0) ngenerations = 4

      dsig  = dsig0

      i_top = 0
      i_atop = 0
      i_wp = 0
      i_wm = 0
      i_lwp = 0
      i_lwm = 0
      i_bfromtop = 0
      i_abfromatop = 0
      i_jet = 0

      do jhep=1,nhep
         id=idhep(jhep)
         if(idhep(jhep).eq.6) i_top = jhep
         if(idhep(jhep).eq.-6) i_atop = jhep
         id=abs(id)
         if(id.eq.5.or.id.eq.24) then
            if(sonofid(6,jhep)) then
               if(idhep(jhep).eq.5) i_bfromtop = jhep
               if(idhep(jhep).eq.-5) i_abfromatop = jhep
               if(idhep(jhep).eq.24) i_wp = jhep
               if(idhep(jhep).eq.-24) i_wm = jhep
            endif
         endif
         if(id.ge.11.and.id.le.14) then
            if(sonofid(24,jhep)) then
               if(idhep(jhep).eq.-11.or.idhep(jhep).eq.-13) i_lwp = jhep
               if(idhep(jhep).eq.11.or.idhep(jhep).eq.13) i_lwm = jhep
               if(idhep(jhep).eq.-12.or.idhep(jhep).eq.-14)i_nuwm = jhep
               if(idhep(jhep).eq.12.or.idhep(jhep).eq.14) i_nuwp = jhep
            endif
         endif


! Select them like this because sometimes the NLO ones have only Born process, and so the 5th entry in phep is junk
         if(whcprg.eq.'NLO') then
            if(jhep.gt.2) then ! If it's a final state parton
               if(abs(idhep(jhep)).lt.6.or.idhep(jhep).eq.21) then
                  i_jet = jhep
               endif
            endif
         elseif(whcprg.eq.'LHE') then
            i_jet = 5
         endif

c for jets, using only final state particles excluding leptons
         if(isthep(jhep).eq.1.and.
     1           (abs(idhep(jhep)).lt.11.or.abs(idhep(jhep)).gt.16)) then
            IsForClustering(jhep) = .true.
         else
            IsForClustering(jhep) = .false.
         endif
      enddo


      id1=idhep(1)
      id2=idhep(2)
      if(id1.eq.21) id1=0
      if(id2.eq.21) id2=0

      p_top=phep(1:4,i_top)
      p_tb=phep(1:4,i_atop)
      p_jet=phep(1:4,i_jet)

      mjets = maxjets
      call buildjets(mjets,j_kt,j_eta,j_rap,j_phi,j_p,jetvec,
     1     isForClustering)

      i_bjet = in_jet(i_bfromtop,jetvec)
      i_abjet = in_jet(i_abfromatop,jetvec)


      call getyetaptmass(p_top,y,eta,pt,mass)
      y_t=y
      call getyetaptmass(p_tb,y,eta,pt,mass)
      y_tbar=y
      call getyetaptmass(p_top+p_tb,y,eta,pt,mass)
      y_ttb=y
      mttbar=mass

      deltay=y_t-y_tbar

      call getyetaptmass(p_jet,y,eta,pt,mass)
  		ptj1 = pt
  		y_j1 = y

c Analysis - make the cuts

      if(whcprg.eq.'NLO') rho=rho_idx

      do jxx = 1,3

      	condition1 = .false.

      	if(jxx.eq.1) then
      		prefix1='-incl'
      		condition1 = .true.
            counter1=counter1+1
        	elseif(jxx.eq.2) then
      		prefix1='-str'
            if((rho.eq.1.and.deltay.lt.0).or.(rho.eq.2.and.deltay.gt.0)) then
      			condition1 = .true.
               counter2=counter2+1
      		endif
      	elseif(jxx.eq.3) then
      		prefix1='-unstr'
            if((rho.eq.1.and.deltay.gt.0).or.(rho.eq.2.and.deltay.lt.0)) then
      			condition1 = .true.
      		endif
C       	elseif(jxx.eq.4) then
C             prefix1='-vstr'
C             if((rho.eq.1.and.deltay.lt.-2).or.(rho.eq.2.and.deltay.gt.2)) then
C                condition1 = .true.
C             endif
C          elseif(jxx.eq.5) then
C             prefix1='-vunstr'
C             if((rho.eq.1.and.deltay.gt.-2).or.(rho.eq.2.and.deltay.lt.2)) then
C                condition1 = .true.
C             endif
C          elseif(jxx.eq.6) then
C       		prefix1='-qqb'
C       		if(rho.gt.2) then
C       			condition1 = .true.
C       		endif
      	endif

      	do lxx=1,2

      		condition2 = .false.

      		if(lxx.eq.1) then
      			prefix2 = '-no-cuts'
      			if(ptj1.gt.0) then 	! Making sure there is a jet
      				condition2 = .true.
      			endif
            elseif(lxx.eq.2) then
               prefix2 = '-gg'
               if(ptj1.gt.0) then
                  if(id1.eq.0.and.id2.eq.0) then
                     condition2 = .true.
                  endif
               endif
C       		elseif(lxx.eq.2) then
C       			prefix2 = '-wa'   ! These cuts don't work as the radiation which is generated could be generated by the old AND new sudakov in the veto algorithm
C       			if(ptj1.gt.0) then
C                   if(sw.eq.0) then
C          				condition2 = .true.
C                   endif
C       			endif
C       		elseif(lxx.eq.3) then
C       			prefix2 = '-coll'
C       			if(ptj1.gt.0) then
C                   if(sw.eq.1) then
C       					condition2 = .true.
C       				endif
C       			endif
      		endif

      		l1=lenocc(prefix1)
      		l2=lenocc(prefix2)

            if(jxx.ne.6) then
      		   if(condition1.and.condition2) then
      			   call filld('pT-j1-2GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
      			   call filld('pT-j1-5GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
      			   call filld('pT-j1-10GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
      			   call filld('pT-j1-50GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
      			   call filld('pT-j1-200GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
      		   endif
            elseif(lxx.eq.1.and.jxx.eq.6) then
               if(condition1.and.condition2) then
                  call filld('pT-j1-2GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
                  call filld('pT-j1-5GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
                  call filld('pT-j1-10GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
                  call filld('pT-j1-50GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
                  call filld('pT-j1-200GeV'//prefix1(1:l1)//prefix2(1:l2),ptj1,dsig)
               endif
            endif
      	enddo
      enddo

C       write(24,*) 'total events:', counter1
C       write(24,*) 'stretched events:', counter2
C       write(24,*) 'stretched percentage:', (100.0 * counter2)/(1.0*counter1)
!      counter1=counter1+1
C       if(powheginput('#dan_flag').eq.1d0) then  ! stretched
C C          rweight = min(ggbornplanar1,ggbornplanar2)/(ggbornplanar1+ggbornplanar2)
C C          if(random().lt.rweight) then
C C             counter2=counter2+1
C             if(id1.eq.0.and.id2.eq.0) then
C                call filld('pT-j1-2GeV',ptj1,dsig)
C                call filld('pT-j1-5GeV',ptj1,dsig)
C                call filld('pT-j1-10GeV',ptj1,dsig)
C                call filld('pT-j1-50GeV',ptj1,dsig)
C                call filld('pT-j1-200GeV',ptj1,dsig)
C             endif
C C          else
C C             if(id1.eq.0.and.id2.eq.0) then
C C                call filld('pT-j1-2GeV',0d0,dsig)
C C                call filld('pT-j1-5GeV',0d0,dsig)
C C                call filld('pT-j1-10GeV',0d0,dsig)
C C                call filld('pT-j1-50GeV',0d0,dsig)
C C                call filld('pT-j1-200GeV',0d0,dsig)
C C             endif
C C          endif
C       elseif(powheginput('#dan_flag').eq.2d0) then  ! unstretched
C C          rweight = max(ggbornplanar1,ggbornplanar2)/(ggbornplanar1+ggbornplanar2)
C C          if(random().lt.rweight) then
C C             counter2=counter2+1
C             if(id1.eq.0.and.id2.eq.0) then
C                call filld('pT-j1-2GeV',ptj1,dsig)
C                call filld('pT-j1-5GeV',ptj1,dsig)
C                call filld('pT-j1-10GeV',ptj1,dsig)
C                call filld('pT-j1-50GeV',ptj1,dsig)
C                call filld('pT-j1-200GeV',ptj1,dsig)
C             endif
C C          else
C C             if(id1.eq.0.and.id2.eq.0) then
C C                call filld('pT-j1-2GeV',0d0,dsig)
C C                call filld('pT-j1-5GeV',0d0,dsig)
C C                call filld('pT-j1-10GeV',0d0,dsig)
C C                call filld('pT-j1-50GeV',0d0,dsig)
C C                call filld('pT-j1-200GeV',0d0,dsig)
C C             endif
C C          endif
C       else
C          if(id1.eq.0.and.id2.eq.0) then
C             call filld('pT-j1-2GeV',ptj1,dsig)
C             call filld('pT-j1-5GeV',ptj1,dsig)
C             call filld('pT-j1-10GeV',ptj1,dsig)
C             call filld('pT-j1-50GeV',ptj1,dsig)
C             call filld('pT-j1-200GeV',ptj1,dsig)
C          endif
C       endif

      end


      function in_jet(i_part,jetvec)
      implicit none
      include 'hepevt.h'
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer in_jet,jetvec(maxtracks),i_part
      integer j
      logical sonofhep
      external sonofhep
      do j=1,nhep
         if(jetvec(j).ne.0) then
            if(sonofhep(i_part,j)) then
               in_jet = jetvec(j)
               return
            endif
         endif
      enddo
      in_jet = 0
      end

      function prodvec2(vec1,vec2)
      implicit none
      real * 8 prodvec2,vec1(4),vec2(4)
      prodvec2=vec1(4)*vec2(4)-vec1(1)*vec2(1)
     1 -vec1(2)*vec2(2)-vec1(3)*vec2(3)
      end

      function deltaphi(azi1,azi2)
      implicit none
      real * 8 pi
      parameter(pi = 3.141592653589793D0)
      real * 8 deltaphi,azi1,azi2
      deltaphi = abs(azi1-azi2)
      deltaphi = min(deltaphi,2*pi - deltaphi)
      end

      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end



      subroutine get_afb_hist(f_idx,b_idx,afb_idx)
      implicit none
      include 'pwhg_bookhist-new.h'
      real*8   f,ef,b,eb
      integer  ixx,f_idx,b_idx,afb_idx
      
      do ixx=1,nbins(f_idx)
         f=yhistarr2(ixx,f_idx)
         ef=errhistarr2(ixx,f_idx)
         b=yhistarr2(ixx,b_idx)
         eb=errhistarr2(ixx,b_idx)
         if((f+b).gt.0d0) then         ! Guard against division by zero.
            yhistarr2(ixx,afb_idx)=(f-b)/(f+b)
            errhistarr2(ixx,afb_idx)=2*sqrt((f*ef)**2+(b*eb)**2)
     1                              /(f+b)**2
         else
            yhistarr2(ixx,afb_idx)=0d0
            errhistarr2(ixx,afb_idx)=0d0
         endif
      enddo

      end


      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.16) then
         islept = .true.
      else
         islept = .false.
      endif
      end

      function phepDot(p_A,p_B)
      implicit none
      real * 8  phepDot
      real * 8  p_A(4),p_B(4)
      phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)
     1       -p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end

c     calculate the separation in the lego plot between the two momenta
c     p1 and p2 in azi and pseudorapidity
      function rsepn_p(p1,p2)
      implicit none
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 rsepn_p,p1(0:3),p2(0:3)
      real * 8 eta1,phi1,eta2,phi2
      real * 8 delphi
      real * 8 pseudorapidity,azi
      external pseudorapidity,azi

      phi1 = azi(p1)   
      phi2 = azi(p2)
      eta1 = pseudorapidity(p1)
      eta2 = pseudorapidity(p2)

      delphi = abs(phi1-phi2)
      if (delphi.gt.pi) then
         delphi = 2*pi-delphi
      endif
      if (delphi.lt.0 .or. delphi.gt.pi) then
         print*,' problem in rsepn. delphi = ',delphi
      endif
      rsepn_p = sqrt( (eta1-eta2)**2 + delphi**2 )
      end

      function azi(p)
      implicit none
      real * 8 pi,pi2
      parameter(pi = 3.141592653589793D0, pi2 = 9.869604401089358D0)
      real * 8 azi,p(0:3)
      azi = atan(p(1)/p(0))
      if (p(0).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end


      function pseudorapidity(p)
      implicit none
      real * 8 p(0:3),pseudorapidity
      real * 8 mod, costh
      mod = sqrt(p(1)**2+p(2)**2+p(3)**2)
      costh = p(3)/mod
      pseudorapidity=0.5*log((1+costh)/(1-costh))
      end


      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,jetvechep,
     1                                               isForClustering)
c     arrays to reconstruct jets
      implicit  none
      include  'hepevt.h'
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer   mjets,jetvechep(maxtracks)
      real * 8  kt(maxjets),eta(maxjets),rap(maxjets),
     1     phi(maxjets),pjet(4,maxjets)
      logical   isForClustering(maxtracks)
      real * 8  ptrack(4,maxtracks),pj(4,maxjets)
      integer   jetvec(maxtracks),itrackhep(maxtracks)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8  r,palg,ptmin,pp,tmp
      logical sonofid
      external sonofid
C - Initialize arrays and counters for output jets
      ptrack = 0
      jetvec = 0
      ntracks=0
      pjet = 0
      pj = 0
      kt = 0
      eta = 0
      rap = 0
      phi = 0
C - Extract final state particles to feed to jet finder
      do j=1,nhep
         if(.not.isForClustering(j)) cycle
         if(ntracks.eq.maxtracks) then
            write(*,*) 'analyze: need to increase maxtracks!'
            write(*,*) 'ntracks: ',ntracks
            call exit(-1)
         endif
         ntracks=ntracks+1
         ptrack(:,ntracks) = phep(1:4,j)
         itrackhep(ntracks)=j
      enddo
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------- C
C - Run FastJet - C
C --------------- C
C - R = 0.7   radius parameter
C - f = 0.75  overlapping fraction
      palg  = -1
      r     = 0.5d0
      ptmin = 0.1d0
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,jetvec)
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pp+pjet(3,j))/(pp-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      jetvechep = 0
      do j=1,ntracks
         jetvechep(itrackhep(j))=jetvec(j)
      enddo
      end

      function sonofid(m,k)
      implicit none
      logical sonofid
      integer m,k
      include  'hepevt.h'
      integer j,kcurr
      integer ngenerations
      common/cngenerations/ngenerations
      kcurr=k
      do j=1,ngenerations
         if(abs(idhep(kcurr)).eq.m) then
            sonofid = .true.
            return
         endif
         kcurr = jmohep(1,kcurr)
         if(kcurr.eq.0) then
            sonofid = .false.
            return
         endif
      enddo
      sonofid=.false.
      end


      function sonofhep(m,k)
      implicit none
      logical sonofhep
      integer m,k
      include  'hepevt.h'
      integer j,kcurr
      integer ngenerations
      common/cngenerations/ngenerations
      kcurr=k
      do j=1,ngenerations
         if(kcurr.eq.m) then
            sonofhep = .true.
            return
         endif
         kcurr = jmohep(1,kcurr)
         if(kcurr.eq.0) then
            sonofhep = .false.
            return
         endif
      enddo
      sonofhep = .false.
      end


