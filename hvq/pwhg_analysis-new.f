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
      integer j,l
      character * 20 prefix
      integer nbins
      parameter (nbins=11)
      real * 8 pT_tt_bins(nbins+1)
      data pT_tt_bins/  0d0, 10d0, 25d0, 50d0,100d0,
     1     150d0,200d0,250d0,300d0,400d0,600d0, 900d0/          
      real * 8 m_tt_bins(nbins+1)
      data m_tt_bins/ 320d0,360d0,400d0,450d0,500d0,
     1     550d0,600d0,650d0,700d0,800d0,900d0,1000d0/          
      character * 2 digit(20)
      data digit/'01','02','03','04','05','06','07','08','09','10',
     1           '11','12','13','14','15','16','17','18','19','20'/
      integer lenocc
      external lenocc

      call inihists

c number of jets with a pT greater than X
      call bookupeqbins('Njets-pt10'  ,1d0,-0.5d0,10.5d0)
      call bookupeqbins('Njets-pt20'  ,1d0,-0.5d0,10.5d0)
      call bookupeqbins('Njets-pt30'  ,1d0,-0.5d0,10.5d0)
      call bookupeqbins('Njets-pt40'  ,1d0,-0.5d0,10.5d0)

c pT of the hardest, second hardest, etc jet
      call bookupeqbins('pT-j1',1d0,10d0,300d0)
      call bookupeqbins('pT-first-emission',1d0,0d0,300d0)
      call bookupeqbins('pT-j2',1d0,10d0,300d0)
      call bookupeqbins('pT-j3',1d0,10d0,300d0)

c observables for the top quark      
      call bookupeqbins('pT-t',10d0,0d0,300d0)
      call bookupeqbins('mass-t',10d0,0d0,300d0)
      call bookupeqbins('mass-zoom-t',1d-1,160d0,180d0)
      call bookupeqbins('rapidity-t',1d-1,-4d0,4d0)

c observables for the antitop
      call bookupeqbins('pT-tbar',10d0,0d0,300d0)
      call bookupeqbins('mass-tbar',10d0,0d0,300d0)
      call bookupeqbins('mass-zoom-tbar',1d-1,160d0,180d0)
      call bookupeqbins('rapidity-tbar',1d-1,-4d0,4d0)

c observables for the ttbar system
      call bookupeqbins('pT-ttbar-all-mttbar',1d0,0d0,300d0)
      call bookupeqbins('invariant-mass-t-tbar',10d0,0d0,1000d0)
      call bookupeqbins('rapidity-ttbar',1d-1,-4d0,4d0)

c rapitdity asymmetry
      call bookupeqbins('yt-ytb',1d-1,-4d0,4d0)

c rapidity of the jets
      call bookupeqbins('y-j1',1d-1,-4d0,4d0)
      call bookupeqbins('y-j2',1d-1,-4d0,4d0)

c hardest, 2nd hardest jet rapidity in the ttbar-rapidity=0 frame
      call bookupeqbins('rapidity-j1-minus-rapidity-ttbar',1d-1,-4d0,4d0)
      call bookupeqbins('rapidity-j2-minus-rapidity-ttbar',1d-1,-4d0,4d0)

c azimuthal angles. phi (the function azi(p)) is defined between -pi and pi
c       call bookupeqbins('phi-t',1d-1,-3.2d0,3.2d0)
c       call bookupeqbins('phi-tbar',1d-1,-3.2d0,3.2d0)
c       call bookupeqbins('phi-j1',1d-1,-3.2d0,3.2d0)

c azimuthal angle of the ttbar system i.e. the top and the anti top considered together (not the same as below - think about it)
c without showering, this should be equivalent to phi-j1
c       call bookupeqbins('phi-t-tbar',1d-1,-3.2d0,3.2d0)

c angular separation of the top and the anti-top
      call bookupeqbins('delta-phi-t-tbar',1d-1,0d0,3.2d0)

      end

      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      integer jpref
      character * 20 prefix(18)
      common/ccccprefix/jpref,prefix
      data whcprg/'NLO   '/
      real * 8  dsig0,dsig
      logical   ini
      data      ini/.true./
      save      ini
      integer   ihep                ! HEPEVT index.
      real * 8 p_top(4),p_tb(4),p_wp(4),p_wm(4),p_lwp(4),p_lwm(4),
     1         p_nuwp(4),p_nuwm(4),p_b(4),p_bb(4),y,eta,pt,mass,
     2         ptzmf(4),plzmf(4)
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer mjets,jetvec(maxtracks)
      logical   isForClustering(maxtracks)
      real * 8 j_kt(maxjets),j_eta(maxjets),j_rap(maxjets),
     1     j_phi(maxjets),j_p(4,maxjets),y_t,y_tbar
      integer j,id,i_top,i_atop,i_bfromtop,i_abfromatop,
     1     i_wp,i_wm,i_lwp,i_lwm,i_nuwp,i_nuwm,i_bjet,i_abjet,jhep,
     1     i_part,njets20,njets30,njets40
      real * 8 mtop,mtb,mwp,mwm,mb,mbb,p_bmax,e_bmax,xb,
     1     p_bbmax,e_bbmax,xbb,ewp,pwp,ewm,pwm,xw,
     2     dy,deta,dphi,dr,cth1,cth2
      integer jcth1,i_jets,ixx,jzz,qxx,jet_index
      integer jet_position(maxjets)
      real * 8 w(4),pb(4),ptb
      real * 8 prodvec2,powheginput
      real * 8 azi,deltaphi
      logical sonofid
      external sonofid
      integer in_jet
      external in_jet,azi,deltaphi
      integer ngenerations,inotfound,iprodrad
      common/cngenerations/ngenerations
      character * 2 digit(20)
      data digit/'01','02','03','04','05','06','07','08','09','10',
     1           '11','12','13','14','15','16','17','18','19','20'/
      integer id1,id2
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

      if(whcprg.eq.'NLO') then
         i_top = 3
         i_atop = 4
         i_wp = 5
         i_wm = 6

         i_lwp = 7
         i_lwm = 9
         i_nuwp = 8
         i_nuwm = 10
         i_bfromtop = 11
         i_abfromatop = 12

         IsForClustering = .false.
         IsForClustering(13) = .true.
         IsForClustering(i_bfromtop) = .true.
         IsForClustering(i_abfromatop) = .true.
C --------------------------------------------- C
C - LHE PARTICLE TOP RECONSTRUCTION: MC TRUTH - C
C --------------------------------------------- C
      else
c Build top MC; find the last top (tbar)
c in the event record, i.e. before decays
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
c for jets, using only final state particles excluding leptons
            if(isthep(jhep).eq.1.and.
     1           (abs(idhep(jhep)).lt.11.or.abs(idhep(jhep)).gt.16)) then
               IsForClustering(jhep) = .true.
            else
               IsForClustering(jhep) = .false.
            endif
         enddo
      endif

      inotfound = 0
      if(i_top.eq.0) then
         write(*,*) 'top not found'
         inotfound = inotfound + 1
      endif
      if(i_atop.eq.0) then
         write(*,*) 'antitop not found'
         inotfound = inotfound + 1
      endif
      if(i_wp.eq.0) then
         write(*,*) 'wp not found'
         inotfound = inotfound + 1
      endif
      if(i_wm.eq.0) then
         write(*,*) 'wm not found'
         inotfound = inotfound + 1
      endif
      if(i_lwp.eq.0) then
         write(*,*) 'lwp not found'
         inotfound = inotfound + 1
      endif
      if(i_lwm.eq.0) then
         write(*,*) 'lwm not found'
         inotfound = inotfound + 1
      endif
      if(i_nuwp.eq.0) then
         write(*,*) 'nuwp not found'
         inotfound = inotfound + 1
      endif
      if(i_nuwm.eq.0) then
         write(*,*) 'nuwm not found'
         inotfound = inotfound + 1
      endif
      if(i_bfromtop.eq.0) then
         write(*,*) 'b from top not found'
         inotfound = inotfound + 1
      endif
      if(i_abfromatop.eq.0) then
         write(*,*) 'bbar from tbar not found'
         inotfound = inotfound + 1
      endif

      if(inotfound.gt.0) return

c 
      if(whcprg.ne.'NLO'.and.whcprg.ne.'LHE') then
c Setup a flag:
c 1 for events with no production radiation in the LHE
c 2 for events with production radiation
c 0 otherwise
         if((phep(1,3)+phep(1,4))**2+(phep(2,3)+phep(2,4))**2.lt.1d-2)
     1        then
            iprodrad=1
c            if(scalup.gt.46) then
c               write(*,*) ' warning: scalup = ',scalup,
c     1              ' in radiation in decay'
c            endif
         else
            iprodrad=2
         endif
      else
         iprodrad = 0
      endif

      id1=idhep(1)
      id2=idhep(2)
      if(id1.eq.21) id1=0
      if(id2.eq.21) id2=0

      if(whcprg.eq.'LHE') then
         if(powheginput('#subprocess').eq.1) then
            if(id1.eq.0.and.id2.eq.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.2) then
            if(id1.gt.0.and.id2.lt.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.3) then
            if(id1.lt.0.and.id2.gt.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.4) then
            if(id1.gt.0.and.id2.eq.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.5) then
            if(id1.lt.0.and.id2.eq.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.6) then
            if(id1.eq.0.and.id2.gt.0) then
               continue
            else
               return
            endif
         elseif(powheginput('#subprocess').eq.7) then
            if(id1.eq.0.and.id2.lt.0) then
               continue
            else
               return
            endif
         endif
      endif
               

      p_top=phep(1:4,i_top)
      p_tb=phep(1:4,i_atop)
      p_wp=phep(1:4,i_wp)
      p_wm=phep(1:4,i_wm)
      p_lwp=phep(1:4,i_lwp)
      p_lwm=phep(1:4,i_lwm)
      p_nuwp=phep(1:4,i_nuwp)
      p_nuwm=phep(1:4,i_nuwm)
      p_b=phep(1:4,i_bfromtop)
      p_bb=phep(1:4,i_abfromatop)

      mjets = maxjets
      call buildjets(mjets,j_kt,j_eta,j_rap,j_phi,j_p,jetvec,
     1     isForClustering)

      i_bjet = in_jet(i_bfromtop,jetvec)
      i_abjet = in_jet(i_abfromatop,jetvec)

c We now have mjets number of jets, however, some of these are b-jets.
c For our analysis, we sometimes want to only see the non b-jets, so this piece of code
c picks out the position of the 3 hardest non-bjets (at most - sometimes less than this if the jets are too soft to be picked up)
c This runs into a problem when we do hadronisation, as the b hadrons are no longer picked up properly, however we can just turn this off for the moment


c jet_postition starts as an array of 1..maxjets, then any b jets are set to zero, then we shift the following entries up
c So if we have 4 jets, and a b-jet as the second hardest jet then jet_position will go
c step 1 -> [1,2,3,4]
c step 2 -> [1,0,3,4]
c step 3 -> [1,3,4]
c Therefore, the 3 hardest non-b jets are located at jets numbers 1,3, and 4
c which can be easily accesed as jet_position(1),jet_position(2), and jet_position_(3)
c so now, the rapidity of say the 2nd hardest non b jet is given by j_rap(jet_position(2))

c first set the position of the b jets to zero
      do i_jets=1,maxjets
         jet_position(i_jets)=i_jets
         if(i_jets.eq.i_bjet.or.i_jets.eq.i_abjet) then
            jet_position(i_jets)=0
         endif
      enddo
c Next, if there are no zeroes (i.e. neither of the b-jets are the hardest jets) then we simply take the first n jets
      if(i_bjet.eq.0.and.i_abjet.eq.0) then
         continue
c if there is only one zero, we only need to shift the indices in the array once
      elseif(i_bjet.eq.0.and.i_abjet.ne.0) then
         jet_index=1
         do jzz=1,maxjets
            if(jet_position(jzz).eq.0) then
               jet_index=jzz
            endif
         enddo
       
         do jzz=jet_index,maxjets-1
            jet_position(jzz)=jet_position(jzz+1)
         enddo
      elseif(i_bjet.ne.0.and.i_abjet.eq.0) then
         jet_index=1
         do jzz=1,maxjets
            if(jet_position(jzz).eq.0) then
               jet_index=jzz
            endif
         enddo
       
         do jzz=jet_index,maxjets-1
            jet_position(jzz)=jet_position(jzz+1)
         enddo
      else
c now the (ususal) case where both the b jets are in the in the list of hardest jets
c We delete the two zeroes
         do qxx=1,2
            jet_index=1
            do jzz=1,maxjets+1-qxx
               if(jet_position(jzz).eq.0) then
                  jet_index=jzz
               endif
            enddo
       
            do jzz=jet_index,maxjets-qxx
               jet_position(jzz)=jet_position(jzz+1)
            enddo
         enddo
      endif

      call getyetaptmass(p_top,y,eta,pt,mass)
      y_t=y
      call getyetaptmass(p_tb,y,eta,pt,mass)
      y_tbar=y

      njets20 = 0
      njets30 = 0
      njets40 = 0
      do j=1,mjets
         if(j_kt(j).gt.20) then
            njets20 = njets20 + 1
         endif
         if(j_kt(j).gt.30) then
            njets30 = njets30 + 1
         endif
         if(j_kt(j).gt.40) then
            njets40 = njets40 + 1
         endif
      enddo

      call filld('Njets-pt10',dble(mjets),dsig)
      call filld('Njets-pt20',dble(njets20),dsig)
      call filld('Njets-pt30',dble(njets30),dsig)
      call filld('Njets-pt40',dble(njets40),dsig)

c transverse momentum of the n-th hardest jet
      call filld('pT-j1',j_kt(jet_position(1)),dsig)
      call filld('pT-j2',j_kt(jet_position(2)),dsig)
      call filld('pT-j3',j_kt(jet_position(3)),dsig)
c this below is the pT of the hardest emission as generated by POWHEG (ie not to do with the jet clusetering)
c Therefore it only gives meaningful results if we don't shower the les houches file
      call filld('pT-first-emission',sqrt(phep(1,5)**2 + phep(2,5)**2),dsig)

c observables for the top quark
      call getyetaptmass(p_top,y,eta,pt,mass)
      call filld('pT-t',pt,dsig)
      call filld('mass-t',mass,dsig)
      call filld('mass-zoom-t',mass,dsig)
      call filld('rapidity-t',y,dsig)

c observables for the anti-top
      call getyetaptmass(p_tb,y,eta,pt,mass)
      call filld('pT-tbar',pt,dsig)
      call filld('mass-tbar',mass,dsig)
      call filld('mass-zoom-tbar',mass,dsig)
      call filld('rapidity-tbar',y,dsig)

c observables for the ttbar system
      call getyetaptmass(p_top+p_tb,y,eta,pt,mass)
      call filld('pT-ttbar-all-mttbar',pt,dsig)
      call filld('invariant-mass-t-tbar',mass,dsig)
      call filld('rapidity-ttbar',y,dsig)

c Rapidity asymmetry
      call filld('yt-ytb',y_t-y_tbar,dsig)

c Rapidity of the jets
      if(j_kt(jet_position(1)).ne.0) then
         call filld('y-j1',j_rap(jet_position(1)),dsig)
         call filld('rapidity-j1-minus-rapidity-ttbar',j_rap(jet_position(1))-y,dsig)
c          call filld('phi-j1',j_phi(jet_position(1)),dsig)
      endif
      if(j_kt(jet_position(2)).ne.0) then
         call filld('y-j2',j_rap(jet_position(2)),dsig)
         call filld('rapidity-j2-minus-rapidity-ttbar',j_rap(jet_position(2))-y,dsig)
      endif

c Angular plots
c       call filld('phi-t',azi(p_top),dsig)
c       call filld('phi-tbar',azi(p_tb),dsig)
c       call filld('phi-t-tbar',azi(p_top+p_tb),dsig) ! should be equivalent to phi_j1 for no showering etc, as the ttbar system only recoils off the hardest jet
      call filld('delta-phi-t-tbar',deltaphi(azi(p_top),azi(p_tb)),dsig)

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


      subroutine yetaptmassplot(p,dsig,prefix)
      implicit none
      real * 8 p(4),dsig
      character *(*) prefix
      real * 8 y,eta,pt,m
      call getyetaptmass(p,y,eta,pt,m)
      call filld(prefix//'_y',y,dsig)
      call filld(prefix//'_eta',eta,dsig)
      call filld(prefix//'_pt',pt,dsig)
      call filld(prefix//'_m',m,dsig)
      end


      subroutine pwhgfinalopshist
      implicit none
      include 'pwhg_bookhist-new.h'
      integer  f_idx,b_idx,x_idx,ixx,jxx,j,k,l
      integer  indexhist
      real * 8 num,numerr,den,denerr

      f_idx=indexhist('dF-dpT')
      b_idx=indexhist('dB-dpT')
      x_idx=indexhist('dAFB-dpT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dpT')
      b_idx=indexhist('dB-dpT')
      x_idx=indexhist('F-gt-pT')
      call integratehist(f_idx,x_idx,1)
      x_idx=indexhist('B-gt-pT')
      call integratehist(b_idx,x_idx,1)
      x_idx=indexhist('F-lt-pT')
      call integratehist(f_idx,x_idx,-1)
      x_idx=indexhist('B-lt-pT')
      call integratehist(b_idx,x_idx,-1)

      f_idx=indexhist('F-gt-pT')
      b_idx=indexhist('B-gt-pT')
      x_idx=indexhist('AFB-gt-pT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('F-lt-pT')
      b_idx=indexhist('B-lt-pT')
      x_idx=indexhist('AFB-lt-pT')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dm_tt')
      b_idx=indexhist('dB-dm_tt')
      x_idx=indexhist('dAFB-dm_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('dF-dm_tt')
      b_idx=indexhist('dB-dm_tt')
      x_idx=indexhist('F-gt-m_tt')
      call integratehist(f_idx,x_idx,1)
      x_idx=indexhist('B-gt-m_tt')
      call integratehist(b_idx,x_idx,1)
      x_idx=indexhist('F-lt-m_tt')
      call integratehist(f_idx,x_idx,-1)
      x_idx=indexhist('B-lt-m_tt')
      call integratehist(b_idx,x_idx,-1)

      f_idx=indexhist('F-gt-m_tt')
      b_idx=indexhist('B-gt-m_tt')
      x_idx=indexhist('AFB-gt-m_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      f_idx=indexhist('F-lt-m_tt')
      b_idx=indexhist('B-lt-m_tt')
      x_idx=indexhist('AFB-lt-m_tt')
      call get_afb_hist(f_idx,b_idx,x_idx)

      end


      subroutine integratehist(integrand,integral,direction)
      implicit none
      include 'pwhg_bookhist-new.h'
      integer  integrand,integral,direction
      integer  indexhist,ixx,jxx
      external indexhist
      if(direction.eq.1) then
        do ixx=1,nbins(integrand)
          do jxx=ixx,nbins(integrand)
            yhistarr2(ixx,integral)=yhistarr2(ixx,integral)
     1         +yhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand))
            errhistarr2(ixx,integral)=errhistarr2(ixx,integral)
     1         +(errhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand)))**2
          enddo
          errhistarr2(ixx,integral)=sqrt(errhistarr2(ixx,integral))
        enddo
      else if(direction.eq.-1) then
        do ixx=1,nbins(integrand)
          do jxx=ixx,1,-1
            yhistarr2(ixx,integral)=yhistarr2(ixx,integral)
     1         +yhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand))
            errhistarr2(ixx,integral)=errhistarr2(ixx,integral)
     1         +(errhistarr2(jxx,integrand)
     1         *(xhistarr(jxx+1,integrand)-xhistarr(jxx,integrand)))**2
          enddo
          errhistarr2(ixx,integral)=sqrt(errhistarr2(ixx,integral))
        enddo
      else
         write(*,*) 'subroutine integratehist: error!'
         write(*,*) '--------------------------------'
         write(*,*) 'Invalid direction specified for histogram given.'
         write(*,*) 'direction=1/-1 only - integral from lowest bin to'
         write(*,*) 'highest bin or highest bin to lowest bin.'
         call exit(-1)
      endif

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
      ptmin = 10d0
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

      subroutine boost2reson4(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(4),pin(4,nm),pout(4,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(4)
      vec(1)=pres(1)/(beta*pres(4))
      vec(2)=pres(2)/(beta*pres(4))
      vec(3)=pres(3)/(beta*pres(4))
      call mboost4(nm,vec,-beta,pin,pout)
      end


      
      subroutine mboost4(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(4,m),vout(4,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
      enddo
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end

      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end
