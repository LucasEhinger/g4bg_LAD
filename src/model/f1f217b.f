CCC   F1F217 version 0.95b  -  Oct. 24, 2018                              CCC
CCC   Collection of subroutines to calculate inclusive cross sections    CCC
CCC   for range of nuclei.  For A > 2 the parameterization is based on   CCC
CCC   by M. E. Christy, T. Gautam, and A Bodek to 12C, 27Al, 56Fe and    CCC
CCC   64Cu.  However, the fit scales relatively well with 'A' and        CCC
CCC   should be good for all nuclei with 10 < A < 80.                    CCC
CCC   Also included is the proton cross section fit and a preliminary    CCC
CCC   deuteron/neutron fit by M. E. Christy, N. Kalantarians, J. Either  CCC
CCC   and W. Melnitchouk (to be published) based on both inclusive       CCC
CCC   deuteron and tagged deuteron data on n/d from BONuS.               CCC
CCC   To get reduced cross sections or structure functions just call     CCC
CCC

      SUBROUTINE RESCSD(w2,q2,eps,doqe,f1d,f2d,fLd,wfn,sigm)
CCCC  Calcultes deuteron cross sections from smeared p+n.             CCCC
CCCC  Requires SQESUB be called first to read in smearing functions.  CCCC
CCCC  This should be one at the lowest level to keep from reading     CCCC
CCCC  multiple times.

      IMPLICIT none

      real*8 w2,q2,eps,t1,x,gamma2
      real*8 sigtp,sigtn,siglp,sigln,sigt,sigl,sigm,m  
      real*8 m2,pi2,alpha,f1d,f2d,fLd,f1dqe,f2dqe,fLdqe
      real*8 xvaln(100),off_mKP_fit,delta,xfr,xfl
      integer i,j,k,ntssot,wfn,drn
      logical doqe, first/.false./ 
c      INCLUDE 'parm.cmn'
      external off_mKP_fit

      drn = 5
      xfr = 0.95
      xfl = 1.E-3
      alpha = 1./137.036
      m = (0.938272+0.939565)/2.0d0  !!! average p,n
      m2 = m*m
      pi2 = 3.14159*3.14159
      call SQESUB(w2,q2,wfn,f2dqe,f1dqe,fLdqe,first) 

      x = q2/(w2-m2+q2)
      gamma2 = 1.+4.*m2*x*x/q2
      t1 = gamma2*f2dqe-2.*x*f1dqe

      if(f2dqe.LT.0.0.OR.f1dqe.LT.0.0.OR.fLdqe.LT.0.0)
     &    write(34,*) w2,q2,f2dqe,f1dqe,fLdqe
        

c      write(6,*) "In rescsd:  ",q2,x,gamma2,w2,doqe

      call SMEARSUB(w2,q2,wfn,f2d,f1d,fLd)
      delta = off_mKP_fit(x,wfn,drn)
      if(x.GT.xfr) delta =  off_mKP_fit(xfr,wfn,drn)
      if(x.LT.xfl) delta =  off_mKP_fit(xfl,wfn,drn)
      if(q2.LT.0.01) delta = 0.0
c      delta = 0.0D0
      f1d = f1d/(1.0D0-delta)
      f2d = f2d/(1.0D0-delta)
      fLd = fLd/(1.0D0-delta)
      if(doqe) then 
       f1d = f1d+f1dqe
       f2d = f2d+f2dqe
       fLd = fLd+fLdqe
      endif

      sigt = 0.3894e3*f1d*pi2*alpha*8.0/abs(w2-m2)
      sigl = 0.3894e3*fLd*pi2*alpha*8.0/abs(w2-m2)/2.*abs(w2-m2+q2)/q2
      sigm = sigt+eps*sigl

c      sigm = eps*sigl

c      write(6,*) "In rescsd:  ",q2,x,gamma2,w2,sigm

      return
      end



      SUBROUTINE SMEARSUB(w2,q2,wfn,f2s,f1s,fLs)
CCCC   Fermi smears structure functions.  Requires subroutne GETFY  CCCC
CCCC   to be initialized first.                                     CCCC

      IMPLICIT none

      real*8 w2,q2,mp,mp2,x,z,wwz,alpha,alpha2,gamma,gamma2
      real*8 y,fy11,fy12,fy2,inc,f2s,f1s,fLs,f2pi,f2ni,xvaln(100)
      real*8 sigtp,sigtn,siglp,sigln,f1pi,f1ni,fLpi,fLni
      real*8 f1i(1000),f2i(1000),vfact1(1000),vfact2(1000)
      real*8 ymin,ymax,pi,pi2,epsd
      integer i,j,nbins,wfn
c      logical first/.true./
      logical firsty/.false./
      real*8 dsimps
      external dsimps

      mp = (0.938272+0.939565)/2.0d0  !!! average of p, n
      mp2 = mp*mp
      alpha = 1/137.03599
      alpha2 = alpha*alpha 
      pi = 3.14159
      pi2 = pi*pi 
      epsd = -2.2E-03
      x = q2/(w2-mp2+q2)

      f2s = 0.0
      f1s = 0.0
      fLs = 0.0
      if(x.LT.0.0D0) return

      gamma2 = (1.+4.*mp2*x*x/q2)   
      gamma = sqrt(gamma2)
      nbins = 60

      if(w2.LT.1.6) nbins = 100

c      nbins = 200   !!! adjust as needed

      if(w2.GE.2.5.AND.Q2.GE.2.5) nbins = 40

CCC///   Fill array for Simpson's rule calculation   ///CCC
       

c      if(x.LT.1.0) then
c        ymin = max(x*(1.0D0-2.0D0*epsd*mp/(w2-mp2)),x)
c      else
c        ymin = x
c      endif 
CCC   Note problem if Q2 = 0  CCC       
      ymin = max(x*(1.0D0-2.0D0*epsd*mp/q2),x)
c      ymin = min(ymin,2.0d0)
      ymax = min(1.0d0+gamma2/2.0d0+epsd/mp,2.0d0)

c      write(6,*) w2,ymin,ymax 

      if(ymin.GE.ymax) return

c      ymax = 2.0d0      

c      write(6,*) ymin,ymax

      if(ymin.EQ.0) write(6,*) "Error ymin = 0" 
      inc = (ymax-ymin)/float(nbins)
      y = ymin-inc

      do i=1,nbins+1      !!!!    First calculate the cross smearing function for each y bin    !!!!
        y = y + inc
        z = x/y

        wwz = mp2 + q2*(1.d0/z - 1.d0)   !!! W^2 at x=z  !!!

        call rescsp(wwz,q2,sigtp,sigLp) 
        call rescsn(wwz,q2,sigtn,sigLn)
  
        f1pi = 2.*z*sigtp/0.3894e3/pi2/alpha/8.0*abs(wwz-mp2)  !!! 2xF1p
        fLpi = 2*q2/abs(wwz-mp2+q2)*sigLp/0.3894e3/pi2/alpha/8.0*
     &        abs(wwz-mp2)                                     
        f1ni = 2.*z*sigtn/0.3894e3/pi2/alpha/8.0*abs(wwz-mp2)  !!! 2xF1n
        fLni =  2*q2/abs(wwz-mp2+q2)*sigLn/0.3894e3/pi2/alpha/8.0*
     &        abs(wwz-mp2)  
        f2pi = (f1pi+fLpi)/(1.+4.*mp2*z*z/q2)
        f2ni = (f1ni+fLni)/(1.+4.*mp2*z*z/q2)
        f2i(i) = f2pi+f2ni  !!!!  Proton + Neutron  !!!
        f1i(i) = f1pi+f1ni  !!!!  Actually 2xF1     !!!

        call getfy(gamma,y,wfn,fy11,fy12,fy2,firsty)
        vfact1(i) = fy11*f1i(i)+fy12*f2i(i)
        vfact2(i) = fy2*f2i(i)
        if(vfact1(i).LT.0.0) vfact1(i) = 0.0
        if(vfact2(i).LT.0.0) vfact2(i) = 0.0

      enddo

      vfact1(1) = 0.0d0
      vfact2(1) = 0.0d0


CCC///    Integrate over y    ///CCC

      f1s = dsimps(vfact1,ymin,ymax,nbins)/2./x  !!! F1d
      f2s = dsimps(vfact2,ymin,ymax,nbins)       !!! F2d
      fLs = gamma2*f2s-2.*x*f1s                  !!! FLd

c      write(6,*) "in smearsub 2:  ",x,ymin,q2,gamma,f1s,f2s,fLs

 2000 format(7f12.4)

      end





      SUBROUTINE RESCSP(w2,q2,sigtp,siglp)
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis

      data xvalp / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /

      do i=1,50
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodp(1,w2,q2,xval1,sigTp)
      call resmodp(2,w2,q2,xvalL,sigLp)
c      call disp(w2,q2,sigTdis,sigLdis)
      if(w2.GT.8.0) then 
c        write(6,*) "resmod: ",w2,q2,sigTp,sigLp
        call disp(w2,q2,sigTdis,sigLdis)
c        write(6,*) "dismod: ",w2,q2,sigtdis,sigLdis
        if(w2.LE.10.0) then
         sigTp = (10.0-w2)*sigTp+(w2-8.0)*sigTdis
         sigLp = (10.0-w2)*sigLp+(w2-8.0)*sigLdis
         sigTp = sigTp/2.0
         sigLp = sigLp/2.0
         else
          sigTp = sigTdis
          sigLp = sigLdis
        endif
      endif

      return
      end



      SUBROUTINE RESCSN(w2,q2,sigtn,sigln)
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xvaln(100),xval1(50),xvalL(50)
      integer i

      data xvaln / 
c     & 0.12300E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
c     & 0.14333E+01,0.13223E+00,0.21779E+00,0.10504E+00,0.19322E+00,
c     & 0.16957E+00,0.38000E+00,0.76000E+01,0.44484E+01,0.22412E+01,
c     & 0.18607E+01,0.13787E-02,0.98347E+04,0.18160E+00,0.26770E+01,
c     & 0.34272E+00,0.18609E-05,0.47226E+01,0.63725E-06,0.23236E-01,
c     & 0.84561E+03,0.29815E+00,0.27082E+01,0.00000E+00,0.10000E+01,
c     & 0.10000E+01,0.20000E+01,0.67680E+01,0.56785E+01,0.43946E+00,
c     & 0.56784E+01,0.20954E+03,0.18496E+00,0.16769E+01,0.15454E+00,
c     & -.83108E+02,0.30574E+01,0.10000E-02,0.93889E+00,-.89302E-02,
c     & -.62215E-01,0.19800E+01,0.45000E+00,0.40853E+01,0.14462E+00,
c     & 0.10046E+01,0.99612E+00,0.99930E+00,0.99357E+00,0.10192E+01,
c     & 0.10000E+01,0.99669E+00,0.10002E+01,0.10000E+01,0.10000E+01,
c     & 0.10094E+01,0.10064E+01,0.16270E+03,0.86884E+01,0.57864E+01,
c     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
c     & 0.30131E+02,0.55479E+01,0.93622E+00,0.00000E+00,0.76450E+02,
c     & 0.66412E-01,0.85428E+01,0.00000E+00,0.00000E+00,0.20000E+01,
c     & 0.10000E+01,0.00000E+00,0.23539E+01,0.10360E-02,0.67452E+00,
c     & 0.00000E+00,0.23312E+03,0.35946E+00,0.23000E+01,-.30155E-01,
c     & 0.59653E+01,-.10616E+01,0.00000E+00,0.00000E+00,0.13868E+04,
c     & 0.29402E+03,0.39655E+00,0.00000E+00,0.26453E+00,0.00000E+00 /

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13223E+00,0.21779E+00,0.10504E+00,0.19322E+00,
     & 0.16957E+00,0.38000E+00,0.77664E+01,0.39601E+01,0.43374E+00,
     & 0.25344E+01,0.28233E+03,0.25363E+04,0.15952E+04,0.71568E+05,
     & 0.90380E-01,0.33761E+03,0.41783E+01,0.17259E+01,0.21868E-01,
     & 0.91452E+03,0.27680E+00,0.26991E+01,0.75512E-02,0.76650E+01,
     & 0.98129E-01,0.47138E+00,0.34484E+01,0.26717E+03,0.17708E+03,
     & 0.38133E+01,0.25583E+03,0.29679E+00,0.20201E+01,0.98158E-01,
     & -.11203E+03,0.11255E+02,0.10730E-03,0.71866E+00,-.54297E-02,
     & -.40653E-01,0.19900E+01,0.40113E+00,0.48000E+01,0.14462E+00,
     & 0.10046E+01,0.99441E+00,0.99961E+00,0.99589E+00,0.10185E+01,
     & 0.10000E+01,0.10017E+01,0.99948E+00,0.10000E+01,0.10000E+01,
     & 0.99720E+00,0.99733E+00,0.14451E+03,0.46369E+01,0.58155E+01,
     & 0.00000E+00,0.63794E+01,0.84321E-01,0.16109E+01,0.00000E+00,
     & 0.36950E+02,0.40077E-01,0.58459E+01,0.00000E+00,0.10358E+03,
     & 0.17872E-03,0.94839E+01,0.00000E+00,0.00000E+00,0.20000E+01,
     & 0.10000E+01,0.00000E+00,0.40315E+01,0.55366E+00,0.46885E+00,
     & 0.00000E+00,0.49373E+04,0.14302E+00,0.36927E+01,-.49413E-01,
     & 0.67573E+01,-.67472E+00,0.00000E+00,0.00000E+00,0.71548E+03,
     & 0.17543E+03,0.35600E+00,0.00000E+00,0.20010E+01,0.00000E+00 /

      do i=1,50
        xval1(i) = xvaln(i)
        xvalL(i) = xvaln(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodn(1,w2,q2,xval1,sigtn)
      call resmodn(2,w2,q2,xvalL,sigLn)

      return
      end




      SUBROUTINE RESMODP(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 



        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       
      if(w2.LT.1.159) sig = 0.0

 1000  format(8f12.5)

      RETURN 
      END 
  


      SUBROUTINE RESMODN(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.939565

c      mp = 0.938272

      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)
      if(w.LT.(mp+mpi)) wdif(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) wdif(2) = 0.0

c      write(6,*) "here"

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+2.*mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LT.(mp+mpi)) xpr(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) xpr(2) = 0.0

      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)*dip

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
c     &      (1.-xb)**(nr_coef(i,2)+0.0*t )
     &      (1.-xpr(1))**(nr_coef(i,2)+nr_coef(i,3)*log(q2+m0))
c     &        *xpr(1)**(xval(41)+xval(42)*t)
     &        *xb**(xval(41)+xval(42)*log(q2+m0)+
     &                  nr_coef(i,4)*q2*sqrt(q2))
c     &        *(1./(1.+q2/nr_coef(i,3)))**nr_coef(i,4)
        enddo

CCC  next line is test CCC nonres3 directory ******

c        sig_nr = sig_nr + nr_coef(i,3)/(1.+q2/nr_coef(i,4))**2*
c     &               (1.-xpr(1))**2

        

      endif


      sig_res = abs(sig_res)
      sig_nr = abs(sig_nr)

      sig = sig_res + sig_nr
       
      if(w.LT.(mp+mpi)) sig = 0.0

 1000  format(8f12.5)

      RETURN 
      END 

      SUBROUTINE DISP(w2,q2,sigt,sigl)
      IMPLICIT none

      real*8 w2,q2,q2t,x,sigt,sigl,f1,f2,fL,r,dr,r1,r2
      Integer i,modt
      character*1 targ
      real*8 mp2,pi,pi2,alpha,t1,t2
      logical goodfit


      targ = 'P'
      modt = 12
      mp2 = 0.938272
      mp2 = mp2*mp2
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.0359

      x = q2/(w2+q2-mp2)

      call f2allm(x,q2,f2)
      call gf2glob(x,q2,targ,modt,f2)
      call r1998(x,q2,r,dr,goodfit)
      if(q2.LE.0.15) then
        q2t = 0.15
        call r1998(x,q2t,r1,dr,goodfit)
        r = r1/0.25*q2
c        write(6,*) r1,r
        call f2allm(x,q2,f2)
      endif


      f1 = f2/2./x/(r+1.0)*(1.0+4.0*mp2*x*x/q2)
      fL = 2.*x*r*f1    

      sigt = 0.3894e3*f1*pi2*alpha*8.0/abs(w2-mp2)
      sigL = r*sigt

      return
      end


      SUBROUTINE SQESUB(w2,q2,wfn,f2s,f1s,fLs,first)
      IMPLICIT none

      real*8 w2,wwz,q2,x,y,z,Z1,A,t1
      real*8 nu,nuel,q4,q,gamma,gamma2,alpha,alpha2
      real*8 epel,w2el,kappa,mp,mp2,md,w2min
      real*8 f1,f2,fy11,fy12,fy2,f2s,f1s,fLs,f1sos,f2sos,fLsos
      real*8 ymin,pi,pi2,GM,GE,GMp,GMn,GEp,GEn,tau,mcor,ecor
      real*8 GM2,GE2,GD,mu_p,mu_n,a1,a2,a3,b1,b2,b3,b4,b5
      integer i,j,k,bin,off,wfn
      logical thend,first,firsty
      real*8 f1d_of,sf_of
      external f1d_of,sf_of

      thend = .false.
      firsty = .false.
      if(first) firsty = .true.
      
      Z1 = 1.
      A = 2.
      mp = (0.938272+0.939565)/2.0d0    !!! average p,n
      mp2 = mp*mp
      md = 1.8756
      alpha = 1/137.03599
      alpha2 = alpha*alpha 
      pi = 3.14159
      pi2 = pi*pi
      mu_p = 2.793D0
      mu_n = -1.913148
    

CCC///   Read in smearing function array  ///CCC
CCC

      nu = (w2+q2-mp2)/2./mp   !!! Fermi smeared !!!

      nuel = q2/2./mp  !!! In nucleon rest frame !!!
      q = dsqrt(q2)
      q4 = q2*q2
      x = q2/(w2+q2-mp2)
      tau = q2/(4*mp2)  !!!  for elastic at this Q^2  !!!
c      ww = mp2 + 2.d0*mp*nu - q2  !!! Fermi smeared  !!!
      gamma2 = (1.+4.*mp2*x*x/q2)
      if(gamma2.LE.6.and.firsty) then   
       gamma = sqrt(gamma2)
       call getfy(gamma,x,wfn,fy11,fy12,fy2,firsty)
      endif

      if(nu.LT.(q2/2.0/md)) then
        f1s = 0.0D0
        f2s = 0.0D0
        fLs = 0.0D0
        return
      endif
   

CCCC    Calculate elastic formfactors at same Q^2  CCCC

        GMp = mu_p/(1.-0.29940*q+5.65622*q2-5.57350*q*q2+      !  with Hall A data  !
     &         6.31574*q2**2-1.77642*q*q2**2+0.30100*q2**3)

        GEp = 1./(1.-0.31453*q+7.15389*q2-14.42166*q*q2+
     &         23.33935*q2**2-14.67906*q*q2**2+4.03066*q2**3)

        GMp = 2.792782     !!!    
     &       * (1.D0 -1.465D0*tau + 1.260D0*tau**2 + 0.262D0*tau**3)
     &       / (1.D0 + 9.627D0*tau + 0.0D0*tau**2 + 0.0D0*tau**3
     &               + 11.179D0*tau**4 + 13.245D0*tau**5)


        a1 = -1.651D0
        a2 = 1.287D0
        a3 = -0.185D0
        b1 = 9.531D0
        b2 = 0.591D0
        b3 = 0.D0
        b4 = 0.D0
        b5 = 4.994D0

        GEp = (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)

        a1 = -2.151D0
        a2 = 4.261D0
        a3 = 0.159D0
        b1 = 8.647D0
        b2 = 0.001D0
        b3 = 5.245D0
        b4 = 82.817D0
        b5 = 14.191D0

        GMp = mu_p
     &      * (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)

c        GEp = GMp/mu_p   !!! test !!!


        GD = (1./(1 + q2/0.71))**2

        GMn = mu_n           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &       / (1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)

c        GMn = mu_n*GD


        GEn = (1.700D0*tau / (1+ 3.300D0*tau)) * GD

        mcor = 1.0-0.09*q2/(0.1+q2)
        ecor = 1.0+0.15*q2/(0.1+q2)
        gmn = mcor*gmn  !!! test
        gen = ecor*gen  !!! test
  
c        GEp = 1.0*GD 
c        GMP = 0.97*GMP  !!! TEST     

c        GEn = 0.0  !!!  TEST  !!!

        GM2 = GMp*GMp + GMn*GMn  !!! Sum p+n  !!!
        GE2 = GEp*GEp + GEn*GEn


        off = 1   !!! CC1
        call offshellqe(w2,q2,wfn,off,GM2,GE2,f1s,f2s)
        fLs = gamma2*f2s-2.*x*f1s

c        call f1f2qe09(Z1,A,q2,w2,f1s,f2s)
c        write(6,*) "sqesub:  ",w2,q2,f1s,f2s,fLs

c        if(fLs.LT.0.0) write(6,2000) w2,q2,f1s,f2s,fLs

        if(f1s.LT.0.0) f1s = 0.0
        if(f2s.LT.0.0) f2s = 0.0
        if(fLs.LT.0.0) fLs = 0.0
 
c       if(w2.LT.1.0) write(6,2000) w2,q2,f1s,f2s,fLs

 2000 format(8f9.4)
c 2000 format(5f9.4,2e11.3)

      end





      SUBROUTINE GETFY(gamma,ypass,wfn,fy11,fy12,fy2,firsty)
      IMPLICIT none
      real*8 gamma,y,ypass,gammav(50),fy,yv(50,1001),yv2(50,1001)
      real*8 yv12(50,1001),fyv11(50,1001),fyv2(50,1001),fyv12(50,1001)
      real*8 t1,fylow,fyhi,fy11,fy12,fy2
      integer i,j,k,bin,bin2,jlow,jhi,wfn
      logical thend,firsty
      COMMON/FYCOM/ yv,yv2,yv12,fyv11,fyv2,fyv12,gammav

      thend = .false.

      y = ypass
      if(ypass.GE.1.9999) y = 1.999  !!! For numerical stability

CCC///   Read in smearing function array  ///CCC

c      firsty = .true.
      i = 0
      if(firsty) then
        if(wfn.EQ.1) then
          open(unit=34,file='f11.WBARELav18',status='old')        
          open(unit=35,file='f2.WBARELav18',status='old')
          open(unit=36,file='f12.WBARELav18',status='old')
        elseif(wfn.EQ.2) then
          open(unit=34,file='f11.WBARELcdbonn',status='old')        
          open(unit=35,file='f2.WBARELcdbonn',status='old')
          open(unit=36,file='f12.WBARELcdbonn',status='old')
        elseif(wfn.EQ.3) then
          open(unit=34,file='f11.WBARELwjc1',status='old')        
          open(unit=35,file='f2.WBARELwjc1',status='old')
          open(unit=36,file='f12.WBARELwjc1',status='old') 
        elseif(wfn.EQ.4) then
          open(unit=34,file='f11.WBARELwjc2',status='old')        
          open(unit=35,file='f2.WBARELwjc2',status='old')
          open(unit=36,file='f12.WBARELwjc2',status='old')
        else
          write(6,*)  "Warning:  wavefunction undefined"
        endif

     
        dowhile(.not.thend)   !!!  First read in smearing function to temp vector!!!
          i = i+1
          j = int(float(i-1)/1001.)+1 
          k = i-(j-1)*1001
          read(34,*,end=999) gammav(j),yv(j,k),fyv11(j,k)
          read(35,*,end=999) gammav(j),yv2(j,k),fyv2(j,k)
          read(36,*,end=999) gammav(j),yv12(j,k),fyv12(j,k)

        enddo
999     thend = .true. 
      endif
      firsty = .false.
      close(34)
      close(35)
      close(36)

      bin = int(y/2.*1000)+1          !!! find y-bin below elastic at each x,Q^2
      jlow = int((gamma-1.0)/0.2)+1   !!! Assumes 0.2 increments in gamma
      jhi = jlow+1
      bin2 = bin+1

      fylow = (fyv11(jlow,bin2)*(y-yv(jlow,bin))+fyv11(jlow,bin)
     &        *(yv(jlow,bin2)-y))/0.002
      fyhi =  (fyv11(jhi,bin2)*(y-yv(jhi,bin))+fyv11(jhi,bin)
     &        *(yv(jhi,bin2)-y))/0.002

      fy11 = (fylow*(gammav(jhi)-gamma)+fyhi*(gamma-gammav(jlow)))/0.2

c      write(6,*) y,fyv11(jlow,bin)

      fylow = (fyv12(jlow,bin2)*(y-yv(jlow,bin))+fyv12(jlow,bin)
     &        *(yv(jlow,bin2)-y))/0.002
      fyhi =  (fyv12(jhi,bin2)*(y-yv(jhi,bin))+fyv12(jhi,bin)
     &        *(yv(jhi,bin2)-y))/0.002

      fy12 = (fylow*(gammav(jhi)-gamma)+fyhi*(gamma-gammav(jlow)))/0.2


      fylow = (fyv2(jlow,bin2)*(y-yv(jlow,bin))+fyv2(jlow,bin)
     &        *(yv(jlow,bin2)-y))/0.002
      fyhi =  (fyv2(jhi,bin2)*(y-yv(jhi,bin))+fyv2(jhi,bin)
     &        *(yv(jhi,bin2)-y))/0.002

      fy2 = (fylow*(gammav(jhi)-gamma)+fyhi*(gamma-gammav(jlow)))/0.2


      if(fy11.LT.0) fy11 = 0.0d0
      if(fy2.LT.0) fy2 = 0.0d0
      if(fy12.LT.0) fy12 = 0.0d0


c      if(y.LT.0.015) then   !!! fix for numerical issues in fy at small x
c        fy11 = 0.0
c        fy12 = 0.0
c        fy2 = 0.0
c      endif

     
 2000 format(5f9.4,2e11.3)

      return
      end




      SUBROUTINE offshellqe(w2,q2,wfn,off,GM2,GE2,f1os,f2os)
CCC    Original code by M.E. Christy with significant borrowing of code from    CCC
CCC    J. Ethier, W. Melnitchouk, et.al.                                        CCC
CCC    Addition of CC1 offshell by J. Ethier (July 1, 2014)                     CCC
CCC    See reference (put references here)                                      CCC
CCC    Oct. 9, fixed bug to keep from returning NAN if a < b                    CCC
c  Nucleon light-cone momentum distribution function in deuteron in
C  weak binding approximation (WBA), as a function of the fraction (y)
C  of deuteron's momentum carried by nucleon, defined as
C		y = p.q/p_D.q
C  (sometimes labeled as y_D) so that in Bjorken limit y -> y0 in [0,1].
C
C  Kulagin, Petti, NPA 765, 126 (2006), Eq. (43)
C  Kahn, WM, Kulagin, PRC 79, 035205 (2009)
C
C  Relativistic kinematics implemented (WM): May 2010
C  Relativistic convolution implemented (cf. AQV): Sep 2010
C  
C  All momenta in MeV.
C ***********************************************************************

      IMPLICIT none

      real*8 x,w2,qsqr,q2,GM2,GE2,FF,pFF,off1,off2,off3,corr
      real*8 q,gamma,gamma2,mp,mp2,md,md2,hc,tau,xth,y,ep,p2
      real*8 pv,pv2,pz,pz2,pt2,xd,w2d,f1,f2,ft,fttilda,f2tilda,flux
      real*8 cc,pvmin,pvmax,a,b,c,rt,prod,ftv,fttildav(10001),f2v
      real*8 f2tildav(10001),f1os,f2os,dsimps,dgauss,pvinc,wfnorm(10)
      real*8 U,W,VS,VT
      external dsimps,dgauss
      integer i,j,nbins,wfn,off
      logical thend,first,firsty
      data wfnorm / 1.0,1.0,1.05,1.02,1.0,1.0,1.0,1.0,1.0,1.0 /

c      write(6,2001) x,q2,gm2,ge2

c      off = 2

      f1os = 0.0D0
      f2os = 0.0D0

      nbins = 1000

c      do i=1,10
c        wfnorm(i) = 1.0D0
c      enddo

c      write(6,*) w2,q2,wfn

      hc = 0.197327D0      !!! convert from GeV -> fm    !!!
      mp = (0.938272+0.939565)/2.0D0    !!! average p,n  !!!
      mp2 = mp*mp
      md = 2.0D0*mp-0.002224575D0
    
      x = q2/(w2+q2-mp2)
      xd = x*mp/md
      w2d = md*md+q2*(1.0D0/xd-1.0D0)

      qsqr = q2*(1.D0 + ( q2 / (4.D0*mp2*x*x)))

      tau = q2/(4.0D0*mp2)  !!!  for elastic at this Q^2  !!!
      gamma2 = (1.+4.0D0*mp2*x*x/q2)   
      gamma = sqrt(gamma2)
      a = dsqrt(1.d0 + w2d/qsqr)
      b = w2d/(2.D0*mp*dsqrt(qsqr))
      prod = b/(a*a - 1.D0)
      c = 1.D0+((a*a-1.0D0)*(1.0-a*a/b/b))
 
      if(c.GE.0.0D0) then
        rt = dsqrt(c)
      else
        return
      endif

      ff = (GE2+tau*GM2)/(1.D0+tau)
      pff = ((sqrt(GE2)-sqrt(GM2))/(1.D0+tau))**2      !!! Pauli FF squared

      pvmin = mp*prod*(1.d0-rt)*sign(1.D0,(a-b))
      pvmax = mp*prod*(1.D0+rt)
      pvinc = (pvmax-pvmin)/float(nbins)

c      write(6,*) "offshellqe:  ",w2,q2,x,rt

      pv = pvmin-pvinc
CCC  Loop over initial nucleon 3-momentum, pv  CCC
      do i=1,nbins+1              !!! integrate over nucleon 3-momentum 
        fttildav(i) = 0.0D0
        f2tildav(i) = 0.0D0
        pv = pv+pvinc
        pv2 = pv*pv              !!! 3-momentum squared of nucleon 
        ep = sqrt(mp2+pv2)       !!! total energy of nucleon
        p2 = (md-ep)**2-pv2      !!! vituality of nucleon
        xth =  q2/(q2-p2+mp2)    !!! x for given virutuality  
        y = x/xth
        pz = (y*mp-md+ep)/gamma  !!! nucleon momentum along q-vector
        pz2 = pz*pz
        pt2 = pv2-pz2            !!! nucleon transverse momentum
        flux = 1.0D0+gamma*pz/mp

c         if(w2.LT.1.0) write(6,2000) w2,q2,f1s,f2s,fLs,f1sos,f2sos,fLsos

        off1 = (mp2-p2)/q2
        off2 = (mp2-p2)/4./mp2
        corr = 1.D0+xth*xth*(4.D0*p2+6.D0*pt2)/q2

        U = 0.D0        ! S-wave
        W = 0.D0        ! D-wave
        VS = 0.D0       ! Singlet P-wave
        VT = 0.D0       ! Triplet P-wave

        IF (wfn.EQ.2) CALL CDBONN(pv/hc,U,W)
        IF (wfn.EQ.3) CALL WJC1(pv/hc,U,W,VS,VT)
              ! wfns normalized to ~105% (5% in V' term)
        IF (wfn.EQ.4) CALL WJC2(pv/hc,U,W,VS,VT)
              ! wfns normalized to ~102%

C...Output wavefunctions in fm^3/2 => MeV^-3/2
           U = U / hc**1.5D0
           W = W / hc**1.5D0
           VS = VS / hc**1.5D0
           VT = VT / hc**1.5D0
           CC = (U**2 + W**2 + VS**2 + VT**2)
        if(wfn.EQ.1) then
          call av18(pv/hc,cc)    
          cc = cc/(hc*hc*hc) 
        endif 

c        write(6,*) "offshellqe:  ", x,pv,p2,xth,pvmin,pvmax

        if(off.EQ.1) then    !!!  CC1 
          ftv = x*(x/y)*(gm2*(1.D0-off1))
          f2v = y*(ff-off2*pff)
          fttilda = ftv + (2.0D0*xth*xth*pt2*f2v/q2)
          fttildav(i) = fttilda*flux*cc*mp/(4.*gamma)
          fttildav(i) = fttildav(i)*2.D0*pv  
          f2tilda = (mp*x*(ff-off2*pff))/(4.d0*gamma**3)
          f2tildav(i) = f2tilda*flux*corr*cc/xth
          f2tildav(i) = f2tildav(i)*2.D0*pv
        else if(off.EQ.2) then    !!!  CC2
          ftv = x*(x/y)*(gm2-off1*(ff-off2*pff))
          f2v = y*ff
          fttilda = ftv + (2.0D0*xth*xth*pt2*f2v/q2)
          fttildav(i) = fttilda*flux*cc*mp/(4.*gamma)
          fttildav(i) = fttildav(i)*2.D0*pv
          f2tilda = (mp*x*ff)/(4.d0*gamma**3) 
          f2tildav(i) = f2tilda*flux*corr*cc/xth
          f2tildav(i) = f2tildav(i)*2.D0*pv
        endif

      enddo

      ft = dsimps(fttildav,pvmin,pvmax,nbins)
      f2os = dsimps(f2tildav,pvmin,pvmax,nbins)

      f1os = ft/2./x

 
      f1os = f1os/wfnorm(wfn)
      f2os = f2os/wfnorm(wfn)


c      if(f1os.LT.0.0.OR.f1os.LT.0.0) 
c     &      write(6,2001) w2,q2,f1os,f2os,gamma2


      if(f1os.LT.0.0D0) f1os = 0.0D0
      if(f2os.LT.0.0D0) f2os = 0.0D0

 2000 format(8f9.4)
 2001 format(5f9.4)

      return 
      end










C ***********************************************************************
        FUNCTION off_mKP_fit (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]), constrained to vanish at x=0
C
C  Defined such that F2d = F2d(conv) + del^off F2d 
C       with off_mKP_fit = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Fit by Eric Christy (6/14/12).
C  Fortran code by W. Melnitchouk (6/17/12).
C
C  wfn: 1 (AV18)
C	2 (CD-Bonn)
C	3 (WJC-1)
C	4 (WJC-2)
C
C  dRN: 0 (0.0%)	[% change in nucleon radius in the deuteron]
C	1 (0.3%)
C	2 (0.6%)
C	3 (0.9%)
C	4 (1.2%)
C	5 (1.5%)
C	6 (1.8%)
C	7 (2.1%)
C	8 (2.4%)
C	9 (2.7%)
C      10 (3.0%)
C
C ***********************************************************************
	IMPLICIT NONE
	REAL*8	off_mKP_fit, x
	INTEGER	wfn, dRN
	REAL*8  p(0:9)

	off_mKP_fit = 0.D0
	IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
	IF (dRN.LT.0 .OR. dRN.GT.10) RETURN
! .......................................................................
	IF (wfn.EQ.1) THEN		! AV18
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02628D0
	    p(1) = 1.20029D0
	    p(2) = 7.49503D0
	    p(3) = 2.01901D0
	    p(4) = 0.00789D0
	    p(5) = 0.46739D0
	    p(6) = 0.73242D0
	    p(7) = 0.00328D0
	    p(8) = 0.87228D0
	    p(9) = 0.06400D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = 0.03638D0
	    p(1) = 0.38307D0
	    p(2) = 8.01156D0
	    p(3) = 2.30992D0
	    p(4) = 0.09027D0
	    p(5) = 0.69521D0
	    p(6) = 0.75973D0
	    p(7) = -0.05098D0
	    p(8) = 1.18963D0
	    p(9) = -0.19192D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02260D0
	    p(1) = 1.45377D0
	    p(2) = 0.50628D0
	    p(3) = 13.92200D0
	    p(4) = 0.03558D0
	    p(5) = 0.75147D0
	    p(6) = 0.86335D0
	    p(7) = -0.01383D0
	    p(8) = 1.04749D0
	    p(9) = 0.42099D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06410D0
	    p(1) = 1.18883D0
	    p(2) = 6.96799D0
	    p(3) = 8.87113D0
	    p(4) = 0.02603D0
	    p(5) = 0.70504D0
	    p(6) = 1.44139D0
	    p(7) = 0.00004D0
	    p(8) = -1.14305D0
	    p(9) = 0.73785D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06237D0
	    p(1) = 2.03192D0
	    p(2) = 4.01755D0
	    p(3) = 6.83741D0
	    p(4) = 0.04701D0
	    p(5) = -0.00457D0
	    p(6) = 1.30967D0
	    p(7) = -0.00996D0
	    p(8) = -0.42418D0
	    p(9) = 0.27524D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.06759D0
	    p(1) = 1.95103D0
	    p(2) = 3.54215D0
	    p(3) = 11.77533D0
	    p(4) = 0.09269D0
	    p(5) = 0.56534D0
	    p(6) = 0.98398D0
	    p(7) = -0.03031D0
	    p(8) = 3.26913D0
	    p(9) = -0.45923D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.07007D0
	    p(1) = 2.30938D0
	    p(2) = 4.94226D0
	    p(3) = 8.95701D0
	    p(4) = 0.06933D0
	    p(5) = 0.07145D0
	    p(6) = 1.94887D0
	    p(7) = -0.01210D0
	    p(8) = 5.92311D0
	    p(9) = 0.14312D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11965D0
	    p(1) = 2.06149D0
	    p(2) = 5.38881D0
	    p(3) = 12.08265D0
	    p(4) = 0.19668D0
	    p(5) = 0.61820D0
	    p(6) = 0.80489D0
	    p(7) = -0.08735D0
	    p(8) = 3.74802D0
	    p(9) = -0.70773D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.14735D0
	    p(1) = 2.27109D0
	    p(2) = 8.23092D0
	    p(3) = 7.31581D0
	    p(4) = 0.11953D0
	    p(5) = 0.67459D0
	    p(6) = 1.59118D0
	    p(7) = -0.02700D0
	    p(8) = 4.52840D0
	    p(9) = -1.77765D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.27194D0
	    p(1) = 2.01340D0
	    p(2) = 10.71380D0
	    p(3) = 8.84886D0
	    p(4) = 0.09345D0
	    p(5) = 0.49802D0
	    p(6) = 1.28523D0
	    p(7) = -0.00474D0
	    p(8) = 0.58703D0
	    p(9) = 0.88354D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.69848D0
	    p(1) = 1.48173D0
	    p(2) = 17.44991D0
	    p(3) = 12.73730D0
	    p(4) = 0.13118D0
	    p(5) = 0.34598D0
	    p(6) = 1.65884D0
	    p(7) = -0.02215D0
	    p(8) = 1.21306D0
	    p(9) = 0.96399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.2) THEN		! CD-Bonn
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02820D0
	    p(1) = 0.85879D0
	    p(2) = 9.48856D0
	    p(3) = 2.18885D0
	    p(4) = 0.00070D0
	    p(5) = -5.61817D0
	    p(6) = 14.80512D0
	    p(7) = 0.00348D0
	    p(8) = -1.30292D0
	    p(9) = -0.73075D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.02996D0
	    p(1) = 0.35717D0
	    p(2) = 6.53843D0
	    p(3) = 3.88389D0
	    p(4) = 0.00758D0
	    p(5) = -16.50399D0
	    p(6) = 77.60083D0
	    p(7) = 0.00320D0
	    p(8) = 0.42334D0
	    p(9) = 0.28545D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03261D0
	    p(1) = 0.91185D0
	    p(2) = 8.49348D0
	    p(3) = 10.19681D0
	    p(4) = 0.01598D0
	    p(5) = 0.83748D0
	    p(6) = 1.55960D0
	    p(7) = 0.00085D0
	    p(8) = -0.63447D0
	    p(9) = 0.65632D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.03034D0
	    p(1) = 1.58677D0
	    p(2) = 3.21753D0
	    p(3) = 11.66572D0
	    p(4) = 0.04999D0
	    p(5) = 0.56688D0
	    p(6) = 0.94941D0
	    p(7) = -0.01453D0
	    p(8) = -0.89157D0
	    p(9) = 0.27160D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.04831D0
	    p(1) = 1.75241D0
	    p(2) = 4.74662D0
	    p(3) = 8.29052D0
	    p(4) = 0.04730D0
	    p(5) = 0.33550D0
	    p(6) = 1.18790D0
	    p(7) = -0.00678D0
	    p(8) = -0.42800D0
	    p(9) = 0.36573D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.09019D0
	    p(1) = 1.22091D0
	    p(2) = 1.30114D0
	    p(3) = 17.58701D0
	    p(4) = 0.08312D0
	    p(5) = 0.66902D0
	    p(6) = 0.60767D0
	    p(7) = -0.02035D0
	    p(8) = 0.95978D0
	    p(9) = 1.11322D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.17060D0
	    p(1) = 1.42115D0
	    p(2) = 7.24672D0
	    p(3) = 5.80680D0
	    p(4) = 0.09200D0
	    p(5) = 0.43367D0
	    p(6) = 1.56378D0
	    p(7) = -0.02338D0
	    p(8) = 0.44968D0
	    p(9) = 0.29678D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11026D0
	    p(1) = 1.85213D0
	    p(2) = 6.74413D0
	    p(3) = 7.74362D0
	    p(4) = 0.08467D0
	    p(5) = 0.24708D0
	    p(6) = 1.12274D0
	    p(7) = -0.01505D0
	    p(8) = 0.44209D0
	    p(9) = 0.36126D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.15291D0
	    p(1) = 1.83333D0
	    p(2) = 7.76495D0
	    p(3) = 7.04783D0
	    p(4) = 0.09206D0
	    p(5) = 0.08655D0
	    p(6) = 1.27460D0
	    p(7) = -0.01659D0
	    p(8) = 0.45536D0
	    p(9) = 0.29407D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.24143D0
	    p(1) = 1.50401D0
	    p(2) = 9.33393D0
	    p(3) = 11.62779D0
	    p(4) = 0.09454D0
	    p(5) = 0.36361D0
	    p(6) = 0.82058D0
	    p(7) = -0.00802D0
	    p(8) = 0.34851D0
	    p(9) = 0.50844D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.22196D0
	    p(1) = 1.87228D0
	    p(2) = 10.18898D0
	    p(3) = 9.21038D0
	    p(4) = 0.11850D0
	    p(5) = 0.34360D0
	    p(6) = 1.28278D0
	    p(7) = -0.01754D0
	    p(8) = 0.54540D0
	    p(9) = 0.53457D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.3) THEN		! WJC-1
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.02322D0
	    p(1) = 0.11213D0
	    p(2) = 3.71079D0
	    p(3) = 5.51496D0
	    p(4) = 0.00877D0
	    p(5) = 0.84639D0
	    p(6) = 0.66227D0
	    p(7) = -0.00621D0
	    p(8) = -0.39896D0
	    p(9) = 0.32012D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.00058D0
	    p(1) = 2.33827D0
	    p(2) = 2.35664D0
	    p(3) = 36.75823D0
	    p(4) = -0.00752D0
	    p(5) = 0.05286D0
	    p(6) = 1.27262D0
	    p(7) = 0.01269D0
	    p(8) = 1.72720D0
	    p(9) = 0.20652D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03373D0
	    p(1) = 0.93858D0
	    p(2) = 0.15704D0
	    p(3) = 10.71630D0
	    p(4) = -0.00235D0
	    p(5) = -0.11937D0
	    p(6) = 0.74925D0
	    p(7) = 0.00452D0
	    p(8) = 2.96830D0
	    p(9) = -2.89070D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.08982D0
	    p(1) = 0.73060D0
	    p(2) = -0.16543D0
	    p(3) = 12.37035D0
	    p(4) = 0.04407D0
	    p(5) = 0.47361D0
	    p(6) = 0.74570D0
	    p(7) = -0.00933D0
	    p(8) = 0.53186D0
	    p(9) = 0.26943D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.11990D0
	    p(1) = 1.19824D0
	    p(2) = 3.06386D0
	    p(3) = 8.55017D0
	    p(4) = 0.05815D0
	    p(5) = 0.06123D0
	    p(6) = 1.45024D0
	    p(7) = -0.01414D0
	    p(8) = 0.48172D0
	    p(9) = 0.25171D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.15292D0
	    p(1) = 1.01991D0
	    p(2) = 1.20661D0
	    p(3) = 13.31860D0
	    p(4) = 0.02571D0
	    p(5) = -1.56438D0
	    p(6) = 2.69042D0
	    p(7) = -0.00000D0
	    p(8) = 0.29759D0
	    p(9) = 0.97967D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.35935D0
	    p(1) = 0.44637D0
	    p(2) = -0.25510D0
	    p(3) = 16.70057D0
	    p(4) = 0.10634D0
	    p(5) = 0.61659D0
	    p(6) = 0.58524D0
	    p(7) = -0.03335D0
	    p(8) = 0.93904D0
	    p(9) = 0.89819D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.97384D0
	    p(1) = -0.24934D0
	    p(2) = -0.61349D0
	    p(3) = 18.43254D0
	    p(4) = 0.18772D0
	    p(5) = 0.49599D0
	    p(6) = 0.61366D0
	    p(7) = -0.08116D0
	    p(8) = 0.87175D0
	    p(9) = 0.24026D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.21641D0
	    p(1) = 1.74710D0
	    p(2) = 5.19387D0
	    p(3) = 10.61285D0
	    p(4) = 0.06655D0
	    p(5) = 0.01300D0
	    p(6) = 0.94503D0
	    p(7) = -0.00642D0
	    p(8) = 0.48859D0
	    p(9) = 0.16331D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.32283D0
	    p(1) = 1.71708D0
	    p(2) = 7.51556D0
	    p(3) = 9.68202D0
	    p(4) = 0.09871D0
	    p(5) = 0.18788D0
	    p(6) = 0.80490D0
	    p(7) = -0.01673D0
	    p(8) = 0.48879D0
	    p(9) = 0.21016D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.32064D0
	    p(1) = 2.07514D0
	    p(2) = 9.34847D0
	    p(3) = 8.17225D0
	    p(4) = 0.10772D0
	    p(5) = 0.50272D0
	    p(6) = 1.30663D0
	    p(7) = -0.01215D0
	    p(8) = 0.59432D0
	    p(9) = 0.65399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.4) THEN		! WJC-2
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.03490D0
	    p(1) = 0.78902D0
	    p(2) = -0.25256D0
	    p(3) = 7.98679D0
	    p(4) = 0.00913D0
	    p(5) = 0.74835D0
	    p(6) = 0.60145D0
	    p(7) = -0.00464D0
	    p(8) = 0.41358D0
	    p(9) = 0.22524D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.01119D0
	    p(1) = 0.50514D0
	    p(2) = 19.35710D0
	    p(3) = 3.32395D0
	    p(4) = 0.00670D0
	    p(5) = 1.38279D0
	    p(6) = 1.24216D0
	    p(7) = 0.00049D0
	    p(8) = 0.38623D0
	    p(9) = 0.23497D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02653D0
	    p(1) = 1.27315D0
	    p(2) = -0.53410D0
	    p(3) = 14.08029D0
	    p(4) = 0.01474D0
	    p(5) = 1.82129D0
	    p(6) = 1.99455D0
	    p(7) = -0.00090D0
	    p(8) = 3.96583D0
	    p(9) = 4.61316D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06301D0
	    p(1) = 1.10373D0
	    p(2) = -0.26356D0
	    p(3) = 15.04038D0
	    p(4) = 0.02428D0
	    p(5) = -0.15349D0
	    p(6) = 3.03168D0
	    p(7) = 0.00127D0
	    p(8) = -0.73818D0
	    p(9) = 0.07474D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06150D0
	    p(1) = 2.15792D0
	    p(2) = 2.18241D0
	    p(3) = 9.84713D0
	    p(4) = 0.03608D0
	    p(5) = -0.13604D0
	    p(6) = 1.12241D0
	    p(7) = -0.00695D0
	    p(8) = -0.35646D0
	    p(9) = 0.31793D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.07179D0
	    p(1) = 1.97917D0
	    p(2) = 3.47662D0
	    p(3) = 10.00224D0
	    p(4) = 0.04587D0
	    p(5) = 0.06416D0
	    p(6) = 1.10677D0
	    p(7) = -0.00391D0
	    p(8) = -0.42677D0
	    p(9) = 0.26619D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.09883D0
	    p(1) = 1.96788D0
	    p(2) = 5.19182D0
	    p(3) = 8.82173D0
	    p(4) = 0.06468D0
	    p(5) = 0.11297D0
	    p(6) = 1.63850D0
	    p(7) = -0.00872D0
	    p(8) = 0.52753D0
	    p(9) = 0.41794D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.14258D0
	    p(1) = 2.00822D0
	    p(2) = 6.23508D0
	    p(3) = 7.81846D0
	    p(4) = 0.07064D0
	    p(5) = -0.05869D0
	    p(6) = 1.24848D0
	    p(7) = -0.01160D0
	    p(8) = 0.48932D0
	    p(9) = 0.22001D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.16184D0
	    p(1) = 2.16963D0
	    p(2) = 7.62378D0
	    p(3) = 7.33369D0
	    p(4) = 0.09197D0
	    p(5) = 0.15692D0
	    p(6) = 1.80734D0
	    p(7) = -0.01561D0
	    p(8) = 0.53224D0
	    p(9) = 0.39357D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.20205D0
	    p(1) = 2.28733D0
	    p(2) = 9.10375D0
	    p(3) = 7.24877D0
	    p(4) = 0.08325D0
	    p(5) = 0.36941D0
	    p(6) = 2.39131D0
	    p(7) = -0.00057D0
	    p(8) = 0.41640D0
	    p(9) = 0.90531D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.95664D0
	    p(1) = 1.11409D0
	    p(2) = 19.00631D0
	    p(3) = 15.97282D0
	    p(4) = 0.15616D0
	    p(5) = 0.40229D0
	    p(6) = 0.85878D0
	    p(7) = -0.03123D0
	    p(8) = 6.75437D0
	    p(9) = -3.83159D0
	  ENDIF

	ENDIF

	off_mKP_fit = -( p(0) * x**p(3) * DEXP(p(1) * x**p(2))
     &                + p(4) * x*DEXP(((x-p(5))/p(6))**2)
     &                + x**0.5D0 * p(7) * DEXP(((x-p(9))/p(8))**2) )

	RETURN
	END

C ***********************************************************************
	SUBROUTINE AV18 (p,rho)
C
C  Deuteron momentum distribution rho = u^2 + w^2, where u and w are
C  S- and D-state wave functions in momentum space, normalized s.t.
C  int dp p^2 rho(p) = 1.
C  Ref: Fantoni and Pandharipande, Nucl. Phys. A427 (1984) 473
C
C  Input file has rho normalized s.t. 4*pi int dp p^2 rho(p) = 1,
C  so that for output, rho needs to be multiplied by 4*pi.
C
C  p in 1/fm, rho in 1/fm^3
C
C.. Uses IMSL interpolation routine DQDVAL.
C.. For compilation on jlabs1:
C.. > use imsl
C.. > f77 -o objectfile file.f -R/site/vni/lib/lib.solaris
C..       -L/site/vni/lib/lib.solaris -limsl -lsocket -lnsl
C.. IMSL decommissioned 8/30/11
C
C ***********************************************************************
        IMPLICIT NONE
        INTEGER	ip,np
        PARAMETER (np=200)
        REAL*8  p,rho
        REAL*8  Parr(np),RHOarr(np),rho_int,dum
        LOGICAL readin /.FALSE./
	REAL*8	pi
	SAVE

C...Value of pi
        pi = 4*DATAN(1.D0)

	rho = 0.D0

        IF (readin) GO TO 123
C...Read data from file
c        OPEN (10,FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/av18.dat',
        OPEN (10,FILE='av18.dat',
     &		FORM='FORMATTED')
	DO ip=1,9
          READ (10,*)
	ENDDO
        DO ip=1,np
	  READ (10,*) Parr(ip), RHOarr(ip)
        ENDDO
        CLOSE (10)
        readin = .TRUE.
c        print *, '... AV18 data read ...'

 123	IF (p.LE.Parr(1)) rho = 4*pi * RHOarr(1)
	IF (p.GT.Parr(1) .AND. p.LE.Parr(np)) THEN
	  CALL Pinterp (Parr,RHOarr,np,p,rho_int,dum,2)
	  rho = 4*pi * rho_int
c     &    rho = 4*pi * DQDVAL(p,np,Parr,RHOarr,.FALSE.)
	ENDIF

c	print *, 'Parr(1),Parr(np)=',Parr(1),Parr(np)
c	print *, 'p,rho=',P, RHO

        RETURN
        END


C **********************************************************************
	SUBROUTINE CDBONN (q,u0,u2)
C
C  Deuteron wave function from CD-Bonn NN potential model.
C  q in 1/fm, u0,u2 in fm^3/2.
C
C  Normalization \int dq q^2 (u0^2+u2^2) = 1.
C
C  Sent by Charlotte Elster, April 8, 2009.
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=95)
        REAL*8  q,u0,u2,
     &		qgrid(nq),uqgrid(nq),wqgrid(nq),weight(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hc
	SAVE

        pi = 4*DATAN(1.D0)
        hc = 197.327D0		! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
c     &		FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/cdbn.qwave',
     &	       FILE='cdbn.qwave',
     &	       STATUS='OLD')
c        READ (10,100)
        READ (10,*)
        READ (10,*)
C...Momentum space [qgrid in MeV, uqgrid in MeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id), weight(id), uqgrid(id), wqgrid(id)
c          write(*,*) qgrid(id), weight(id), uqgrid(id), wqgrid(id)
          qgrid(id) = qgrid(id) / hc            ! MeV => 1/fm
          uqgrid(id) = uqgrid(id) * hc**1.5D0   ! MeV^-3/2 => fm^3/2
          wqgrid(id) = wqgrid(id) * hc**1.5D0
        ENDDO
 100    FORMAT (2(1X/))
 101    FORMAT (3X,D13.6,3X,D20.13,3X,D20.13)
c	PRINT *, '... CD-Bonn model read...'
	init = .TRUE.

C...Evaluate wave function
c 999	u0 = DQDVAL (q,nq,qgrid,uqgrid,.FALSE.)
c	u2 = DQDVAL (q,nq,qgrid,wqgrid,.FALSE.)
 999	CALL Pinterp (qgrid,uqgrid,nq,q,u0,dum,2)
	CALL Pinterp (qgrid,wqgrid,nq,q,u2,dum,2)

        RETURN
        END


C **********************************************************************
	SUBROUTINE WJC1 (q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-1 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~105%
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,DQDVAL,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hcM,hcG
	SAVE

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
c     &	      FILE='/u/home/wmelnitc/Work/EMC/D/Wfn/wjc-1.dat',
     &	      FILE='wjc-1.dat',
     &	      STATUS='OLD')

C...Momentum space [qgrid in MeV, ugrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
c	PRINT *, '... WJC-1 model read...'
	init = .TRUE.

C...Evaluate wavefunction
c 999	u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
c	w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
c	vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
c	vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)
 999	CALL Pinterp (qgrid,ugrid,nq,q,u,dum,1)
	CALL Pinterp (qgrid,wgrid,nq,q,w,dum,1)
	CALL Pinterp (qgrid,vtgrid,nq,q,vt,dum,1)
	CALL Pinterp (qgrid,vsgrid,nq,q,vs,dum,1)

        RETURN
        END


C **********************************************************************
	SUBROUTINE WJC2 (q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-2 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~102%
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),dum
        LOGICAL init /.FALSE./
        REAL*8  pi,hcM,hcG
	SAVE

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (init) GO TO 999     ! Data already read
C...Read data from file
	OPEN (10, FORM='FORMATTED',
     &	      FILE='wjc-2.dat',
     &	      STATUS='OLD')

C...Momentum space [qgrid in MeV, ugrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
c	PRINT *, '... WJC-2 model read...'
	init = .TRUE.

C...Evaluate wavefunction
c 999	u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
c	w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
c	vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
c	vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)
 999	CALL Pinterp (qgrid,ugrid,nq,q,u,dum,1)
	CALL Pinterp (qgrid,wgrid,nq,q,w,dum,1)
	CALL Pinterp (qgrid,vtgrid,nq,q,vt,dum,1)
	CALL Pinterp (qgrid,vsgrid,nq,q,vs,dum,1)

        RETURN
        END

**************************************************************
***       File contains various interpolation codes        ***
**************************************************************
***      Code obtained from Alberto Accardi, Dec. 2008     ***
**************************************************************
*
* II  Polynomial function interpolation of given order
      subroutine pinterp(xa,ya,n,x,y,dy,order)
*     programmer: Alberto Accardi
*     date: 2/05/01
*
*  A. COMMENTARY
*
*     Performs an interpolation using a polynomial function
*     interpolation at a given order: given an x, it uses "order" points 
*     to its left and "order" to its right to perform the interpolation
*
*     xa(*) = (DP) array with tabulated abscissae (any dimension)
*     ya(*) = (DP) array with tabulated function  (any dimension)
*     n     = (I)  number of tabulated points  
*     x     = (DP) abscissa at which compute the interpolated function
*     y     = (DP) value of the function at x
*     dy    = (DP) estimated error (usually larger than real error)
*     order = (I)  order of the interpolation (see intro)  
*                  If order = 0 performs a linear interpolation
*                  between the nearest neighbours lattice point
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      integer n, order

      double precision xa(*), ya(*), x, y, dy, tempx(n)
     :     , x1(2*order), y1(2*order), xmax, xmin, ymax, ymin  

      integer i, nlow, nmin

*
*  C. ACTION
*

      do i = 1, n
         tempx(i) = xa(i)
      end do
      call hunt(tempx,n,x,nlow)

      if (order.ge.1) then
         if (nlow.lt.order) then
            nmin = 0
         else if (nlow.le.n-order) then
            nmin = nlow-order
         else
            nmin = n-2*order
         end if
         do i = 1, 2*order
            x1(i) = xa(nmin+i) 
            y1(i) = ya(nmin+i) 
         end do
         call polintnum(x1,y1,2*order,x,y,dy)
      else
         ymax = ya(nlow+1)
         ymin = ya(nlow)
         xmax = xa(nlow+1)
         xmin = xa(nlow)
         y = ymin + (ymax-ymin)/(xmax-xmin) * (x-xmin)
      end if

      return
      end


************************************************************************
*
* III search in an ordered table 
      SUBROUTINE hunt(xx,n,x,jlo)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given an array xx(1:n) and given a value x, returns a value j
*     suchthat x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
*     either decreasing or increasing. j=0 or j=n is returned to
*     indicate that x is out of range.
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER jlo,n

      double precision x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd

*
*  C. ACTION
*

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo 
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END


************************************************************************
*
* IV  Polynomial interpolation and extrapolation
      SUBROUTINE polintnum(xa,ya,n,x,y,dy)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given arrays xa and ya of length n, and given a value x, this
*     routine returns a value y and an error estimate dy. If P(x) is the
*     polynomial of degree N-1 such that P(xa_i) = ya_i, i=1,...,n
*     then the returned value y = P(x).
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER n,NMAX
      double precision dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

*
*  C. ACTION
*

      if (n.gt.nmax) then
         print*, 'ERROR(polintnum): order larger than max', n,'>', nmax 
         stop
      end if
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
c          if(den.eq.0.)pause 'failure in polintnum'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      
      SUBROUTINE F1F2QE17(Z, A, qsq, wsq, xvalc,F1, F2)

C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS  
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c Based on the earlier code F1F2QE09 by P. Bosted.
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invariant mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus


      IMPLICIT NONE     

      REAL*8 Z, A, avgN, F1, F2, wsq, qsq, R
      REAL*8 amp/0.938272/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      real*8 xvalc(40),mcor,ecor
      integer IA, izz, izzmin, izp, izznom, izdif
      
      real*8 a1,a2,a3,b1,b2,b3,b4,b5,gd


! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

CCCCCCCCCCCCCCCCCCC  New FormFactors  CCCCCCCCCCCCCCCCCCCCCCCCC


        a1 = -1.651D0
        a2 = 1.287D0
        a3 = -0.185D0
        b1 = 9.531D0
        b2 = 0.591D0
        b3 = 0.D0
        b4 = 0.D0
        b5 = 4.994D0

        GEp = (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)


        a1 = -2.151D0
        a2 = 4.261D0
        a3 = 0.159D0
        b1 = 8.647D0
        b2 = 0.001D0
        b3 = 5.245D0
        b4 = 82.817D0
        b5 = 14.191D0

        GMp = RMUP
     &      * (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)

        GD = (1./(1 + qsq/0.71))**2
        GMn = -1.913148           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &  /(1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)

*        GMn = -1.913148*GD

        GEn = (1.700D0*tau/(1+ 3.300D0*tau))*GD

        mcor = 1.0-0.09*qsq/(0.1+qsq)
        ecor = 1.0+0.15*qsq/(0.1+qsq)
        GMn = mcor*GMn
        GEn = ecor*GEn


! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
! changed 4/09
      if(IA.eq.3) kf=0.115
      if(iA.eq.3) Es=0.015
c      if(IA.eq.3) kf=0.193
c      if(iA.eq.3) Es=0.016 
! changed 4/09
      if(IA.gt.3) kf=0.193
      if(iA.gt.3) Es=0.016
      if(IA.gt.7) kf=0.228
      if(iA.gt.7) Es=0.020 
c changed 5/09
      if(iA.gt.7) Es=0.0165
c Test MEC  10/13
c         if(iA.gt.7) Es=0.015
c      if(IA.gt.7) kf = 0.228
c      if(iA.gt.7) Es= 0.015

c Test MEC  10/13
         if(iA.gt.7) Es=0.015
         if(IA.gt.7) kf=0.228
c

c
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
! changed 5/09 
      if(iA.gt.55) Es=0.018 


! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)

      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)

      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip)) / kf

CCMEC - testing below  CCCC

      FY = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip)) / kf

      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      if(F1.LT.0.0) F1 = 0.
      if(nu.gt.0. .and.f1.gt.0.) then
        R = (F2 / nu) / (F1 / amp) * (1. + nu**2 / qsq) - 1.0
      else
        r = 0.4/qsq
      endif
      
      return
      end

      


                                                                        
      SUBROUTINE GSMEARING(Z, A, W2, Q2, xval, F1, F2, FL)

CCC   Returns Fermi smeared structure functions.  Smearing function is a Gaussian.  CCC
CCC   Note:  not tested for nuclei with A < 12 or A > 64                            CCC      

      implicit none
      real*8 Z,A,q2,w2,f1,f2,fL,xval(40)
      real*8 nu,x,mp,mp2,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta
      real*8 pf,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,rc,off_mKP_fit
      real*8 dw2dpf,r,zt,at
      real*8 xxp(100),fytot,fytot2,norm

      real*8 emcfac,emcfacL,xfr,xfl
      logical goodfit
      INTEGER ISM,drn,wfn
      external off_mKP_fit

      drn = 5
      xfr = 0.95
      xfl = 1.E-3
      wfn = 2

      mp = 0.9382727
      mp2 = mp*mp
      pi = 3.141593
      x = q2/(q2+w2-mp2)
      nu = (w2+q2-mp2)/2./mp      
      qv = sqrt(nu**2 + q2)

      If(A.EQ.3) then
c         deltae = 0.014
         deltae = 0.015
         kf = 0.115
      ELSEIf(A.GE.4) then
c         deltae = 0.014
         deltae = 0.016
         kf = 0.193
      elseif(A.GE.10.) then
        deltae = 0.015  !!! energy shift !!!
        kf = 0.228      !!! fermi momentum  !!!
      elseif(A.GE.27.) then
        deltae = 0.017
        kf = 0.238
      elseif(A.GE.56.) then
        deltae = 0.023
        kf = 0.241
      endif
c      es=deltae

      norm = 20.471
      norm = norm*2.0
c      norm = 1.0

      f1p = 0.0D0
      f1n = 0.0D0
      f2p = 0.0D0
      f2n = 0.0D0
      fLp = 0.0D0
      fLn = 0.0D0

      fytot = 0.0D0
      fytot2 = 0.0D0


! adjust pf to give right width based on kf
      pf = 0.5 * kf 
! assume this is 2 * pf * qv
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + mp) 

      DO ism = 1,99

CCC   
        xxp(ism) = -3.0+6.0*float(ism-1)/98.0

        fyuse = 1.0D0/norm*exp(-0.5*xxp(ism)*xxp(ism))    !!! Gaussian !!!

CCC  Next is from f1f209 CCC

        WSQP = W2 + XXp(ISM) * PF * DW2DPF - es * dw2des

CCC

        fytot = fytot+fyuse

c           write(6,2000) w2,q2,ism,xxp(ism),fyuse, wsqp, fytot


        IF(WSQP.GT. 1.159) THEN
          xt = q2/(q2+Wsqp-mp2)

          offshell = 1.0D0   !!!  Currently include all modification effects in EMC factors  !!!
CCC   Next is medium modification factor  CCC


          emcfac = xval(26)*(1.0-xval(27)*xt*xt)**xval(28)*
     &               (1.0-xval(29)*exp(-1.0*xval(30)*xt))

          emcfacL = (1.0 + xval(35)*xt + 
     &          xval(36)*xt*xt)/(1.0+xval(31)*log(1.0+xval(32)*q2)) 


c          at = 1.0
c          zt = 1.0
c          call f1f2in09(zt,at,q2,wsqp,xval,f1pp,f2pp,r)  !!! must turn off threading and recompile  !!!
c          fLpp =  (1.+4.*xt*xt*mp2/q2)*f2pp-2.0*xt*f1pp
c          at = 1.0
c          zt = 0.0
c          call f1f2in09(zt,at,q2,wsqp,xval,f1nn,f2nn,r)  !!! must turn off threading and recompile  !!!
c          fLnn =  (1.+4.*xt*xt*mp2/q2)*f2nn-2.0*xt*f1nn

          call sf(wsqp,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)
          f1pp = f1pp*emcfac*offshell
          f1nn = f1nn*emcfac*offshell
          fLpp = fLpp*emcfac*emcfacL*offshell
          fLnn = fLnn*emcfac*emcfacL*offshell
          f2pp = (2.*xt*f1pp+fLpp)/(1.+4.*xt*xt*mp2/q2)
          f2nn = (2.*xt*f1nn+fLnn)/(1.+4.*xt*xt*mp2/q2)

          F1p = F1p + F1pp * Fyuse
          F1n = F1n + F1nn * Fyuse
          F2p = F2p + F2pp * Fyuse
          F2n = F2n + F2nn * Fyuse      
          FLp = FLp + FLpp * Fyuse
          FLn = FLn + FLnn * Fyuse

         ENDIF

      ENDDO

      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)


c      write(6,*) fytot,fytot2

 2000 format(2f7.3,1i4,4f10.4)



      RETURN                                                            
      END                                          

      
      
      SUBROUTINE MEC2017(z,a,w2,q2,xvalc,f1corr)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine for Transverse Enhancement new the QE and Delta due to meson   CCC
CCC   exchange currents and isobar excitations in the medium.  This is assumed  CCC
CCC   to be due to quasi-deuteron 2-body currents.  Shape is a distorted        CCC
CCC   Gaussian in W^2 with a cut-off near the 2-body threshold at x=2.          CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,mp/0.938272/,mp2,w,nu
      real*8 a1,b1,c1,t1,dw2,dw22,emc,w2min,damp
      integer i
      real*8 x, f1corr, xvalc(40)

      mp2 = mp*mp
      f1corr = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu )

      if(a.lt.2.5) return

      a1 = A*q2**2*xvalc(1)*exp(-1.0*q2/xvalc(2))/
     &                                (xvalc(3)+q2)**xvalc(4)

      b1 = xvalc(5)+xvalc(6)*q2

c        c1 = xvalc(7)/(1.+xvalc(8)*sqrt(q2))

      c1 = xvalc(33)+xvalc(34)*q2

c         c1 = 0.290   !!! Test  !!!


      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2+q2*(1.-1./2.2)-1.0*mp2

      if(dw2.LT.0.0) dw2 = 0.0
      f1corr = a1*(exp(-1.*t1)*sqrt(dw2))

      f1corr = f1corr+xvalc(38)*exp(-1.*((w2-1.25)/xvalc(39))**2.0)*
     &               q2*exp(-1.*xvalc(40)*q2)


c       f1corr = f1corr

       if(f1corr.LE.1.0E-9 ) f1corr=0.0


      return
      end
      

      SUBROUTINE SF(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)
CCCC   Converts reduced cross sections to structure functions for protons and neutrons  CCCCC 

      IMPLICIT none

      real*8 w2,q2,x,sigtp,siglp,sigtn,sigln,f1p,f2p,fLp
      real*8 f1n,f2n,fLn,pi,pi2,alpha,mp,mp2
      Integer i
 

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      x = q2/(q2+w2-mp2)

   
      call rescsp(w2,q2,sigTp,sigLp)
      call rescsn(w2,q2,sigTn,sigLn)

      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f1n = sigTn/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLn = sigLn*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)
      f2n = (2.*x*f1n+fLn)/(1.+4.*mp2*x*x/q2)

      return
      end


      
      SUBROUTINE SFCROSS(w2,q2,A,Z,opt,sigt,sigl,f1,f2,fL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine to return reduced cross sections and structure functions.      CCC
CCC   opt:  0 (total), 1 (QE only), 2 (inelastic only), 3 (inelastic + MEC)     CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      IMPLICIT none
 
      real*8 w2,q2,A,Z,sigt,sigl,f1,f2,fl
      real*8 x,alpha,pi,pi2,mp,mp2,xvalc(40)
      integer lstart,i, opt
      real*8 xvalt(160) / 
     & 0.40249E+00,0.70000E+01,0.33195E+00,0.66463E+01,0.80000E+00,
     & 0.14605E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.97866E+00,0.96776E+00,0.10363E+01,0.99254E+00,0.98301E+00,
     & 0.10000E+01,0.97634E+00,0.10055E+01,0.98454E+00,0.99203E+00,
     & 0.99209E+00,0.99553E+00,0.10000E+01,0.10125E+01,0.10060E+01,
     & 0.10664E+01,0.42389E+00,0.14867E+01,0.13360E+00,0.47378E+02,
     & 0.90026E+01,0.29186E-02,0.18698E+00,0.49604E-01,-.16739E+00,
     & 0.16172E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+00,
     & 0.35826E+00,0.60000E+01,0.34200E+00,0.67220E+01,0.80000E+00,
     & 0.11683E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.10109E+01,0.10945E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10009E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.99738E+00,0.10000E+01,0.10086E+01,0.10000E+01,
     & 0.10780E+01,0.44725E+00,0.13428E+01,0.13512E+00,0.44000E+02,
     & 0.99966E+01,0.93603E-02,0.20825E+00,0.62977E-01,-.21674E+00,
     & 0.16388E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+00,
     & 0.35274E+00,0.60000E+01,0.33300E+00,0.66500E+01,0.80000E+00,
     & 0.13096E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.10061E+01,0.99196E+00,0.10134E+01,0.10157E+01,0.10000E+01,
     & 0.10000E+01,0.97920E+00,0.10080E+01,0.98964E+00,0.10000E+01,
     & 0.10000E+01,0.99570E+00,0.10000E+01,0.10252E+01,0.10000E+01,
     & 0.10780E+01,0.80008E-01,0.91568E+01,0.41462E+00,0.44000E+02,
     & 0.20200E-01,0.99597E+00,0.27976E+00,0.22938E-07,-.23486E+00,
     & 0.16106E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+00,
     & 0.35269E+00,0.60000E+01,0.33300E+00,0.66500E+01,0.80000E+00,
     & 0.13100E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10031E+01,
     & 0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     & 0.10780E+01,0.83622E-01,0.91410E+01,0.40578E+00,0.44000E+02,
     & 0.51500E-02,0.99373E+00,0.27989E+00,0.26265E-07,-.25044E+00,
     & 0.15755E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+00 /   
      
      LOGICAL GOODFIT/.true./  

      mp = .938272
      mp2 = mp*mp 

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi
      x = q2/(q2+w2-mp2)
      
      if(A.GE.3.) then
         lstart = 0  
      elseif(A.GE.27.) then
        lstart = 40
      elseif(A.GE.56.) then   
        lstart = 80
      elseif(A.GE.64.) then
        lstart = 120
      endif 
      do i=1,40
         xvalc(i) = xvalt(i+lstart)
      enddo   
      if(A.LE.4) xvalt(4) = 6.0
      
c      write(6,*) xvalc
      call csfit(w2,q2,A,Z,xvalc,opt,sigt,sigL)

c      write(6,*) w2,q2,sigt,sigL
      
      f1 = sigt
      fL = sigl*2.*x
      sigt = 0.3894e3*8.0d0*pi2*alpha/abs(w2-mp2)*sigt
      sigL =  0.3894e3*8.0d0*pi2*alpha/abs(w2-mp2)*sigL
      
      return
      end

      
      SUBROUTINE CSFIT(w2,q2,A,Z,xvalc,opt,sigt,sigL)
      IMPLICIT none

      real*8 e,ep,th,q2,w2,x,cs,flux,kappa,sin2,tan2,csmod
      real*8 f1,fl,f1qe,flqe,r,rqe,sigt,sigl,f1mec,sigm,f2,f2qe,f2mec
      real*8 alpha,pi,pi2,mp,mp2,res,veff,foc,Z,A,xvalc(40)
      real*8 psip,fy1,fy2,int1,int2,rat,f1t,f2t,fLt
      integer i,j,k,ntot,opt
      LOGICAL GOODFIT/.true./  
      character*40 filename
      mp = .938272
      mp2 = mp*mp

      alpha = 1./137.036
      pi = 3.141593
      pi2 = pi*pi

      x = q2/abs(w2-mp2+q2)
 

      int1 = 0.0D0
      int2 = 0.0D0
      rat = 1.0D0

CCC  NEXT bit only needed if fitting scaling function CCC
      do i=1,120        
        psip = -2.0+0.06*i
       FY1 = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
       FY2 = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip))
       int1 = int1+fy1  
       int2 = int2+fy2
      enddo
      rat= int1/int2
CCC

      call gsmearing(Z,A,w2,q2,xvalc,f1,f2,fL)

      r = fL/2.0D0/x/f1

c      write(6,*) w2,q2,f1,f2,fL
      
      call f1f2qe17(Z,A,q2,w2,xvalc,f1qe,f2qe)

c      write(6,*) "17:  ", w2,q2,f1qe,f2qe     
c      call f1f2qe16(Z,A,q2,w2,xvalc,f1qe,f2qe)
c      write(6,*) "16:  ", w2,q2,f1qe,f2qe
            
      call MEC2017(Z,A,w2,q2,xvalc,f1mec)
c      write(6,*) "17:  ", w2,q2,f1mec     
c      call MEC2016(Z,A,w2,q2,xvalc,f1mec)
c      write(6,*) "16:  ", w2,q2,f1mec
      
      
      f1qe = f1qe*rat
      f2qe = f2qe*rat
c      write(6,*) "!!!  ",f1,f1qe,f1mec
      if(opt.EQ.1) then       !!! QE only
         f1 = 0.0
         f2 = 0.0
         fL = 0.0
         f1mec = 0.0
      elseif(opt.EQ.2) then   !!! inelastic only
         f1qe = 0.0
         f2qe = 0.0
         f1mec = 0.0
      elseif(opt.EQ.3) then   !!! inelastic + MEC
         f1qe = 0.0
         f2qe = 0.0
      endif
      
      f1 = f1 + f1qe + f1mec
      f2mec = 2.*x*f1mec/(1.+4.*x*x*mp2/q2) 
      f2 = f2 + f2qe + f2mec
      fL = (1.+4.*x*x*mp2/q2)*f2-2.0*x*f1
c      if(F1.LE.0.0) write(6,*) "F1 < 0"
      if(F1.LE.0.0) F1=0.0
      if(F2.LE.0.0) F2=0.0
      if(FL.LE.0.0) FL=0.0
      sigt = f1
      sigl = fL/2./x


 2000 format(6f10.4)
          
      return

      end
      



      






