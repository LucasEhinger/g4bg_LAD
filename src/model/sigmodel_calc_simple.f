        FUNCTION sigmodel_calc(e1pass,e2pass,thpass,zpass,apass,mpass,xflag,model)
C       +______________________________________________________________________________
c	
C       Calculate cross section using Peter's F1F209.f routine
c       
c       ARGUMENTS:
c       
c       E1:		-	Incident energy in GeV.
c       E2:		- Scattered energy in GeV.
c       TH:		- Scattering angle in Degrees.
c       A:		- 'A' of nucleus.
c       Z:		- Number of protons in nucleus.
c       M_TGT:	- Mass of target nucleus in GeV/c2.
c       M_REC:	- Mass of recoiling nucleon in GeV/c2.
c       E_SEP:	- Separation energy for target nucleus in GeV/c2.
c       SIG  :	- Calculated cross section in nb/(MeV-ster).// should be nb/(GeV-ster)
C       ______________________________________________________________________________

        implicit none
	include 'math_physics.inc'

C       Declare arguments.
        real sigmodel_calc
	real*8		e1,e2,th,a,z,m_tgt
	real*4          e1pass,e2pass,thpass,mpass
	real*4          sig,factpass,sig_dis_pass,sig_qe_pass
	integer         zpass,apass
	logical		modelhack

C       Declare locals.

	real*8	 	sig_qe,sig_dis,y,normfac,fact
        real*8		thr,cs,sn,tn,elastic_peak
	real*8          Q2,nu,WSQ, x
	real*8          F1,F2,FL,W1,W2,sigmott,r
	real*8          W1p,W1n,W1D,W2p,W2n,W2D
	real*8          inelastic_it
	integer         xflag !flag for which xsec to calculate 1=both 2=QE only 3=DIS only
	logical         first
        integer         sigdis_model,model  !DIS_MODEL defined in TARG
        real*8          eps,f1d,f2d,fLd,f2dqe,f1dqe,fLdqe,sigt,sigl
        real*8          sigtp,siglp
        real*8          kappa,flux
        integer         wfn,opt
        logical         doqe,dfirst/.false./
        real*8          pi2,alp

	save

        real*8 emc_func_xem
	external emc_func_xem

	real*8 emc_func_slac
	external emc_func_slac

	data first/.true./

	e1=dble(e1pass)
	e2=dble(e2pass)
	th=dble(thpass)
	a=dble(apass)
	z=dble(zpass)
	m_tgt=dble(mpass)
        sigdis_model=model
        alp = 1./137.036
        pi2 = 3.14159*3.14159


	sig =0.0
	sig_qe=0.0
	sig_dis=0.0

C       If at or beyond elastic peak, return ZERO!

	thr = th*d_r
	cs = cos(thr/2.)
	sn = sin(thr/2.)
	tn = tan(thr/2.)
	elastic_peak = e1/(1.+2.*e1*sn**2/m_tgt)
      	if (e2.ge.elastic_peak) then
       	   sig = 0.0
       	   return
       	endif
c	write(6,*) 'got to 1'

	Q2 = 4.*e1*e2*sn**2
	nu=e1-e2
	WSQ = -Q2 + m_p**2 + 2.0*m_p*nu 
        x = Q2/2/m_p/nu
        eps = 1.0/(1+2.*(1+Q2/(4*m_p**2*x**2))*tn**2)
        kappa = abs((wsq-m_p**2))/2./m_p
        flux = alp*kappa/(2.*pi2*Q2)*e2/e1/(1.-eps)
        wfn=2

	F1=0
	F2=0
	r=0
	if((xflag.eq.1).or.(xflag.eq.3)) then
c----------------------------------------------------------------
c       
c       do inelastic stuff
c	   call F1F2IN09(Z, A, Q2, WSQ, F1, F2, r)
C Use old Bodek fit + SLAC EMC fit for now, b/c F1F2IN09 doesn't like large Q2,W2
          if(sigdis_model .eq. 1) then
	     if(wsq.gt.1.1664) then
	      call ineft(Q2,sqrt(wsq),W1p,W2p,dble(1.0))
	      call ineft(Q2,sqrt(wsq),W1D,W2D,dble(2.0))

	      W1n=2.0*W1D-W1p
	      W2n=2.0*W2D-W2p

	      W1=Z*W1p+(A-Z)*W1n
	      W2=Z*W2p+(A-Z)*W2n
              sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
              sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2)
              sig_dis = sig_dis*emc_func_slac(x, A)
             endif
          endif
  
          if(sigdis_model .eq. 2) then
C       Mott cross section
c          call gsmearing(Z,A,WSQ,Q2,F1,F2)
            if(A.gt.2)then
              opt=3
              call SFCROSS(WSQ,Q2,A,Z,opt,sigt,sigl,f1,f2,fL)
              sig_dis = 1000*flux*(sigt+eps*sigl)
            endif
            if(A.EQ.1.0)then
              call RESCSP(wsq,q2,sigtp,siglp)
              sig_dis = 1000*flux*(sigtp+eps*siglp)
            endif
C       Convert F1,F2 to W1,W2
c           W1 = F1/m_p
c           W2 = F2/nu
c          sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
c          sig_dis = 1d3*sigmott*(W2+2.0*W1*tn**2)
          endif 

          if(sigdis_model .eq. 3) then
             if((Z.eq.1.) .and. (A.eq.2.)) then
                doqe = .false.
                wfn=2
                call RESCSD(WSQ,Q2,eps,doqe,f1d,f2d,fLd,wfn,sig_dis)
                sig_dis=1000*flux*sig_dis

             else 
               write(6,*) 'Wrong DIS model !! This model is for '//
     >         'Deuterium only'
               sig_dis=0.0
             endif
          endif
CDG apply "iteration" correction (used for XEM analysis)
CDG DO not use this for more "generic" stuff.
CDG        sig_dis = sig_dis*inelastic_it(x,A)
c        else
c             W1=0.0
c             W2=0.0
	endif


	if((xflag.eq.1).or.(xflag.eq.2)) then
           if((Z.eq.1.) .and. (A.eq.2.)) then
               dfirst = .false.
               call SQESUB(WSQ,Q2,wfn,f2dqe,f1dqe,fLdqe,dfirst)
               sigt = 0.3894e3*f1dqe*pi2*alp*8.0/abs(wsq-m_p**2)
               sigl = 0.3894e3*fLdqe*pi2*alp*8.0/abs(wsq-m_p**2)/2.*abs(wsq-m_p**2+q2)/q2
               sig_qe = sigt+eps*sigl
               sig_qe=1000*flux*sig_qe

           else
               opt=1
               call SFCROSS(WSQ,Q2,A,Z,opt,sigt,sigl,f1,f2,fL)
               sig_qe = sigt+eps*sigl
               sig_qe = 1000*flux*sig_qe
           endif
c	   call F1F2QE09(Z, A, Q2, WSQ, F1, F2)
C       Convert F1,F2 to W1,W2
c	   W1 = F1/m_p
c	   W2 = F2/nu
C       Mott cross section
c	   sigmott=(19732.0/(2.0*137.0388*e1*sn**2))**2*cs**2/1.d6
c	   sig_qe = 1d3*sigmott*(W2+2.0*W1*tn**2)
C Temp test - DJG May 23, 2013
c	   sig_qe=sig_qe/0.8
	endif


	sig = sig_qe + sig_dis !sig is already real*4


	sig_qe_pass = sig_qe ! pass back as real*4
	sig_dis_pass = sig_dis

	return
	end


	real*8 function inelastic_it(x,A)
C DJG: Correction to inelastic cross section to F1F209 from XEM
C DJG: 40 degree data. Just a simple one-pass iteration for use
C DJG: to check our model dependence.
	real*8 A,x
	real*8 p1,p2,p3,p4
	real*8 x1,x2,xit


	inelastic_it = 1.0 ! set to 1 by default

	if(A.gt.1.5.and.A.lt.2.5) then
	   x1=0.3172d0
	   x2=0.7275d0
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=1.2394d0 - 2.1805d0*xit + 5.6853d0*xit**2
     >     -4.3908d0*xit**3
	endif

	if(A.gt.2.5 .and. A.lt.3.5) then
	   x1=0.3172
	   x2=1.055
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=2.9235d0 - 16.075d0*xit + 46.426d0*xit**2
     >     -56.779d0*xit**3 + 25.007d0*xit**4
	endif

	if(A.gt.3.5 .and. A.lt.4.5) then
	   x1=0.3172
	   x2=1.0927
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=1.505d0 - 4.8103d0*xit + 15.221d0*xit**2
     >     -19.713d0*xit**3 + 8.9693d0*xit**4
	endif

	if(A.gt.8.5 .and. A.lt.9.5) then
	   x1=0.6478
	   x2=1.0929
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.89139d0 + 0.44009d0*xit + -0.44163d0*xit**2
	endif

	if(A.gt.11.5 .and. A.lt.13.5) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.88964d0 + 0.39884d0*xit - 0.36051d0*xit**2
	endif

	if(A.gt.61.0 .and. A.lt.66.0) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.77646d0 + 0.90481d0*xit - 0.83815d0*xit**2
	endif


	if(A.gt.196.0 .and. A.lt.198.0) then
	   x1=0.3172
	   x2=1.005
	   if(x.lt.x1) then
	      xit=x1
	   elseif(x.gt.x2) then
	      xit=x2
	   elseif(x.ge.x1 .and. x.le.x2) then
	      xit=x
	   endif
	   inelastic_it=0.70439d0 + 1.0510d0*xit - 0.91679d0*xit**2
	endif

	return
	end

	   
c-------------------------------------------------------------------------------------------
	real*8 function emc_func_slac(x,A)
	real*8 x,A,atemp
	real*8 alpha,C

	atemp = A
!	if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these 2...
!	   atemp = 12
!	endif

	alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
	1    -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
	2    +775.767*x**7 - 205.872*x**8

	C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
	
	emc_func_slac = C*atemp**alpha
	return 
	end      
	   

	
