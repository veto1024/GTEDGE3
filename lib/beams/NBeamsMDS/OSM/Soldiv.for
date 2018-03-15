      PROGRAM SOLDIV
C		SOL-DIVERTOR MODEL MAIN PROGRAM
	 USE MSIMSL 
	INCLUDE 'SOLDIV.FI' 
													
 	parameter (nit=9,npsi0=22,nmesh=25,nsol=51,nn=30)
	dimension soln(nit),solv(nit),sole(nit),
	1          gam25con(nn),gamQ25con(nn),gamorbl(nn),gamQorbl(nn),
	2	      gamsum(nn), gamQsum(nn), gamiller(nn),gamQmiller(nn),
	3	      fheat(nsol),fpart(nsol),gamvorbl(nn),vloss(npsi0,nit),
     4          eloss(npsi0,nit),lo(npsi0,nit),floss(npsi0,nit),   
     5          v0m(nit,nit),v0p(nit,nit),xke0m(nit,nit),xke0p(nit,nit),
	6	      angle0(nit),angleD(nit),epminhat(npsi0),vlosshat(npsi0),
     7		  epmin(npsi0,nit),emin(npsi0,nit),
     8	      cos0(npsi0),xkeDM(nit,nit),
     9		  xkeDp(nit,nit),flos(nit),elos(nit),vlos(nit),
	1	      numloss(npsi0),eminhat(npsi0),flosshat(npsi0),
	2		  elosshat(npsi0),xmomiol(nmesh),forbl2(nmesh),	
     4	      xmorbl2(nmesh),eorbl2(nmesh)
      real nbdep1(51),nbdep2(51),nbdep3(51),floss,vloss,eloss,ephi,
	1     forbl,xmorbl,eorbl,ffrac1(51),retcur(51),
	2     ffrac2(51),ffrac3(51),lost1(51), lost2(51),lost3(51),
     3     tally(51),erex0(51),ephi0(51),erextot(51),fpsi0(51),
	4     nbdepkept1(50),nbdepkept2(50),nbdepkept3(50),nbdeplost1(50),
     5     nbdeplost2(50),nbdeplost3(50),NBIlost(51),NBIreturn   
	integer rholength, psinum,nmass,npsi,n0,ioptFIOL

	CALL SOLDATA(nbdep1,nbdep2,nbdep3,erextot,fpsi0)

      xxx=epot(24)

	fuelplout = fuelplin
	thetain = theta
	N1 = 0
	DELN = 0.02
	CALL GEOMETRY(N1)

C		INITIALIZE
	AP = 39.44*RMAJOR*AMINOR*SQRT((1.+ELONG**2)/2.)
c	ap= 39.44*rmajor*aminor*sqrt(elong)
	VP = 19.72*RMAJOR*(AMINOR**2)*ELONG
	ITERTDL = 0
	yyt = fuelplout
	iyyt=ioptq95 
	IMARFE = 0
	TISOL = TSOL
	x = yltebarx
	x = fheate
	x = xlpm(20)
	TESOL = TSOL
	TIPED = TPED
	TEPED = TPED
	XNBAR = XNPED
	TIBAR = TPED
	TEBAR = TPED
	
	xxx = xlne(20)

	XLNBAR = 0.05
	XLTIBAR = 0.06
	XLTEBAR = 0.05
	VSEP = 1.E6
	TSPL = TSEP
	TPL = TSEP
	TPF = TSEP
	TDIV = TSEP
	TEMPSEP = TSEP
	SVIOND = 1.E-14
	SVTOTD = 1.E-14
	DEL = 1.E-2
	DELN = 2.*DEL
	DELNV = DELN
	DELEA = DELN/7.
	DQIN = FLUXHEAT*XLPERP/DELEA
	DELNIN = FLUXPART*XLPERP/DELN
	GAMOUTSPL = 1E20
	GAMOUTPF  = 1E20
	GAMOUTPL  = 1E20
	GAMZERO = 1E20
	XNAV = 1.E20
	XNEL = XNAV
 	CSD = SQRT(2.*XK*TD/XMASS)
	CALL DISTR
	AAA = xnsepex
	wall = 0.609


c	 constants
  	RLOSSIOL = 1.0
c	 Toggle for neutral beam ion orbit loss (ioptFIOL: 1=yes,0=no)
	ioptFIOL = 1
 	ioptionorb = 1
	ioptxtran = 0
 	xm = 3.343e-27
	eq = 1.6e-19
      Icur = plasmacur
	Btor = Bphi
	radwall=   0.90
    	radminor = aminor*sqrt(0.5*(1+elong**2))
	rminorw = wall*sqrt(0.5*(1+elong**2))
	dr = 0.1
	RM = rmajor
	xmu0 =1.257					 
      zmass = 3.34e-27
	delnaf = 0.02
	frholength = 50
	nbi = 75000
55	OPEN(unit=121,FILE='orbloss.TXT',STATUS='UNKNOWN')
c	do loop over launch surface radii,
c      nmass = 1 (fast ions), nmass = 2 (deuterium), nmass = 3 (carbon)

	do 2222 nmass = 1,3

	rholength = 24
	if(nmass.eq.1)then
	rholength = frholength
	endif

	DO 2221 n0 = 1, rholength

	delna0 = delna
	do 3 m = 1,rholength+1
	erex0(m) = erex(m)
3	continue

	if(nmass.eq.1)then
	   delna0 = delnaf*aminor
	   do 4 m = 1,rholength+1
	      erex0(m) = erextot(m)
4	   continue
	lost1(n0) = 0.0
      lost2(n0) = 0.0
	lost3(n0) = 0.0
	tally(n0) = 0.0
      endif 
 
	xdr =  (rholength+1-n0)
	rminor0	= radminor - xdr*delna0

c	set loss surface radius	
	rminorD = radminor
	if(nmass.eq.2)then
	rminorD = rminorw
	endif

	nD = rholength+1
c	xloss input
	thetax = 4.712
	deltathetx = 0.15			

	deltarx = radminor*deltathetx
	rminorx = radminor - deltarx
	nx = 25
	do 5 m=1,rholength
	if(rhor(25-m).lt.rminorx) goto 5
     	nx=25-m
5	continue 

	xpi = 3.1416
	xy = xmu0*Icur/6.283


c	Calculate electrostatic potential from exp Erad	in edge
c	Calculate minimum x-loss energy as av Erad between FS & Sep
	ephi0(rholength+1) = 0.0
	ersum = erex0(rholength+1)

	do 10 n =1,rholength
	m = rholength+1-n
	ephi0(m) = ephi0(m+1) + radminor*delna0*0.5*(erex0(m)+erex0(m+1))
      ersum =	 ersum + erex0(m)
	erav(m) = ersum/(n+1)
	if(nmass.eq.2)then
	ephi(m) = ephi0(m)
	endif
 10	continue
    
      if(nmass.eq.2)then
	ephi0(rholength) = 0.0
	endif

c	4-26-11 formulation

	xxI = xmu0*Icur/(6.1832*(radminor**2))	
	fphiD = 1./sqrt(1. + (xxI**2)*(rminorD**2)/(Btor**2))
 	fphi0 = 1./sqrt(1. + (xxI**2)*(rminor0**2)/(Btor**2))

c 	fphix = 1./sqrt(1. + (xxI**2)*(rxloss**2)/(Btor**2))

	Tion = xti(n0)
c	 dummy variable tion for fast ion scenario (no need for ep)
	if(nmass.eq.1)then
	Tion = 1
	endif

c	calc loss for different direction cosines thru 200, for launch surfance n0
c     only do one calculation for fast ions - all oriented at the same psi0
c      22 different psi0 locations for deuterium and carbon
	psi0 = 0.01
	psinum = 22      
	if(nmass.eq.1)then
	psi0 = fpsi0(n0)
	psinum = 1 
	endif
			
	do 200 npsi = 1,psinum
	if(npsi.eq.2) psi0 = 0.1
	if(npsi.gt.2.and.npsi.lt.11) psi0 = psi0 + 0.1
	if(npsi.eq.11) psi0 = 0.99

	if(npsi.eq.12) psi0 = -0.01
	if(npsi.eq.13) psi0 = -0.1
	if(npsi.gt.13.and.npsi.lt.22) psi0 = psi0 - 0.1
	if(npsi.eq.22) psi0 = -0.99

	Wtrap = rmajor*erex0(nD)/(1.0+psi0**2) 

	theta0 = 0.0
	thetaD = 0.0
	

	do 99 n = 1,9 
	emin(npsi,n) = 0.0
	do  90 m = 1,9 	
	h0 = 1. + (rminor0/RM)*cos(theta0)
 	hd = 1. + (rminorD/RM)*cos(thetaD)

    	delphi = ephi0(n0) - ephi0(nD)
	xpot = 2.*(eq/xm)*delphi
      xflux= xxI*(eq/(1.*xm*hd*fphiD))*(rminor0**2-rminorD**2)
c****charge and mass enter as e/m in xpot and xflux***************	
	a =(h0/hd)*(1.-psi0**2)-(((h0/hd)*(fphi0/fphid)*psi0)**2) - 1.
	bterm = 2.*xflux*(h0*fphi0)*psi0
	c = (xflux**2) - xpot
	abc = 4.*a*c/(bterm**2)
	if(abc.gt.1.0) goto 80
	v0p(n,m) = -1.*(bterm/(2.*a))*(1.0+sqrt(1.-4.*a*c/(bterm**2)))
	v0m(n,m) = -1.*(bterm/(2.*a))*(1.0-sqrt(1.-4.*a*c/(bterm**2)))
	goto 85
80	continue	
    	v0p(n,m) = 0.0
	v0m(n,m) = 0.0
c***but min reduced energy depends on mass**********************
85	xke0m(n,m) = 0.5*xm*(v0m(n,m)**2)/eq
	xke0p(n,m) = 0.5*xm*(v0p(n,m)**2)/eq
	xkeDm(n,m) = xke0m(n,m) + (ephi0(n0)-ephi0(nD))
	xkeDp(n,m) = xke0p(n,m) + (ephi0(n0)-ephi0(nD))
	thetaD = thetaD + 3.14159/4.
90	continue
	thetaD = 0.0
	theta0 = theta0 + 3.14159/4.
99	continue
	angle0(1) = 0.0
	angleD(1) = 0.0
	do 150 n=2,9
	angle0(n) = angle0(n-1) + 3.14159/4.
	angled(n) = angled(n-1) + 3.14159/4.
150	continue 

c	Minimum Loss Energy for each Theta0
c	Assumes that for each launch theta0 and npsi, the thetaD with the smallest  
c	minimum Loss Energy defines the lower energy of the loss cone at that theta0
c	emin(npsi,theta0)
      if(nmass.eq.1) goto 1560 
	if(npsi.ge.12) goto 1560
	
	
	Do 155 n=1,8
	emin(npsi,n) = 1.0e9
	
	do 154 m=1,8
	if(v0m(n,m).lt.0.0) goto 154 
	if(emin(npsi,n).gt.xke0m(n,m).and.v0m(n,m).gt.0.0) then
		emin(npsi,n) = xke0m(n,m)
		LO(npsi,n) = m
	endif
154	continue
155	continue
 	  
	       
	Do 156 n=1,8
 	epmin(npsi,n) = emin(npsi,n)/(Tion)
156	continue

	goto 160
c	co-current ions
1560	Do 158 n=1,8
	emin(npsi,n) = 1.0e9
c 	if(v0p(n,m).lt.0.0) emin(npsi,n) = xke0p(n,1)
	do 157 m=1,8
	if(v0m(n,m).gt.0.0) goto 157 
	if(emin(npsi,n).gt.xke0p(n,m).and.v0p(n,m).gt.0.0) then 
		emin(npsi,n) = xke0p(n,m)
		LO(npsi,n) = m
	endif
157	continue
158	continue

      if(nmass.eq.1)then
      do 500 n=1,8
	    if(emin(npsi,n)<nbi)then
	       lost1(n0) = 1.0+lost1(n0)
	    endif
	    if(emin(npsi,n)<nbi/2)then
	       lost2(n0) = 1.0+lost2(n0)
	    endif
	    if(emin(npsi,n)<nbi/3)then
	       lost3(n0) = 1.0+lost3(n0)
	    endif
	tally(n0) = 1.0+tally(n0)
500	continue
      endif	 

	Do 159 n=1,8
	epmin(npsi,n) = emin(npsi,n)/(Tion)
159	continue

160  	if(nmass.gt.1)then
c	Evaluate Ion, Velocity and Energy Orbit Loss, for each pitch angle & theta0          	
      do 161 n = 1,8 
	A = 1.5
	X = EPMIN(npsi,n)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	floss(npsi,n) = (value1 - value)/value1

	A = 2.0
	X = EPMIN(npsi,n)
	VALUE = GAMI(A,X)
	VALUEv= GAMMA(A)
	vloss(npsi,n) = (valuev - value)/valuev

	A = 2.5
	EVALUE = GAMI(A,X)
	EVALUE1= GAMMA(A)
	eloss(npsi,n) = (evalue1 - evalue)/evalue1
161	continue

	avfloss = 0.
	avvloss = 0
	aveloss = 0
	do 165 n = 1,8
	avfloss = avfloss + floss(npsi,n)/8.
	avvloss = avvloss + vloss(npsi,n)/8. 
	aveloss = aveloss + eloss(npsi,n)/8.
165	continue
	

c     ************2nd option  3/14/12****************************************************
c	assumes that for each r0 and npsi, the ions will sample every theta0 on the flux surface 
c	and will be lost from the theta0 with the lowest emin(r0,theta0,npsi) 

c	3/23/12**************I think the 1st option is better***************************
	eminhat(npsi) = 1.e4
	numloss(npsi) = 8
	do 1599 n=1,8
	if(emin(npsi,n).lt.eminhat(npsi)) then
	eminhat(npsi) = emin(npsi,n)
	numloss(npsi) = n
c	numloss is the index of theta0 for which emin(npsi,n) is least
	endif
	epminhat(npsi) = eminhat(npsi)/Tion   
1599	continue
1600   do 1610 n = 1,8 
	A = 1.5
	X = EPMINhat(npsi)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	flosshat(npsi) = (value1 - value)/value1

	A = 2.0
	X = EPMINhat(npsi)
	VALUE = GAMI(A,X)
	VALUEv= GAMMA(A)
	vlosshat(npsi) = (valuev - value)/valuev

	A = 2.5
	EVALUE = GAMI(A,X)
	EVALUE1= GAMMA(A)
	elosshat(npsi) = (evalue1 - evalue)/evalue1
1610	continue

c*****************end 2nd option*************************************************
	endif

200	continue

		
c	This completes calculation of Loss Fractions for each Theta0 for all pitch angles


c	Integrate over cosine of pitch angle to Calculate Loss Fractions for each Theta0 
	cos0(1) = 0.975
	cos0(2) = 0.90
	do 215 npsi = 3,10																		  
	cos0(npsi) = cos0(npsi-1) - 0.1
215	continue
	cos0(11) = 0.025 

	cos0(12) = -0.975
	cos0(13) = -0.90
	do 216 npsi = 14,21
	cos0(npsi) = cos0(npsi-1) + 0.1
216	continue
	cos0(22) = -0.025 

      Do 225 n=1,8 

	Flos(n) = 0.0
	Vlos(n) = 0.0
 
	Elos(n) = 0.0
	Flos(n) = 0.05*(floss(1,n)+floss(12,n))
	Vlos(n) = 0.05*(vloss(1,n)*cos0(1)+vloss(12,n)*cos0(12))

 	 
	Elos(n) = 0.05*(eloss(1,n)+eloss(12,n))
	Do 220 npsi= 2,10
	Flos(n) = Flos(n) + 0.1*(floss(npsi,n)+floss(npsi+11,n))
	Vlos(n) = Vlos(n) + 0.1*(vloss(npsi,n)*cos0(npsi) + 
     1	vloss(npsi+11,n)*cos0(npsi+11))


	Elos(n) = Elos(n) + 0.1*(eloss(npsi,n)+eloss(npsi+11,n))
220	continue
	Flos(n) = Flos(n) + 0.05*(floss(11,n)+floss(22,n))
	Vlos(n) = Vlos(n) + 0.05*(vloss(11,n)*cos0(11) + 
     1	vloss(22,n)*cos0(22)) 

	Elos(n) = Elos(n) + 0.05*(eloss(11,n)+eloss(22,n))

	Flos(n) = 0.5*Flos(n)
	Vlos(n) = 0.5*Vlos(n)
	Elos(n) = 0.5*Elos(n)
	
225	continue

c	calculate av loss from flux surface assuming poloidal symmetry of fluxes across FS	 
	flosav = 0.0
	vlosav = 0.0
	elosav = 0.0
	do 226 n = 1,8
	flosav = flosav + flos(n)/8.
	vlosav = vlosav + vlos(n)/8.
	elosav = elosav + elos(n)/8.
226	continue

	if(nmass.gt.1)then
	if(nmass.eq.2) then
	yy1(n0) = RLOSSIOL*vlosav
	endif

	FORBL(N0) = RLOSSIOL*FLOSAV
	XMORBL(N0) = RLOSSIOL*VLOSAV
c  ************yy1 is filling in for xmorbl, which has transfer problems**************
	
	EORBL(N0) = RLOSSIOL*ELOSAV

	if(nmass.eq.3) then 
	yy2(n0) = RLOSSIOL*vlosav 
 	endif
c************option 2 integral over cosine****************************
	
	Flos2 = 0.0
	Vlos2 = 0.0
	Elos2 = 0.0
	Flos2 = 0.05*(flosshat(1)+flosshat(12))
	Vlos2 = 0.05*(vlosshat(1)*cos0(1)+vlosshat(12)*cos0(12))	 
	Elos2 = 0.05*(elosshat(1)+elosshat(12))
	Do 2200 npsi= 2,10
	Flos2 = Flos2 + 0.1*(flosshat(npsi)+flosshat(npsi+11))
	Vlos2 = Vlos2 + 0.1*(vlosshat(npsi)*cos0(npsi) + 
     1	vlosshat(npsi+11)*cos0(npsi+11))
	Elos2 = Elos2 + 0.1*(elosshat(npsi)+elosshat(npsi+11))
2200	continue
	Flos2 = Flos2 + 0.05*(flosshat(11)+flosshat(22))
	Vlos(n) = Vlos(n) + 0.05*(vlosshat(11)*cos0(11) + 
     1	vlosshat(22)*cos0(22)) 
	Elos2 = Elos2 + 0.05*(elosshat(11)+elosshat(22))

	Flos2 = 0.5*Flos2
	Vlos2 = 0.5*Vlos2						  
	Elos2 = 0.5*Elos2
	FORBL2(N0) = RLOSSIOL*FLOS2
 	XMORBL2(N0) = RLOSSIOL*VLOS2
	EORBL2(N0) = RLOSSIOL*ELOS2
	endif
c*******************************************************************

      if(nmass.eq.1)then
      ffrac1(n0) = lost1(n0)/tally(n0)
      ffrac2(n0) = lost2(n0)/tally(n0)
      ffrac3(n0) = lost3(n0)/tally(n0)
  	endif

2221 	continue

2222	CONTINUE
c	Reduce NBI source rate by fast ion loss [ions/cm^3-s]
c     If no NBI IOL, then just scaled deposition profiles
	do 501 n0=1,50
	nbdepkept1(n0) = 0.624E25*pbeam/(nbi)*nbdep1(n0)
 	nbdepkept2(n0) = 0.624E25*pbeam/(nbi/2)*nbdep2(n0)
	nbdepkept3(n0) = 0.624E25*pbeam/(nbi/3)*nbdep3(n0)
	nbdeplost1(n0) = 0
	nbdeplost2(n0) = 0
	nbdeplost3(n0) = 0
	if(ioptFIOL.eq.1)then   
	nbdepkept1(n0) = 0.624E25*pbeam/nbi*(1-ffrac1(n0))*nbdep1(n0)
	nbdeplost1(n0) = 0.624E25*pbeam/nbi*ffrac1(n0)*nbdep1(n0)
	nbdepkept2(n0) = 0.624E25*pbeam/(nbi/2)*(1-ffrac2(n0))*nbdep2(n0)
	nbdeplost2(n0) = 0.624E25*pbeam/(nbi/2)*ffrac2(n0)*nbdep2(n0)
	nbdepkept3(n0) = 0.624E25*pbeam/(nbi/3)*(1-ffrac3(n0))*nbdep3(n0)
	nbdeplost3(n0) = 0.624E25*pbeam/(nbi/3)*ffrac3(n0)*nbdep3(n0)
	endif
501	continue		
c     Subtract all fast ions that are lost from total source rate to
c     core before attenuation (calculated in fuel.for)
	do 502 n0=1,50
	spellet = spelletreal + nbdepkept1(n0)*0.02+nbdepkept2(n0)*0.02
	1       +nbdepkept3(n0)*0.02
	fuelmp = fuelmpreal + nbdeplost1(n0)*0.02+nbdeplost2(n0)*0.02
	1       +nbdeplost3(n0)*0.02
	NBIlost(n0) = nbdeplost1(n0)*aminor/VP+nbdeplost2(n0)
	1              *aminor/VP +nbdeplost3(n0)*aminor/VP
	retcur(n0) =  xk*NBIlost(n0)
502	continue

	write(*,*) xk
c	 Project lost particles from NBI onto shorter rho vector
c      Note this loss continues into the plasma farther than rho=0.864
	do 503 n0 = 1,25
	if(n0.le.3)then
	NBIreturn(n0) = NBIlost(44)
	endif
	if(n0.le.7.and.n0.gt.3)then
	NBIreturn(n0) = NBIlost(45)
	endif
	if(n0.le.10.and.n0.gt.7)then
	NBIreturn(n0) = NBIlost(46)
	endif
	if(n0.le.14.and.n0.gt.10)then
	NBIreturn(n0) = NBIlost(47)
	endif
	if(n0.le.17.and.n0.gt.14)then
	NBIreturn(n0) = NBIlost(48)
	endif
	if(n0.le.21.and.n0.gt.17)then
	NBIreturn(n0) = NBIlost(49)
	endif
	if(n0.gt.21)then
	NBIreturn(n0) = NBIlost(50)
	endif
503	continue  

c	 --------------------------------------------------------------------


1043	format (I3, 3f10.3,e10.3,3f10.3,e10.3)
1044	format (I3, 2e10.3,2x,I3,2x,5e10.3)
1046  format(I2,2x,f10.2,2x,f10.2,2x,3f10.2,2x,f10.2,2x,2e12.3,2x,
     12e12.3)
	write(121,'(1x,20A)') '  n0    forblav   morblav  eorblav  Vparav 
     1    forbl2   xmorbl2   eorbl2   Vpar2   ' 						 
	xm = 3.343e-27
	do 2223 n = 1,24
	xmv = 1.13*sqrt(2.*xk*xti(n)/xm)*yy1(n)
	xmv2 = 1.13*sqrt(2.*xk*xti(n)/xm)*xmorbl2(n)
	write(121,1043) n, forbl(n), xmorbl(n), eorbl(n), xmv, forbl2(n),
	1 xmorbl2(n),eorbl2(n),xmv2
	xmomiol(n) = xmv
2223	continue
	write (121,'(1x20A)')' n0   intrinD   intrinC'
	do 2224 n = 1, 24
	yy1(n)=  1.13*sqrt(2.*xk*xti(n)/xm)*yy1(n)
 	yy2(n)=  1.13*sqrt(2.*xk*xti(n)/(6.*xm))*yy2(n)
	write(121,1044) n, yy1(n), yy2(n)
2224	continue 
c	this ends loop on n0 for calculating forbl(n0), morbl(n0),eorbl(n0)
      write (121,'(1x20A)')'n0        lost1       tally        ffrac1
	1 ffrac2    ffrac3    nbdep1    nbdeplost1   nbdepkept1  NBIreturn'
	do 2225 n = 1,50
	write(121,1046) n,lost1(n),tally(n),ffrac1(n),ffrac2(n),ffrac3(n),
	1               nbdep1(n),nbdeplost1(n),nbdepkept1(n),retcur(n)
2225	continue

c     	call IOL(erextot,fpsi0,nbdep,floss,vloss,eloss,lo)     

c     --------------------------------------------------------------------

C	SOL/DIVERTOR PLASMA + PARTICLE/HEAT FLUX + CORE	SOLUTION 
25  	CALL DIVERT3

C	SOLUTION CONVERGED		
C	IF(ITSTOP.EQ.1) GOTO 350

 
C	DIVERTOR/SOL PLASMA-NEUTRAL & CORE PART/HEAT FLUX INTO SOL
C     CONVERGED.

C		EVALUATE THE MARFE DENSITY LIMIT																		
C		EVALUATE RADIAL INSTABILITY GROWTH RATES IN TRANSPORT BARRIER
50	CALL DENLIM 

C	ADJUST CONFINEMENT FOR MARFE (H89=HN=1,C89=C89N=0 IF MI>1)
C	IF(JJOPTPED.NE.10)  GOTO 75
C	XMARAVCAL = 1./(((1./DENLIMXPT(14)) + (1./DENLIMMP(14)))/2.)
C     XMARFCAL = XNBAR/XMARAVCAL 
C     	IF(XMARFCAL.LT.1.005) GOTO 75
C	IMARFE = IMARFE + 1
C	IF(IMARFE.GT.1) GOTO 75
C 	H89 = 1.0 
C	HCONF = 1.0
C	HN = HRAT
C	C89 = 0.0
C	C89N = 0.0
C	GOTO 25


C	CALL LH	
C		EVALUATE THE DISRUPTION DENSITY LIMIT

75	IF(JIT.EQ.1) GOTO 100



      CALL DISRUPT
C		EVALUATE EDGE TRANSPORT COEFFICIENTS
c	CALL TRANSCOEF 
C		EVALUATE IMPURITY CONTAINMENT IN DIVERTOR
c	CALL IMPURITY

100	CALL EDITDIV

C		EVALUATE THERMAL STABILITY OF DIVERTOR
	CALL DIVSTAB

 	IF(IOPTEDGE.EQ.1) CALL EDGECALC 


C*************CALCULATE IOL FLUXES INTO SOL*************************
c	**********integrate the loss to SOL octants********************

	xgamn =   gamion(n0,1)*(1.-forbl(n0))
	xgamv = 1.0
	xgame = gamheati(n0)*(1.-eorbl(n0))
	xgame = 1.0
	xgamn = 1.0
	do 230 n = 1,8
	Soln(n) = 0.0
	Solv(n) = 0.0
	Sole(n) = 0.0
	do 230 np  = 1,8
	areafactor = rhor(n0)/rhor(25)
	do 229 npsi = 1,22 
	if(LO(npsi,np).eq.n)  then
	  soln(n) = soln(n) +  xgamn*floss(npsi,np)*areafactor
  	  solv(n) = solv(n) +  xgamv*vloss(npsi,np)*areafactor
	  sole(n) = sole(n) +  xgame*eloss(npsi,np)*areafactor
	endif
229	continue
230	continue


c	 ***********Miller Eq & 1D divertor solution*********************************
c 	IF(IOPTEDGE.EQ.1) CALL DIVSOL


	do 233 n = 1,30
	
      gam25con(n) = fpart(n+10)*(1.- Forbl(24))
	gamQ25con(n) = fheat(n+10)*(1.- Eorbl(24))
	if(n.ge.22.and.n.le.25) gamorbl(n) =  7.5*gamion(21,1)*Forbl(21)
	if(n.ge.22.and.n.le.25) gamQorbl(n) =  7.5*gamheati(19)*Eorbl(19) 	
	if(n.ge.22.and.n.le.25) gamVorbl(n) =  7.5*gamion(21,1)*XMorbl(21)
     2	*xmas(1)

      if(n.ge.6.and.n.le.9) 
     1gamorbl(n) =  7.5*(fpart(n+10)*Forbl(24) -
     2	gamion(21,1)*Forbl(21))
	if(n.ge.6.and.n.le.9) 
     1gamQorbl(n) =  7.5*(fheat(n+10)*Eorbl(24) 
     2	- gamheati(19)*Eorbl(19))
	  
	if(n.ge.6.and.n.le.9) 
	1 gamVorbl(n) =  7.5*(fpart(n+10)*XMorbl(24) -
	2	gamion(21,1)*XMorbl(21))*xmas(1)
    	gamsum(n) = gam25con(n) + gamorbl(n)
	   gamQsum(n) = gamQ25con(n) + gamQorbl(n)
	  gamiller(n) = fpart(n+10)
	  gamQmiller(n) = fheat(n+10) 
233   continue
c***********distribution of fluxes into SOL**************
	OPEN(unit=132,FILE='distrsol.TXT',STATUS='UNKNOWN') 
c	********ion orbit loss fluxes***********************
	write (132,1061) 
	write (132,1007) (soln(n),n=1,8)
	write (132,1008) (solv(n),n=1,8)
	write (132,1009) (sole(n),n=1,8)
	x=areafactor
	write (132,'(1x,35A)')'mesh  gamcon   gamorb  gamsum    gamiller '
	do 234 n = 1, 30 
	write (132,10110) n,gam25con(n),gamorbl(n),gamsum(n),gamiller(n)
234	continue 
	write (132,'(1x,35A)')'mesh  gamQcon  gamQorb  gamQsum gamQmiller'
 	do 235 n = 1, 30 
	write (132,10110) n,gamQ25con(n),gamQorbl(n),gamQsum(n),
     1	gamQmiller(n)   
235	continue  
	write (132,'(1x,35A)')'	mesh  par_mom '
	do 236 n=1,30
	write (132,10110) n, gamVorbl(n)
236   continue

c     ----------Formatting Options----------------
1000	format(1x,"theta0 =",f10.3) 
1001  format(1x,"theta0=",9f9.3)
1002  format(1x,10e9.3) 
1003  format(6x, "floss",6x, "eloss")
1004	format(1x,2f10.3) 
1005	format(1x,"Flos =", f6.3,3x,"Vlos =",f6.3,3x,"Elos =",f6.3) 
1006	format(1x,"Flos2 =", f6.3,3x,"Vlos2 =",f6.3,3x,"Elos2 =",f6.3) 
1007	format(1x,"soln =", 8e9.3)
1008	format(1x,"solv =", 8e9.3)
1009	format(1x,"sole =", 8e9.3)
1010  format (9e9.3)
10110 format (i4,6e9.3) 
1061  format(1x, "ion orbit loss particle, momentum and energy")
1015	format(1x,"psi0 =", f6.3,3x,"Wtrap =",f9.2,3x,"Tion =",f8.3,
	1		3x,"r0 =",f6.3,3x,"rho =",f6.3)
1016	format(1x,"emin=",8f10.2)
1020	format(1x,"flos =", 8f6.3,1x,"av=",f6.3)
1021	format(1x,"vlos =", 8f6.3,1x,"av=",f6.3)


1022	format(1x,"elos =", 8f6.3,1x,"av=",f6.3)
1023	format(1x,"option2 flos=",f6.3,1x,"vlos=",f6.3,"elos=",f6.3)
1030	format(1x,"floss=", 8f6.3,1x,"sum=",f6.3)
1031	format(1x,"eloss=", 8f6.3,1x,"sum=",f6.3)
1040 	format(1x,"LO(npsi=1,22, m= 1,8)",23I3)
1041	format(f4.2,8i4)
1042	format (i4) 
	
	CLOSE(121,STATUS='UNKNOWN')
	
	CLOSE(132,STATUS='UNKNOWN')

C***************END ION ORBIT LOSS CALCULATION************************ 
1051	CONTINUE

C		SOLVE EDGE NEUTRAL & ION DISTRIBUTIONS

	yy=delx
	yz=delxreal 	 
	GOTO 400
	a1=f(1)

C		TERMINATION EDIT
300	WRITE (6,'(1x,50A)')'TEMPERATURE BELOW 1 EV OR COOLFRAC GT UNITY.
     2  TERMINATED'

305	RECYCLE =  CURSEP*(1.-ALPHASEP)/FLUXPART
	FUELPL = FUELPLIN + FUELPLOUT 
	WRITE(6,118) FUELPF,FUELPL,FUELMP,SPELLET
 	WRITE (6,105) NZIMP1,FZ1,NZIMP2,FZ2,ZEFFC
 	WRITE (6,113) IZINTRIN,FZINTRIN,IZINJECT,FZINJECT,ZEFF
118	FORMAT(1X,'PFFUEL=',E8.3,1X,'PLFUEL=',E8.3,1X,'SOLFUEL=',E8.3,1X,
     2       'PELFUEL=',E8.3)
105	FORMAT (1X,'CORE IMP#1 Z & CONC=',I3.0,F6.4,2X,'IMP#2 Z & CONC=',
     2          I3.0,F6.4,2X,'CORE ZEFF=',F6.4)
113   FORMAT (1X,' DIV IMP#1 Z & CONC=',I3.0,F6.4,2X,'IMP#2 Z & CONC='
     2       ,I3.0,F6.4,2X,'DIV  ZEFF=',F6.4)
120	FORMAT (1X,'POWFRACS: CORERAD='F5.3,1X,'DIVRAD=',F5.3,1X,'ATWALL='
     2,F5.3,1X,'DPLAS=',F5.3,1X,'TOTAL=',F5.3) 

	WRITE(6,120) FRACRAD,FRACRADIV,FRACATOM,FRACPLASDP,FRACTOT
117	FORMAT(1X,'HEATFLUX=',E9.3,1X,'PARTFLUX=',E9.3,1X,'FRACHEATSOL=',E
     29.3,1X,'RECYCLE=',F8.3)
	WRITE(6,117) FLUXHEAT,FLUXPART,FRACSOL,RECYCLE
103   FORMAT (1X,'CENTRAL NE=',E10.3,1X,'EDGE NE=',E10.3,1X,'AVGNE=',
     2	 	E10.3,1X,'ALPHAN=',F4.2)
104	FORMAT (1X,'CENTRAL TE=',E10.3,1X,'EDGE TE=',E10.3,1X,'AVGTE=',
     2		 E10.3,1X,'ALPHAT=',F4.2)
	WRITE (6,103) XN0,XNPED,XNAV,ALPHAN
	WRITE (6,104) T0,TPED,TAV,ALPHAT
	WRITE(6,'(1X,14A)') 'PLASMA DENSITY'  
 	WRITE(6,1011) XND,XNDIV,XNSEP,XNBAR,XNPED
	WRITE(6,'(1X,20A)') 'PLASMA TEMPERATURE'
	WRITE(6,1011) TD, TDIV,TSEP,TBAR,TPED 
	WRITE(6,'(1X,15A)') 'NEUTRAL DENSITY'
	WRITE(6,1011) XNOD,XNODIV,XNOSOL,XNOBAR,XNOPED
	WRITE (6,121) RPARTT,RPARTN,RPARTP,RPARTTC,RPARTC,itern
	WRITE (6,123) FLUXNEUTIN,FLUXIONIN,FLUXPART
121	FORMAT (1X,'PARTICLE CONVERGE: TOTS0LDIV=',E10.5,1X,'ATOMSOLDIV=',
     2E10.5,1X,'IONSOLDIV=',E10.5,1X,'TOTCHAM=',E10.5,1X,'IONCORE=',E10.
     35,'#ITER=',I4) 
123	FORMAT (1X,'FLUXNEUTIN=',E9.4,1X,'FLUXIONIN=',E9.4,1X,
	2        'FLUXPART=',E9.4)

1011	FORMAT(7E10.3)
      GOTO 400

350	WRITE (6,'(1X,27A)') 'ITERATION TERMINATED AT 100'  
	GOTO 100
375	WRITE (6,'(1X,32A)') 'ITERATION STOPPED; X > 1 5 TIMES' 		
400	STOP
	END

      