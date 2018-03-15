      subroutine IOL(erextot,fpsi0,nbdep,floss,vloss,eloss,lo)
	include 'soldiv.fi' 
	
	parameter (nit=9,npsi0=22,nmesh=25,nsol=51,nn=30)
 	dimension v0m(nit,nit),v0p(nit,nit),xke0m(nit,nit),xke0p(nit,nit),
	1		  angle0(nit),angleD(nit),epminhat(npsi0),vlosshat(npsi0),
     2		  eloss(npsi0,nit),epmin(npsi0,nit),emin(npsi0,nit),
     3	      cos0(npsi0),keDm(nit,nit),xkeDM(nit,nit),Hmiller(nsol),
     4		  xkeDp(nit,nit),flos(nit),elos(nit),vlos(nit),
	5	      numloss(npsi0),eminhat(npsi0),flosshat(npsi0),
	6		  elosshat(npsi0),floss2(npsi0),vloss2(npsi0),
     7          eloss2(npsi0),xmomiol(nmesh),forbl2(nmesh),	
     8	      xmorbl2(nmesh),eorbl2(nmesh),lo(npsi0,nit),
	9          floss(npsi0,nit),vloss(npsi0,nit),ffrac(50)
      real nbdep(51),nbsr, nbsrold,nbdepkept(50), nbdeplost,
	1     erextot(51),fpsi0(51),lost(51),tally(51),ffrac,erex0(51),
	2     ephi0(51),forbl,xmorbl,eorbl,floss,vloss,eloss,ephi
	integer rholength, psinum,nmass,npsi,n0

c	 constants
 	RLOSSIOL = 1.0
 	ioptionorb = 1
	ioptxtran = 0
 	xm = 3.343e-27
	eq = 1.6e-19
	nbi = 74000
	nbsr = 6.54e20
      Icur = plasmacur
	Btor = Bphi
	radwall=   0.90
    	radminor = aminor*sqrt(0.5*(1+elong**2))
	dr = 0.1
	RM = rmajor
	xmu0 =1.257					 
	zmass = 3.34e-27
	delnaf = 0.02
	frholength = 50
55	OPEN(unit=121,FILE='orbloss.TXT',STATUS='UNKNOWN')
      open(unit=550,file='bugging.txt',status='old')
      tallytot = 0.0
	losttot = 0.0

c	do loop over launch surface radii,
c      nmass = 1 (fast ions), nmass = 2 (deuterium), nmass = 3 (carbon)
	do 2222 nmass = 1,3

	rholength = 24
	if(nmass.eq.1)then
	rholength = frholength
	endif

	DO 2221 n0 = 1, rholength

	delna0 = delna
	do 3 m = 1,rholength
	erex0(m) = erex(m)
3	continue
 
	if(nmass.eq.1)then
	   delna0 = delnaf*aminor
	   do 4 m = 1,rholength
	      erex0(m) = erextot(m)
4	   continue
          lost(n0) = 0.0
	   tally(n0) = 0.0
      endif 
 
	xdr =  (rholength+1-n0)
	rminor0	= radminor - xdr*delna0

c	set loss surface radius	
	rminorD = radminor
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
 10	continue

	if(nmass.eq.2)then
	do 503 n= 1,24
	ephi(n) = ephi0(n)
503	continue
	endif


c	4-26-11 formulation

	xxI = xmu0*Icur/(6.1832*(radminor**2))	
	fphiD = 1./sqrt(1. + (xxI**2)*(rminorD**2)/(Btor**2))
 	fphi0 = 1./sqrt(1. + (xxI**2)*(rminor0**2)/(Btor**2))
 	fphix = 1./sqrt(1. + (xxI**2)*(rxloss**2)/(Btor**2))

	Tion = xti(n0)
c	 dummy variable tion for fast ion scenario
	if(nmass.eq.1)then
	Tion = 1
	endif

c	calc loss for different direction cosines thru 200, for launch surfance n0
c     only do one calculation for fast ions - all oriented at the same psi0
      if(nmass.eq.1)then
	psi0 = fpsi0(n0)
	psinum = 1 
	endif

c      22 different psi0 locations for deuterium and carbon	
	if(nmass.gt.1)then
	psi0 = 0.01
	psinum = 22
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
	       lost(n0) = 1.0+lost(n0)
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
      ffrac(n0) = lost(n0)/tally(n0)
	losttot = losttot + lost(n0)
	tallytot = tallytot + tally(n0)
  	endif

	write(*,*) nmass,n0

2221 	continue

2222	CONTINUE

	ff = losttot/tallytot
c	Reduce NBI source rate by fast ion loss [ions/s]
	do 501 n0=1,50			   
	nbdepkept(n0) = 0.624E25*pbeam/nbi*(1-ffrac(n0))*nbdep(n0)
	nbdeplost(n0) = 0.624E25*pbeam/nbi*ffrac(n0)*nbdep(n0)
501	continue		
c     Subtract all fast ions that are lost from total source rate to
c     core before attenuation (calculated in fuel.for)

	do 502 n0=1,50
	nbsr = nbsr - nbdeplost(n0)*0.02
502	continue


c	 --------------------------------------------------------------------


1043	format (I3, 3f10.4,e10.3,3f10.3,e10.3)
1044	format (I3, 2e10.3,2x,I3,2x,5e10.3)
1045  format(I4,2f10.3)
1046  format(I3,2x,4f6.3,2x,2e12.3)
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
	write (121,'(1x20A)')' n0   intrinD   intrinC '
	do 2224 n = 1, 24
	yy1(n)=  1.13*sqrt(2.*xk*xti(n)/xm)*yy1(n)
 	yy2(n)=  1.13*sqrt(2.*xk*xti(n)/(6.*xm))*yy2(n)
	write(121,1044) n, yy1(n), yy2(n)
2224	continue 
c	this ends loop on n0 for calculating forbl(n0), morbl(n0),eorbl(n0)
      write (121,'(1x20A)')' n0    lost   tally    ffrac   nbdep
	1nbdeplost     nbdepkept'
	do 2225 n = 1,50
	write(121,1046) n,lost(n),tally(n),ffrac(n),nbdep(n),
	1                nbdeplost(n),nbdepkept(n)
2225	continue
	write(121,'(1x,20A)') 'losttot   tallytot   ff'
      write(121,1045)	losttot, tallytot, ff


	return
	end
