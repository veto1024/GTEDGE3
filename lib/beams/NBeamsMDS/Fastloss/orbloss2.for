	subroutine orbloss (erex,delna,ephi,xti,erav,rhor,deltathetx,
     2			ioptionorb,ioptxtran,zbar,ATNUM,exne)

c    include 'soldiv.fi'
	
      parameter (nit=9,npsi0=22,nmesh=25)
	dimension v0m(nit,nit),v0p(nit,nit),xke0m(nit,nit),xke0p(nit,nit),
	1			angle0(nit),angleD(nit),tfloss(npsi0,nit), ffloss(npsi0,nit),
     3			teloss(npsi0,nit),feloss(nspi0,nit),epmin(npsi0,nit),
     2			emin(npsi0,nit),cos0(npsi0),erex(nmesh),ephi(nmesh),
     3			xkeDm(nit,nit),xkeDp(nit,nit),xti(nmesh),tflos(nit),telos(nit),
     4			fflos(nit),felos(nit),erav(nit),rhor(nmesh)


	ioptionorb = 1
	ioptxtran = 1
 	xm = 3.343e-27
	eq = 1.6e-19
      Icur = 1.5
	Btor = -1.98
	radwall=   0.88
	radminor = 0.84
	rmajor = 1.75
	dr = 0.1
	RM = 1.75
	xmu0 =1.257
	zmass = 3.34e-27

	rminor0	= 0.725
	n0 = 1
	rminorD = 0.84
	nD = 25
c	xloss input
	thetax = 4.712
	deltathetx = 0.15			

	deltarx = radminor*deltathetx
	rminorx = radminor - deltarx
	nx = 25
	do 5 m=1,24
	if(rhor(25-m).lt.rminorx) goto 5
     	nx=25-m
5	continue 

	xpi = 3.1416
	xx = xmu0*Icur/(6.283*radminor)
	xy = xmu0*Icur/6.283


c	Calculate electrostatic potential from exp Erad
c	Calculate minimum x-loss energy as av Erad between FS & Sep
	ephi(25) = 0.0
	ersum = erex(25)
	do 10 n =1,24
	m = 25-n
	ephi(m) = ephi(m+1) + radminor*delna*0.5*(erex(m)+erex(m+1))
	ersum =	 ersum + erex(m)
	erav(m) = ersum/(n+1)
10	continue 


c	4-26-11 formulation
55	OPEN(121,FILE='orbloss.TXT',STATUS='UNKNOWN')
    open(122,file=orbloss2.txt',status='unknown') 
 	xk = 1.6e-19
	


	xxI = xmu0*Icur/(6.1832*(radminor**2))	
	fphiD = 1./sqrt(1. + (xxI**2)*(rminorD**2)/(Btor**2))
 	fphi0 = 1./sqrt(1. + (xxI**2)*(rminor0**2)/(Btor**2))
	fphix = 1./sqrt(1. + (xxI**2)*(rxloss**2)/(Btor**2))

	Tion = xti(n0)
	
	psi0 = 0.01
	
	do 200 npsi = 1,22
	if(npsi.eq.2) psi0 = 0.1
	if(npsi.gt.2.and.npsi.lt.11) psi0 = psi0 + 0.1
	if(npsi.eq.11) psi0 = 0.99

	if(npsi.eq.12) psi0 = -0.01
	if(npsi.eq.13) psi0 = -0.1
	if(npsi.gt.13.and.npsi.lt.22) psi0 = psi0 - 0.1
	if(npsi.eq.22) psi0 = -0.99

	Wtrap = rmajor*erex(nD)/(1.0+psi0**2) 
 
	theta0 = 0.0
	thetaD = 0.0
	
	

	do 100 n = 1,9 	
	do  90 m = 1,9 	
	h0 = 1. + (rminor0/RM)*cos(theta0)
 	hd = 1. + (rminorD/RM)*cos(thetaD)
	
    	delphi = ephi(n0) - ephi(nD)
	xpot = 2.*(eq/xm)*delphi
      xflux= xxI*(eq/(1.*xm*hd*fphiD))*(rminor0**2-rminorD**2)

	a =(h0/hd)*(1.-psi0**2)-(((h0/hd)*(fphi0/fphid)*psi0)**2) - 1.
	b = 2.*xflux*(h0*fphi0)*psi0
	c = (xflux**2) - xpot
	abc = 4.*a*c/(b**2)
	if(abc.gt.1.0) goto 80
	v0p(n,m) = -1.*(b/(2.*a))*(1.0+sqrt(1.-4.*a*c/(b**2)))
	v0m(n,m) = -1.*(b/(2.*a))*(1.0-sqrt(1.-4.*a*c/(b**2)))
	goto 85
80	continue	
    	v0p(n,m) = 0.0
	v0m(n,m) = 0.0

85	xke0m(n,m) = 0.5*xm*(v0m(n,m)**2)/xk
	xke0p(n,m) = 0.5*xm*(v0p(n,m)**2)/xk
	xkeDm(n,m) = xke0m(n,m) + (ephi(n0)-ephi(nD))
	xkeDp(n,m) = xke0p(n,m) + (ephi(n0)-ephi(nD))
	thetaD = thetaD + 3.14159/4.
90	continue
	thetaD = 0.0
	theta0 = theta0 + 3.14159/4.
100	continue
	angle0(1) = 0.0
	angleD(1) = 0.0
	do 150 n=2,9
	angle0(n) = angle0(n-1) + 3.14159/4.
	angled(n) = angled(n-1) + 3.14159/4.
150	continue 

c	Minimum Loss Energy for each Theta0
c	Assumes minimum Loss Energy ThetaD for each Theta0
	if(npsi.ge.12) goto 1560
	
	
	Do 155 n=1,8
	emin(npsi,n) = 1.0e9
c	if(v0m(n,m).gt.0.0)	emin(npsi,n) = xke0m(n,1)
	do 154 m=1,8
	if(v0m(n,m).lt.0.0) goto 154 
	if(emin(npsi,n).gt.xke0m(n,m).and.v0m(n,m).gt.0.0) 
	1	emin(npsi,n) = xke0m(n,m)
154	continue
155	continue
	Do 156 n=1,8
	epmin(npsi,n) = emin(npsi,n)/(Tion)
156	continue
	goto 160
c	ctr-current ions
1560	Do 158 n=1,8
	emin(npsi,n) = 1.0e9
c 	if(v0p(n,m).lt.0.0) emin(npsi,n) = xke0p(n,1)
	do 157 m=1,8
	if(v0m(n,m).gt.0.0) goto 157 
	if(emin(npsi,n).gt.xke0p(n,m).and.v0p(n,m).gt.0.0) 
	1	emin(npsi,n) = xke0p(n,m)
157	continue
158	continue
	Do 159 n=1,8
	epmin(npsi,n) = emin(npsi,n)/(Tion)
159	continue

c	Evaluate Fast and Thermal Ion and Energy Orbit Loss, for each pitch angle & theta0          	
160   do 161 n = 1,8 
	A = 1.5
	X = EPMIN(npsi,n)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	tfloss(npsi,n) = (value1 - value)/value1
	
	fvalue = (2/3)*(enbi**(3/2)-epmin(npsi,n)**(3/2))
	fvalue1 = (2/3)*enbi**^3/2)
	ffloss(npsi,n) =  (fvalue-value1+value)/(fvalue1-value1)
	
		

	A = 2.5
	EVALUE = GAMI(A,X)
	EVALUE1= GAMMA(A)
	teloss(npsi,n) = (evalue1 - evalue)/evalue1
	
	fevalue = (2/5)*(enbi**(5/2)-epmin(npsi,n)**(5/2))
	fevalue1 = (2/5)*enbi**^5/2)
	feloss(npsi,n) = (fevalue-evalue1+evalue)/(fevalue1-evalue1)
	
161	continue

	avtfloss = 0.
	avteloss = 0
	avffloss = 0.
	avfeloss = 0
	
	do 165 n = 1,8
	avtfloss = avtfloss + tfloss(npsi,n)/8.
	avteloss = avteloss + teloss(npsi,n)/8.
	avffloss = avffloss + ffloss(npsi,n)/8.
	avfeloss = avfeloss + feloss(npsi,n)/8.
165	continue

	write(121,1015) psi0,Wtrap,Tion,rminor0,rminor0/radminor 
	write(121,1020) (tfloss(npsi,n)/2.,n=1,8),avtfloss/2.
 	write(121,1021) (teloss(npsi,n)/2.,n=1,8),avteloss/2. 

	write(121,'(1x,20A)') 'launch energy, minus sign' 
	write (121,1001) (angle0(n),n=1,9) 
	do 175, m=1,9
 	write (121,1002) angled(m),(xke0m(n,m),n=1,9)
175	continue
	write(121,'(1x,20A)') 'launch energy, plus sign'
	write (121,1001) (angle0(n),n=1,9) 
	do 185, m=1,9
 	write (121,1002) angled(m),(xke0p(n,m),n=1,9)
185	continue	
c	write(121,'(1x,20A)') 'escape energy, minus sign'
c	write (121,1001) (angle0(n),n=1,9) 
c	do 195, m=1,9
c	write (121,1002) angled(m),(xkeDm(n,m),n=1,9)
c195	continue
	write(121,'(1x,20A)') 'velocity, minus sign' 
	write (121,1001) (angle0(n),n=1,9) 
	do 190, m=1,9
	write (121,1002) angled(m),(v0m(n,m),n=1,9)
190	continue
	write(121,'(1x,20A)') 'velocity, plus sign'
	write (121,1001) (angle0(n),n=1,9) 
	do 195, m=1,9
	write (121,1002) angled(m),(v0p(n,m),n=1,9)
195	continue
	
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
	tFlos(n) = 0.0
	tElos(n) = 0.0
	fFlos(n) = 0.0
	fElos(n) = 0.0
	
	tFlos(n) = 0.05*(tfloss(1,n)+tfloss(12,n))
	tElos(n) = 0.05*(teloss(1,n)+teloss(12,n))
	fFlos(n) = 0.05*(ffloss(1,n)+ffloss(12,n))
	fElos(n) = 0.05*(feloss(1,n)+feloss(12,n))
	
	Do 220 npsi= 2,10
	tFlos(n) = tFlos(n) + 0.1*(tfloss(npsi,n)+tfloss(npsi+11,n))
	tElos(n) = tElos(n) + 0.1*(teloss(npsi,n)+teloss(npsi+11,n))
	fFlos(n) = fFlos(n) + 0.1*(ffloss(npsi,n)+ffloss(npsi+11,n))
	fElos(n) = fElos(n) + 0.1*(feloss(npsi,n)+feloss(npsi+11,n))
220	continue
	tFlos(n) = tFlos(n) + 0.05*(tfloss(11,n)+tfloss(22,n))
	tElos(n) = tElos(n) + 0.05*(teloss(11,n)+teloss(22,n))
	fFlos(n) = fFlos(n) + 0.05*(ffloss(11,n)+ffloss(22,n))
	fElos(n) = fElos(n) + 0.05*(feloss(11,n)+feloss(22,n))

	tFlos(n) = 0.5*tFlos(n)
	tElos(n) = 0.5*tElos(n)
	fFlos(n) = 0.5*fFlos(n)
	fElos(n) = 0.5*fElos(n)
	
225	continue
	write(121,'(1x,35A)') '        0.0   1/4   2/4   3/4   4/4   5/4         
     1 6/4   7/4   ' 
	write(121,1020) (tflos(n),n=1,8)
	write(121,1021) (telos(n),n=1,8) 
	tflosav = 0.0
	telosav = 0.0
	do 226 n = 1,8
	tflosav = tflosav + tflos(n)/8.
	telosav = telosav + telos(n)/8.
226	continue
	write(121,1005) tFlosav,tElosav
c	write(121,'(1x,35A)') '       psi0     0.0   1/4   2/4   3/4   4/4            
c     1  5/4   6/4   7/4    av  ' 
c	psi0 = 0.1 
c	do 250 npsi = 1,22
c	if(npsi.eq.2) psi0 = 0.1
c	if(npsi.gt.2.and.npsi.lt.11) psi0 = psi0 + 0.1
c	if(npsi.eq.11) psi0 = 0.99

c	if(npsi.eq.12) psi0 = -0.01
c	if(npsi.eq.13) psi0 = -0.1
c	if(npsi.gt.13.and.npsi.lt.22) psi0 = psi0 - 0.1
c	if(npsi.eq.22) psi0 = -0.99
	
c  	write(121,1030) psio,(floss(npsi,n),n=1,8),avfloss(npsi)
c250	continue
c	write(121,'(1x,35A)') '        psi0    0.0   1/4   2/4   3/4   4/4            
c    1  5/4    6/4   7/4   av  ' 
c	psi0 = 0.1
c	do 255 n = 1,22
c	if(npsi.eq.2) psi0 = 0.1
c	if(npsi.gt.2.and.npsi.lt.11) psi0 = psi0 + 0.1
c	if(npsi.eq.11) psi0 = 0.99

c	if(npsi.eq.12) psi0 = -0.01
c	if(npsi.eq.13) psi0 = -0.1
c	if(npsi.gt.13.and.npsi.lt.22) psi0 = psi0 - 0.1
c	if(npsi.eq.22) psi0 = -0.99

c	write(121,1031) psi0, (eloss(npsi,n),n=1,8),aveloss(npsi)
c255	continue
 
1000	format(1x,"theta0 =",f10.3) 
1001  format(1x,"theta0=",9f9.3)
1002  format(1x,10e9.3) 
1003  format(6x, "tfloss",6x, "teloss")
1004	format(1x,2f10.3) 
1005	format(1x,"tFlos =", f6.3,3x,"tElos =",f6.3) 
1010  format (9e9.3) 
1015	format(1x,"psi0 =", f6.3,3x,"Wtrap =",f8.3,3x,"Tion =",f8.3,
	1		3x,"r0 =",f6.3,3x,"rho =",f6.3)
1020	format(1x,"tflos =", 8f6.3,1x,"av=",f6.3)
1021	format(1x,"telos =", 8f6.3,1x,"av=",f6.3)
1030	format(1x,"tfloss=", 8f6.3,1x,"sum=",f6.3)
1031	format(1x,"teloss=", 8f6.3,1x,"sum=",f6.3)
 
 	CLOSE(121,STATUS='UNKNOWN')







400 	return
	END

