	PROGRAM SOLDIV
C		SOL-DIVERTOR MODEL MAIN PROGRAM
		USE MSIMSL 
	INCLUDE 'SOLDIV.FI' 
	
	
	

	parameter (nit=9,npsi0=22,nmesh=25)
	dimension v0m(nit,nit),v0p(nit,nit),xke0m(nit,nit),xke0p(nit,nit),
	1			angle0(nit),angleD(nit),floss(npsi0,nit),
     3			eloss(npsi0,nit),epmin(npsi0,nit),emin(npsi0,nit),
     2			cos0(npsi0),keDm(nit,nit),xkeDM(nit,nit),
     3			xkeDp(nit,nit),flos(nit),elos(nit)

	CALL SOLDATA

	 
      xxx=epot(24)

	fuelplout = fuelplin
	thetain = theta
	xy = thetain
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
	x = ap
	TESOL = TSOL
	TIPED = TPED
	TEPED = TPED
	XNBAR = XNPED
	TIBAR = TPED
	TEBAR = TPED
	hgy = xnsol
	xxx = hconf
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
C		SOL/DIVERTOR PLASMA + PARTICLE/HEAT FLUX + CORE	SOLUTION 
25	CALL DIVERT3
C		SOLUTION CONVERGED		
C	IF(ITSTOP.EQ.1) GOTO 350
C	IF(ITERTD.GE.4) GOTO 300 
40	CONTINUE
 
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

c		GOTO 1051
C**********ION ORBIT LOSS CALCULATION*********************************
 	ioptionorb = 1
	ioptxtran = 0
 	xm = 3.343e-27
	eq = 1.6e-19
      Icur = plasmacur
	Btor = Bphi
	radwall=   0.88

	radminor = aminor*sqrt(0.5*(1+elong**2))
	
	dr = 0.1
	RM = rmajor
	xmu0 =1.257					 
	zmass = 3.34e-27

	rminor0	= 0.725
	n0 = 1
	rminorD = radminor
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
	
	

	do 99 n = 1,9 	
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
99	continue
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

c	Evaluate Ion and Energy Orbit Loss, for each pitch angle & theta0          	
160   do 161 n = 1,8 
	A = 1.5
	X = EPMIN(npsi,n)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	floss(npsi,n) = (value1 - value)/value1

	A = 2.5
	EVALUE = GAMI(A,X)
	EVALUE1= GAMMA(A)
	eloss(npsi,n) = (evalue1 - evalue)/evalue1
161	continue

	avfloss = 0.
	aveloss = 0
	do 165 n = 1,8
	avfloss = avfloss + floss(npsi,n)/8.
	aveloss = aveloss + eloss(npsi,n)/8.
165	continue

	write(121,1015) psi0,Wtrap,Tion,rminor0,rminor0/radminor 
	write(121,1020) (floss(npsi,n)/2.,n=1,8),avfloss/2.
 	write(121,1021) (eloss(npsi,n)/2.,n=1,8),aveloss/2. 

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
	Flos(n) = 0.0
	Elos(n) = 0.0
	Flos(n) = 0.05*(floss(1,n)+floss(12,n))
	Elos(n) = 0.05*(eloss(1,n)+eloss(12,n))
	Do 220 npsi= 2,10
	Flos(n) = Flos(n) + 0.1*(floss(npsi,n)+floss(npsi+11,n))
	Elos(n) = Elos(n) + 0.1*(eloss(npsi,n)+eloss(npsi+11,n))
220	continue
	Flos(n) = Flos(n) + 0.05*(floss(11,n)+floss(22,n))
	Elos(n) = Elos(n) + 0.05*(eloss(11,n)+eloss(22,n))

	Flos(n) = 0.5*Flos(n)
	Elos(n) = 0.5*Elos(n)
	
225	continue
	write(121,'(1x,35A)') '        0.0   1/4   2/4   3/4   4/4   5/4         
     1 6/4   7/4   ' 
	write(121,1020) (flos(n),n=1,8)
	write(121,1021) (elos(n),n=1,8) 
	flosav = 0.0
	elosav = 0.0
	do 226 n = 1,8
	flosav = flosav + flos(n)/8.
	elosav = elosav + elos(n)/8.
226	continue
	write(121,1005) Flosav,Elosav
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
1003  format(6x, "floss",6x, "eloss")
1004	format(1x,2f10.3) 
1005	format(1x,"Flos =", f6.3,3x,"Elos =",f6.3) 
1010  format (9e9.3) 
1015	format(1x,"psi0 =", f6.3,3x,"Wtrap =",f8.3,3x,"Tion =",f8.3,
	1		3x,"r0 =",f6.3,3x,"rho =",f6.3)
1020	format(1x,"flos =", 8f6.3,1x,"av=",f6.3)
1021	format(1x,"elos =", 8f6.3,1x,"av=",f6.3)
1030	format(1x,"floss=", 8f6.3,1x,"sum=",f6.3)
1031	format(1x,"eloss=", 8f6.3,1x,"sum=",f6.3)
 
 	CLOSE(121,STATUS='UNKNOWN')

C***************END ION ORBIT LOSS CALCULATION************************ 
1051	CONTINUE

C		SOLVE EDGE NEUTRAL & ION DISTRIBUTIONS
	IF(IOPTEDGE.EQ.1) CALL EDGECALC 
c	IF(IOPTEDGE.EQ.1) CALL DIVSOL	
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