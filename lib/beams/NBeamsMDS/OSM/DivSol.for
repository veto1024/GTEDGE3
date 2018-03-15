	SUBROUTINE DIVSOL
	INCLUDE 'SOLDIV.FI'
c	  Solves 1D Sol-Div density, temp., flow & potential eqs
c	  5/6/08
	parameter (NCAP=51,jq=2,jp=6)
      dimension AA(NCAP),dpsi(NCAP),AALPHA(NCAP),AACOEF(NCAP),
     1		SORC(NCAP),QHEAT(NCAP),SRAD(NCAP),SAT(NCAP),SOL(NCAP),
     2		XT(NCAP),YN(NCAP),DVDPSI(NCAP),gpart(ncap),difp(ncap),
     3		xnupart(ncap),xnumom(ncap),DPRESS(NCAP),VPAR(NCAP),
     4		GAMF(NCAP),AMCOEF(NCAP),APCOEF(NCAP),XNT32(JQ),BIGQ(NCAP),
     5		BIGAM(NCAP),SPART(NCAP),xnov(ncap),xnovcold(ncap),
     6		bigMOM(ncap),vee(ncap),ETAM(NCAP), PRESZ(NCAP),ynm(ncap),
     7		VISCZ(NCAP),CONVZ(NCAP),VTHERM(NCAP),YNOLD(NCAP),
     8		xtold(ncap),ESPOT(NCAP),XJ(NCAP),DELJB(NCAP),DOTRZ(NCAP),
     9		ZLNM(NCAP),ZLTM(NCAP),ZLPM(NCAP),ZLVM(NCAP),vdrz(ncap),
     1		DELSUM(NCAP),VEXB(NCAP),GAMErXB(NCAP),YNZ(NCAP),
     2		BIGAMZ(NCAP),BIGMOMZ(NCAP),YNZM(NCAP),VEZ(NCAP),
     3		YNA(NCAP),BIGAMA(NCAP),VEA(NCAP),SA(NCAP),FZCARBON(NCAP),
     4    	FZARGON(NCAP),GAMDB(NCAP),VEXBZ(NCAP),VEXBA(NCAP),
     5        vep(ncap),GAMErXBZ(NCAP),GAMErXBA(NCAP),vdrdb(ncap),
     6        VDRIFT(NCAP),dothetz(ncap),vdrdbcarb(ncap),vdrzcarb(ncap),
     7	    deljcarb(ncap),gamdbcarb(ncap),vdRzarg(ncap),
     8	    deljarg(ncap),gamdbarg(ncap),vdrdbarg(ncap),XZEFF(NCAP),
	9        ZCARB(NCAP),DELTEMP(NCAP),ER(NCAP),AXCOEF(NCAP),
     1	    EPAR(NCAP),VDRIFTZ(NCAP),VDRIFTA(NCAP),ASOR(NCAP),
     2	    BIGA(NCAP),ALF(NCAP),XJB(NCAP),SUMLN(NCAP),SUMJB(NCAP),
     3	    SUMJ(NCAP),SUMCUR(NCAP),GRADBCUR(NCAP),THET(NCAP),
	4        XJPS(NCAP),xpirthet2(ncap),vparz(ncap),xpirthet2z(ncap),
	5	    SR(NCAP),vdia(ncap),vdiaz(ncap),vdiar(ncap),angle(ncap),
     6	    bigH(ncap),cheatx(ncap),delrthet(ncap),zn(ncap),zt(ncap),
	7	    xpirthet(ncap),hrmet(ncap),hthetmet(ncap), xpirphi(ncap),
	8	    delthet(ncap),convphi(ncap),convthet(ncap),veepar(ncap),
     9	    veephi(ncap),veethet(ncap),xpirthet1(ncap),
	1	    heatnet(ncap),tno(ncap),vezphi(ncap),vezthet(ncap),
     2	    vezpar(ncap),xpirphiz(ncap),xpirthet1z(ncap),
     3	    xpirthetz(ncap),coef1(ncap),coef2(ncap),coef3(ncap),		
     4	    delrm(ncap),delrp(ncap),bigmomB(ncap),bigqB(ncap),
     5        bigamB(ncap)
	
      double precision TEDP,XLZ2DP,DLDT2DP,Tdbl,xlzdbl,dlzdbl
	real iter

c	INITIAL STAGNATION POINTS & FLUXES 
	nstagp = 26
	nstagq = 26
	NSTAG =  26
	Gstag = 0.0
	Qstag = 0.0
c	Thomsom location and measured values
	nThom = 30
	xnsep = xnsepex
	tsep = 0.5*(tsepexe + tsepexi)
	tsep = tsepexe
	bigMom(nstag) = 1.5*4.*xk*tsep*xnsep

C	SET IMPURITY = 0 TO USE INPUT FZ FOR CARB AND ARG IN RADIATION CALCULATION, 1 TO USE CALC FZ
	IMPURITY = 0

C*********************OPTIONS***************OPTIONS*****************************
C	SET IGRADB = 0 TO TURN OFF GRADB & CURV DRIFTS, 1 TO TURN THEM ON
	IGRADB = 1	
C	SET IOPTEXB = 0 TO TURN OFF EXB DRIFTS, 1 TO TURN THEM ON
	ioptexb = 1
C	SET IOPTDIA = 0 TO TURN OFF DIA DRIFTS, 1 TO TURN THEM ON
	IOPTDIA = 1
C	SET ITRANS = 0 TO TURN OFF TRANSPORT LOSS TERMS, 1 TO TURN THEM ON
C	ITRANS = 0.  THIS DRIVES NEGATIVE SOLUTION, BOHM IS TOO MUCH TRANSPORT 
C	SET IB = 0 TO HAVE B IN NORMAL (-) DIRECTION, =1 TO REVERSE B DIRECTION (+)
	IB = 0
c******************************************************************************
c	SET PARALLEL DIFFUSIVITY MULTIPLIER
	REDKAPPA =1.50
c******************************************************************************
c	set heat and particle perp influx multipliers
	cpart =0.35
	cheat = 0.9
	
	 
	fluxheat = cheat*fluxheat
	fluxpart = cpart*fluxpart
c*************************************************************
c	ANOMALOUS PARALLEL MOMENTUM DRAGS
	XNUMOMDIV = 0.0E3
	XNUMOMSOL = 0.0E2
C*************************************************************
c**********Poloidal Asymmetry**********************************
c	if iflux = 1, uses Miller equil dist, iflux = 0 uniform dist--qperp and gamperp
	iflux=1
c	if isol = 1, uses Miller equil dist, isol = 0 uniform dist--deln,delT, delE, delheat
	isol = 1
c	********************************************************
c	symetric (isym=0) or nonsymetric (isym=1) divertor
	isym = 1 
c***************neutral multiplier, used to adjust solution***********
c	no = nocalc/cn
	cn = 9.5
	 
c**********************************************************************

C*****************************DIVERTOR FLUX SURFACE EXPANSION**********
c	epdin is expansion inner divertor width relative to <sol>
c	epxin is expansion inner xpt width relative to <sol>
	EPDin = 3.0
	epxin = 3.0
c	epdout is expansion outer divertor width relative to <sol>
c	epxout is expansion outer xpt width relative to <sol>
	epdout = 3.0
	epxout = 3.0
c******************radiation multiplier***********************************
c	frad multiplier in divertor and sol srad
	fradin = 0.1
	fradsol = 0.1
	fradout = 0.1
c************************************************************************	
c	Rinout is ratio of neutral density in inner to outer divertor
	Rinout= 1.15
C***************************************************************** 
c	TEMPERATURE DENSITY OPTION
C*****	IOPTTHOM = 0 FOR DIRECT SOLUTION FROM ITERATED T(1), N(1)
C*****	 		 = 1 FOR SOLUTION NORMALIZED TO T(30)=TSEP N(30)=NSEP
	ioptthom = 0

	xtout= 10.
	xtin = 10. 


C	VBIASIN IS THE BIAS OF THE INNER DIVERTOR RELATIVE TO THE OUTER DIVERTOR
	VBIASIN = 0.0

	ITRANS = 0
C	ALB IS ANOMALY FACTOR OVER BOHM TRANSPORT
	ALB = 1.0
c	DTRAN IS THE RADIAL TRANSPORT COEFFICIENT
	DTRAN = 0.0
c	DELTIN IS THE TEMP GRAD SCALE LENGTH TO BE USED TO CALCULATE ER & SET SOL WIDTH
	DELTIN = 0.01
	XJIN = 1.E4	
c		xlpar1 and xlpar2 are the parallel distances from top of plasma to divertor
c	    plate on the outboard and inboard, respectively.
c	   	put 10 mesh intervals in each divertor channel & 2 mesh intervals of same size
c		in lower SOL

c	************symmetric*********************
	divdepth2=divdepth1
	xlpar2 = xlpar1 
	xlperp2 = xlperp1
c	*******************************************
c	divdepth is the projection of the divertor leg on the poloidal plane.
c	dellt is the projection of the thickness of the recycling region at the
c	divertor target onto the poloidal plane.  these quantities used in 2-pt model
	sinal = abs(zx-zsep1)/divdepth1
	xl1 = (xlpar1-xlperp1)-4.*dellt/betag
	xl2 = (xlpar2-xlperp2)-4.*dellt/betag

	dpsi(1) = dellt/betag
	dpsi(2) = 3.*dellt/betag
	do 5 n = 3, 10
	dpsi(n) = xl2/7.5
5	continue
	dpsi(ncap-1) = dellt/betag
	dpsi(ncap-2) = 3.*dellt/betag
	do 6 n = 3,10
	m = ncap - n
	dpsi(m) = xl1/7.5
6	continue

c		divide rest of SOL into 31 equal intervals
	zlsol = (xlpar1+xlpar2) - 0.5*(dpsi(10)+dpsi(41))
	do 7 n = 11,40
	dpsi(n) = zlsol/30.
7	continue
	xlperpcheck = 0.5*(dpsi(10)+dpsi(41))
	xlparcheck = 0.
	do 8 n = 1,ncap-1
	xlparcheck = xlparcheck + dpsi(n)
8	continue
	do 9 n = 11, 40
	xlperpcheck = xlperpcheck + dpsi(n)
9	continue	
C		cells with heat & part fluxes from plasma

C	DELT CALC AS DELT=CHIBOHM*Tsep*nsep/heatflux
	CONDHEAT = FLUXHEAT - 3.*XK*FLUXPART*TSEP 
	DELT = (5.*alb*XK*XNSEP*(TSEP**2)/(32.*Bfield*FLUXHEAT))
c	DBOHM = TSEP/16*B 
c	CHIBOHM = (5./32.)*ALB*TSEP/B
C	DELN CALC AS DELN =NSEP*TSEP/16*B*FLUXPART 
	DELN = XNSEP*TSEP/(16.*BFIELD*FLUXPART)
	DELE = 2.*DELT/7.
c********************common gradient scale length************
	delt = deltin
	dele = deltin
	deln = deltin
c***************************************************************
c	goto 1111
c	shafranov shift parameter--shifted circles
	betap = 0.5 
	xinduct = 0.5  
	dR0dr= -1.*(aminor/rmajor)*(betap + 0.5*xinduct)
c	evaluate integrals over theta
	dx = 6.28/30.
	x = -0.5*dx
	xint = 0.0
	xnum = 0.0
	xden = 0.0
	abarm = 0.0	
	d23 = 0.0
	x24 = 0.0
	xmas(1) = 3.34e-27
	
	do 230 n = 1,30
	elong0 = elong
	triang0 = triang
c	if(n.gt.15) elong0 = 2.32
c	if(n.lt.11) then
c		elong0 = elongtop
c		triang0= triangtop
c	endif
c	elong0 = 1.0
c	triang0 = 0.0
c	dR0dr = 0.0
	x = x + dx 
c	integrals in G(r,theta)
	xx = asin(triang0)
	angle(n) = x
	xint = xint + dx
	zz = x + xx*sin(x)
	yx = rmajor + aminor*cos(zz)
	ys = sqrt((cos(zz)**2)+((elong0*sin(x))**2)) 
	yc = sqrt(((sin(zz)**2)*((1+xx*cos(x))**2)+((elong0*cos(x))**2)))
	yr = 1.+dR0dr*cos(x)
	yr0 = 1.+dR0dr
	zy = xx*sin(x)
	yr = cos(zy)+dR0dr*cos(x)+xx*cos(x)*sin(x)*sin(zz)
	yc = sqrt((sin(zz)**2)*((1.+xx*cos(x))**2) + (elong0*cos(x))**2)
	ftheta = yr/yc
	delrthet(n) = (yc/elong0)/yr
	delthet(n) = aminor*ys
	xden = xden + dx*yx*ys
	xnum = xnum + dx*(yr/yc)*yx*ys*elong0
	abarm = abarm + dx*yx*ys
230	continue
c	evaluate H(r,theta)

	sem = (xnum/xden)/yr0
	abarm = (abarm/6.28)*(aminor/rmajor)
	se = sem
	do 232 n = 1,30
	bigH(n) =(xnum/xden)*delrthet(n)
232	continue
	do 235 n = 11, 40
		if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
		angthet(n)= angle(m)	
235	continue	


	delheat = delE
	delpart = deln

c	coef to account for change in radial width with location
C****************radial widths**************************************
c	coef3 is the ratio of local ribbon width to the average ribbon width
	do 240 n = 1,ncap-1
c	isol = 0 for uniform ribbon width
	if(isol.eq.0) then
		coef3(n) = 1.0
		goto 240
	endif
			if(n.ge.11.and.n.le.33) m = 34 - n
			if(n.ge.34.and.n.le.40) m = 64 - n
	
	if(n.lt.10)	coef3(n) = epdin
	if(n.eq.10) coef3(n) =epxin
	if(n.gt.10.and.n.lt.41) coef3(n) = 1./bigH(m)
	if(n.eq.41) coef3(n) = epxout
	if(n.gt.41) coef3(n) = epdout

240	continue
c	coef2 and coef3 are the coefficients for integrating the heat and
c	particle fluxes taking into account flux surface expansion/contraction

c	isol = 0 for uniform ribbon width
c	if(isol.eq.0) then
c		coef2(n) = 1.0
c		coef3(n) = 1.0
c		goto 246
c	endif
c	isol = 1 for Miller del/delav and expansion in divertor

	delrm(1) = epdin
	do 241 n = 2, ncap-1
	dm = dpsi(n)/(dpsi(n)+dpsi(n-1))
	delrm(n) = dm*coef3(n-1) + (1. - dm)*coef3(n)
241	continue
	delrp(50) = epdout
	do 242 n = 1, ncap-2
	dp = dpsi(n)/(dpsi(n)+dpsi(n+1))
	delrp(n) = dp*coef3(n+1) + (1. - dp)*coef3(n)
242	continue 
	do 245 n = 1, ncap-1
	coef1(n) = delrm(n)/delrp(n)
	coef2(n) = coef3(n)/delrp(n)
	goto 245


	
     	if(n.lt.10) then 
		coef1(n) = 1.
		coef2(n) = 1.
		
	endif
	if(n.eq.10) then
		coef1(n) = (epdin+epxin)/((1./bigH(23))+epxin)
		coef2(n) = (2.*epxin)/((1./bigH(23))+epxin)
	    
	endif
	if(n.eq.11) then
		coef1(n) = ((1./bigH(m))+epxin)/((1./bigH(m-1))+(1./bigH(m)))
		coef2(n) =	(2./bigH(m))/((1./bigH(m-1))+(1./bigH(m))) 
		
	endif	
	if(n.ge.12.and.n.lt.33) then
		coef1(n) = ((1./bigH(m))+(1./bigH(m+1)))/
     1     		((1./bigH(m))+(1./bigH(m-1)))
	    coef2(n) = (2./bigH(m))/((1./bigH(m))+(1./bigH(m-1)))
		
	endif
	if(n.eq.33) then
		coef1(n) = ((1./bigH(m))+(1./bigH(2)))/
     1		((1./bigH(30))+(1./bigH(1)))
		coef2(n) = (2./bigH(m))/((1./bigH(30))+(1./bigH(1)))
 		
	endif
	if(n.eq.34) then
		coef1(n) = ((1./bigH(m))+(1./bigH(1)))/
     1			((1./bigH(m))+(1./bigH(m-1)))
		coef2(n) = (2./bigH(m))/((1./bigH(m))+(1./bigH(m-1)))
		
	endif 	
	if(n.ge.35.and.n.lt.40) then
		coef1(n) = ((1./bigH(m))+(1./bigH(m+1)))/
     1			   ((1./bigH(m))+(1./bigH(m-1)))	
     		coef2(n) = (2./bigH(m))/((1./bigH(m))+(1./bigH(m-1)))
		
	endif
	if(n.eq.40) then
		coef1(n) = ((1./bigH(m))+(1./bigH(m+1)))/
     1		(epxout + (1./bigH(m)))
		coef2(n) = (2./bigH(m))/(epxout + (1./bigH(m)))	
		
      endif
	if(n.eq.41) then
		coef1(n) = (epxout +(1./bigH(24)))/(epxout+epdout)
		coef2(n) = (2.*epxout)/(epxout+epdout)
		
      endif
	if(n.eq.42) then
		coef1(n) = (epxout+epdout)/(epdout+epdout)
		coef2(n) = (2.*epdout)/(epdout+epdout)
		
	endif
	if(n.gt.42) then
		coef1(n) = 1.0
		coef2(n) = 1.0
		
	endif  
245	continue
246	continue

c****set Uniform(iflux=0) or Miller Eq (iflux=1) Qperp and Gamperp distributions**************
	do 10 n = 11, 40
	cheatx(n) = 1.0
	if(n.ge.11.and.n.le.33) m = 34 - n
	if(n.ge.34.and.n.le.40) m = 64 - n

	if(iflux.eq.1) then	 		
		cheatx(n) = bigH(m)
		angthet(n)= angle(m)	
	endif
c********flux surface expansion/contraction*************
c	if(isol.eq.1) then 
c		delheat = delheat/bigH(m)
c		delpart = delpart/bigH(m)
c	endif
	C3 = coef3(n)
	qheat(n) = cheatx(n)*fluxheat/(delheat*c3)
	gpart(n) = cheatx(n)*fluxpart/(delheat*c3)

10	continue
	            													
c*********************end variable width setup***************************

 
c		goto 914
1111	continue


c*********************neutral distribution*********************************

c******************************
c	adjusts neutral level	
	xnod=xnod/cn
	xncoldd = xncoldd/cn 	
	xnosolpxpt = xnosolxpt/cn
	xnocoldsolxpt = xnocoldsolxpt/cn
	xnosol=xnosol/cn
	xnocoldsol=xnocoldsol/cn
c***************************************
c	sets neutral distribution in div-sol
	xnov(1) = xnod*rinout
	xnov(ncap-1)= xnod
  	xnovcold(1) = xncoldd*rinout
	xnovcold(ncap-1) =xncoldd

	aldiv = log(xnod/xnosolxpt)/9.5
	cldiv = log(xncoldd/xnocoldsolxpt)/9.5
	
	do 11 n = 2,10
	
	xnov(n) = (xnod*rinout)*exp(-1.*(n-1)*aldiv)	
	xnov(ncap-n) = xnov(n)/rinout

	xnovcold(n) = (xncoldd*rinout)*exp(-1.*(n-1)*cldiv)
	xnovcold(ncap-n) = xnovcold(n)/rinout
11	continue
C	PUMP OPTION FOR LOWER OUTER DIVERTOR
C	REMOVES 5% NEUTRALS IN OUTER RECYCLING REGION
C	IPUMP = 1 PUMP ON; = 0 PUMP OFF
C*************SET BELOW ALSO BEFORE ATCOOL***********************
	IPUMP = 0
	IF(IPUMP.EQ.1) THEN
		XNOV(NCAP-1) = 0.95*XNOV(NCAP-1)
C		XNOV(NCAP-2) = 0.5*XNOV(NCAP-2)
		XNOVCOLD(NCAP-1) = 0.95*XNOVCOLD(NCAP-1)
C		XNOVCOLD(NCAP-2) = 0.5*XNOVCOLD(NCAP-2)
   	ENDIF
c******************************************************************
C	sol distribution
	alsol = log(xnosolxpt/xnosol)/9.
	clsol = log(xnocoldsolxpt/xnocoldsol)/9.
	do 12 n = 11,25
	xnov(n) = xnosolxpt*rinout*exp(-1.*(n-9)*alsol)
	xnov(ncap-n) = xnov(n)/rinout
	xnovcold(n) = xnocoldsolxpt*rinout*exp(-1.*(n-9)*clsol)
	xnovcold(ncap-n) = xnovcold(n)/rinout
12	continue
C		initial conditions
	XNSOL = XNSEP
	ITER = 0
13	xnd1=xnd
	ITER = ITER + 1	
	xnd2=xnd
	td1 = td
	td2 = td
c	TD2 = 13.0
c		atoMic physics frequencies & initial densities
c			recycling region inner leg
	xnupart(1) = xnov(1)*sviond - xnd*svrecd
	xnumom(1) = xnovcold(1)*svatd + xnd*svrecd
	yn(1) = xnd2
	XT(1) = TD2

c			divertor channel inner leg
	do 15 n = 2,9
	yn(n) = xnd - (xnd - xnsep)*(n-1)/9.
	XT(N) = TD - (td - tsep)*(n-1)/9.
	yn(51-n) = yn(n)
	xt(51-n) = xt(n)
	xnupart(n) = xnov(n)*sviondiv - yn(n)*svrecdiv
 	xnumom(n) = xnovcold(n)*svatdiv + yn(n)*svrecdiv
15	continue
c			SOL	
	xnupart(10) = xnov(10)*svionsol
	xnumom(10) = xnovcold(10)*svatsol
	yn(10) = xnsep
	xt(10) = tsep
	do 16 n =  11,40
	xnupart(n) = xnov(n)*svionsol
	xnumom(n) = xnovcold(n)*svatsol
	yn(n) = xnsep
	XT(N) = TSep
16	continue
	xnupart(41) = xnov(41)*svionsol
	xnumom(41) = xnovcold(41)*svatsol
	yn(41) = xnsep
	xt(41) = tsep   
c			divertor channel outer leg
	do 17 n = 42,49
	
	xnupart(n) = xnov(n)*sviondiv - yn(n)*svrecdiv
  	xnumom(n) = xnovcold(n)*svatdiv + yn(n)*svrecdiv
17	continue
c			recycling region outer leg
	xnupart(50) = xnov(50)*sviond - xnd*svrecd
 	xnumom(50) = xnovcold(50)*svatd + xnd*svrecd   
	yn(50) = xnd1
	XT(50) = TD1

	ZEFF=(1.+(IZINJECT**2)*FZINJECT + 4.*FHE +(IZINTRIN**2)*FZINTRIN)/
     2      (1.+ IZINJECT*FZINJECT + 2.*FHE + IZINTRIN*FZINTRIN) 
	do 180 n = 2,ncap-1
	zn(n) = 0.5*(yn(n-1)+yn(n))
	zt(n) = 0.5*(xt(n-1)+xt(n))
180	continue
	zn(1) = yn(1)
	zn(ncap) = yn(ncap-1)
	zt(1) = xt(1)
	zt(ncap) = xt(ncap-1)		
c	the entire subsequent calculation it iterated back to #1000 10 times


 	iconverge = 0 

1000	CONTINUE

c***************************************************************************
c***************************************************************************
C	ITERATE ON HEAT BALANCE AND TEMPERATURE DISTRIBUTION THROUGH 211
	ITERT = 0
	Do 211 j = 1, 50
	
	
	DO 18 N = 1, NCAP-1
	YNOLD(N) = YN(N)
	XTOLD(N) = XT(N)
18	CONTINUE
	ITERT = ITERT + 1
C			COULOMB LOGARITHM,PARALLEL CONDUCTIVITY & VISCOSITY
	XNN = 0.5*(XNSEP+XND)
	XTT = 0.5*(TSEP+TD)	
 	Y = 25.3 - 1.15*LOG10(1.E-6*XNN) + 2.3*LOG10(XTT)
 	IF(XTT.LT.50.) Y = 23.4-1.15*LOG10(1.E-6*XNN)+3.45*LOG10(XTT)
	
	 
	XKAPPA = (1.e4/(ZEFF*Y))*REDKAPPA
	
	xmas1 = 3.343e-27
 	ETAPAR = XMAS1*XKAPPA														   
	ETAPAR =4.8e-6/(ZEFF*Y)

C		EQ 13.13 FUSION PLASMA PHYSICS
	XK = 1.6E-19 
	X = SECEL
	XPI = 3.14159
	
	XMASEL = 9.108E-31
C		assuming TI=TE at plate
	GAMSHEATHIN = 2. + 2./(1.-X) +
     1	 0.5*LOG((((1-X)**2)*(XMAS(1)/XMASEL))/(4.*XPI))
c		subtract convective flux
C	gamsheathin = gamsheathin - 5.       
	GAMSHEATHOUT = GAMSHEATHIN 


C		SETUP

	SUMHEAT1 = 0.0
	SUMRAD1 = 0.0
	SUMAT1 = 0.0
	SUMHEAT2 = 0.0
 	SUMRAD2 = 0.0
	SUMAT2 = 0.0
	SUMSHEATH = 0.0
	SUMPART1 = 0.0
	SUMPART2 = 0.0 
	sumheatsep = 0.0
	chsep = 0.0
	sumpartsep = 0.0
	
	IF(IMPURITY.GT.0) THEN
	FZINTRIN = YNZ(1)/YN(1)
	FZINJECT = YNA(1)/YN(1)
	ENDIF
	ZARG = 10.
	xxx= XNOv(1)
	xxy=XNovCOLD(1)
	tno(1) = todd
	CALL ATCOOL2(XT(1),YN(1),SATZ)
	zcarb(1) = 4.0
	if(xt(1).gt.100.) zcarb(1) = 6.0
	XZEFF(1) = (YN(1)+(Zcarb(1)**2)*YNZ(1)+(ZARG**2)*YNA(1))/
	1	(YN(1)+Zcarb(1)*YNZ(1)+ZARG*YNA(1))
 	
	SAT(1) = SATZ
	SRAD(1) = SRADZ*fradin
	SPART(1) = XNUPARTZ*YN(1) + GPART(1)
	xnumom(1) = xnumomz	+ xnumomdiv
	xnupart(1)= xnupartz		
	tnodiv = 0.5*(tpl+tpf)
	IF(IMPURITY.GT.0) THEN
 	FZINTRIN = YNZ(2)/YN(2)
	FZINJECT = YNA(2)/YN(2)
	ENDIF 
	zcarb(2) = 4.0
	if(xt(2).gt.100.) zcarb(2) = 6.0
 	xxx= xnov(2)
	xxy=xnovcold(2)
	tno(2) = tnodiv
	CALL ATCOOL2(XT(2),YN(2),SATZ)
	XZEFF(2) = (YN(2)+(Zcarb(2)**2)*YNZ(2)+(ZARG**2)*YNA(2))/
	1	(YN(2)+Zcarb(2)*YNZ(2)+ZARG*YNA(2))
 
	SAT(2) = SATZ
	SRAD(2) = SRADZ*fradin
	SPART(2) = XNUPARTZ*YN(2) + GPART(2)
	xnumom(2) = xnumomz + xnumomdiv
	xnupart(2)= xnupartz
      DO 25 N =3,NCAP-2
	
	
		xxx = xnov(n)
		xxy = xnovcold(n)
	
	tnoz = 0.5*(tpl+tpf)
	
C	REGIONS 10-41 ARE SOL
	
	tno(n) = tnoz

C	REGIONS 10 & 41 HAVE X-PT NEUTRAL DENSITIES
c	IF(N.EQ.10.OR.N.EQ.41) THEN
c		XNOZ = XNOSOLXPT
c		XNOCOLDZ = XNOCOLDSOLXPT
c		tnoz = tpf
c	ENDIF
	IF(IMPURITY.GT.0) THEN
 	FZINTRIN = YNZ(N)/YN(N)
	FZINJECT = YNA(N)/YN(N)
	ENDIF 
	zcarb(n) = 4.0
	if(xt(n).gt.100.) zcarb(n) = 6.0
	
		xxx=xnov(n)
		xxy=xnovcold(n)
	tno(n) = tnoz
	
	if(n.ne.49) goto 1900
	xrtyu=0.9
1900	CALL ATCOOL2(XT(N),YN(N),SATZ)
	XZEFF(N) = (YN(N)+(Zcarb(n)**2)*YNZ(N)+(ZARG**2)*YNA(N))/
	1	(YN(N)+Zcarb(n)*YNZ(N)+ZARG*YNA(N))
 	
	SAT(N) = SATZ
	if(n.le.10) fradiv = fradin
	if(n.gt.10.and.n.lt.41) fradiv = fradsol
	if(n.ge.41) fradiv = fradout
	SRAD(N) = SRADZ*fradiv
	SPART(N) = XNUPARTZ*YN(N) + GPART(N)
	xnumom(n) = xnumomz
	if(n.lt.10.or.n.gt.41) xnumom(n) = xnumom(n) + xnumomdiv
	if(n.ge.10.and.n.le.41) xnumom(n) = xnumom(n) + xnumomsol

	xnupart(n) = xnupartz
25	CONTINUE 	
	N = NCAP
	IF(IMPURITY.GT.0) THEN
      FZINTRIN = YNZ(N-1)/YN(N-1)
	FZINJECT = YNA(N-1)/YN(N-1)
	ENDIF
C	********PUMP OPTION  outer divertor********************************************
	IF(IPUMP.EQ.1) THEN
		XNOD = 0.95*XNOD
		XNODCOLDD = 0.95*XNCOLDD
	ENDIF 
	zcarb(n-1) = 4.0
	if(xt(n-1).gt.100.) zcarb(n-1) = 6.0
	tno(n-1) = todd
	CALL ATCOOL2(XT(N-1),YN(N-1),SATZ)
	XZEFF(N-1)=(YN(N-1)+(Zcarb(n-1)**2)*YNZ(N-1)+(ZARG**2)*YNA(N-1))/
	1	(YN(N-1)+Zcarb(n-1)*YNZ(N-1)+ZARG*YNA(N-1))
 
	IF(IPUMP.EQ.1) THEN
		XNOD = XNOD/0.95
	 	XNCOLDD = XNCOLDD/0.95
	ENDIF
	SAT(N-1) = SATZ
	SRAD(N-1) = SRADZ*fradiv
	SPART(N-1) = XNUPARTZ*YN(N-1) + GPART(N-1)
	xnumom(n-1) = xnumomz + xnumomdiv
	xnupart(n-1) = xnupartz
	
 	

C*********************INTEGRAL BALANCES****************************
	sumpart2 = 0.0
	sumpart1 = 0.0
	sumheat2 = 0.0
	sumheat1 = 0.0
	sumrad1 = 0.0
	sumrad2 = 0.0
	sumat2 = 0.0
	sumat1 = 0.0

	DPERP = 0.0
	SBD = 1.0
	IF(IB.EQ.1) SBD = 0.0
c	heat & particle balance terms
c***note that the outer SOL-DIV is #1 and the inner is #2************
	IF(ITRANS.EQ.1) DPERP = ALB*XT(1)/(16.*BFIELD)

	SUMPART2=SUMPART2 + (SPART(1)-DPERP*YN(1)/
	1((DELN*coef3(1))**2))*DPSI(1)*coef3(1)
	1  +(VDRIFT(1)*YN(1)/(DELN*coef3(1)))*DPSI(1) +
	1	SBD*SEXB*DPSI(1)*coef3(1)
	SUMHEAT2 = SUMHEAT2 + (QHEAT(1))*DPSI(1)*coef3(1) 
	SUMRAD2 = SUMRAD2 + (SRAD(1))*DPSI(1)*coef3(1)
	SUMAT2 = SUMAT2 + (SAT(1))*DPSI(1)*coef3(1)


	IF(ITRANS.EQ.1) DPERP = ALB*XT(2)/(16.*BFIELD)
 	SUMHEAT2 = SUMHEAT2 + (QHEAT(2))*DPSI(2)*coef3(2)	
 	SUMRAD2 = SUMRAD2 + (SRAD(2))*DPSI(2)*coef3(2)
	SUMAT2 = SUMAT2 + (SAT(2))*DPSI(2)*coef3(2)
	SUMPART2 = SUMPART2+(SPART(2)-DPERP*YN(2)/(DELN*coef3(2))**2)*
     1	DPSI(2)*coef3(2) 
	1  +(VDRIFT(2)*YN(2)/(DELN*coef3(2))+ SBD*SEXB)*DPSI(2)*coef3(2)
	do 30 n = 3, ncap-2
	  
	if(n.lt.nstagq) then
     		SUMHEAT2 = SUMHEAT2 + (QHEAT(N))*DPSI(N)*coef3(n)
		SUMRAD2= SUMRAD2 + (Srad(N))*DPSI(N)*coef3(n)  
		SUMAT2= SUMAT2 + (SAT(N))*DPSI(N)*coef3(n)
	endif	 
	if(n.ge.nstagq) then
		SUMHEAT1 = SUMHEAT1 + (QHEAT(N))*DPSI(N)*coef3(n)
		SUMRAD1= SUMRAD1 + (Srad(N))*DPSI(N)*coef3(n)  
		SUMAT1= SUMAT1 + (SAT(N))*DPSI(N)*coef3(n)
	endif

	if(n.ge.nstagp) goto 22
	IF(ITRANS.EQ.1) DPERP = ALB*XT(N)/(16.*BFIELD)
	SUMPART2 = SUMPART2 + SPART(N)*DPSI(N)*coef3(n) 
	IF(N.LE.10) SUMPART2 = SUMPART2
     1  +(VDRIFT(N)*YN(N)/(DELN*coef3(n)))*DPSI(N)*coef3(n)
     2	+ SBD*SEXB*DPSI(N)-DPERP*YN(N)/
     3	((DELN*coef3(n))**2)*DPSI(N)*coef3(n)
c	IF(N.GT.10.AND.N.LT.41)	SUMPART2 = SUMPART2
c     1  -	ABS(VEXB(N)*YN(N)/(coef3(n)*DELN)*DPSI(N)*coef3(n)
	IF(N.GE.41)	SUMPART2 = SUMPART2
     1  +(VDRIFT(N)*YN(N)/(DELN*coef3(n)))*DPSI(N)*coef3(n)
     2	-DPERP*YN(N)/((DELN*coef3(n))**2)*DPSI(N)*coef3(n)
	IF(N.GT.41) SUMPART2 = SUMPART2 + (1.- SBD)*SEXB*DPSI(N)*coef3(n) 
	GOTO 23


22	SUMPART1 = SUMPART1 + SPART(N)*DPSI(N)*coef3(n)
	IF(N.LE.10) SUMPART1 = SUMPART1
     1  +(VDRIFT(N)*YN(N)/(DELN*coef3(n)))*DPSI(N)*coef3(n)
     2	+ SBD*SEXB*DPSI(N)-DPERP*YN(N)/
     3  ((DELN*coef3(n))**2)*DPSI(N)*coef3(n)
c	IF(N.GT.10.AND.N.LT.41)	SUMPART1 = SUMPART1
c     1  -	ABS(VEXB(N)*YN(N)/(DELN/BigH(m)))*DPSI(N)*coef3(n)
	IF(ITRANS.EQ.1) dperp = ALB*xt(n)/(16.*bfield)
	IF(N.GE.41)	SUMPART1 = SUMPART1
     1  +(VDRIFT(N-1)*YN(N)/(DELN*coef3(n)))*DPSI(N)*coef3(n)
     2	-DPERP*YN(N)/((DELN*coef3(n))**2)*DPSI(N)*coef3(n)	
	IF(N.GT.41) SUMPART2 = SUMPART2 + (1.- SBD)*SEXB*DPSI(N)*coef3(n) 

23	CONTINUE
	if (n.eq.30) then
		sumheatsep = sumheat1 + sumheat2-sumrad1-sumrad2-sumat1-sumat2
		sumpartsep = sumpart1 + sumpart2
	endif
	

30	CONTINUE

	N  = NCAP
	IF(ITRANS.EQ.1) dperp = ALB*xt(n)/(16.*bfield)
	SUMHEAT1 = SUMHEAT1 + (QHEAT(N-1))*DPSI(N-1)*coef3(n-1)
 	SUMRAD1 = SUMRAD1 + (SRAD(N-1))*DPSI(N-1)*coef3(n-1)
	SUMAT1 = SUMAT1 + (SAT(N-1))*DPSI(N-1)*coef3(n-1)  
	SUMPART1 = SUMPART1 + SPART(N-1)*DPSI(N-1)*coef3(n-1) 
	1  +	(VDRIFT(N-1)*YN(NCAP-1)/
	3  (DELN*coef3(n-1)))*DPSI(NCAP-1)*coef3(n-1)
	2 -(DPERP*YN(N-1)/((DELN*coef3(n-1))**2))*DPSI(NCAP-1)*coef3(n-1)
     3	+ (1.- SBD)*SEXB*DPSI(N-1)*coef3(n-1)	

35	CONTINUE
	sumheatnet = 0.0 
	do 3599 n = 1,ncap-1
	heatnet(n) = qheat(n)-srad(n)-sat(n)
	sumheatnet = sumheatnet + heatnet(n)*dpsi(n)*coef3(n)
3599	continue 
	chsep = 0.0
	do 3600 n = 2,30
		chsep = chsep -0.5*(bigq(n)+bigq(n-1))*dpsi(n-1)*coef3(n-1) 
3600	continue 
	delheat1 = sumheat1-sumrad1-sumat1
	delheat2 = sumheat2-sumrad2-sumat2
	delheatot = delheat1 + delheat2
	delpart1 = sumpart1
	delpart2 = sumpart2
	delpartot = delpart1 + delpart2
	delmom1 = sumom1
	delmom2 = sumom2

	tq = sumheat1+sumheat2-sumrad1-sumrad2-sumat1-sumat2
	tm = sumom1 + sumom2
	sumom = tm
	tp =sumpart1 + sumpart2	
	tqsep = sumheatsep
	tpsep = sumpartsep
	tmsep = sumomsep

c**************NEW ITERATION 7/18/2010***********************
c	first 5 iterations on a symmetric divertor
c****heat outfluxes balance net heat sources***************
	if(itert.eq.1) bigqB(1) = 0.5*delheatot 
	if(itert.gt.1) bigqB(1) = -1.0*delheatot+bigQB(51)
	if(bigqB(1).gt.0.0) bigqB(1) = -0.5*delheatot
	bigq(1) = bigqB(1)/delrm(1)
	
c******particle outfluxes balance net particle sources*******
	if(itert.eq.1) bigamB(1) = 0.5*delpartot
	if(itert.gt.1) bigamB(1) = -1.0*delpartot+bigamB(51)
	if(bigamB(1).gt.0.0) bigamB(1) = -0.5*delpartot
	bigam(1) = bigamB(1)/delrm(1)
c*********sheath boundary condition****************************
	
      zt(1) = bigQ(1)/(xk*gamsheathin*bigam(1))
 	csin=sqrt(2.*xk*zt(1)/xmas1)
	zn(1) = -1.0*bigam(1)/csin

c********sheath boundary condition on flow velocity************
	bigMom(1) = 4.*xk*zt(1)*zn(1)

c************************************************************
	 	 
C	SOLVE FOR THE PARTICLE,HEAT & MOMENTUM FLUXES
	
C****solve for product of bigam & bigq with coef3***************
	bigamB(1) = delrm(1)*bigam(1)
	bigQB(1) = delrm(1)*bigQ(1)	
3999	do 4000 n=1,ncap-1
	
	VDRIFT(N) = (VEXB(N)+VDRDB(N)+vdia(n))
	if(n.lt.10.or.n.gt.41) vdrift(n) = -1.*(ABS(vdrdb(n))+
	1  ABS(vexb(n))+ABS(VDIA(N)))

	c3 = coef3(n)
	BIGAMB(N+1) = BIGAMB(N) +
     1	c3*DPSI(N)*(SPART(N) + (VDRIFT(N)*YN(N)/(c3*DELN)))
	IF(N.LT.10) BIGAM(N+1) = BIGAM(N+1)	+ SBD*SEXB*DPSI(N)*c3
	IF(N.GT.41) BIGAM(N+1)=BIGAM(N+1) + (1.-SBD)*SEXB*DPSI(N)*c3 
    
     	bigQB(n+1) = bigQB(n) + c3*dpsi(n)*(qheat(n)-sat(n)-srad(n))	
4000	continue
	if(bigamB(51).lt.0.0) bigamB(51) = 0.5*delpartot
	if(bigqB(51).lt.0.0) bigqB(51) = 0.5*delheatot
c*******construct bigam and bigQ*********************************
	bigam(1) = bigamB(1)/delrm(1)
	bigQ(1) = bigQB(1)/delrm(1)
	bigam(ncap) = bigamB(ncap)/delrp(ncap-1)
	bigQ(ncap) = bigQB(ncap)/delrp(ncap-1)
	do 4001 n = 1,ncap-2	
	bigam(n+1) = bigamB(n+1)/delrp(n)
	bigQ(n+1) = bigQB(n+1)/delrp(n)
4001	continue

c*****************MOMENTUM INTEGRALS**************************************
	 
C	coupled first-order momentum flux-velocity formulation 5/28/08
	 sumom1 = 0.0
	 SUMOM2 = 0.0

	DO 159 N = 1, NCAP-1
	
	IF(N.LT.NSTAGP)
     1	SUMOM2 = SUMOM2 - 0.5*DPSI(N)*XMAS1*XNUMOM(N)*
 	2				(BIGAM(N)+BIGAM(N+1))*coef3(n)
	IF(N.GE.NSTAGP)
     1	SUMOM1 = SUMOM1 - 0.5*DPSI(N)*XMAS1*XNUMOM(N)*
 	2				(BIGAM(N)+BIGAM(N+1))*coef3(n)
	if(n.eq.30) sumomsep = sumom1+sumom2 
159	CONTINUE
	tmsep = sumomsep
	tm = sumom1 + sumom2 

	bigmomB(1) = bigmom(1)*delrm(1)
	do 4090 n = 1,ncap-1
	bigmomB(n+1) = bigmomB(n) - 0.5*DPSI(N)*XMAS1*XNUMOM(N)*
 	2				(BIGAM(N)+BIGAM(N+1))*coef3(n)
4090	continue
		do 4091 n = 1,ncap-1
	bigmom(n) = bigmomB(n)/delrm(n)
4091	continue

	bigmom(51) = bigmomB(51)/delrp(50)
c************evaluate outer density & temp from sheath conditions*********
 
	zt(51)= bigQ(51)/(xk*gamsheathout*bigam(51))
	csout=sqrt(2.*xk*zt(51)/xmas1)
	zn(51)= bigam(51)/csout
	bigMomout = 4.*xk*zt(51)*zn(51)	


C**********************TEMPERATURE************************************
c	coupled first-order heat flux-temperature formulation 5/19/08

c*******new solution for T 6/17/2010***********************************
                                                                                                                                                             	SOL(1) = (zT(1))**3.5
     	VTHERM(1) = SQRT(2.*xk*zT(1)/XMAS1)
	xtold(1) = zt(1)
	DO 6155 N = 2,Nstag
	xtold(n) = zt(n)
	SOL(N) = SOL(N-1) 
	1	- (7./(2.*XKAPPA))*0.5*(bigQ(N)+bigQ(n-1))*dpsi(n-1)
c	2	+ (7./(2.*XKAPPA))*1.25*xk*(bigam(N)*zt(n)
c     3	+bigam(n-1)*zt(n-1))*dpsi(n-1)
	
	zT(N) =(SOL(N))**(2./7.)
6155	CONTINUE
	ztstag2=zt(nstag)
c	goto 7500
	mstag = ncap - nstag
	SOL(NCAP) = (ZT(51))**3.5
	DO 6165 M = 1,MSTAG
	N = NCAP-M
	SOL(N) = SOL(N+1) 
	1	+ (7./(2.*XKAPPA))*0.5*(bigQ(N)+bigQ(n+1))*dpsi(n)
c	2	- (7./(2.*XKAPPA))*1.25*xk*(bigam(N)*zt(n)
c     3	+bigam(n+1)*zt(n+1))*dpsi(n)
	zT(N) =(SOL(N))**(2./7.)
6165	CONTINUE  	
	ztstag1 = zt(nstag)
	zt(nstag) = 0.5*(ztstag1 + ztstag2)
	sol(nstag) = (zt(nstag))**3.5



c************************************************************************

	

7500	CONTINUE

C*********************GRADB & CURVATURE DRIFTS*******************************
	BDIR = 1.0
	IF(IB.EQ.1) BDIR = -1.0 
C	ASSIGN THETA VALUE TO EACH SOL LOCATION
	xpi = 3.14159 
	DO 290 N = 11,25
	THET(N) = XPI/2. + (25.-N+0.5)*XPI/15.
290	CONTINUE
	THET(11) = THET(12)
	DO 291 N = 26,40
	THET(N) = XPI/2. - (N-26.+0.5)*XPI/15.
291	CONTINUE
	THET (40) = THET(39)	 
C	PARALLEL PLASMA CURRENT
C	XJ(1) = -1.0*EQ*YN(1)*VTHERM(1)*(1.  
C     1 -0.5*SQRT(XMAS1/(XPI*XMASEL))*EXP(ESPOT(1)/XT(1))*
C     2	(1.-SELEC**2))

c	no drifts if igradb = 0
	IF(IGRADB.EQ.0) GOTO 311
	ALPHA1 = ATAN(ABS((RSEP1-RX)/(ZX-ZSEP1)))	
	ALPHA2 = ATAN(ABS((RSEP2-RX)/(ZX-ZSEP2)))
	xjb(1) = 0.0
	DO 300 N = 1, NCAP-1
C	NR dot DOWNWARD UNIT VECTOR
	IF(N.LT.10) then
		DOTRZ(N) = -1.*COS(ALPHA1)
		dothetz(n) = -1.*sin(alpha1)
	endif
	if(n.eq.10) then
		DOTRZ(N) = 	sqrt(15./16.)
		dothetz(n)= -0.25  
	 endif	 

c	IF(N.GE.10.AND.N.LT.23) then
c		DOTRZ(N) = 0.0
c		dothetz(n) = -1.0
c	endif
c	if(n.eq.23) then
c		dothetz(n) = -0.75
c		dotrz(n)   = -1.*sqrt(7./16.)
c	endif
c	if(n.eq.24) then
c		dothetz(n) = -0.5
c		dotrz(n) = -1.*sqrt(3./4.)
c	endif
	IF(N.EQ.25) then
		DOTRZ(N) = 	-1.*sqrt(15./16.)
		dothetz(n)= -0.25
	endif 
	if(n.eq.26)	then
		dotrz(n) = dotrz(25)
		dothetz(n) = -1.0*dothetz(25)		
	endif
	IF(N.GT.26.AND.N.LT.41) then
	
		gn = 2.*(n-26.)+1.
	
		DOTRZ(N) = -1.*SIN(thet(N))
		dothetz(n) = cos(thet(N))
c	symmetrized
		dotrz(ncap-n) = dotrz(n)
		dothetz(ncap-n)	= -1.*dothetz(n)
	

	endif
	if(n.eq.41) then
		 DOTRZ(N) = sqrt(15./16.)
		dothetz(n)= 0.25   
	endif	 
	IF(N.GT.41) then
		DOTRZ(N) = -1.*COS(ALPHA1)
		dothetz(n) = sin(alpha1)
	endif
	epsilon = 1.0
	if(n.gt.10.and.n.lt.41) then
		if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
		epsilon = 1./bigH(m)
	endif 

	

 
	ZLNM(N) = 1./(epsilon*DELN)
	ZLTM(N) = 1./(epsilon*DELN)
	ZLPM(N) = 1./(epsilon*DELN)
	
	
	
	
 


		
	eq = 1.6e-19
300	CONTINUE

	xjb(1) = 0.
	do 301 n =1,ncap-1
	vdrz(n) =3.*bdir*xk*xt(n)/(eq*abs(rmajor*bfield))
	vdrdb(n) = vdrz(n)*dotrz(n)	
	deljb(n) = 0.0
	if(n.gt.10.and.n.lt.41) deljb(n) =      		
     1	 2.*yn(n)*eq*vdrz(n)*dotrz(n)*zlpm(n)*dpsi(n)
	xjb(n+1) = xjb(n) + deljb(n)
301	continue


c************************************************************
c	7/31/08 adjust gradB current to integrate to 0 around SOL	
	do 305 n = 11,40
	deljb(n) =  deljb(n) -  xjb(42)/30.  
305	continue
	xjb(1) = 0.0
	do 306 n = 2,ncap
	xjb(n) = xjb(n-1) + deljb(n-1)
306	continue
c*************************************************************
C	PARALLEL PARTICLE FLOW	& CURRENT
 	DO 310 N = 1, 50
	gamdb(n+1) = yn(n)*vdrz(n)*dothetz(n)*betag
	BIGAM(N+1) = BIGAM(N+1) + GAMDB(N)
C	XJ(N+1) = XJ(N+1) + 2.*EQ*GAMDB(N)
c	radial drift velocity
	vDrdb(n) = vdrz(n)*dotrz(n)
310	CONTINUE
	
C*****************ELECTROSTATIC POTENTIAL*****************************
311	SELEC = 0.0
	EPOTBIASIN = 0.0
	EPOTBIASOUT = 0.0
	xpi =  3.14159

	GOTO 325
C	******NEW FEBS SOLUTION OF DELJ=0 & ELECTRON MOM EQ*******7/24/08
c	ESPOT(1) = ESPOTBIASIN -1.0*XT(1)*LOG(0.5*SQRT(XMAS1/(XPI*XMASEL))
c	1  *(1.-SELEC))
c    2	/(1.0 - XJ(1)/(YN(1)*EQ*SQRT(2.*XK*XT(1)/XMAS1)))
c	ESPOT(NCAP) = ESPOTBIASOUT -1.0*XT(NCAP-1)*
c     1	LOG(0.5*SQRT(XMAS1/(XPI*XMASEL))*(1.-SELEC))
c     2 /(1.0 - XJ(NCAP)/(YN(NCAP-1)*EQ*SQRT(2.*XK*XT(NCAP-1)/XMAS1)))
c	DO 370 N =2,NCAP-1
c      APCOEF(N) = -1.0*(XT(N)**1.5)/DPSI(N)
c	AXCOEF(N)=((XT(N)**1.5)/DPSI(N) +(XT(N-1)**1.5)/DPSI(N-1))
c	AMCOEF(N) = -1.*(XT(N-1)**1.5)/DPSI(N-1)
c	ASOR(N) = -1.71*((XT(N)**1.5)*(XT(N+1)-XT(N))/
c     1	(0.5*(DPSI(N)+DPSI(N+1)))	
c	2	-(XT(N-1)**1.5)*(XT(N)-XT(N-1))/(0.5*(DPSI(N)+DPSI(N-1))))

c	ASOR(N) = ASOR(N) - ((XT(N)**2.5)*(YN(N+1)-YN(N))/
c     1		(0.5*YN(N)*(DPSI(N)+DPSI(N+1)))-
c	2 (XT(N)**2.5)*(YN(N)-YN(N-1))/(0.5*YN(N-1)*(DPSI(N)+DPSI(N-1))))
c	IF(IGRADB.EQ.1) 
c     1 ASOR(N) = ASOR(N) + 3.*(VDRZ(N)*DPSI(N)+VDRZ(N-1)*DPSI(N-1))*
c     2	DOTRZ(N)/(RMAJOR*BFIELD*SIGM0)
c370	CONTINUE 
c	ASOR(2) = ASOR(2) - AMCOEF(2)*ESPOT(1)
c	ASOR(NCAP-1) = ASOR(NCAP-1) - APCOEF(NCAP-1)*ESPOT(NCAP)  
c	BIGA(2) = APCOEF(2)/AXCOEF(2)
c	ALF(2) = ASOR(2)/AXCOEF(2)
c	DO 375 N =3,NCAP-1
c	BIGA(N) = APCOEF(N)/(AXCOEF(N) - AMCOEF(N)*BIGA(N-1))
c      ALF(N)= (ASOR(N)-AMCOEF(N)*ALF(N-1))/
c	1	(AXCOEF(N) - AMCOEF(N)*BIGA(N-1))
c375	CONTINUE
c	ESPOT(NCAP-1) = ALF(NCAP-1)
c	DO 380 M = 1,NCAP-3
c	N = NCAP-1-M
c	ESPOT(N) = BIGA(N)*ESPOT(N+1) + ALF(N)
c380   CONTINUE		 
C	6/23/08 NOTES
325	ESPOTBIAS = 0.
C	ESPOT(1) = ESPOTBIAS -1.0*XT(1)*LOG(0.5*SQRT(XMAS1/(XPI*XMASEL))*
C	1  (1.-SELEC))
C     2	/(1.0 - ABS(XJ(1))/(YN(1)*EQ*SQRT(2.*XK*XT(1)/XMAS1)))

c****************************
	eta0 = 0.1
c*****************************
C	INTEGRALs IN POTENTIAL EQ
	
	sumlng = 0.0
	DO 328 N = 1,NCAP-1 
c	log density
	sumlng = sumlng + 0.5*(zt(n)+zt(n+1))*log(zn(n+1)/zn(n))
328	continue
c	DLN INTEGRAL
	SUMJB(1) = 0.0
	do 329 n = 2,ncap
 

	IF(N.EQ.2) 	DELSUM(N) =	 XT(1)*LOG(zN(2)/zN(1))
	IF(N.GT.2.AND.N.LT.NCAP) DELSUM(N) = 
     1	XT(N-1)*LOG(zN(N+1)/zN(N))
	IF(N.EQ.NCAP) DELSUM(N) =  
     1	XT(NCAP-1)*LOG(zN(NCAP)/zN(NCAP-1)) 
	IF(N.EQ.2) SUMLN(N) = DELSUM(N)
	IF(N.GT.2) SUMLN(N) = SUMLN(N-1) + DELSUM(N)

C	GRAD_B DRIFT CURRENT INTEGRAL IN POTENTIAL EQ (RAD & PARALLEL DRIFTS)
	IF(IGRADB.EQ.0) GOTO 330
	
	IF(N.LE.10) SUMJB(N) = SUMJB(N-1)
	2 + ETA0*DPSI(N-1)*EQ*(GAMDB(N)+GAMDB(N-1))/(XT(N-1)**1.5)

	if(n.gt.10.and.n.lt.41)
   	1SUMJB(N) = SUMJB(N-1) + eta0*DPSI(N-1)*0.5*(XJB(N)+XJB(N-1))/
 	1			(XT(N-1)**1.5)
	2 + ETA0*DPSI(N-1)*EQ*(GAMDB(N)+GAMDB(N-1))/(XT(N-1)**1.5)

	IF(N.GE.41) SUMJB(N) = SUMJB(N-1)
     2 + ETA0*DPSI(N-1)*EQ*(GAMDB(N)+GAMDB(N-1))/(XT(N-1)**1.5)
329	continue
C	ETA DPSI/T^3/2 INTEGRAL IN POTENTIAL EQ
330	SumCUR(1) = 0.0 
	DO 331 N = 2,NCAP 
	sumCUR(N) = sumCUR(N-1) +  eta0*DPSI(N-1)/(XT(N-1)**1.5)
331	CONTINUE 
C	TOTAL CURRENT DUE TO GRAD_B DRIFTS
	GRADBCUR(1)	= 0.0
	DO 332 N = 2, NCAP
	GRADBCUR(N) = XJB(N) + 2.*EQ*GAMDB(N)
332	CONTINUE 
   	 
C*****************SOLUTION FOR JIN*****************************************
c	goto 418 
411	yjin = xjin
	CONST =1.71*(zT(NCAP)-zT(1)) + SUMLNG + SUMJB(NCAP)
	3		-(VBIASOUT-VBIASIN)
	eq = 1.6e-19
	CSIN = SQRT(2.*XK*zT(1)/XMAS1)
	CSOUT = SQRT(2.*XK*zT(NCAP)/XMAS1)
	vx = 0.5*sqrt(xmas1/(xpi*xmasel))*(1-selec)
412	xx = (XJIN)/(zN(1)*EQ*CSIN)
c	if(xx.le.-1.0) xx = -.95 
	xin = zt(1)*log(vx/(1. + xx))
c	if(itert.eq.1) then
c	sqout = sqrt(2.*xk*zt(ncap)/xmas1)
c	tp = sumpart1 + sumpart2
c     zn(ncap) = (tp - zn(1)*sqin)/sqout
c	endif
	zjout = (XJIN)+XJB(NCAP) 
	vb = (1.-(zjout/(zN(NCAP)*EQ*CSout)) )
	xout = zt(ncap)*log(vx/vb)
	SUMOUT = SUMCUR(NCAP)
	xjinold = xjin
 
c	XJIN = (CONST-xout+xin)/SUMOUT 
      xjin = 20000
 	DELXJIN = ABS((XJIN-XJINOLD)/XJIN)
	delxjin = 0.01
	IF(DELXJIN.LT. 0.05) goto 418
	GOTO 412

418	CONTINUE
	
c	xjin = yjin + 0.5*(xjin-yjin)

C**********************CURRENT*********************************
	xj(1) = xjin
	DO 419 N = 2,NCAP
	XJ(N) = XJIN + XJB(N) + 2.*EQ*GAMDB(N) 
419	CONTINUE
C************************potential******************************
	SUMJ(1) = 0.0
	ESPOT(1) = XIN + vbiasin
	do 340 n = 2, ncap

    
	IF(N.NE.NCAP)	
     1	ESPOT(N) = ESPOT(1) + 1.71*(zT(N)-zT(1)) + 
     1				SUMLN(N) - SUMJB(N) - XJIN*SUMCUR(N)
     	IF(N.EQ.NCAP) 
     1		ESPOT(N) = ESPOT(1) + 1.71*(zT(N)-XT(1)) + SUMLN(N)
     2		- SUMJB(N) - XJIN*SUMCUR(N) 
	 
340	CONTINUE


	

C************************************************************************

C**************************EXB*********************************************
420	 CONTINUE
c	parallel electric field 
	 do 421 n = 1, ncap-1
	EPAR(N) = -1.*((ESPOT(N+1)-ESPOT(N))/(DPSI(N)))
421	continue
	if(ioptexb.eq.0) goto 435
C	EXB DRIFTDUE TO PARALLEL ELECTRIC FIELD
	DO 425 N = 1, NCAP-1
	VEXB(N) = -1.*BDIR*epar(n)/BFIELD
425	CONTINUE
	SEXB = 0.0
	IF(IB.EQ.0) THEN
	DO 426 N = 42, NCAP-1
	epd = coef3(n)
	SEXB = SEXB + ABS(VEXB(N)*YN(N)/(EPD*DELN))*dpsi(n)/DIVDEPTH1
426	CONTINUE
	ENDIF
	GOTO 428
	DO 427 N = 1, 9
	epd = coef3(n)
	SEXB = SEXB + ABS(VEXB(N)*YN(N)/(epd*DELN))*dpsi(n)/DIVDEPTH2
427	CONTINUE
428	CONTINUE


c	a fraction fexb of the ions drifting from outer leg reach inner leg
	fexb = 1.0
	sexb = fexb*sexb
	
c	EXB DRIFT DUE TO RADIAL ELECTRIC FIELD, BOHM TRANSPORT
	DO 431 N = 11,40
c	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
c	1	(FLUXHEAT/(XK*XT(N)))
     	delt = deltin

	if(n.ge.11.and.n.le.33) m = 34 - n
	if(n.ge.34.and.n.le.40) m = 64 - n

	GAMErXB(N) = BDIR*BETAG*(YN(N)*espot(n)/BFIELD)/(DELT*coef3(n))
	BIGAM(N) = BIGAM(N) + 0.5*(GAMErXB(N)+GAMErXB(N-1))
	ER(N) = (ESPOT(N)/(DELT/bigH(m)))
431	CONTINUE
	DO 432 N= 1,10
c	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
c	1	(FLUXHEAT/(XK*XT(N))) 
	delt = deltin
C	TRANSITION FROM OUTWARD Er IN REGION 10 TO INWARD Er IN REGION 1
	TRANS = -1. + 2.*(N-1)/9.
	epd = coef3(n)
	GAMErXB(N) = TRANS*BDIR*BETAG*(YN(N)*espot(N)/BFIELD)/(EPD*DELT)
	GAMErXB(NCAP-N) = TRANS*BDIR*BETAG*(YN(NCAP-N)*espot(NCAP-N)/
     1	BFIELD)/(EPD*DELT)
C	ER(N) = (3.*XT(N)/DELT)*TRANS
	ER(N) = (ESPOT(N)/(EPD*DELT))*TRANS
C	ER(NCAP-N) = (3.*XT(NCAP-N)/DELT)*TRANS 
	ER(NCAP-N) = (ESPOT(NCAP-N)/(EPD*DELT))*TRANS 
432	CONTINUE 
	DO 433 N=2,10
	BIGAM(N+1) = BIGAM(N) + 0.5*(GAMErXB(N)+GAMErXB(N-1))
c	write(*,*) n,bigam(n+1),bigam(n)
	BIGAM(NCAP-N) = BIGAM(NCAP-N) + 0.5*(GAMErXB(N)+GAMErXB(NCAP-N-1))
433	CONTINUE	 

C**************************************************************************
435	continue
	DO 436 N = 1,10	
	epd = coef3(n)
	TRANS = -1. + 2.*(N-1)/9.
C	ER(N) = (3.*XT(N)/DELT)*TRANS*BDIR
	ER(N) = (ESPOT(N)/(EPD*DELT))*TRANS
C 	ER(NCAP-N) = (3.*XT(NCAP-N)/DELT)*TRANS*BDIR
 	ER(NCAP-N) = (ESPOT(NCAP-N)/(EPD*DELT))*TRANS 

 436	CONTINUE
	DO 437, N=11,40
	if(n.ge.11.and.n.le.33) m = 34 - n
	if(n.ge.34.and.n.le.40) m = 64 - n

C	ER(N) =  (3.*XT(N)/DELT)
	ER(N) = (ESPOT(N)/(DELT*coef3(n)))
 437	CONTINUE
C********************************end ExB **********************************

c*****************diamagnetic drift*****************************************
	if(ioptdia.eq.0) goto 441
c	used only for particle source/sink, not for current
	vdia(1) = -1.0*bdir*xk*(yn(2)*xt(2)-yn(1)*xt(1))/
	1	(yn(1)*eq*betag*bfield*dpsi(1))
	do 440 n = 2, ncap-2
	vdia(n) = -1.0*bdir*xk*(yn(n+1)*xt(n+1)-yn(n)*xt(n))/
	1	(yn(n)*eq*betag*bfield*dpsi(n))
440	continue
      vdia(ncap-1) = -1.0*bdir*xk*(yn(ncap-1)*xt(ncap-1)-
     2  yn(ncap-2)*xt(ncap-2))/(yn(ncap-1)*eq*betag*bfield*dpsi(ncap-1))
441	continue
C**********************DENSITY & VELOCITY*********************************

c	new solution of momentum and particle balance eqs for n and v 6/16/08
	csin = sqrt(2.*xk*zt(1)/xmas1)
	csout = sqrt(2.*xk*zt(ncap)/xmas1)
c	write(*,*) bigam
	do 350 n = 2,ncap-1

	bq = -bigmom(n)
	aq = 2.*xk*zt(n)
	cq = xmas1*(bigam(n)**2)
	SR(N) = 4.*aq*cq/(bq**2)
c	if(cq.gt.bq*bq/(4.*aq)) goto 350

	zn(n) = (-1.*bq/(2.*aq))*(1+sqrt(1.-4.*aq*cq/(bq**2)))
	ynm(n)=	(-1.*bq/(2.*aq))*(1-sqrt(1.-4.*aq*cq/(bq**2)))
	if(n.eq.ncap) zn(n) = (-1.*bq/(2.*aq))
 	if(sr(n).ge.1.0) zn(n) = (-1.*bq/(2.*aq)) 
350	continue

	VEE(1) = -1.0*SQRT(2.*XK*zT(1)/XMAS1)
	do 355 n = 2,ncap
	vee(n) = bigam(n)/zn(n)
355	continue

c	average density and temp in mesh interval
	do 359 n = 1,ncap-1  
	yn(n) = 0.5*(zn(n)+zn(n+1))
	xt(n) = 0.5*(zt(n)+zt(n+1)) 
359	continue
c************************************************
	

360	continue
C*********************CONVERGENCE*****************************
	write(6,555) itert,zt(1),zt(51),zn(1),zn(51),zn(30),zt(30)
555		format (i5,6e10.4)

C	CHECK CONVERGENCE

	  if(itert.eq.1) goto 213
	IFCON = 0
	EPCON = 0.02

c	iterate on particle and energy fluxes out
	gamout = bigam(51)
	qout = bigQ(51)
	x = (gamout-gamoutold)/gamout
	if(abs(x).gt.epcon) goto 210
	x = (qout-qoutold)/qout
	if(abs(x).gt.epcon) goto 210

	goto 225
210	con = 0.0
	if(itert.gt.1) con = 0.9
c	bigam(1) = con*gaminold+(1.-con)*gamin
	bigam(51)= con*gamoutold+(1.-con)*gamout
c	bigQ(1) = con*qinold+(1.-con)*qin
	bigQ(51)= con*qoutold+(1.-con)*qout
213	gamoutold = bigam(51)
	qoutold = bigQ(51)

211	continue 
C	211 end of do loop for D n,T, gamma.  Particle and heat fluxes converged.
225	continue 
	

C******THE CALCULATION IS NOW CONVERGED WITH FIXED INPUT IMPURITY DISTRIBUTIONS
C******NEXT CARBON AND ARGON DISTRIBUTIONS ARE CALCULATED	USING THESE D DISTRIBUTIONS
C******THEN THE CALCULATION IS LOOPED BACK THROUGH THE ABOVE DEUTERIUM CALCULATION AND
C******THE IMPURITY CALC 10 TIMES ON THE "GOTO 1000 IF IMPURITY.LT.10" STATEMENT AT THE END

C*************CARBON IMPURITIES*****************************************

C	SPUTTERING YIELDS CARBON (FPP p 324)
C	TRANPORT MODEL KIELHACKER NUCLEAR FUSION, 31, 535, 1991
	ETF1=447.
	ETH1=30.
	Q1 = 0.1
	ETFZ=5687.
	ETHZ=42.
	QZ=1.5
C	ENERGY = T + 3T SHEATH ACCELERATION
	EPIN1 = 4.*XT(1)/ETF1
	EPINZ = 4.*XT(1)/ETFZ
	DELTA1IN = ETH1/XT(1)
	GIN1 =(1.-DELTA1IN**.6667)*((1-DELTA1IN)**2)
	DELTAZIN = ETHZ/XT(1)
	GINZ =(1.-DELTAZIN**.6667)*((1-DELTAZIN)**2)
	SN1IN =3.441*SQRT(EPIN1)*LOG(EPIN1+2.718)/(1.+6.355*SQRT(EPIN1)
	1		+EPIN1*(6.882*SQRT(EPIN1)-1.708))
	SNZIN =3.441*SQRT(EPINZ)*LOG(EPINZ+2.718)/(1.+6.355*SQRT(EPINZ)
     1		+EPINZ*(6.882*SQRT(EPINZ)-1.708)) 
	YIZIN = Q1*SN1IN*GIN1
	YZZIN = QZ*SNZIN*GINZ
	EPOUT1 = 4.*XT(NCAP-1)/ETF1
	EPOUTZ = 4.*XT(NCAP-1)/ETFZ
	DELTA1OUT = ETH1/XT(NCAP-1)
	GOUT1 =(1.-DELTA1OUT**.6667)*((1-DELTA1OUT)**2)
	DELTAZOUT = ETHZ/XT(NCAP-1)
	GOUTZ =(1.-DELTAZOUT**.6667)*((1-DELTAZOUT)**2)
	SN1OUT =3.441*SQRT(EPOUT1)*LOG(EPOUT1+2.718)/(1.+6.355*SQRT(EPOUT1)
	1		+EPOUT1*(6.882*SQRT(EPOUT1)-1.708)))
	SNZOUT =3.441*SQRT(EPOUTZ)*LOG(EPOUTZ+2.718)/(1.+6.355*SQRT(EPOUTZ)
     1		+EPOUTZ*(6.882*SQRT(EPOUTZ)-1.708))) 
	YIZOUT = Q1*SN1OUT*GOUT1
	YZZOUT = QZ*SNZOUT*GOUTZ
	yizin =  0.0001
	yizout = 0.0001
	yzzout = 0.0001
	yzzin  = 0.0001

C	PARTICLE  FLUX
	XMASZ = XMAS(2)
	IF(IMPURITY.EQ.0) THEN
	DO 504 N = 1,NCAP-1
	YNZ(N) = FZINTRIN*YN(N)
504	CONTINUE
	ENDIF

	CSIIN = SQRT(2.*XK*XT(1)/XMAS1)
	CSZIN = SQRT(2.*XK*XT(1)/XMASZ)
 	CSIOUT = SQRT(2.*XK*XT(NCAP-1)/XMAS1)
	CSZOUT = SQRT(2.*XK*XT(NCAP-1)/XMASZ)

c	*******beginning of 575 do loop TO ITERATE CARBON SOLUTION 10 TIMES***
	itertz = 0
	DO 575 MM = 1, 10
505	itertz = itertz+1
	INTER = 0
506	INTER = INTER + 1
c	PRINT *,"MM,ITERt=",MM,ITERt
	SZINI = YN(1)*csiin*yizin/(dpsi(1)+DPSI(2))
	SZINZ =	ynz(1)*cszin*yzzin/(DPSI(1)+DPSI(2))
	SZOUTI= YN(NCAP-1)*csiout*yizout/(dpsi(ncap-1)+DPSI(NCAP-2))	   	
	SZOUTZ= ynz(ncap-1)*cszout*yzzout/(DPSI(NCAP-1)+DPSI(NCAP-2))
	sumz2 =  (szinI)*(dpsi(1)+DPSI(2))
	sumz1 =  (szoutI)*(dpsi(ncap-1)+DPSI(NCAP-2))
	XLOSS = 0.
	IF(MM.EQ.1) GOTO 511
	do 510 n = 1,ncap-1
	epsilon = epd
	if(n.gt.10.and.n.lt.41) then
		if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
		epsilon = 1./bigH(m)
	endif 

	DPERPZ = 0.
	
	IF(ITRANS.EQ.1) dperpz = ALB*xt(n)/(16.*bfield)
	XLOSS = XLOSS +	 dpsi(n)*dperpz*ynz(n)/((epsilon*deln)**2)
	if(n.gt.25) goto 509 
C**********WHEN THE DPERP TRANSPORT IS ON, THESE SUMS GO NEGATIVE ON 2ND ITER,
C****LEADING TO NEGATIVE YNZ(1)*********DON'T USE ITRANS=1********************	 	
	sumz2 = sumz2 - dpsi(n)*dperpz*ynz(n)/((epsilon*deln)**2)
	goto 510
509	sumz1 = sumz1 - dpsi(n)*dperpz*ynz(n)/((epsilon*deln)**2)
510	continue
511	CONTINUE
 	gamma =  (SZINI+SZOUTI)/(XLOSS-SZINZ-SZOUTZ)

	SUM = SZINZ+SZINI+SZOUTZ+SZOUTI 
C 	IF(INTER.EQ.1) GOTO 506
c*********************************
	RCARB = 0.99
c***********************************
	YNZ(1) = SUMZ2/(CSZIN*(1-RCARB))
	YNZ(NCAP-1) = SUMZ1/(CSZOUT*(1-RCARB))
	VEZ(1) = -1.*CSZIN	
	VEZ(NCAP) = CSZOUT
	bigamz(1) = -1.*YNz(1)*cszin*(1.-rcarb)
	bigamz(ncap) =YNz(NCAP-1)*cszout*(1.-rcarb) 
	CRES = (0.457/(1.077+ZEFF) + 0.29*ZEFF)
	CE2 = 1.5*(1.-0.6934/1.3167**ZEFF)
	EP0 = 8.854E-12
	ZEE = ZCARB(1)


C	NEW INTEGRATION 7/3/08
	DO 520 N = 1,NCAP-2

	IF(ITRANS.EQ.1) DPERPZ = ALB*XT(N)/(16.*BFIELD)
	EPSILON = epd
	if(n.gt.10.and.n.lt.41) then
		if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
		epsilon = 1./bigH(m)
	endif 

	
	BIGAMZ(N+1) = BIGAMZ(N) - DPSI(N)*DPERPZ*YNZ(N)/
	1 ((epsilon*DELN)**2)

C	RADIAL EXB and gradB DRIFTs
	VDRIFTZ(N) = VEXBZ(n)+vdrdbcarb(n)+vdiaz(n)
	IF(N.Le.10.or.n.ge.41) 
     1	VDRIFTZ(n) = -1.*(abs(VDRDBCARB(N)) +abs(vexbz(n))
	2	+ abs(vdiaz(n)))
	BIGAMZ(N+1) = BIGAMZ(N+1)
     1  +(VDRIFTZ(N))*YNz(n)/(epsilon*DELN)*DPSI(n)	
	IF(N.LT.10)	BIGAMZ(N+1) = BIGAMZ(N+1) + SBD*SEXBZ*DPSI(n)
	IF(N.GT.41)	BIGAMZ(N+1) = BIGAMZ(N+1) + (1.-SBD)*SEXBZ*DPSI(n)
C	PARALLEL ErXB FLUX
	BIGAMZ(N+1) = BIGAMZ(N+1) + GAMErXBZ(N)
c	parallel gradB flux
	BIGAMz(N+1) = BIGAMz(N+1) + GAMDBcarb(N) 
C	SPUTTERED SOURCE DISTRIBUTED OVER FIRST 2 MESH
	IF(MM.GT.1.AND.N.EQ.1) BIGAMZ(2) = BIGAMZ(2) + DPSI(1)*SZINI
	IF(MM.GT.1.AND.N.EQ.2) BIGAMZ(3) = BIGAMZ(3) + DPSI(2)*SZINI
	IF(MM.GT.1.AND.N.EQ.NCAP-2) BIGAMZ(NCAP-1) = BIGAMZ(NCAP-1) + 
     1	DPSI(NCAP-2)*SZOUTI 
	IF(MM.GT.1.AND.N.EQ.NCAP-1) BIGAMZ(NCAP) = BIGAMZ(NCAP) + 
     1	DPSI(NCAP-1)*SZOUTI
	Z0 = (YNZ(N)/YN(N))*(ZCARB(N)**2)
 	CI1 = ((1.+0.24*Z0)*(1.+.93*Z0))/((1.+2.65*Z0)*(1.+.285*Z0))
	CI2 =1.56*(1+1.414*Z0)*(1.+0.52*Z0)/
	1	((1.+2.65*Z0)*(1.+0.285*Z0)*(Z0+SQRT(0.5*(1.+XMAS1/XMASZ))))
      RES = CRES/(1.9E3*XT(N)**1.5)
	IF(N.NE.1.AND.N.NE.NCAP-1) DT =0.5*XK*(XT(N+1)-XT(N-1))
	IF(N.EQ.1) DT = 0.5*XK*(XT(2)-XT(1))
	IF(N.EQ.NCAP-1) DT = 0.5*XK*(XT(NCAP-1)-XT(NCAP-2))
C			COULOMB LOGARITHM
	Y = SQRT(YN(N))*(ATNUM(1)**2)*ATNUM(2) 
	X = (EP0/EQ)**1.5	
	COULOG(2,1) = LOG(12.*3.1416*(XT(N)**1.5)*X/Y)
C			REDUCED MASS
	XMR = XMAS(1)*XMAS(2)/(XMAS(1)+XMAS(2))
	XMR12 = XMAS(1)*(1.+XMAS(1)/XMAS(2))
	C1 = 1./((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5)) 
	XNUC(2,1)=3.34*(COULOG(2,1)*((ATNUM(1)*ATNUM(2))**2)*1.E-6*YN(N))
     2			/(C1*SQRT(XMR12*1E3)*(XT(N)**1.5))
	XNUZI = XNUC(2,1)
	XNUZI =	YN(N)*(ATNUM(2)**2)*(ATNUM(1)**2)*(EQ**2.5)*COULOG(2,1)
     2		  / (47.25*SQRT(XMR)*(EP0**2)*(XT(N)**1.5))

	IF(N.LT.NCAP-2) VEZ(N+1) = 2.*BIGAMZ(N+1)/(YNZ(N+1)+YNZ(N+2))
	IF(N.EQ.NCAP-2) VEZ(N+1) = BIGAMZ(N+1)/YNZ(N+1)
	IF(N.NE.1.AND.N.NE.NCAP-1) DT =0.5*XK*(XT(N+1)-XT(N-1))
	IF(N.EQ.1) DT = 0.5*XK*(XT(2)-XT(1))
	IF(N.EQ.NCAP-1) DT = 0.5*XK*(XT(NCAP-1)-XT(NCAP-2))
 
	YNZ(N+1) = (XT(N)/XT(N+1))*YNZ(N) + 
     2		XMASZ*(VEZ(N)*BIGAMZ(N)-VEZ(N+1)*BIGAMZ(N+1))/XT(N+1)
     2		-ynz(n)*zcarb(n)*eq*(espot(n+1)-espot(n))/xt(n+1) +
	3		(YNZ(N)/XT(N+1))*(ZCARB(N)**2)*(CE2/ZEFF+CI2)*DT	 -
     3(YNZ(N)/XT(N+1))*(ZCARB(N)**2/XZEFF(N))*EQ*RES*0.5*(XJ(N+1)+XJ(N))
     2*DPSI(N)+CI1*XMASZ*XNUZI*0.5*((YNZ(N)/YN(N))*(BIGAM(N)+BIGAM(N+1))
     5	 - (BIGAMZ(N)+BIGAMZ(N+1)))*dpsi(n)/XT(N+1)
	VEZ(N+1) = 2.*BIGAMZ(N+1)/(YNZ(N+1)+YNZ(N))	 	
 
 
520	CONTINUE
	n  = ncap
	bigamz(n) = bigamz(n-1) +(VDRIFTZ(N-1))*YNz(n-1)/
	1	(EPD*DELN)*DPSI(n-1)
     1	 - DPSI(N)*DPERPZ*YNZ(N)/((EPD*DELN)**2)
	2	 + DPSI(N-1)*SZOUTI
	3	+  GAMErXBZ(N-1) 	+ GAMDBcarb(N-1) 
	4	+ (1.-SBD)*SEXBZ*DPSI(n-1)

C*********************GRADB & CURVATURE DRIFTS CARBON*******************
	
	IF(IGRADB.EQ.0) GOTO 1420 
	BDIR = 1.0
	IF(IB.EQ.1) BDIR = -1.0
C	PARALLEL PLASMA CURRENT

 	DO 1400 N = 1, NCAP-1
	epsilon = epd
		if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
	if(n.ge.11.and.n.le.40) epsilon = 1.0/bigH(m)

	ZLNM(N) = 1./(DELN*epsilon)
	ZLTM(N) = 1./(DELN*epsilon)
	ZLPM(N) = 1./(DELN*epsilon)
	
	
	
 
C	CHANGE IN PARALLEL CURRENT DUE TO RADIAL GRADB & CURV DRIFT CURRENTS

	vdrzcarb(n) =  BDIR*3.*XK*XT(N)/(zcarb(N)*eq*abs(RMAJOR*BFIELD))
	DELJcarb(N)=2.*ynz(n)*zcarb(N)*eq*dpsi(n)*vdrzcarb(n)*
     1	zlpm(n)*dotrz(n)
c	radial drift velocity
	vDrdbcarb(n) = vdrzcarb(n)*dotrz(n)
  	XJ(N+1) = XJ(N+1) + DELJcarb(N)
1400	CONTINUE

C	PARALLEL PARTICLE FLOW
 	DO 1410 N = 1, 50
	gamdbcarb(n) = ynz(n)*vdrzcarb(n)*dothetz(n)*betag
C	bigamz(n+1) = bigamz(n+1) +gamdbcarb(n)
1410	CONTINUE
1420	continue
C**************************EXB CARBON*********************************************

	if(ioptexb.eq.0) goto 535
C	EXB DRIFTDUE TO PARALLEL ELECTRIC FIELD
	DO 525 N = 1, NCAP-1
	VEXBZ(N) = vexb(n)
525	CONTINUE
	SEXBZ = 0.0
	IF(BDIR.EQ.0) THEN
	DO 526 N = 42, NCAP-1
	SEXBZ = SEXBZ + ABS(VEXBZ(N)*YNz(N)/(EPD*DELN))*dpsi(n)/DIVDEPTH1
526	CONTINUE
	ENDIF
	GOTO 528
	DO 527 N = 1, 9
	SEXBZ = SEXBZ + ABS(VEXBZ(N)*YNz(N)/(EPD*DELN))*dpsi(n)/DIVDEPTH1
527	CONTINUE
528	CONTINUE
c	a fraction fexb of the ions drifting from outer leg reach inner leg
	fexb = 1.0
	sexbZ = fexb*sexbZ
		
c	EXB DRIFT DUE TO RADIAL ELECTRIC FIELD, BOHM TRANSPORT
	DO 531 N = 11,40
	FON = XNOV(N)/YN(N)
	CALL CXRCEFITS(6,XT(N),fon,celz,dlzzdt,ZAV)
	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
 	1	(FLUXHEAT/(XK*XT(N))-2.5*FLUXPART)
     	delt = deltin
	if(n.ge.11.and.n.le.33) m = 34 - n
	if(n.ge.34.and.n.le.40) m = 64 - n
  	 epsilon = 1./bigH(m)
	GAMErXBZ(N) = BDIR*BETAG*(YNZ(N)*ESPOT(N)/BFIELD)/(epsilon*DELT)
C	BIGAMZ(N) = BIGAMZ(N) + 0.5*(GAMErXBZ(N)+GAMErXBZ(N-1))
531	CONTINUE
	DO 532 N= 1,10
	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
 	1	(FLUXHEAT/(XK*XT(N))-2.5*FLUXPART)
	delt = 	deltin
	GAMErXBZ(N) = -1.*BDIR*BETAG*(YNZ(N)*ESPOT(N)/BFIELD)/(epd*DELT)
	GAMErXBZ(NCAP-N)=-1.*BDIR*BETAG*(YNZ(NCAP-N)*ESPOT(NCAP-N)/
     1	BFIELD)/(epd*DELT)
532	CONTINUE 

C	DO 533 N=2,10
c	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
c 	1	(FLUXHEAT/(XK*XT(N))-2.5*FLUXPART)	
C	BIGAMZ(N) = BIGAMZ(N) + 0.5*(GAMErXBZ(N)+GAMErXBZ(N-1))
C	BIGAMZ(NCAP-N)=BIGAMZ(NCAP-N)+0.5*(GAMErXBZ(N)+GAMErXBZ(NCAP-N-1))
C533	CONTINUE	 

C**************************************************************************
535	continue
c*****************diamagnetic drift carbon*********************************
	if(ioptdia.eq.0) goto 541
c	used only for particle source/sink, not for current
	vdiaz(1) = -0.5*bdir*xk*(ynz(2)*xt(2)-ynz(1)*xt(1))/
	1	(ynz(1)*eq*betag*bfield)
	do 540 n = 2, ncap-2
	vdiaz(n) = -0.5*bdir*xk*(ynz(n+1)*xt(n+1)-ynz(n)*xt(n))/
	1	(ynz(n)*eq*betag*bfield)
540	continue
      vdiaz(ncap-1) = -0.5*bdir*xk*(ynz(ncap-1)*xt(ncap-1)-
     2	ynz(ncap-2)*xt(ncap-2))/(ynz(ncap-1)*eq*betag*bfield)
 
541	continue	
575	continue
c	*****************end of 575 do loop*********************************


C**********************DENSITY & VELOCITY*********************************
		

	DO 580 N = 1, NCAP-1
	FZCARBON(N) = YNZ(N)/YN(N)
580	CONTINUE	
C***************************************************************************

C********************ARGON IMPURITIES****************************************
	XMASA = AZINJECT*1.673E-27
	ATNUM(3) = 10.
	IF(IMPURITY.EQ.O) THEN
	 DO 604 N = 1,NCAP-1
	YNA(N) = FZINJECT*YN(N)
604	CONTINUE
	ENDIF
	CSIIN = SQRT(2.*XK*XT(1)/XMAS1)
	CSAIN = SQRT(2.*XK*XT(1)/XMASA)
 	CSIOUT = SQRT(2.*XK*XT(NCAP-1)/XMAS1)
	CSAOUT = SQRT(2.*XK*XT(NCAP-1)/XMASA)

c***********************************************************************5555
c	*******beginning of 675 do loop TO ITERATE Argon SOLUTION 10 TIMES***
	iterta = 0
	DO 675 MM = 1, 10
605	iterta = iterta+1
	INTER = 0
606	INTER = INTER + 1
c	PRINT *,"MM,ITERt=",MM,ITERt
C	ARGON SOURCE DISTRIBUTION IN DIVERTOR CHANNELS
	divl1 = 0.0
	divl2 = 0.0
	do 607 n = 1,9
	divl2 = divl2 + dpsi(n)
	divl1 = divl1 + dpsi(ncap-n)
607	continue
C	1 TORR-L/s = 6.7E+19 #/s		
	sumA2 =  1.E19
	sumA1 =  1.E19
	DO 608 N = 1,9
	SA(N) = SUMA2/divl2
	SA(NCAP-N) = SUMA1/divl1
608	CONTINUE
 
	XLOSS = 0.
	IF(MM.EQ.1) GOTO 611
	do 610 n = 1,ncap-1

	DPERPa = 0.
	if(n.le.10) epsilon = epd
	if(n.ge.41) epsilon = epd
	if(n.gt.10.and.n.lt.41) then
		if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
		epsilon = 1./bigH(m)
	endif
	IF(ITRANS.EQ.1) dperpa = ALB*xt(n)/(16.*bfield)
	XLOSS = XLOSS +	 dpsi(n)*dperpa*yna(n)/((epsilon*deln)**2)
	

	if(n.gt.25) goto 609  	
	sumA2 = sumA2 + SA(N)*DPSI(N) 
	SUMA2 = SUMA2 - DPSI(N)*DTRAN*YNA(N)/((EPSILON*DELN)**2)
	goto 610
609	sumA1 = sumA1 + SA(N)*DPSI(N) 
	1	- DPSI(N)*DTRAN*YNA(N)/((EPsilon*DELN)**2)

610	continue
611	continue 
	SUMA = SUMA1 + SUMA2





c*********************************
	RARG = 0.99
c***********************************
	YNa(1) = SUMa2/(CSaIN*(1-RARg))
	YNa(NCAP-1) = SUMa1/(CSaOUT*(1-RARg))
	VEa(1) = -1.*CSaIN	
	VEa(NCAP) = CSaOUT
	bigama(1) = -1.*YNa(1)*csain*(1.-rarg)
	bigama(ncap) =YNa(NCAP-1)*csaout*(1.-rarg) 
	CRES = (0.457/(1.077+ZEFF) + 0.29*xZEFF(n))
	CE2 = 1.5*(1.-0.6934/1.3167**xZEFF(n))
	EP0 = 8.854E-12
	ZEE = Zarg


C	NEW INTEGRATION 7/3/08
	ARGCORE = 0.0
	DO 620 N = 1,NCAP-2

	IF(ITRANS.EQ.1) DPERPa = ALB*XT(N)/(16.*BFIELD)
	EPSILON = epd
	
	if(n.gt.10.and.n.lt.41) then
		if(n.le.33) m = 34 - n
		if(n.ge.34) m = 64 - n
		epsilon = 1.0/bigH(m)
	endif

	BIGAMa(N+1) = BIGAMa(N) - DPSI(N)*DPERPa*YNa(N)/
	1 ((EPSILON*DELN)**2)
	
C	RADIAL EXB and gradB DRIFTs
	VDRIFTa(N) = VEXBa(n)+vdrdbarg(n)+vdiar(n)
	IF(N.GT.10.AND.N.LT.41) ARGCORE = ARGCORE + 
     1	VDRIFTA(N)*DPSI(N)*YNA(N)*epsilon*DELN
	IF(N.LT.10.or.n.gt.41) 
     1	VDRIFTA(n) = -1.*(abs(VDRDBARg(N)) +abs(vexba(n)) 
     1	+abs(vdiar(n)))
	BIGAMa(N+1) = BIGAMa(N+1)
     1  +(VDRIFTa(N))*YNa(n)/(epsilon*DELN)*DPSI(n)	
	IF(N.LT.10)	BIGAMa(N+1) = BIGAMa(N+1) + SBD*SEXBa*DPSI(n)
	IF(N.GT.41)	BIGAMa(N+1) = BIGAMa(N+1) + (1.-SBD)*SEXBa*DPSI(n)

C	PARALLEL ErXB FLUX
	BIGAMa(N+1) = BIGAMa(N+1) + GAMErXBa(N)
c	parallel gradB flux
	BIGAMa(N+1) = BIGAMa(N+1) + GAMDBarg(N) 
C	source distribution
	bigama(n+1) = bigama(n+1) + sa(n)*dpsi(n) 
	Z0 = (YNa(N)/YN(N))*(ZARg**2)
 	CI1 = ((1.+0.24*Z0)*(1.+.93*Z0))/((1.+2.65*Z0)*(1.+.285*Z0))
	CI2 =1.56*(1+1.414*Z0)*(1.+0.52*Z0)/
	1	((1.+2.65*Z0)*(1.+0.285*Z0)*(Z0+SQRT(0.5*(1.+XMAS1/XMASa))))
      RES = CRES/(1.9E3*XT(N)**1.5)
	IF(N.NE.1.AND.N.NE.NCAP-1) DT =0.5*XK*(XT(N+1)-XT(N-1))
	IF(N.EQ.1) DT = 0.5*XK*(XT(2)-XT(1))
	IF(N.EQ.NCAP-1) DT = 0.5*XK*(XT(NCAP-1)-XT(NCAP-2))
C			COULOMB LOGARITHM
	Y = SQRT(YN(N))*(ATNUM(1)**2)*zarg 
	X = (EP0/EQ)**1.5	
	COULOG(2,1) = LOG(12.*3.1416*(XT(N)**1.5)*X/Y)
C			REDUCED MASS
	XMR = XMAS(1)*XMASa/(XMAS(1)+XMASa)
	XMR12 = XMAS(1)*(1.+XMAS(1)/XMASa)
	C1 = 1./((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5)) 
	XNUC(2,1)=3.34*(COULOG(2,1)*((ATNUM(1)*zarg)**2)*1.E-6*YN(N))
     2			/(C1*SQRT(XMR12*1E3)*(XT(N)**1.5))
	XNUaI = XNUC(2,1)
	XNUaI =	YN(N)*(zarg**2)*(ATNUM(1)**2)*(EQ**2.5)*COULOG(2,1)
     2		  / (47.25*SQRT(XMR)*(EP0**2)*(XT(N)**1.5))

c	IF(N.LT.NCAP-2) VEa(N+1) = 2.*BIGAMa(N+1)/(YNa(N+1)+YNa(N+2))
c	IF(N.EQ.NCAP-2) VEa(N+1) = BIGAMa(N+1)/YNa(N+1)
	IF(N.NE.1.AND.N.NE.NCAP-1) DT =0.5*XK*(XT(N+1)-XT(N-1))
	IF(N.EQ.1) DT = 0.5*XK*(XT(2)-XT(1))
	IF(N.EQ.NCAP-1) DT = 0.5*XK*(XT(NCAP-1)-XT(NCAP-2))
 
	YNa(N+1) = (XT(N)/XT(N+1))*YNa(N) + 
     2		XMASa*(VEa(N)*BIGAMa(N)-VEa(N+1)*BIGAMa(N+1))/XT(N+1)
     2		-ynz(n)*zcarb(n)*eq*(espot(n+1)-espot(n))/xt(n+1) +
	3		(YNa(N)/XT(N+1))*(Zarg**2)*(CE2/xZEFF(n)+CI2)*DT	 -
     3(YNa(N)/XT(N+1))*(ZARg**2/XZEFF(N))*EQ*RES*0.5*(XJ(N+1)+XJ(N))
     2*DPSI(N)+CI1*XMASa*XNUaI*0.5*((YNa(N)/YN(N))*(BIGAM(N)+BIGAM(N+1))
     5	 - (BIGAMa(N)+BIGAMa(N+1)))*dpsi(n)/XT(N+1)
	VEa(N+1) = 2.*BIGAMa(N+1)/(YNa(N+1)+YNa(N))	 	
 
 
620	CONTINUE

	n  = ncap
	bigama(n) = bigama(n-1) +(VDRIFTa(N-1))*YNa(n-1)/
	1	(EPD*DELN)*DPSI(n-1)
     1	 - DPSI(N)*DPERPa*YNa(N)/((EPd*DELN)**2)
	2	 + DPSI(N-1)*Sa(n-1)
	3	+  GAMErXBa(N-1) 	+ GAMDBarg(N-1) 
     4	+ (1.-SBD)*SEXBa*DPSI(n-1)
C*********************GRADB & CURVATURE DRIFTS ARGON*******************
	
	IF(IGRADB.EQ.0) GOTO 1620 
	BDIR = 1.0
	IF(IB.EQ.1) BDIR = -1.0

C	PARALLEL PLASMA CURRENT

 	DO 1600 N = 1, NCAP-1
	epsilon = epd
	 	if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
	if(n.gt.10.and.n.lt.41) epsilon = 1.0/bigH(m)
	ZLNM(N) = 1./(epsilon*DELN)
	ZLTM(N) = 1./(epsilon*DELN)
	ZLPM(N) = 1./(epsilon*DELN)
	
	
C	CHANGE IN PARALLEL CURRENT DUE TO RADIAL GRADB & CURV DRIFT CURRENTS

	vdrzarg(n) =  BDIR*3.*XK*XT(N)/(zarg*eq*abs(RMAJOR*BFIELD))
	DELJarg(N)=2.*yna(n)*zarg*eq*dpsi(n)*vdrzarg(n)*
     1	zlpm(n)*dotrz(n)
c	radial drift velocity
	vDrdbarg(n) = vdrzarg(n)*dotrz(n)
  	XJ(N+1) = XJ(N+1) + DELJarg(N)
1600	CONTINUE

C	PARALLEL PARTICLE FLOW
 	DO 1610 N = 1, 50
	gamdbarg(n) = yna(n)*vdrzarg(n)*dothetz(n)*betag
	bigama(n+1) = bigama(n+1) +gamdbarg(n)
1610	CONTINUE
1620	continue
C**************************EXB CARBON*********************************************

	if(ioptexb.eq.0) goto 635
C	EXB DRIFTDUE TO PARALLEL ELECTRIC FIELD
	DO 625 N = 1, NCAP-1
	VEXBa(N) = vexb(n)
625	CONTINUE
	SEXBa = 0.0
	IF(BDIR.EQ.0) THEN
	DO 626 N = 42, NCAP-1
	SEXBa = SEXBa + ABS(VEXBa(N)*YNa(N)/(EPD*DELN))*dpsi(n)/DIVDEPTH1
626	CONTINUE
	ENDIF
	GOTO 628
	DO 627 N = 1, 9
	SEXBa = SEXBa + ABS(VEXBa(N)*YNa(N)/(EPD*DELN))*dpsi(n)/DIVDEPTH1
627	CONTINUE
628	CONTINUE


c	a fraction fexb of the ions drifting from outer leg reach inner leg
	fexb = 1.0
	sexba = fexb*sexba
		
c	EXB DRIFT DUE TO RADIAL ELECTRIC FIELD, 
	DO 631 N = 11,40
	FON = XNOV(N)/YN(N)
	CALL CXRCEFITS(6,XT(N),fon,celz,dlzzdt,ZAV)
C	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
C	1	(FLUXHEAT/(XK*XT(N))-2.5*FLUXPART)
     	delt = deltin  
			if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n

	GAMErXBz(N) = BDIR*BETAG*(YNa(N)*ESPOT(N)/BFIELD)/(DELT/bigH(m))
C	BIGAMa(N) = BIGAMa(N) + 0.5*(GAMErXBa(N)+GAMErXBa(N-1))
631	CONTINUE
	DO 632 N= 1,10
C	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
C	1	(FLUXHEAT/(XK*XT(N))-2.5*FLUXPART)
	delt = 	deltin
	GAMErXBa(N) = -1.*BDIR*BETAG*(YNa(N)*ESPOT(N)/BFIELD)/(epd*DELT)
	GAMErXBa(NCAP-N)=-1.*BDIR*BETAG*(3.*YNa(NCAP-N)*XT(NCAP-N)/
     1	BFIELD)/(epd*DELT)
632	CONTINUE 
C	DO 633 N=2,10
c	DELT = (5.*ALB*YN(N)*XT(N)/(32.*BFIELD))/
c 	1	(FLUXHEAT/(XK*XT(N))-2.5*FLUXPART)	
C	BIGAMa(N) = BIGAMa(N) + 0.5*(GAMErXBa(N)+GAMErXBa(N-1))
C	BIGAMa(NCAP-N)=BIGAMa(NCAP-N)+0.5*(GAMErXBa(N)+GAMErXBa(NCAP-N-1))
C633	CONTINUE	 

C**************************************************************************
635	continue

c*****************diamagnetic drift argon***********************************
	if(ioptdia.eq.0) goto 641
c	used only for particle source/sink, not for current
	vdiar(1) = -0.5*bdir*xk*(yna(2)*xt(2)-yna(1)*xt(1))/
	1	(yna(1)*eq*betag*bfield)
	do 640 n = 2, ncap-2
	vdiar(n) = -0.5*bdir*xk*(yna(n+1)*xt(n+1)-yna(n)*xt(n))/
	1	(yna(n)*eq*betag*bfield)
640	continue
      vdiar(ncap-1) = -0.5*bdir*xk*(yna(ncap-1)*xt(ncap-1)-
     2	yna(ncap-2)*xt(ncap-2))/(yna(ncap-1)*eq*betag*bfield)
641	continue	
675	continue
c	*****************end of 675 do loop*********************************


C**********************DENSITY & VELOCITY*********************************
		

	DO 680 N = 1, NCAP-1
	FZARgON(N) = YNa(N)/YN(N)
680	CONTINUE	
C***************************************************************************


 	


C**************************IMPURITY ITERATION**************************************
C	SET IMPURITY SWITCH ON FOR USE OF IMPURITY DISTRIBUTIONS IN RADIATION CALCULATION
	IMPURITY = IMPURITY + 1

C	ITERATE	BACK THROUGH RADIATION HEAT BALANCE WITH CALCULATED IMPURITY DISTRIBUTIONS
	IF(IMPURITY.LT.5) GOTO 1000
C**************************ASYMMETRY ITERATION**************************************
	iconverge = iconverge + 1 

c*****isym=0 symetric divertor, =1 nonsymetric divertor 	
	IF(iconverge.eq.1.and.isym.eq.0) goto 825 
c****iflux=0 for uniform q, = 1 for Miller eq q**************** 
	if(iconverge.ge.1.and.iflux.eq.1) goto 825
			   
	

	goto 1000
825	continue 

C***********************************************************************************

c*************moments of the poloidal potential distribution**************
c	poloidal average
	xnum = 0.0
	xdenom = 0.0
	xsin = 0.0
	xcos = 0.0
	do 850 n = 11,41
	xnum = xnum + 0.5*(dpsi(n)+dpsi(n-1))*espot(n)
	xdenom = xdenom + 0.5*(dpsi(n)+dpsi(n-1))
	xsin =	xsin + sin(angthet(n))*0.5*(dpsi(n)+dpsi(n-1))*espot(n)
	xcos =	xcos + cos(angthet(n))*0.5*(dpsi(n)+dpsi(n-1))*espot(n)
850	continue
	epot0 = xnum/xdenom	  
	epotsin = xsin/(epot0*xdenom)
	epotcos = xcos/(epot0*xdenom)
c***********************************************************************

c	*****************compute potential profile from Erad_exp***********************
c	epot0 = 0.01
	epot(25) = epot0
c	epot(25) = 177.0
	
c	do 8500 n = 1,25
c	erex(n) = 0.0
c8500	continue	 
	do 851 nn = 1,24
	n=25-nn
	epot(n) = epot(n+1) + 0.5*delna*(erex(n+1)+erex(n))
851	continue	  
	do 852 n= 1,25
	xlphim(n) = erex(n)/epot(n)
c	epota(n) = epot(n)
c	epot(n) = epot(n)/xte(n)
852	continue
c***********************************************************************************

c*********************calculate poloidal velocity**********************************
	ioptpot = 2
      do 860 n = 1,25
	call poloidal(n)
860	continue
c***********************************************************************************

c****************radial flux of toroidal momentum***********************************
c	goto 914
		deltheta = 0.2094 
      do 910 n = 11, 40
	
c	 elong0 = elong
c	 triang0 = triang
c	elong0 = 1.0
c	triang0 = 0.0
c	if(n.lt.11) then
c		elong0 = elongtop
c		triang0= triangtop
c	endif
	  
	x = angthet(n)
	xx = asin(triang0) 	
	zz = x + xx*sin(x)
c	yx is R
	yx = rmajor + aminor*cos(zz)
c	ys is h_theta
	ys = aminor*sqrt((cos(zz)**2)+((elong0*sin(x))**2)) 
c	yt is num 2nd term in pi_rtheta
	yt = aminor*sin(zz)*(1.+xx*cos(x))
	yc = sqrt(((sin(zz)**2)*((1+xx*cos(x))**2)+((elong0*cos(x))**2)))
	gyro = eq*abs(bfield)/xmas(1)
	gyroz = 0.5*gyro
	eta4 = 0.5*xk*(zn(n)*zt(n)+zn(n+1)*zT(n+1))/gyro
	eta4z = 2.0*eta4
	sign = bfield/abs(bfield)
	vparl = 0.5*(vee(n)+vee(n+1))
 	veepar(n) = vparl
	vezpar(n) = 0.5*(vez(n)+vez(n+1))
	ratio = 1./sqrt(1.+betag**2)
	veephi(n) = sign*veepar(n)*ratio
	vezphi(n) = sign*vezpar(n)*ratio
c	sign = 1.0
	
	delvphi = (sign*vee(n)-sign*vee(n+1))*ratio
	delvphiz =(sign*vez(n)-sign*vez(n+1))*ratio
c	if(n.gt.32) angthet(n+1) = angthet(n+1) - 6.28 

	delvthet = vee(n)-vee(n+1)
	delvthetz = vez(n)-vez(n+1)
	dvphidthet = delvphi/deltheta
	dvphidthetz = delvphiz/deltheta
	dvthetdthet = delvthet/deltheta
	dvthetdthetz= delvthetz/deltheta
	yr = cos(zy)+dR0dr*cos(x)+xx*cos(x)*sin(x)*sin(zz)

	delrthet(n) = (yc/elong0)/yr  	
	
	
   	xpirphi(n) = (-1.*eta4/ys)*(dvphidthet + veephi(n)*yt/yx)
	xpirphiz(n) =(-1.*eta4z/ys)*(dvphidthetz+vezphi(n)*yt/yx)
	xpirthet1(n) = 0.5*xpirphi(n)
	xpirthet1z(n)= 0.5*xpirphiz(n)
	xpirthet2(n) = (-0.5*eta4/ys)*betag*dvthetdthet	
	xpirthet2z(n)= (-0.5*eta4z/ys)*betag*dvthetdthetz 
	

910	continue
	vdriftsum = 0.0 
	do 911 n = 11,40
	gyro = eq*abs(bfield)/xmas(1)
	eta4 = 0.5*xk*(zn(n)*zt(n)+zn(n+1)*zT(n+1))/gyro
	sign = bfield/abs(bfield)

	if(n.eq.40) delrthet(n+1) = delrthet(11)
	dlnr = log(delrthet(n)) -	log(delrthet(n+1))
	ys = aminor*sqrt((cos(zz)**2)+((elong0*sin(x))**2))
	veethet(n) = betag*0.5*(vee(n)+vee(n+1))
	vezthet(n) = betag*0.5*(vez(n)+vez(n+1))
	
	xpirthet2(n)=xpirthet2(n)-(0.5*eta4/ys)*veethet(n)*(dlnr/deltheta)
	xpirthet2z(n)=xpirthet2z(n)-
     1	(0.5*eta4z/ys)*vezthet(n)*(dlnr/deltheta)
	xpirthet(n) = xpirthet1(n) + xpirthet2(n)
	xpirthetz(n) = xpirthet1z(n) + xpirthet2z(n)
	vradial = vrad1(25) 
	convphi(n) = 0.5*(zn(n)+zn(n+1))*xmas(1)*veephi(n)*vradial
	convthet(n) = 0.5*(zn(n)+zn(n+1))*xmas(1)*veethet(n)*vradial
	vdriftsum = vdriftsum + vdrift(n)						   
911	continue 
c**********temporary for calculating k=1, delta=0*************
914	do 915 n = 11, 40
	if(n.ge.11.and.n.le.33) m = 34 - n
		if(n.ge.34.and.n.le.40) m = 64 - n
	  	hrmet(n) = 1./delrthet(m)
	    hthetmet(n) = delthet(m)
915	continue
920	continue
c***************************************************************
5000	format(I2,7e10.3)
	OPEN(123,FILE='DIVSOL.TXT',STATUS='UNKNOWN') 

	write(123,'(1x,35A)') 'psi D dens   C dens   Arg dens    VEL_D 
     1 VEL_C     VEL_A     temp  '
     	do 755 n = 1,51
	write (123,5000) n,yn(n),ynz(n),yna(n),VEE(n),VEZ(n),VEA(n),
     1	xt(n)
755   continue
	
	Write(123,'(1x,35A)') 'psi potential   erad   EPAR   CURRENT
     1FZ_CARB  FA_ARG '
      do 760 n = 1,51
	write (123,5000) n,espot(n),er(n),EPAR(N),xj(n),
     1	fzcarbon(n),fzargon(n)
760   continue
	Write(123,'(1x,35A)') 'psi GAM_D   GAMEXB_D  GAMGRADB_D	  BIGMOM
	1  RAD   q/<q>  '
      do 765 n = 1,51
	write (123,5000) n,BIGAM(N),GAMErXB(N),GAMDB(N),BIGMOM(N),SRAD(N),
	1	cheatx(n)
765   continue
	Write(123,'(1x,35A)') 'psi GAM_C   GAMEXB_C  GAMGRADB_C  GAM_Ar
     1GAMEXB_Ar  GAMGRADB_Ar '
      do 770 n = 1,51
	write (123,5000) n,BIGAMZ(N),GAMErXBZ(N),GAMDBCARB(N),BIGAMA(N),
	1	GAMErXBA(N),GAMDBARG(N)
770   continue
	Write(123,'(1x,35A)') 'psi VrEXB_D  VrDRDB_D  VrEXB_C  VrGRB_C
     1VrEXB_Ar  VrGR_Ar   VrdiaD  '
      do 775 n = 1,51
	write (123,5000) n,VEXB(N),VDRDB(N),VEXBz(N),VDRDBCARB(N),
	1	VEXBA(N),VDRDBARG(N),vdia(n)
775   continue
	WRITE(123,'(1X,35A)')'PSI  CURRENT    XJIN   GRADBCUR   RGBCUR 
     1PGBCUR    DELJB   bigQ  '
	currgradb = 0.0
	DO 780 N = 1,NCAP
	 currgradb = currgradb + deljb(n)
      WRITE (123,5000) N, XJ(N),XJIN,GRADBCUR(N),XJB(N),2.*EQ*GAMDB(N),
	1					DELJB(N),bigQ(N)
780	CONTINUE
	WRITE(123,'(1X,35A)')'PSI    neuts    sat     srad       T 
     1   n       bigG     bigQ  '

	 DO 790 N = 1,NCAP
	 WRITE (123,5000) N, Xnov(N),sat(n),srad(N),zT(N),zN(N),
	1					bigam(N),bigQ(N)
790	CONTINUE
		
	WRITE(123,'(1X,35A)')'PSI   theta    hr       hthet   '

	 DO 795 N = 1,NCAP
	 WRITE (123,5000) N, angthet(N), hrmet(n), hthetmet(n)
795	CONTINUE
	WRITE(123,'(1X,35A)')'PSI   vtor    pirphi    convphi 	total
	1vpar  xpirthet1  '

	 DO 800 N = 1,NCAP
	
	
	 WRITE (123,5000) N, veephi(n),xpirphi(n),convphi(n),
     1	xpirphi(n)+convphi(n),veepar(n),xpirthet1(n)
800	CONTINUE
	WRITE(123,'(1X,35A)')'PSI   theta    vperp  pirthet    convthet 
	1total   pirthet2 '

	 DO 805 N = 1,NCAP

 

	 WRITE (123,5000) N,angthet(N),veethet(n),xpirthet(n),convthet(n),
	1 xpirthet(n)+convthet(n),xpirthet2(n)	

805	CONTINUE 
	write(123,'(1x,35a)')'psi  vperpz  vparz    vezphi    piphiz  
     1pithetz   '	 
	do 810 n=1,ncap
	write (123,5000) n,vezthet(n),vezpar(n),vezphi(n),xpirphiz(n),
     1	xpirthetz(n)
810	continue
	write(123,'(1x,35a)')'psi    VEE      VEZ    press  '
	do 815 n=1,ncap
	if(n.eq.1) znz = ynz(1)
	if(n.eq.ncap) znz = ynz(ncap-1)
	if(n.gt.1.and.n.lt.ncap) znz = 0.5*(ynz(n-1)+ynz(n))
	press = (2.*zn(n)+(1.+zcarb(n))*znz)*xk*zt(n)
	write (123,5000) n, vee(n), vez(n), press
815	continue		    													  


900	FORMAT(1X,8A,E7.3)
901	format(1x,'nstagp=',i3,1x, 'nstagq=',i3,1x,'tncap=',f6.2) 
	write (123,901)nstagp,nstagq,ztncap

	OPEN(124,FILE='vpol.TXT',STATUS='UNKNOWN') 
902	format(1x,'epot0=',e10.4,1x, 'epotcos=',f8.4,1x, 'epotsin=',f8.4)
	write (124,902) epot0,epotcos,epotsin
	
903	write (124,'(1x,35a)')'nn   erex     epot      vths1      vths2   
     1  vthex' 
	do 905 n = 1, 25
	write(124,5000) n,erex(n),epot(n),vths1(n),vths2(n),vthexp(n)
905	continue
	RETURN
	END

	
