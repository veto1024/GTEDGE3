	SUBROUTINE EDGECALC
	INCLUDE 'SOLDIV.FI'
	parameter (ml=26,jq=2,jp=6)
	dimension ynold(ml),yniold(ml),
     1			velthet1old(ml),ratdrag(ml),beamdot(ml),
	2			errthet(ml),errtor(ml),cxcool(ml),coolion(ml),
	3			radcool(ml),zbar(ml),qie(ml),tiold(ml),s0v(jp,ml),
	4			delvelthet(ml),vphizterm(ml),erterm(ml),
	5			velthet2old(ml),ererad(ml),eradaold(ml),
	6			atcool(ml),gamconde(ml),gamcondi(ml),teold(ml),
	7			vtor2old(ml),pressd(ml),veld(ml),ssv(jp),bbv(jp,jp),
	8			b0v(jp,jp,ml),CHITHI(ML),CHITHE(ML)


	double precision TEDP,XLZ2DP,DLDT2DP,Tdbl,xlzdbl,dlzdbl
c	bfield = bphi
C	CALCULATES NEUTRAL & ION DENSITY, TEMP AND VELOCITY DISTRIBUTIONS IN EDGE
C	SET TOROIDAL ROTATION BC	
	
	

1591	do 3, n=1,25
	omegt(n,1) = torv(n)/rmajor
	omegt(n,4) = torv(n)/rmajor
3	continue 
	OMEGT(25,1) = OMEGTSEP1	
	OMEGT(25,4) = OMEGTSEP2
	OMEGT(25,5) = OMEGTSEP2C
	OMEGT(25,2) = OMEGTSEP1C
	OMEGT(25,3) = OMEGTSEP1S
 	OMEGT(25,6) = OMEGTSEP2S

c	**********************************
	ioptshear = 0 
c	sheareal = sheare
c	sheare = 1.0	
	ioptequil = 0
c	**********************************
	ATNUM(1) = ZION
	ATNUM(2) = ZIMP

	ynpedz = xnped
	ynbarz = xnbar
	if(ioptsoln.eq.1) then
	ynbarz = 0.5*(xnpedex+xnsepex)
	ynpedz = xnpedex
	endif
	YNI(25,1) = xnsol
	yni(25,2) = fracz*xnsol
	if(joptedped.eq.1) then
		yni(25,1) = xnsep
		yni(25,2) = fracz*xnsep
	endif
	
	TEED(25) = tsepexe
	TIED(25) = tsepexi
	ti(25) = tsepexi
	tel(25) = tsepexe
	GAMION(25,1) = enh*FLUXPART
	gamheat(25) = fluxheat
	 
	gamheate(25) = fheate*gamheat(25)
	gamheati(25) =(1.-fheate)*gamheat(25)
	R2 = 0.01
c	gamion(25,2) = r2*gamion(25,1)
	gamion(25,2) = 0.0
	xni(1) = xnsol
	xni(2) = fracz*xnsol
	tep = tsep
	temp(1) = tsep
	temp(2) = tsep
	XLNA(1)	= YLNBARX
 	XLTA(1) = YLTIBARX
	XLNA(2)	= YLNBARX
	XLTA(2) = YLTIBARX
	XLVA(1) = YLVBARX
	XLVA(2) = YLVBARX
	
	CNEUT = 1.
	IF(IOPTELN.EQ.0) CNEUT = 0.
	TED = teed(25)
	TID = tied(25)		
	IF(TED.LT.1.E-1) TED = 1.05E-1 
	IF(TED.GT.1E3) TED = .95E3
	IF(TID.LT.1.E-1) TID = 1.05E-1 
 	IF(TID.GT.1E3) TID = .95E3
     	YND = xnsol
 	IF(YND.GT.1E22) YND = 0.95E22
	IF(YND.LT.1E16) YND = 1.1E16
	TND = tied(25)
	IF(TND.GE.1000) TND = 995.
	CALL INTERP(TED,TID,TND,YND)
	SVEL = SEL(1)
	SVELN= SELN(1)
	SVCX = SCX(1)
	SVATA(25)= (SEL(1) + CNEUT*SELN(1)*YNO(25)/YNI(25,1))+SCX(1)
	SVIONA(25) = SION(1) 
	SVREC = RECOM(1)
	XNUIONI(25) = YNO(25)*SVIONA(25)
	VCOLD = SQRT(XK*TSPL/XMASS)
	ynocold(25) = 0.5*(gamouteff/vcold + coldno(25))
      XNUATI(25) =  YNOCOLD(25)*SVATA(25)

	xni(1) = xnsol
	xni(2) = fracz*xnsol
C	***************TEMPORARY*************************
 
	DELTAN = (XNPEDEX-XNSOL)/24.
	DELTATE = (TPEDEXE-TSEPEXE)/24.
	DELTATI = (TPEDEXI-TSEPEXI)/24.
c	***********************
	rhor(25) = 1.0
	if(ioptran.eq.1) then
	chiion(25) = chixpi
	chiel(25) = chixpe
	endif
 	DO 5, NN=1,24
	N = 25-NN
	rhor(n) = rhor(n+1) - delna/(aminor*SQRT(0.5*(1.+ELONG**2)))
	if(ioptran.eq.1) then
	chiion(n)=chixpi
	if(rhor(n).lt.pedrhoti) chiion(n)=chitop 
	chiel(n)=chixpe
	if(rhor(n).lt.pedrhote) chiel(n)=chetop
	endif
	sheare(n) = 1.0
	shearfac(n) = 0.0
 	xlnm(n) = 1./ylnbarx
	xlnim(n)= 1./ylnbarx
	xlnzm(n)= 1./ylnbarx
 	xltim(n)= 1./yltibarx 
	xltem(n)= 1./yltebarx 
	xlvm(n) = -2.*(torv(n)-torv(n-1))/(delna*(torv(n)+torv(n-1)))
	if(ioptxlvm.eq.1) xlvm(n) = 1./xlv1
	xlvm2(n) = xlvm(n)
	vtor1old(n) = torv(n)
	YNI(N,1) = YNI(N+1,1) + DELTAN
	yni(n,2) = fracz*yni(n,1)
	TEED(N) = TEED(N+1) + DELTATE
	TIED(N) = TIED(N+1) + DELTATI
	tel(n) = teed(n)
	ti(n) = tied(n)
5	CONTINUE
	
	sheare(25) = 1.0
	shearfac(25) = 0.0
     	xlnm(25) = 1./ylnbarx 
	xlnim(25)= 1./ylnbarx
	xlnzm(25)= 1./ylnbarx
 	xltim(25)= 1./yltibarx 
	xltem(25)= 1./yltebarx
	xlvm(25) = -2.*(torv(25)-torv(24))/(delna*(torv(25)+torv(24)))
	if(ioptxlvm.eq.1) xlvm(25) = 1./xlv1 
	xlvm(1) = xlvm(2)
	xlvm2(25) = xlvm(25)
	xlvm2(1) = xlvm(1)
	vtor1old(25) = torv(25)
C	*****************************************
	
c	"Outer" iteration converging velocity calculations & density calculations 
c		through 300
c	if(ioptapproach.eq.5) then
c	first converges the solution using vphi1=vphi2= exp vphi2, vthet2=vthetexp 
c	iconverge = 1
c	ioptapproach = 1
c	endif 
6	kk=0
7	do 300 jt = 1,101
	mjt = 0 

	CALL NEUTDIST
4999	continue 
c	**********control of edge calculation***********************************
c	ioptapproach = 1: 1) fit nudrags to vphi1=vphi2= vphi2exp using vth1=vth2=vth2exp, and then 
c	calculates vth1 and vth2 to cal gyro and inertial drags, cal atomic drag and construct
c	anom drag as dif between fit drag and (gyro+atomic+inert drags). 2) cal vphi1 & vphi2 as check
c	and cal vth1 & vth2 to recal gyro and inert drags, then recal anom drag.	
	if(ioptapproach.eq.1) goto 5000
c	ioptapproach = 2: 1) first does part 1 of ioptapproach=1 cal, then 2) repeats part 1 fitting 
c	drags to vphi2=vphi2exp and vphi1=vphi2exp+difcal using vth2=vth2exp and vth1=vth2exp+difcal.
c	3) step 2 above
	if(ioptapproach.eq.2) then
	ioptap = 2
	ioptapproach = 1
	endif
	if(ioptapproach.eq.5) then
	ioptap = 5
	ioptapproach = 1
	endif
	ntap = 0	 
c	************end control section****************************************
	
c	*****"Inner" iteration converging ion and neutral distributions through 200******
 
5000 	DO 200 IT = 1,100
	XNUIONI(25) = YNO(25)*SVIONA(25)
      XNUATI(25) =  YNOCOLD(25)*SVATA(25)
	XNUDRAGATOMIC(25) = 0.5*(XNUIONI(25)+xnuioni(24)) +	
     1	coldno(25)*svata(25)
     
c	update velocity gradient scale length
	kk = kk+1
	goto 9
	if(ioptxlvm.eq.1) goto 9
	if(ioptxlvm.eq.3) goto 9

 	if(it.gt.1) then
	do 8, n=2,25
	if(ioptxlvm.eq.2) xlvm(n) = -2.*(vphiex1(n)-vphiex1(n-1))/
	1		(delna*(vphiex1(n)+vphiex1(n-1)))
c	if(ioptxlvm.eq.3) xlvm(n) = -2.*(vphiex2(n)-vphiex2(n-1))/
c	1		(delna*(vphiex2(n)+vphiex2(n-1)))
	if(ioptxlvm.eq.4) xlvm(n) = -2.*(vtor1(n)-vtor1(n-1))/
	1		(delna*(vtor1(n)+vtor1(n-1)))
	if(ioptxlvm.eq.5) xlvm(n) = -2.*(omegt(n,1)-omegt(n-1,1))/
     1		(delna*(omegt(n,1)+omegt(n-1,1)))

8	continue
	xlvm(1) = xlvm(2)
	endif 

9	continue
	do 10 n = 1,25
	velthet1old(n) = velthet1(n)
	velthet2old(n) = velthet2(n)
	eradaold(n) = erada(n)
	
10	continue
 	 
	do 25 n = 1,25
	ynold(n) = yno(n)
	yniold(n) = yni(n,1)
25	continue
c	impurity radiation & average charge
	TDBL = 0.5*(Tel(24)+tel(25))
 	eTDBL = 0.5*(Tel(24)+tel(25))

	IZ1 = izintrin
	IZ2 = atnum(1)
	if(iz1.eq.4.or.iz1.eq.6.or.iz1.eq.74) then
	fon = yno(25)/(0.5*(yni(24,1)+yni(25,1)))
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,eTdbl,fon,exlzdbl,edlzdbl,ZAV)
	zbar2(25) = zav
c*******correct input separatrix density to electron density*******
c	yni(25,1)=yni(25,1)/(1.+fracz*zbar2(25))
c	yni(25,2)=fracz*yni(25,1)
	zne(25) = yni(25,1) + zbar2(25)*yni(25,2)
	goto 28 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
28	XLradZ(25)= eXLZDBL
	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(25) = XLradZ(25)*1.e-13
	xni(1) = yni(25,1)
	xni(2) = yni(25,2)
	vrad1(25) = gamion(25,1)/yni(25,1)
 
	temp(1) = ti(25)
	temp(2) = ti(25)
	tep     = tel(25)
	call param(25)
	xnuc12(25) = xnuc(1,2)
	xnuc21(25) = xnuc(2,1)
	delma = delna
c	*********************	
      cncmult = 1.0 
	cetaimult = 0.0
	cetgmult = 1.0
	cedwmult = 0.0
	radmultedge = 1.0
c******************************



c*****note******** 
	ntorque = 0

	call edgerotran(25,ntorque)
    	thetw1(25) = thetw(1)
	thetw2(25) = thetw(2) 
	ratdrag(25) = xnudragatomic(25)/ynudrag1(25)

      

	
	if(chiion(25).le.0.0) chiion(25) = 1.0
	if(chiel(25).le.0.0) chiel(25) = 1.0
 
c	inverse temperature gradient scale lengths	

	xltim(25) = (gamheati(25)/(yni(25,1)*xk*ti(25))-2.5*vrad1(25))/
     1 	(chiion(25))
c	xltim(25) = 1./ylti
	yne = yni(25,1)*atnum(1)+yni(25,2)*zbar2(25)
	vrade=(atnum(1)*yni(25,1)*vrad1(25)+zbar2(25)*yni(25,2)*vrad2(25))
     1	/yne
	xltem(25) = (gamheate(25)/(yne*xk*tel(25))-2.5*vrade)
     1		/(chiel(25))
c	xltem(25) = 1./ylte 
c	***************************infer chi from exp***********************
c	if(ioptran.eq.1) then
c		xltim(25) = 1./yltibarx
c		xltem(25) = 1./yltebarx 
c	endif
c	********************************************************************
	xltzm(25) = xltim(25) 
c	inverse pressure and ion density gradient scale lengths
c	if(it.gt.1) diff(25) = 0.5*(diff(24)+diff(25)) 
	xlpm(25) = (vrad1(25)-vpinchi(25))/diff(25)
	if(xlpm(25).lt.0.0) xlpm(25) = 0.0
c	if(it.gt.1) xlpm(25) = xlpm(24)
 	c10=1.0
c	if(kk.eq.3) c10 = 0.1 
c	if(kk.eq.4) c10 = 0.2 
c	if(kk.eq.5) c10 = 0.3 
c	if(kk.eq.6) c10 = 0.5 
c	if(kk.eq.7) c10 = 0.7 
c	if(kk.ge.8) c10 = 1.0

	xlnmold = xlnm(25)
c	xlnm(25) = xlpm(25)-c10*xltim(25)
c	xlnm(25) = xlpm(25)-c10/ylti  
c	if(kk.gt.3) xlnm(25) = 0.5*(xlnm(25)+xlnmold)
	if(xlnm(25).lt.0.0) xlnm(25) = 0.0			
c	**********skip***********************************
	goto 35
	xnumi = dzz(25)*(vrad1(25)-(vpinchi(25)+(dii(25)-diz(25))*
     1	xltim(25))) +
     1  diz(25)*(vrad2(25)-(vpinchz(25)+(dzz(25)-dzi(25))*xltzm(25)))
	xnumz = dii(25)*(vrad2(25)-(vpinchz(25)+(dzz(25)-dzz(25))*
     1	xltim(25))) +
     1   dzi(25)*(vrad1(25)-(vpinchi(25)+(dii(25)-diz(25))*xltzm(25)))
      denom = dii(25)*dzz(25) - diz(25)*dzi(25)			
	xlnim(25) = xnumi/denom
	xlnzm(25) = xnumz/denom
c	****************************************************
35	continue	

	DO 100 NN = 1,24
c	if(it.gt.1) yno(n)=0.5*(yno(n)+ynold(n))
	cedwmult = 1.0		
	MIT = 0	
	N = 25 - NN
	j=n
	delma = delna

	XNI(1) = YNI(N,1)
c	YNIOLD(n) = YNI(N,1)
	XNI(2) = YNI(N,2)
	TEMP(1) = TI(N)
	TEMP(2) = TI(N)
	TEP = TEl(N)
	
	XNUIONI(N) = YNO(N)*SVIONA(N)
	ynocold(n) = 0.5*(coldno(n)+coldno(n+1))
	XNUATI(N) =  YNOCOLD(N)*SVATA(N)
	
	CALL PARAM(N)
	xnuc12(n) = xnuc(1,2)
	xnuc21(n) = xnuc(2,1)
     
c	integrate particle and heat fluxes inward from separatrix
c	***********skip*******************
	goto 45
c	old density pedestal formulation
	thetint = 0.5
	DENS(N) = thetint*YNI(N+1,1) + (1.-thetint)*YNI(N,1)
	
	GAMION(N,1) = GAMION(N+1,1) - dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n)) 
45	continue
c	********************************
c	new ion & plasma density and ion & electron temp formulation
	dens(n) = 0.5*(yni(n+1,1)+yni(n,1) )
	tele = 0.5*(tel(n)+tel(n+1))
	tiav = 0.5*(ti(n)+ti(n+1))
	rz = (atnum(1)**2/xmas(1) + zbar2(n)**2/xmas(2))*1.67e-27  
	cequil = 6.32e-14*rz/(tele**1.5) 
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)

	qie(n) = cequil*yne*(tiav-tele)
	cxcool(n)=1.5*dens(n)*tiav*xk*xnuati(n)
c	if(n.eq.24) cxcool(n) = 1.5*dens(n)*tiav*xk*xnuati(n)
c**********************
	cmulteq = 0.0

c	if(kk.eq.2) cmulteq = 0.1
c	if(kk.eq.3) cmulteq = 0.2
c	if(kk.eq.4) cmulteq = 0.5
c	if(kk.ge.5) cmulteq = 1.0
c**********************
c	if(ioptequil.eq.1) cmulteq = 0.5
c	if(ioptequil.eq.2) cmulteq = 1.0
c*******************************************
c	delma = radmultedge*delna
c*******************************************
	
	
c	electron heat flux
	EIONi = 17.5
	IF(dens(n).LE.1.E21) 
     2     EIONi = 17.5 + (5.+37.5/Tel(n))*LOG10(1.E21/dens(n))
    	IF(dens(n).GT.1.E21)
     2    EIONi = (30.6 - 16.4*EXP(-5.E19/dens(n)))*
     3              EXP(5.45/(Tel(n)*EXP((dens(n)/1.37E20)**0.26)))
	TDBL = 0.5*(Tel(n)+tel(n+1))
	eTDBL = 0.5*(Tel(n)+tel(n+1))
	IZ1 = izintrin
	IZ2 = atnum(1)
	if(iz1.eq.4.or.iz1.eq.6.or.iz1.eq.74) then
	fon = yno(n)/(0.5*(yni(n,1)+yni(n+1,1)))
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,eTdbl,fon,exlzdbl,edlzdbl,ZAV)
	zbar2(n) = zav
	goto 55 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
55	XLradZ(n)= eXLZDBL
	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(n) = XLradZ(n)*1.e-13
	coolion(n) = xk*eioni*yne*xnuioni(n)
	radcool(n) = 0.5*(yni(n,2)+yni(n+1,2))*yne*xlradz(n)
	radcool(n) = radmultedge*radcool(n)

	ntorque = 0
	CALL EDGEROTRAN(N,ntorque)
	thetw1(n) = thetw(1)
	thetw2(n) = thetw(2)
c	particle fluxes on ions & impurities
 	thetint = 0.5
c	DENS(N) = thetint*YNI(N+1,1) + (1.-thetint)*YNI(N,1)
	GAMION(N,1) = GAMION(N+1,1) - dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n)) 
	GAMION(N,2) = GAMION(N+1,2) 
     1		-	DENS(N)*XNUIONz(N)*DELMA*(1.+fracz*zbar2(n))
	

   	beamdot(n) =  dens(n)*xnuionb(n)
	
	gamheati(n) = gamheati(n+1) + delma*(cxcool(n) + cmulteq*qie(n))
 	gamheati(n) = gamheati(n) - fionb(n)*qnb(n)*delma	 
	gamheate(n)=gamheate(n+1)+delma*
     1					(coolion(n)+radcool(n)-cmulteq*qie(n))
	gamheate(n) = gamheate(n) -	(1.-fionb(n))*qnb(n)*delma
c	time dependent
	gamheati(n) = gamheati(n)-dlnw_dt*1.5*yni(n,1)*xk*ti(n)*delma
	gamheate(n) = gamheate(n)-dlnw_dt*1.5*yne*xk*tel(n)*delma
	gamion(n,1) = gamion(n,1)-dln_dt*yni(n,1)*delma
     	
c*******************

	delma = delna 
c*******************
c	assume fracheate, instead of calculating equilibration
	gamheat(n) = gamheate(n) + gamheati(n)
c	gamheate(n) = fheate*gamheat(n)
c	gamheati(n) = gamheat(n)- gamheate(n)
c	solve poloidal velocities & calculate coefficient for advancing density inward

60	continue
c	new ion & impurity density, ion & electron temp formulation
c	calculate inverse temperature gradient scale lengths
	if(chiion(n).le.0.0) chiion(n) = 1.0
	if(chiel(n).le.0.0) chiel(n) = 1.0
c	use chi chang-hinton os + chi eta-i 
	xltim(n)=(gamheati(n)/(yni(n,1)*xk*ti(n))-2.5*vrad1(n))/
	1			(chiion(n))
c	if(rhor(n).gt.pedrhoti) xltim(n) = 1./ylti 
c	if(rhor(n).le.pedrhoti) xltim(n) = 1./xltitop 

c	if(xltim(n).lt.0.0) xltim(n) = 0.0
	

	xltzm(n) = xltim(n) 
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	vrade=(atnum(1)*yni(n,1)*vrad1(n)+zbar2(n)*yni(n,2)*vrad2(n))/yne
	xltem(n) = (gamheate(n)/(yne*xk*tel(n))-2.5*vrade)
	1	/(chiel(n))
c	if(rhor(n).gt.pedrhote) xltem(n) = 1./ylte
c	if(rhor(n).le.pedrhote) xltem(n) = 1./xltetop 
c	if(xltem(n).lt.0.0) xltem(n) = 0.0
c	***************************infer chi from exp***********************
c	if(ioptran.eq.1) then
c		if(rhor(n).gt.pedrhoti) xltim(n) = 1./yltibarx
c		if(rhor(n).gt.pedrhote) xltem(n) = 1./yltebarx 
c	endif
c	********************************************************************



c	**************skip****************************** 	
	goto 65
c	calculate ion & impurity inverse density gradient scale lengths 
     	xnumi = dzz(n)*(vrad1(n)-(vpinchi(n)+(dii(n)-diz(n))*xltim(n))) +
     1	    diz(n)*(vrad2(n)-(vpinchz(n)+(dzz(n)-dzi(n))*xltzm(n)))
	xnumz = dii(n)*(vrad2(n)-(vpinchz(n)+(dzz(n)-dzz(n))*xltim(n))) +
     1	    dzi(n)*(vrad1(n)-(vpinchi(n)+(dii(n)-diz(n))*xltzm(n)))
      denom = dii(n)*dzz(n) - diz(n)*dzi(n)			
	xlnim(n) = xnumi/denom
	xlnzm(n) = xnumz/denom
c	advance ion & impurity densities
	c6 = 0.0	
	if(n.lt.25) yni(n,1) = yni(n+1,1)*
	1	(1.+delma*0.5*(c6*xlnim(n)+(1.-c6)*xlnim(n+1)))/
     2	(1.-delma*0.5*(c6*xlnim(n)+(1.-c6)*xlnim(n+1)))
	if(n.eq.25) yni(n,1)=yni(n+1,1)/(1.-delma*xlnim(n))
	if(yni(n,1).le.0.0) yni(n,1) = 1.e19 
	if(n.lt.25) yni(n,2) = yni(n+1,2)*
 	1	(1.+delma*0.5*(c6*xlnzm(n)+(1.-c6)*xlnzm(n+1)))/
     2	(1.-delma*0.5*(c6*xlnzm(n)+(1.-c6)*xlnzm(n+1)))
	if(n.eq.25) yni(n,2)=yni(n+1,2)/(1.-delma*xlnzm(n))
	if(yni(n,2).le.0.0) yni(n,2) = 1.e17
c	88888888888888
	yni(n,2) = fracz*yni(n,1)
 	zne(n) = yni(n,1) + zbar2(n)*yni(n,2)
	frazimp(n) = yni(n,2)/yni(n,1)
c	******************************************************
c	inverse pressure and ion density gradient scale lengths
c	old fz = const density pedestal formulation 
65	continue 

	xlpm(n) = (vrad1(n)-vpinchi(n))/diff(n)	
	if(xlpm(n).lt.0.0) xlpm(n) = xltim(n)
	
	xlnmold = xlnm(n) 
	xlnm(n) = xlpm(n)-c10*xltim(n) 
c	if(rhor(n).gt.pedrhoti) xlnm(n) = xlpm(n)-c10/ylti 
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpm(n)-c10/xltim(n) 
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpm(n)-c10/xltitop  
c	if(kk.gt.3) xlnm(n) = 0.5*(xlnm(n)+xlnmold)
	if(xlnm(n).lt.0.0) xlnm(n) = 0.0
c	integrate ion density inward from separatrix
	c6 = 0.5		
	yni(n,1) = yni(n+1,1)*
	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
     2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
c	yni(n,1)=yni(n+1,1)*exp(delna*(c6*xlnm(n)+(1.-c6)*xlnm(n+1))) 
c	if(n.eq.24) yni(n,1)=yni(n+1,1)/(1.-delma*xlnm(n))
	if(yni(n,1).le.0.0) yni(n,1) = 1.e19 
	if(yni(n,1).gt.2.e20) yni(n,1) = 2.e20
	yni(n,2) = fracz*yni(n,1)
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	zne(n) = yne 
	vrad1(n) = gamion(n,1)/yni(n,1)
	vrad2(n) = gamion(n,2)/yni(n,2)
c	************edit terms***************************	
c	ratdrag(n) = xnudragatomic(n)/ynudrag1(n)
c	delvelthet(n) = xmas(1)*(xnuc12(n)+ynudrag1(n))*
c	1	(velthet1(n)-velthet2(n))/(fp*eq*bthet)
c	vphizterm(n) =xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
c	erterm(n) =xmas(1)*(xnuc12(n)+ynudrag1(n))*
c	1	((erada(n)/bthet)+velthet1(n)/fp)/(eq*bthet)
c	***************************************************
 
c	integrate ion & electron temperatures inward from separatrix
	tiold(n) = ti(n)
	teold(n) = tel(n)
	ti(n) = ti(n+1)*
     1	(1.+delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))/
     2	(1.-delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))
c	if(n.eq.24) ti(n) = ti(n+1)/(1.-delma*xltim(n))	
	
      if(ti(n).lt.0.0) ti(n) = 100.
c	if(ti(n).gt.1.e3) ti(n) = 1.e3
	tel(n) = tel(n+1)*
     1	(1.+delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))/
     2	(1.-delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))
c	if(n.eq.24) te(n) = te(n+1)/(1.-delma*xltem(n))	
	scv = tel(n)	
      if(tel(n).lt.0.0) tel(n) = 100.
c	if(tel(n).gt.1e3) tel(n) = 1.e3
C	INTEGRATE TOROIDAL ROTATION FREQ INWARD FROM SEPARATRIX

c	add anomalous convective velocity to calculated radial velocity
c	vanom = 0.0 
c	if(j.lt.23) vanom =0.0*(24-j)
c	vrad1(j) = vrad1(j)+ vanom
c	vrad2(j) = vrad2(j)+ vanom 
	goto 1104
	if(n.eq.24) call torotate(25)
	do 1045 mm = 1,6 
	s0v(mm,25) = sv(mm,25)
1045  continue
	CALL TOROTATE(N)
	
      DO 1055 mm = 1,6
	s0v(mm,j) = sV(mm,j)
c	SV(mm,j) = SV(mm,j)+SV(mm,j+1)
	DO 1050 kk = 1,6
	b0v(mm,kk,j) = bv(mm,kk,j) 
c	xsv = (AV(mm,kk,j)+AV(mm,kk,j+1))/delma - BV(mm,kk,j+1)
	
c	SV(mm,j) = SV(mm,j) + xsv*OMEGT(j+1,kk)
	SV(mm,j) = SV(mm,j) + (AV(mm,kk,j)/delma)*OMEGT(j+1,kk)

c	BV(mm,kk,j) = BV(mm,kk,j) +(AV(mm,kk,j+1)+AV(mm,kk,j))/DELMA
	BV(mm,kk,j) = BV(mm,kk,j) +(AV(mm,kk,j)/DELMA)
 

	BBV(mm,kk) = BV(mm,kk,j)
c	write(6,2111) mm,kk,bbv(mm,kk) 
1050	CONTINUE 
c	remove anomalous convective velocity from radial velocity
	vrad1(j) = vrad1(j)- vanom
	vrad2(j) = vrad2(j)- vanom 
	
	SSV(mm) = SV(mm,j)
	
1055  CONTINUE
2111	format (I2,I2,6e10.3)

	CALL LSLRG(6,BBV,6,SSV,1,SSV)
	DO 1060 kk = 1,6
	OMEGT(j,kk) = SSV(kk)
1060	CONTINUE	
1104	continue
	x = vtor1(n)
	y = vtor2(n)  
 
100	CONTINUE

		
	call neutdist
	
c	converge on main ion density and temp
	do 150 n = 1,24	
 	IF(ABS(YNIOLD(n)/YNI(N,1)-1.).GT.0.02) MIT = MIT + 1
	IF(ABS(tIOLD(n)/tI(n)-1.).GT.0.02) MIT = MIT + 1
 	if(abs(teold(n)/tel(n)-1.).gt.02) mit = mit+1 

	IF(YNI(N,1).LE.0.0) YNI(N,1) = 1.E19
	IF(YNI(N,2).LE.0.0) YNI(N,2) = 1.E19*fracz

c	if(it.gt.1) yni(n,1) = 0.5*(yni(n,1)+yniold(n))
	
150	continue 
	iF(MIT.EQ.0) GOTO 225

	nj = 50
	
	if(iconverge.eq.0.and.it.eq.nj) then
	write (6,199) nj
	jwarn = 1
	goto 225
	endif
199	format(1x,'not converged on 200 loop after iterations=',I3.0) 
200	CONTINUE
c	************approach control check*************************
225	if(ioptap.eq.2) then
	ioptapproach = 2
	ioptap = 0
	goto 5000
	endif 
	if(ioptap.eq.5) then
	if(ntap.eq.0) then
	ioptapproach = 2
	ntap = ntap + 1
	goto 5000
	endif
	ioptapproach = 5
	ioptap = 0
	goto 5000
	endif
250	continue	
c	converge on poloidal & toroidal velocities ***************turned off*******
c	ioptshear = 1 
c	if(n.gt.15) sheare = sheareal
c************convergence on velocities & Er turned off--trouble converging 3/18/04

c	goto 374
c*********************************************************************************
c	do 275 n = 1,25
c	errthet(n) = abs((velthet1old(n)/velthet1(n))-1.)
c	if(errthet(n).gt.0.05) mjt = mjt + 1
c	errthet(n) = abs((velthet2old(n)/velthet2(n))-1.)
c 	if(errthet(n).gt.0.05) mjt = mjt + 1 
c	ererad(n) = abs(eradaold(n)/erada(n))
c	if(ererad(n).gt.0.05) mjt = mjt + 1 
c	velthet1(n) = 0.5*(velthet1(n)+velthet1old(n))
c	velthet2(n) = 0.5*(velthet2(n)+velthet2old(n))
	
c	errtor(n) = abs((vtor1old(n)/vphiex1(n))-1.)
c	if(errtor(n).gt.0.05) mjt = mjt + 1
	
275	continue

	if(mjt.ne.0) goto 298
	 
	goto 350

298	do 299 n=1,25	
      vtor1old(n) = vphiex1(n)
299	continue
			
300	continue		 
	
350	continue
c**********************************
	goto 700
c*********************************
c	introduce shear into density profile calculation******turned off*******
c  note 3/5/04  much better results to te & ti w/o doing this, when weak equil used.
374	kl = 0
375	continue
	do 400 n = 1,25
c	diff(n) = diff(n)/(sheare(n)**1.5)
	xlpm(n) = (vrad1(n)-vpinchi(n))/(diff(n)) 
c	if(kl.eq.3) c10 = 0.1 
c	if(kl.eq.4) c10 = 0.2 
c	if(kl.eq.5) c10 = 0.3 
c	if(kl.eq.6) c10 = 0.5 
c	if(kl.eq.7) c10 = 0.7 
c	if(kl.ge.8) c10 = 1.0
c	if(xlpm(n).lt.0.0) xlpm(n) = 0.0
	if(jjcon.eq.1) c10=1.0 
	xlnm(n) = xlpm(n)-c10*xltim(n)
c	if(rhor(n).gt.pedrhoti) xlnm(n) = xlpm(n)-c10/ylti  
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpm(n)-c10/xltim(n)  
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpm(n)-c10/xltitop 
	if(xlnm(n).le.0.0) xlnm(n) = 0.0
400	continue
	kl = kl + 1
	if(kl.gt.50) goto 700	
c	integrate ion density inward from separatrix
	c6 = 0.5
	do 450 nmesh = 1,24
	n = 25-nmesh
	yniold(n) = yni(n,1)		
	yni(n,1) = yni(n+1,1)*
	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
     2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
c	yni(n,1)=yni(n+1,1)*exp(delna*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
c	yni(n,1) = 0.5*(yni(n,1)+yniold(n))
c	if(n.eq.24) yni(n,1)=yni(n+1,1)/(1.-delma*xlnm(n))
	if(yni(n,1).le.0.0) yni(n,1) = 1.e19 
	if(yni(n,1).gt.1.e20) yni(n,1) = 1.e20
450	continue
	
	do 451 n = 1,24
	yni(n,2) = fracz*yni(n,1)
	zne(n) = yni(n,1) + zbar2(n)*yni(n,2)
	yne = zne(n)
451	continue		 
c	recalculate neutral distribution
	call neutdist
	thetint = 0.5

	temp(1) = ti(25)
	temp(2) = ti(25)
	tep = tel(25) 
	xni(1) = yni(25,1)
	xni(2) = yni(25,2)
	call param(25)
	xnuc12(25) = xnuc(1,2)
	c12 = 1.0
	 if(ynudrag1(25).le.0.0) ynudrag1(25) = 0.0
 
	diff(25) =  xmas(1)*xk*ti(25)*xnuc12(25)*
     1			((c12*ynudrag1(25)/xnuc12(25))+1.-atnum(1)/zbar2(25))/
     2			((eq*atnum(1)*bthet)**2)
	vrad1(25) = gamion(25,1)/yni(25,1)
 
	nm = 25
	ntorque = 0
	call edgerotran (nm,ntorque)
	
c	goto 434
c	vpinchi(25) = (-1.*xmtor(1)/yni(25,1)  + 
c     1	 xmas(1)*ynudrag1(25)*((erada(25)/bthet)+velthet1(25)/fp) +
c     2	 xmas(1)*xnuc12(25)*(velthet1(25)-velthet2(25))/fp)/
c     3	 (eq*atnum(1)*bthet)
c	if(ioptpinchi.eq.5) then
c	vpinchi(25) = (-1.*xmtor(1)/yni(25,1)  + 
c     1	xmas(1)*(ynudrag1(25)+xnuc12(25))*
c     1	((erada(25)/bthet)+velthet1(25)/fp) -
c     2	 xmas(1)*xnuc12(25)*torv(25))/
c     3	 (eq*atnum(1)*bthet)
c	 
c	diff(n) =  xmas(1)*xk*ti(25)*xnuc12(25)*
c	1			((ynudrag1(25)/xnuc12(25))+1.)/
c     2			((eq*atnum(1)*bthet)**2)
c	endif
c434	continue 
c	recalculate outward ion flux and heat fluxes
	do 675 nn = 1,25
	n = 26-nn
	j=n
	delma = delna

	temp(1) = ti(n)
	temp(2) = ti(n)
	tep = tel(n) 
     	xni(1) = yni(n,1)
	xni(2) = yni(n,2)
	call param(n)
	xnuc12(n) = xnuc(1,2)
	thetint = 0.5

	DENS(N) = thetint*YNI(N+1,1) + (1.-thetint)*YNI(N,1)
	XNUIONI(N) = YNO(N)*SVIONA(N)
	ynocold(n) = 0.5*(coldno(n)+coldno(n+1))
	XNUATI(N) =  YNOCOLD(N)*SVATA(N)
	GAMION(N,1) = GAMION(N+1,1) - dens(n)*xnuionb(n)*delna	
     1			-	DENS(N)*XNUIONI(N)*DELNA*(1.+fracz*zbar2(n))


	ntorque = 0
	call edgerotran(n,ntorque)
c	original formulation eliminating vphi in all terms 
c	goto 433
c	if(ioptpinchi.eq.2) then
c	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1)  + 
c     1	 xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp) +
c    2	 xmas(1)*xnuc12(n)*(velthet1(n)-velthet2(n))/fp)/
c     3	 (eq*atnum(1)*bthet)
c		 
c	diff(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
c
c	1			((ynudrag1(n)/xnuc12(n))+1.-atnum(1)/zbar2(n))/
c     2			((eq*atnum(1)*bthet)**2)
c	endif
c	use measured vphi-z
c	if(ioptpinchi.eq.5) then
c	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1) - eq*atnum(1)*ephia + 
c     1 xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)+velthet1(n)/fp) 
c     2	 - xmas(1)*xnuc12(n)*torv(n))/
c     3	 (eq*atnum(1)*bthet)
c	diff(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
c     1			((ynudrag1(n)/xnuc12(n))+1.)/
c     2			((eq*atnum(1)*bthet)**2)
 
c	endif
c433	continue
	     
      vrad1(n) = gamion(n,1)/yni(n,1)
c	goto 675
	tiav = 0.5*(ti(n)+ti(n+1)) 
	tele = 0.5*(tel(n)+tel(n+1))
	cxcool(n)=1.5*dens(n)*tiav*xnuati(n)*xk
	

	EIONi = 17.5
	IF(dens(n).LE.1.E21) 
     2     EIONi = 17.5 + (5.+37.5/Tele)*LOG10(1.E21/dens(n))
    	IF(dens(n).GT.1.E21)
     2    EIONi = (30.6 - 16.4*EXP(-5.E19/dens(n)))*
     3              EXP(5.45/(Tele*EXP((dens(n)/1.37E20)**0.26)))
	TDBL = 0.5*(Tel(n)+tel(n+1))
	eTDBL = 0.5*(Tel(n)+tel(n+1))
	IZ1 = izintrin


	IZ2 = atnum(1)
	if(iz1.eq.4.or.iz1.eq.6.or.iz1.eq.74) then
	fon = yno(n)/(0.5*(yni(n,1)+yni(n+1,1)))
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,eTdbl,fon,exlzdbl,edlzdbl,ZAV)
	zbar2(n) = zav
	goto 455 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
455	XLradZ(n)= eXLZDBL
	 
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(n) = XLradZ(n)*1.e-13
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)

	coolion(n) = xk*eioni*yne*xnuioni(n)
	radcool(n) = 0.5*(yni(n,2)+yni(n+1,2))*yne*xlradz(n)
	radcool(n) = radmultedge*radcool(n)
	atcool(n) = cxcool(n)+coolion(n)

c	equilibrated heat flux calculation
	
	rz = (atnum(1)**2/xmas(1) + zbar2(n)**2/xmas(2))*1.67e-27  
	cequil = 6.32e-14*rz/(tele**1.5) 
	
	qie(n) = cequil*yne*(tiav-tele)
c**********************
c	cmulteq = 0.0

c	if(kl.eq.2) cmulteq = 0.1
c	if(kl.eq.3) cmulteq = 0.2
c	if(kl.eq.4) cmulteq = 0.5
c	if(kl.ge.5) cmulteq = 1.0

c*****************************
c********************************
	beamdot(n) =  dens(n)*xnuionb(n)
	
	gamheati(n) = gamheati(n+1) + delma*(cxcool(n) + cmulteq*qie(n))
	gamheati(n) = gamheati(n) - fionb(n)*qnb(n)*delma	
	gamheate(n)=gamheate(n+1)+delma*
  	1					(coolion(n)+radcool(n)-cmulteq*qie(n))
	gamheate(n) = gamheate(n) - (1.-fionb(n))*qnb(n)*delma
c	time dependent
 	gamheati(n) = gamheati(n)-dlnw_dt*1.5*yni(n,1)*xk*ti(n)*delma
	gamheate(n) = gamheate(n)-dlnw_dt*1.5*yne*xk*tel(n)*delma
	gamion(n,1) = gamion(n,1)-dln_dt*yni(n,1)*delma

	delma = delna
c**************************************************************
c	assume fracheate, instead of calculating equilibration
	gamheat(n) = gamheate(n) + gamheati(n)
c	gamheate(n) = fheate*gamheat(n)
c	gamheati(n) = gamheat(n)- gamheate(n)

c	calculate inverse temperature gradient scale lengths
	 
	xltim(n)=(gamheati(n)/(yni(n,1)*xk*ti(n))-2.5*vrad1(n))/
	1			(chiion(n))
c	if(rhor(n).gt.pedrhoti) xltim(n) = 1./ylti
c	if(rhor(n).le.pedrhoti) xltim(n) = 1./xltitop 
	gamcondi(n)= (gamheati(n)/(yni(n,1)*xk*ti(n))-2.5*vrad1(n))
	xltzm(n) = xltim(n) 
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	vrade=(atnum(1)*yni(n,1)*vrad1(n)+zbar2(n)*yni(n,2)*vrad2(n))/yne
	xltem(n) = (gamheate(n)/(yne*xk*tel(n))-2.5*vrade)
     1	/(chiel(n))
c	if(rhor(n).gt.pedrhote) xltem(n) = 1./ylte
c	if(rhor(n).le.pedrhote) xltem(n) = 1./xltetop  
	gamconde(n) = (gamheate(n)/(yne*xk*tel(n))-2.5*vrade)
c	***************************infer chi from exp***********************
c	if(ioptran.eq.1) then
c		if(rhor(n).gt.pedrhoti) xltim(n) = 1./yltibarx
c		if(rhor(n).gt.pedrhote) xltem(n) = 1./yltebarx 
c	endif
c	********************************************************************

 
c	integrate ion & electron temperatures inward from separatrix
	tiold(n) = ti(n)
	teold(n) = tel(n)
	ti(n) = ti(n+1)*
     1	(1.+delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))/
     2	(1.-delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))
	ti(n)=0.5*(ti(n)+tiold(n))
c	if(n.eq.24) ti(n) = ti(n+1)/(1.-delma*xltim(n))		
      if(ti(n).lt.0.0) ti(n) = 100.
	tel(n) = tel(n+1)*
     1	(1.+delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))/
     2	(1.-delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))
	tel(n)=0.5*(tel(n)+teold(n)) 
c	if(n.eq.24) te(n) = te(n+1)/(1.-delma*xltem(n))		
      if(tel(n).lt.0.0) tel(n) = 100.
C	INTEGRATE TOROIDAL ROTATION FREQ INWARD FROM SEPARATRIX
	call torotate(n)
	     
 	DO 1155 mm = 1,6
	DO 1150 kk = 1,6
	xsv = (AV(mm,kk,j)+AV(mm,kk,j-1))/delma - BV(mm,kk,j)
	SV(mm,j-1) = SV(mm,j-1)+SV(mm,j) + xsv*OMEGT(j,kk)/DELMA
	BV(mm,kk,j-1) = BV(mm,kk,j-1) +(AV(mm,kk,j-1)+AV(mm,kk,j))/DELMA
1150	CONTINUE 
	BBV(mm,kk) = BV(mm,kk,j-1)
	SSV(mm) = SV(mm,j-1) 
1155  CONTINUE
	if(nm.eq.1) goto 1204 
	CALL LSLRG(6,BBV,6,SSV,1,SSV)
	DO 1160 kk = 1,6
	OMEGT(j,kk) = SSV(kk)
1160	CONTINUE	
1204	CONTINUE      
	goto 675			  
C	NEOCLASSICAL CHI FOR IONS
c	simple neoclassical chi 
	ep = aminor*SQRT((1.+ELONG**2)/2.)/rmajor 
	ep = rhor(n)*aminor/rmajor
	bfield = abs(bphi) 
	OMI =EQ*BFIELD/XMAS(1)
	CSOUND = SQRT(XK*TEL(N)/XMAS(1))
 	rhot = csound/omi
 
	CHINC(n) = ((RHOTi(1)*bfield/bthet)**2)*xnuc12(n)*(EP**0.5) 
C	CHANG-HINTON CHI
	ALFA = XNi(2)*(zbar2(n)**2)/(yni(n,1)*(ATNUM(1)**2))
	qa = ep*bfield/bthet
	XMUii =(xnuc(1,1)*Q95*RMAJOR/(vth(1)*(EP**1.5)))*(1.+1.54*alfa)
	
	dp = 0.
600	G1 = (1. + 1.5*((EP**2)+ep*dp)+.375*(ep**3)*dp)/(1.+.5*ep*dp)
	G2 =SQRT(1.-(EP**2))*(1.+0.5*ep*dp)/(1.+(dp/ep)*(sqrt(1.-ep**2)
	1	-1))
	A1 =(0.66*(1.+1.54*ALFA)+(1.88*SQRT(EP)-1.54*EP)*(1.+3.75*ALFA))/
	1	(1.+1.03*SQRT(XMUii)+0.31*XMUii)
	A2 =0.59*XMUii*EP*(1.+1.33*ALFA*(1.+0.6*ALFA)/(1.+1.79*ALFA))/
     1	(1.+0.74*XMUii*(EP**1.5))	 
	
	betap = 2.*yni(n,1)*xk*ti(n)/((bthet**2)/(2.*1.257e-6)) 
	CHICH = CHINC(n)*(xnuc(1,1)/xnuc(1,2))*(A1*G1+A2*(G1-G2))
	if(dp.eq.0) then
		chich0 = chich
		dp = -1.*ep*(betap+0.5*log(1.65+0.89*(qa-1.)))
  		goto 600
	endif
	
	
c	orbit squeezing
		sheare(n) = 1.0
	if(n.lt.22) goto 625
	dEdr = (erada(n+1)-erada(n))/delma 
	
	para = 1./(bthet*vth(1))
c	if(ioptshear.eq.0) sheare(n) =	1.-(rhoti(1)/abs(fp))*dlnEdr*para
	sheare(n) =	1.-(rhoti(1)/abs(fp))*dEdr*para
	if(n.eq.24) sheare(25) = sheare(24)
	if(abs(sheare(n)).lt.1.0) sheare(n) = 1.0
 
625	chinc(n) = chinc(n)/(abs(sheare(n))**1.5)
	chichos(n) = chich/(abs(sheare(n))**1.5)
	

C	ITG-MODE CHI FOR IONS
	
c	ETAI(n) = XLTIM(N)/XLNM(N)
	
	CHIETAI(N) = 1.25*((CSOUND**2)*RHOT/OMI)*
     1	SQRT(XLTIM(N)/RMAJOR)
650	chiion(n) = cncmult*chichos(n) + cetaimult*chietai(n)
	if(ioptran.eq.1) then
		chiion(n) = chixpi
		if(rhor(n).lt.pedrhoti) chiion(n) = chitop  
	endif 
C	ETG-MODE CHI FOR ELECTRONS
	EMASS = 9.1E-31
	CLIGHT = 3.E8 
c	ETAE(n) =  XLTEM(N)/XLNM(N)
c	****exp etae**************
	etae(n) = 1.43
c	**************************
	CSE = SQRT(2.*XK*TEL(N)/EMASS)
	eplasfreq = 56.4*sqrt(YNI(N,1))
	CHIEETG(N)=0.13*((CLIGHT/EPLASFREQ)**2)*CSE*SHEARM*ETAE(n)*
	1		(1.+ETAE(n))/(Q95*RMAJOR)
C	TRAPPED ELECTRON MODE W/INTERPOLATION TO COLLISIONLESS DRIFT MODE
C		CHI FOR ELECTRONS 
	CSE = SQRT(2.*XK*TEL(N)/EMASS)
	RHOS = CSOUND/OMII(1)
	RHOTE = 3.37E-6*SQRT(TEL(N)/BFIELD)
	CHIEDW(N) = 2.5*(EP**1.5)*(CSOUND**2)*(RHOS**2)*XLNM(N)*XLTEM(N)/
	1	(XNUEI*(1.+0.1/XNUEIAST(N)))
	if(chiedw(n).gt.5.0) cedwmult = 0. 
C	do not use for collisionless regime because of previous results
c	IF(XNUEISTAR.LT.1)	CHIEDW(N) = 0.0
	chiel(n) = cetgmult*chieetg(n) + cedwmult*chiedw(n)
	if(ioptran.eq.1) then
		chiel(n) = chixpe
		if(rhor(n).lt.pedrhote) chiel(n) = chetop  
	endif
675	continue
c	iterate ion and neutral calculation to convergence
	do 680 n = 1,24
	xx = abs(yniold(n)/yni(n,1)-1.)
	if(abs(xx).gt.0.01) goto 375
	yy = abs(tiold(n)/ti(n) -1.)
	if(abs(yy).gt.0.01) goto 375
	zz = abs(teold(n)/tel(n)-1.)
	if(abs(zz).gt.0.01) goto 375  
680	continue
	if(kl.lt.10) goto 375
700	continue
c	now converge velocities, densities & temperatures

	if(iconverge.ne.1) goto 740
	iconverge = 0
	ioptapproach = 5
	goto 6


740	continue 
	OPEN(121,FILE='pedestal.TXT',STATUS='UNKNOWN') 

	
c	6/24/05 now solve for vtheta and vphi on fixed temperature & density profiles,
c		   using 1) neoclassical visc for poloidal calc and 2) both exper and gyro 
c		   drag frequencies.

c	TEMPORARY ****************************
	GOTO 1610 
C	***************************************
	do 1580 n=1,25
	vphia(1) = vphiex1(n)
  	vphia(2) = vphiex2(n)
	vtor1old(n) = vphiex1(n)
	vtor2old(n) = vphiex2(n)
	velthet1old(n) = velthet1(n)
	velthet2old(n) = velthet2(n)

1580	continue
	epv  = 0.05
	epd = 1.0
c	iterate vtheta-vphi solutions thru 1605
	do 1605 nit = 1,1 
	nn = 0
	do  1600 n = 1, 25
	j = 26 - n
	if(nit.gt.1) then
	vphia(1) = torvel(j,1)
	vphia(2) = torvel(j,2)

	endif
c
c	goto 1467
c	poloidal velocities 
      call poloidal(j)
	vtheta(1) = epd*vtheta(1)*VTH(1)*fp+(1.-epd)*velthet1old(n)
	vtheta(2) = epd*vtheta(2)*VTH(2)*fp+(1.-epd)*velthet2old(n)

	vpol(1,j) = vtheta(1)
	vpol(2,j) = vtheta(2)
 	velthet1(j) = vtheta(1)
	velthet2(j) = vtheta(2)
	if(ioptvisc.eq.1) then
c	gyroviscous drag
	XNUDRAGvis1(j)=XNUDRAG(1)
	XNUDRAGvis2(j)=XNUDRAG(2)
	endif
c	gyroviscous
	xnudragyro1(j) = xnudrag(1) 
      xnudragyro2(j) = xnudrag(2) 
	 
c	inertial terms
1467	AM = (aminor*SQRT(0.5*(1.+ELONG**2)))
	ep = am*rhor(j)/rmajor 
	ep = rhor(j)*aminor/rmajor
 	rminor=am*rhor(j)
c	xnuinert1(n) = (vrad1(n)/rmajor)*(1.-rmajor*xlvm(n)) -
c	1	0.5*(ep*velthet1(n)/rmajor)*thetinert(1)
c	xnuinert2(n) = (vrad2(n)/rmajor)*(1.-rmajor*xlvm(n)) -
c    	1	0.5*(ep*velthet2(n)/rmajor)*thetinert(2)
c	convection terms 9/30/05
c	xnuinert1(n)=(1./rminor)+(xnuioni(n)+xnuionb(n))/vrad1(n)-xlvm(n)
c	xnuinert2(n)=0.0	
c	anomalous drag
c	xnudraganom1(j) = xnudtot1(j)-(xnuioni(j) + xnuati(j) +xnuionb(j))
c	1			-xnudragyro1(j)-xnuinert1(j)
c	xnudraganom2(j) = xnudtot2(j)-xnudragyro2(j)-xnuinert2(j)

 
c	pressure gradients
	do 1590 k=1,2
 	xz = atnum(1)
	if(k.eq.2) xz = zbar2(j)
	PRESS(k)=-1.*(ti(j)/(xz*BTHET))*xlpm(j)
1590	continue
c	nbi momentum input
	xmtor(1) = xmomtor1(j)  
	XMTOR(2) = XMOMTOR2(j)   
c	density and mom trans freqs 
	xni(1) = yni(j,1)
	xni(2) = yni(j,2)
	xnuc(1,2) = xnuc12(j)
	xnuc(2,1) = xnuc21(j)
c	exper inferred drag
	ynud1 = ynudrag1(j)/xnuc(1,2)
	ynud2 = ynudrag2(j)/xnuc(2,1)
c	gyro + atomic + anomalous drag
	xnudragatomic(j) = xnuioni(j) + xnuati(j) + xnuionb(j) 
c     	xnudraganom1(j) = 0.0
c	xnudraganom2(j) = 0.0
 	ynud1=(xnudragyro1(j)+xnudragatomic(j)+xnudraganom1(j)+
     1	xnuinert1(n))/xnuc(1,2)
	ynud2 = (xnudragyro2(j)+xnudraganom2(j)+xnuinert2(j))/xnuc(2,1)
c	ynud1= xnudtot1(n)/xnuc(1,2)
c	ynud2= xnudtot2(n)/xnuc(2,1)
 

	ynudrag1(j) = ynud1*xnuc(1,2)
	ynudrag2(j) = ynud2*xnuc(2,1)
	bv1 = ynud1
	bv2 = ynud2
c	toroidal velocities
	y1 = xmtor(1) 
	y2 = xmtor(2) 
 	pressd(n) = press(1)-press(2)
	veld(j) = (velthet1(j)-velthet2(j))/fp
      brack(j) = ((velthet1(j)-velthet2(j))/fp - (press(1) -press(2)))
	
c	brack(j) = 0.0 
     	vtor1(j) = ((y1+y2)/(xni(1)*xmas(1)*xnuc(1,2))+
     1  ynud2*brack(j))/(ynud1+ynud2)
	vtor1(j) = vtor1(j) + (atnum(1)*eq*(bthet*gamion(j,1) +
     2		xni(1)*ephia*atnum(1))+zbar2(j)*eq*xni(2)*ephia)/
     3		(xni(1)*xmas(1)*xnuc(1,2)*(ynud1 + ynud2)) 
c	vtor2(n) = vtor1(n) - 
c	1	(velthet1(n)-velthet2(n))/fp + (press(1) -press(2))
	vtor2(j) = vtor1(j) - brack(j)
c	vtor2(n)= brack(n)/ynud2 +
c	1  (y2+xni(2)*zbar2(n)*eq*ephia)/(xni(2)*xmas(2)*xnuc(2,1)*ynud2)
c	vtor2(n) = (1.+ynud1)*vtor1(n) -
c    1	(atnum(1)*eq*(bthet*gamion(n,1)+xni(1)*ephia) + y1)/
c   2	(xni(1)*xmas(1)*xnuc(1,2))
	vtor2(j) = vtor1(j)/(1.+ynud2) + 
	1(atnum(2)*eq*xni(2)*ephia + y2)/(xni(1)*xmas(1)*xnuc(1,2))
     
     	
					
     	vtor1(j) = epd*vtor1(j)+(1.-epd)*vtor1old(j)
	vtor2(j) = epd*vtor2(j)+(1.-epd)*vtor2old(j)

c	erfb1(n) = bthet*(vtor1(n) - vtheta(1)/fp + press(1))
	ioptpress = 1
c	if(ioptpress.eq.1) press(2) = bpcarb(n)/bthet 
	erfb2(j) = bthet*(vtor2(j) - vtheta(2)/fp + press(2))
	eradfb2(j) = bthet*(torv(j) - vthexp(j)/fp + press(2))


1600	continue 
c	integrate toroidal velocities inward from separatrix
	torvel(25,1) = torv(25)
	torvel(25,2) = torv(25)
	veltor(25,1) = torv(25)
	veltor(25,2) = torv(25)
	nn = 0.0
	do 1603, nm = 1,24
	j = 25-nm
	call rotate(j)
	epad = 1.0
	torvel(j,1) = epad*torvel(j,1) + (1.-epad)*vtor1old(j)
	torvel(j,2) = epad*torvel(j,2) + (1.-epad)*vtor2old(j)
c	check convergence of vtor
 	xx1 = abs(torvel(j,1)/vtor1old(j)-1.0) 
	xx2 = abs(torvel(j,2)/vtor2old(j)-1.0)
	if(xx1.gt.epv) nn = nn+1
	if(xx2.gt.epv) nn = nn+1
	vtor1old(j) = torvel(j,1)
	vtor2old(j) = torvel(j,2)
	velthet1old(j)=velthet1(j)
	velthet2old(j)=velthet2(j)
 
1603	continue      
     
      
	if(nn.eq.0) goto 1610
	do 1604, nm = 1,25
	j = 26-nm
C	INTEGRATE TOROIDAL ROTATION FREQ INWARD FROM SEPARATRIX
	CALL TOROTATE(j)
	DO 1555 mm = 1,6
	DO 1550 kk = 1,6
	xsv = (AV(mm,kk,j)+AV(mm,kk,j-1))/delma - BV(mm,kk,j)
	SV(mm,j-1) = SV(mm,j-1)+SV(mm,j) + xsv*OMEGT(j,kk)/DELMA
	BV(mm,kk,j-1) = BV(mm,kk,j-1) +(AV(mm,kk,j-1)+AV(mm,kk,j))/DELMA
1550	CONTINUE 
	BBV(mm,kk) = BV(mm,kk,j-1)
	SSV(mm) = SV(mm,j-1) 
1555  CONTINUE
	if(nm.eq.1) goto 1604 
	CALL LSLRG(6,BBV,6,SSV,1,SSV)
	DO 1560 kk = 1,6
	OMEGT(j,kk) = SSV(kk)
1560	CONTINUE	
1604	continue 

1605	continue
	write(6,'(1x,25A)')'vtheta-vphi iteration not converged' 
	write(121,'(1x,25A)')'vtheta-vphi iteration not converged' 
1610	continue 
c	goto 1594 
C	INFER CHI's FROM EXP N,T,GSCL & CALC. Q, GAM
C	DO 1650 N= 1,25
C	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
C	XCHIE(N) = EXLTE(N)*((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
C     1					2.5*GAMEL/EXNE(N))	 
C	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))
C	XCHII(N) = EXLTI(N)*((GAMHEATI(N)/(XNION*XK*XTI(N)))-
C    1					2.5*GAMION(N,1)/XNION)
C1650	CONTINUE     	 
c	inference of exp chi's*************************************
c	new 12/5/05******initialize to input exp data**************
	 XMPROT = 1.673E-27
	EQ = 1.6E-19
	BTHET = EP*abs(BPHI)/QSAFE
	FP = BTHET/BPHI	
	XK = 1.6E-19 
	XMASSELECTRON = 9.11E-31
	XMAS(1) = AION*XMPROT
	XMAS(2) = AIMP*XMPROT  
	ATNUM(1) = ZION
	ATNUM(2) = ZIMP	  
	do 1350 n= 1,25
	zne(n) = exne(n)
	if(zbar2(n).eq.0.0) zbar2(n) = atnum(2) 
	yni(n,1) = exne(n)/(atnum(1)+fracz*zbar2(n)) 
	yni(n,2) = fracz*yni(n,1)
	ti(n) = xti(n)
	tel(n) = xte(n)
	teed(n) = tel(n)
	tied(n) = ti(n)
	gamion(n,2) = 0.0
	xlnm(n) = 1./xlne(n)
 	xlnim(n)= 1./xlne(n)
	xlnzm(n)= 1./xlne(n)
 	xltim(n)= 1./exlti(n)
	xltem(n)= 1./exlte(n)
	xlvm(n) = 1./exlv(n)
	xlpm(n) = xlnm(n) + xltim(n) 
	xnudrag(1) = xnudtot1(n)
	diff(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
     1			((xnudrag(1)/xnuc12(n))+1.-atnum(1)/zbar2(n))/
     2			((eq*atnum(1)*bthet)**2) 		
1350	continue
	rhor(25) = 1.0 	 
	DO 1355, NN=1,24
 	N = 25-NN
	rhor(n) = rhor(n+1) - delna/(aminor*SQRT(0.5*(1.+ELONG**2)))
1355	continue
c	calculate neutral density distribution
	call neutdist  
c	calculate atomic physics  	
	do 1360 n = 1,25
	ted = xte(n)
	tid = xti(n)
	tnd = xti(n)
	ynd = yni(n,1) 
	CALL INTERP(TED,TID,TND,YND)
	SVEL = SEL(1)
	SVELN= SELN(1)
	SVCX = SCX(1)
	SVATA(n)= (SEL(1) + CNEUT*SELN(1)*YNO(n)/YNI(n,1))+SCX(1)
	SVIONA(n) = SION(1) 
	SVREC = RECOM(1)
	XNUIONI(n) = YNO(n)*SVIONA(n)
	VCOLD = SQRT(XK*TSPL/XMASS)
	ynocold(n) = coldno(n)
	if(n.eq.25 )ynocold(25) = 0.5*(gamouteff/vcold + coldno(25))
      XNUATI(n) =  YNOCOLD(n)*SVATA(n)
		EIONe = 17.5
	IF(ynd.LE.1.E21) 
     2     EIONe = 17.5 + (5.+37.5/Ted)*LOG10(1.E21/ynd)
    	IF(ynd.GT.1.E21)
     2    EIONe = (30.6 - 16.4*EXP(-5.E19/ynd))*
     3              EXP(5.45/(Ted*EXP((ynd/1.37E20)**0.26)))	
	alphael(n) = xnuioni(n)*(2.5*(xnu-1.)+xnu*eione/ted)
	alphaion(n) = 1.5*xnuati(n)*(xnu-1.)
1360	continue
c	calculate heat and particle flux distributions at separatrix
	gamheat(25) = fluxheat
	gamheate(25) = fheate*gamheat(25)
	gamheati(25) =(1.-fheate)*gamheat(25)
      
	GAMION(25,1) = enh*FLUXPART
	xnsol = yni(25,1)
c	calculate radiation at separatrix
		TDBL = 0.5*(Tel(24)+tel(25))
  	eTDBL = 0.5*(Tel(24)+tel(25))

	IZ1 = izintrin
	IZ2 = atnum(1)
	if(iz1.eq.4.or.iz1.eq.6.or.iz1.eq.74) then
	fon = yno(25)/(0.5*(yni(24,1)+yni(25,1)))
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,eTdbl,fon,exlzdbl,edlzdbl,ZAV)
	zbar2(25) = zav
	zne(25) = yni(25,1) + zbar2(25)*yni(25,2)
	goto 1368 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
1368	XLradZ(25)= eXLZDBL

	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(25) = XLradZ(25)*1.e-13
	dlzdte = edlzdbl*1.e-13
	alphael(25) = alphael(25) + yni(25,2)*(xlradz(25)/etdbl-dlzdte)
	xni(1) = yni(25,1)
	xni(2) = yni(25,2)
	vrad1(25) = gamion(25,1)/yni(25,1)
	do 1575 j = 1,24
	n = 25-j
c	new ion & plasma density and ion & electron temp formulation
	dens(n) = 0.5*(yni(n+1,1)+yni(n,1) )
	tele = 0.5*(tel(n)+tel(n+1))
	tiav = 0.5*(ti(n)+ti(n+1))
	rz = (atnum(1)**2/xmas(1) + zbar2(n)**2/xmas(2))*1.67e-27  
	cequil = 6.32e-14*rz/(tele**1.5) 
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)

	qie(n) = cequil*yne*(tiav-tele)
	cxcool(n)=1.5*dens(n)*tiav*xk*xnuati(n)
c	if(n.eq.24) cxcool(n) = 1.5*dens(n)*tiav*xk*xnuati(n)
c**********************
	cmulteq = 1.0

c	electron heat flux
	EIONi = 17.5
	IF(dens(n).LE.1.E21) 
     2     EIONi = 17.5 + (5.+37.5/Tel(n))*LOG10(1.E21/dens(n))
    	IF(dens(n).GT.1.E21)
     2    EIONi = (30.6 - 16.4*EXP(-5.E19/dens(n)))*
     3              EXP(5.45/(Tel(n)*EXP((dens(n)/1.37E20)**0.26)))
	TDBL = 0.5*(Tel(n)+tel(n+1))
	eTDBL = 0.5*(Tel(n)+tel(n+1))
	IZ1 = izintrin
	IZ2 = atnum(1)
	if(iz1.eq.4.or.iz1.eq.6.or.iz1.eq.74) then
	fon = yno(n)/(0.5*(yni(n,1)+yni(n+1,1)))
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,eTdbl,fon,exlzdbl,edlzdbl,ZAV)
	zbar2(n) = zav
	goto 1565 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
1565	XLradZ(n)= eXLZDBL
	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(n) = XLradZ(n)*1.e-13
	dlzdte = edlzdbl*1.e-13
 	alphael(n) = alphael(n) + yni(n,2)*(xlradz(n)/tel(n)-dlzdte) 
	coolion(n) = xk*eioni*yne*xnuioni(n)
	radcool(n) = 0.5*(yni(n,2)+yni(n+1,2))*yne*xlradz(n)
	radcool(n) = radmultedge*radcool(n)

	call edgerotran(n,10)
c	particle fluxes on ions & impurities
 	thetint = 0.5
c	DENS(N) = thetint*YNI(N+1,1) + (1.-thetint)*YNI(N,1)
	delma = delna 
	GAMION(N,1) = GAMION(N+1,1) - dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n)) 
	GAMION(N,2) = GAMION(N+1,2) 
     1		-	DENS(N)*XNUIONz(N)*DELMA*(1.+fracz*zbar2(n))
	

   	beamdot(n) =  dens(n)*xnuionb(n)
	
	gamheati(n) = gamheati(n+1) + delma*(cxcool(n) + cmulteq*qie(n))
 	gamheati(n) = gamheati(n) - fionb(n)*qnb(n)*delma
		 
	gamheate(n)=gamheate(n+1)+delma*
     1					(coolion(n)+radcool(n)-cmulteq*qie(n))
	gamheate(n) = gamheate(n) -	(1.-fionb(n))*qnb(n)*delma
	x=gamheate(n)
	y=gamheati(n)
	z=gamion(n,1)
c	time dependent
	gamheati(n) = gamheati(n)-dlnw_dt*1.5*yni(n,1)*xk*ti(n)*delma
	gamheate(n) = gamheate(n)-dlnw_dt*1.5*yne*xk*tel(n)*delma
	gamion(n,1) = gamion(n,1)-dln_dt*yni(n,1)*delma
	x=gamheate(n)
	y=gamheati(n)
	z=gamion(n,1)


c	poloidal rotation
	call poloidal(n)
	xnudragyro1(n) = xnudrag(1)
	xnudragyro2(n) = xnudrag(2)
	velthet1(n) = vtheta(1)*VTH(1)*fp
	velthet2(n) = vtheta(2)*VTH(2)*fp
 
1575  continue     	
c*******************	
	call poloidal(25)
	xnudragyro1(25) = xnudrag(1)
	xnudragyro2(25) = xnudrag(2)
 	velthet1(25) = vtheta(1)*VTH(1)*fp
 	velthet2(25) = vtheta(2)*VTH(2)*fp

C	INFER CHI's FROM EXP N,T,GSCL & CALC. Q, GAM
	do 1589 n = 1,25
 	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
	XCHIE(N) = EXLTE(N)*((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
     1					2.5*GAMEL/EXNE(N))	 
	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))
	XCHII(N) = EXLTI(N)*((GAMHEATI(N)/(XNION*XK*XTI(N)))-
     1					2.5*GAMION(N,1)/XNION)
c	convective fractions
	econve(n) = 2.5*gamel*xk*xte(n)/gamheate(n) 
	econvi(n) = 2.5*(gamion(n,1)+gamion(n,2))*xk*xti(n)/gamheati(n) 
c	calculate theoretical chi's

1589	continue 
	do 7100 nmesh =1,25 
	n= nmesh
	EQ = 1.6E-19
	XK = 1.6E-19
	EP0 = 8.854E-12
	AM = (aminor*SQRT(0.5*(1.+ELONG**2)))
 	ep = aminor*rhor(n)/rmajor 
	fc = (((1.-ep)**2)/sqrt(1.-(ep**2)))/(1.+1.46*sqrt(ep)+0.2*ep)

	do 107 k=1,2
	xni2=yni(n,2)
	if(k.eq.1) xni2=yni(n,1)
	do 107 j=1,2 
	z2 = atnum(k)
	z1 = atnum(j)
	if(k.eq.2) z2=zbar2(n)
	if(j.eq.2) z1=zbar2(n)
	Y = SQRT(XNI2)*(z2**2)*z1 
 	X = (EP0/EQ)**1.5	
	COULOG(J,K) = LOG(12.*3.1416*(xti(n)**1.5)*X/Y)
107	continue
c		BRAGINSKI COLLISION FREQUENCIES
	XMR11 = XMAS(1)*(1.+XMAS(1)/XMAS(1))
	XMR12 = XMAS(1)*(1.+XMAS(1)/XMAS(2))
 	XMR21 = XMAS(2)*(1.+XMAS(2)/XMAS(1))
	XMR22 = XMAS(2)*(1.+XMAS(2)/XMAS(2))

	XNI1 = EXNE(N)/(1.+ ZBAR2(N)*FRACZ)
	XNI2 = FRACZ*XNI1
	C1 = 1./((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5)) 
	XNUC(1,1) = 3.34*(COULOG(1,1)*(ATNUM(1)**4)*1.E-6*XNI1)/
	2			(C1*SQRT(XMR11*1E3)*(XTI(N)**1.5))
	XNUC(1,2)=3.34*(COULOG(1,2)*((ATNUM(1)*ZBAR2(N))**2)*1.E-6*XNI2)
	2			/(C1*SQRT(XMR12*1E3)*(XTI(N)**1.5))
 	XNUC(2,1)=3.34*(COULOG(2,1)*((ATNUM(1)*ZBAR2(N))**2)*1.E-6*XNI1)
     2			/(C1*SQRT(XMR21*1E3)*(XTI(N)**1.5))
	XNUC(2,2) = 3.34*(COULOG(2,2)*(ZBAR2(N)**4)*1.E-6*XNI2)/
	2			(C1*SQRT(XMR22*1E3)*(XTI(N)**1.5))
	xnuc12(n) = xnuc(1,2) 



	XVTHI = SQRT(2.*XK*XTI(N)/XMAS(1))
	XVTHIMP = SQRT(2.*XK*XTI(N)/XMAS(2))
	EMASS = 9.1E-31
	xvthe = sqrt(2.*xk*xte(n)/emass)
	CLIGHT = 3.E8 
	CSE = SQRT(2.*XK*XTE(N)/EMASS) 
	XNELECTRON = XNI(1)*(ATNUM(1)**2) + XNI(2)*(zbar2(n)**2) 
	XNUEI =EXNE(N)/(6.4E14*((1.E-3*XTE(N))**1.5)) 
	XNUEIAST(N) = XNUEI*ABS(QSAFE)*RMAJOR/XVTHE

	ep = aminor*rhor(n)/rmajor 
	bfield = abs(bphi) 
	OMI =EQ*BFIELD/XMAS(1)
	CSOUND = SQRT(XK*XTE(N)/XMAS(1))
 	rhot = csound/omi
C	NEOCLASSICAL CHI FOR IONS
c	simple neoclassical chi 
 
	CHINC(n) = ((RHOTi(1)*bfield/bthet)**2)*XNUc(1,2)*(EP**0.5) 
C	CHANG-HINTON CHI
	ALFA = XNi2*(zbar2(n)**2)/(XNi1*(ATNUM(1)**2))
	qa = ep*bfield/bthet
	XMUii =(XNUc(1,1)*Q95*RMAJOR/(XvthI*(EP**1.5)))*(1.+1.54*alfa)
	
	dp = 0.
6100	G1 = (1. + 1.5*((EP**2)+ep*dp)+.375*(ep**3)*dp)/(1.+.5*ep*dp)
	G2 =SQRT(1.-(EP**2))*(1.+0.5*ep*dp)/(1.+(dp/ep)*(sqrt(1.-ep**2)
	1	-1))
	A1 =(0.66*(1.+1.54*ALFA)+(1.88*SQRT(EP)-1.54*EP)*(1.+3.75*ALFA))/
	1	(1.+1.03*SQRT(XMUii)+0.31*XMUii)
	A2 =0.59*XMUii*EP*(1.+1.33*ALFA*(1.+0.6*ALFA)/(1.+1.79*ALFA))/
     1	(1.+0.74*XMUii*(EP**1.5))	 
	
	betap = 2.*Xni1*xk*Xti(n)/((bthet**2)/(2.*1.257e-6)) 
	CHICH = CHINC(n)*(xnuc(1,1)/xnuc(1,2))*(A1*G1+A2*(G1-G2))
	if(dp.eq.0) then
		chich0 = chich
		dp = -1.*ep*(betap+0.5*log(1.65+0.89*(qa-1.)))
  		goto 6100
	endif
	
	
c	orbit squeezing
	
c	if(n.eq.25) then
c	dlnEdr = log(erada(n)/erada(n-1))/delma 
c	goto 6110
c	endif
	sheare(n) = 1.0
c	if(n.lt.10) goto 6125
	dEdr = (erada(n+1)-erada(n))/delma 
	
6110	para = 1./(bthet*XvthI)
c	if(ioptshear.eq.0) sheare(n) =	1.-(rhoti(1)/abs(fp))*dlnEdr*para
	sheare(n) =	1.-(rhot/abs(fp))*dEdr*para
	if(n.eq.24) sheare(25) = sheare(24)
	if(abs(sheare(n)).lt.1.0) sheare(n) = 1.0
c6125	chinc(n) = chinc(n)/(abs(sheare(n))**1.5)
	sheare(n) = 1.0 
	chichos(n) = chich/(abs(sheare(n))**1.5)
	
     
      dngscl = XLNE(N)
	dtIgscl = EXLTI(N)
	DTEGSCL = EXLTE(N)

	etai(n) = dngscl/dtigscl
	etae(n) = dngscl/dtegscl

C	ITG-MODE CHI FOR IONS
	
c	ETAI(n) = XLTIM(NMESH)/XLNM(NMESH)
c	ETAi(n) =  ylnbarx/yltibarx
c 	if(rhor(n).lt.pedrhon) etai(n) = ylntop/xltitop
	CHIET = ((CSOUND**2)*RHOT/OMI)*
     1	SQRT(1./(EXLTI(N)*RMAJOR))
	chietair(n) = 1.25*chiet
	tauw = xti(n)/(zeff*xte(n))
	fweiland = 2.64*(1.-.6667*exlti(n)/xlne(n))*tauw/fc -
     1           (exlti(n)/rmajor)*((((2.-rmajor/xlne(n))/fc)**2)/8. +
     2		   20.*(tauw**2)/9.)
      chietaiw(n) = chiet*fweiland 	  
	chiion(n) = cncmult*chichos(n)+cetaimult*chietai(n)
	if(ioptran.eq.1) then
	CHIION(N) = XCHII(N)
C	chiion(n) = chixpi
C		if(rhor(n).lt.pedrhoti) chiion(n) = chitop  
 	endif 
   	
C	ETG-MODE CHI FOR ELECTRONS
	
c	
c	****exp etae**************
c	etae(n) = 1.43
c	ETAe(n) = XLTeM(NMESH)/XLNM(NMESH)

c	ETAE(n) =  ylnbarx/yltebarx
c	if(rhor(n).lt.pedrhon) etae(n) =  ylntop/xltetop

c	**************************
	etae(n) = xlne(n)/exlte(n)
	shearm = am*rhor(n)*ssi95
	eplasfreq = 56.4*sqrt(EXNE(N))
	CHIEETG(N)=0.13*((Clight/EPLASFREQ)**2)*xvthe*SHEARM*ETAE(n)*
	1		(1.+ETAE(n))/(Q95*RMAJOR)
c	if(chieetg(n).gt.5.0) chieetg(n) = 1.0 
C	TRAPPED ELECTRON MODE W/INTERPOLATION TO COLLISIONLESS DRIFT MODE
C		CHI FOR ELECTRONS 
	
	RHOS = CSOUND/OMII(1)
	RHOTE = 3.37E-6*SQRT(XTE(N)/BFIELD)
	ylnylte = dngscl*dtegscl
c	if(rhor(n).lt.pedrhon) ylnylte = ylntop*xltetop
C	if(xlnm(n).ne.0.0)
c	wesson
      CHIEDW(N) = 2.5*(EP**1.5)*(CSOUND**2)*(RHOS**2)/(ylnylte*
     1	XNUEI*(1.+0.1/XNUEIAST(N)))
c	Kalupin NF,45,468(2005)
	zz = (xnuei/ep)/(Xte(n)/(bfield*dngscl)) 
c	chiedw(n) = sqrt(2./ep)*etae(n)*xnuei*(rhos**2)/(1.+zz**2) 
c	if(chiedw(n).gt.5.0) chiedw(n) = 1.0 
C	do not use for collisionless regime because of previous results
c	IF(XNUEISTAR.LT.1)	CHIEDW(N) = 0.0
	chiel(n) = cetgmult*chieetg(n) + cedwmult*chiedw(n)
	if(ioptran.eq.1) then
	CHIEL(N) = XCHIE(N)
C		chiel(n) = chixpe
C		if(rhor(n).lt.pedrhote) chiel(n) = chetop  
	endif
C	RESISTIVE BALLOONING MODE CHI FOR ELECTRONS
	res = (2.8e-8)*zeff/(((1.e-3)*Xte(n))**1.5)
	res = sqrt(9.1e-31*xk)*zbar2(n)*17./
	1		(6.*sqrt(3.)*3.14*((8.85e-12)**2)*(xte(n)**1.5))
	beta =2.*EXne(n)*xk*Xte(n)*(2.*1.257e-6)/((bfield**2))
	XLPIM = (1./EXLTI(N))+(1./XLNE(N))

	alpha = 2.*1.257e-6*rmajor*(q95**2)*Xni1*xk*Xti(n)*xlpIm/
	1		(bfield**2)
	shearq = 1.0 
      chierb(n) = 10.0*res*sqrt(beta*3.343/9.1e-4)*((alpha/shearq)**1.5) 
c	gudzar pf,5, 3712,1993
	
	pi = 3.14159    
	chierb(n)= ((2.*q95*rhote)**2)*xnuei*rmajor/dngscl
c	paleoclassical for electrons (NF,45,1120,2005)
	AM = (aminor*SQRT(0.5*(1.+ELONG**2)))
	ep = aminor*rhor(n)/rmajor 
	fc = (((1.-ep)**2)/sqrt(1.-(ep**2)))/(1.+1.46*sqrt(ep)+0.2*ep)

	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(eXNe(n)) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(xte(n)**1.5)*X/Yz)
	ome = 56.4*sqrt(exne(n))
	xlame = 1.2e16*(xte(n)**2)*17./(exne(n)*zeff*couloge)
	xnuestar = (rmajor*q95)/((ep**1.5)*xlame)
	

	xmunu = ((zeff+1.414 - log(2.414))/
	1		(zeff*(1 +sqrt(xnuestar)+xnuestar)))*(1.-fc)/fc
	const=(1.414+zeff)/(1.414+(13./4.)*zeff) + 	xmunu  
      		
      bigL = xlame
	vthe = sqrt(2.*xk*xte(n)/xme)
	xnue = (vthe/xlame)
	dele = (3.e8/ome)
	xlmax = rmajor*sqrt(3.1416*q95)/sqrt(dele*ssi95/aminor)
	if(xlmax.lt.xlame) bigL = xlmax 
	chiepale(n)=1.5*(1. + bigL/(3.1416*q95*rmajor))*const 
	1			*xnue*(dele**2)
	chiepale(n)=1.5*(1. + bigL/(3.1416*q95*rmajor))*const
	1			*(1.4e3*zeff/xte(n)**1.5)*(couloge/17.)
	xx = xnue*(dele**2)
	xy = (1.4e3*zeff/xte(n)**1.5)*(couloge/17.)
c	drift alfven (CPP,38,118,1998) 
	exlpe = 1./((1./xlne(n))+(1./exlte(n))) 
	chigb = (rhot**2)*csound/exlpe
	xmuu = (exlpe/(q95*rmajor))*sqrt(xmas(1)*xte(n)/(xme*xti(n)))
	beta = 2.*1.257e-6*exne(n)*xk*xte(n)/(bfield**2)
	betan =sqrt(xmas(1)/xme)*beta*q95*rmajor/exlpe
	xnun =((xmas(1)/xme)**0.25)*sqrt(q95*rmajor*exlpe)/xlame
	chiperp = (1./((1.+betan**2)**3)) + xnun**2
	chiperp = chiperp/(1.+betan**2+xnun**1.333)
	chiperp = sqrt(chiperp)
	chida(n) = chigb*chiperp/sqrt(xmuu)
c	trapped electron mode drift wave
	ftr = 1.-fc
	chitem = (ep*ftr/xnue)*((csound*rhot)**2)/(xlne(n)*exlte(n))
	xxx= 0.1*cse/(xnue*q95*rmajor*(ep**1.5))
	fwesson = 2.5/(1. + xxx)
	chitemw(n) = chitem*fwesson
	yyy= csound*ep/(xlne(n)*xnue)**2.
	fkadom = 1./(1.+fkadom)
	chitemk(n) = chitem*fkadom

c	Bateman(PPCF,48,A93,2006) ETG
	r = rho(n)*aminor	
	xles = q95*rhote*rmajor/exlte(n)	
c	if(xles.ge.dele) goto 1805
	chietg0 = 0.06*(dele**2)*vthe/sqrt(exlte(n)*rmajor)
c	goto 1810
1805	chietg0 = 0.06*(q95**2)*((rmajor/exlte(n))**1.5)*
	1  (rhote/(eq*bfield))*xk*xte(n)*
     2	((1./exlte(n))-1.88*(ssi95/(q95*rmajor))*(1.+xte(n)/xti(n)))
1810	continue
c	etg cutoff
	cutoffetg = 0.8*rmajor/xlne(n)
	xw = (1.+zeff*xte(n)/xti(n))*(1.33+1.91*r*ssi95/q95)*(1.-1.15*ep)
	if(xw.gt.cutoffetg) cutoffetg = xw
	xa = rmajor/exlte(n)
	if(xa.le.cutoffetg) then
	chietg0 = 0.0
	chieetg(n) = 0.0
	endif 
	chieetg(n) = chietg0
c	chi due to thermal instabilities

c	atomic physics	
 	ted = xte(n)
	tid = xti(n)
	tnd = xti(n)
	ynd = exne(n) 
	CALL INTERP(TED,TID,TND,YND)
	SVEL = SEL(1)
	SVELN= SELN(1)
	SVCX = SCX(1)
	SVATA(n)= (SEL(1) + CNEUT*SELN(1)*YNO(n)/YNI(n,1))+SCX(1)
	SVIONA(n) = SION(1) 
	SVREC = RECOM(1)
	XNUIONI(n) = YNO(n)*SVIONA(n)
	VCOLD = SQRT(XK*TSPL/XMASS)
	ynocold(n) = coldno(n)
	if(n.eq.25 )ynocold(25) = 0.5*(gamouteff/vcold + coldno(25))
      XNUATI(n) =  YNOCOLD(n)*SVATA(n)
		EIONe = 17.5
	IF(ynd.LE.1.E21) 
     2     EIONe = 17.5 + (5.+37.5/Ted)*LOG10(1.E21/ynd)
    	IF(ynd.GT.1.E21)
     2    EIONe = (30.6 - 16.4*EXP(-5.E19/ynd))*
     3              EXP(5.45/(Ted*EXP((ynd/1.37E20)**0.26)))	
	alphael(n) = xnuioni(n)*(2.5*(xnu-1.) + xnu*eione/ted) +
	1			2.5*xnuioni(n)*(xnu-1.)
	alphaion(n) = 1.5*xnuati(n)*(xnu-1.) + 2.5*xnuioni(n)*(xnu-1.)
	
c	radiation
	eTDBL = xte(n)
	IZ1 = izintrin
	IZ2 = atnum(1)
	if(iz1.eq.4.or.iz1.eq.6.or.iz1.eq.74) then
	fon = yno(n)/exne(n)
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,eTdbl,fon,exlzdbl,edlzdbl,ZAV)
	zbar2(n) = zav
	goto 1965 
	endif

	CALL cefits (IZ1, eTDBL, eXLZDBL, 1, eDLZDBL)
1965	XLradZ(n)= eXLZDBL
	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(n) = XLradZ(n)/1.6e-12
	dlzdte = edlzdbl/1.6e-12
 	alphael(n) = alphael(n) + 
     1		cfzinttb*exne(n)*1.e-6*(xlradz(n)/xte(n)-dlzdte) 
	chii0 = chichos(n)
	chie0 = chiepale(n)
 
	omegthi(n) =0.667*((alphaion(n)-2.5*xnu*vrad1(n)/exlti(n)-
	1	chii0*xnu/(exlti(n)**2)) - chii0*(xkrtherm**2))
	omegthe(n) =0.667*((alphael(n)-2.5*xnu*vrad1(n)/exlte(n)-
	1	chie0*xnu/(exlte(n)**2)) - chie0*(xkrtherm**2))

	delchii(n) = omegthi(n)/(xkrtherm**2)
 	delchie(n) = omegthe(n)/(xkrtherm**2)

c	delchii(n)=0.667*((alphaion(n)-2.5*xnu*vrad1(n)/exlti(n)-
c	1	2.5*chii0*xnu/exlti(n))/(xkrtherm**2) - chii0)
c	delchie(n)=0.667*((alphael(n)-2.5*xnu*vrad1(n)/exlte(n)-
c    1	2.5*chie0*xnu/exlte(n))/(xkrtherm**2) - chie0)
 
		 






C	THERMAL INSTABILITIES

	XLINSTABI = 1.
	XLINSTABE = 1.
C	XLINSTAB IS THE RADIAL CORRELATION LENGTH OF INSTABILITY IN CM
	CHI0ION = CHINC(N)
	CHI0EL  = CHIEPALE(N)
 
C	CHI0ION & CHI0EL ARE THE BACKGROUND CHI'S IN M^2/S	 
	XKRMI2= 1.E-4
	XKRME2= 1.E-4

	CHITHI(N) = .667*(XKRMI2*(ALPHAION(N) - CHI0ION*XNU/(EXLTI(N)**2))
	1	- 2.5*XNU*VRAD1(N)/EXLTI(N) - CHI0ION)
	CHITHE(N) = .667*(XKRME2*(ALPHAEL(N) - CHI0EL*XNU/(EXLTE(N)**2))-
 	1	2.5*XNU*VRAD1(N)/EXLTE(N) - CHI0EL)

	

7100	continue 
	goto 1594
c	integrate ion & electron temps & ion dens inward from separatrix
	
	ti(25) = xti(25)
	tel(25) = xte(25)
c	yni(25,1)=EXNE(25)/(ATNUM(1)+FRACZ*ZBAR2(25))
 	zne(25) = exne(25)
	do 1593 j = 1, 24
	n = 25-j
	c6 = 1.0
	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
	xltem(N) = (((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
     1					2.5*GAMEL/EXNE(N)))/xchie(n)	 
	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))
	xltim(N) = (((GAMHEATI(N)/(XNION*XK*XTI(N)))-
     1					2.5*GAMION(N,1)/XNION))/xchii(n)
 
		
c	yni(n,1) = yni(n+1,1)*
c	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
c    2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
c	yni(n,2) = fracz*yni(n,1)
c	zne(n) = yni(n,1) + zbar2(n)*yni(n,2)

	ti(n) = ti(n+1)*
     1	(1.+delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))/
     2	(1.-delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))
c	if(n.eq.24) ti(n) = ti(n+1)/(1.-delma*xltim(n))	
	
      if(ti(n).lt.0.0) ti(n) = 100.
c	if(ti(n).gt.1.e3) ti(n) = 1.e3
	tel(n) = tel(n+1)*
     1	(1.+delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))/
     2	(1.-delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))
1593	continue
1594  continue 
c	goto 700
c	******************end heat transport investigation************	


	 
	OPEN(122,FILE='pedestalplot.TXT',STATUS='UNKNOWN') 
	if(jwarn.eq.1) write(121,199)

	write(121,'(1x,35A)') '   rho	     ne         ti         te     
     1  erad    OMEGT01'
	write(122,'(1x,35A)') '   rho	     ne         ti         te  ' 
c	sepdif = xnsepex - xnesepreal
     	do 750 n = 1,25
	zne(n) = zne(n) - sepdif
	write (121,1000) rhor(n),zne(n),ti(n),tel(n),erada(n),OMEGT(n,1)
	write (122,1000) rhor(n),zne(n),ti(n),tel(n)
	sion(n) = zne(n)*0.5*(ynuioni(n) + ynuioni(n+1))
c	if(n.ne.25) xlnm(n+1) = log(yni(n,1)/yni(n+1,1))/delna
750	continue
	
	write(121,'(1x,35A)') '   rho	   gamheate    gamheati   gamion     
     1 econve    econvi'
     	do 755 n = 1,25
	write (121,1000) rhor(n),gamheate(n),gamheati(n),gamion(n,1),
     1econve(n),econvi(n)
755	continue 
      write(121,'(1x,35A)') '   rho	    nudrag1*    nuatom1   nuinert1       
     1nugyro1  nurip1  nuanom1  '
     	do 760 n = 1,25
	write (121,1000) rhor(n),xnudtot1(n),xnuati(n)+xnuioni(n)+xnuionb(
     1n),xnuinert1(n),xnudragyro1(n),xnudragrip(n,1),xnudraganom1(n)
760	continue 
	write(121,'(1x,35A)') '   rho	    nudrag2*    nuatom2   nuinert2       
     1nugyro2  nurip2  nuanom2  '
	xcb = 0.0	 
     	do 761 n = 1,25
	write (121,1000) rhor(n),xnudtot2(n), xcb,
     1xnuinert2(n),xnudragyro2(n),xnudragrip(n,2),xnudraganom2(n)
761	continue 
 
	write(121,'(1x,35A)') '   rho	     xno      coldno    cxcool       
     1radcool  ioncool   qnbe  '
 
	do 765 n = 1,25
	write (121,1000) rhor(n),yno(n),coldno(n),cxcool(n),
     1radcool(n),coolion(n), qnbe(n)
765	continue 
	write(121,'(1x,35A)') '   rho	     lp        lti       ln    
     1 lte	   lv    qnbi '
   
 
	do 770 n = 1,25
	write (121,1000) rhor(n),1./xlpm(n),1./xltim(n),1./xlnm(n),
     1					1./	xltem(n),1./xlvm(n),qnbi(n)
     
          
770	continue 

	write(121,'(1x,35A)') '  rho     chiion     chiel     diff    vrad 
	1/D    vpin/D   NUEI* ' 
	do 771 n = 1,25
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	eradterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)))/
     1			(eq*bthet)
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((velthet1(n)/fp)))/
     1			(eq*bthet)

 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp5 = zbeam-ephia/bthet+vphiterm+eradterm+vtheterm

	write(121,1000) rhor(n),chiion(n),chiel(n),diff(n),
     1	vrad1(n)/diff(n),vp5/diff(n),XNUEIAST(N)
771	continue 
	WRITE(121,'(1X,35A)') '	 rho     xchi-ion    xchi-el  DIFF   CHITH
     1I    CHITHE     qie '
	do 7710 n = 1,25
	write(121,1000) rhor(n),xchii(n),xchie(n),DIFF(N),CHITHI(N),
     1	CHITHE(N),qie(n)

7710	continue	  
	write(121,'(1x,35A)') ' chiiexp    chiineo  chichos   chiida  chie
     1tair(n)  chietaiw(n)  chithermi  '					 
	do 772 n = 1,25
	write(121,1000) xchii(n),chinc(n),chichos(n), chida(n),chietair(n)
	1				,chietaiw(n),delchii(n) 

772	continue 
	write(121,'(1x,35A)') '  chieexp   chiepale  chierb    chitemw, 
     2chitemk  chieetg  chitherme  '   
	do 773, n= 1,25
	write(121,1000) xchie(n),chiepale(n),chierb(n),chitemw(n),
     1				chitemk(n),chieetg(n),delchie(n)

773	continue     	 
	write(121,'(1x,35A)') ' alphaion  chii(k+L)  omegi	 alphael  chie
	1(k+L)  omege  '
	do 774, n=1,25
	chii0 = chichos(n)
	chie0 = chiepale(n)

	write(121,1000) alphaion(n),-1.*chii0*((xkredge**2)+1./
	1				(exlti(n)**2)),omegthi(n),alphael(n),
     2 -1.*chie0*((xkredge**2)+1./(exlte(n)**2)),omegthe(n),
	3  -2.5*vrad1(n)*xnu/exlti(n)
774	continue       
	write(121,'(1x,35A)') '   rho     cosion    cosimp    sinion  
	1 sinimp	  eta-i    eta-e  '   
 
	do 775 n = 1,25
	write (121,1000) rhor(n),cosion(n),cosimp(n),sinion(n),sinimp(n),
     1  etai(n),etae(n) 	
     
775	continue 
	write(121,'(1x,35A)') '   rho     thetion     thetimp     G  
	1 erad6	   eradfb  '   
 
	do 780 n = 1,25
	write (121,1000) rhor(n),thetwid1(n),thetwid2(n),gy(n),
     1	erad6(n), eradfb(n)
    
     
780	continue 
	write(121,'(1x,35A)') '   rho      beam      ephi      drag  
	1 fric   vpinch2  '   

	
 
	do 785 n = 1,25
	ntorque = 1
c	 call edgerotran(n,ntorque)
    	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
	zdrag =	(xmas(1)*ynudrag1(n)*((erada(n)/bthet)+velthet1(n)/fp))/
	1	 (eq*bthet)
	zfrict = (xmas(1)*xnuc12(n)*(velthet1(n)-velthet2(n))/fp)/
     1	 (eq*bthet)
	vp2 = zbeam - ephia/bthet+zdrag+zfrict 
	write (121,1000) rhor(n),zbeam,-1.*ephia/bthet,zdrag, zfrict,vp2
 
    
     
785	continue 


	write(121,'(1x,35A)') '   rho    beam&ephi  phiterm   eradterm   v
	1theterm   vpinch5  ' 
	do 786 n = 1,25
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	eradterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)))/
	1			(eq*bthet)
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((velthet1(n)/fp)))/
     1			(eq*bthet)

 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp5 = zbeam-ephia/bthet+vphiterm+eradterm+vtheterm
     	write (121,1000) rhor(n),zbeam-1.*ephia/bthet,vphiterm,eradterm,
	1				vtheterm,vp5
786	continue 

      write(121,'(1x,35A)') '   rho    beam&ephi  vphiterm   vtheterm
     1vpinch0  ' 
	do 787 n = 1,25
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*
	1	(velthet1(n)-velthet2(n))/fp)/(eq*bthet)


 	vphiterm = 1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp0 = zbeam-ephia/bthet+vphiterm+vtheterm
     	write (121,1000) rhor(n),zbeam-1.*ephia/bthet,vphiterm,
	1				vtheterm,vp0
787	continue 

	
	write(121,'(1x,35A)') '   rho      Bvphi     -Bvthet  press  
	1  eradfb2  eradex  '   
 
	do 790 n = 1,25
	press(2) = -1.*ti(n)*xlpm(n)/zbar2(n)	
c	press(2) = bpcarb(n)/bthet
c	erfb2(n) = bthet*(vtor2(n) - vtheta(2)/fp + press(2))
c	eradfb2(n)=bthet*torv(n)+bfield*vthexp(n)+bthet*press(2)

	eradfb2(n)=bthet*vtor2(n)-bfield*velthet2(n)+bthet*press(2)
	write (121,1000) rhor(n),vtor2(n)*bthet,-1.*velthet2(n)*bfield,
     1	bthet*press(2),eradfb2(n),erex(n)
c   	write (121,1000) rhor(n),torv(n)*bthet,1.*vthexp(n)*bfield,
c     1	bthet*press(2),eradfb2(n),erex(n)


    
     
790	continue 
	write(121,'(1x,35A)') '   rho     vthet1     vthet2      vthexp  
	1eradnew  eradexp'
      do 795 n = 1,25
	write (121,1000) rhor(n),velthet1(n), velthet2(n),vthexp(n),
     1eradfbnew(n),erex(n)    
795	continue
	write(121,'(1x,35A)') '   rho     vtheory     erfb1      erad      
	1 eradfb2   erexp     erfb2  '
	do 800 n = 1,25
	write (121,1000) rhor(n),vtheory(n), erfb1(n),erada(n)
     1,	eradfb2(n),erex(n),erfb2(n)
800	continue 
	write (121,'(1x,35a)') '  omegi0    omegic    omegis    omefz0  
     1omegzc    omegzs  '
	do 805 n = 1,25 
      write (121,1000) omegt(n,1),omegt(n,2),omegt(n,3),omegt(n,4),
     1       omegt(n,5), omegt(n,6)		      
805	continue 
	WRITE (121,'(1X,35A)')'  VTH1      VTHTSIN  VTHTCOS  VTHTCOSZ   VT
     1HTSINZ    VTH2     vth2exp'
      DO 810, N=1,25
	WRITE (121,1000) Vpol(1,N),VTHTSIN(1,N),VTHTCOS(1,N),
     1   VTHTCOS(2,N),VTHTSIN(2,N),Vpol(2,N),vthexp(n)
810	CONTINUE 
    	write(121,'(1x,35a)')'  SV1      SV2         SV3       SV4    
     1 SV5      SV6    brack'
      DO 815, N=1,25
	WRITE(121,1000) S0V(1,N),S0V(2,N),S0V(3,N),S0V(4,N),S0V(5,N),
     1	S0V(6,N),brack(n)
815	CONTINUE
	WRITE(121,'(1X,35A)')'  rho      vphi1      vphi2   vphi2ex     vp
     1hi1ex, nustar-iz  nustar-zi'

	DO 820, N=1,25
		do 816 j=1,2
 	VTH(J) = SQRT(2.*XK*TI(N)/XMAS(J))
816	continue 
      xxxiz = xnuc12(n)*rmajor*qsafe/vth(1)		  	 
	xxxzi = xnuc21(n)*rmajor*qsafe/vth(2)	 

	WRITE(121,1000)rhor(n),vtor1(N),vtor2(N),torv(N),vphiex1(n), 
     1	xxxiz,xxxzi
     	
820	CONTINUE

      WRITE(121,'(1X,35A)')'  BV11      BV12   	 BV13      BV14
     1 BV15    bv16    '
	
      DO 825, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) B0V(1,1,N),B0V(1,2,N),B0V(1,3,N),B0V(1,4,N),
     1	B0V(1,5,n),b0v(1,6,n)
825	CONTINUE 
	 WRITE(121,'(1X,35A)')'  rho       vpinch   	 xlpm      ne     
     1no     omega '
     	goto 1559
	do 827 jj = 1,1
      do 826 nn = 1,25
	n = 26 - nn
c	if(ioptpinchi.eq.5) then
	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1)  + 
     1 xmas(1)*(ynudrag1(n)+xnuc12(n))*((erfb2(n)/bthet)+velthet1(n)/fp) 
     2	 - xmas(1)*xnuc12(n)*torv(n))/
     3	 (eq*atnum(1)*bthet)
	 
	diff(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
	1			((ynudrag1(n)/xnuc12(n))+1.)/
     2			((eq*atnum(1)*bthet)**2)
	xlpm(n) = (vrad1(n)-vpinchi(n))/diff(n) 
c	endif
	if(xlpm(n).lt.0.0) xlpm(n) = 0.0
	
	 xlnmold = xlnm(n) 
	xlnm(n) = xlpm(n)-c10*xltim(n) 
c	if(rhor(n).gt.pedrhoti) xlnm(n) = xlpm(n)-c10/ylti 
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpm(n)-c10/xltim(n) 
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpm(n)-c10/xltitop  
c	if(kk.gt.3) xlnm(n) = 0.5*(xlnm(n)+xlnmold)
	if(xlnm(n).lt.0.0) xlnm(n) = 0.0
c	integrate ion density inward from separatrix
	c6 = 0.5
	if(n.eq.25) goto 826		
	yni(n,1) = yni(n+1,1)*
	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
     2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
c	yni(n,1)=yni(n+1,1)*exp(delna*(c6*xlnm(n)+(1.-c6)*xlnm(n+1))) 
c	if(n.eq.24) yni(n,1)=yni(n+1,1)/(1.-delma*xlnm(n))
	if(yni(n,1).le.0.0) yni(n,1) = 1.e19 
	if(yni(n,1).gt.2.e20) yni(n,1) = 2.e20
	yni(n,2) = fracz*yni(n,1)
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	zne(n) = yne
826	continue
	call neutdist
827	continue
1559	continue			 
	DO 830, N=1,25

	toromega = torv(n)/rmajor 
	WRITE(121,1000) rhor(N),vpinchi(N),xlpm(N),zne(n),yno(n),toromega
830	CONTINUE 
 	
	WRITE(121,'(1X,35A)')'  aV11      aV12   	 aV13     av14   av15
	1	av16  '      
	
      DO 835, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) aV(1,1,N),aV(1,2,N),aV(1,3,N),av(1,4,n),av(1,5,n),
	1	av(1,6,n)
835	CONTINUE 
	WRITE(121,'(1X,35A)')'  sinion   cosion   vthtsin     vthtcos   
     1vphisin   vphicos  '	        
	
      DO 840, N=1,25
 	WRITE(121,1000) sinion(n),cosion(n),vthtsin(1,n),vthtcos(1,n),
	1	vphisin(1,n),vphicos(1,n)
840	CONTINUE 

	WRITE(121,'(1X,35A)')'  rho      torvel1    torvel2   vphi2ex   vp  
     1hi1ex   veltor1    veltor2'

	DO 845, N=1,25
		do 841 j=1,2
 	VTH(J) = SQRT(2.*XK*TI(N)/XMAS(J))
841	continue 
       

	WRITE(121,1000)rhor(n),torvel(N,1),torvel(N,2),torv(N),vphiex1(n), 
     1	veltor(n,1),veltor(n,2)
     	
845	CONTINUE
	WRITE(121,'(1X,35A)')'  rho      chiphi1    chiphi2   nudragp1     
     1nudragp2, nuinert1    nuinert2'

	DO 850, N=1,25
		do 846 j=1,2
 	VTH(J) = SQRT(2.*XK*TI(N)/XMAS(J))
846	continue 
       

	WRITE(121,1000)xlpm(n),chiphi(N,1),chiphi(N,2),xnudragp(N,1),
     1	xnudragp(n,2), ynuinert(n,1),ynuinert(n,2)
     	
850	CONTINUE
	WRITE(121,'(1X,35A)')'   RHO      A1        A2        SPOLA        
	1SPOLB        S2        S1  '
	DO 855 N=1,25
	WRITE(121,1000)RHOR(N),APOL1(N),APOL2(N),SPOLA(N,2),SPOLB(N,2),
     1	SPOL2(N),SPOL1(N)
855	CONTINUE      

	write(121,1005) chixpi,chixpe,chitop,chetop
 
	write(121,1006) ylti,ylte,xltitop,xltetop
1000	format(7e10.3)
1005	format(1x,'chixpi=',f5.2,1x,'chixpe=',f5.2,1x,'chitop=',f5.2,1x,
     1	'chetop=',f5.2) 
1006	format(1x,'ylti=',f6.3,1x,'ylte=',f6.3,1x,'xltitop=',f6.3,
	1	1x,'xltetop=',f6.3) 

	
	RETURN
	end










