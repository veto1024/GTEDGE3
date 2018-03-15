	SUBROUTINE EDGECALC1.for
	INCLUDE 'SOLDIV.FI'
	parameter (ml=26,jq=2)
	dimension ynold(ml),yniold(ml),
     1			velthet1old(ml),ratdrag(ml),beamdot(ml),
	2			errthet(ml),errtor(ml),cxcool(ml),coolion(ml),
	3			radcool(ml),zbar(ml),qie(ml),tiold(ml),
	4			delvelthet(ml),vphizterm(ml),erterm(ml),
	5			velthet2old(ml),ererad(ml),eradaold(ml),
	6			atcool(ml),gamconde(ml),gamcondi(ml),teold(ml),
	7			zne(ml)

	double precision TEDP,XLZ2DP,DLDT2DP,Tdbl,xlzdbl,dlzdbl
c	bfield = bphi
C	CALCULATES NEUTRAL & ION DENSITY, TEMP AND VELOCITY DISTRIBUTIONS IN EDGE
C	SET TOROIDAL ROTATION BC	
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
	chiion(25) = chixpi
	chiel(25) = chixpe
 	DO 5, NN=1,24
	N = 25-NN
	rhor(n) = rhor(n+1) - delna/(aminor*SQRT(0.5*(1.+ELONG**2)))
	chiion(n)=chixpi
	if(rhor(n).lt.pedrhoti) chiion(n)=chitop 
	chiel(n)=chixpe
	if(rhor(n).lt.pedrhote) chiel(n)=chetop
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
	

c	********************************
	
c	"Inner" iteration converging ion and neutral distributions through 200
c	********************************8	 
 	DO 200 IT = 1,100
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

	IZ1 = atnum(2)
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
	
	temp(1) = ti(25)
	temp(2) = ti(25)
	tep     = tel(25)
	call param(25)
	xnuc12(25) = xnuc(1,2)
	xnuc21(25) = xnuc(2,1)
	delma = delna
c*****note******** 
	ntorque = 0

	call edgerotran(25,ntorque)
    	thetw1(25) = thetw(1)
	thetw2(25) = thetw(2) 
	ratdrag(25) = xnudragatomic(25)/ynudrag1(25)

      
c	*********************	
      cncmult = 2.0 
	cetaimult = 1.0
	cetgmult = 1.0
	cedwmult = 1.0
	radmultedge = 1.0
c******************************
	
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
	IZ1 = atnum(2)
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
	CALL TOROTATE(N)
	DO 1055 I = 1,6
	DO 1050 J = 1,6
	SV(I) = SV(I) - AV(I,J)*OMEGT(N+1,J)/DELMA
	BV(I,J) = BV(I,J) - AV(I,J)/DELMA
1050	CONTINUE 
1055  CONTINUE
	CALL LSLRG(6,BV,6,SV,1,SV)
	DO 1060 I = 1,6
	OMEGT(N,I) = SV(I)
1060	CONTINUE	
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
	iF(MIT.EQ.0) GOTO 250


	if(iconverge.eq.0.and.it.eq.15) then
	write (6,199)
	jwarn = 1
	goto 740
	endif
199	format(1x,'not converged on 200 loop') 
200	CONTINUE


250	continue	
c	converge on poloidal & toroidal velocities ***************turned off*******
c	ioptshear = 1 
c	if(n.gt.15) sheare = sheareal
c************convergence on velocities & Er turned off--trouble converging 3/18/04

c	goto 374
c*********************************************************************************
	do 275 n = 1,25
	errthet(n) = abs((velthet1old(n)/velthet1(n))-1.)
	if(errthet(n).gt.0.05) mjt = mjt + 1
	errthet(n) = abs((velthet2old(n)/velthet2(n))-1.)
 	if(errthet(n).gt.0.05) mjt = mjt + 1 
	ererad(n) = abs(eradaold(n)/erada(n))
c	if(ererad(n).gt.0.05) mjt = mjt + 1 
	velthet1(n) = 0.5*(velthet1(n)+velthet1old(n))
	velthet2(n) = 0.5*(velthet2(n)+velthet2old(n))
	
	errtor(n) = abs((vtor1old(n)/vphiex1(n))-1.)
c	if(errtor(n).gt.0.05) mjt = mjt + 1
	
275	continue

	if(mjt.ne.0) goto 298
	 
	goto 350

298	do 299 n=1,25	
      vtor1old(n) = vphiex1(n)
299	continue
			
300	continue		 
	
350	continue
c	goto 700
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
	if(xlpm(n).lt.0.0) xlpm(n) = 0.0
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
	do 675 nn = 1,24
	n = 25-nn
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
	IZ1 = atnum(2)


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
	DO 1155 I = 1,6
	DO 1150 J = 1,6
	SV(I) = SV(I) - AV(I,J)*OMEGT(N+1,J)/DELMA
	BV(I,J) = BV(I,J) - AV(I,J)/DELMA
1150	CONTINUE 
1155  CONTINUE
	CALL LSLRG(6,BV,6,SV,1,SV)
	DO 1160 I = 1,6
	OMEGT(N,I) = SV(I)
1160	CONTINUE	     
 
	goto 675			  
C	NEOCLASSICAL CHI FOR IONS
c	simple neoclassical chi 
	ep = aminor*SQRT((1.+ELONG**2)/2.)/rmajor 
	bfield = abs(bphi) 
	OMI =EQ*BFIELD/XMAS(1)
	CSOUND = SQRT(XK*TEL(N)/XMAS(1))
 	rhot = csound/omi
 
	CHINC(n) = ((RHOTi(1)*bfield/bthet)**2)*XNUc(1,2)*(EP**0.5) 
C	CHANG-HINTON CHI
	ALFA = XNi(2)*(zbar2(n)**2)/(XNi(1)*(ATNUM(1)**2))
	qa = ep*bfield/bthet
	XMUii =(XNUc(1,1)*Q95*RMAJOR/(vth(1)*(EP**1.5)))*(1.+1.54*alfa)
	
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
	
625	chinc(n) = chinc(n)/(abs(sheare(n))**1.5)
	chichos(n) = chich/(abs(sheare(n))**1.5)
	

C	ITG-MODE CHI FOR IONS
	
	ETAI = XLTIM(N)/XLNM(N)
	
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
	ETAE =  XLTEM(N)/XLNM(N)
c	****exp etae**************
	etae = 1.43
c	**************************
	CSE = SQRT(2.*XK*TEL(N)/EMASS)
	eplasfreq = 56.4*sqrt(YNI(N,1))
	CHIEETG(N)=0.13*((CLIGHT/EPLASFREQ)**2)*CSE*SHEARM*ETAE*(1.+ETAE)/
	1			(Q95*RMAJOR)
C	TRAPPED ELECTRON MODE W/INTERPOLATION TO COLLISIONLESS DRIFT MODE
C		CHI FOR ELECTRONS 
	CSE = SQRT(2.*XK*TEL(N)/EMASS)
	RHOS = CSOUND/OMII(1)
	RHOTE = 3.37E-6*SQRT(TEL(N)/BFIELD)
	CHIEDW(N) = 2.5*(EP**1.5)*(CSOUND**2)*(RHOS**2)*XLNM(N)*XLTEM(N)/
	1	(XNUEI*(1.+0.1/XNUEISTAR))
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
	OPEN(122,FILE='pedestalplot.TXT',STATUS='UNKNOWN') 
	if(jwarn.eq.1) write(121,199)

	write(121,'(1x,35A)') '   rho	     ne         ti         te     
     1  erad    OMEGT01'
	write(122,'(1x,35A)') '   rho	     ne         ti         te  ' 
	sepdif = xnsepex - xnesepreal
     	do 750 n = 1,25
	zne(n) = zne(n) - sepdif
	write (121,1000) rhor(n),zne(n),ti(n),tel(n),erada(n),OMEGT(n,1)
	write (122,1000) rhor(n),zne(n),ti(n),tel(n)
	sion(n) = zne(n)*0.5*(ynuioni(n) + ynuioni(n+1))
c	if(n.ne.25) xlnm(n+1) = log(yni(n,1)/yni(n+1,1))/delna
750	continue
	
	write(121,'(1x,35A)') '   rho	   gamheate    gamheati   gamion     
     1  vrad    vpinchi'
     	do 755 n = 1,25
	write (121,1000) rhor(n),gamheate(n),gamheati(n),gamion(n,1),
     1vrad1(n),vpinchi(n)
755	continue 
      write(121,'(1x,35A)') '   rho	    nudrag*    nuati    nuioni       
     1nuvsic1  nuvisc2  '
     	do 760 n = 1,25
	write (121,1000) rhor(n),ynudrag1(n),xnuati(n),xnuioni(n),
     1xnudragvis1(n),xnudragvis2(n)
760	continue 
	write(121,'(1x,35A)') '   rho	     xno      coldno    cxcool       
     1radcool  atcool'
 
	do 765 n = 1,25
	write (121,1000) rhor(n),yno(n),coldno(n),cxcool(n),
     1radcool(n),atcool(n)
765	continue 
	write(121,'(1x,35A)') '   rho	     lpm        ltim       lnm    
     1 ltem	   lvm'
   
 
	do 770 n = 1,25
	write (121,1000) rhor(n),xlpm(n),xltim(n),xlnm(n),xltem(n),xlvm(n)
     
          
770	continue 

	write(121,'(1x,35A)') '  rho     chiion     chiel    diff' 
	do 771 n = 1,25
	write(121,1000) rhor(n),chiion(n),chiel(n),diff(n)
771	continue 


	write(121,'(1x,35A)') '   rho     cosion    cosimp    sinion  
	1 sinimp	  diff  '   
 
	do 775 n = 1,25
	write (121,1000) rhor(n),cosion(n),cosimp(n),sinion(n),sinimp(n),
     1  diff(n) 	
     
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


	write(121,'(1x,35A)') '   rho      beam      ephi      vphiterm  
	1dragfric vpinch5  ' 
	do 786 n = 1,25
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	zdragfric =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)+
     1	velthet1(n)/fp))/(eq*bthet)
 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp5 = zbeam-ephia/bthet+vphiterm+zdragfric
     	write (121,1000) rhor(n),zbeam,-1.*ephia/bthet,vphiterm,zdragfric
	1,vp5
786	continue 

      
	
	write(121,'(1x,35A)') '   rho      vphi     -vthet/fp  -ti*lpm/B  
	1  eradfb2/b eradfb1/b  '   
 
	do 790 n = 1,25
	write (121,1000) rhor(n),torv(n),-1.*velthet2(n)/fp,
     1	-1.*ti(n)*xlpm(n)/(zbar2(n)*bthet),eradfb2(n)/bthet,
     2	eradfb1(n)/bthet

    
     
790	continue 
	write(121,'(1x,35A)') '   rho     vthet1     vthet2      vthexp  
	1eradnew  eradexp'
      do 795 n = 1,25
	write (121,1000) rhor(n),velthet1(n), velthet2(n),vthexp(n),
     1eradfbnew(n),erex(n)    
795	continue
	write(121,'(1x,35A)') '   rho     vtheory     vthet2     erexp      
	1 eradfb2   erad  '
	do 800 n = 1,25
	write (121,1000) rhor(n),vtheory(n), velthet2(n),erex(n)
     1,	eradfb2(n),erada(n)
800	continue 
	write (121,'(1x,35a)') '  omegi0    omegic    omegis    omefz0  
     1omegzc    omegzs  '
	do 805 n = 1,25 
      write (121,1000) omegt(n,1),omegt(n,2),omegt(n,3),omegt(n,4),
     1       omegt(n,5), omegt(n,6)		      
805	continue 
	WRITE (121,'(1X,35A)')'  VTH1      VTHTSIN  VTHTCOS  VTHTCOSZ   VT
     1HTSINZ  VTH2   '
      DO 810, N=1,25
	WRITE (121,1000) Vpol(1,N),VTHTSIN(1,N),VTHTCOS(1,N),
     1   VTHTCOS(2,N),VTHTSIN(2,N),Vpol(2,N)
810	CONTINUE 
    	write(121,'(1x,35a)')'  SV1      SV2         SV3       SV4    
     1 SV5      SV6    '
      DO 815, N=1,25
	WRITE(121,1000) SV1(N),SV2(N),SV3(N),SV4(N),SV5(N),SV6(N)
815	CONTINUE
	WRITE(121,'(1X,35A)')'  xlncm1    xlnsm1   xlncm2    xlnsm2   nu
     1star-iz   nustar-zi'

	DO 820, N=1,25
		do 816 j=1,2
 	VTH(J) = SQRT(2.*XK*TI(N)/XMAS(J))
816	continue 
      xxxiz = xnuc12(n)*rmajor*qsafe/vth(1)		  	 
	xxxzi = xnuc21(n)*rmajor*qsafe/vth(2)	 

	WRITE(121,1000)xlncm1(N),xlnsm1(N),xlncm2(N),xlnsm2(N),xxxiz,xxxzi
     	
820	CONTINUE

      WRITE(121,'(1X,35A)')'  BV11      BV22   	 BV23      BV24
     1 BV25    bv33   '
	
      DO 825, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) BV11(N),BV22(N),BV23(N),BV24(N),BV25(N),bv33(n)
825	CONTINUE 
	 WRITE(121,'(1X,35A)')'  BV12      BV13   	 BV14      BV32
     1 BV21    bv31   '
	
      DO 830, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) BV12(N),BV13(N),BV14(N),BV32(N),BV21(N),bv31(n)
830	CONTINUE 
 	
	WRITE(121,'(1X,35A)')'  aV11      aV12   	 aV13     cos0  sino'      
	
      DO 835, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) aV11(N),aV12(N),aV13(N),cosdeno(1),sindeno(2)
835	CONTINUE 
	WRITE(121,'(1X,35A)')'  sinion   vthtsin cosion+1    vthtcos
     1	  '      
	
      DO 840, N=1,25
 	WRITE(121,1000) sinion(n),vthtsin(1,n),cosion(n)+1.,vthtcos(1,n)
840	CONTINUE 



	write(121,1005) chixpi,chixpe,chitop,chetop
 
	write(121,1006) ylti,ylte,xltitop,xltetop
1000	format(6e10.3)
1005	format(1x,'chixpi=',f5.2,1x,'chixpe=',f5.2,1x,'chitop=',f5.2,1x,
     1	'chetop=',f5.2) 
1006	format(1x,'ylti=',f6.3,1x,'ylte=',f6.3,1x,'xltitop=',f6.3,
	1	1x,'xltetop=',f6.3) 

	y = qnbi(n)
	y = qnbe(n)
	y = f(1)
	y = f(2) 
	y = eradfb(n)
c	y = thetwid(n)
	y = sheare(n)
	y = rhot1(n)
	Y  = FLUXINNEUT(N)
	y = cosden(1)
	Y  = BEAMSUM(N)
	Y = PRESSUM(N)
	Y = VPOLSUM(N)
	y = sinpot
	y = cospot
	Y = DRAGSUM(N)
	Y = ERADA(N)
	Y = YNUDRAG1(N)
	Y = YNUDRAG2(N)
	y = pinert(n)
	y = sinden(1)
	y = chieetg(n)
	y = vrad1(n)
	y = vrad2(n)
	y = vpin(1,n)
	y = dif(n,1,2)
	y = acoef(n,1)
	y = ynuati(n)
	y = xnuionb(n)
	y = vtor1(n)
	y = vtor2(n)
	y = vpinchi(n)
	y = diff(n)
	y = coefp(n)
	y = ratnu(n)
	y=xmomtor1(n)
	y=xmomtor2(n)

	y = press(1)
	y = vphivth1
	y = xfrac
	y = fricmom(n)
	y = tormom(n)
	y = eph(m)
	y = dragmom(n)
	y = xnudraggyro(n)
	y = xmtor(n)
	y = gamouteff
	nnn=niter
	y = xnsep
	y= press1(n)
	y= press2(n)
	y = xnsol
	y = ynuioni(n)
	y = torv(n)
	y = vphiex1(n)
	y = vphiex2(n)
	x = ioptvdif
	x = ioptvphi
	x = ioptvthet
	y = ioptpinch
	y = elong
	x = ioptxvlm
	x = ioptapproach
	y = sinion(n)
	y = cosion(n)
 	y = phimom(n)
	y = omegt(n,1)
	y = omegt(n,2)

	RETURN

	END