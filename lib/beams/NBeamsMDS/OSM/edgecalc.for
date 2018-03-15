	SUBROUTINE EDGECALC
	INCLUDE 'SOLDIV.FI'
	parameter (ml=26,jq=2,jp=6)
	dimension ynold(ml),yniold(ml),	xnuatim(ml),zkHS1(ml),zkHS2(ml),
     1			velthet1old(ml),ratdrag(ml),beamdot(ml),
	2			errthet(ml),errtor1(ml),cxcool(ml),coolion(ml),
	3			radcool(ml),zbar(ml),qie(ml),tiold(ml),s0v(jp,ml),
	4			delvelthet(ml),vphizterm(ml),erterm(ml),
	5			velthet2old(ml),ererad(ml),eradaold(ml),epota(ml),
	6			atcool(ml),gamconde(ml),gamcondi(ml),teold(ml),
	7			vtor2old(ml),pressd(ml),veld(ml),ssv(jp),bbv(jp,jp),
	8			b0v(jp,jp,ml),CHITHI(ML),CHITHE(ML),chich8(ml),
	9			vphical(jq),ynud(jq),errtor2(ml),bnud1(ml),bnud2(ml),
	1			exd21(ml),exd22(ml),exd23(ml),exd51(ml),exd52(ml),
     2			exd53(ml),vpinch2(ml),vpinch5(ml),ssion(ml),ypi(ml),
     3			vpinch3(ml),exd31(ml),exd11(ml),exd41(ml),vpinch4(ml),
	4			heatvisc(ml),heatin(ml),eta0(ml),eta4(ml),qcondi(ml),	 
	5			znec(ml),rhorn(ml),qconde(ml),vth1cor(ml),
     6			eradnm(ml,ml),Eminx(ml,ml),Xtranp(ml,ml), 
	7			Xtrane(ml,ml),xlossn(ml),xlosse(ml),Wminx(ml),
	8			xnsource(ml),sinknx(ml),xesource(ml),gamionorb(ml),
	9			sinkex(ml),Wminx2(ml),dtheta(ml),dthet(ml),E1(ml,ml),
	1			seplossn(ml),seplosse(ml),sepnx(ml),sepex(ml),
	2			fraclossn(ml,ml),fraclosse(ml,ml),Exloss(ml),FX(ml),
	3			EX(ml),gamionorbX(ml),gamheatiorb(ml),
     4			gamheatiorbX(ml),xchiiorbX(ml),xchiiorbx15(ml),
     5			XXN(ML,ML),XXE(ML,ML),GAMIONx(mL), gamheatix(ll),
	6			xXCHII(ml),xXCHII15(ml),xXCHIIorb(ml),xXCHIIorb15(ml),
	7			qcond25(ml),qcond15(ml),qcond25orb(ml),qcond15orb(ml),
	8			qcond25orbx(ml),qcond15orbx(ml),qcond25xtran(ml),
	9			qcond15xtran(ml),qcond25xtranorb(ml),xgamionorb(ml),
     1			qcond15xtranorb(ml),xheatXtranorb(ml),gamin(ml),
	2			gamout(ml),xtransink(ml),xtransource(ml),cdxtran(ml), 
     3		cdionorb(ml),xebin(ml),xebout(ml),velthet3(ml),xnet(ml),
	4			gxg(ml,ml),dxdrg(ml,ml),xchieorb(ml),xnuc120(ml),
	5            snion(ml), xLnmin(ml),xmomiol(ml),torq1(ml),torq2(ml),
     6			xddrag1(ml),xddrag2(ml),XNUDZERO1(ML),XNUDZERO2(ML),
     7			DELV1(ML),DELV0(ML),xnudtot1s(ml),xnudtot2s(ml),
     8			GAMIONORBCUR(ML), zkHS(jq), vpol61(ml),vpol62(ml),
     9		drive1(ml),drive2(ml) ,drive1v(ml),drive2v(ml),fop(ml),
	2			xlpmom(ml)
  
	DIMENSION EMIN(ML,ML),EPMIN(ML,ML)

	double precision TEDP,XLZ2DP,DLDT2DP,Tdbl,xlzdbl,dlzdbl
	real NBIreturn
c************************************************************************
c	ipolopt = 1 uses vthet2 = vexp, =0 calculates vthet2
	ipolopt = 1
c************************************************************************
	
	SE = SEGEOM
c**************INTERPRETIVE MODE GO TO 666**************************	
	goto 666  
C*******************************************************************

c	bfield = bphi
C	CALCULATES NEUTRAL & ION DENSITY, TEMP AND VELOCITY DISTRIBUTIONS IN EDGE
C	SET TOROIDAL ROTATION BC	
	
	jdrag = 1
c	jdrag=1 causes xnud to be inferred from exp vphi

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
c	use ioptapproach=0 for inferring chi
	ioptapproach = 0
c	ioptped = 0, do pedestal profile calc; = 1, use input pedestal profiles
	ioptped = 0
	if(ioptped.eq.1) goto 1616 
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
	rhorn(25) = 1.0
	if(ioptran.eq.1) then
	chiion(25) = chixpi
	chiel(25) = chixpe
	endif
 	DO 5, NN=1,24
	N = 25-NN
	rhor(n) = rhor(n+1) - delna/(aminor*SE)
	rhorn(n) = rhorn(n+1) - delna/aminor
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
	if(ioptxlvm.eq.3) xlvm(n) = -2.*(vphiex2(n)-vphiex2(n-1))/
	1		(delna*(vphiex2(n)+vphiex2(n-1)))
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
	radmultedge = qzmultcore
	radmultedge = 1.
c******************************

	xs = xnav

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
c	if(it.gt.1) diffA(25) = 0.5*(diffA(24)+diffA(25)) 
	xlpm(25) = (vrad1(25)-vpinchi(25))/diffA(25)
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
	
	GAMION(N,1) = GAMION(N+1,1) -dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n)) 
45	continue
c	********************************
c	new ion & plasma density and ion & electron temp formulation
	dens(n) = 0.5*(yni(n+1,1)+yni(n,1) )
	tele = 0.5*(tel(n)+tel(n+1))
	tiav = 0.5*(ti(n)+ti(n+1))
	rz = (atnum(1)**2/xmas(1) + zbar2(n)**2/xmas(2))*1.67e-27  
	cequil = 6.32e-14*rz 
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)
c	eq 4.90 fpp
	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(yne) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(tele**1.5)*X/Yz)

	cequil = 7.9e-42*couloge*zeff/xmas(1) 

	qie(n) = cequil*yne*(tiav-tele)/(tele**1.5)

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
c	zbar2(n) = 6.0
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
	GAMION(N,1) = GAMION(N+1,1)-dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n)) 
 	GAMION(N,2) = GAMION(N+1,2) 
     1		-	DENS(N)*XNUIONz(N)*DELMA*(1.+fracz*zbar2(n))
c	fraction of heating to ions
	fb1 = 0.75
 	fb2 = 0.15
	fb3 = 0.10

	ecrit = 19.*tel(n)
	xc= sqrt(1.e3*eb/ecrit)
	xx = atan((2.*xc-1.)/1.732)
	yy = log((xc**2+2.*xc+1.)/(xc**2-xc+1.))
	fion1 = 2.*(((xx+.5236)/1.732)-yy/6.)/(xc**2) 
	xc= sqrt(1.e3*(eb/2.)/ecrit)
	xx = atan((2.*xc-1.)/1.732)
	yy = log((xc**2+2.*xc+1.)/(xc**2-xc+1.))
	fion2 = 2.*(((xx+.5236)/1.732)-yy/6.)/(xc**2)
	xc= sqrt(1.e3*(eb/3.)/ecrit)
 	xx = atan((2.*xc-1.)/1.732)
	yy = log((xc**2+2.*xc+1.)/(xc**2-xc+1.))
	fion3 = 2.*(((xx+.5236)/1.732)-yy/6.)/(xc**2) 
	fionb(n) = fb1*fion1 + fb2*fion2 + fb3*fion3
	qnbi(n) = fionb(n)*qnb(n)
	qnbe(n) = (1.-fionb(n))*qnb(n) 
	

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

	xlpm(n) = (vrad1(n)-vpinchi(n))/diffA(n)	
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
c	diffA(n) = diff%(n)/(sheare(n)**1.5)
	xlpm(n) = (vrad1(n)-vpinchi(n))/(diffA(n)) 
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
 
	diffA(25) =  xmas(1)*xk*ti(25)*xnuc12(25)*
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
c	diffA(n) =  xmas(1)*xk*ti(25)*xnuc12(25)*
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
	GAMION(N,1) = GAMION(N+1,1) -dens(n)*xnuionb(n)*delna	
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
c	diffa(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
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
c	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
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
c	zbar2(n) = 6.0
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
	cequil = 6.32e-14*rz 
c	eq 4.90 fpp
	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(yne) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(tele**1.5)*X/Yz)
 
	cequil = 7.9e-42*couloge*zeff/xmas(1) 
	
	qie(n) = cequil*yne*(tiav-tele)/(tele**1.5)
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
	1(1.+delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))/
	2(1.-delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))
	ti(n)=0.5*(ti(n)+tiold(n))
c		if(n.eq.24) ti(n) = ti(n+1)/(1.-delma*xltim(n))		
	if(ti(n).lt.0.0) ti(n) = 100.
	tel(n) = tel(n+1)*
	1(1.+delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))/
	2(1.-delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))
	tel(n)=0.5*(tel(n)+teold(n)) 
c		if(n.eq.24) te(n) = te(n+1)/(1.-delma*xltem(n))		
      if(tel(n).lt.0.0) tel(n) = 100.
C		INTEGRATE TOROIDAL ROTATION FREQ INWARD FROM SEPARATRIX
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
	ep = aminor*SE 
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
	CHICH(n) = CHINC(n)*(xnuc(1,1)/xnuc(1,2))*(A1*G1+A2*(G1-G2))
	if(dp.eq.0) then
		chich0 = chich(n)
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
	chichos(n) = chich(n)/(abs(sheare(n))**1.5)
	

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
c	GOTO 1610 
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
1467	AM = (aminor*SE)
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
c	nbi momentum input & iol momentum input
	  
	xmtor(1) = xmomtor1(j)
c	IOL momentum input
c*****yy1 is filling in for the iol delta-V term*****************
	xmomiol(j) = yni(j,1)*xmas(1)*yy1(j) 
	xmtor(1) = xmtor(1) + xmomiol(j)  
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
c	goto 7200 
1616	continue 
c	*************************************************************************
	
C	INFER CHI's FROM EXP N,T,GSCL & CALC. Q, GAM
C	DO 1650 N= 1,25
C	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
C	XCHIE(N) = EXLTE(N)*((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
C     1					2.5*GAMEL/EXNE(N))	 
C	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))
C	XCHII(N) = EXLTI(N)*((GAMHEATI(N)/(XNION*XK*XTI(N)))-
C    1					2.5*GAMION(N,1)/XNION)
C1650	CONTINUE 


C	****************************************************************************
C	INTERPRETIVE MODE, USES EXP. Ne, Te, Ti TO CALCULATE HEAT & PARTICLE FLUXES
C	IN ORDER TO INFER CHI'S.  THEN CALCULATES ROTATION VELOCITIES.    	 
C	******************************************************************************
c	goto 3700
c	new 12/5/05******initialize to input exp data**************

666	continue
c	*****************compute potential profile from Erad_exp***********************
      se = SQRT(0.5*(1.+ELONG**2))
	epot(25) = 600.0 
	do 667 nn = 1,24
	n=25-nn
	epot(n) = epot(n+1) + 0.5*(delna/se)*(erex(n+1)+erex(n))
667	continue	  
	do 668 n= 1,25
	xlphim(n) = erex(n)/epot(n)
	epota(n) = epot(n)
	epot(n) = epot(n)/xte(n)
668	continue
c	***********************************************************************************

	QSAFE = Q95
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
C	set plasma n & T & gsl to input experimental values	  
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
	XNI(1) = YNI(N,1)
	XNI(2) = YNI(N,2)
	TEMP(1) = TI(N)
	TEMP(2) = TI(N)
	TEP = TEl(N)
      CALL PARAM(N)
	xnuc12(n) = xnuc(1,2)
	xnuc21(n) = xnuc(2,1)
	xnuc11(n) = xnuc(1,1)
	xnuc22(n) = xnuc(2,2)
c	xnudrag(1) = xnudtot1(n)
c	diffa(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
c    1			((xnudrag(1)/xnuc12(n))+1.-atnum(1)/zbar2(n))/
c   2			((eq*atnum(1)*bthet)**2) 		
1350	continue
	SE = sqrt((1.+elong**2)/2)
	rhor(25) = 1.0 
	rhorn(25) = 0.0	 
	DO 1355, NN=1,24
 	N = 25-NN
	rhor(n) = rhor(n+1) - delna/(aminor*SE)
	rhorn(n) = rhorn(n+1) - delna/(se)
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
	xnuatim(n) = yno(n)*svata(n)
		EIONe = 17.5
	IF(ynd.LE.1.E21) 
     2     EIONe = 17.5 + (5.+37.5/Ted)*LOG10(1.E21/ynd)
    	IF(ynd.GT.1.E21)
     2    EIONe = (30.6 - 16.4*EXP(-5.E19/ynd))*
     3              EXP(5.45/(Ted*EXP((ynd/1.37E20)**0.26)))	
	alphael(n) = xnuioni(n)*(2.5*(xnu-1.)+xnu*eione/ted)
	alphaion(n) = 1.5*xnuati(n)*(xnu-1.)
1360	continue

C***************X-LOSS AND X-TRANSPORT CALCULATION***************	
	
c***************************************************************************
C	X-Loss and X-Transport
c**************************************************************************
	psix = 0.5
      rhotopx = 1. - 0.5*deltathetx  
	
c	X -transport sinks and sources
c	sinks
	do 7989 n = 1,25
	sinknx(n) = 0.0
	sinkex(n) = 0.0

	if(rhor(n).gt.rhotopx) goto 7986
		Wminx(n) = 1.e5
		nx = n+1
	goto 7989
7986	delr = delna*SQRT(0.5*(1.+ELONG**2))
	se = SQRT(0.5*(1.+ELONG**2))
 	delrho = delna/aminor*se

	xniexp = exne(n)/(atnum(1)+fracz*zbar(n))
	rr = rhor(n)*aminor*SQRT(0.5*(1.+ELONG**2))
c	Minimum W to gradB drift delrho/2 out of node n
	Wminx2(n)=(rmajor*abs(erex(n))/deltathetx)*(delrho/rhor(n))/
     2			(1.+psix**2)	 
	A = 1.5
	X = WMINX(n)/xti(n)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	fns = (value1 - value)/value1
	sinknx(n) = xniexp*abs(erex(n)/bphi)*fns/(6.28*rr)
	A = 2.5
	EVALUE = GAMI(A,X)
	EVALUE1= GAMMA(A) 
	ens = (value1 - value)/value1
	sinkex(n) =  xk*xti(n)*xniexp*abs(erex(n)/bphi)*ens/(6.28*rr)

7989	continue

      
796	continue 
c	Erad integral	
c	do 798 n = nx,25
c	erint = erex(n)
c	do 797 m = n+1,25
c	erint = erint + erex(m)
c	eradnm(n,m) = erint/(m-n)
c	Eminx(n,m) = 2.*(1.- (rhor(n)/rhor(m)))*rmajor*abs(eradnm(n,m))/
c	1             ((1.+psix**2)*deltathetx*xti(n))
c	Minimum energy to gradB drift across separatrix
c	Wminx(n) = eminx(n,25)*xti(n) 
c	A = 1.5
c	X = EMINX(n,m)
c	VALUE = GAMI(A,X)
c	VALUE1= GAMMA(A)
c	Xtranp(n,m) = (value1 - value)/value1

c	A = 2.5
c	EVALUE = GAMI(A,X)
c	EVALUE1= GAMMA(A)
c	Xtrane(n,m) = (evalue1 - evalue)/evalue1

c797	continue 
c	xlossn(n) = xtranp(n,25)
c	xlosse(n) = xtrane(n,25)
c798	continue
	delrho = delna*SQRT(0.5*(1.+ELONG**2))
 
	Eminx(25,25) = (delrho/rhor(25))*rmajor*abs(erex(25))/
	2				((1.+psix**2)*deltathetx*xti(25))
	Wminx(25) = Eminx(25,25)*xti(25) 
 	A = 1.5
	X = EMINX(25,25)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	Xtranp(25,25) = (value1 - value)/value1

	A = 2.5
	EVALUE = GAMI(A,X)
	EVALUE1= GAMMA(A)
	Xtrane(25,25) = (evalue1 - evalue)/evalue1

	xlossn(25) = xtranp(25,25)
	xlosse(25) = xtrane(25,25)

c	X-Loss across separatrix
	do 7987 n = 1,25
	xniexp = exne(n)/(atnum(1)+fracz*zbar(n))
	rr = rhor(n)*aminor*SQRT(0.5*(1.+ELONG**2))
	partloss = (xniexp*abs(erex(n))/(abs(bphi)))*
     2			xlossn(n)*delna
	energyloss = xti(n)*(xniexp*abs(erex(n))/(abs(bphi)))*
     2			xlosse(n)*delna*xk
	xlossion = xlossion + partloss
	xlosspow = xlosspow + energyloss
7987	continue

c	sources
	do 7991 n = 2, 25
	xnsource(n) = 0.0
	xesource(n) = 0.0
	do 7990 m = 1,n-1
	xniexp = exne(m)/(atnum(1)+fracz*zbar(m))
	rr = rhor(m)*aminor*SQRT(0.5*(1.+ELONG**2)) 

	xnsource(n) = xnsource(n) + 
     1			(xniexp*abs(erex(m))/(abs(bphi)))*
     2			(xtranp(m,n)-xtranp(m,n-1))/(6.28*rr)
	xesource(n) = xesource(n) + 
	1				xk*xti(m)*(xniexp*abs(erex(m))/(abs(bphi)))*
     2			(xtrane(m,n)-xtrane(m,n-1))/(6.28*rr)
7990  continue
7991	continue	

c**************************end X-loss*********************************************
c**************tally the sum X-Loss treatment**************************************

c*******calculate cumulative dthet ExB drift for successive radial increments*************************
	OPEN(129,FILE='XLOSS.TXT',STATUS='UNKNOWN') 
c	write(129,6601)
c	delr is radial mesh increment in effective circular model, aminor = 0.6m is actual minor radius
	delr = delna
	se = SQRT(0.5*(1.+ELONG**2))
	delrho = delna/aminor*se

	 
	W = 0.0
	do 6600 nw = 1,200
	W = W + 0.1e3
	write (129,6603) W
      do 6500 n = nx,24
	dthet(n) = rmajor*erex(n)*(delrho/rhor(n))/(W*(1.+psix**2))
      do 6400 m = n+1,24
	dthet(m) = dthet(m-1) +
	2			 rmajor*erex(m)*(delrho/rhor(m))/(W*(1.+psix**2))
6400	continue
c	write(129,6604) n
	write(129,6602) n,(dthet(m),m=nx,24)
	do 6450 m = 1,25
	dthet(m) = 0.
6450  continue	 
6500	continue 
6600	continue
c6601	format(129, '(1x,35A)') '  17    18    19    20    21    22   23
c     2 24   25   '
6602	format(1x,I2,1x,20f8.3)
6603  format(1x,"W =",e8.3)
6604  format(1x,"launch node =",I2)

c	sources and sinks
c	input from 6400 calculation
	do 6800 n = 1,25
	do 6800 m = 1,25
	E1(n,m) = 0.
6800	continue

c******************************************done*********************************************************
c	these dthet must now be examined visually from xloss.txt to determine the minimum energies required
c	to gradB drift from n to m radial mesh points, which are then entered by hand below
c******************************************************************************************************* 


c	min W for entry in mesh n to gradB drift out mesh m	***************

c	goto 23302
	goto 9876 
c	shot 11897@1525 L-mode, deltathetx = 0.15 
	
     
   
	  EMIN(12,13) = 0.5
	  EMIN(12,14) = 0.9 
		EMIN(12,15) = 1.3
		EMIN(12,16) = 1.7
		EMIN(12,17) = 2.2
		EMIN(12,18) = 2.6
		EMIN(12,19) = 3.1
		EMIN(12,20) = 3.7
		EMIN(12,21) = 4.2
		EMIN(12,22) = 4.8 
		EMIN(12,23) = 5.5
		EMIN(12,24) = 6.2
		EMIN(12,25) = 6.9
		
		EMIN(13,14) = 0.5 
		EMIN(13,15) = 0.9
		EMIN(13,16) = 1.3
		EMIN(13,17) = 1.8
		EMIN(13,18) = 2.2
		EMIN(13,19) = 2.7
		EMIN(13,20) = 3.3
		EMIN(13,21) = 3.8
		EMIN(13,22) = 4.4 
		EMIN(13,23) = 5.1
		EMIN(13,24) = 5.8
		EMIN(13,25) = 6.5

		EMIN(14,15) = 0.5
		EMIN(14,16) = 0.9
		EMIN(14,17) = 1.4
		EMIN(14,18) = 1.8
		EMIN(14,19) = 2.3
		EMIN(14,20) = 2.9
		EMIN(14,21) = 3.4
		EMIN(14,22) = 4.0 
		EMIN(14,23) = 4.7
		EMIN(14,24) = 5.4
		EMIN(14,25) = 6.1

		EMIN(15,16) = 0.5
		EMIN(15,17) = 0.9
		EMIN(15,18) = 1.4
		EMIN(15,19) = 1.9
		EMIN(15,20) = 2.4
		EMIN(15,21) = 3.0
		EMIN(15,22) = 3.6 
		EMIN(15,23) = 4.3
		EMIN(15,24) = 5.0
		EMIN(15,25) = 5.7

		EMIN(16,17) = 0.5
		EMIN(16,18) = 1.0
		EMIN(16,19) = 1.5
		EMIN(16,20) = 2.0
		EMIN(16,21) = 2.6
		EMIN(16,22) = 3.2
		EMIN(16,23) = 3.8
		EMIN(16,24) = 4.5
		EMIN(16,25) = 5.3

		  EMIN(17,18) = 0.5
		  EMIN(17,19) = 1.0
		  EMIN(17,20) = 1.6
		  EMIN(17,21) = 2.1
		  EMIN(17,22) = 2.7
		  EMIN(17,23) = 3.4
		  EMIN(17,24) = 4.1
		  EMIN(17,25) = 4.8


		  EMIN(18,19) = 0.6
		  EMIN(18,20) = 1.1
		  EMIN(18,21) = 1.6
		  EMIN(18,22) = 2.3
		  EMIN(18,23) = 2.9
		  EMIN(18,24) = 3.6
		  EMIN(18,25) = 4.3

		  EMIN(19,20) = 0.6
		  EMIN(19,21) = 1.1
		  EMIN(19,22) = 1.8 
		  EMIN(19,23) = 2.4
		  EMIN(19,24) = 3.1
		  EMIN(19,25) = 3.8

		  EMIN(20,21) = 0.7
		  EMIN(20,22) = 1.2 
		  EMIN(20,23) = 1.7
		  EMIN(20,24) = 2.6
		  EMIN(20,25) = 3.3

		  EMIN(21,22) = 0.7
		  EMIN(21,23) = 1.3
		  EMIN(21,24) = 2.0
		  EMIN(21,25) = 2.7

		  EMIN(22,23) = 0.7
		  EMIN(22,24) = 1.4
		  EMIN(22,25) = 2.1

		  EMIN(23,24) = 0.7
		  EMIN(23,25) = 1.5

		  EMIN(24,25) = 0.8
C	goto 7899



c	shot 11897@2140 H-mode, deltathetx = 0.15
	 
	
9876  continue   
   
	  EMIN(12,13) = 1.1
	  EMIN(12,14) = 1.7 
		EMIN(12,15) = 2.2
		EMIN(12,16) = 2.2
		EMIN(12,17) = 2.2
		EMIN(12,18) = 2.2
		EMIN(12,19) = 2.2
		EMIN(12,20) = 1.e3
		EMIN(12,21) = 1.e3
		EMIN(12,22) = 1.e3 
		EMIN(12,23) = 1.E3
		EMIN(12,24) = 1.E3
		EMIN(12,25) = 1.E3
		
		EMIN(13,14) = 0.7 
		EMIN(13,15) = 1.2
		EMIN(13,16) = 1.2
		EMIN(13,17) = 1.2
		EMIN(13,18) = 1.2
		EMIN(13,19) = 1.e3
		EMIN(13,20) = 1.e3
		EMIN(13,21) = 1.e3
		EMIN(13,22) = 1.e3 
		EMIN(13,23) = 1.E3
		EMIN(13,24) = 1.E3
		EMIN(13,25) = 1.E3

		EMIN(14,15) = 0.5
		EMIN(14,16) = 0.6
		EMIN(14,17) = 0.6
		EMIN(14,18) = 1.e3
		EMIN(14,19) = 1.e3
		EMIN(14,20) = 1.e3
		EMIN(14,21) = 1.e3
		EMIN(14,22) = 1.e3 
		EMIN(14,23) = 1.E3
		EMIN(14,24) = 1.E3
		EMIN(14,25) = 1.E3

		EMIN(15,16) = 0.5
		EMIN(15,17) = 1.e3
		EMIN(15,18) = 1.e3
		EMIN(15,19) = 1.e3
		EMIN(15,20) = 1.e3
		EMIN(15,21) = 1.e3
		EMIN(15,22) = 1.e3 
		EMIN(15,23) = 1.e3
		EMIN(15,24) = 1.e3
		EMIN(15,25) = 1.e3

		EMIN(16,17) = 0.12
		EMIN(16,18) = 0.9
		EMIN(16,19) = 1.7
		EMIN(16,20) = 2.8
		EMIN(16,21) = 4.1
		EMIN(16,22) = 5.6
		EMIN(16,23) = 7.2
		EMIN(16,24) = 8.9
		EMIN(16,25) = 10.5

		  EMIN(17,18) = 0.6
		  EMIN(17,19) = 1.5
		  EMIN(17,20) = 2.6
		  EMIN(17,21) = 3.9
		  EMIN(17,22) = 5.4
		  EMIN(17,23) = 7.0
		  EMIN(17,24) = 8.6
		  EMIN(17,25) = 10.2


		  EMIN(18,19) = 0.9
		  EMIN(18,20) = 2.0
		  EMIN(18,21) = 3.3
		  EMIN(18,22) = 4.8
		  EMIN(18,23) = 6.4
		  EMIN(18,24) = 8.1
		  EMIN(18,25) = 9.6

		  EMIN(19,20) = 1.2
		  EMIN(19,21) = 2.5
		  EMIN(19,22) = 4.0 
		  EMIN(19,23) = 5.6
		  EMIN(19,24) = 7.2
		  EMIN(19,25) = 8.8

		  EMIN(20,21) = 1.4
		  EMIN(20,22) = 2.9 
		  EMIN(20,23) = 4.5
		  EMIN(20,24) = 6.1
		  EMIN(20,25) = 7.7

		  EMIN(21,22) = 1.5
		  EMIN(21,23) = 3.1
		  EMIN(21,24) = 4.8
		  EMIN(21,25) = 6.4

		  EMIN(22,23) = 1.6
		  EMIN(22,24) = 3.3
		  EMIN(22,25) = 4.9

		  EMIN(23,24) = 1.7
		  EMIN(23,25) = 3.3

		  EMIN(24,25) = 1.6
C	goto 7899

c	*c	shot 123301, deltathetx = 0.15  
	EMIN(12,13) = 1.3
	EMIN(12,14) = 2.3 
	EMIN(12,15) = 2.8
	EMIN(12,16) = 2.9
	EMIN(12,17) = 2.9
	EMIN(12,18) = 2.9
	EMIN(12,19) = 1.E3
	EMIN(12,20) = 1.E3
	EMIN(12,21) = 1.E3
	EMIN(12,22) = 1.E3 
	EMIN(12,23) = 1.E3
	EMIN(12,24) = 1.E3
	EMIN(12,25) = 1.E3
	
	EMIN(13,14) = 1.0 
	EMIN(13,15) = 1.6
	EMIN(13,16) = 1.6
	EMIN(13,17) = 1.6
	EMIN(13,18) = 1.E3
	EMIN(13,19) = 1.E3
	EMIN(13,20) = 1.E3
	EMIN(13,21) = 1.E3
	EMIN(13,22) = 1.E3 
	EMIN(13,23) = 1.E3
	EMIN(13,24) = 1.E3
	EMIN(13,25) = 1.E3

	EMIN(14,15) = 0.6
	EMIN(14,16) = 0.7
	EMIN(14,17) = 16.0
	EMIN(14,18) = 1.E3
	EMIN(14,19) = 1.E3
	EMIN(14,20) = 1.E3
	EMIN(14,21) = 1.E3
	EMIN(14,22) = 1.E3 
	EMIN(14,23) = 1.E3
	EMIN(14,24) = 1.E3
	EMIN(14,25) = 1.E3

	EMIN(15,16) = 0.1
	EMIN(15,17) = 18.2
	EMIN(15,18) = 18.2
	EMIN(15,19) = 18.2
	EMIN(15,20) = 18.2
	EMIN(15,21) = 18.2
	EMIN(15,22) = 18.2 
	EMIN(15,23) = 18.2
	EMIN(15,24) = 18.2
	EMIN(15,25) = 18.2

	EMIN(16,17) = 0.7
	EMIN(16,18) = 2.3
	EMIN(16,19) = 4.8
	EMIN(16,20) = 7.4
	EMIN(16,21) = 9.0
	EMIN(16,22) = 9.6 
	EMIN(16,23) = 9.6
	EMIN(16,24) = 9.6
	EMIN(16,25) = 9.6

	  EMIN(17,18) = 1.7
	  EMIN(17,19) = 4.1
	  EMIN(17,20) = 6.7
	  EMIN(17,21) = 9.0
	  EMIN(17,22) = 9.0
	  EMIN(17,23) = 9.0
	  EMIN(17,24) = 9.0
	  EMIN(17,25) = 9.0


	  EMIN(18,19) = 2.5
	  EMIN(18,20) = 5.1
	  EMIN(18,21) = 6.7
	  EMIN(18,22) = 7.3
	  EMIN(18,23) = 7.3
	  EMIN(18,24) = 7.3
	  EMIN(18,25) = 7.3

	  EMIN(19,20) = 2.6
	  EMIN(19,21) = 4.3
	  EMIN(19,22) = 4.9 
	  EMIN(19,23) = 4.9
	  EMIN(19,24) = 4.9
	  EMIN(19,25) = 1.E3

	  EMIN(20,21) = 1.7
	  EMIN(20,22) = 1.7 
	  EMIN(20,23) = 2.3
	  EMIN(20,24) = 1.E3
	  EMIN(20,25) = 1.E3

	  EMIN(21,22) = 0.6 
	  EMIN(21,23) = 1.E3
	  EMIN(21,24) = 1.E3
	  EMIN(21,25) = 1.E3

	  EMIN(22,23) = 0.9
	  EMIN(22,24) = 2.8
	  EMIN(22,25) = 5.2

	  EMIN(23,24) = 2.0
	  EMIN(23,25) = 4.4

	  EMIN(24,25) = 2.5

C	goto 7899
c	shot 123302, deltathetx = 0.15
	 
	
23302 continue     
      EMIN(12,13) = 0.5
	EMIN(12,14) = 1.4 
	EMIN(12,15) = 2.5
	EMIN(12,16) = 4.0
	EMIN(12,17) = 5.8
	EMIN(12,18) = 8.0
	EMIN(12,19) = 10.4
	EMIN(12,20) = 13.2
	EMIN(12,21) = 16.2
	EMIN(12,22) = 19.4 
	EMIN(12,23) = 1.E3
	EMIN(12,24) = 1.E3
	EMIN(12,25) = 1.E3
	
	EMIN(13,14) = 0.9 
	EMIN(13,15) = 2.0
	EMIN(13,16) = 3.5
	EMIN(13,17) = 5.3
	EMIN(13,18) = 7.5
	EMIN(13,19) = 9.9
	EMIN(13,20) = 12.7
	EMIN(13,21) = 15.8
	EMIN(13,22) = 18.9 
	EMIN(13,23) = 1.E3
	EMIN(13,24) = 1.E3
	EMIN(13,25) = 1.E3

	EMIN(14,15) = 1.1
	EMIN(14,16) = 2.7
	EMIN(14,17) = 4.5
	EMIN(14,18) = 6.6
	EMIN(14,19) = 9.1
	EMIN(14,20) = 11.9
	EMIN(14,21) = 14.9
	EMIN(14,22) = 18.0 
	EMIN(14,23) = 1.E3
	EMIN(14,24) = 1.E3
	EMIN(14,25) = 1.E3

	EMIN(15,16) = 1.6
	EMIN(15,17) = 3.4
	EMIN(15,18) = 5.5
	EMIN(15,19) = 8.0
	EMIN(15,20) = 10.8
	EMIN(15,21) = 13.8
	EMIN(15,22) = 16.9 
	EMIN(15,23) = 1.e3
	EMIN(15,24) = 1.e3
	EMIN(15,25) = 1.e3

	EMIN(16,17) = 1.9
	EMIN(16,18) = 4.0
	EMIN(16,19) = 6.5
	EMIN(16,20) = 9.2
	EMIN(16,21) = 12.3
	EMIN(16,22) = 15.4
	EMIN(16,23) = 18.6
	EMIN(16,24) = 1.e3
	EMIN(16,25) = 1.e3

	  EMIN(17,18) = 2.2
	  EMIN(17,19) = 4.7
	  EMIN(17,20) = 7.4
	  EMIN(17,21) = 10.5
	  EMIN(17,22) = 13.6
	  EMIN(17,23) = 16.8
	  EMIN(17,24) = 19.7
	  EMIN(17,25) = 1.e3


	  EMIN(18,19) = 2.5
	  EMIN(18,20) = 5.3
	  EMIN(18,21) = 8.3
	  EMIN(18,22) = 11.5
	  EMIN(18,23) = 14.6
	  EMIN(18,24) = 17.5
	  EMIN(18,25) = 20.0

	  EMIN(19,20) = 2.8
	  EMIN(19,21) = 5.9
	  EMIN(19,22) = 9.0 
	  EMIN(19,23) = 12.2
	  EMIN(19,24) = 15.0
	  EMIN(19,25) = 17.6

	  EMIN(20,21) = 3.1
	  EMIN(20,22) = 6.2 
	  EMIN(20,23) = 9.4
	  EMIN(20,24) = 12.2
	  EMIN(20,25) = 14.8

	  EMIN(21,22) = 3.2
	  EMIN(21,23) = 6.3
	  EMIN(21,24) = 9.2
	  EMIN(21,25) = 11.7

	  EMIN(22,23) = 3.2
	  EMIN(22,24) = 6.1
	  EMIN(22,25) = 8.6

	  EMIN(23,24) = 2.9
	  EMIN(23,25) = 5.4

	  EMIN(24,25) = 2.6
c	update 123303 deltathetx = 0.15  3/26/12
	EMIN(12,13) = 0.15
	EMIN(12,14) = 0.20 
	EMIN(12,15) = 1.0
	EMIN(12,16) = 2.3
	EMIN(12,17) = 4.0
	EMIN(12,18) = 6.4
	EMIN(12,19) = 9.2
	EMIN(12,20) = 12.5
	EMIN(12,21) = 15.9
	EMIN(12,22) = 19.2 
	EMIN(12,23) = 1.E3
	EMIN(12,24) = 1.E3
	EMIN(12,25) = 1.E3
	
	EMIN(13,14) = 0.3 
	EMIN(13,15) = 1.2
	EMIN(13,16) = 2.4
	EMIN(13,17) = 4.2
	EMIN(13,18) = 6.5
	EMIN(13,19) = 9.4
	EMIN(13,20) = 12.6
	EMIN(13,21) = 16.0
	EMIN(13,22) = 19.3 
	EMIN(13,23) = 1.E3
	EMIN(13,24) = 1.E3
	EMIN(13,25) = 1.E3

	EMIN(14,15) = 0.8
	EMIN(14,16) = 2.1
	EMIN(14,17) = 3.9
	EMIN(14,18) = 6.2
	EMIN(14,19) = 9.1
	EMIN(14,20) = 12.3
	EMIN(14,21) = 15.7
	EMIN(14,22) = 19.0 
	EMIN(14,23) = 1.E3
	EMIN(14,24) = 1.E3
	EMIN(14,25) = 1.E3

	EMIN(15,16) = 1.3
	EMIN(15,17) = 3.1
	EMIN(15,18) = 5.4
	EMIN(15,19) = 8.3
	EMIN(15,20) = 11.5
	EMIN(15,21) = 14.4
	EMIN(15,22) = 18.2 
	EMIN(15,23) = 1.e3
	EMIN(15,24) = 1.e3
	EMIN(15,25) = 1.e3

	EMIN(16,17) = 1.8
	EMIN(16,18) = 4.2
	EMIN(16,19) = 7.0
	EMIN(16,20) = 10.2
	EMIN(16,21) = 13.7
	EMIN(16,22) = 16.9
	EMIN(16,23) = 19.4
	EMIN(16,24) = 1.e3
	EMIN(16,25) = 1.e3

	  EMIN(17,18) = 2.4
	  EMIN(17,19) = 5.2
	  EMIN(17,20) = 8.5
	  EMIN(17,21) = 11.9
	  EMIN(17,22) = 15.2
	  EMIN(17,23) = 17.7
	  EMIN(17,24) = 19.2
	  EMIN(17,25) = 19.7


	  EMIN(18,19) = 2.9
	  EMIN(18,20) = 6.2
	  EMIN(18,21) = 9.6
	  EMIN(18,22) = 12.8
	  EMIN(18,23) = 15.3
	  EMIN(18,24) = 16.8
	  EMIN(18,25) = 17.4

	  EMIN(19,20) = 3.3
	  EMIN(19,21) = 6.7
	  EMIN(19,22) = 10.0 
	  EMIN(19,23) = 12.5
	  EMIN(19,24) = 14.0
	  EMIN(19,25) = 14.6

	  EMIN(20,21) = 3.5
	  EMIN(20,22) = 6.7 
	  EMIN(20,23) = 9.2
	  EMIN(20,24) = 10.7
	  EMIN(20,25) = 11.3

	  EMIN(21,22) = 3.3
	  EMIN(21,23) = 5.8
	  EMIN(21,24) = 7.3
	  EMIN(21,25) = 7.9

	  EMIN(22,23) = 2.5
	  EMIN(22,24) = 4.0
	  EMIN(22,25) = 4.6

	  EMIN(23,24) = 1.5
	  EMIN(23,25) = 2.2

	  EMIN(24,25) = 0.6

  	goto 7899

	EMIN(12,13) = 0.5
	EMIN(12,14) = 1.4 
	EMIN(12,15) = 2.5
	EMIN(12,16) = 4.0
	EMIN(12,17) = 5.8
	EMIN(12,18) = 8.0
	EMIN(12,19) = 10.4
	EMIN(12,20) = 13.2
	EMIN(12,21) = 16.2
	EMIN(12,22) = 19.4 
	EMIN(12,23) = 1.E3
	EMIN(12,24) = 1.E3
	EMIN(12,25) = 1.E3
	
	EMIN(13,14) = 0.9 
	EMIN(13,15) = 2.0
	EMIN(13,16) = 3.5
	EMIN(13,17) = 5.3
	EMIN(13,18) = 7.5
	EMIN(13,19) = 9.9
	EMIN(13,20) = 12.7
	EMIN(13,21) = 15.8
	EMIN(13,22) = 18.9 
	EMIN(13,23) = 1.E3
	EMIN(13,24) = 1.E3
	EMIN(13,25) = 1.E3

	EMIN(14,15) = 1.1
	EMIN(14,16) = 2.7
	EMIN(14,17) = 4.5
	EMIN(14,18) = 6.6
	EMIN(14,19) = 9.1
	EMIN(14,20) = 11.9
	EMIN(14,21) = 14.9
	EMIN(14,22) = 18.0 
	EMIN(14,23) = 1.E3
	EMIN(14,24) = 1.E3
	EMIN(14,25) = 1.E3

	EMIN(15,16) = 1.6
	EMIN(15,17) = 3.4
	EMIN(15,18) = 5.5
	EMIN(15,19) = 8.0
	EMIN(15,20) = 10.8
	EMIN(15,21) = 13.8
	EMIN(15,22) = 16.9 
	EMIN(15,23) = 1.e3
	EMIN(15,24) = 1.e3
	EMIN(15,25) = 1.e3

	EMIN(16,17) = 1.9
	EMIN(16,18) = 4.0
	EMIN(16,19) = 6.5
	EMIN(16,20) = 9.2
	EMIN(16,21) = 12.3
	EMIN(16,22) = 15.4
	EMIN(16,23) = 18.6
	EMIN(16,24) = 1.e3
	EMIN(16,25) = 1.e3

	  EMIN(17,18) = 2.2
	  EMIN(17,19) = 4.7
	  EMIN(17,20) = 7.4
	  EMIN(17,21) = 10.5
	  EMIN(17,22) = 13.6
	  EMIN(17,23) = 16.8
	  EMIN(17,24) = 19.7
	  EMIN(17,25) = 1.e3


	  EMIN(18,19) = 2.5
	  EMIN(18,20) = 5.3
	  EMIN(18,21) = 8.3
	  EMIN(18,22) = 11.5
	  EMIN(18,23) = 14.6
	  EMIN(18,24) = 17.5
	  EMIN(18,25) = 20.0

	  EMIN(19,20) = 2.8
	  EMIN(19,21) = 5.9
	  EMIN(19,22) = 9.0 
	  EMIN(19,23) = 12.2
	  EMIN(19,24) = 15.0
	  EMIN(19,25) = 17.6

	  EMIN(20,21) = 3.1
	  EMIN(20,22) = 6.2 
	  EMIN(20,23) = 9.4
	  EMIN(20,24) = 12.2
	  EMIN(20,25) = 14.8

	  EMIN(21,22) = 3.2
	  EMIN(21,23) = 6.3
	  EMIN(21,24) = 9.2
	  EMIN(21,25) = 11.7

	  EMIN(22,23) = 3.2
	  EMIN(22,24) = 6.1
	  EMIN(22,25) = 8.6

	  EMIN(23,24) = 2.9
	  EMIN(23,25) = 5.4

	  EMIN(24,25) = 2.6

	goto 7899
c	shot 123302, deltathetx = 0.20
	EMIN(8,9)   = 0.8
	EMIN(8,10)  = 1.3
	EMIN(8,11)  = 1.4
	EMIN(8,12)  = 1.4
	EMIN(8,13)  = 1.4
	EMIN(8,14)  = 1.4
	EMIN(8,15)  = 1.E3
	EMIN(8,16)  = 1.E3
	EMIN(8,17)  = 1.E3
	EMIN(8,18)  = 1.E3
	EMIN(8,19)  = 1.E3
	EMIN(8,20)  = 1.E3
	EMIN(8,21)  = 1.E3
	EMIN(8,22)  = 1.E3
	EMIN(8,23)  = 1.E3
	EMIN(8,24)  = 1.E3
	EMIN(8,25)  = 1.E3

	EMIN(9,10)  = 0.5
	EMIN(9,11)  = 0.7
	EMIN(9,12)  = 0.7
	EMIN(9,13)  = 0.7
	EMIN(9,14)  = 1.E3
	EMIN(9,15)  = 1.E3
	EMIN(9,16)  = 1.E3
	EMIN(9,17)  = 1.E3
	EMIN(9,18)  = 1.E3
	EMIN(9,19)  = 1.E3
	EMIN(9,20)  = 1.E3
	EMIN(9,21)  = 1.E3
	EMIN(9,22)  = 1.E3
	EMIN(9,23)  = 1.E3
	EMIN(9,24)  = 1.E3
	EMIN(9,25)  = 1.E3

	EMIN(10,11)  = 0.2
	EMIN(10,12)  = 0.2
	EMIN(10,13)  = 1.E3
	EMIN(10,14)  = 1.E3
	EMIN(10,15)  = 1.E3
	EMIN(10,16)  = 1.E3
	EMIN(10,17)  = 1.E3
	EMIN(10,18)  = 1.E3
	EMIN(10,19)  = 1.E3
	EMIN(10,20)  = 1.E3
	EMIN(10,21)  = 1.E3
	EMIN(10,22)  = 1.E3
	EMIN(10,23)  = 1.E3
	EMIN(10,24)  = 1.E3
	EMIN(10,25)  = 1.E3

	EMIN(10,12)  = 0.2
 	EMIN(11,13)  = 0.5
	EMIN(11,14)  = 1.2
	EMIN(11,15)  = 2.0
	EMIN(11,16)  = 3.1
	EMIN(11,17)  = 4.5
	EMIN(11,18)  = 6.1
	EMIN(11,19)  = 7.9
	EMIN(11,20)  = 10.0
	EMIN(11,21)  = 12.3
	EMIN(11,22)  = 14.7
	EMIN(11,23)  = 17.0
	EMIN(11,24)  = 19.2
	EMIN(11,25)  = 1.E3
	  
	EMIN(12,13) = 0.4
	EMIN(12,14) = 1.1 
	EMIN(12,15) = 1.9
	EMIN(12,16) = 3.0
	EMIN(12,17) = 4.4
	EMIN(12,18) = 6.0
	EMIN(12,19) = 7.8
	EMIN(12,20) = 9.9
	EMIN(12,21) = 12.2
	EMIN(12,22) = 14.6 
	EMIN(12,23) = 16.9
	EMIN(12,24) = 19.1
	EMIN(12,25) = 1.E3
	
	EMIN(13,14) = 0.7 
	EMIN(13,15) = 1.5
	EMIN(13,16) = 2.7
	EMIN(13,17) = 4.0
	EMIN(13,18) = 5.6
	EMIN(13,19) = 7.5
	EMIN(13,20) = 9.6
	EMIN(13,21) = 11.8
	EMIN(13,22) = 14.2 
	EMIN(13,23) = 16.7
	EMIN(13,24) = 18.7
	EMIN(13,25) = 1.E3

	EMIN(14,15) = 0.9
	EMIN(14,16) = 2.0
	EMIN(14,17) = 3.4
	EMIN(14,18) = 5.0
	EMIN(14,19) = 6.8
	EMIN(14,20) = 8.9
	EMIN(14,21) = 11.2
	EMIN(14,22) = 13.5 
	EMIN(14,23) = 15.9
	EMIN(14,24) = 18.1
	EMIN(14,25) = 20.0

	EMIN(15,16) = 1.2
	EMIN(15,17) = 2.5
	EMIN(15,18) = 4.2
	EMIN(15,19) = 6.0
	EMIN(15,20) = 8.1
	EMIN(15,21) = 10.4
	EMIN(15,22) = 12.7 
	EMIN(15,23) = 15.1
	EMIN(15,24) = 17.2
	EMIN(15,25) = 19.1

	EMIN(16,17) = 1.4
	EMIN(16,18) = 3.0
	EMIN(16,19) = 4.9
	EMIN(16,20) = 7.0
	EMIN(16,21) = 9.2
	EMIN(16,22) = 11.6
	EMIN(16,23) = 14.0
	EMIN(16,24) = 16.1
	EMIN(16,25) = 18.0

	  EMIN(17,18) = 1.7
	  EMIN(17,19) = 3.5
	  EMIN(17,20) = 5.6
	  EMIN(17,21) = 7.9
	  EMIN(17,22) = 10.2
	  EMIN(17,23) = 12.6
	  EMIN(17,24) = 14.7
	  EMIN(17,25) = 16.6


	  EMIN(18,19) = 1.8
	  EMIN(18,20) = 4.0
	  EMIN(18,21) = 6.3
	  EMIN(18,22) = 8.6
	  EMIN(18,23) = 11.0
	  EMIN(18,24) = 13.1
	  EMIN(18,25) = 15.0

	  EMIN(19,20) = 1.9
	  EMIN(19,21) = 4.4
	  EMIN(19,22) = 7.0 
	  EMIN(19,23) = 9.1
	  EMIN(19,24) = 11.3
	  EMIN(19,25) = 13.2

	  EMIN(20,21) = 2.1
	  EMIN(20,22) = 4.7 
	  EMIN(20,23) = 7.1
	  EMIN(20,24) = 9.2
	  EMIN(20,25) = 11.1

	  EMIN(21,22) = 2.3
	  EMIN(21,23) = 4.8
	  EMIN(21,24) = 6.9
	  EMIN(21,25) = 8.8

	  EMIN(22,23) = 2.4
	  EMIN(22,24) = 4.6
	  EMIN(22,25) = 6.5

	  EMIN(23,24) = 2.2
	  EMIN(23,25) = 4.1

	  EMIN(24,25) = 2.0


c	123301 rmp

c	deltathetx = 0.1	 
	E1(17,18) = 2.5e3
	E1(17,19) = 6.2e3
	E1(17,20) = 7.7e3
 	E1(17,21) = 10.0e3
	E1(17,22) = 10.0e3
	E1(17,23) = 10.0e3
	E1(17,24) = 10.9e3
	E1(17,25) = 13.4e3
	 

	E1(18,19) = 3.7e3
	E1(18,20) = 7.6e3
	E1(18,21) = 10.0e3
	E1(18,22) = 10.0e3
	E1(18,23) = 10.0e3
	E1(18,24) = 10.0e3
	E1(18,25) = 10.9e3
	
	
      E1(19,20) = 3.7e3
	E1(19,21) = 6.4e3
	E1(19,22) = 7.3e3
	E1(19,23) = 7.3e3
	E1(19,24) = 7.3e3
	E1(19,25) = 10.0e5
	

	E1(20,21) = 2.5e3
	E1(20,22) = 2.5e3
	E1(20,23) = 3.4e3
	E1(20,24) = 10.0e5
	E1(20,25) = 10.0e5


	E1(21,22) = 1.3e3
	E1(21,23) = 10.0e5
	E1(21,24) = 10.0e5
	E1(21,25) = 10.0e5


	E1(22,23) = 1.3e3
	E1(22,24) = 4.1e3
	E1(22,25) = 7.7e3
	

	E1(23,24) = 2.9e3
	E1(23,25) = 6.5e3
	
	E1(24,25) = 3.7e3
c*****************this is where the Emin (m,mp) go ****************************

C	X-LOSS FRACTION, CUMULATIVE WITH RADIUS******old way**********************
	do 5310 n = 1, 25
	Exloss(n) = 100.0e3
5310	continue
c	delthetx = 0.1 
	Exloss(17) = 13.4e3
	Exloss(18) = 10.9e3
	Exloss(19) = 10.9e3
	Exloss(20) = 10.9e3
	Exloss(21) = 10.9e3
	Exloss(22) = 7.7e3
	Exloss(23) = 6.5e3
	Exloss(24) = 3.7e3

c	delthetx = 0.15

c	exloss(16) = 9.6e3
c	exloss(17) = 9.0e3
c	exloss(18) = 7.3e3
c	exloss(19) = 7.3e3
c	exloss(20) = 7.3e3
c	exloss(21) = 7.3e3
c	exloss(22) = 5.2e3
c	exloss(23) = 4.4e3
c	exloss(24) = 2.5e3

c	delthetx = 0.2
c	Exloss(8)  = 7.8e3
c	Exloss(9)  = 7.8e3
c	Exloss(10) = 7.8e3
c	Exloss(11) = 7.8e3
c	Exloss(12) = 7.8e3
c	Exloss(13) = 7.8e3
c	Exloss(14) = 7.8e3
c	Exloss(15) = 7.8e3
c	Exloss(16) = 7.2e3
c	Exloss(17) = 6.7e3
c	Exloss(18) = 5.5e3
c	exloss(19) = 5.5e3
c	exloss(20) = 5.5e3
c	exloss(21) = 5.5e3
c	Exloss(22) = 3.7e3
c	Exloss(23) = 3.3e3
c	Exloss(24) = 1.9e3
c	calculate cumulative xloss fractions as a function of radius.
c	assumption: all particles at a given radius are swept through xloss
c	region, so that all particles with energy > Exloss and lost.

	do 5315 n = 1,25
	if(n.lt.nx) then
	fx(n) = 0.0
	ex(n) = 0.0 
	go to 5315 
	endif
	 
      exmin = Exloss(n-1)/xti(n)
      A = 1.5
	X = ExMIN
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	FX(n) =  (value1 - value)/value1
	
	A = 2.5
	EVALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	EX(n) = (evalue1 - evalue)/evalue1

	exmin = Exloss(n)/xti(n)
	A = 1.5
	X = ExMIN
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	FX(n) =  (value1 - value)/value1  - fx(n) 
	fx(n) = fx(n-1)+fx(n)
	A = 2.5
	EVALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	EX(n) = (evalue1 - evalue)/evalue1 - ex(n)
	ex(n) = ex(n-1) + ex(n)
5315  continue	 
5301	write (129, '(1x,35A)') '  rho     EXLOSS      IONs      ENERGY  '    
5302	format(f6.3,e10.3,2f10.3) 
	do 5320 n = 1, 25
	write(129,5302) rhor(n),Exloss(n), fx(n),ex(n)
5320  continue
	 
c	min W to gradB drift beyond 2nd mesh interval
	

	do 6900 n = 1,25
	sinknx(n) = 0.
	sinkex(n) = 0.
	xesource(n) = 0.
	xnsource(n) = 0.
6900  continue
	sumsinkn = 0.0
	sumsinke = 0.0
	sumsepn = 0.0
	sumsepe = 0.0
	fraclosseold = 0.0
	fraclossnold = 0.0
c	total ion (#/s) and ion energy (Watts) input rates to X-loss region
c	at each radial mesh point 
	do 6910 n = nx,24
	xniexp = exne(n)/(atnum(1)+fracz*zbar(n))
	rr = rhor(n)*aminor*SQRT(0.5*(1.+ELONG**2))
	sinknx(n) =	(xniexp*abs(erex(n)/bphi))*
	2			(6.28*rmajor*delrho)
	sinkex(n) =	xk*xti(m)*(xniexp*abs(erex(n)/bphi))*
	2			(6.28*rmajor*delrho)
	

c	X-transport fractions
	do 6905 m = n+1,25  
      exmin = E1(n,m)/xti(n)
	A = 1.5
	X = ExMIN
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	fraclossn(n,m) =  (value1 - value)/value1
	A = 2.5
	EVALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	fraclosse(n,m) = (evalue1 - evalue)/evalue1	
6905  continue	 
c	fraction gradB drifting across separatrix
	exmin = E1(n,25)/xti(n)	
	A = 1.5
	X = ExMIN
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	fracseplossn =  (value1 - value)/value1
	A = 2.5
	EVALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	fracseplosse = (evalue1 - evalue)/evalue1
	  
c	# part & amount energy grad B loss across separatrix
     	seplossn(n) = sinknx(n)*fraclossn(n,25)
	seplosse(n) = sinkex(n)*fraclosse(n,25)
	sumsepn = sumsepn + sepnx(n)
	sumsepe = sumsepe + sepex(n)

c	# part & amount of energy gradB drift down into next region      
	sinkex(n) = sinkex(n)*fraclosse(n,n+1)
	sinknx(n) = sinknx(n)*fraclossn(n,n+1)
c	fraclosseold = fraclosse
c	fraclossnold = fraclossn
c	# part participating in X-transport or X-loss
	sumsinkn = sumsinkn + sinknx(n)
	sumsinke = sumsinke + sinkex(n)

c	# part & amount energy X-loss
c     	sepnx(n) = sinknx(n)*fraclossn
c    	sepex(n) = sinkex(n)*fraclosse
	
	sumxtransn = sumsinkn - sumsepn
	sumxtranse = sumsinke - sumsepe
6910	continue
c	sources into other mesh intervals & X-loss
	sumsourcen = 0.0
	sumsourcee = 0.0
	do 6920 m = nx + 1,25
	do 6915 n = nx, m-1 
	xniexp = exne(n)/(atnum(1)+fracz*zbar(n))
 	rr = rhor(n)*aminor*SQRT(0.5*(1.+ELONG**2))
	xnn =	(xniexp*abs(erex(n)/bphi))/(6.28*rr)
	xee =	xk*xti(n)*(xniexp*abs(erex(n)/bphi))/(6.28*rr)
	A = 1.5
	X = E1(n,m)/xti(n)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	flossn1 =  (value1 - value)/value1
	if(m.eq.25) then
		flossn2 = 0.0
		goto 6911
	endif
	X = E1(n,m+1)/xti(n)
	VALUE = GAMI(A,X) 
	flossn2 =  (value1 - value)/value1
6911	xnsource(m) = xnsource(n) + xnn*(flossn1 - flossn2) 

	A = 2.5
	X = E1(n,m)/xti(n)
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	flosse1 =  (value1 - value)/value1
	if(m.eq.25) then
		flosse2 = 0.0
		goto 6912
	endif
	X = E1(n,m+1)/xti(n)
	VALUE = GAMI(A,X) 
	flosse2 =  (value1 - value)/value1
6912	xesource(m) = xesource(n) + xee*(flosse1 - flosse2) 
	sumsourcen = sumsourcen + xnsource(m)
	sumsourcee = sumsourcee + xesource(m) 
6915	continue
6920  continue
	rr = rhor(25)*aminor*SQRT(0.5*(1.+ELONG**2))
 
	xlossnsum = xnsource(25)*6.28*rr
	xlossesum = xesource(25)*6.28*rr

C	*****************X-TRANSPORT***7/27/2011*************************
c	minimum energy EMIN(m,m') for transport from mesh m to mesh m'>m (keV)
c	must be created by hand from previous XLOSS calculation 7/27/2011
C	1.E3 INDICATES THAT EMIN > 20 KEV, EMIN THE SAME FOR SEVERAL MESH INDICATES
C	THAT IONS PASS THROUGH THESE MESH 
7899	 continue
c    IOPXTRAN = 1 to include X-transport
	ioptxtran = 1 
C	**********************************



c	************calculate transport sinks and sources**********************
c***************new way ***************************************************
c	epmin in keV, xti in eV
	do 3510 m = nx,24
	do 3505 mp = m+1,25
	epmin(m,mp) = emin(m,mp)/(1.e-3*xti(m))
3505	continue
3510	continue
c	*********calculate particle & energy transport from mesh m to mesh mp********
c*****************old way**********************************************************
	do 3520 m = nx,24
	do 3515 mp = m+1,25	
C	PARTICLES		
	A = 1.5
	X =  epmin(m,mp)
	VALUE = GAMIC(A,X)
	VALUE1= GAMMA(A)
	value2 = gami(A,X)
	flossn1 =  (value1-value2)/value1



	Y =  epmin(m,mp+1)
	VALUE = GAMIC(A,Y)
	value2 = gami(A,Y)
	flossn2 =  (value1-value2)/value1
	XXN(m,mp) = flossn1 - flossn2
		if(mp.eq.25) XXN(m,mp) = flossn1
 
c	if(abs(xxn(m,mp)).lt.1.e-2) xxn(m,mp) = 0.0

C	ENERGY
	A = 2.5
	VALUE = GAMIC(A,X)
	VALUE1= GAMMA(A)
	value2 = gami(A,X)
	flosse1 =  (value1-value2)/value1
	value2 = gami(A,Y)
	flosse2 =  (value1-value2)/value1
	XXE(m,mp) = flossE1 - flossE2
	if(mp.eq.25) XXE(m,mp) = flosse1

c	if(abs(xxe(m,mp)).lt.1.e-2) xxe(m,mp) = 0.0
3515	continue
3520	continue
c	****************calculate particle and energy sinks at each mesh *************
	do 3524 m = 1,25
	sinknx(m) = 0.0
	sinkex(m) = 0.0
3524 	continue
	do 3530 m = nx,24
	xniexp = exne(m)/(atnum(1)+fracz*zbar(m))
 	sinknx(m) =	(xniexp*abs(erex(m)/bphi))
	sinkex(m) =	xk*xti(m)*(xniexp*abs(erex(m)/bphi))
	sinknsum = 0.0
	sinkesum = 0.0
	do 3525 mp =m+1,25
	sinknsum = sinknsum + XXN(m,mp)
	sinkesum = sinkesum + XXE(m,mp)
3525	continue
	sinknx(m) = sinknx(m)*sinknsum 
	sinkex(m) = sinkex(m)*sinkesum
	xl = 6.28*rhor(m)*aminor*sqrt(0.5*(1.+elong**2))
	sinknx(m) = sinknx(m)/xl
	sinkex(m) = sinkex(m)/xl
3530	continue
c	**********calculate particle and energy sources at each mesh**************
	do 3534 m = 1,25
	xnsource(m) = 0.0
	xesource(m) = 0.0
3534	continue 
	do 3545 m = nx+1,25 
	do 3540 mp = nx,m-1
	xniexp = exne(mp)/(atnum(1)+fracz*zbar(mp))
	xnsource(m) = xnsource(m) +	(xniexp*abs(erex(mp)/bphi))*XXN(mp,m)
	xesource(m) = xesource(m) + xk*xti(mp)*(xniexp*abs(erex(mp)/bphi))
	2				*XXE(mp,m)
	xl = 6.28*rhor(m)*aminor*sqrt(0.5*(1.+elong**2))
 
	xnsource(m) = xnsource(m)/xl
	xesource(m) = xesource(m)/xl
3540  continue
3545	continue 

c***********************new way 11/10/11******happy birthday USMC**************
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c	ExB drift ion flux into x-region from plasma
	do 3550 n = nx,25
	xniexp = exne(n)/(atnum(1)+fracz*zbar(n)) 
	gamin(n) =  0.5*xniexp*abs(erex(n)/bphi)
c	loss averaged over flux surface
	xl = 6.28*rhor(n)*aminor*sqrt(0.5*(1.+elong**2))
 	xtransink(n) = gamin(n)/xl
3550	continue

c	ExB drift out of x-region into plasma
	do 3560 n = nx+1,25
	if(n.lt.25) then
	 delr = (rhor(n+1)-rhor(n))*aminor*sqrt(0.5*(1.+elong**2))
	 depdr =	 (epmin(nx,n+1)-epmin(nx,n))/delr
	endif
	gamout(n) = 0.5*gamin(nx)*sqrt(epmin(nx,n))*exp(-1.0*epmin(nx,n))
	1			*depdr*delr 
	do 3555 np= nx,n-1
	depdr = (epmin(np,n+1)-epmin(np,n))/delr
	gamout(n) = gamout(n) + gamin(np)*sqrt(epmin(np,n))*
     1			exp(-1.0*epmin(np,n))*depdr*delr 
3555	continue
	gamout(n) = gamout(n)/GAMMA(1.5) 
	xlg = 6.28*rhor(n)*aminor*sqrt(0.5*(1.+elong**2))
	xtransource(n) = gamout(n)/xlg 
3560  continue

C******************Xloss fractions*******************************************
	do 8315 n = 1,25
	if(n.lt.nx) then
	fx(n) = 0.0
	ex(n) = 0.0 
	go to 8315 
	endif
	 
      exmin = 1e3*Emin(n,25)/xti(n)
      A = 1.5
	X = ExMIN
	VALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	FX(n) =  (value1 - value)/value1
	
	A = 2.5
	EVALUE = GAMI(A,X)
	VALUE1= GAMMA(A)
	EX(n) = (evalue1 - evalue)/evalue1

	
8315  continue	 
8301	write (129, '(1x,35A)') '  rho   EMIN(n,25)    FX       EX  '    
8302	format(f6.3,e10.3,2f10.3) 
	do 8320 n = 1, 25
	write(129,5302) rhor(n),Emin(n,25), fx(n),ex(n)
8320	continue 

c***********************current density xtran*********************************

c	x-trans current density, assuming the ion drifts from one radius
c	to another is uncompensated by an electron flow, and spreading 
c	the ion flow over the flux surface to get an average current density
	

C	COMPENSATING INWARD CURENT DENSITY FOR XTRAN
	do 5550 mp = nx,25
	
	xniexp = exne(mp)/(atnum(1)+fracz*zbar(mp))
 	xebin(mp) =	(xniexp*abs(erex(mp)/bphi))
5550  continue
	
	do 5555 mp = nx,24
	xebout(mp) = 0.0
	do 5554 mpp= nx,mp-1
	delr = (rhor(mp+1)-rhor(mp))*aminor*sqrt(0.5*(1.+elong**2))
	X = epmin(mpp,mp)
	DXDR = (epmin(mpp,mp)-epmin(mpp,mp-1))/delr
	GX =2.*sqrt(X/3.1416)*exp(-1.*X)*DXDR 
	gxg(mpp,mp) = gx
	dxdrg(mpp,mp) = dxdr
	
5554	continue
	
5555	continue
	do 5558 mp = nx,25
	xebout(mp) = 0.0
	do 5557 mpp = nx, mp-1
	xebout(mp) = xebout(mp) + 0.5*(xebin(mpp)*gxg(mpp,mp) +
     1            xebin(mpp+1)*gxg(mpp+1,mp))*delr
5557	continue
      xnet(mp) = 0.5*(xebin(mp)+xebin(mp+1)) - xebout(mp)
5558	continue	      	 

	do 5560 m = nx+1,25
	cdxtran(m) = 0.0
	do 5559 mp = nx, m
	
	cdxtran(m) = cdxtran(m) - (eq/6.28)*xnet(mp)*delr 
5559  continue
c	reduction of the # of ions that can exb drift into x-region
	cdxtran(m) = 0.5*cdxtran(m)
5560	continue
	z =delr
	
	


	
c	radr =  rhor(m)*aminor*sqrt(0.5*(1.+elong**2))
c	A = 1.5
c	X =  epmin(m,mp)
c	VALUE = GAMIC(A,X)
c	VALUE1= GAMMA(A)
c	value2 = gami(A,X)
c	flossn1 =  (value1-value2)/value1
c	cdxtran(mp) = cdxtran(mp) + eq*exbin*delr*flossn1/(6.28*radr)
c5555	continue
c5560  continue
c	xtrans current density, assuming that electron flow compensates
c	the ions that return to the plasma in order to maintain charge 
c	neutrality and only the ions escaping across the separatrix
c	constitute an outward current 
c	cdxtran(nx-1) = 0.0
c	do 5565 m = nx,24
c	xniexp = exne(m)/(atnum(1)+fracz*zbar(m))
c	exbin =	(xniexp*abs(erex(m)/bphi))
c	delr = (rhor(m+1)-rhor(m))*aminor*sqrt(0.5*(1.+elong**2))
c	radr =  rhor(m)*aminor*sqrt(0.5*(1.+elong**2))
c	A = 1.5
c	X =  epmin(m,25)
c	VALUE = GAMIC(A,X)
c	VALUE1= GAMMA(A)
c	value2 = gami(A,X)
c	flossn1 =  (value1-value2)/value1
c 	cdxtran(m) = cdxtran(m-1) + eq*exbin*delr*flossn1/(6.28*radr)
c5565	continue
 
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



C**************************END OF X-LOSS & X-TRANSPORT CALCULATION*************
c*******will calculate below the effect on total particle and heat fluxes**********


	 
C	**************************************************************
c	calculate heat and particle flux distributions 
C	***************************************************************
	gamheat(25) = fluxheat
	gamheate(25) = fheate*gamheat(25)
	gamheati(25) =(1.-fheate)*gamheat(25)

      call edgerotran(25,10)
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
c	zbar2(25) = 6.0 
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

	ssion(25) = yni(25,1)*XNUIONI(25)*(1.+fracz*zbar2(25))
c    1			+  yni(25,1)*xnuionb(25)

	do 1575 j = 1,24
	n = 25-j
c	new ion & plasma density and ion & electron temp formulation
	dens(n) = 0.5*(yni(n+1,1)+yni(n,1) )
	tele = 0.5*(tel(n)+tel(n+1))
	tiav = 0.5*(ti(n)+ti(n+1))
	rz = (atnum(1)**2/xmas(1) + zbar2(n)**2/xmas(2))*1.67e-27  
	cequil = 6.32e-14*rz 
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)
c	eq 4.90 Fusion Plasma Physics for equilibration
	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(yne) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(tele**1.5)*X/Yz)
 
	cequil = 7.9e-42*couloge*zeff/xmas(1) 

	qie(n) = cequil*yne*(tiav-tele)/(tele**1.5)
	cxcool(n)=1.5*dens(n)*tiav*xk*xnuati(n)
c	if(n.eq.24) cxcool(n) = 1.5*dens(n)*tiav*xk*xnuati(n)
c**********************
	cmulteq = 1.0

c	electron and ion heat flux and particle flux calculations
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
c	zbar2(n) = 6.0 
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
	radmultedge = 1.0
	radcool(n) = radmultedge*radcool(n)
	sourceion(n) = yne*xnuioni(n)

	
c	particle fluxes of ions & impurities
c 	thetint = 0.5
c	DENS(N) = thetint*YNI(N+1,1) + (1.-thetint)*YNI(N,1)
	delma = delna 
	GAMION(N,1) = GAMION(N+1,1) -dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n))


	GAMION(N,2) = GAMION(N+1,2) 
     1		-	DENS(N)*XNUIONz(N)*DELMA*(1.+fracz*zbar2(n))
	ssion(n) =  yni(N,1)*XNUIONI(N)*(1.+fracz*zbar2(n))
c	1			+ yni(n,1)*xnuionb(n)	
	call edgerotran(n,10)
   	beamdot(n) =  dens(n)*xnuionb(n)
c	*************heat fluxes******************************************	

	gamheati(n) = gamheati(n+1) + delma*(cxcool(n) + cmulteq*qie(n))
 	gamheati(n) = gamheati(n) - fionb(n)*qnb(n)*delma

c	*****************************************************************
	gamheate(n)=gamheate(n+1)+delma*
     1					(coolion(n)+radcool(n)-cmulteq*qie(n))
	gamheate(n) = gamheate(n) -	(1.-fionb(n))*qnb(n)*delma
c	goto 1572



c	*****************x-transport corrections************************	
c	*****************old way****************************************
C	ion flux reduced by the part contributed by X-transport
	do 1940 n = 1,nx
	gamionX(n) = gamion(n,1)
	gamheatiX(n) = gamheati(n)
1940	continue
 
	do 1941 n = nx,24
	
	gamionX(n+1) = gamionX(n) +dens(n)*xnuionb(n)*delma	
     1			+	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n))
     2 + 0.5*(xnsource(n)+xnsource(n+1)-sinknx(n)- sinknx(n+1))*delma
		 
	gamheatiX(n+1) = gamheatiX(n) - delma*(cxcool(n) - cmulteq*qie(n))
     1   			+ fionb(n)*qnb(n)*delma
	2 + 0.5*(xesource(n)+xesource(n+1)-sinkex(n)-sinkex(n+1))*delma
1941	continue
c	*********new way 11/10/11****************************************
C	ion flux reduced by the part contributed by X-transport
	do 1950 n = 1,nx
	gamionX(n) = gamion(n,1)
	gamheatiX(n) = gamheati(n)
1950	continue
 
	do 1951 n = nx,24
	
	gamionX(n+1) = gamionX(n) +dens(n)*xnuionb(n)*delma	
     1			+	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n))
     2 + (xtransource(n)-xtransink(n))*delma
		 
	gamheatiX(n+1) = gamheatiX(n) - delma*(cxcool(n) - cmulteq*qie(n))
     1   			+ fionb(n)*qnb(n)*delma
	2 + 0.5*(xesource(n)+xesource(n+1)-sinkex(n)-sinkex(n+1))*delma
1951	continue
c***********************************************************************
c	time dependent corrections
	dnped_dt = 0.5*(dlnn_dt(n)*yni(n,1) +	dlnn_dt(n+1)*yni(n+1,1))
	une = 0.25*(yni(n,1)+yni(n+1,1))*atnum(1)*xk*(tel(n)+tel(n+1)) +
     1		0.25*(yni(n,2)+yni(n+1,2))*zbar2(n)*xk*(tel(n)+tel(n+1))
	dweped_dt= 0.5*(dlnwe_dt(n)+dlnwe_dt(n+1))*1.5*une
	uni = 0.25*(yni(n,1)+yni(n+1,1))*xk*(ti(n)+ti(n+1))
	dwiped_dt= 0.5*(dlnwi_dt(n)+dlnwi_dt(n+1))*1.5*uni  
	gamheati(n) = gamheati(n)+dwiped_dt*delma
	gamheate(n) = gamheate(n)+dweped_dt*delma
	gamion(n,1) = gamion(n,1)+dnped_dt*delma
 	gamionX(n)  = gamionX(n)+dnped_dt*delma
 	gamheat(n) = gamheati(n)+gamheate(n)
1575	continue 

c	adjustment to surface flux boundary conditions
c	to account for time-dependent n and nT in pedestal
1572  continue
c	goto 1573
	ye=0.
	yi=0.
	yn=0.
	dlnpeped_dt =dlnpped_dt
	dlnpiped_dt =dlnpped_dt
	do 1576 n = 1,24
	une = 0.25*(yni(n,1)+yni(n+1,1))*atnum(1)*(tel(n)+tel(n+1))*xk +
     1	  0.25*(yni(n,2)+yni(n+1,2))*zbar2(n)*(tel(n)+tel(n+1))*xk 
	uni = 0.25*(yni(n,1)+yni(n+1,1))*(ti(n)+ti(n+1))*xk 
c	set time derivatives to input profile values data20
	
	ye = ye + 0.5*(dlnwe_dt(n)+dlnwe_dt(n+1))*1.5*une*delma 
	yi = yi + 0.5*(dlnwi_dt(n)+dlnwi_dt(n+1))*1.5*uni*delma
	yn = yn + 0.5*(dlnn_dt(n)+dlnn_dt(n+1))*yni(n,1)*delma
1576  continue
c	adjust heat & particle flux for change in boundary conditions
c	***5/17/07***assume gam constant between ELMs ***no adjustment************
c	*******6/22/07  adjust edge fluxes************************
	do 1577 n= 1,25 
	gamheati(n) = gamheati(n) - yi
	gamheate(n) = gamheate(n) - ye
	gamion(n,1) = gamion(n,1) - yn
	gamionX(n)  = gamionX(n)  - yn
1577	continue	 






c	*****************CALCULATE XNUDRAG ***************************************
c	**************************************************************************
c	goto 1579
c**************above temporary********************************
	do 1578 n = 1,25

c	***************pressure gradient************************
c*********calculated from input gradients*******************	
      xlp = 1./xlpm(n) 

	 
	do 1908 j=1,2
	xz = atnum(1)

	if(j.eq.2) xz = zbar2(n)
	PRESS(J)=-1.*(ti(n)/(xz*BTHET))*xlpm(n)
1908	continue
	press1(n) = press(1)
	press2(n) = press(2) 
c	***********************************************
	ioptapproach = 0
c	iolopt = 1 for radial D flux reduced by ion orbit loss (including return current) 
      IOLOPT = 1
	ziol = 1.0

	if(iolopt.eq.1) ziol = 1. - 2.*forbL(n)
 	open(500,file='bugging.txt',status='old')
	vrad1(n) = gamion(n,1)*ziol/yni(n,1)
	vrad2(n) = gamion(n,2)/yni(n,2)

C	******calculate nudrag-exp & poloidal rotation & nudragyro & vphical****** 
c	ntorque=20 makes edgerotran calculate only momentum input & return w/o poloidal calc	
	call edgerotran(n,20)
c	inferred nudrag	
	xnuc(1,1) = xnuc11(n)
	xnuc(1,2) = xnuc12(n)
	xnuc(2,1) = xnuc21(n)
	xnuc(2,2) = xnuc22(n)
	xni(1) = yni(n,1)
	xni(2) = yni(n,2)
	xmtor(1) = xmomtor1(n)
c		xmv = 1.13*sqrt(2.*xk*xti(n)/xm)*xmorbl(n)
c*******includes iol intrinsic rotation in xnud inference************************
c		xmomiol(n) = yni(n,1)*xmas(1)*yy1(n) 
c		xmtor(1) = xmtor(1) + xmomiol(j)  

	xmtor(2) = xmomtor2(n) 
c	ziol = 1.0
c	if(iolopt.eq.1) ziol = 1. - 2.*forbL(n)
	ioptnu = 0 
c*******if ioptnu = 1, calculate nudrag w input iol correction*********************
c*******if ioptnu = 0, nudrag is calculated using the calculated iol correction****
	y11 = xmtor(1)+atnum(1)*eq*(xni(1)*ephia+bthet*gamion(n,1)*ziol)
 	y1 = y11/(xni(1)*xmas(1)*xnuc(1,2))
	torq1(n) = y11
c*****ion fractions for 100% loss**********need to input for every new shot********	
	if(ioptnu.eq.1) then
c********shot123302 @ 2600ms**********************************************
		fop(1) =  .020
		fop(2) =  .023
		fop(3) =  .026
		fop(4) =	 .029
		fop(5) =	 .032
		fop(6) =	 .036	
		fop(7) =	 .041
		fop(8) =	 .045
		fop(9) =	 .052
		fop(10) = .059
		fop(11) = .066
		fop(12) = .074
		fop(13) = .085
		fop(14) = .097
		fop(15) = .111
		fop(16) = .134
		fop(17) = .149
		fop(18) = .174
		fop(19) = .198
		fop(20) = .237
		fop(21) = .289
		fop(22) = .373
		fop(23) = .525
		fop(24) = .806
		fop(25) = .000
c***********************************************************
		zp100 = 1. - 2.*fop(n)	
		y11 = xmtor(1)+atnum(1)*eq*(xni(1)*ephia+bthet*gamion(n,1)*zp100)
 		y1 = y11/(xni(1)*xmas(1)*xnuc(1,2))	
	  endif
c********************************************************************************
		y22 = xmtor(2) + zbar2(n)*eq*(xni(2)*ephia + bthet*gamion(n,2))
		y2 = y22/(xni(2)*xmas(2)*xnuc(2,1))
	xx = xni(1)*xmas(1)/(xni(1)*xmas(1)+xni(2)*xmas(2)) 

	
	torq2(n) = y22

c	4-29-13 new approx************************************************************
	xddrag1(n) = y11/(xni(1)*xmas(1)*torv(n))
	xddrag2(n) = y22/(xni(2)*xmas(2)*torv(n))
	ndrag=3
		 
c	approximation #0  nudrag1=nudrag2, zero-order approx (vj-vk)_0
	if(ndrag.eq.0) then
	xnud0 = xx*xnuc(1,2)*(y1+y2)/torv(n)
      brack0 = (xnuc(1,2)*y1-xnud0*torv(n))/(xnuc(1,2)+xnud0)
	xnudtot1(n) = xnud0
	xnudtot2(n) = xnud0
	vtor1(n) = torv(n) + brack0
	vtor2(n) = torv(n)
	endif

c	approximation #1  vphi1=vphi2, zero order approx nudj_1 & nudk_1
	if(ndrag.eq.1) then 
	xnudtot1(n) = xnuc(1,2)*y1/torv(n)
	bnud1(n) = xnudtot1(n)/xnuc12(n)
	xnudtot2(n) = xnuc(2,1)*y2/torv(n)
	bnud2(n) = xnudtot2(n)/xnuc21(n)
     	vtor1(n) = torv(n)
	vtor2(n) = torv(n)
	endif
cXXXXXXXX****THIS IS NUDRAG INFERENCE FROM EXP USING PERT THEORY****XXXXXXXXXXXXX
c	approximation #2  nudrag1=nudrag2, first-order approx nud_2 and (vj-vk)_2
	if(ndrag.eq.2) then
	xnud0 = xx*xnuc(1,2)*(y1+y2)/torv(n)
	
      brack0 = (xnuc(1,2)*y1-xnud0*torv(n))/(xnuc(1,2)+xnud0)
	xnud2p = xnud0/(1.+xx*brack0/torv(n))
	xnudtot1(n) = xnud2p
	xnudtot2(n) = xnud2p
	brack2 = (xnuc(1,2)*y1-xnud2p*torv(n))/(xnuc(1,2)+xnud2p)
 	vtor1(n) = torv(n) + brack2
	vtor2(n) = torv(n)
	endif	

c	approximation #3, extend approx #1 to use (vj-vk)_0  APRIL 2013
	if(ndrag.eq.3) then

 	ZR = 1.0
	DELV1(N) = 0.
	XNUDTOT1(N) = (Y11+Y22)/(((xni(1)*xmas(1)+
     2	ZR*xni(2)*xmas(2))*TORV(N))	+ XNI(1)*XMAS(1)*DELV1(N))
	XNUDZERO(N) = XNUDTOT1(N)
	XNUDTOT2(N) = XNUDTOT1(N)
	Q1 = 1. - XNUC(1,2)/(XNUC(1,2)+XNUDTOT1(N))
	Q2 = 1. - XNUC(2,1)/(XNUC(2,1)+XNUDTOT2(N))
	RACK = 1.- XNUC(2,1)*XNUC(1,2)/((XNUC(1,2)+
     2		XNUDTOT1(N))*(XNUC(2,1)+XNUDTOT2(N)))
	DELV1(N) = (Q2*Y11/(XNI(1)*XMAS(1)*(XNUC(1,2)+XNUDTOT1(N))) -
     2			Q1*Y22/(XNI(2)*XMAS(2)*(XNUC(2,1)+XNUDTOT2(N))))/RACK

      DELV0(N) = DELV1(N)

 	XNUDTOT1(N) = (Y11+Y22)/(((xni(1)*xmas(1)+
     2	ZR*xni(2)*xmas(2))*TORV(N)) + XNI(1)*XMAS(1)*DELV1(N))
	XNUDTOT2(N) = XNUDTOT1(N)
	Q1 = 1. - XNUC(1,2)/(XNUC(1,2)+XNUDTOT1(N))
	Q2 = 1. - XNUC(2,1)/(XNUC(2,1)+XNUDTOT2(N))
	BRACK = 1.- XNUC(2,1)*XNUC(1,2)/((XNUC(1,2)+
     2		XNUDTOT1(N))*(XNUC(2,1)+XNUDTOT2(N)))
	DELV1(N) = (Q2*Y11/(XNI(1)*XMAS(1)*(XNUC(1,2)+XNUDTOT1(N))) -
     2			Q1*Y22/(XNI(2)*XMAS(2)*(XNUC(2,1)+XNUDTOT2(N))))/RACK
	XNUDZERO(N) = XNUDTOT1(N)
C	separate drag frequencies for the two species from individ mom balances
	XNUDTOT1s(N) = (Y11-XNI(1)*XMAS(1)*XNUC(1,2)*DELV1(N))/
	2				(XNI(1)*XMAS(1)*(TORV(N)+DELV1(N)))
	XNUDTOT2s(N) = (Y22+XNI(2)*XMAS(2)*XNUC(2,1)*DELV1(N))/
	2				(XNI(2)*XMAS(2)*TORV(N))


c	then uses sep tor eq to evaluate xnud_main and xnud_imp from respect tor eq using (vphi_main - vphi_imp)_0  
C	xnudtot1(n) = xnuc(1,2)*((1.+y1/torv(n))/(1.+brack0/torv(n)) - 1.)
	
C	xnudtot2(n) = xnuc(2,1)*(y2/torv(n)+brack0/torv(n))
c	XNUDTOT2(N) = (Y22 + XNI(2)*XMAS(2)*XNUC(2,1)*BRACK0)/
c     1				(XNI(2)*XMAS(2)*TORV(N))	 
c	bnud2(n) = xnudtot2(n)/xnuc21(n)
c	XNUDTOT1(N) = ((xni(1)*xmas(1)+xni(2)*xmas(2))*XNUDZERO(N)-
c     1	XNI(2)*XMAS(2)*XNUDTOT2(N))/(XNI(1)*XMAS(1))
c	xnudtot1(n) = xnudzero(n) 
c	bnud1(n) = xnudtot1(n)/xnuc12(n)
 
c	then sets vphi_main = vphi_impex + (vphi_main - vphi_imp)_0 
	vtor1(n) = torv(n) + delv1(n)
	vtor2(n) = torv(n)
	endif

c	approximation #4, extend approx #1 to use (vj-vk)_1
	if(ndrag.eq.4) then 
	xnudtot1(n) = xnuc(1,2)*y1/torv(n)
	xnudtot2(n) = xnuc(2,1)*y2/torv(n)
	Rack = (y1-xnudtot1(n)*torv(n)/xnuc(1,2))/
	1							(1.+xnudtot1(n)/xnuc(1,2))
	xnudtot1(n) = xnuc(1,2)*((torv(n)+y1)/(torv(n)+rack)-1.)
	xnudtot2(n) = xnuc(2,1)*((torv(n)+rack+y1)/torv(n)-1.) 
	vtor1(n) = torv(n) + rack
	vtor2(n) = torv(n) 
	endif

c	new approx 4-29-13
c	xnudtot1(n) = xddrag1(n)
c	xnudtot2(n) = xddrag2(n)
CXXXXXXXXXXXX*********END NUDRAG INFERENCE FROM EXP*************XXXXXXXXXXXXXXXXXXXXXXX
	
c	**********************6/25/07 toroidal inference*********************************
c	this option solves the 2 tor mom eqs for a common nudrag and vphi_main, using measured vphi_imp
c	resulting quadratic eq has 2 solutions for nudrag and 2 corresponding solutions for vphi_main
	 xx = xni(1)*xmas(1)/(xni(1)*xmas(1)+xni(2)*xmas(2)) 
	y1 = xmtor(1) + atnum(1)*eq*(xni(1)*ephia + bthet*gamion(n,1))
 	y2 = xmtor(2) + zbar2(n)*eq*(xni(2)*ephia + bthet*gamion(n,2))
	gkj = 10.0
	a =	xni(2)*xmas(2)*torv(n)*gkj
	b = xni(2)*xmas(2)*torv(n)*gkj - (xmtor(1) + xmtor(2)) + 
	1	xni(1)*xmas(1)*(xnuc(1,2)*torv(n)) + y1
	c = -1.*xnuc(1,2)*(xmtor(1)+xmtor(2))
	xnudinf1(n) = -0.5*(b/a)*(1.+sqrt(1.-4.*a*c/(b**2)))*xnuc(1,2) 
	xnudinf2(n) = -0.5*(b/a)*(1.-sqrt(1.-4.*a*c/(b**2)))*xnuc(1,2) 
	vftor2(n) = torv(n)
	vftor1(n) = torv(n)/(1.+xnudinf1(n)/xnuc(1,2)) + 
     1	y1/(xni(2)*xmas(2)*(xnuc(2,1)+xnudinf1(n)))
	vftor0(n) = torv(n)/(1.+xnudinf2(n)/xnuc(1,2)) + 
     1	y1/(xni(2)*xmas(2)*(xnuc(2,1)+xnudinf2(n)))
c	**********end toroidal********************************************************
	 
	xx99 = xnustar(1,2)
	yy99 = xnuatomstar(1)


c	poloidal rotation & gyroviscosity
	vphia(1) = vtor1(n)
	vphia(2) = vtor2(n)
c	temp fix ***********%%%%%%%%%%&&&&&&&&&&&***********
C	if(n.eq.19) then
C	xx=44.
C	endif
C	call poloidal(n) 
C	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
C	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
C	velthet1(n) = vtheta(1)*vth(1)*fp
C	velthet2(n) = vtheta(2)*vth(2)*fp
C	xnudragyro1(n) = xnudrag(1)
C	xnudragyro2(n) = xnudrag(2)
	ynudrag1(n) = xnudtot1(n)
	ynudrag2(n) = xnudtot2(n)
	VT1 = VELTHET1(N)
	VT2 = VELTHET2(N)
c	*************repeat poloidal calc w/ v_imp = vthexp*****added 6/28/07*************
c	ipolopt =1
c	call poloidal(n)
c	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
c	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
c	vvelthet1(n) = vtheta(1)*vth(1)*fp
c	vvelthet2(n) = vtheta(2)*vth(2)*fp
c	xxnudragyro1(n) = xnudrag(1)
c	xxnudragyro2(n) = xnudrag(2)
c937	ipolopt = 0
 

c	*******************6/26/07 poloidal inference ******************************************
c	this option solves the 2 pol mom bal eqs for f_main and vpol_main, using the collisional f_imp
c	and the measured vpol_imp

c	6/26/07 version neglects poloidal asymmetries
		
	ep = rhor(n)*aminor/rmajor
	xnu11star = (qedge(n)*rmajor*xnuc(1,1))/vth(1)
	xnu22star = (qedge(n)*rmajor*xnuc(2,2))/vth(2)
	f2 = xnu22star/((1.+xnu22star)*((ep**1.5)+xnu22star))
	f1 = xnu11star/((1.+xnu11star)*((ep**1.5)+xnu11star))  
	fpolimp_shaing(n) = f2
	vpol_imp(n) = vthexp(n)
	xnuionstar1 = xnuioni(n)*qedge(n)*rmajor/vth(1)
	xnuion2 = 0.0
	xnuionstar2 = xnuion2*qedge(n)*rmajor/vth(2)
	xnu12star =  xnuc(1,2)*qedge(n)*rmajor/vth(1)
	xnu21star =  xnuc(2,1)*qedge(n)*rmajor/vth(2)
 	
	rho_main = vth(1)/(eq*bphi/xmas(1))
	rho_imp = vth(2)/(atnum(2)*eq*bphi/xmas(2))
	yy = qedge(n)*rmajor*vrad1(n)/(rho_main*vpol_imp(n))
	yz = qedge(n)*rmajor*vrad2(n)/(rho_imp*vpol_imp(n))
	
c*	********calculate the H-S Kj factors for the heat flux terms************************************* 
	 
	alph = (atnum(2)**2)*xni(2)/xni(1) 
	ym00iB = 0.53 + alph
	ym00iP = 3.54
	Di =2.23 + 5.32*alph +2.4*(alph**2)
	ym00iPS =(3.02 + 4.25*alph)/Di
	yk01iB = 0.71 + alph
	yk01iP = 10.63
	yk01iPS = (12.42 +20.13*alph)/Di
	ynustari = xnu11star/(ep**1.5)

	ym00zB = 0.53 + 1./alph
	ym00zP = 3.54
	Dz =2.23 + 5.32/alph +2.4/(alph**2)
	ym00zPS =(3.02 + 4.25/alph)/Dz
	yk01zB = 0.71 + 1./alph
	yk01zP = 10.63
	yk01zPS = (12.42 +20.13/alph)/Dz 
	ynustarz = xnu22star/(ep**1.5)

	yg = 1.46*sqrt(ep)-0.46*(ep**1.5)
	yg = yg/(1. - yg)

	ym00i = yg*ym00iB/((1.+2.92*ynustari*ym00iB/ym00iP)*
	1		(1.+ynustari*ym00iP*(ep**1.5)/(6.*ym00iPS)))
	ym00z = yg*ym00zB/((1.+2.92*ynustarz*ym00zB/ym00zP)*
     1		(1.+ynustarz*ym00zP*(ep**1.5)/(6.*ym00zPS)))
	yk01i = yg*yk01iB/((1.+2.92*ynustari*yk01iB/yk01iP)*
 	1		(1.+ynustari*yk01iP*(ep**1.5)/(6.*yk01iPS)))
 	yk01z = yg*yk01zB/((1.+2.92*ynustarz*yk01zB/yk01zP)*
     1		(1.+ynustarz*yk01zP*(ep**1.5)/(6.*yk01zPS)))
	ym01i = 2.5*ym00i - yk01i 
	ym01z = 2.5*ym00z - yk01z
	
      xki = ym01i/ym00i
	xkz = ym01z/ym00z
	zkHS1(n) = xki
	zkHS2(n) = xkz
c	****************************************************************************************************	 

c	*****extend the vrad source terms to included the heat flux viscosity contribution***********
	cvisc = 1. 
	s1 = eq*vrad1(n)*bphi/xmas(1)
c     	s1 = s1 + cvisc*(vth(1)*qedge(n)*f1/(rmajor))*
c     2	(xki*xk*xti(n)/(exlti(n)*eq*bphi))
	s2 = atnum(2)*eq*vrad2(n)*bphi/xmas(2)
c	s2 = s2	+ cvisc*(vth(2)*qedge(n)*f2/(rmajor))*
c     2	(xkz*xk*xti(n)/(exlti(n)*eq*bphi))
c	************************************************************************************************

c	option #1  uses collision fpol_imp=1/nukk for imp,& solves for fpol_main and vpol_main
c	modified 7/1/07 to use shaing fol_imp 

      xx = cvisc*(qedge(n)*vth(2)/(rmajor*xnuc(2,1)))*f2
	vpol_main(n)=vpol_imp(n)*(1.+xnuion2/xnuc(2,1)+xx) + s2/xnuc(2,1)
	
	fpol_main(n)= (xnuc(1,2)*(vpol_imp(n)-vpol_main(n))
	1     - xnuioni(n)*vpol_main(n) - s1)/
	2		((cvisc*qedge(n)*vth(1)*vpol_main(n))/rmajor)	
	fpol_col(n) = 1./xnu11star

	fpol_shaing(n)= xnu11star/((1.+xnu11star)*((ep**1.5)+
     1	xnu11star))
	fpol_imp(n) = f2	 
c	option #2  uses same fpol for both species, solves for fpol and vpol_main
 
	a = (cvisc**2)*vth(1)*vth(2)*(qedge(n)**2)*vthexp(n)/
	1		(xnuc(2,1)*(rmajor**2))
	b =cvisc*(vth(1)*qedge(n)/rmajor)*((1.+xnion2/xnuc(2,1))*vthexp(n)
	1  + s2/xnuc(2,1)) 
	2+cvisc*(vth(2)*qedge(n)/(rmajor*xnuc(2,1)))*(xnuc(1,2)+xnuioni(n))
	3	*vthexp(n)
	c = (xnuc(1,2)+xnuioni(n))*((1.+xnion2/xnuc(2,1))*vthexp(n)
	1	+ s2/xnuc(2,1))
	2	-xnuc(1,2)*vthexp(n) + s1
	fpol1(n) = -0.5*(b/a)*(1.+sqrt(1.-4.*a*c/(b**2)))
	fpol2(n) = -0.5*(b/a)*(1.-sqrt(1.-4.*a*c/(b**2)))
	vpol3(n) = (xnuc(1,2)*vthexp(n)-s1)/
	1 (xnuc(1,2) + xnuioni(n) + cvisc*vth(1)*qedge(n)*fpol1(n)/rmajor)
	vpol4(n) = (xnuc(1,2)*vthexp(n)-s1)/
	1 (xnuc(1,2) + xnuioni(n) + cvisc*vth(1)*qedge(n)*fpol2(n)/rmajor)

c	option #3  uses fpol_shaing and solves for both vpols
	f1 = fpol_shaing(n)
	xnu22star = (qedge(n)*rmajor*xnuc(2,2))/vth(2)
	f2 = xnu22star/((1.+xnu22star)*((ep**1.5)+xnu22star))
	fpolimp_shaing(n) = f2
	a1 = xnuc(1,2)+xnuioni(n)+cvisc*(vth(1)*qedge(n)*f1/rmajor)
	a2 = xnuc(2,2)+xnuion2   +cvisc*(vth(2)*qedge(n)*f2/rmajor) 
	b1 = -1.*xnuc(1,2)
	b2 = -1.*xnuc(2,1)
	s1 = -1.*s1
	s2 = -1.*s2
	vpol5(n) = (a2*s1-b1*s2)/(a1*a2-b1*b2)
	vpol6(n) = (s2-b2*vpol5(n))/a2



c	**************end poloidal inference********************************
c*************calculate Vtheta from eq (6) of structure paper --spring 2013********
	vth(1) = sqrt(2.*xk*ti(n)/xmas(1)) 
	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
	se = sqrt((1.+ elong**2)/2.)
	ep = se*rhor(n)*aminor/rmajor
	xnu11star = (qedge(n)*rmajor*xnuc(1,1))/vth(1)
	xnu22star = (qedge(n)*rmajor*xnuc(2,2))/vth(2)
	f2 = xnu22star/((1.+xnu22star)*((ep**1.5)+xnu22star))
	f1 = xnu11star/((1.+xnu11star)*((ep**1.5)+xnu11star))  
 
      xnuvisc1 = (qedge(n)**2)*f1
	xnuvisc2 = (qedge(n)**2)*f2
	xnuions1 = xnuioni(n)*qedge(n)*rmajor/vth(1)
	xnucxs1 = xnuatim(n)*qedge(n)*rmajor/vth(1)
	xnuions2 = xnuionz(1)*qedge(n)*rmajor/vth(2)
 	xnucxs2 = xnuatim(n)*qedge(n)*rmajor/vth(2)

	xnutheta1= xnuvisc1 + xnuions1 + xnucxs1 +
     2	 xnuc(1,2)*qedge(n)*rmajor/vth(1)
	xnutheta2= xnuvisc2 + xnuions2 + xnucx2 +
	2     xnuc(2,1)*qedge(n)*rmajor/vth(2)
	zv = 1.0

	drive1v(n) = xnuvisc1*xki*ti(n)/(3.*(bphi**2)*exlti(n)) 
	drive1(n) = zv*drive1v(n) +
     2(eq*gamion(n,1)*ziol*qedge(n)*rmajor)/(yni(n,1)*xmas(1)*vth(1))
	drive2v(n) = xnuvisc2*xkz*ti(n)/(3.*zbar2(n)*(bphi**2)*exlti(n))  
      drive2(n) = zv*drive2v(n)+zbar2(n)*eq*gamion(n,2)*qedge(n)*rmajor/
     3(yni(n,2)*xmas(2)*vth(2))
     
	xnu12star = (qedge(n)*rmajor*xnuc(1,2))/vth(1)
	xnu21star = (qedge(n)*rmajor*xnuc(2,1))/vth(2)
 	bracc = 1.-xnuc12star*xnuc21star/(xnutheta1*xnutheta2)
	2 
	vpol61(n) = -1.*bphi*(drive1(n) + 
     2		xnuc12star*drive2(n)/xnutheta2)/bracc
	vpol62(n) = -1.*bphi*(drive2(n) + 
     2    	xnuc21star*drive1(n)/xnutheta1)/bracc



1578	continue 

c**************temporary input nudrag*************************************

c	 ynudrag1(1) =.311E+05
c	 ynudrag1(2) =.320E+05
c	 ynudrag1(3) =.427E+03
c	 ynudrag1(4) =.450E+03
c	 ynudrag1(5) =.480E+03
c	 ynudrag1(6) =.513E+03
c	 ynudrag1(7) =.549E+03
c	 ynudrag1(8) =.591E+03
c	 ynudrag1(9) =.642E+03
c	 ynudrag1(10) =.695E+03
c	 ynudrag1(11) =.761E+03
c	 ynudrag1(12) =.837E+03
c	 ynudrag1(13) =.928E+03
c	 ynudrag1(14) =.105E+04
c	 ynudrag1(15) =.120E+04
c	 ynudrag1(16) =.138E+04
c	 ynudrag1(17) =.163E+04
c	 ynudrag1(18) =.193E+04
c	 ynudrag1(19) =.238E+04
c	 ynudrag1(20) =.299E+04
c	 ynudrag1(21) =.360e+04
c	 ynudrag1(22) =.447E+04
c	 ynudrag1(23) =.559E+04
c	 ynudrag1(24) =.107E+05
c	 ynudrag1(25) =.247E+06

c*************************************************************************





	

C	***************INFER CHI's FROM EXP N,T,GSCL & CALC. Q, GAM*****************
1579	do 1589 n = 1,25
 	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
	XCHIE(N) = EXLTE(N)*((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
     1					1.5*GAMEL/EXNE(N))	 
	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))

	XCHII(N) = EXLTI(N)*((GAMHEATI(N)/(XNION*XK*XTI(N)))-
     1					2.5*GAMION(N,1)/XNION)
	XCHII15(N) = EXLTI(N)*((GAMHEATI(N)/(XNION*XK*XTI(N)))-
     1					1.5*GAMION(N,1)/XNION)
	xchieorb(n) =  EXLTE(N)*((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
     1					1.5*GAMEL*(1.-forbl(n))/EXNE(N))	
c	chii w/ion orbit loss
	orbf = 1. - forbl(n)
	orbe = 1. - eorbl(n)
	XCHIIorb(N) = EXLTI(N)*(((GAMHEATI(N)*orbe)/(XNION*XK*XTI(N)))-
     1					2.5*GAMION(N,1)*orbf/XNION)			
	XCHIIorb15(N) = EXLTI(N)*(((GAMHEATI(N)*orbe)/(XNION*XK*XTI(N)))-
     1					1.5*GAMION(N,1)*orbf/XNION)

c    chii w/ion orbit and X-LOSS
	orbfx = orbf*(1. - fx(n))
	orbex = orbe*(1. - ex(n))
	XCHIIorbX(N) = EXLTI(N)*(((GAMHEATI(N)*orbex)/(XNION*XK*XTI(N)))-
     1					2.5*GAMION(N,1)*orbfx/XNION)
	XCHIIorbX15(N) = EXLTI(N)*(((GAMHEATI(N)*orbex)/(XNION*XK*XTI(N)))
     1				-	1.5*GAMION(N,1)*orbfx/XNION)

c	ion particle fluxes reduced by ion orbit and X loss
	  gamionorb(n) = gamion(n,1)*(1.-forbl(n)) 
	  cdionorb(n) = eq*gamion(n,1)*forbl(n)
	  gamheatiorb(n) = gamheati(n)*(1. - eorbl(n))
	  gamionorbX(n) =  gamion(n,1)*(1.-forbl(n))*(1. - fx(n)) 
	  gamheatiorbX(n) = gamheati(n)*(1. - eorbl(n))*(1. - ex(n))

c	chii with X-transport
	xXCHII(N) = EXLTI(N)*((GAMHEATIx(N)/(XNION*XK*XTI(N)))-
     1					2.5*GAMIONx(N)/XNION)
	xXCHII15(N) = EXLTI(N)*((GAMHEATIx(N)/(XNION*XK*XTI(N)))-
     1					1.5*GAMIONx(N)/XNION)
	xXCHIIorb(N) = EXLTI(N)*((GAMHEATIx(N)*orbe/(XNION*XK*XTI(N)))-
     1					2.5*GAMIONx(N)*orbf/XNION)
	xXCHIIorb15(N) = EXLTI(N)*((GAMHEATIx(N)*orbe/(XNION*XK*XTI(N)))-
     1					1.5*GAMIONx(N)*orbf/XNION)


c	conductive heat fluxes
	qcond25(n) = gamheati(n) - 2.5*gamion(n,1)*xk*xti(n) 
	qcond15(n) = gamheati(n) - 1.5*gamion(n,1)*xk*xti(n)
	qcond25orb(n) =  gamheati(n)*orbe - 2.5*gamion(n,1)*orbf*xk*xti(n) 
	qcond15orb(n) =  gamheati(n)*orbe - 1.5*gamion(n,1)*orbf*xk*xti(n) 
	qcond25orbx(n)=gamheati(n)*orbex - 2.5*gamion(n,1)*orbfx*xk*xti(n) 
	qcond15orbx(n)=gamheati(n)*orbex - 1.5*gamion(n,1)*orbfx*xk*xti(n) 
	qcond25xtran(n) =  gamheatix(n) - 2.5*gamionx(n)*xk*xti(n) 
	qcond15xtran(n) =  gamheatix(n) - 1.5*gamionx(n)*xk*xti(n) 
	qcond25xtranorb(n)=gamheatix(n)*orbe-2.5*gamionx(n)*orbf*xk*xti(n) 
	qcond15xtranorb(n)=gamheatix(n)*orbe-1.5*gamionx(n)*orbf*xk*xti(n) 

c	xtrans ion particle fluxes reduced by ion orbit loss
	xgamionorb(n) = gamionx(n)*orbf 
	xheatXtranorb(n) =gamheatix(n)*orbe
c	***********Pinch Velocities & Diffusion Coefficients*********************
	ipolopt = 1
c	 ipolopt = 1 uses vthet2 = vexp, =0 calculates vthet2

    	call poloidal(n) 
	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
	velthet1(n) = vtheta(1)*vth(1)*fp
	velthet2(n) = vtheta(2)*vth(2)*fp
 	xnudragyro1(n) = xnudrag(1)
	xnudragyro2(n) = xnudrag(2)
      erada(n) = erex(n)
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
	zdrag =	(xmas(1)*ynudrag1(n)*((erada(n)/bthet)+velthet2(n)/fp))/
	1	 (eq*bthet)
	zfrict = 0.0
	vp2 = zbeam - ephia/bthet+zdrag+zfrict 
c	if(vp2.gt.0.0) vp2 = 0.0
	vpinch2(n) = vp2
		zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
	zdrag =	(xmas(1)*ynudrag1(n)*((erada(n)/bthet)+velthet2(n)/fp))/
	1	 (eq*bthet)
	zfrict = (xmas(1)*(xnuc12(n)+ynudrag1(n))*
	1		(velthet1(n)-velthet2(n))/fp)/(eq*bthet)
	vp3 = zbeam - ephia/bthet+zdrag+zfrict
c	if(vp3.gt.0.0) vp3 = 0.0 
	vpinch3(n) = vp3
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	eradterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)))/
	1			(eq*bthet)
	delvt2 = velthet2(n)-vthexp(n)							
c*********correction to calculated D vel to add difference between exp
c****and calculated C velocity to calculated D velocity***************
	vt1 = velthet1(n)-delvt2
	vth1cor(n) = vt1
c**********************************************************************
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*
	1((velthet1(n)/fp)))/ (eq*bthet)
		
	zfrict = 0.
	zdrag = eradterm+vtheterm
 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp4 = zbeam-ephia/bthet+vphiterm+eradterm+vtheterm+zfrict
c	if(vp4.gt.0.0) vp4 = 0.0
	vpinch4(n) = vp4
	xnuc120(n) = xnuc12(n)
c******vp4 is evaluated using VphiCexp and Vthet1=Vthet3 from rad mom bal as of 2/21/12
c	this should be preferred Vpinch evaluation***************************************
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	eradterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)))/
	1			(eq*bthet)
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((velthet1(n)/fp)))/
     1			(eq*bthet)

	zfrict = (xmas(1)*(xnuc12(n)+ynudrag1(n))*
 	1		(velthet1(n)-velthet2(n))/fp)/(eq*bthet)
	 
c*******************************************************************
	zdrag = eradterm+vtheterm
 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)

	vp5 = zbeam-ephia/bthet+vphiterm+eradterm+vtheterm+zfrict

c*********code is set up to use vp5 and vpinch, course of least effort

	vpinch5(n) = vp5
c************************************************************************
	EXD51(N) = ((GAMION(N,1)/XNION)-VPINCH5(N))*
     1	xlne(n)*exlti(n)/(XLNE(N) + EXLTI(N))
	 EXD41(N) = ((GAMION(N,1)/XNION)-VPINCH4(N))*
	1	xlne(n)*exlti(n)/(XLNE(N) + EXLTI(N))
c*****EXC41 evaluated using preferred vp4 eval with vthet1 from rad mom bal 2/21/12
	 EXD53(N) = ((GAMION(N,1)/XNION))*(XLNE(N))
	   EXD21(N) = ((GAMION(N,1)/XNION)-VPINCH2(N))*
     1	xlne(n)*exlti(n)/(XLNE(N) + EXLTI(N))
	 EXD31(N) = ((GAMION(N,1)/XNION)-VPINCH3(N))*
     1	xlne(n)*exlti(n)/(XLNE(N) + EXLTI(N))

	 EXD22(N) = ((GAMION(N,1)/XNION)-VPINCH2(N))*(XLNE(N))
	 EXD23(N) = ((GAMION(N,1)/XNION))*(XLNE(N))
 	
	exd11(n) = xmas(1)*(xnuc12(n)+ynudrag1(n))*xti(n)/
	1			(eq*(bthet**2))


c	************************************************************


	

c	sensitivity studies
c	without convection but with spatially varying heat flux
	XCHIE2(N) = EXLTE(N)*(GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))
 	XCHII2(N) = EXLTI(N)*(GAMHEATI(N)/(XNION*XK*XTI(N)))
c	without convection and with flat heat flux
	XCHIE3(N) = EXLTE(N)*(GAMHEATE(25)/(EXNE(N)*XK*XTE(N)))	
	XCHII3(N) = EXLTI(N)*(GAMHEATI(25)/(XNION*XK*XTI(N)))
c	with ion orbit loss
 	
	XCHII2(N) = (EXLTI(N)/(XNION*XK*XTI(N)))*
     1	(GAMHEATI(N)*(1.-EORBL(N))-	
     2	2.5*GAMION(N,1)*(1.-FORBL(N))*XK*XTI(N))
 	XCHII3(N) = (EXLTI(N)/(XNION*XK*XTI(N)))*
     1	(GAMHEATI(N)*(1.-EORBL(N))-	
     2	1.5*GAMION(N,1)*(1.-FORBL(N))*XK*XTI(N))


c	convective fractions
	 econve(n) = 2.5*gamel*xk*xte(n)/gamheate(n) 
	econvi(n) = 2.5*(gamion(n,1)+gamion(n,2))*xk*xti(n)/gamheati(n) 
c	calculate theoretical chi's

1589	continue 
	do 7100 nmesh =1,25 
	n= 26-nmesh
	EQ = 1.6E-19
	XK = 1.6E-19
	EP0 = 8.854E-12
	AM = (aminor*SE)
 	ep = aminor*rhor(n)/rmajor 
c	untrapped (fc) and trappled (ftr) particle fractions, with orbit squeezing
	fc = (((1.-ep)**2)/sqrt(1.-(ep**2)))/(1.+1.46*sqrt(ep)+0.2*ep)
	ftr = 1. - fc
	if(shearho(n).eq.0.0) shearho(n) = 1.e-6
	ss =aminor*shearho(n)*rhor(n)/qedge(n)
	ep = aminor*rhor(n)/rmajor 
	bfield = abs(bphi) 
	OMI =EQ*BFIELD/XMAS(1)
	CSOUND = SQRT(XK*XTE(N)/XMAS(1))
 	rhot = csound/omi
 	XVTHI = SQRT(2.*XK*XTI(N)/XMAS(1))
	para = 1./(bthet*XvthI)
	sheare(n) =	abs(1.-(rhot/abs(fp))*dEdredge(n)*para)
	ftros = ftr*sqrt(abs(sheare(n)))
	fcos = 1. - ftros
c	orbit squeezing correction to trapped particle fraction********************
	ftr = ftros
	fc = fcos 
c	****************************************************************************
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
	xnuc21(n) = xnuc(2,1)
	xnuc22(n) = xnuc(2,2)
	xnuc11(n) = xnuc(1,1)

	XVTHI = SQRT(2.*XK*XTI(N)/XMAS(1))
	vthi = xvthi
	XVTHIMP = SQRT(2.*XK*XTI(N)/XMAS(2))
	EMASS = 9.1E-31
	xvthe = sqrt(2.*xk*xte(n)/emass)
	CLIGHT = 3.E8 
	CSE = SQRT(2.*XK*XTE(N)/EMASS) 
	XNELECTRON = XNI(1)*(ATNUM(1)**2) + XNI(2)*(zbar2(n)**2) 
	xnelectron = exne(n)
	XNUEI =EXNE(N)/(6.4E14*((1.E-3*XTE(N))**1.5)) 
	XNUEIAST(N) = XNUEI*ABS(QSAFE)*RMAJOR/XVTHE

	
C	NEOCLASSICAL CHI FOR IONS
c	simple neoclassical chi 
 
	CHINC(n) = (((RHOT)*bfield/bthet)**2)*XNUc(1,2)*(EP**0.5) 
C	CHANG-HINTON CHI
	ALFA = XNi2*(zbar2(n)**2)/(XNi1*(ATNUM(1)**2))
	qa = qedge(25)
	XMUii=(XNUc(1,1)*Qedge(n)*RMAJOR/(XvthI*(EP**1.5)))*(1.+1.54*alfa)
	
	dp = 0.
6100	G1 = (1. + 1.5*((EP**2)+ep*dp)+.375*(ep**3)*dp)/(1.+.5*ep*dp)
	G2 =SQRT(1.-(EP**2))*(1.+0.5*ep*dp)/(1.+(dp/ep)*(sqrt(1.-ep**2)
	1	-1))
	A1 =(0.66*(1.+1.54*ALFA)+(1.88*SQRT(EP)-1.54*EP)*(1.+3.75*ALFA))/
	1	(1.+1.03*SQRT(XMUii)+0.31*XMUii)
	A2 =0.59*XMUii*EP*(1.+1.33*ALFA*(1.+0.6*ALFA)/(1.+1.79*ALFA))/
     1	(1.+0.74*XMUii*(EP**1.5))	 
	
	betap = 2.*eXne(n)*xk*Xti(n)/((bthet**2)/(2.*1.257e-6)) 
	CHICH(n) = CHINC(n)*(xnuc(1,1)/xnuc(1,2))*(A1*G1+A2*(G1-G2))
	chich8(n) = chich(n)
	if(dp.eq.0) then
		chich0 = chich(n)
		dp = -1.*ep*(betap+0.5*log(1.65+0.89*(qa-1.)))
  		goto 6100
	endif
	
	
c	orbit squeezing
	

	

6110	continue 
c	if(ioptshear.eq.0) sheare(n) =	1.-(rhoti(1)/abs(fp))*dlnEdr*para


	
c6125	chinc(n) = chinc(n)/(abs(sheare(n))**1.5)
 
	chichos(n) = chich(n)/(abs(sheare(n))**1.5)
	
     
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
	chiet = rhot*xte(n)/abs(bfield)
	tauw = xti(n)/(zeff*xte(n))
	xx = ((1./exlti(n)) - (0.6667/xlne(n)))*(tauw/(rmajor*(1.-ftr))) -
     2	(1./8.)*((2/rmajor)-(1./xlne(n)))**2/((1.-ftr)**2) -
     3	(20./9.)*(tauw**2/(rmajor**2))
 	fweiland = 3.333*sqrt(xx)
    	fromanelli = 1.25*sqrt(1./(exlti(n)*rmajor))
	Ci = 0.014
	Ci1 = 0.00014  
	fhorton1 = Ci1*(qedge(n)**2)*(rmajor**1.5)/(exlti(n)**2.5)
	fhorton2 = Ci*qedge(n)*(rmajor**0.5)/(exlti(n)**1.5)
	if(ioptran.eq.1) then
	CHIION(N) = XCHII(N)
C	chiion(n) = chixpi
C		if(rhor(n).lt.pedrhoti) chiion(n) = chitop  
 	endif 
c	Weiland ITG as given by Kalupin, et al., NF,45,468,2005.

	tau =xti(n)/(zeff*xte(n))
	chiix = xte(n)*rhot/(bfield) 
	
       
c	etai threshold  F. Jenko et al., PoP,7,1904,2000 & 8,4096,2001
	rlcrit = 0.8*rmajor/xlne(n)
	xlcrit=(1.+xti(n)/(zeff*xte(n)))*
	1		(1.33+1.91*ss/qedge(n))*(1.-1.15*ep)
      if(xlcrit.gt.rlcrit) rlcrit = xlcrit
	xz = rmajor/exlti(n)
	crititg(n) = xz/rlcrit
	chietair0(n) = 0.0
	chietaiw0(n) = 0.0
	chietaih0(n) = 0.0
	if(xz.gt.rlcrit) then
		chietaiw0(n) = chiet*fweiland 
		chietair0(n) = chiet*fromanelli
		chietaih0(n) = chiet*fhorton1

	endif
c     ExB shear suppression (T.S.Hahm & K. H. Burrell, PoP,2,1648,1995)		     	 	 
     	omexb(n) = abs((dEdredge(n)/(aminor*bfield)) -
	1	(erex(n)/(rhor(n)*aminor))*(1.- ss))


	ysitg(n) = omexb(n)*sqrt(rmajor*exlti(n))/(0.3*csound)
	fsitg(n) = 1./(1. + (ysitg(n)**2))

	gamaxitg(n) = 0.3*csound/sqrt(exlti(n)*rmajor)

 	chietair1(n) = fsitg(n)*chietair0(n)
	chietaiw1(n) = fsitg(n)*chietaiw0(n)
	chietaih1(n) = fsitg(n)*chietaih0(n)
	chietai2h1(n) = fsitg(n)*chietai2h0(n) 

c	magnetic shear suppression ala Zolotukhin, 
	chietair2(n) = chietair1(n)/(ss**1.8)
   	chietaiw2(n) = chietaiw1(n)/(ss**1.8)
	chietaih2(n) = chietaih1(n)/(ss**1.8)
	chietai2h2(n) = chietai2h1(n)/(ss**1.8)

c	WEILAND ITG/TEM MODEL
C	critical etai
	tau = zeff*xte(n)/xti(n)
	RHOTS = CSOUND/OMII(1)
	RHOTE = 3.37E-6*SQRT(XTE(N)/BFIELD)
	xkrhoi = 0.3
	xkrhoe = 0.3
	dome = -1.*(xkrhoe/rhots)*xte(n)/(xlne(n)*bphi)

	domi = (xkrhoi/rhots)*xti(n)/(xlne(n)*bphi)

	
	epn = 2.*xlne(n)/rmajor 
	bracket = 5./3. -tau/4. + tau/(4.*epn) -
	1	(10./3.+ tau/4.-10./(9.*tau))*epn + 
     2    (5./3.+	tau/4.-10./(9.*tau))*(epn**2)
	etaithresh(n) = 0.6667 - 0.5*tau + epn*(tau/4.+10./(9.*tau)) +
     1	tau/(4.*epn) - (xkrhoi**2)*bracket/(2.*epn)	 
c	etaethresh(n) = 0.6667 - 0.5*tau + epn*(tau/4.+10./(9.*tau)) +
c     1	tau/(4.*epn) - (xkrhoe**2)*bracket/(2.*epn)	 
c	if(etaicriti(n).lt.0.0) etaicriti(n) = 0.0 
c	if(etaicrite(n).lt.0.0) etaicrite(n) = 0.0
 	gamitg(n)=dome*sqrt(epn/tau)*sqrt(etai(n)-etaithresh(n))/
	1										(1.+xkrhoi**2)
	gamtem(n)= gamitg(n)
 	xx =  etai(n)-etaithresh(n)
	if(xx.lt.0.0) then
	gamitg(n) = 0.0
	gamtem(n) = 0.0
	endif
	omreali(n) = 0.5*dome*(1. - epn*(1.+10./(3.*tau)) -
     1	(xkrhoi**2)*(1. + (1.+etai(n))/tau - epn -5.*epn/(3.*tau)))
	omreale(n) = 0.5*dome*(1. - epn*(1.+10./(3.*tau)) -
     1	(xkrhoe**2)*(1. + (1.+etai(n))/tau - epn -5.*epn/(3.*tau)))
	omdi = 3.*(xkrhoi/rhots)*xti(n)/(rmajor*bphi)
	omde = -3.*(xkrhoi/rhots)*xti(n)/(rmajor*bphi)
	chiitgx(n) = (etai(n)-0.667-10.*epn/(9.*tau))/etai(n)
	chiitgx(n) = chiitgx(n)*(gamitg(n)**3)*((rhots/xkrhoi)**2)
	chiitgx(n) = chiitgx(n)/((omreali(n)-5.*omdi/3.)**2+gamitg(n)**2)

	ghati = gamitg(n)/dome
	ghate = gamtem(n)/dome

	orhati = omreali(n)/dome 
	orhate = omreale(n)/dome
      capNi = (orhati**2 - ghati**2-(10./3.)*orhati*epn +
     1	(5./3.)*(epn**2))**2+ 4.*(ghati**2)*((orhati- (5./3.)*epn)**2)
	capNe = (orhate**2 - ghate**2-(10./3.)*orhate*epn +
     1	(5./3.)*(epn**2))**2+ 4.*(ghate**2)*((orhate- (5./3.)*epn)**2)


	ohati2 = orhati**2 + ghati**2 
	ohate2 = orhate**2 + ghate**2 

      delti =	ohati2*(ohati2*(epn-1.) + 
	1		orhati*epn*(14./3.-2.*etae(n)-(10./3.)*epn) +
     2		(5./3.)*(epn**2)*(2.*etae(n)-11./3.+(7./3.)*epn) -
     3		(5./(3.*tau))*(epn**2)*(1.+etae(n)-(5./3.)*epn)  )	+
     4		(50./(9.*tau))*orhati*(epn**3)*(1.-epn) -
     5		(25./(9.*tau))*(epn**4)*(7./3.-etae(n)-5.*epn/3.) 
	delti = delti/capNi 
	delte =	ohate2*(ohate2*(epn-1.) + 
     1		orhate*epn*(14./3.-2.*etae(n)-(10./3.)*epn) +
     2		(5./3.)*(epn**2)*(3.*etae(n)-8./3.+(2./3.)*epn)) -
     4		(50./9.)*orhate*(epn**3)*(1.-epn) +
     5		(25./9.)*(epn**4)*(7./3.-etae(n)-5.*epn/3.) 
	delte = delti/capNe 
	

c	orbit squeezing corrected trap particle fraction ftros
	bracketi=etai(n)-.667-(1.- ftros)*10.*epn/(9.*tau)-
	1	2.*ftros*delti/3.
	brackete=etae(n)-.667-2.*delte/3.
	chietaiw9(n) = (rhots**2)*(gamitg(n)**3)/((xkrhoi**2)*etai(n))
	chietaiw9(n) = chietaiw9(n)*bracketi/
	1		((omreali(n)-5.*omdi/3.)**2 + gamitg(n)**2)
	chietaiw91(n) = chietaiw9(n)*fsitg(n)/(ss**1.8) 
	chitemw9(n)=(rhots**2)*(gamtem(n)**3)*ftros/((xkrhoe**2)*etae(n))
	chitemw9(n) = chitemw9(n)*brackete/
	1	((omreale(n)-5.*omde/3.)**2 + gamtem(n)**2)

	chitemw91(n) = chitemw9(n)*fsitg(n)/(ss**1.8)
			
C	ETG-MODE CHI FOR ELECTRONS
	
c	
c	****exp etae**************
c	etae(n) = 1.43
c	ETAe(n) = XLTeM(NMESH)/XLNM(NMESH)

c	ETAE(n) =  ylnbarx/yltebarx
c	if(rhor(n).lt.pedrhon) etae(n) =  ylntop/xltetop

c	**************************
	eq = 1.6e-19 
	RHOS = CSOUND/OMII(1)
	RHOTE = 3.37E-6*SQRT(XTE(N)/BFIELD)
	vthe = sqrt(2.*xk*xte(n)/xme) 
	etae(n) = xlne(n)/exlte(n)
	
	eplasfreq = 56.4*sqrt(EXNE(N))
	ss =aminor*shearho(n)*rminor/qedge(n)
	shearm = ss
	ome = 56.4*sqrt(exne(n))
	dele = (3.e8/ome)
	chiwesson =0.13*((Clight/EPLASFREQ)**2)*xvthe*SS*ETAE(n)*
	1		(1.+ETAE(n))/(Qedge(n)*RMAJOR)

c	Bateman, Kritz et al. model of ETG MODE Transport(PPCF,48,A93,2006)
c	Quasi-linear coefficient W.Horton, P. Zhu, G.T.Joang, et al., PoP,7,1494,2000
	chihort = 0.06*((Clight/eplasfreq)**2)*xvthe/sqrt(rmajor*exlte(n))
	xlce =qedge(n)*rhote*rmajor/exlte(n)

	if(xlce.ge.dele) then
c	fitted numerical Cse = 0.03 from Horton,et al.  PoP,7,1494,2000	for Tore Supra
	Cse = 0.03 
	coefe =  Cse*(qedge(n)**2)*((rmajor/exlte(n))**1.5)*
     1		  ((rhote**2)*vthe/(xk*xte(n)))
	squarebracket = ((xk*xte(n)/exlte(n)) -
     2		  1.88*(ss*xk*xte(n)/(qedge(n)*rmajor))*
     3		  (1.+zeff*xte(n)/xti(n)))
	chihort = coefe*squarebracket 
	endif 
c	etg threshold  F. Jenko et al., PoP,7,1904,2000 & 8,4096,2001
	rlcrit = 0.8*rmajor/xlne(n)
	xlcrit=(1.+zeff*xte(n)/xti(n))*(1.33+1.91*ss/qedge(n))*(1.-1.5*ep)
      if(xlcrit.gt.rlcrit) rlcrit = xlcrit
	xx = rmajor/exlte(n)
	critetg(n) = xx/rlcrit
	chieetg(n) = 0.0
	chieetg1(n) = 0.0

c	if(xx.gt.rlcrit) then
		chieetg(n) = chihort
		chieetg1(n) = chiwesson
c	endif
	     	 	 
     	
c	maximum growth rate etg modes (W. Horton, P. Zhu, G. T. Hoang, et al.,
c	PoP, 7, 1494,2000.)

	gamaxetg(n) = vthe*1.414/sqrt(exlte(n)*rmajor)
c	if(omexb(n).gt.gamaxetg(n)) then
c	chieetg(n) = 0.0
c	chieetg1(n) = 0.0
c	endif

c     shear suppression (Hamaguchi & Horton,)	

		ysetg(n) = rhor(n)*aminor*sqrt(xmas(1)/(xk*xte(n)))*
	1			   abs(omexb(n)/(fp*ss))	
		fsetg(n) = 1./(1. + ysetg(n)**2) 
					
c	if(ys(n).gt.1.0) then
c	chieetg(n) = 0.0
c	chieetg1(n) = 0.0 
c      endif
c	partial shear suppression Hamaguchi
c	Wesson with & w/o Hamaguchi shear suppression
	chieetg0(n) = chieetg1(n)
	chieetg1(n) = fsetg(n)*chieetg1(n)
c	Horton with & w/o Hamaguchi shear suppression	 
	chieetg2(n) = chieetg(n) 
	chieetg(n) = fsetg(n)*chieetg(n)
	
 
      
C	TRAPPED ELECTRON MODE W/INTERPOLATION TO COLLISIONLESS DRIFT MODE
C		CHI FOR ELECTRONS 
	
	
	ylnylte = dngscl*dtegscl
c	if(rhor(n).lt.pedrhon) ylnylte = ylntop*xltetop
C	if(xlnm(n).ne.0.0)
c	wesson
      CHIEDW(N) = 2.5*(EP**1.5)*(CSOUND**2)*(RHOS**2)/(ylnylte*
     1	XNUEI*(1.+0.1/XNUEIAST(N)))
c	partial suppression
	 
	chiedw(n) = fsitg(n)*chiedw(n) 
c	Kalupin NF,45,468(2005)
	zz = xnuei/(ep*dome)
	chiedw0(n) = ftr*etae(n)*(xnuei/ep)*(rhos**2)/(1.+zz**2)
	chiedw1(n) = fsitg(n)*chiedw0(n)
	chiedw2(n) = chiedw1(n)/ss**1.8 
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
	AM = (aminor*SE)
	ep = aminor*rhor(n)/rmajor 
	fc = (((1.-ep)**2)/sqrt(1.-(ep**2)))/(1.+1.46*sqrt(ep)+0.2*ep)

	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(eXNe(n)) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(xte(n)**1.5)*X/Yz)
c	electron plasma frequency
	ome = 56.4*sqrt(exne(n))
c	ion & electron cyclotron frequencies
	omcyci = 0.96e8*zeff*bfield/2.
	omcyce = 1.76e11*bfield
	xlame = 1.2e16*(xte(n)**2)*17./(exne(n)*zeff*couloge)
	xnuestar = (rmajor*qedge(n))/((ep**1.5)*xlame)
	starnue(n) = xnuestar
	abar = aminor*elong/SE
 
	res = (9e9/1.96)*1.15e-14*zeff*couloge/(xte(n)**1.5) 
     
	beta =2.*EXne(n)*xk*Xte(n)*(2.*1.257e-6)/((bfield**2))
	XLPIM = (1./EXLTI(N))+(1./XLNE(N))

	alpha = 2.*1.257e-6*rmajor*(q95**2)*Xni1*xk*Xti(n)*xlpIm/
	1		(bfield**2)
	shearq = ss 
      chierb(n) = 10.0*res*sqrt(beta*3.343/9.1e-4)*((alpha/shearq)**1.5) 
c	gudzar pf,5, 3712,1993, Kalupin, et al., NF,45,468,2005.
	
	pi = 3.14159    
	chierb0(n)= ((2.*qedge(n)*rhote)**2)*xnuei*rmajor/xlne(n)

c	shear suppression of chierb
  		
	xlrb=6.28*qedge(n)*rmajor*sqrt(2.*exne(n)*(eq**2)*res*rhot/
     1		(xme*ome*sqrt(2.*rmajor*xlne(n))))
	taudrb = (xlrb**2)/chierb0(n)
 
      ysdrb(n) = 	omexb(n)*taudrb
	fsdrb(n) = 1./(1.+(ysdrb(n)**2))	 
	chierb(n) = fsdrb(n)*chierb0(n)	

c	paleoclassical for electrons (NF,45,1120,2005)
c	corrected 2/13/12

	xmunu = ((zeff+1.414 - log(2.414))/
	1 (zeff*(1 +sqrt(xnuestar)+xnuestar)))*sqrt(sheare(n))*(1.-fc)/fc
	const=(1.414+zeff)/(1.414+(13./4.)*zeff) + 	xmunu 
	dlndhatdnu = xmunu*(1.+0.5/sqrt(xnuestar))/const
	dlndhatdnu = dlndhatdnu/(1 +sqrt(xnuestar)+xnuestar) 
      dhat =(6.28**2)*rmajor*((1.+elong**2)/(2.*elong))*1.5*zeff*
     1		1.4e3*(couloge/17.)*const			
      bigL = xlame
	
	vthe = sqrt(2.*xk*xte(n)/xme)
	xnue = (vthe/xlame)
	dele = (3.e8/ome)
	xlmax = rmajor*qedge(n)*sqrt((3.1416*abar)/(dele*shearho(n)))

 	if(xlmax.lt.xlame) bigL = xlmax 
	bigel = (3.1416*qedge(n)*rmajor)

	bigM(n) = (1./bigel)/((1./xlmax)+(1./xlame))
c	bigM(25) = 0.
c	chiepale(n)=1.5*(1. + bigL/(3.1416*qedge(n)*rmajor))*const 
c	1			*xnue*(dele**2)
c	chiepale(n)=1.5*(1. + bigL/(3.1416*qedge(n)*rmajor))*const
c	1			*(1.4e3*zeff/xte(n)**1.5)*(couloge/17.)
c	correction 2/13/12
	const2 = ((1.414+zeff)/(1.414+(13.*zeff/4.)) + (1.+0.533/zeff))/
     1		 ((1 +sqrt(xnuestar)+1.65*(1.+0.533*zeff)*xnuestar)*
     2		 (1.+rmajor*qedge(n)/xlame))
      const2 = const2*(1-fc)/fc 
		
	chiepale(n)=1.5*(1. + bigM(n))*const
     1			*(1.4e3*zeff/xte(n)**1.5)*(couloge/17.)
 
c		electron neoclassical
	xnue = vthe/xlame
	rhote = 3e-6*sqrt(xte(n))/bfield 
	chieneo(n)= xnue*(rhote**2)*(qedge(n)**2)/(ep**1.5)
 
	xx = xnue*(dele**2)
	xy = (1.4e3*zeff/xte(n)**1.5)*(couloge/17.)
c	paleoclassical flux
	gnu(n)=(0.5*sqrt(xnuestar)+xnuestar)/(1.+sqrt(xnuestar)+xnuestar)
	gnu = 0.	
	 
	dep = 0.	
     	coefy(n) = (xk*dhat*exne(n)/sqrt(xte(n)))*
	1	(1.0-aminor*rhor(n)*(1./xlne(n))*(1.-xnuestar*dlndhatdnu) +
     2	0.5*aminor*rhor(n)*(1./exlte(n))*(1.-4.*xnuestar*dlndhatdnu))
	surfA(n) = (6.28**2)*rmajor*aminor*rhor(n)*sqrt(0.5*(1.+elong**2))
c	calculated surface power flow in electrons 
	epower(n)=(gamheate(n)-2.5*xk*zeff*gamion(n,1)*xte(n))*surfA(n)
c	uses calculated electron power flow for paleo power flow at separatrix 
	epowpcsep2= epower(25)

c	uses Callen's estimated paleo power flow at separatrix
	epowpcsep3 = -1.*coefy(25)


	if(n.lt.25) then
		sum = 0.0
		do 2660 np = n,24
		sum = sum + (1.+0.5*(bigM(np+1)+bigM(np)))*
     1		(coefy(np+1)-coefy(np))
2660		continue
	endif
	VpA = (6.28**2)*rmajor*aminor*rhor(n)*(1.+elong**2)/(2.*elong)
	epowpc2 = epowpcsep2 + sum 
	fluxpal(n) = epowpc2/VpA
      chiepale2(n) = fluxpal(n)/(exne(n)*xk*xte(n)/exlte(n)) 
	epowpc3 = epowpcsep3 + sum 
	fluxpal(n) = epowpc3/VpA
      chiepale3(n) = fluxpal(n)/(exne(n)*xk*xte(n)/exlte(n))

		if(n.eq.1) then
c	integrate paleo power flow out from rho = .86
	cpalfrac = 1.0	 
	epower86 = (gamheate(1)-2.5*xk*zeff*gamion(1,1)*xte(1))*surfA(1)
	epower86 = gamheate(1)*(1.-econve(1))*surfA(1)
 
	VpA = (6.28**2)*rmajor*aminor*rhor(n)*(1.+elong**2)/(2.*elong)
	fluxpal(1) = cpalfrac*epower86/vpA
	chiepale4(1) = fluxpal(1)/(exne(1)*xk*xte(1)/exlte(1))
	sum = 0.0
      do 2665 np = 1,24
	sum = sum + (1.+0.5*(bigM(np+1)+bigM(np)))*
     1		(coefy(np+1)-coefy(np))	
	VpA = (6.28**2)*rmajor*aminor*rhor(np+1)*(1.+elong**2)/(2.*elong)
	fluxpal(np+1) = (cpalfrac*epower86 - sum)/VpA
	chiepale4(np+1) = fluxpal(np+1)/
	1					(exne(np+1)*xk*xte(np+1)/exlte(np+1))
2665  continue
		endif

c	drift alfven (CPP,38,118,1998) 
	exlpe = 1./((1./xlne(n))+(1./exlte(n))) 
	chigb = (rhot**2)*csound/exlpe
	xmuu =(exlpe/(qedge(n)*rmajor))*
     1		sqrt(xmas(1)*xte(n)/(xme*xti(n)))
	beta = 2.*1.257e-6*exne(n)*xk*xte(n)/(bfield**2)
	betan =sqrt(xmas(1)/xme)*beta*qedge(n)*rmajor/exlpe
	xnun =((xmas(1)/xme)**0.25)*sqrt(qedge(n)*rmajor*exlpe)/xlame
	chiperp = (1./((1.+betan**2)**3)) + xnun**2
	chiperp = chiperp/(1.+betan**2+xnun**1.333)
	chiperp = sqrt(chiperp)
	chida0(n) = chigb*chiperp/sqrt(xmuu)
c	partial shear suppression 
	valfv = bfield/sqrt(1.257e-6*exne(n)*xmas(1)) 
	
	taualfv = (xlalfv**2)/chida0(n)
	ysda = omexb(n)*taualfv
	fsda(n) = 1./(1.+(ysda**2))
	chida(n) = fsda(n)*chida0(n)
c	stochastic magnetic field electrons

c	numerical calculation of dm by todd evans 7/26/06
	chiemag(n) = vthe*dmag(n)
	 
c	trapped electron mode drift wave

	chitem = (ep*ftr/xnue)*((csound*rhot)**2)/(xlne(n)*exlte(n))
	xxz = 0.1*cse/(xnue*qedge(n)*rmajor*(ep**1.5))
	fwesson = 2.5/(1. + xxz)
	chitemw0(n) = chitem*fwesson
	yyy= (csound*ep/(xlne(n)*xnue))**2.
	fkadom = 1./(1.+yyy)
	chitemk0(n) = chitem*fkadom
c	partial suppression	by ExB and magnetic shear
	chitemw(n) = fsitg(n)*chitemw0(n)/ss**1.8
	chitemk(n) = fsitg(n)*chitemk0(n)/ss**1.8

c	Bateman(PPCF,48,A93,2006) ETG
c	r = rho(n)*aminor	
c	xles = qedge(n)*rhote*rmajor/exlte(n)	
c	if(xles.ge.dele) goto 1805
c	chietg0 = 0.06*(dele**2)*vthe/sqrt(exlte(n)*rmajor)
	
c	goto 1810
c1805	chietg0 = 0.06*(qedge(n)**2)*((rmajor/exlte(n))**1.5)*
c	1  (clight*rhote/(eq*bfield))*xk*xte(n)*
c     2((1./exlte(n))-1.88*(ssi95/(qedge(n)*rmajor))*(1.+xte(n)/xti(n)))
c1810	continue
c	etg cutoff
c	cutoffetg = 0.8*rmajor/xlne(n)
c	xw = (1.+zeff*xte(n)/xti(n))*(1.33+1.91*r*shearho(n))*
c	1		(1.-1.15*ep)
c	if(xw.gt.cutoffetg) cutoffetg = xw
c	xa = rmajor/exlte(n)
c	if(xa.le.cutoffetg) then
c	chietg0 = 0.0
c	chieetg(n) = 0.0
c	endif 

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
	xnuatim(n)=	 YNO(n)*SVATA(n)

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
c	zbar2(n) = 6.0
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


c	chie0 = 0.1
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
	CHI0ION = CHICHos(N)
	CHI0EL  = CHIEPALE(N)
c 	chi0ion = 0.1
c	chi0el  = 0.1
C	CHI0ION & CHI0EL ARE THE BACKGROUND CHI'S IN M^2/S	 
	XKRMI2= 2.E-2
	XKRME2= 2.E-2

	CHITHI(N) = .667*(XKRMI2*(ALPHAION(N) - CHI0ION*XNU/(EXLTI(N)**2))
	1	- 2.5*XNU*VRAD1(N)/EXLTI(N) - CHI0ION)
	CHITHE(N) = .667*(XKRME2*(ALPHAEL(N) - CHI0EL*XNU/(EXLTE(N)**2))-
 	1	2.5*XNU*VRAD1(N)/EXLTE(N) - CHI0EL)

	

7100	continue 
C	**************END INFERRENCE OF CHI'S*******************************
7200 	continue
c***********************************
	goto 1594
c*************************************
c	integrate ion & electron temps & ion dens inward from separatrix
	
c	ti(25) = xti(25)
c	tel(25) = xte(25)
	yni(25,1)=EXNE(25)/(ATNUM(1)+FRACZ*ZBAR2(25))
	ypi(25) = yni(25,1)*xk*xti(25)
 	zne(25) = exne(25)
	znec(25) = exne(25)

	do 1593 j = 1, 24
	n = 25-j
	c6 = 1.0
c	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
c	xltem(N) = (((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
c    1					2.5*GAMEL/EXNE(N)))/xchie(n)	 
c	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))
c	xltim(N) = (((GAMHEATI(N)/(XNION*XK*XTI(N)))-
c     1					2.5*GAMION(N,1)/XNION))/xchii(n)
C 	xpon = -0.5*delna*(((1./exlti(n))+(vpinch4(n)-vrad1(n))/exd21(n))+
C	1		((1./exlti(n+1))+(vpinch4(n+1)-vrad1(n+1))/exd21(n+1)))
c	yni(n,1) = yni(n+1,1)*exp(xpon) 
	xpon = -0.5*delna*((vpinch4(n)-vrad1(n))/exd11(n) +
     1		(vpinch4(n+1)-vrad1(n+1))/exd11(n+1))
	ypi(n) = ypi(n+1)*exp(xpon)
	yni(n,1) = ypi(n)/(xk*xti(n))

c	yni(n,1) = yni(n+1,1)*(1.- delna*((1./exlti(n+1)) +
c     2		   (vpinch2(n+1)-vrad1(n+1))/exd11(n+1)))
     		
c	yni(n,1) = yni(n+1,1)*
c	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
c    2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
	yni(n,2) = fracz*yni(n,1)
	 
	znec(n) = yni(n,1) + zbar2(n)*yni(n,2)

c	ti(n) = ti(n+1)*
c     1	(1.+delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))/
c     2	(1.-delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))
c	if(n.eq.24) ti(n) = ti(n+1)/(1.-delma*xltim(n))	
c	x = xnudragyro1(n)
c	 y = xnudragyro2(n)
 
c      if(ti(n).lt.0.0) ti(n) = 100.
c	if(ti(n).gt.1.e3) ti(n) = 1.e3
c	tel(n) = tel(n+1)*
c     1	(1.+delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))/
c     2	(1.-delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))
1593	continue
	 
1594  continue 
c	goto 700
c	******************end heat transport investigation************	
	goto 3700
c	***************CALCULATE ROTATION*****************************

	do 3600 n = 1,25
	vphia(2) = torv(n)
  	vphia(1) = vphiex1(n)
	xni(1) =  yni(n,1)
	xni(2) =  yni(n,2)
	press(1) = -1.*(ti(n)/(atnum(1)*BTHET))*xlpm(n)
	press(2) = -1.*(ti(n)/(zbar2(n)*BTHET))*xlpm(n)

	xnuc(1,2) =	xnuc12(n)
	xnuc(2,1) = xnuc21(n)
c	ynud(1) =	xnud1(n)
c	ynud(2) =	xnud2(n)
	xmtor(1) =	xmomtor1(n)
	xmtor(2) =	xmomtor2(n)
c	if(ioptvthet.eq.0) then
c		velthet2(n) = vthexp(n)
c		velthet1(n) = velthet2(n) + delvelthet(n)
c	endif
C	
	call edgerotran(n,0)

C	**********INFER MOMENTUM TRANSPORT FREQUENCY******************* 
C	**********CALCULATE TOROIDAL & POLOIDAL ROTATION***************
c	goto 6799 
c	xnum = xni(1)*eq*atnum(1)*bthet*vrad1(n) + xmtor(1) + 
c	1		 xmtor(2) + (xni(1)*atnum(1)+xni(2)*zbar2(n))*eq*ephia
	bracket = (velthet1(n)-velthet2(n))/fp + 
     1			ti(n)*xlpm(n)*(1.-atnum(1)/zbar2(n))/(atnum(1)*bthet)
	xnum = (xmtor(1)/xni(1)+eq*ephia+eq*bthet*vrad1(n))/xmas(1) -
     1		xnuc12(n)*bracket	 
      
c	xdenom = (xni(1)*xmas(1)+xni(2)*xmas(2))*torv(n) + 
c    1		xni(1)*xmas(1)*bracket  
	xdenom = torv(n) + bracket
	xnud1(n) = xnum/xdenom 
	xnudtot1(n) = xnud1(n)
	vphical(1) =((xmtor(1)/xni(1)+eq*ephia+eq*bthet*vrad1(n))/xmas(1)
     1			+ xnuc21(n))/(xnuc12(n)+xnud1(n))	 
	xnum2 = (xmtor(2)/xni(2)+zbar2(n)*(eq*ephia+eq*bthet*vrad2(n)))/
     1	xmas(2) + xnuc21(n)*bracket
 	xnud2(n) = xnum2/torv(n)
	xnudtot2(n) = xnud2(n)
6799	continue 
	vphical(2) = xnum2/xnud2(n)
	vphical(1) = vphical(2)*(1. + xnud2(n)/xnuc21(n)) -
     1	 (xmtor(2)+xni(2)*zbar2(n)*eq*(ephia+bthet*vrad2(n)))/
     2	 (xni(1)*xmas(1)*xnuc12(n)) 	 
c	vphical(1) = ((xmtor(1)+xmtor(2))/(xni(1)*xmas(1)*xnuc(1,2))+
c	1  ynud(2)*((velthet1(n)-velthet2(n))/fp - (press(1) -press(2))))/
c	2	(ynud(1)+ynud(2))
c	 vphical(1) = vphical(1) + (atnum(1)*eq*(bthet*gamion(n,1) +
c     2		xni(1)*ephia)+zbar2(n)*eq*xni(2)*ephia)/
c     3		(xni(1)*xmas(1)*xnuc(1,2)*(ynud(1) + ynud(2))) 
c	vphical(2) = (1.+ynud(1))*vphical(1) -
c     1	(atnum(1)*eq*(bthet*gamion(n,1)+xni(1)*ephia) + xmtor(1))/
c     2	(xni(1)*xmas(1)*xnuc(1,2))	
 

	vtor1(n) = vphical(1)
	vtor2(n) = vphical(2)
	vphia(1) = vphical(1)
	vphia(2) = vphical(2)
3990	ioptapproach = 5
c	call edgerotran(n,0)
c	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
c	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
c	vtheta(1) = vtheta(1)*VTH(1)*fp
c	vtheta(2) = vtheta(2)*VTH(2)*fp
c	velthet1(n) = vtheta(1)
c	velthet2(n) = vtheta(2)
c	xnudragyro1(n) = xnudrag(1)
c	xnudragyro2(n) = xnudrag(2)

3600	continue

 
3700	continue
	OPEN(unit=121,FILE='pedestal.TXT',STATUS='UNKNOWN') 
	 
   	OPEN(122,FILE='pedestalplot.TXT',STATUS='UNKNOWN') 
	if(jwarn.eq.1) write(121,1999)
1999	format(1x,'not converged on 2000 loop after iterations=',I3.0) 
1030  format (i4)
	do 7440  n = nx,24
	write(121,1030) n 
	write(121,1000)  (epmin(n,m),m=n+1,25)
7440  continue
	write(121,'(1x,35A)') '   rho	     ne         ti         te     
     1  vrad    vpinch5  gamconve  '
	write(122,'(1x,35A)') '   rho	     ne         ti         te  ' 
c	sepdif = xnsepex - xnesepreal
     	do 750 n = 1,25
	gamconve = econve(n)*gamheate(n)
	zne(n) = zne(n) - sepdif
	write (121,1000) rhor(n),zne(n),ti(n),tel(n),vrad1(n),vpinch5(n),
	1   gamconve	
	write (122,1000) rhor(n),zne(n),ti(n),tel(n)
	sion(n) = zne(n)*0.5*(ynuioni(n) + ynuioni(n+1))
c	if(n.ne.25) xlnm(n+1) = log(yni(n,1)/yni(n+1,1))/delna
750	continue

	write(121,'(1x,35A)')'   rho    gamion  gamionorb gamionorbCUR gam
	1ionX   xgamionorb  xtransink  xtransourc'
	do 751 n = 1,25
	GAMIONORBCUR(N) = GAMION(N,1)*(1. - 2.*FORBL(N))
	write(121,1000) rhor(n),gamion(n,1), gamionorb(n),gamionorbCUR(n),
	1   gamionX(n), xgamionorb(n),xtransink(n),xtransource(n)
 751  continue
	Write(121,'(1x,35A)')'   rho     gamheati   heatiorb heatiorbX
	1heatiXtran heatiXorb ' 
     	do 752 n = 1, 25
	write(121,1000)rhor(n),gamheati(n),gamheatiorb(n),gamheatiorbX(n),
	1  gamheatiX(n),xheatXtranorb(n)
 752  continue
	write(121,'(1x,35A)')'   rho      chii    chiiorb   chiiorbX    ch
     1ii15  chiiorb15  chiiorbX15  '	 
	do 753 n = 1, 25
	write(121,1000) rhor(n),xchii(n), xchiiorb(n), xchiiorbX(n),
	1xchii15(n),xchiiorb15(n),xchiiorbX15(n)
753   continue	        
	write(121,'(1x,35A)')'   rho      xxchii   xxchiiorb  xxchii15
	1xxchiiorb15  chii  chii15'	 
	do 754 n = 1, 25
	write(121,1000) rhor(n),xxchii(n), xxchiiorb(n), xxchii15(n),
     1	xxchiiorb15(n),xchii(n),xchii15(n)
754   continue 
	write(121,'(1x,35A)')'   rho       qi25     qi15    qi25orb    qi1
	25orb  qi25orbx  qi15orbx       qi25xtran  qi15xtran  '
	do 7545 n = 1, 25
	write(121,1000) rhor(n),qcond25(n),qcond15(n),qcond25orb(n),
     2  qcond15orb(n),qcond25orbx(n),qcond15orbx(n)
     
     ,
7545	continue     	 
	write(121,'(1x,35A)')'   rho       qi25     qi15    qi25xtran  qi1
	25xtran q25xtranorb q15xtranorb  '
	do 7546 n = 1, 25
	write(121,1000) rhor(n),qcond25(n),qcond15(n),qcond25xtran(n),
     2	qcond15xtran(n),qcond25xtranorb(n),qcond15xtranorb(n)
7546	continue     	 


	write(121,'(1x,35A)') '   rho	   gamheate    gamheati   gamion     
     1heatinerti heatvisci gamconvi  '
     	do 755 n = 1,25
	gamconde(n) = (1.-econve(n))*gamheate(n)
	gamconvi = 2.5*xk*gamion(n,1)*xti(n)
	heatin(n)= gamion(n,1)*0.5*xmas(1)*(vtor1(n)**2+velthet1(n)**2)
	call param(n)
	vthi = sqrt(2.*xti(n)*xk/xmas(1))
	eta0(n) = yni(n,1)*xmas(1)*vthi*qedge(n)*rmajor*f(1)
	eta4(n) = yni(n,1)*xmas(1)*xti(n)*xk/(eq*bfield)
	heatvisc(n) = vtor1(n)*(fp*eta0(n)*vrad1(n)-
	1			  0.5*eta4(n)*(4.0*vtor1(n)+velthet1(n))) -
     1				0.5*velthet1(n)*(eta0(n)*vrad1(n)+
     2				eta4(n)*(vtor1(n)+0.5*velthet1(n)))
c	Vsin multiplier	
      vees = 0.0 
      heatvisc(n) = (vees/rmajor)*heatvisc(n)	 
	write (121,1000) rhor(n),gamheate(n),gamheati(n),
     1	gamion(n,1), heatin(n),heatvisc(n), gamconvi
755	continue
     	write(121,'(1x,35A)') '   rho		eta0        eta4	  qcondi 
	1qconde  chii1    chii2  chii3   '
	do 756 n = 1,25
	gamconvi = 2.5*xk*gamion(n,1)*xti(n)
	qcondi(n) = gamheati(n)-gamconvi-heatin(n)-heatvisc(n)
	Q1 = gamheati(n)
	Q2 = gamheati(n) - gamconvi
	Q3 = qcondi(n)
	chii1 = Q1*exlti(n)/(yni(n,1)*xk*ti(n))
	chii2 = Q2*exlti(n)/(yni(n,1)*xk*ti(n))
	chii3 = Q3*exlti(n)/(yni(n,1)*xk*ti(n))
	write(121,1000) rhor(n),eta0(n),eta4(n),qcondi(n),gamconde(n), 
     1	chii1,chii2,	chii3
756	continue 

	write(121,'(1x,35A)') '   rho      Qheati	Qheatiorb   Qconvi25 
	2Qconvi15 Qconvi25orb Qconvi15orb  '
	
	do 757 n= 1,25
	Qi = gamheati(n)
	Qiorb= gamheati(n)*(1. - eorbl(n))
	Qconvi25 = 2.5*xk*gamion(n,1)*xti(n)
 	Qconvi15 = 1.5*xk*gamion(n,1)*xti(n)
	Qconvi25orb = 2.5*xk*gamion(n,1)*xti(n)*(1.- forbl(n))
	Qconvi15orb = 1.5*xk*gamion(n,1)*xti(n)*(1.- forbl(n))
	write(121,1000) rhor(n),Qi,Qiorb,qconvi25,qconvi15,qconvi25orb,
     2	qconvi15orb
757	continue
	write(121,'(1x,35A)')'   rho     chii25    chii25orb chii15   
     2chii15orb '
	do 758 n = 1,25 
	gamionorb(n) = gamion(n,1)*(1-forbl(n))
	Qi = gamheati(n)
	Qiorb= gamheati(n)*(1. - eorbl(n))
	Qconvi25 = 2.5*xk*gamion(n,1)*xti(n)
 	Qconvi15 = 1.5*xk*gamion(n,1)*xti(n)
	Qconvi25orb = 2.5*xk*gamion(n,1)*xti(n)*(1.- forbl(n))
	Qconvi15orb = 1.5*xk*gamion(n,1)*xti(n)*(1.- forbl(n))
	xx = exlti(n)/(yni(n,1)*xk*xti(n))
	chii25 = xx*(Qi-Qconvi25)
	chii25orb = xx*(Qiorb-Qconvi25orb)
	chii15 = xx*(Qi-Qconvi15)
 	chii15orb = xx*(Qiorb-Qconvi15orb)
	write(121,1000) rhor(n),chii25,chii25orb,chii15,chii15orb
758	continue
	write(121,'(1x,35A)')'  rho       chii25    chiiorb25   chii15  ch
	1iiorb15 gamionorb cdxtran '
	do 759 n = 1, 25
	write(121,1000) rhor(n),xchii(n),xchiiorb(n),xchii15(n),
     1	xchiiorb15(n),gamionorb(n),cdxtran(n)
759	continue 
	write(121,'(1x,35A)')'   rho     Eionorb    Extran     Emom    
	2Evpol   Epress     Etot    Eexp   Etot2   '
	do 7595 n=1,25
	xniexp = exne(n)/(atnum(1)+fracz*zbar2(n))
	xnzexp = xniexp*fracz
	rhomassd = xniexp*xmas(1)*xnudtot1(n) + xnzexp*xmas(2)*xnudtot2(n)
	Eionorb = -1.*cdionorb(n)*((fp*bphi)**2)/rhomassd
C     cdionorb is outward ion orb current dens, cdxtran is compensating inward xtran cur den
	Extran	= +1.*cdxtran(n)*((fp*bphi)**2)/rhomassd
	Emom    = (xmomtor1(n)+xmomtor2(n))*abs(fp*bphi)/rhomassd
	velthet3(n) = fp*vtor1(n)-erex(n)/bphi - ti(n)*xlpm(n)/bphi
	Evpol   = -1.*(xniexp*xmas(1)*xnudtot1(n)*velthet3(n)*bphi +
     2	xnzexp*xmas(2)*xnudtot2(n)*velthet2(n)*bphi)
	Evpol = Evpol/rhomassd 
	Epress  = -1.*xlpm(n)*ti(n)*(xniexp*xmas(1)*xnudtot1(n)+
     2xnzexp*xmas(2)*xnudtot2(n)/zbar2(n))/rhomassd
c 	Extran = 0.0
c	Eionorb = 0.0
	Etot = Eionorb+Extran+Emom+Evpol+Epress
	Etot2 = Emom+Evpol+Epress
      write(121,1000) rhor(n),Eionorb,Extran,Emom,Evpol,Epress,Etot,
     2	erex(n), etot2	 
7595	continue     
	write(121,'(1x,35A)')'   rho    -1.*cdionorb  cdxtran '
	do 7611 n = 1,25
	write(121,1000) rhor(n),-1.*cdionorb(n),cdxtran(n)
7611	continue
 	write(121,'(1x,35A)')'   rho    velthet1  velthet3 '  
	do 7612 n = 1,25
	write(121, 1000)  rhor(n),velthet1(n),velthet3(n) 
7612	continue
 
	write(121,'(1x,35A)') '   rho	    nudrag1*    nuatom1   momphi1       
     1nugyro1  nurip1  xdrag1sep  '
     	do 760 n = 1,25
	write (121,1000) rhor(n),xnudtot1(n),xnuati(n)+xnuioni(n)+xnuionb(
     1n),xmomtor1(n),xnudragyro1(n),xnudragrip(n,1),XNUDtot1s(n)
	dragfreq(n) = xnudtot1(n) 
760	continue
	write(121,'(1x,35A)') '   rho	    nudrag2*    nuatom2   momphi2       
     1nugyro2  nurip2  xdrag2sep  '
	xcb = 0.0	 
     	do 761 n = 1,25
	write (121,1000) rhor(n),xnudtot2(n), xcb,
     1xmomtor2(n),xnudragyro2(n),xnudragrip(n,2),XNUDtot2s(n)
761	continue

	write(121,'(1x,35A)')'  rho      vphi1     vphi2     vphi2exp 
	    
     2vphi1ex    DELV0   DELV1    vphi1sep   vphi2sep '
c*********Eq 4 from Structure... paper April, 2013*******************
	do 7615 n = 1,25
	vphi1p=(torq1(n)-yni(n,1)*xmas(1)*xnuc12(n)*delv1(n))/
     2		(yni(n,1)*xmas(1)*xnudtot1s(n))
      vphi2p=(torq2(n)+yni(n,2)*xmas(2)*xnuc21(n)*delv1(n))/
     2		(yni(n,2)*xmas(2)*xnudtot2s(n))
c   this is vphi1p & vphi2p of p86 notes calculated w/sep mom bal and sep nud1 & nud2
	torq1(n) = torq1(n)/(yni(n,1)*xmas(1)*(xnuc12(n)+xnudtot1(n)))
	torq2(n) = torq2(n)/(yni(n,2)*xmas(2)*(xnuc21(n)+xnudtot2(n)))
	xx = (xnuc12(n)+xnudtot1(n))*(xnuc21(n)+xnudtot2(n))
	brackx = 1. -  xnuc12(n)*xnuc21(n)/xx
	vph1=(torq1(n)+xnuc12(n)*torq2(n)/(xnuc12(n)+xnudtot1(n)))/brackx
	vph2=(torq2(n)+xnuc21(n)*torq1(n)/(xnuc21(n)+xnudtot2(n)))/brackx
	vphi1ex = torv(n) + delv1(n)
c  this is vphi1 & vphi2 on p70 of notes calc w/ nud1=nud2=nud0
	write(121, 1000) rhor(n),vph1,vph2,torv(n),vphi1ex,DELV0(N),
     2	DELV1(N),vphi1p, vphi2p
7615	continue
 	 	write(121,'(1x,35A)')'  rho      eradcur    eradrot   eradpres    
     2eradtot  eradex  '
c*********Eq 4 from Structure... paper June, 2013*******************
	do 76155 n = 1,25
	torq1(n) = torq1(n)*((xnuc12(n)+xnudtot1(n))/
	2	(xnuc12(n)+xnudzero(n)))
	torq2(n) = torq2(n)*((xnuc21(n)+xnudtot2(n))/
	2	(xnuc21(n)+xnudzero(n)))
	xx = (xnuc12(n)+xnudzero(n))*(xnuc21(n)+xnudzero(n))
	brackx = 1. -  xnuc12(n)*xnuc21(n)/xx
	vph1=(torq1(n)+xnuc12(n)*torq2(n)/(xnuc12(n)+xnudzero(n)))/brackx
	vph2=(torq2(n)+xnuc21(n)*torq1(n)/(xnuc21(n)+xnudzero(n)))/brackx
	vphi1ex = torv(n) + delv1(n)
	
c***************calculate Erad from Eq 7 of structure paper--June, 2013******
	xjk = yni(n,1)*xmas(1)/(yni(n,2)*xmas(2))
	rotj = (vpol61(n)*bphi-vph1*bthet)/(1. + (1./xjk))
	rotk = (vpol62(n)*bphi-vph2*bthet)/(1. + xjk)
	eradrot = -1.*(rotj+rotk)
c****resistivity****************
	xloglambda = 18.0
	xme = 9.1e-31
	vnum = sqrt(xme)*sqrt(xk)*xloglambda
	Vdenom = 32.*sqrt(3.14)*((8.85e-12)**2)*((2.*xte(n))**1.5)
	zeff = (yni(n,1)+(zbar2(n)**2)*yni(n,2))/
	2          (yni(n,1)+zbar2(n)*yni(n,2))
	vnum = (1.03e-4)*zeff*xloglambda
	vdenom = xte(n)**1.5
	etares  = vnum/vdenom
c*******************************
	eradcur = -1.*eq*etares*(gamion(n,1)*forbl(n))
	pressjk = xk*(yni(n,1)+yni(n,2))*ti(n)
	eradpres= -1.*pressjk*xlpm(n)/(eq*(yni(n,1)+zbar2(n)*yni(n,2)))
	eradtot = eradcur + eradrot + eradpres

	write(121,1000) rhor(n),eradcur,eradrot,eradpres,eradtot,erex(n)
76155	continue


	write(121,'(1x,35A)') '   rho	     xno      coldno    cxcool       
     1radcool  ioncool   qnbe  '
 
	do 765 n = 1,25
	write (121,1000) rhor(n),yno(n),coldno(n),cxcool(n),
     1radcool(n),coolion(n), qnbe(n)
765	continue 
	write(121,'(1x,35A)') '   rho	 -(dp/dr)/p -(dTi/dr)/T -(dn/dr)/n    
     1-(dTe/dr)/T xnuion qnbi '
   
	do 770 n=1,25
	xlnm(n) = 1./xlne(n)
	write (121,1000) rhor(n),xlpm(n),xltim(n),xlnm(n),
     1						xltem(n),xnuioni(n)+xnuionb(n),qnbi(n)
     
          
770	continue 
		write(121,'(1x,35A)') '   rho	      exlti     exlte    exlne
     1 sion       zbar  rhon  '
        
	do 7701 n=1,25
	xlnm(n) = 1./xlne(n)
	write (121,1000) rhor(n),exlti(n),exlte(n),xlne(n),ssion(n),
     1   zbar2(n),rhorn(n)          
7701	continue 
C*****************gradient scale lengths*********************************\
	write(121,'(1x,35A)') '  rho       Ln  	    snion(n)   forb       
     2gamion  1/r term   S term   nujk term  '
	do 7702 n = 1, 24
	snion(n) = yni(n,1)*(xnuioni(n)*(1.+fracz*zbar2(n))+xnuionb(n))
	rr = rhor(n)*aminor*sqrt(0.5*(1.+ elong**2))
	denumb = yni(n,1)*velthet1(n)
	denom = (1.-2.*forbl(n))*gamion(n,1)
      denom = 1.
	xLnmin1 = (1./rr)*(1.+denumb/denom)
	xLnmin2 = -1.*snion(n)/denom
	xLnmin3 = -1.*yni(n,1)*xnuc12(n)/denom
	xlnmin(n) = xLnmin1+xLnmin2+xLnmin3
	xlength =  1./xlnmin(n)
	write (121,1000) rhor(n),xlength, snion(n),forbl(n),gamion(n,1),
	2                 xLnmin1,xLnmin2,xLnmin3 
     				  	   
7702  continue
	write(121,'(1x,35A)') '  rho       Ln  	    exlne     exlti   for       
     2b      gamion    xlpjmin    xlnjmin  '
	do 7703 n=1,25 
	alpha = fracz*(zbar2(n)**2)
	const1 = 0.47 + 0.35/(0.66+alpha)
	const2 = 0.30 + 0.41/(0.58+alpha)
	rr = rhor(n)*aminor*sqrt(0.5*(1.+ elong**2))
	ep = rr/rmajor
	const3 = ((eq*bthet)**2)/
	2		 (2.*(ep**2)*yni(n,1)*xmas(1)*xnuc12(n)*xk*ti(n))
	xlpjmin = (const3*gamion(n,1)*(1.-forbl(n)) + 
	2			(const2/exlti(n)))/const1
	xlnmin(n) = xlpjmin-(1./exlti(n))
	xlength = 1/xlnmin(n) 
	write (121,1000) rhor(n),xlength, xlne(n),  exlti(n),forbl(n),
	2         gamion(n,1), xLpjmin,xLnmin(n)

7703	continue
c**************************************************************************
	write(121,'(1x,35A)') '  rho     xchii2     xchie2    diffii    di 
	1ffiz    diffzz   diffzi ' 
	do 771 n = 1,25
	
	write(121,1000) rhor(n),xchii2(n),xchie2(n),diffii(n),
     1	diffiz(n),diffzz(n),diffzi(n)                                       
771	continue 
	WRITE(121,'(1X,35A)') '	 rho       q         shear     bigM      
     1Y     xnuei*     qie '
	do 7710 n = 1,25
	write(121,1000) rhor(n),qedge(n),shearho(n),bigM(N),coefy(N),
     1	starnue(N),qie(n)

7710	continue
  
	write(121,'(1x,35A)') ' chiiexp56    neo1   w/6      da15
     1w/ExB   thermin  '		  
	write(121,'(1x,35A)') ' cond q      neoch   w/orbsq  driftAlfen
     1da15&28   ions31  '					 
	do 772 n = 1,25
	
	write(121,1000) xchii(n),chich(n),chichos(n),chida0(n),chida(n),
     1	chithi(n)

772	continue 

	write(121,'(1x,35A)')' chiiexp56  itg8      itg8&23    itg8&24   i
	1tg20    itg20&23  itg20&24  elect31  '
	write(121,'(1x,35A)')' orb15  q     Romanelli w/ExB      wExB&dqdr K
	1alupin  w/ExB     wExB&dqdr therminstab  '
 
	do 7721 n = 1, 25
	write(121,1000) xchii3(n),chietair0(n),chietair1(n),chietair2(n),
     1             chietaiw0(n),chietaiw1(n),chietaiw2(n),chithe(n)
7721	continue      	 
	write(121,'(1x,35A)')' chiiexpe60 itg9a    it9a&23   itg9a&24  itg  
	19b    itg9b&23   itg9b&24 '
	write(121,'(1x,35A)')' chii orb25 horton#1  w/ExB    wExB&dqdr ho  
	1rt#2   wExB     wExB&dqdr '
 
	do 7722 n = 1, 25
	write(121,1000) xchii2(n),chietaih0(n),chietaih1(n),chietaih2(n),
	1			chietai2h0(n),chietai2h1(n),chietai2h2(n)
7722	continue   
	write(121,'(1x,35A)') ' chieexp56 chiepal32 chiepal35  temwesson w    
     2ExB&dqdr tem43     tem43&47      '   
	write(121,'(1x,35A)') ' conduct q formula   powerbal   wesson  w/E   
     2xB&dqdr kadomtsev  wExB&dqdr    '   
	do 773, n= 1,25
	write(121,1000) xchie(n),chiepale(n),chiepale4(n),chitemw0(n),
     1				chitemw(n),chitemk0(n),chitemk(n)
773	continue
	write(121,'(1x,35A)')'  chieorb   drb48   drb48&50  drb48&51  et 
	1g39    etg39&41    etg42   etg42&41  '
	write(121,'(1x,35A)')'  total Q     resbal  w/ExB     wExB&dqdr we  
	1sson w/Hamaguchi  Horton   w/Hamaguchi  '
 
	do 7731, n = 1,25
	ss =aminor*shearho(n)*rminor/qedge(n)
	chierb2= chierb(n)/ss**1.8
	write(121,1000) xchieorb(n),chierb0(n),chierb(n),chierb2, 	
     1	chieetg0(n),chieetg1(n),chieetg2(n),chieetg(n)
7731	continue  
   
	write(121,'(1x,35A)') ' omexb     fsitg      fsetg      fsdrb    
     1fsda   etgRLcrit   itgRLcrit ysitg' 
	write(121,'(1x,35A)') ' Eq21       Eq22       Eq40       Eq49    
     1Eq28     Eq38		Eq7      Eq22' 
 
	do 774, n=1,25
	write(121,1000) omexb(n),fsitg(n),fsetg(n),fsdrb(n),fsda(n),
	1	critetg(n),crititg(n),ysitg(n)
774	continue
	write(121,'(1x,35A)')' etai/thresh etai   etaithresh  gamitg      
	1gamtem   omreal   itg12b    tem44	'

	do 7745, n=1,25
	write(121,1000) etai(n)/etaithresh(n),etai(n),
     1  etaithresh(n),gamitg(n),gamtem(n),omreali(n),chietaiw9(n),
     2chitemw9(n)
7745	continue
	write(121,'(1x,35A)')'  itg12b  itg12b&23  itg12b&24  tem44
     1tem44&46  tem44&47  itg12a '
      do 7746 n=1,25
	chiiweil0 = chietaiw9(n)
	chieweil0 = chitemw9(n)
	chiiweil1 = fsitg(n)*chiiweil0
	chieweil1 = fsitg(n)*chieweil0
	ss =aminor*shearho(n)*rhor(n)/qedge(n)
	chiiweil2 = chiiweil1/(ss**1.8)
	chieweil2 = chieweil1/(ss**1.8)
	write(121,1000) chiiweil0,chiiweil1,chiiweil2,chieweil0,chieweil1,
     1				chieweil2,chiitgx(n)
7746  continue   		 		  	 
c	write(121,'(1x,35A)') ' alphaion  chii(k+L)  omegi	 alphael  chie
c	1(k+L)  eflux  fluxpal'
c	do 774, n=1,25
c	chii0 = chichos(n)
c	chie0 = chiepale(n)

c	write(121,1000) alphaion(n),-1.*chii0*((xkredge**2)+1./
c	1				(exlti(n)**2)),omegthi(n),alphael(n),
c    2 -1.*chie0*((xkredge**2)+1./(exlte(n)**2)),epower(n)/surfA(n),
c	3  fluxpal(n)
c774	continue       
	write(121,'(1x,35A)') '   rho     chiedw0    chiedw1    chiedw2  
	1 eta-i    eta-e  '   
 
	do 775 n = 1,25
	write (121,1000) rhor(n),chiedw0(n),chiedw1(n),chiedw2(n),
     1  etai(n),etae(n) 	
     
775	continue 
	write(121,'(1x,35A)') '   rho     thetion     thetimp     G   VTHE 
	1TA1calc VTH1cor  vthet2 '   
 
	do 780 n = 1,25
	
	write (121,1000) rhor(n),thetwid1(n),thetwid2(n),gy(n),
     1	velthet1(n), vth1cor(n), velthet2(n)	 
780	continue

	write(121,'(1x,35A)') '   rho      beam      ephi      drag  
	1 fric   vpinch2  '
	do 787 n = 1,25
	zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
	 zdrag =	(xmas(1)*ynudrag1(n)*((erada(n)/bthet)+velthet2(n)/fp))/
	1	 (eq*bthet)
	zfrict = 0.0
	vp2 = zbeam - ephia/bthet+zdrag+zfrict 
	xnuc12(n) = xnuc120(n)
	write (121,1000) rhor(n),zbeam,-1.*ephia/bthet,zdrag, zfrict,vp2
787	continue 

	write(121,'(1x,35A)') '   rho      beam      ephi      drag  
	1 fric   vpinch3  '   

	
 
	do 785 n = 1,25
		zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
	zdrag =	(xmas(1)*ynudrag1(n)*((erada(n)/bthet)+velthet2(n)/fp))/
     1	 (eq*bthet)
	zfrict = (xmas(1)*(xnuc12(n)+ynudrag1(n))*
	1		(velthet1(n)-velthet2(n))/fp)/(eq*bthet)
	vp3 = zbeam - ephia/bthet+zdrag+zfrict 
	 
     	write (121,1000) rhor(n),zbeam,-1.*ephia/bthet,zdrag, zfrict,vp3
 
    
     
785	continue 
		write(121,'(1x,35A)') '   rho    beam&ephi  erad      vpol   
	1vphi    vpinch4  ' 
	do 7860 n = 1,25
		zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	 eradterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)))/
	1			(eq*bthet)
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((velthet1(n)/fp)))/
     1			(eq*bthet)
	zfrict = 0.
	zdrag = eradterm+vtheterm
 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp4 = zbeam-ephia/bthet+vphiterm+eradterm+vtheterm+zfrict
 

     	write (121,1000) rhor(n),zbeam-1.*ephia/bthet,eradterm, vtheterm,
     1			vphiterm,vp4
7860	continue 


	write(121,'(1x,35A)') '   rho    beam&ephi    erad      vpol   
	1vphi    vpinch5  ' 
	do 786 n = 1,25
		zbeam =  -1.0*xmomtor1(n)/(yni(n,1)*eq*bthet) 
 	 eradterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)))/
	1			(eq*bthet)
	vtheterm =	(xmas(1)*(ynudrag1(n)+xnuc12(n))*((velthet2(n)/fp)))/
     1			(eq*bthet)
	zfrict = (xmas(1)*(xnuc12(n)+ynudrag1(n))*
 	1		(velthet1(n)-velthet2(n))/fp)/(eq*bthet)

	polterm = zfrict + vtheterm 
 	vphiterm = -1.*xmas(1)*xnuc12(n)*torv(n)/(eq*bthet)
	vp5 = zbeam-ephia/bthet+vphiterm+eradterm+vtheterm+zfrict

     	write (121,1000) rhor(n),zbeam-1.*ephia/bthet,eradterm,polterm,
     1			vphiterm,vp5
786	continue 

		write(121,'(1x,35A)') '   RHO      EXD11     EXD22     EXD31    
     1EXD41    EXD51    EXD53   '
     
	DO 7711 N = 1,25
	WRITE(121,1000) RHOR(N),EXD11(N),EXD22(N),EXD31(N),EXD41(N),
	1	EXD51(N),EXD53(N)
7711  CONTINUE	
	write(121,'(1x,35A)') '   rho     vpinch2    vpinch3   vpinch4
	1 vpinch5  vtherm'
	do 7712 n= 1, 25
 	 VTHerm = SQRT(2.*XK*TI(N)/XMAS(1))
	write(121,1000)rhor(n),vpinch2(n),vpinch3(n),vpinch4(n),vpinch5(n)
	1,vtherm
7712	continue	
	
	write(121,'(1x,35A)') '   rho      Bvphi     -Bvthet  press  
	1  eradfb2  eradex  '   
 
	do 790 n = 1,25
c	press(2) = -1.*ti(n)*xlpm(n)/zbar2(n)	
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
	1eradnew  eradexp   vthet3'
      do 795 n = 1,25
	write (121,1000) rhor(n),velthet1(n), velthet2(n),vthexp(n),
     1eradfbnew(n),erex(n),velthet3(n)    
795	continue
	write(121,'(1x,35A)') '   rho     vpol61     vpol62      vthexp  
	1   vthet3'
      do 7951 n = 1,25
	write (121,1000) rhor(n),vpol61(n), vpol62(n),vthexp(n),
     1velthet3(n)    
7951	continue

	write(121,'(1x,35A)') '   rho     vtheory     Wminx2     Wminx      
	1eradav   erexp    ephi  '
	do 7966 n = 1,25
	write (121,1000) rhor(n),vtheory(n), Wminx2(n), Wminx(n),erav(n),
     2	erex(n),ephi(n)
7966	continue 

	write (121,'(1x,35a)') '   rho     xlossn    xlosse    partloss   
     2 energyloss  	'
		xlossion = 0.0
	 xlosspow = 0.0 
	do 799 n = 1,25
	write (121,1000) rhor(n),xlossn(n),xlosse(n),partloss,energyloss
799	continue
 	write (121,1023) xlossion,xlosspow
1023	format (5x,"ion X-loss =",1x, e10.3,5x, "power X-loss =",1x,e10.3)	
      
	write (121,'(1x,35a)') '  omegi0    omegic    omegis    omefz0  
     1omegzc    omegzs  '
	do 805 n = 1,25 
      write (121,1000) omegt(n,1),omegt(n,2),omegt(n,3),omegt(n,4),
     1       omegt(n,5), omegt(n,6)		      
805	continue 
	WRITE (121,'(1X,35A)')'  rho    vphiSIN      VphiCOS     nsin     
	1  nCOS  '
      DO 806, N=1,25
	WRITE (121,1000) rhor(n), vphisin(n,1),vphicos(n,1),sinion(n),
	1				  cosion(n)
806	CONTINUE 

	WRITE (121,'(1X,35A)')'  VTH1      VTHTSIN  VTHTCOS  VTHTCOSZ   VT
     1HTSINZ    VTH2     vth2exp'
      DO 810, N=1,25
	WRITE (121,1000) Vpol(1,N),VTHTSIN(1,N),VTHTCOS(1,N),
     1   VTHTCOS(2,N),VTHTSIN(2,N),Vpol(2,N),vthexp(n)
810	CONTINUE 
	WRITE(121,'(1X,35A)')'  rho      vphi1      vphi2   vphi1intr   vp
     1hi2intr, vphi1cor  vphi2cor  '
	do 812 n= 1,24
	x = vtor1(n)-yy1(n)
	y = vtor2(n)-yy2(n)
	write (121,1000) rhor(n),vtor1(n),vtor2(n),yy1(n),YY2(n),x,y
812	continue
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
     1 BV15    bv16    xnuiI  '
	
      DO 825, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) B0V(1,1,N),B0V(1,2,N),B0V(1,3,N),B0V(1,4,N),
     1	B0V(1,5,n),b0v(1,6,n),xnuc12(n)
825	CONTINUE 
	  
	WRITE(121,'(1X,35A)')'  rho       potential  vt1HS      vt2HS 	      
     1vt1SS    vt2SS    vtexp     vt1HS* '
 	do 8260 n = 1,25 
	write(121,1000) rhor(n),epota(n),vths1(n),vths2(n),velthet1(n),
	1	velthet2(n),vthexp(n),vpolhs(n)
8260	continue
	
c     	goto 1559
	do 827 jj = 1,1
      do 826 nn = 1,25
	n = 26 - nn
c	if(ioptpinchi.eq.5) then
c	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1)  + 
c    1 xmas(1)*(ynudrag1(n)+xnuc12(n))*((erfb2(n)/bthet)+velthet1(n)/fp) 
c     2	 - xmas(1)*xnuc12(n)*torv(n))/
c     3	 (eq*atnum(1)*bthet)
	vpinchi(n) = vpinch4(n) 
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
	1			((ynudrag1(n)/xnuc12(n))+1.-atnum(1)/zbar2(n))/
     2			((eq*atnum(1)*bthet)**2)
	xlpmom(n) = (vrad1(n)-vpinchi(n))/diffA(n) 
c	endif
c	if(xlpmom(n).lt.0.0) xlpmom(n) = 0.0
	c10 = 1.0
	 xlnmold = xlnm(n) 
	xlnm(n) = xlpmom(n)-c10*xltim(n) 
c	if(rhor(n).gt.pedrhoti) xlnm(n) = xlpmom(n)-c10/ylti 
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpmom(n)-c10/xltim(n) 
c	if(rhor(n).le.pedrhoti) xlnm(n) = xlpmom(n)-c10/xltitop  
c	if(kk.gt.3) xlnm(n) = 0.5*(xlnm(n)+xlnmold)
c	if(xlnm(n).lt.0.0) xlnm(n) = 0.0
c	integrate ion density inward from separatrix
	c6 = 0.5
	if(n.eq.25) goto 826		
	yni(n,1) = yni(n+1,1)*
	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
     2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
c	yni(n,1)=yni(n+1,1)*exp(delna*(c6*xlnm(n)+(1.-c6)*xlnm(n+1))) 
c	if(n.eq.24) yni(n,1)=yni(n+1,1)/(1.-delma*xlnm(n))
	if(yni(n,1).le.0.0) yni(n,1) = 1.e19 
c	if(yni(n,1).gt.2.e20) yni(n,1) = 2.e20
	yni(n,2) = fracz*yni(n,1)
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	znec(n) = yne
826	continue
	write(121,'(1x,35A)')'  rho       nicalc     necalc    neexp   vpi
	1nch     vrad         D       xnu12     '
	do 8261 n = 1,25
	write(121,1000) rhor(n), yni(n,1), znec(n), exne(n),vpinchi(n),
     1						 vrad1(n), exd11(n), xnuc12(n)	 
8261	continue	
	call neutdist
827	continue
	 WRITE(121,'(1X,35A)')'  rho       vpinch     vrad1    xlpmom    
	2 xlpmex   Diff  '
	 do 8255 n = 1, 25
	 write(121,1000) rhor(n), vpinch4(n), vrad1(n),xlpmom(n),xlpm(n),
	2 diffA(n)
8255  Continue	 

1559	continue
c	*******************temporary*********************************
	goto 4099	
			 
	DO 828, N=1,25

	toromega = torv(n)/rmajor 
	WRITE(121,1000) rhor(N),vpinchi(N),xlpm(N),zne(n),yno(n),toromega
828	CONTINUE 
c	nudrags from ndrag=3 pert analysis, first nu0, assuming delvphi=0, 
c	then estimate delvphi_0 and recalculate nudrag 1&2.  
c	nugyros from calculate vpol, vph_main = vphi_impex + delphi_0
	write (121,'(1x,35A)')' nudrag_main nudrag_imp nudrag0  vphi_imp  
	1vphi_main  nugyro_main nugyro_imp  '
	do 829 n = 1, 25
	write(121,1000) xnudtot1(n),xnudtot2(n),xnudzero(n),vtor2(n),
     1		vtor1(n),xnudragyro1(n),xnudragyro2(n)
829	continue 
c	calculate vthetas using	vphis above, then recal using also vthet2exp
	write (121,'(1x,35A)')'vthexp_impvthcal_impvthcal_main vpcal_imp
     1vpcal_main vpolHS_main vpolHS_imp'
	do 830 n = 1, 25
	write (121,1000) vthexp(n),velthet2(n),velthet1(n),vpcal2(n),
	1	vpcal1(n),vths1(n),vths2(n)
830	continue
c	vphis & nudrags from 6/25/07 quadratic eq sol for vphi_impex	 
 	 write (121,'(1x,35A)') '  vtkdg1  vtkdg2 vphi_imp  vphi_main1
     1vphi_main2  xnuatom  epota'
      do 831 n = 1,25
	write (121,1000) vtkdg1(n),vtkdg2(n),vftor2(n),vftor1(n),
     1	vftor0(n), xnuati(n)+xnuioni(n)+xnuionb(n),epota(n)	
831	continue
c	6/26/07 insert above solving pol mom eq using f_imp collisional and V_impex
	write (121,'(1x,35A)')' fpol_main  vpol_main  vpol_imp fpol_col  f
     1pol_shaing fimp_shaing '
	do 832 n = 1,25
	write (121,1000) fpol_main(n),vpol_main(n),vpol_imp(n),fpol_col(n)
	1,	fpol_shaing(n),fpol_imp(n)
832	continue
c	option #2 solution quadratic eqs for fpol_main and vpol_main
	write(121,'(1x,35A)')'   fpol+    fpol-   fpol_col    fpol_shaing  
     1vpol+    vpol-     vpolimp  '
      do 833 n = 1, 25
	write (121,1000) fpol1(n),fpol2(n),fpol_col(n),fpol_shaing(n),
     1		vpol3(n),vpol4(n),vpol_imp(n)
833	continue  
c	option #3 solve for vj and vk using fshaing for each j and k
	write (121,'(1x,35A)')'  vpol_ion   vpol_imp  fshaing_imp  press1
	1   press2  '
	do 834 n = 1, 25
	write (121,1000) vpol5(n),vpol6(n),fpolimp_shaing(n),press1(n),
	1	press2(n)
834	continue   				 
	WRITE(121,'(1X,35A)')'  cospot   sinpot    potential   ephi/T   
	1av15 	av16  '      
	
      DO 835, N=1,25
	toromega = torv(n)/rmajor 
	WRITE(121,1000) cpot(n),spot(n),estatpot(n),estatpot(n)/xte(n),
     1	av(1,5,n),av(1,6,n)
835	CONTINUE 
	WRITE(121,'(1X,35A)')'  sinion   cosion   vthtsin     vthtcos   
     1vphisin   vphicos  '	        
	
      DO 840, N=1,25
 	WRITE(121,1000) sinion(n),cosion(n),vthtsin(1,n),vthtcos(1,n),
	1	vphisin(1,n),vphicos(1,n)
840	CONTINUE 

	WRITE(121,'(1X,35A)')'  rho      torvel1    torvel2   vphi2ex   vp  
     1hi1ex   veltor1    veltor2  vtherm1  '

	DO 845, N=1,25
		do 841 j=1,2
 	VTH(J) = SQRT(2.*XK*TI(N)/XMAS(J))
841	continue 
       

	WRITE(121,1000)rhor(n),torvel(N,1),torvel(N,2),torv(N),vphiex1(n), 
     1	veltor(n,1),veltor(n,2), vth(1)
     	
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
1000	format(9e10.3)
1005	format(1x,'chixpi=',f5.2,1x,'chixpe=',f5.2,1x,'chitop=',f5.2,1x,
     1	'chetop=',f5.2) 
1006	format(1x,'ylti=',f6.3,1x,'ylte=',f6.3,1x,'xltitop=',f6.3,
	1	1x,'xltetop=',f6.3) 

C**********************************************************************
c**********************************************************************					
C	CALCULATION USING INFERRED CHI's & XNUD's



		 
c 	ATNUM(1) = ZION
c	ATNUM(2) = ZIMP

c	ynpedz = xnped
c	ynbarz = xnbar
c	if(ioptsoln.eq.1) then
c	ynbarz = 0.5*(xnpedex+xnsepex)
c	ynpedz = xnpedex
c	endif
c	YNI(25,1) = xnsol
c	yni(25,2) = fracz*xnsol
c	if(joptedped.eq.1) then
c		yni(25,1) = xnsep
c		yni(25,2) = fracz*xnsep
c	endif
c	
c	TEED(25) = tsepexe
c	TIED(25) = tsepexi
c	ti(25) = tsepexi
c	tel(25) = tsepexe
c	GAMION(25,1) = enh*FLUXPART
c	gamheat(25) = fluxheat
	 
c	gamheate(25) = fheate*gamheat(25)
c	gamheati(25) =(1.-fheate)*gamheat(25)
c	gamion(25,2) = 0.0
c	xni(1) = xnsol
c	xni(2) = fracz*xnsol
c	tep = tsep
c	temp(1) = tsep
c	temp(2) = tsep
c	XLNA(1)	= YLNBARX
c	XLTA(1) = YLTIBARX
c	XLNA(2)	= YLNBARX
c	XLTA(2) = YLTIBARX
c	XLVA(1) = YLVBARX
c	XLVA(2) = YLVBARX
	
c	CNEUT = 1.
c	IF(IOPTELN.EQ.0) CNEUT = 0.
c	TED = teed(25)
c	TID = tied(25)		
c	IF(TED.LT.1.E-1) TED = 1.05E-1 
c	IF(TED.GT.1E3) TED = .95E3
c	IF(TID.LT.1.E-1) TID = 1.05E-1 
c	IF(TID.GT.1E3) TID = .95E3
c     	YND = xnsol
c	IF(YND.GT.1E22) YND = 0.95E22
c	IF(YND.LT.1E16) YND = 1.1E16
c	TND = tied(25)
c	IF(TND.GE.1000) TND = 995.
c	CALL INTERP(TED,TID,TND,YND)
c	SVEL = SEL(1)
c	SVELN= SELN(1)
c	SVCX = SCX(1)
c	SVATA(25)= (SEL(1) + CNEUT*SELN(1)*YNO(25)/YNI(25,1))+SCX(1)
c	SVIONA(25) = SION(1) 
c	SVREC = RECOM(1)
c	XNUIONI(25) = YNO(25)*SVIONA(25)
c	VCOLD = SQRT(XK*TSPL/XMASS)
c	ynocold(25) = 0.5*(gamouteff/vcold + coldno(25))
c     XNUATI(25) =  YNOCOLD(25)*SVATA(25)

c	xni(1) = xnsol
c	xni(2) = fracz*xnsol
C	***************TEMPORARY*************************
 
c	DELTAN = (XNPEDEX-XNSOL)/24.
c	DELTATE = (TPEDEXE-TSEPEXE)/24.
c	DELTATI = (TPEDEXI-TSEPEXI)/24.
c	***********************
c	rhor(25) = 1.0
	chiion(25) = xchii(25)
	chiel(25) =  xchie(25)
	
 	DO 1003, NN=1,24
	N = 25-NN
	chiion(n)= xchii(n)
	chiel(n)= xchie(n)
	 	
	vtor1old(n) = torv(n)
	yni(n,2) = fracz*yni(n,1)
	TEED(N) = TEl(N)
	TIED(N) = TI(N) 
	
1003	CONTINUE
	
	   	
	vtor1old(25) = torv(25)
C	*****************************************
	
c	"Outer" iteration converging velocity calculations & density calculations 
c		through 3000
1004	kk=0
1007	do 3000 jt = 1,1
	 mjt = 0 

	CALL NEUTDIST

c	***********end control section****************************************
	
c	*****"Inner" iteration converging ion and neutral distributions through 2000******
 
1500 	DO 2000 IT = 1,1
	XNUIONI(25) = YNO(25)*SVIONA(25)
      XNUATI(25) =  YNOCOLD(25)*SVATA(25)
	XNUDRAGATOMIC(25) = 0.5*(XNUIONI(25)+xnuioni(24)) +	
     1	coldno(25)*svata(25)
     
c	update velocity gradient scale length
	kk = kk+1
	
 	
	xlvm(n) = -2.*(vtor1(n)-vtor1(n-1))/
	1		(delna*(vtor1(n)+vtor1(n-1)))
c	if(ioptxlvm.eq.5) xlvm(n) = -2.*(omegt(n,1)-omegt(n-1,1))/
c    1		(delna*(omegt(n,1)+omegt(n-1,1)))
     	xlvm(1) = xlvm(2)
	do 1010 n = 1,25
	velthet1old(n) = velthet1(n)
	velthet2old(n) = velthet2(n)
	eradaold(n) = erada(n)
1010	continue
 	do 1025 n = 1,25
	ynold(n) = yno(n)
	yniold(n) = yni(n,1)
1025	continue
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
	goto 1028 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
1028	XLradZ(25)= eXLZDBL
	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(25) = XLradZ(25)*1.e-13
	xni(1) = yni(25,1)
	xni(2) = yni(25,2)
	vrad1(25) = gamion(25,1)/yni(25,1)
	vphia(2) = torv(25)
	vphia(1) = vphiex1(25)
 
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
	radmultedge = qzmultcore
	radmultedge = 1.
c******************************
	xs = xnav
c*****note******** 
c	ntorque = 0
c	vphia(1) = vtor1(25)
c	vphia(2) = vtor2(25)
c	call edgerotran(25,20)
c	call poloidal(25) 
c
c	vth(1) = sqrt(2.*xk*ti(25)/xmas(1))
c	vth(2) = sqrt(2.*xk*ti(25)/xmas(2))
c
c	vtheta(1) = vtheta(1)*VTH(1)*fp
c	vtheta(2) = vtheta(2)*VTH(2)*fp
c	velthet1(25) = vtheta(1)
c	velthet2(25) = vtheta(2)
	do 808 j=1,2
 	xz = atnum(1)
	if(j.eq.2) xz = zbar2(25)
	PRESS(J)=-1.*(ti(25)/(xz*BTHET))*xlpm(25)
808	continue

	eradfb(25) = bthet*(torv(25) - vtheta(2)/fp + press(2))

	vpinchi(25) = -1.*xmomtor1(25)/yni(25,1)  +  eq*ephia
	vpinchi(25) = vpinchi(25) +  xmas(1)*
     1 (xnudtot1(25)+xnuc12(25))*(eradfb(25)/bthet+velthet1(25)/fp) 
     	vpinchi(25) = vpinchi(25) - xmas(1)*xnuc12(25)*torv(25)
     	vpinchi(25) = vpinchi(25)/(eq*atnum(1)*bthet)
	 
	diffA(25) =  xmas(1)*xk*ti(25)*xnuc12(25)*
     1			((xnudtot1(25)/xnuc12(25))+1.-zbar2(25)/atnum(1))/
     2			((eq*atnum(1)*bthet)**2)

    	thetw1(25) = thetw(1)
	thetw2(25) = thetw(2) 
	ratdrag(25) = xnudragatomic(25)/ynudrag1(25)
     	
	if(chiion(25).le.0.0) chiion(25) = 1.0
	if(chiel(25).le.0.0) chiel(25) = 1.0
 
c	inverse temperature gradient scale lengths	

	xltim(25) = (gamheati(25)/(yni(25,1)*xk*ti(25))-2.5*vrad1(25))/
     1 	(chiion(25))

	yne = yni(25,1)*atnum(1)+yni(25,2)*zbar2(25)
	vrade=(atnum(1)*yni(25,1)*vrad1(25)+zbar2(25)*yni(25,2)*vrad2(25))
     1	/yne
	xltem(25) = (gamheate(25)/(yne*xk*tel(25))-2.5*vrade)
     1		/(chiel(25))

	xltzm(25) = xltim(25) 
c	inverse pressure and ion density gradient scale lengths
c	if(it.gt.1) diff%(25) = 0.5*(diffA(24)+diffA(25)) 
	xlpm(25) = (vrad1(25)-vpinchi(25))/diffA(25)
c	if(xlpm(25).lt.0.0) xlpm(25) = 0.0
c	if(it.gt.1) xlpm(25) = xlpm(24)
 	c10=1.0

	xlnmold = xlnm(25)
	xlnm(25) = xlpm(25)-c10*xltim(25)
c	if(xlnm(25).lt.0.0) xlnm(25) = 0.0			

c	88888888888temp fix	****************************88888888888
c	goto 1301
		
	DO 1001 NN = 1,24
	cedwmult = 1.0		
	MIT = 0	
	N = 25 - NN
	j=n
	delma = delna

	XNI(1) = YNI(N,1)
	XNI(2) = YNI(N,2)
	TEMP(1) = TI(N)
	TEMP(2) = TI(N)
	TEP = TEl(N)
	vphia(2) = torv(n)
	vphia(1) = vphiex1(n)
	XNUIONI(N) = YNO(N)*SVIONA(N)
	ynocold(n) = 0.5*(coldno(n)+coldno(n+1))
	XNUATI(N) =  YNOCOLD(N)*SVATA(N)
	
	CALL PARAM(N)
	xnuc12(n) = xnuc(1,2)
	xnuc21(n) = xnuc(2,1)
     
c	integrate particle and heat fluxes inward from separatrix
c	new ion & plasma density and ion & electron temp formulation
	dens(n) = 0.5*(yni(n+1,1)+yni(n,1) )
	tele = 0.5*(tel(n)+tel(n+1))
	tiav = 0.5*(ti(n)+ti(n+1))
	rz = (atnum(1)**2/xmas(1) + zbar2(n)**2/xmas(2))*1.67e-27  
	cequil = 6.32e-14*rz 
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)
c	eq 4.90 fpp
	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(yne) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(tele**1.5)*X/Yz)

	cequil = 7.9e-42*couloge*zeff/xmas(1) 

	qie(n) = cequil*yne*(tiav-tele)/(tele**1.5)

	cxcool(n)=1.5*dens(n)*tiav*xk*xnuati(n)
c**********************
c	cmulteq = 0.0

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
c	zbar2(n) = 6.0
	goto 1054 
	endif

	CALL cefits (IZ1, TDBL, XLZDBL, 1, DLZDBL)
1054	XLradZ(n)= eXLZDBL
	
C		CONVERT FROM ERG-CM3/EV-S TO J-M3/S
	XLradZ(n) = XLradZ(n)*1.e-13
	coolion(n) = xk*eioni*yne*xnuioni(n)
	radcool(n) = 0.5*(yni(n,2)+yni(n+1,2))*yne*xlradz(n)
	radcool(n) = radmultedge*radcool(n)

	ntorque = 0
	CALL EDGEROTRAN(N,20)
	vphia(1) = vtor1(n)
	vphia(2) = vtor2(n) 

c**************temp fix********************888
c	call poloidal(n)
c	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
c	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
c 
c	vtheta(1) = vtheta(1)*VTH(1)*fp
c	vtheta(2) = vtheta(2)*VTH(2)*fp
c  	velthet1(n) = vtheta(1)
c	velthet2(n) = vtheta(2)
c	thetw1(n) = thetw(1)
c	thetw2(n) = thetw(2)
c	particle fluxes on ions & impurities
 	thetint = 0.5
	DENS(N) = thetint*YNI(N+1,1) + (1.-thetint)*YNI(N,1)
	GAMION(N,1) = GAMION(N+1,1) -dens(n)*xnuionb(n)*delma	
     1			-	DENS(N)*XNUIONI(N)*DELMA*(1.+fracz*zbar2(n)) 
	GAMION(N,2) = GAMION(N+1,2) 
     1		-	DENS(N)*XNUIONz(N)*DELMA*(1.+fracz*zbar2(n))
	vrad1(n) = gamion(n,1)/yni(n,1)
	vrad2(n) = gamion(n,2)/yni(n,2)
	do 908 j=1,2
 	xz = atnum(1)
	if(j.eq.2) xz = zbar2(n)
	PRESS(J)=-1.*(ti(n)/(xz*BTHET))*xlpm(n)
908	continue

	eradfb(n) = bthet*(torv(n) - vtheta(2)/fp + press(2))
 
c	fraction of heating to ions
	fb1 = 0.75
 	fb2 = 0.15
	fb3 = 0.10

	ecrit = 19.*tel(n)
	xc= sqrt(1.e3*eb/ecrit)
	xx = atan((2.*xc-1.)/1.732)
	yy = log((xc**2+2.*xc+1.)/(xc**2-xc+1.))
	fion1 = 2.*(((xx+.5236)/1.732)-yy/6.)/(xc**2) 
	xc= sqrt(1.e3*(eb/2.)/ecrit)
	xx = atan((2.*xc-1.)/1.732)
	yy = log((xc**2+2.*xc+1.)/(xc**2-xc+1.))
	fion2 = 2.*(((xx+.5236)/1.732)-yy/6.)/(xc**2)
	xc= sqrt(1.e3*(eb/3.)/ecrit)
 	xx = atan((2.*xc-1.)/1.732)
	yy = log((xc**2+2.*xc+1.)/(xc**2-xc+1.))
	fion3 = 2.*(((xx+.5236)/1.732)-yy/6.)/(xc**2) 
	fionb(n) = fb1*fion1 + fb2*fion2 + fb3*fion3
	qnbi(n) = fionb(n)*qnb(n)
	qnbe(n) = (1.-fionb(n))*qnb(n) 
	

   	beamdot(n) =  dens(n)*xnuionb(n)
	
	gamheati(n) = gamheati(n+1) + delma*(cxcool(n) + cmulteq*qie(n))
 	gamheati(n) = gamheati(n) - fionb(n)*qnb(n)*delma	 
	gamheate(n) = gamheate(n+1)+delma*
     1					(coolion(n)+radcool(n)-cmulteq*qie(n))
	gamheate(n) = gamheate(n) -	(1.-fionb(n))*qnb(n)*delma
c	time dependent
c	gamheati(n) = gamheati(n)-dlnwi_dt(n)*1.5*yni(n,1)*xk*ti(n)*delma
c	gamheate(n) = gamheate(n)-dlnwe_dt(n)*1.5*yne*xk*tel(n)*delma
c	gamion(n,1) = gamion(n,1)-dlnn_dt(n)*yni(n,1)*delma
1001  continue
c	adjustment to surface flux boundary conditions
c	to account for time-dependent n and nT in pedestal
	ye=0.
	yi=0.
	yn=0.
	dlnpeped_dt =dlnpped_dt-dlnw_dt
	dlnpiped_dt =dlnpped_dt-dlnw_dt
c******changed 6-22-07 to subtract global dlnn_dt and dlnw_dt
	do 1776 n = 1,24
	yne = 0.5*(yni(n,1)+yni(n+1,1))*atnum(1) +
     1		0.5*(yni(n,2)+yni(n+1,2))*zbar2(n)
c	set time derivatives to input profile values data20
	dlnpeped_dt = dlnwe_dt(n)- dlnw_dt
	ye = ye + dlnpeped_dt*1.5*yne*xk*tel(n)*delma
	dlnpiped_dt = dlnwi_dt(n)- dlnw_dt
	yi = yi + dlnpiped_dt*1.5*yni(n,1)*xk*ti(n)*delma
	dlnnped_dt = dlnn_dt(n)-dln_dt
	yn = yn + dlnnped_dt*yni(n,1)*delma
1776  continue
c*****changed 6-22-07 end ***************

c	adjust heat & particle flux for change in boundary conditions
	do 1777 n = 1,25
	gamheati(n) = gamheati(n) - yi
	gamheate(n) = gamheate(n) - ye
	gamion(n,1) = gamion(n,1) - yn
1777	continue	 
	do 1013 n=1,25 
	vrad1(n) = gamion(n,1)/yni(n,1)
c*******************
	delma = delna 
c*******************
c	assume fracheate, instead of calculating equilibration
	gamheat(n) = gamheate(n) + gamheati(n)

c	new ion & impurity density, ion & electron temp formulation
c	calculate inverse temperature gradient scale lengths
	if(chiion(n).le.0.0) chiion(n) = 1.0
	if(chiel(n).le.0.0) chiel(n) = 1.0

	xltim(n)=(gamheati(n)/(yni(n,1)*xk*ti(n))-2.5*vrad1(n))/
	1			(chiion(n))
	xltzm(n) = xltim(n) 
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	vrade=(atnum(1)*yni(n,1)*vrad1(n)+zbar2(n)*yni(n,2)*vrad2(n))/yne
	xltem(n) = (gamheate(n)/(yne*xk*tel(n))-2.5*vrade)
	1	/(chiel(n))

	vpinchi(n) = -1.*xmomtor1(n)/yni(n,1)  +  eq*ephia 
	vpinchi(n) = vpinchi(n) +
     1 xmas(1)*(xnudtot1(n)+xnuc12(n))*(eradfb(n)/bthet+velthet1(n)/fp)
      vpinchi(n) = vpinchi(n) - xmas(1)*xnuc12(n)*torv(n)
     	vpinchi(n) = vpinchi(n)/(eq*atnum(1)*bthet)
	 
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
     1			((xnudtot1(n)/xnuc12(n))+1.-zbar2(n)/atnum(1))/
     2			((eq*atnum(1)*bthet)**2)

	xlpm(n) = (vrad1(n)-vpinchi(n))/diffA(n)	
c	if(xlpm(n).lt.0.0) xlpm(n) = xltim(n)
	xlnmold = xlnm(n) 
	xlnm(n) = xlpm(n)-c10*xltim(n) 
c	if(xlnm(n).lt.0.0) xlnm(n) = 0.0
1013	continue 
c	integrate ion density inward from separatrix
	c6 = 0.5
	yni(25,1) = exne(25)/(atnum(1)+fracz*zbar2(25)) 
	yni(25,2) = fracz*yni(25,1)
	ti(25) = xti(25)
	tel(25) = xte(25)
 	do 1014 nm=1,24
	n = 25 - nm	
	yni(n,1) = yni(n+1,1)*
	1	(1.+delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))/
     2	(1.-delma*0.5*(c6*xlnm(n)+(1.-c6)*xlnm(n+1)))
	if(yni(n,1).le.0.0) yni(n,1) = 1.e19 
	if(yni(n,1).gt.2.e20) yni(n,1) = 2.e20
	yni(n,2) = fracz*yni(n,1)
	yne = yni(n,1)*atnum(1)+yni(n,2)*zbar2(n)
	zne(n) = yne 
	vrad1(n) = gamion(n,1)/yni(n,1)
	vrad2(n) = gamion(n,2)/yni(n,2)

c	integrate ion & electron temperatures inward from separatrix
	tiold(n) = ti(n)
	teold(n) = tel(n)
	ti(n) = ti(n+1)*
     1	(1.+delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))/
     2	(1.-delna*0.5*(c6*xltim(n)+(1.-c6)*xltim(n+1)))
	
      if(ti(n).lt.0.0) ti(n) = 100.
	tel(n) = tel(n+1)*
     1	(1.+delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))/
     2	(1.-delna*0.5*(c6*xltem(n)+(1.-c6)*xltem(n+1)))
	if(tel(n).lt.0.0) tel(n) = 100.

C	INTEGRATE TOROIDAL ROTATION FREQ INWARD FROM SEPARATRIX

c	goto 1164
c	if(n.eq.24) call torotate(25)
c	do 1145 mm = 1,6 
c	s0v(mm,25) = sv(mm,25)
c1145  continue
c	CALL TOROTATE(N)
	
c      DO 1154 mm = 1,6
c	s0v(mm,j) = sV(mm,j)
c	DO 1149 kk = 1,6
c	b0v(mm,kk,j) = bv(mm,kk,j) 
c	SV(mm,j) = SV(mm,j) + (AV(mm,kk,j)/delma)*OMEGT(j+1,kk)
c	BV(mm,kk,j) = BV(mm,kk,j) +(AV(mm,kk,j)/DELMA)
c	BBV(mm,kk) = BV(mm,kk,j)
c1149	CONTINUE 
c	SSV(mm) = SV(mm,j)
c1154  CONTINUE
c	CALL LSLRG(6,BBV,6,SSV,1,SSV)
c	DO 1159 kk = 1,6
c	OMEGT(j,kk) = SSV(kk)
c1159	CONTINUE	
c1164	continue
	
 
1014	CONTINUE
		
	call neutdist
	
c	converge on main ion density and temp
	do 1170 n = 1,24	
 	IF(ABS(YNIOLD(n)/YNI(N,1)-1.).GT.0.02) MIT = MIT + 1
	IF(ABS(tIOLD(n)/tI(n)-1.).GT.0.02) MIT = MIT + 1
 	if(abs(teold(n)/tel(n)-1.).gt.02) mit = mit+1 

	IF(YNI(N,1).LE.0.0) YNI(N,1) = 1.E19
	IF(YNI(N,2).LE.0.0) YNI(N,2) = 1.E19*fracz
1170	continue 
	iF(MIT.EQ.0) GOTO 1225

	nj = 50
	
	if(iconverge.eq.0.and.it.eq.nj) then
	write (6,1199) nj
	jwarn = 1
	goto 1225
	endif
1199	format(1x,'not converged on 2000 loop after iterations=',I3.0) 
2000	CONTINUE

1225	continue	
c	converge on poloidal velocities ********************

c*********************************************************************************
	do 1275 n = 1,25
	errthet(n) = abs((velthet1old(n)/velthet1(n))-1.)
	if(errthet(n).gt.0.05) mjt = mjt + 1
	errthet(n) = abs((velthet2old(n)/velthet2(n))-1.)
 	if(errthet(n).gt.0.05) mjt = mjt + 1 
	

	velthet1(n) = 0.5*(velthet1(n)+velthet1old(n))
	velthet2(n) = 0.5*(velthet2(n)+velthet2old(n))
	velthet1old(n) = velthet1(n)
	velthet2old(n) = velthet2(n)
 
 
1275	continue

	if(mjt.ne.0) goto 1298
	 
	goto 1301

1298	continue
			
3000	continue		 
	
1301	continue
	OPEN(unit=123,FILE='pedestal2.TXT',STATUS='UNKNOWN') 
	if(jwarn.eq.1) write(123,199)

	write(123,'(1x,35A)') '   rho	     ne         ti         te     
     1  vrad    vpinch  -(dp/dr)/p '
	
c	sepdif = xnsepex - xnesepreal
     	do 3750 n = 1,25
	gamconve = 2.5*xk*zeff*gamion(n,1)*xte(n)
	zne(n) = zne(n) - sepdif
	write (123,1000) rhor(n),zne(n),ti(n),tel(n),vrad1(n),vpinchi(n),
	1	xlpm(n)
	
3750	continue
	write(123,'(1x,35A)') '   rho	    nudrag1*    nuatom1   nudrag2*       
     1nugyro1  nugyro2  gamion  '
     	do 3760 n = 1,25
	write (123,1000) rhor(n),xnudtot1(n),xnuati(n)+xnuioni(n)+xnuionb(
     1n),xnudtot2(n),xnudragyro1(n),xnudragyro2(n),gamion(n,1)
	dragfreq(n) = xnudtot1(n) 
3760	continue
	write(123,'(1x,35A)') '   rho	 -(dp/dr)/p -(dTi/dr)/T -(dn/dr)/n    
     1-(dTe/dr)/T  xnuion  qnbi '
   
	do 3770 n=1,25
	write (123,1000) rhor(n),xlpm(n),xltim(n),xlnm(n),
     1	xltem(n),ynuion(n),qnbi(n)
3770	continue 

	write(123,'(1x,35A)') '   rho     cosion    cosimp    sinion  
	1 sinimp	  beta-1    beta-2  '   
 
	do 3775 n = 1,25
	ynud(1) = xnudtot1(n)/xnuc12(n)
	ynud(2) = xnudtot2(n)/xnuc21(n)

	write (123,1000) rhor(n),cosion(n),cosimp(n),sinion(n),sinimp(n),
     1  ynud(1),ynud(2) 	
3775	continue 
	write(123,'(1x,35A)') '   rho     vthet1     vthet2      vthexp  
	1eradfb  eradexp'
      do 3795 n = 1,25
	write (123,1000) rhor(n),velthet1(n), velthet2(n),vthexp(n),
     1eradfb(n),erex(n)    
3795	continue
	WRITE(123,'(1X,35A)')'  rho      vphi1      vphi2   vphi2ex     vp
     1hi1ex, nustar-iz  nustar-zi'

	DO 3820, N=1,25
		do 3816 j=1,2
 	VTH(J) = SQRT(2.*XK*TI(N)/XMAS(J))
3816	continue 
      xxxiz = xnuc12(n)*rmajor*qsafe/vth(1)		  	 
	xxxzi = xnuc21(n)*rmajor*qsafe/vth(2)	 

	WRITE(123,1000)rhor(n),vtor1(N),vtor2(N),torv(N),vphiex1(n), 
     1	xxxiz,xxxzi
     	
3820	CONTINUE
	write(123,'(1x,35A)') '   rho	   gamheate    gamheati   gamion     
     1 econve    econvi  gamconvi  '
     	do 3755 n = 1,25
	gamconvi = 2.5*xk*gamion(n,1)*xti(n)

	write (123,1000) rhor(n),gamheate(n),gamheati(n),gamion(n,1),
     1econve(n),econvi(n),gamconvi
3755	continue 
	
4099	RETURN

	END
