	SUBROUTINE EDGEROTRAN(nmesh,ntorque)
	INCLUDE 'SOLDIV.FI'
	PARAMETER (JQ=2) 
	DIMENSION fluxpav(JQ), y(JQ), ynud(jQ),ynuconv(jQ),vphical3(jQ),
	2			vpinch(jQ),VPHICAL(JQ),VPHICAL2(JQ),grad(jq)

	if(ioptapproach.eq.0) then
c	erad= eradex, vphi1=vphi2=vphi2ex, vthet1&vthet2 calculated;
c	infer nudrag to match vphiex. calc velocities using exp nudrag
	iopterad = 0
	ioptdrag = 2
	ioptvisc = 0
	ioptpinchi =2
	ioptvphi = 2
	ioptvdif = 2
	ioptxvlm =3
	ioptvthet = 5
	ioptexpdata = 1
	endif

	
     	if(ioptapproach.eq.1) then
c	erad= eradex, vphi1=vphi2=vphi2ex, vthet1=vthet2=vthet2ex;
c	infer nudrag to match vphiex. calc velocities using exp nudrag
	iopterad = 0
	ioptdrag = 2
	ioptvisc = 0
	ioptpinchi =2
	ioptvphi = 3
	ioptvdif = 2
	ioptxvlm =3
	ioptvthet = 0
	ioptexpdata = 1
	endif
	
	if(ioptapproach.eq.2) then
c	erad=eradex, vphi2=vphi2ex, vphi1=vphi2ex+difvphical, vthet2=vthet2ex,
c	vthet1=vthet2ex+difvthetcal; infer nudrag to match vphiex. calc vel using nudrag  
	iopterad = 0
	ioptdrag = 2
	ioptvisc = 0
	ioptpinchi = 2
	ioptvphi = 2
	ioptvdif = 2
	ioptxvlm =3
	ioptvthet = 1
	ioptexpdata = 1		 
	endif


	if(ioptapproach.eq.4) then
c	use exp erad, vphi and vtheta; calc velocities using gyro nudrag

	iopterad = 0
	ioptdrag = 2
	ioptvisc = 1
	ioptpinchi =2
	ioptvphi = 5
	ioptvdif = 2
	ioptxvlm =3
	ioptvthet = 0
	ioptexpdata = 1
	endif

	if(ioptapproach.eq.5) then
c	uses all calculated quantities
  	iopterad = 2
	ioptdrag = 2 
	ioptpinchi = 5
	ioptvphi = 5
	ioptvdif = 1
	ioptxvlm =3
	ioptvthet = 5
	ioptexpdata = 0
	endif

     
      

	if(ioptapproach.eq.3) then
	iopterad = 2
	ioptvthet = 5
	ioptdrag = 2
	ioptpinchi = 2
	ioptvphi =5
	ioptvdif = 2
	ioptxvlm =3	   
	endif

	
	
	

 	n = nmesh
	atnum2 = atnum(2)
	atnum(2) = zbar2(nmesh)	
c	if(nmesh.eq.25) goto 44
C	PHI MOMENTUM DEPOSITION IN PEDESTAL
	if(nmesh.eq.25) then
	unatten1 = 1.0 
	unatten2 = 1.0 
	unatten3 = 1.0 
	endif
	IF(NMESH.EQ.25) rRMINOR = AMINOR + deln
c	rminor is plasma radius for effective cylinder, aminor is actual plasma radius (0.6m)
c	bfield = abs(bphi)
	se = sqrt(0.5*(1.+elong**2))
	XMBEAM = ABEAM*1.673E-27 
	fb1 = 0.75
	fb2 = 0.15
	fb3 = 0.10

c	torque
	TORQUE1 = SQRT(2.*XMBEAM/(XK*EB*1.E3))*RTAN*fb1*PBEAM*1.E6
	TORQUE2 = SQRT(2.*XMBEAM/(XK*(EB/2.)*1.E3))*RTAN*fb2*PBEAM*1.E6 
	TORQUE3 = SQRT(2.*XMBEAM/(XK*(EB/3.)*1.E3))*RTAN*fb3*PBEAM*1.E6
	GAMMA = 0.5
	ZEFF = (XNi(1)*(ATNUM(1)**2)+XNi(2)*(zbar2(n)**2))/
	1		(XNi(1)*ATNUM(1)+XNi(2)*zbar2(n))
	 
c	attenuation	
	XLAM1 = 5.5E17*EB/(ABEAM*0.5*DENS(NMESH)*(ZEFF**GAMMA))
	ATTEN1 = 1. - EXP(-1.*(DELMA/COS(ALPHAIN))/XLAM1)
	XLAM2 = 5.5E17*(EB/2.)/(ABEAM*0.5*DENS(NMESH)*(ZEFF**GAMMA))
	ATTEN2 = 1. - EXP(-1.*(DELMA/COS(ALPHAIN))/XLAM2)
 	XLAM3 = 5.5E17*(EB/3.)/(ABEAM*0.5*DENS(NMESH)*(ZEFF**GAMMA))
	ATTEN3 = 1. - EXP(-1.*(DELMA/COS(ALPHAIN))/XLAM3)

c	toroidal momentum input
	xmphi1 = unatten1*atten1*torque1/rmajor
	xmphi2 = unatten2*atten2*torque2/rmajor
 	xmphi3 = unatten3*atten3*torque3/rmajor

 
	x = cos(alphain)
	PI = 3.1416
	rRMINOR = rRMINOR - DELMA			
	VOLM = 4.*(PI**2)*RMAJOR*rRMINOR*DELMA*se
	PEDMPHIDEN1 = XMPHI1/VOLM 
	PEDMPHIDEN2 = XMPHI2/VOLM 
	PEDMPHIDEN3 = XMPHI3/VOLM 
	YZY = PEDMPHIDEN1 + pedmphiden2 + pedmphiden3
	phimom(n) = yzy
c	neutral beam particle & energy sources
	xpartdot1 = (unatten1*atten1*((fb1*pbeam*1.e6)/(xk*eb*1.e3)))/volm
	xpartdot2 = (unatten2*atten2*((fb2*pbeam*1.e6)/
	1			(xk*(eb/2.)*1.e3)))/volm
	xpartdot3 = (unatten3*atten3*((fb3*pbeam*1.e6)/
	1			(xk*(eb/3.)*1.e3)))/volm
 	xpartdot = xpartdot1+xpartdot2+xpartdot3
	xnuionb(nmesh) = xpartdot/dens(nmesh)
	qnb(nmesh) = (xpartdot1+xpartdot2/2.+xpartdot3/3.)*xk*1.e3*eb

c	fraction of heating to ions
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
c	unattenuated beam remaining after crossing region n
	unatten1 = unatten1*(1.-atten1)
	unatten2 = unatten2*(1.-atten2)
	unatten3 = unatten3*(1.-atten3)

400	format(1x,7e10.3) 
c	momentum distributed proportion nZ^2	
 	XX = (zbar2(n)**2)*XNi(2)/
	1		((zbar2(n)**2)*XNi(2)+(ATNUM(1)**2)*XNi(1))
c	momentum distributed equally to all ions	 
	xx = xni(2)/(xni(1)+xni(2))
	XMTOR(1) = (1.-XX)*YZY
	XMTOR(2) = XX*YZY
c	include anomalous torque
	xmtor(1) = xmtor(1) + anomom(n)*(1.-xx)
	xmtor(2) = xmtor(2) + anomom(n)*xx 
	xmomtor1(n) = xmtor(1)
	XMOMTOR2(N) = XMTOR(2)
	xmompol(n,1)=(abs(bthet*bfield)/(bfield**2+bthet**2))*xmomtor1(n)
	xmompol(n,2)=(abs(bthet*bfield)/(bfield**2+bthet**2))*xmomtor2(n)

	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
 	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))

	if(ntorque.eq.10) goto 756 
	if(ntorque.eq.20) goto 756

C	POLOIDAL VELOCITIES
44	R2 = 0.01
c	sputtering yield for C in the 100-1000eV range is about .02,
c	allow for half of this to be trapped in divertor , R2 = 0.01
c	xnuatI(2) = xnocoldav*1.e-13/VTH(2)
c	carbon-hydrogen charge exchange

c	radial velocities
c	fluxpartav = fluxpart - 0.5*delped*xni(1)*xnuionI(1) 
	xnuhat1 = Aminor*xnuatI(1)/VTH(1)
	xnuhat2 = Aminor*xnuatI(2)/VTH(2)
c	x1 = atnum(1)*eq*bfield/(xmas(1)*VTH(1))
c	y1 = qsafe*rmajor/VTH(1)
c	vr1 = fluxpartav/xnbarx
c	vr1hat = x1*y1*vr1
c	x2 = atnum(2)*eq*bfield/(xmas(2)*VTH(2))
c	y2 = qsafe*rmajor/VTH(2)
c	vr2 = r2*fluxpartav/(fracz*xnbarx)
c	vr2hat = x2*y2*vr2
c	fluxpav(1) = fluxpartav
c	fluxpav(2) = r2*fluxpartav	
c		momentum convection
c	xnudconv = ((vr1+fracz*(xmas(2)/xmas(1))*vr2)/
c	1	(1.+	fracz*(xmas(2)/xmas(1))))/delped
595	FORMAT(1X,'OK')
596	format(7e10.3)	
c	 write(6,596) vr1hat,vr2hat
c	***************pressure gradient************************
	
      if(xlnm(nmesh).le.0.0) xlnm(nmesh) = 0.0
	if(xlpm(nmesh).le.0.0) xlpm(nmesh) = xltim(nmesh)			
c	if(xlnm(nmesh).le.0.0) 	then
c	xlnm(nmesh) = 1./xlntop
c	if(rhor(nmesh).gt.pedrhon) xlnm(nmesh) = 1./yln
c	endif
c	if(xlpm(nmesh).le.0.0) then
c	xx = 1./xltitop
c	if(rhor(nmesh).gt.pedrhon) xx = 1./ylti
c	xlpm(nmesh) = 1./yln + 1./xx
c	endif 
	xlp = 1./xlpm(nmesh) 

	do 908 j=1,2
 	xz = atnum(1)
	if(j.eq.2) xz = zbar2(n)
	PRESS(J)=-1.*(ti(nmesh)/(xz*BTHET))*xlpm(nmesh)
908	continue
	press1(nmesh) = press(1)
	press2(nmesh) = press(2) 
	if(nmesh.eq.25) presvth1 = press(1)/vth(1)
c	***********************************************
c	***********toroidal velocities***************************
	if(ioptvphi.eq.1) then
c	use input constant value at all mp	 
	vphia(1) = vphiApedx
	vphia(2) = vphiApedx
 	endif
	ynud2 = ynudrag2(n)/xnuc(2,1)
	if(ioptvphi.eq.2) then
c	use exp value for vphi2 and calculate vphi1 from difference algorithm	 
	vphia(2) = torv(nmesh)
	if(ioptvdif.eq.1)
c	toroidal mom bal
     1	vphia(1) = xnuc(1,2)*vphia(2)/(xnuc(1,2)+ynudrag1(n)) +
	1	(atnum(1)*eq*(ephia+vrad1(n)*bthet)+(xmtor(1)/xni(1)))/
     2	(xmas(1)*(xnuc(1,2)+ynudrag1(n)))
	if(ioptvdif.eq.2)
c	radial mom bal	 
     1	vphia(1) = vphia(2)*(1.+ynud2) -
	2  (xmtor(2)+xni(2)*zbar2(n)*eq*ephia)/(xni(1)*xmas(1)*xnuc(1,2))
 
c     	vphia(1) = vphia(2) + (velthet1(n)-velthet2(n))/fp + 
c	1			(ti(n)/bthet)*xlpm(n)*(1.- 1./zbar2(n))
	endif
	if(ioptvphi.eq.3) then
c	use measured carbon vphi for carbon and deuterium
	vphia(1) = torv(nmesh)
	vphia(2) = torv(nmesh)
	endif

	if(ioptvphi.eq.4) then
c	use calculated values for vphi1 and vphi2
	vphia(1) = vtor1(nmesh)
	vphia(2) = vtor2(nmesh)
	if(it.eq.1) then
	vphia(1) = vphiApedx
	vphia(2) = vphiApedx
	endif
	endif
	if(ioptvphi.eq.5) then
	vphia(1) = torv(nmesh)
 	vphia(2) = torv(nmesh)

c	use new calculated values for vphi1 and vphi2
c	vphia(1) = rmajor*omegt(nmesh,1)
c	vphia(2) = rmajor*omegt(nmesh,4)
	if(it.lt.1) then
	vphia(1) = vphiApedx
	vphia(2) = vphiApedx
	endif
	
	endif	 
	
	
c	************toroidal mom balance NUDRAG*************
	epd = 1.0 
	oldtot1 = xnudtot1(n)
	oldtot2 = xnudtot2(n) 
	if(ioptdrag.eq.2) then 
c	****vphi dif from rad mom bal *******************************
c	bracket =  ti(n)*xlpm(n)*(1.-atnum(1)/zbar2(n))/(atnum(1)*bthet)
c     1			(velthet1(n)-velthet2(n))/fp
c	****vphi dif from tor mom bal ********************************
c	bracket = (atnum(1)*eq*(bthet*vrad1(n)+ephia) + xmtor(1))/
c	1		(xni(1)*xmas(1)*(xnuc(1,2)+xnudtot1(n))) -
c	2		(xnudtot1(n)*vtor1(n))/(xnuc(1,2)+xnudtot1(n))
c	****vphi dif = 0 **********************************************
c	if(ioptvphi.eq.5) bracket = 0.0 
c	bracket = velthet1(n)/fp + erada(n)/bthet + 
c     1	ti(n)*xlpm(n)/(atnum(1)*bthet)
c	******vphi dif from tor mom bal, using exp vphi_I**************
	bracket = 0.0
c*********************************set nudrag option***********************
	jdrag = 1
	epd = 1.0
c***************************************************************************
	mn = 0

635	if(jdrag.eq.1) then
c	calc single nudrag from combined tor mom eqs
	bracket = 0.0	 
	xnum = xni(1)*eq*atnum(1)*bthet*vrad1(n) + xmtor(1) + 
	1		 xmtor(2) + (xni(1)*atnum(1)+xni(2)*zbar2(n))*eq*ephia 
	xdenom = (xni(1)*xmas(1)+xni(2)*xmas(2))*torv(n) + 
     1		xni(1)*xmas(1)*bracket 
      xnud0 = (xnum/xdenom) 
c	correct for vph1 not equal vphi2
	bracket = (xnuc12(n)*y1-xnud0*torv(n))/(xnuc12(n)-xnud0)
	xz = (xnuc12(n)*y1-xnud0*torv(n))/((xnuc12(n)+xnud0)*torv(n))
	xnud1(n) = xnud0/(1.+(xni(1)*xmas(1)/
	1		(xni(1)*xmas(1)+xni(2)*xmas(2)))) 
	xnud2(n) = xnud1(n)
	endif
	
	if(jdrag.eq.2) then
c	calc individual nudrags from respective tor mom eqs
	xnum = xni(1)*eq*atnum(1)*bthet*vrad1(n) + xmtor(1) + 
     1	xni(1)*atnum(1)*eq*ephia - xni(1)*xmas(1)*xnuc(1,2)*bracket
	xdenom = xni(1)*xmas(1)*(bracket + torv(n))
	xnud1(n) = epd*xnum/xdenom + (1.-epd)*oldtot1
	xnum = xmtor(2) + xni(2)*zbar2(n)*eq*ephia + 
     1			xni(1)*xmas(1)*xnuc(1,2)*bracket
	xnud2(n) = epd*xnum/(xni(2)*xmas(2)*vphia(2)) + (1.-epd)*oldtot2 
	endif
	if(jdrag.eq.3) then
c	use input nudrag
	xnud1(n) = dragfreq(n)
	xnud2(n) = dragfreq(n)
	endif
c	calculate inferred experimental deuterium toroidal rotation
c	from toroidal momentum balance
	vphiex2(n) = torv(nmesh)
c	vphiex1(n) = vphiex2(n)
	bracket =  (xnuc12(n)*y1-xnud1(n)*torv(n))/(xnuc12(n)-xnud1(n))
	vphiex1(n) = vphiex2(n) + bracket
	
	mn = mn+1
	if(mn.lt.2) goto 635
	brack(n) = bracket

	brackrad = (velthet1(n)-velthet2(n))/fp + 
     1			ti(n)*xlpm(n)*(1.-atnum(1)/zbar2(n))/(atnum(1)*bthet)
	
	xnudragatomic(nmesh) = xnuioni(n) + xnuati(n) +xnuionb(n)
 
	delb = delbsep*exp((rhor(n)-1.0)*aminor/0.009) 
	xnudragrip(n,1) = 1.25*24.*(delb**2)*vth(1)/rmajor 
	xnudragrip(n,2) = 1.25*24.*(delb**2)*vth(2)/rmajor 

	xnudtot1(n) = xnud1(n) 
	xnudtot2(n) = xnud2(n)
	endif
c	***************************************************

	if(nmesh.eq.25) vphivth1 = vphia(1)/vth(1)

c	************************************************ 


c	*********poloidal velocities**************************
c***********
	if(ioptvthet.eq.0) then
	velthet1(nmesh) = vthexp(nmesh)
	velthet2(nmesh) = vthexp(nmesh) 
	vtheta(1) = velthet1(nmesh)
	vtheta(2) = velthet2(nmesh)
	goto 752
	endif
750	continue 
	if(ntorque.eq.1) then
c	do not do vtheta calc, but update xnudrag for xlpm if ntorque=1

	vtheta(1) = velthet1(nmesh)
	vtheta(2) = velthet2(nmesh)
	thetw(1) = thetw1(nmesh)
	thetw(2) = thetw2(nmesh) 
	do 751 j=1,2
	zz = ATNUM(J)
 	if(j.eq.2) 	zz = ZBAR2(N)
	temp(j) = ti(n)
	grad(j) = rhor(n)*AMINOR*(xlpm(n)+xlvm(n))

	XNUDRAG(J)=GRAD(J)*THETw(J)*TEMP(J)/
     2			(2.*(RMAJOR**2)*zz*(BPHI))
	
	
 
751	continue
	
c	goto 752
	endif 
c*************
	 
c	call poloidal(nmesh)
	
c	vtheta(1) = vtheta(1)*VTH(1)*fp
c	vtheta(2) = vtheta(2)*VTH(2)*fp
	
c	vtheory(n) = vtheta(2)
c	thetw1(nmesh) = thetw(1)
c	thetw2(nmesh) = thetw(2)
c	if(ioptvthet.eq.1) then
c		difvtheta = vtheta(1) -vtheta(2)
c		vtheta(2) = vthexp(n)
c	
c		vtheta(1) = vtheta(2) + difvtheta
c	endif 
c	velthet1(nmesh) = vtheta(1)
c	velthet2(nmesh) = vtheta(2)
c	vpol(1,nmesh) = vtheta(1)
c 	vpol(2,nmesh) = vtheta(2)

c	gyroviscous
 	xnudragyro1(n) = xnudrag(1) 
      xnudragyro2(n) = xnudrag(2) 
	AM = (aminor*SQRT(0.5*(1.+ELONG**2)))
	ep = am*rhor(n)/rmajor 
 	rminor=am*rhor(n)
 
c	inertial terms
c	xnuinert1(n) = (vrad1(n)/rmajor)*(1.-rmajor*xlvm(n)) -
c	1	0.5*(ep*velthet1(n)/rmajor)*thetinert(1) 
c	xnuinert2(n) = (vrad2(n)/rmajor)*(1.-rmajor*xlvm(n)) -
c    	1	0.5*(ep*velthet2(n)/rmajor)*thetinert(2)
c	convection terms 9/30/05
c	xnuinert1(n)=(1./rminor)+(xnuioni(n)+xnuionb(n))/vrad1(n)-xlvm(n)
c	xnuinert2(n)=0.0
	
c*****xnuinert calculated in poloidal*******************   	
c	anomalous drag
	xnudraganom1(n) = xnudtot1(n)-(xnuioni(n) + xnuati(n) +xnuionb(n))
	1			-xnudragyro1(n)-xnuinert1(n)-xnudragrip(n,1)
	xnudraganom2(n) = xnudtot2(n)-xnudragyro2(n)-xnuinert2(n)-
     1			xnudragrip(n,2)

c	if(xnudraganom1(n).lt.0.0) xnudraganom1(n) = 0.0
c	if(xnudraganom2(n).lt.0.0) xnudraganom2(n) = 0.0


752	continue
	if(ioptdrag.eq.2) then
 	xnudrag(1) = xnud1(n)

	xnudrag(2) = xnud2(n)
c	gyroviscous
c 	xnudragyro1(n) = xnudrag(1) 
c      xnudragyro2(n) = xnudrag(2) 

	XNUDRAGvis1(nmesh)=XNUDRAG(1)
	XNUDRAGvis2(nmesh)=XNUDRAG(2)
c	xnudrag(1) = xnud1(n) + xnuioni(n) + xnuati(n) + xnuionb(n)
	endif
 	do 753 j=1,2
	if(xnudrag(j).le.0.0.and.ioptzerodrag.eq.0) xnudrag(j) = 0.0
753	continue
		

	
	
c	convective momentum transfer frequency
c	ynuconv(1) = vr1/delped
c	ynuconv(2) = vr2/delped 
	
	YNUDRAG1(NMESH) = XNUDRAG(1)
 	YNUDRAG2(NMESH) = XNUDRAG(2)
 
	do 233 j = 1,2
c	fluxpav(j) = 0.5*(gamion(nmesh,j)+gamion(nmesh+1,j))
c	y(j) = atnum(j)*eq*bthet*fluxpav(j) + xmtor(j)  
     
c	k=2
c	if(j.eq.2) k=1

c		update nudrag to include atomic physics
c	ynuioni(nmesh) = xnuioni(nmesh)
c	ynuati(nmesh)  = xnuati(nmesh)
c	average no in mesh intervals to get no at mesh points
c	if(kk.gt.1) then
c     if(nmesh.ne.1) 
c     1		ynuioni(nmesh) = 0.5*(xnuioni(nmesh)+ xnuioni(nmesh-1))

c   	ynuati(nmesh)  = 0.5*(xnuati(nmesh) + xnuati(nmesh-1))
c	endif
c	if(j.eq.1) xnudrag(j)=xnudrag(j) + coldno(nmesh)*svata(nmesh) + 
c     1	0.5*(ynuionI(nmesh) + ynuionI(nmesh+1))
	xnudrag(j) = xnudrag(j) + xnudraganom(j)
	ynud(j) = xnudrag(j)/xnuc(j,k) 
	 
233	continue
	ynud(1) = xnud1(n)/xnuc(1,2)
	ynud(2) = xnud2(n)/xnuc(2,1)
	 
	if(ioptdrag.eq.1) then
	YNUDRAG1(NMESH) = XNUDRAG(1)+ xnudraganom(1)
	YNUDRAG2(NMESH) = XNUDRAG(2)+ xnudraganom(2)
	endif
	ratnu(n) = xnudrag(1)/xnuc(1,2)

 	xnudragatomic(nmesh) = xnuioni(n) + xnuati(n) + xnuionb(n)
	fac1 = (1.+ynud(1))*(1.+ynud(2)) - 1. 
	vphical1old = vphical(1)
	vphical2old = vphical(2)
	if(it.eq.1) then
	vphical1old = torv(n)
	vphical2old = torv(n)
	endif
	
	 

	vphical(1) = ((xmtor(1)+xmtor(2))/(xni(1)*xmas(1)*xnuc(1,2))+
	1  ynud(2)*((velthet1(n)-velthet2(n))/fp - (press(1) -press(2))))/
	2	(ynud(1)+ynud(2))
	 vphical(1) = vphical(1) + (atnum(1)*eq*(bthet*gamion(n,1) +
     2		xni(1)*ephia)+zbar2(n)*eq*xni(2)*ephia)/
     3		(xni(1)*xmas(1)*xnuc(1,2)*(ynud(1) + ynud(2))) 
	vphical(2) = (1.+ynud(1))*vphical(1) -
     1	(atnum(1)*eq*(bthet*gamion(n,1)+xni(1)*ephia) + xmtor(1))/
     2	(xni(1)*xmas(1)*xnuc(1,2))
     
      y1 = xmtor(1) + atnum(1)*eq*(xni(1)*ephia + bthet*gamion(n,1))
	y1 = y1/(xni(1)*xmas(1)*xnuc(1,2))	
	y2 = xmtor(2) + zbar2(n)*eq*(xni(2)*ephia + bthet*gamion(n,2))
	y2 = y2/(xni(2)*xmas(2)*xnuc(2,1))
	vphical(1) = ((1.+ynud(2))*y1+y2)/((1.+ynud(1))*(1.+ynud(2))-1.)
	vphical(2) = (vphical(1)+y2)/(1.+ynud(2)) 

c	do 235 j=1,2
c	k=2
c	if(j.eq.2) k=1
c	vphical(j) = ((1.+ynud(k))*y(j) + y(k))/ 
c   1				(xni(j)*xmas(j)*xnuc(j,k)*fac1)
c	vphical2(j) = eradcomb/bthet - press(j) + vtheta(j)/fp
c	vphical3(j) = y(j)/(xni(j)*xmas(j)*xnuc(j,k)*ynud(j))
c235	continue
c	vphidif = (ynud(2)*y(1) - ynud(1)*y(2))/
c	1			(xni(2)*xmas(2)*xnuc(2,1)*fac1)
c     vphical3(2) = vphical3(1) + y(2)/(xni(2)*xmas(2)*xnuc(2,1)) 
c	write(6,596) vphical3(1),vphical3(2),y(1),y(2),vphipedx/VTH(2)
	adj = 1.0 
	vtor1(nmesh) = adj*vphical(1) + (1.-adj)*vphical1old
	vtor2(nmesh) = adj*vphical(2) + (1.-adj)*vphical2old
45	continue 

c	radial electric field--summed over species 
	xnum = 0.
	denom = 0.
	do 236 j=1,2
	xnum = xnum + xmtor(j) + xni(j)*xmas(j)*xnudrag(j)*(press(j) -
     1		vtheta(j)/fp)
      denom = denom + xni(j)*xmas(j)*xnudrag(j)
236	continue
c	write(6,596) xnudrag(1),xnudrag(2),xnudraggyro
	
	fluxsum = (0.5*(gamion(nmesh,1)+gamion(nmesh+1,1))*atnum(1) +
	1		  0.5*(gamion(nmesh,2)+gamion(nmesh+1,2))*atnum(2)	)*
	1			eq*bthet
      beamsum(NMESH) = xmtor(1) + xmtor(2)
	pressum(NMESH) = xni(1)*xmas(1)*xnudrag(1)*press(1) +
     1	 	   xni(2)*xmas(2)*xnudrag(2)*press(2) 
	vpolsum(NMESH) = xni(1)*xmas(1)*xnudrag(1)*vtheta(1)/fp +
     1	 	   xni(2)*xmas(2)*xnudrag(2)*vtheta(2)/fp 
	dragsum(nmesh) = xni(1)*xmas(1)*xnudrag(1) +
     1	 	   xni(2)*xmas(2)*xnudrag(2) 

C	xnum = xnum + fluxsum 
	eradcomb = bthet*xnum/denom		 
	ERADA(NMESH) = 	((BEAMSUM(NMESH)+PRESSUM(NMESH)-VPOLSUM(NMESH))/
	1				 DRAGSUM(NMESH))*BTHET
C		electric field from force balance
c	eradfb = bthet*(vphipedx - vtheta(1)/fp + press(1))
	
	eradfb(n) = bthet*(vphia(2) - vtheta(2)/fp + press(2))
c	eradfb(n) = 0.5*eradfb(n) + 0.5* bthet*(vphia(1) - vtheta(1)/fp +
c     1	 press(1))

	eradfbnew(n) = bthet*(rmajor*omegt(n,2) - vpol(2,n)/fp + press(2))


	eradfb1(n) = bthet*(vtor1(n) - vtheta(1)/fp + press(1))
	eradfb2(n) = bthet*(vtor2(n) - vtheta(2)/fp + press(2))
	erad3(n) = (xmtor(1)+xmtor(2) - xni(2)*xmas(2)*xnudrag(2)*torv(n))
	1			/(xni(1)*xmas(1)*xnudrag(1)) +
	2		   vtheta(1)/fp + ti(n)*xlpm(n)/bthet
	 
c
	erad6(nmesh) = erada(n)
	if(iopterad.eq.2) erada(nmesh) = eradfb2(n)
	if(iopterad.eq.3) erada(nmesh) = eradfb1(n)
	if(iopterad.eq.0) erada(nmesh) = erex(n)
	if(iopterad.eq.4) erada(nmesh) = erad6(n)
	
C	ERADA(NMESH) = ERADCOMB 
	
	
1005	format(1x,'chixpi=',f5.2,1x,'chixpe=',f5.2,1x,'chitop=',f5.2,1x,
     1	'chetop=',f5.2) 
C	********************COEFFICIENTS FOR EDGE ION DISTRIBUTION CALC********

	n = nmesh
	xni(1) = yni(nMESH,1)
	xni(2) = yni(nmesh,2)

	diffii(n) = xmas(1)*xk*ti(1)*(xnuc(1,2)+xnud1(n))/
     2			((eq*atnum(1)*bthet)**2)
	diffiz(n) =	xmas(1)*xk*ti(1)*xnuc(1,2)/((eq*atnum(1)*bthet)**2)
	diffiz(n) = diffiz(n)/zbar2(n)
	diffzz(n) = xmas(2)*xk*ti(n)*(xnuc(2,1)+xnud2(n))/
	2			((eq*zbar2(n)*bthet)**2)
	diffzi(n) =	xmas(2)*xk*ti(1)*xnuc(2,1)/((eq*atnum(1)*bthet)**2)
 	diffzi(n) = diffzi(n)/zbar2(n)

c	goto 250
c	old density formulation
	
	
 	xlt1(n) = xlt11
c	if(n.ge.24) XLT1(N) = yltibarx	


	xltmi(n) = 1./xlt1(n)
	xk = 1.6e-19 

c	vpinchi(n) = (-1.*xmtor(1)/xni(1)  - 
c	2	xmas(1)*xnuc(1,2)*torv(n) +	 xmas(1)*(xnuc(1,2)+ynudrag1(n))*
c	3  ((erada(n)/bthet)+vtheta(1)/fp))/(eq*atnum(1)*bthet)

	if(ioptpinchi.eq.1) then
c	formulation using vphi to evaluate drag term  
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc(1,2)*(1.-atnum(1)/zbar2(n))/
     2			((eq*atnum(1)*bthet)**2)
	

	vpinchi(n) = (-1.*xmtor(1)/xni(1) + 
	2	xmas(1)*ynudrag1(n)*vphia(1) + xmas(1)*xnuc(1,2)*
	3  (vtheta(1)-vtheta(2))/fp)/(eq*atnum(1)*bthet)
	endif
	if(ioptpinchi.eq.10) then
c	formulation using vphi to evaluate drag term  
 	diffA(n) =  xmas(1)*xk*ti(n)*xnuc(1,2)*(1.-atnum(1)/zbar2(n))/
     2			((eq*atnum(1)*bthet)**2)
	diffA(n) = diffA(n)*(1.+ynudrag1(n)/xnuc(1,2))

	vpinchi(n) = (-1.*xmtor(1)/xni(1)   + 
	2	xmas(1)*ynudrag1(n)*vphia(2) + xmas(1)*xnuc(1,2)*
	3  (1.+ynudrag1(n)/xnuc(1,2))*(vtheta(1)-vtheta(2))/fp)/
     4	(eq*atnum(1)*bthet)
	endif

	if(ioptpinchi.eq.0) then
c	formulation retaining vphi2 and vthet2 direct input 
	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1)	 + 
     1	 xmas(1)*ynudrag1(n)*torv(n) +
     2 xmas(1)*(xnuc12(n)+ynudrag1(n))*(velthet1(n)-velthet2(n))/fp)/
     3	 (eq*atnum(1)*bthet)
	x1 = -1.*xmomtor1(n)/yni(n,1)
	x2 = xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp) 
	x3 = xmas(1)*xnuc12(n)*(velthet1(n)-velthet2(n))/fp	 

	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*(1.-atnum(1)/zbar2(n))*
	1			((ynudrag1(n)/xnuc12(n))+1.)/
     2			((eq*atnum(1)*bthet)**2)
	endif


	if(ioptpinchi.eq.2) then
c	original formulation eliminating vphi in all terms 
	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1)	 + 
     1	 xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp) +
     2	 xmas(1)*xnuc12(n)*(velthet1(n)-velthet2(n))/fp)/
     3	 (eq*atnum(1)*bthet)
	x1 = -1.*xmomtor1(n)/yni(n,1)
	x2 = xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp) 
	x3 = xmas(1)*xnuc12(n)*(velthet1(n)-velthet2(n))/fp	 

	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*(1.-atnum(1)/zbar2(n))*
	1			((ynudrag1(n)/xnuc12(n))+1.)/
     2			((eq*atnum(1)*bthet)**2)
c	orbit squeezing	
c	param = 1.e3*eradx/(bthet*cs(1))
c	sheare =	1.-(rhot(1)/abs(fp))*dlnedr*param

	endif
c	use measured vphi-z
	if(ioptpinchi.eq.5) then
	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1)  + 
     1 xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)+velthet1(n)/fp) 
     2	 - xmas(1)*xnuc12(n)*torv(n))/
     3	 (eq*atnum(1)*bthet)
	 
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
	1			((ynudrag1(n)/xnuc12(n))+1.)/
     2			((eq*atnum(1)*bthet)**2)
	endif

c	use calculated vphi-i & vphi-z
	if(ioptpinchi.eq.6) then
	vpinchi(n) = (-1.*xmomtor1(n)/yni(n,1) + 
     1 xmas(1)*(ynudrag1(n)+xnuc12(n))*((erada(n)/bthet)+velthet1(n)/fp) 
     2	 - xmas(1)*xnuc12(n)*omegt(n,1)*rmajor)/
     3	 (eq*atnum(1)*bthet)
	 
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc12(n)*
	1			((ynudrag1(n)/xnuc12(n))+1.)/
     2			((eq*atnum(1)*bthet)**2)
	endif



c	original formulation
 	if(ioptpinchi.eq.3) then 
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc(1,2)*
     1			((xnudrag(1)/xnuc(1,2))+1.)/
     2			((eq*atnum(1)*bthet)**2)
 	vpinchi(n) = (-1.*xmtor(1)/xni(1)  + 
     1	 xmas(1)*ynudrag1(n)*(vphia(1)-ti(n)*xlpm(n)/bthet) +
     2	 xmas(1)*xnuc(1,2)*(vtheta(1)-vtheta(2))/fp)/
     3	 (eq*atnum(1)*bthet)
	endif

	if(ioptpinchi.eq.4) then
c	replace er/b + vthet/fp with vphi + )dp/dr)/neB
 	vpinchi(n) = (-1.*xmtor(1)/xni(1)  + 
     1	 xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp) +
     2	 xmas(1)*xnuc(1,2)*(vtheta(1)-vtheta(2))/fp)/
     3	 (eq*atnum(1)*bthet)
	diffA(n) =  xmas(1)*xk*ti(n)*xnuc(1,2)*
	1			((xnudrag(1)/xnuc(1,2))+1.-atnum(1)/zbar2(n))/
     2			((eq*atnum(1)*bthet)**2)
	endif



	
	tormom(n) = -1.*xmtor(1)/xni(1)/(eq*atnum(1)*bthet)
 	eph(n)	  = - eq*atnum(1)*ephia/(eq*atnum(1)*bthet)
	dragmom(n)= xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp)
     1			/(eq*atnum(1)*bthet)
	fricmom(n)=	xmas(1)*(xnuc(1,2)*(vtheta(1) -	vtheta(2))/fp)
	1			/(eq*atnum(1)*bthet)

	

	coefp(n) = vpinchi(n) + diffA(n)/xlt1(n)
	ratnu(n) = xnudrag(1)/xnuc(1,2)
250	continue
	goto 300
c	new ion & impurity density and ion & electron temp formulation
	vpinchi(n) = (-1.*xmtor(1)/xni(1) - eq*atnum(1)*ephia + 
     2			xmas(1)*ynudrag1(n)*((erada(n)/bthet)+vtheta(1)/fp) +
     3			xmas(1)*(xnuc(1,2)*(vtheta(1) -	vtheta(2))/fp))
     4			/(eq*atnum(1)*bthet)
 	vpinchz(n) = (-1.*xmtor(2)/xni(2) - eq*atnum(2)*ephia + 
     2			xmas(2)*ynudrag2(n)*((erada(n)/bthet)+vtheta(2)/fp) +
	3			xmas(2)*(xnuc(2,1)*(vtheta(2) -	vtheta(1))/fp))
	4			/(eq*atnum(2)*bthet)
	ddi = xk*(xmas(1)*ti(n)*xnuc(1,2)/((atnum(1)*eq*bthet)**2))
	dii(n) = ddi*(1.+ynudrag1(n)/xnuc(1,2))
	diz(n) = ddi*atnum(1)/atnum(2)
	ddz = xk*(xmas(2)*ti(n)*xnuc(2,1)/((atnum(2)*eq*bthet)**2))
 	
      dzz(n) = ddz*(1.+ynudrag2(n)/xnuc(2,1))
	dzi(n) = ddz*atnum(2)/atnum(1)
300	continue
c************vpinch < 0 required*******************
	if(vpinchi(n).gt.0.0)  vpinchi(n) = 0.0
	
c**************************************************
c	if(ioptran.eq.1) goto 700 
C	INFER CHI's FROM EXP N,T,GSCL & CALC. Q, GAM
	
	GAMEL = ATNUM(1)*GAMION(N,1) + ZBAR2(N)*GAMION(N,2)
	XCHIE(N) = EXLTE(N)*((GAMHEATE(N)/(EXNE(N)*XK*XTE(N)))-
     1					2.5*GAMEL/EXNE(N))	 
	XNION = EXNE(N)/(ATNUM(1)+FRACZ*ZBAR2(N))
	XCHII(N) = EXLTI(N)*((GAMHEATI(N)/(XNION*XK*XTI(N)))-
     1					2.5*GAMION(N,1)/XNION)
 


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
	XVTHI = SQRT(2.*XK*XTI(N)/XMAS(1))
	XVTHIMP = SQRT(2.*XK*XTI(N)/XMAS(2))
	EMASS = 9.1E-31
	xvthe = sqrt(2.*xk*xte(n)/emass)
	CLIGHT = 3.E8 
	CSE = SQRT(2.*XK*XTE(NMESH)/EMASS)
 


	

	XNELECTRON = XNI(1)*(ATNUM(1)**2) + XNI(2)*(zbar2(nmesh)**2) 
	XNUEI =EXNE(N)/(6.4E14*((1.E-3*XTE(N))**1.5)) 
	XNUEIAST(N) = XNUEI*ABS(QSAFE)*RMAJOR/XVTHE

	ep = aminor*SQRT((1.+ELONG**2)/2.)/rmajor 
	bfield = abs(bphi) 
	OMI =EQ*BFIELD/XMAS(1)
	CSOUND = SQRT(XK*XTE(NMESH)/XMAS(1))
 	rhot = csound/omi
C	NEOCLASSICAL CHI FOR IONS
c	simple neoclassical chi 
 
	CHINC(n) = ((RHOTi(1)*bfield/bthet)**2)*XNUc(1,2)*(EP**0.5) 
C	CHANG-HINTON CHI
	ALFA = XNi2*(zbar2(n)**2)/(XNi1*(ATNUM(1)**2))
	qa = ep*bfield/bthet
	XMUii =(XNUc(1,1)*Q95*RMAJOR/(XvthI*(EP**1.5)))*(1.+1.54*alfa)
	
	dp = 0.
100	G1 = (1. + 1.5*((EP**2)+ep*dp)+.375*(ep**3)*dp)/(1.+.5*ep*dp)
	G2 =SQRT(1.-(EP**2))*(1.+0.5*ep*dp)/(1.+(dp/ep)*(sqrt(1.-ep**2)
	1	-1))
	A1 =(0.66*(1.+1.54*ALFA)+(1.88*SQRT(EP)-1.54*EP)*(1.+3.75*ALFA))/
	1	(1.+1.03*SQRT(XMUii)+0.31*XMUii)
	A2 =0.59*XMUii*EP*(1.+1.33*ALFA*(1.+0.6*ALFA)/(1.+1.79*ALFA))/
     1	(1.+0.74*XMUii*(EP**1.5))	 
	
	betap = 2.*Xni1*xk*Xti(n)/((bthet**2)/(2.*1.257e-6)) 
	CHICH(n) = CHINC(n)*(xnuc(1,1)/xnuc(1,2))*(A1*G1+A2*(G1-G2))
	if(dp.eq.0) then
		chich0 = chich(n)
		dp = -1.*ep*(betap+0.5*log(1.65+0.89*(qa-1.)))
  		goto 100
	endif
	
	
c	orbit squeezing
	
c	if(n.eq.25) then
c	dlnEdr = log(erada(n)/erada(n-1))/delma 
c	goto 600
c	endif
	sheare(n) = 1.0
	if(n.lt.22) goto 625
	dEdr = (erada(n+1)-erada(n))/delma 
	
600	para = 1./(bthet*XvthI)
c	if(ioptshear.eq.0) sheare(n) =	1.-(rhoti(1)/abs(fp))*dlnEdr*para
	sheare(n) =	1.-(rhoti(1)/abs(fp))*dEdr*para
	if(n.eq.24) sheare(25) = sheare(24)
	if(abs(sheare(n)).lt.1.0) sheare(n) = 1.0
625	chinc(n) = chinc(n)/(abs(sheare(n))**1.5)
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
	CHIETAI(N) = 1.25*((CSOUND**2)*RHOT/OMI)*
     1	SQRT(1./(EXLTI(NMESH)*RMAJOR))
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
	
	eplasfreq = 56.4*sqrt(EXNE(NMESH))
	CHIEETG(N)=0.13*((Cse/EPLASFREQ)**2)*xvthe*SHEARM*ETAE(n)*
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
c	paleoclassical for electrons
	AM = (aminor*SQRT(0.5*(1.+ELONG**2)))
 
	ep = am*rhor(n)/rmajor 
 
	EQ = 1.6E-19
 	XK = 1.6E-19
	EP0 = 8.854E-12
	xme = 9.1e-31 
	Yz = SQRT(eXNe(n)) 
 	X = (EP0/EQ)**1.5	
	COULOGe = LOG(12.*3.1416*(xte(n)**1.5)*X/Yz)
c	ome = 56.4*sqrt(exne(n))
	ome = 1.76e11*bfield 
 
	xlame = 1.2e16*(xte(n)**2)*17./(exne(n)*zeff*couloge)
	xnuestar = (rmajor*q95)/((ep**1.5)*xlame)
	fc = (((1.-ep)**2)/sqrt(1.-(ep**2)))/(1.+1.46*sqrt(ep)+0.2*ep)
	xmunu = ((zeff+1.414 - log(2.414))/
	1		(zeff*(1 +sqrt(xnuestar)+xnuestar)))*(1.-fc)/fc
	const=(1.414+zeff)/(1.414+(13./4.)*zeff) + 	xmunu  
      		
      bigL = xlame
	vthe = sqrt(2.*xk*xte(n)/xme)
	xnue = (vthe/xlame)
	dele = (vthe/ome)
	xlmax = rmajor*sqrt(3.1416*q95)/sqrt((vthe/ome)*ssi95/aminor)
	if(xlmax.lt.xlame) bigL = xlmax 
	chiepale(n)=1.5*(1. + bigL/(3.1416*q95*rmajor))*(const + xmunu)
	1			*xnue*(dele**2)
700	continue 
	atnum(2) = atnum2
756	continue 
	RETURN	 
	END

	
		
