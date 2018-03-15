	SUBROUTINE POLOIDAL(nmesh)
	INCLUDE 'soldiv.FI'
C  SOLVES THE 6 COUPLED EQUATIONS FOR ION & IMPURITY NORMALIZED
C  POLOIDAL VELOCITIES AND SIN & COS COMPONENTS OF DENSITY ASYM.
	PARAMETER (Jq=4,jk=2,jv=201,nm=26) 
	DIMENSION A(Jq), B(Jq), S(Jq),C(Jq,Jq),gam(jq),vrhat(jk),
	2	vthetsin(jk),vthetcos(jk),hnuat(jk),hnuion(jk),sc1(jk),
	3	sc2(jk),cc1(jk),cc2(jk),vc1(jk),vc2(jk),vc3(jk),vc4(jk),
	4	vs1(jk),vs2(jk),vs3(jk),vs4(jk),dela(jk),dels(jk),
	5	XV(JK),YV(JK),FRIC(JK),grad(jk), VPHIHAT(JK),
	6	vt(jv),ss(jk),xmomhat(jk),hnucol(jk),XNUATSTAR(JK),
	7	COSPHI(JK),SINPHI(JK), hnuatom(jk),hnuatomstar(jk),
	8	vrhatT(jk),presshat(jk),aa(jk),chk(jk),sq(jk),shs(jk),
	9	ahs(jk),zk(jk),ym00(jk),erexhat(nm)

298	format(1x,'OK')
301	format(7e11.4)

c************SOL potential option****************************************
c		ioptpot = 1 local phi asymm determined by electron asymm (Max-Boltz)
c				= 2 local phi asym determined by SOL asym
	if(ioptpot.ne.1) ioptpot = 2
c************************************************************************
c	ipolopt = 1 uses vthet2 = vexp, =0 calculates vthet2
c************************************************************************
c	ipolopt = 1
		n = nmesh
c	initialize asymmetries to zero
	vtheta(1) = 0.
	vtheta(2) = 0. 
	cosden(1)=0.
	cosden(2)=0.
	sinden(1)=0.
	sinden(2)=0.
	vphicos(n,1)=0.
	vphicos(n,2)=0.
	vphisin(n,1)=0.
	vphisin(n,2)=0.
	vphic1=0.
	vphics=0.
	vhpis1=0.
	vphis2=0.
	vthetcos(1)=0.
	vthetcos(2) = 0.
	vthetsin(1)=0.
	vthetsin(2)=0

	xk = 1.6e-19 
	nit = 0
	
	iopteq = 2
	C9 = 0.0
	c8 = 0.0
	c7 = 0.0
	epol = 1.e-3
c	****set ERAD TO EXP & ephi/Te to calculated value FROM EXP****
	RADE = EREX(N)
	phite =	 epot(n)
c	if(abs(phite).lt.0.3) phite = 0.3*phite/abs(phite)
	estatpot(n) = epot(n)*xte(n)
c	**************************************************************
	xne = xni(1)*atnum(1)+xni(2)*zbar2(n)  
	GAM(1) = XNi(1)*ATNUM(1)/(XNi(1)*ATNUM(1)+XNi(2)*zbar2(n))
 	GAM(2) = XNi(2)*zbar2(n)/(XNi(1)*ATNUM(1)+XNi(2)*zbar2(n))
	if(ioptedge.ne.0) xne = xni(1)*atnum(1) + xni(2)*zbar2(n)	 
		
	SIGNB = fp/abs(fp)
 	 
C	POLOIDAL ROTATION EQS
			
c	couple to toroidal rotation
c	if(it.gt.4) then 
c	vphia(1) = rmajor*omegt(n,1)	
c	vphia(2) = rmajor*omegt(n,4) 
c	endif
c	5/03 formulation
	
c	goto 1003
	if(ioptedge.eq.1) 	then
  	bfield = bphi
	x1 = atnum(1)*eq*bfield/(xmas(1)*VTH(1))
	y1 = qsafe*rmajor/vth(1)
	z1 = xmompol(n,1)/(yni(n,1)*xmas(1)*vth(1))
	vr1 = gamion(n,1)/yni(nmesh,1)
c********************change************************
c	vrad1(nmesh) = vr1
c	if(vr1.gt.5.0) vrad1(nmesh) = 5.0
c	if(vr1.gt.10.) vrad1(nmesh) = 10.
	vrhat(1) = x1*y1*vr1
c	*********temp*****************
c	vrhat(1) = 0.0

	xmomhat(1) = z1*y1
	x2 = atnum(2)*eq*bfield/(xmas(2)*VTH(2))
 	y2 = qsafe*rmajor/vth(2)
	z2 = xmompol(n,2)/(yni(n,2)*xmas(2)*vth(2))
	

	R2 = 0.01
c	sputtering yield for C in the 100-1000eV range is about .02,
c	allow for half othis to be trapped in divertor , R2 = 0.01
	vr2 = gamion(nmesh,2)/(yni(nmesh,2))

	vrad2(nmesh) = vr2
c	vrad2(nmesh) = vrad1(nmesh)

	vrhat(2) = x2*y2*vr2
	xmomhat(2) = z2*y2

	
 	endif


C	******calculate temp grad source term & add to vrhat********************************************
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

 	xnustar(1,1) = xnu11star
	xnustar(1,2) = xnu12star
	xnustar(2,1) = xnu21star
	xnustar(2,2) = xnu22star


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
c	********************HEAT FLUX source term*******************************************************************

	vrhatT(1)=vrhat(1)+ qedge(n)*f1*xki*xti(n)/(exlti(n)*bphi*vth(1)) 
	vrhatT(2)=vrhat(2)+ qedge(n)*f2*xkz*xti(n)/(exlti(n)*bphi*vth(2)) 
	vphihat(1) = vphia(1)/vth(1)
 	vphihat(2) = vphia(2)/vth(2)
	RIE = TE(N)/TI(N)
c	****************************************************************************************************	 
	DO 100 NN=1,100
 	epd =0.5

	

	vt(nn) =0.0 
 	VTHET1 = VTHETa(1)
	VTHET2 = VTHETa(2)

c	************** 
	zc =1.
c  ****************change**********************
		if(ioptedge.eq.1) then

		hnuat(1) = rminor*coldno(n)*svata(n)/VTH(1)
		avion = xnuioni(n)+ xnuionb(n)
	1							
		hnuion(1) =rminor*avion/VTH(1)
	
	
		hnuat(2) = 0.
		hnuion(2)=0.0

	

	
	hnuatom(1) = hnuat(1) + hnuion(1)
	hnuatom(2) = hnuat(2) + hnuion(2)

C
c	HNUATOM(1) = 0.
c	HNUION(1) = 0.
c	HNUAT(1) = 0.
 

	XNUATSTAR(1) = HNUAT(1)*QSAFE*RMAJOR/RMINOR
	XNUATSTAR(2) = HNUAT(2)*QSAFE*RMAJOR/RMINOR
	XNUATomSTAR(1) = HNUATom(1)*QSAFE*RMAJOR/RMINOR
	XNUATomSTAR(2) = HNUATom(2)*QSAFE*RMAJOR/RMINOR
	
	  endif 

     	 	
c	*******************************************
	

c	**************	
1003	nlimit = 2
	if(ioptedge.eq.3) nlimit = 1
c	Rutherford friction factor C1		
	cf = 0.35 
	cf = 1.0
1004	DO 1005 J=1,2
	 	K=2
	IF(J.EQ.2) K=1
	zcc = zc



      IF(IOPTPOL.EQ.1) GO TO 666	
c	***change sign of f(j) terms  7/2/07******************

	A(J)= 1.*FP*(QSAFE**2)*F(J)*(1.0+(zcc*(2./3.)*cosden(j)+COSPOT+
	1		c7*dela(j))) 
     2   +	cf*FP*XNUSTAR(J,K)*(1.+c9*fric(j)) - 
     3	QSAFE*ep*(VPHIa(J)/VTH(J))*(SINDEN(J) + SINPOT)
      B(J)=-1.*cf*FP*XNUSTAR(J,K)*SQRT(XMAS(J)/
     1	XMAS(K))*(1.+c9*fric(k))
	s(j) = xmomhat(j) - vrhatT(j) 
  	S(J) = s(j) - QSAFE*ep*(((VPHIa(J)/VTH(J))**2)*SINPOT + 
	1	(VPHIa(J)/VTH(J))*((PRESS(J)/VTH(J))*SINPOT + 
     2	0.5*(VPHIa(J)/VTH(J))*SINDEN(J))) 
     3+   (QSAFE**2)*F(J)*FP*(((VPHIa(J)+PRESS(J))/VTH(J))*
     4	COSPOT+c7*dels(j))
	SPOLA(N,J) = S(J)  
	S(J) = S(J) - 0.25*QSAFE*EP*sinden(j) 
	SPOLB(N,J) = S(J)
	zz = atnum(j)
	if(j.eq.2) zz = zbar2(n)
	s(j) = s(j)	- 0.25*QSAFE*EP*zz*(teP/temp(j))*
     2	((1.+cosden(j))*sinpot - sinden(j)*cospot)
	nit = nit + 1
	
c	if(ioptedge.eq.1) ioptedge =100 
	IF (IOPTEDGE.EQ.1) THEN
	  a(j) = a(j) + qsafe*fp*hnuat(j)/ep
	  xy = (xne/xni(j))*(cosdene+cosdeno(j))-(cosden(j)+cosdeno(j)) 
       xz = (xne/xni(j))*(sindene+sindeno(j))-(sinden(j)+sindeno(j)) 
        yb = (1.+cosden(j))*xy+sinden(j)*xz
	  a(j) = a(j) + 0.5*fp*ep*hnuion(j)*yb
	  yx = (vphia(j)/VTH(j))*ep*(cosdene+cosdeno(j))
	  yz = 0.667*qsafe*f(j)*(sindene+sindeno(j))
	  s(j) = s(j) - (xne/xni(j))*hnuion(j)*qsafe*(yx-yz)
	endif
c	GOTO 1005
C	NEW FORMULATION 8/31/07 calc VPH_S,C from rad mom bal***********, 
c	******************************************
666	ECH = atnum(1)
c************SOL potential option****************************************
c		ioptpot = 1 local phi asymm determined by electron asymm (Max-Boltz)
c				= 2 local phi asym determined by SOL asym
c************************************************************************
	 
	presshat(j) = press(j)/vth(j)
	IF(J.EQ.2) ECH = ZBAR2(N)
	if(ioptpot.eq.2) phite = epot(n)/xte(n)
		PHIT = ECH*RIE*PHITE
			
	A(J) =0.5*FP*EP*QSAFE*hnuion(j)*(xne/xni(j))*((1. + cosden(j))*
     1		(cosdene+cosdeno(j))+sinden(j)*(sindene+sindeno(j)))
     2		+(1.+(2./3.)*cosden(j))*(QSAFE**2)*FP*f(j)
     3	+ XNUSTAR(J,K)*FP + XNUATomSTAR(J)*FP
	4	-qsafe*vphihat(j)*ep*sinden(j)
	aa(j) = a(j)
      B(J) = -1.0*XNUSTAR(J,K)*FP*SQRT(XMAS(J)/XMAS(K))

	S(J) = -1.0*VRHATT(j) 
	1		- 0.25*EP*QSAFE*sinden(j)     
	
	
	S(J) = S(J) - VPHIHAT(J)*EP*QSAFE*(0.5*VPHIHAT(J)*sinden(j))
      S(J) = S(J)	- (xne/xni(j))*hnuion(j)*
	1	((-2./3.)*(qsafe**2)*f(j)*(sindene+sindeno(j)) +
	2	 qsafe*ep*vphihat(j)*(cosdene+cosdeno(j))/fp)
	

	if(ioptpot.eq.1) then
c	pre-2/09 version
	S(J) = S(J) + (QSAFE**2)*F(J)*FP*(erex(n)/(vth(j)*bthet))*cospot
	s(j) = s(j) + 0.25*ep*qsafe*ech*(Te(n)/ti(n))*
     1		 (sinden(j)*cosdene - (1.+cosden(j))*sindene)
 	s(j) = s(j) - (xlphim(n)*te(n)/(bthet*vth(j)))*
	1 ((qsafe**2)*f(j)*fp*cosdene + qsafe*ep*vphihat(j)*sindene)
	endif
	if(ioppot.eq.2) then
c	SOL potential effect
	cospot = epotcos
	sinpot = epotsin
	a(j) = a(j) + (2.+(2./3.)*cosden(j))*(qsafe**2)*fp*f(j)
	s (j) = s(j) + 0.25*ep*qsafe*ech*(Te(n)/ti(n))*
     1		 (sinden(j)*cospot - (1.+cosden(j))*sinpot)	     
 	s(j) = s(j) + (erexhat(n)/(bthet*vth(j)))*
	1 ((qsafe**2)*f(j)*fp*cospot - qsafe*ep*vphihat(j)*sinpot)
	s(j) = s(j) + 4.0*(qsafe**2)*fp*vphihat(j)
	endif
 





	sq(j) = s (j)
	  if(j.eq.2) then
	ss20 = -1.0*VRHATT(j)
	ss21 = 	- 0.25*EP*QSAFE*sinden(j) 
	ss22 = 	0.25*ep*qsafe*ech*(Te(n)/ti(n))*
     1		 (sinden(j)*cosdene - (1.+cosden(j))*sindene)
	ss23 = -1.0*(QSAFE**2)*F(J)*FP*(erex(n)/(vth(j)*bthet))
	ss24 = -1.0*(VPHIHAT(J)**2)*EP*QSAFE*0.5*sinden(j)
	ss25 = -1.0*(xne/xni(j))*hnuion(j)*
     1	((-2./3.)*(qsafe**2)*f(j)*(sindene+sindeno(j)) +
	2	 qsafe*ep*vphihat(j)*(cosdene+cosdeno(j))/fp)
 
	ss26 = -1.0*(xlphim(n)*te(n)/(bthet*vth(j)))*
     1 ((qsafe**2)*f(j)*fp*cosdene + qsafe*ep*vphihat(j)*sindene)
	  endif
c	****source and coef for hirshman-sigmar********************
	ym00(1) = ym00i
	ym00(2) = ym00z
	zk(1) = xki
	zk(2) = xkz 
	shs(j) = s(j) + (QSAFE**2)*F(J)*FP*(erex(n)/(vth(j)*bthet))
	shs(j) = shs(j) + vrhatt(j) -
	1	vrhat(j)+qedge(n)*zk(j)*xki*xti(n)/(exlti(n)*bphi*vth(j))
	ahs(j) = a(j) - (1.+(2./3.)*cosden(j))*(QSAFE**2)*FP*f(j)
	xvisc = ym00i*xnustar(j,j)
	if(j.eq.2) xvisc = ym00z*xnustar(j,j)
	ahs(j) = ahs(j) + (qsafe**2)*fp*xvisc
c	************************************************************
c	********Kim-Diamond-Groebner********************************

c	write(6,301) a(j),b(j),s(j),f(j),xnustar(j,k),cs(1),cs(2)

1005	CONTINUE 
	APol1(n) = A(1)
 	Apol2(n) = A(2)
	Bpol1(n) = B(1)
	Bpol2(n) = B(2)
	Spol1(n) = S(1)
	Spol2(n) = S(2) 

1006	continue

c	*****************************trial vhs ***********************

	vpolhs(n) = (shs(1)-b(1)*vthexp(n))/ahs(1)
	
c	**************************************************************      
	  	
	
      VTHETA(1) = (A(2)*S(1)-B(1)*S(2))/(A(2)*A(1) - B(2)*B(1))
	if(ioptedge.ne.3)
     1	VTHETA(2) = (A(1)*S(2)-B(2)*S(1))/(A(1)*A(2) - B(1)*B(2))
c	if(abs(vtheta(1)).gt.1.) goto 100
c	if(abs(vtheta(2)).gt.1.) goto 100
	vth(1) = sqrt(2.*xk*ti(n)/xmas(1))
	vth(2) = sqrt(2.*xk*ti(n)/xmas(2))
	vt1 = vtheta(1)*fp*vth(1)
	vt2 = vtheta(2)*fp*vth(2)
	if(ipolopt.eq.1) then
	vtheta(2) = vthexp(n)/(fp*vth(2))
	vtheta(1) = (s(1)-b(1)*vtheta(2))/a(1)
c	goto 150
	ENDIF
c	goto 150
c	**********solution w/o poloidal asymmetries************	
	   if(nn.eq.1) then
	vpcal1(n) = vt1
	vpcal2(n) = vt2
	   endif
c	*********************************************************

c	*********Hirshman-Sigmar*********************************
		if(nn.eq.1) then
	VTHS1(n) = (Ahs(2)*Shs(1)-B(1)*Shs(2))/
	1			(Ahs(2)*Ahs(1) - B(2)*B(1))
	VTHS2(n) = (Ahs(1)*Shs(2)-B(2)*Shs(1))/
	1			(Ahs(1)*Ahs(2) - B(1)*B(2))
    	vths1(n) = vths1(n)*fp*vth(1)
	vths2(n) = vths2(n)*fp*vth(2)
	vpolhs(n)= (Shs(1)-b(1)*vthexp(n)/(fp*vth(1)))/Ahs(1)
	vpolhs(n)= vpolhs(n)*fp*vth(1)
		endif
	goto 150
c	*********Kim-Diamond-Groebner*****************************
	vtkdg1(n) = -1.0*xki*xti(n)/(exlti(n)*atnum(1)*bphi)
	vtkdg2(n) = -1.0*xti(n)*(((xki+1.5*xkz)/exlti(n)) -
	1			xlpm(n)*(1. - atnum(1)/zbar2(n)))/(atnum(1)*bphi)
c	**************************************************************
1015	CONTINUE

c	if(nmesh.eq.25.and.nn.gt.50) epd = 0.1
	IVTHET = IVTHETREAL
	if(nn.gt.1) then
 	VTHETa(1) = (1.-epd)*VTHETa(1)+epd*VTHET1
	VTHETa(2) = (1.-epd)*VTHETa(2)+epd*VTHET2
	endif
c	***************************************************
	xx = vtheta(1)*vth(1)*fp
c	if(abs(xx).gt.1.e4) vtheta(1) = 1.e4*vtheta(1)/abs(xx)
	yy = vtheta(2)*vth(2)*fp
c	if(abs(xx).gt.1.e4) vtheta(2) = 1.e4*vtheta(2)/abs(xx)
c	****************************************************
	
c	****new feature****iterate the density & velocity asymmetries to convergence
c	on each iterate on the poloidal velocities*********
c	do 90 jj = 1, 50 
	PS1 = VPHIS1
	PS2 = VPHIS2
	PC1 = VPHIC1
	PC2 = VPHIC2
	SPOTS = SINPOT
	CPOTC = COSPOT
C	DENSITY ASYMMETRIES	
	IF(IOPTPOL.EQ.1) GOTO 570
C	********ORIGINAL FORMULATION***********************************	 
c	cosine moment species 1		
	C(1,1)=-1.*F(1)*(QSAFE**2)*FP*VTHETa(1)*((1./3.)-c7*vc1(1))+c7* 	
	1vc2(1)-0.5*(EP**2)*cf*FP*XNUSTAR(1,2)*VTHETa(2)*SQRT(XMAS(1)
     2	/XMAS(2))
	C(1,2) = -0.5*QSAFE*EP*((FP**2)*(VTHETa(1)**2)-0.5) +
	2  0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(1) + c7*vc3(1)
	C(1,3) = 0.5*FP*(EP**2)*cf*XNUSTAR(1,2)*VTHETa(1)
	C(1,4) = 0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(2)
	S(1) = 0.5*F(1)*(QSAFE**2)*FP*(VTHETa(1)*(1.+COSPOT+c8*cc1(1))-
	1	((VPHIa(1)+PRESS(1))/VTH(1))*COSPOT+c8*cc2(1))	
	s(1) = s(1) - QSAFE*EP*xv(1)*c9	+ c7*vc4(1)
	if(ioptedge.eq.1) then
	  c(1,1) = c(1,1)+0.5*ep*qsafe*hnuion(1)*fp*vtheta(1)
	  c(1,2) = c(1,2)-0.5*ep*qsafe*hnuat(1)*hnuion(1)
	  xy = fp*vtheta(1)*cosdeno(1)-hnuion(1)*sindeno(1)
	  xz = cosdeno(1)*(1.+xne/xni(1))+(xne/xni(1))*cosdene
	  yz = (xne/xni(1))*(sindene+sindeno(1))
	  s(1) = s(1)-0.5*ep*qsafe*hnuat(1)*xy-qsafe*hnuion(1)*fp*
     1         (0.5*vtheta(1)*xz +0.333*(qsafe/ep)*f(1)*yz)
      endif	   
c	sine moment species 1
	C(2,1) = 0.5*QSAFE*EP*((FP**2)*(VTHETa(1)**2)-0.5) -
	1	0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(1) +c7*vs3(1)	
	C(2,2) = -1.*(QSAFE**2)*F(1)*FP*VTHETa(1)*((1./3.)+vs1(1))+vs2(1) 
	1  -	0.5*(EP**2)*cf*FP*XNUSTAR(1,2)*VTHETa(2)*SQRT(XMAS(1)/XMAS(2))
	C(2,4) = 0.5*(EP**2)*cf*FP*XNUSTAR(1,2)*VTHETa(1)
 	C(2,3) = -0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(2)
	S(2) = 0.5*(QSAFE**2)*F(1)*FP*(VTHETa(1)*(sinpot+c8*sc1(1))-
	1	((VPHIa(1)+PRESS(1))/VTH(1))*SINPOT+c8*sc2(1))
	S(2) = S(2) - 0.5*QSAFE*EP*(((VPHIa(1)/VTH(1))**2)+c9*yv(1))
	2  + 0.5*QSAFE*(FP**2)*EP*(VTHETa(1)**2) + c7*vs4(1)
	if(ioptedge.eq.1) then
	  c(2,1) = c(2,1) + 0.5*ep*qsafe*hnuion(1)*hnuat(1)
	  c(2,2) = c(2,2) + 0.5*ep*qsafe*hnuion(1)*fp*vtheta(1)
	 yz = (xne/xni(1))*(cosdene+cosdeno(1))
	  xy = fp*vtheta(1)*sindeno(1)+hnuion(1)*cosdeno(1)
	  xz = sindeno(1)*(1.+xne/xni(1))+(xne/xni(1))*sindene
	 								   
	  s(2) = s(2)-0.5*ep*qsafe*hnuat(1)*xy-qsafe*hnuion(1)*fp*
     1         (0.5*vtheta(1)*xz +0.333*(qsafe/ep)*f(1)*yz)
	endif
c	cosine moment species 2
	C(3,3) = -1.*F(2)*FP*(QSAFE**2)*VTHETa(2)*((1./3.)-vc1(2))+c7*
	1  vc2(2) - 0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(1)*SQRT(XMAS(2)
	2	/XMAS(1))
	C(3,4) = -0.5*QSAFE*EP*((FP**2)*(VTHETa(2)**2)-0.5) +
	2	0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(2) +c7*vc3(2)			
	C(3,1) = 0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(2)
	C(3,2) = 0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(1)
	S(3) = 0.5*F(2)*(QSAFE**2)*FP*(VTHETa(2)*(1.+COSPOT+c8*cc1(2))-
	1	((VPHIa(2)+PRESS(2))/VTH(2))*COSPOT+c8*cc2(2))	
	s(3) = s(3) - QSAFE*EP*xv(2)*c9	 +c7*vc4(2)
	if(ioptedge.eq.1) then
	  c(3,3) = c(3,3)+0.5*ep*qsafe*hnuion(2)*fp*vtheta(2)
	  c(3,4) = c(3,4)-0.5*ep*qsafe*hnuat(2)*hnuion(2)
	  xy = fp*vtheta(2)*cosdeno(2)-hnuion(2)*sindeno(2)
	  xz = cosdeno(2)*(1.+xne/xni(2))+(xne/xni(2))*cosdene
	  yz = (xne/xni(2))*(sindene+sindeno(2))
	  s(3) = s(3)-0.5*ep*qsafe*hnuat(2)*xy-qsafe*hnuion(2)*fp*
     1         (0.5*vtheta(2)*xz +0.333*(qsafe/ep)*f(2)*yz)
      endif	   

c	sine moment species 2
	C(4,3) = 0.5*QSAFE*EP*((FP**2)*(VTHETa(2)**2)-0.5) -
	1	0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(2) +c7*vs3(2)
	C(4,4) = -1.*(QSAFE**2)*F(2)*FP*VTHETa(2)*((1./3.)+vs1(2))+c7*
 	1  vs2(2) - 0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(1)*SQRT(XMAS(2)/
	2	XMAS(1))
	C(4,2) = 0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(2)
 	C(4,1) = -0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(1)
	S(4) = 0.5*(QSAFE**2)*F(2)*FP*(VTHETa(2)*(sinpot+c8*sc1(2))-
	1	(((VPHIa(2)+PRESS(2))/VTH(2))*SINPOT+c8*sc2(2))) 
	S(4) = S(4) -0.5*QSAFE*EP*(((VPHIa(2)/vth(2))**2)+C9*YV(2)) - 
     2	0.5*QSAFE*(FP**2)*EP*(VTHETa(2)**2) + c7*vs4(2)
	if(ioptedge.eq.1) then
	  c(4,3) = c(4,3) + 0.5*ep*qsafe*hnuion(2)*hnuat(2)
	  c(2,2) = c(2,2) + 0.5*ep*qsafe*hnuion(2)*fp*vtheta(2)
	  xy = fp*vtheta(2)*sindeno(2)+hnuion(2)*cosdeno(2)
	  xz = sindeno(2)*(1.+xne/xni(2))+(xne/xni(2))*sindene
	  yz = (xne/xni(2))*(cosdene+cosdeno(2))
	  s(4) = s(4)-0.5*ep*qsafe*hnuat(2)*xy-qsafe*hnuion(2)*fp*
     1         (0.5*vtheta(2)*xz +0.333*(qsafe/ep)*f(2)*yz)
 	endif
 
	xx = (0.5*qsafe*ep)

	do 3098 il = 1,4
	s(il) = s(il)/xx
	do 3098	jm = 1,4
	c(il,jm) = c(il,jm)/xx
3098	continue 
	GOTO 575

570	CONTINUE
C	********NEW FORMULATION 8/31/07*******VPHI_S,C,VTHET_S,C from rad mom bal*********
	VPHIS1 = VPHISIN(1,N)
	VPHIS2 = VPHISIN(2,N)
	VPHIC1 = VPHICOS(1,N)
	VPHIC2 = VPHICOS(2,N) 
	PHIT1 = RIE*PHITE
	PHIT2 =ZBAR2(N)*RIE*PHITE
	erhat1 = erex(n)/(vth(1)*bthet)
	erhat2 = erex(n)/(vth(2)*bthet)

	XNE = ATNUM(1)*XNI(1) + ZBAR2(N)*XNI(2)
C	*******SIN MOM, SPECIES 1***************************************	 

	C(1,1) = -0.25*QSAFE + 0.5*QSAFE*(FP**2)*(VTHETA(1)**2)
	1	-0.25*qsafe*atnum(1)*gam(1)*te(n)/ti(n)
     2	-0.333*(qsafe**2)*(f(1)/ep)*HNUION(1)*XNE/XNI(1) 
	C(1,2) =0.5*EP*VRHAT(1) - 0.5*EP*XNUSTAR(1,2)*FP*SQRT(XMAS(1)/
     1       	XMAS(2))*VTHETA(2) + 0.5*EP*fp*XNUATomSTAR(1)*VTHETA(1) -
     2	((QSAFE**2)*F(1)*FP*VTHETA(1)/EP)* 
	3	((vtheta(1)/3.)+(presshat(1)/2.)-(erhat1/2.)-
     4	(xlphim(n)*te(n)/(2.*vth(1)*bthet)))
     5	+0.5*qsafe*fp*vtheta(1)*HNUION(1)*XNE/XNI(1)
	C(1,3)=(-1./3.)*(QSAFE**2)*(F(1)/EP)*HNUION(1)*XNE/XNI(1)*gam(2) -
	2		0.25*qsafe*(te(n)/ti(n))*atnum(1)*gam(2)
	C(1,4) =  0.5*EP*XNUSTAR(1,2)*FP*VTHETA(1) +
	1		0.5*QSAFE*FP*HNUION(1)*vtheta(1)*XNE/XNI(1)*gam(2)	-
	2	0.5*(qsafe**2)*(f(1)*fp/ep)*xlphim(n)*te(n)*zbar2(n)/
     3	(vth(1)*bthet) 
	S(1) =  - 0.5*QSAFE*(VPHIHAT(1)**2) 
     1	- 0.5*qsafe*(FP**2)*(VTHETA(1)**2)
     3	+ (qsafe*HNUION(1)*XNE/XNI(1))*
	4	(-0.5*(FP**2)*vtheta(1)*SINDENO(1) + 
     5	0.333*(QSAFE*F(1)/EP)*COSDENO(1))
C	******SIN MOM, SPECIES 2*****************************************     	 
	C(3,3) = -0.25*QSAFE + 0.5*QSAFE*(FP**2)*(VTHETA(2)**2)
	1 		 -0.25*qsafe*zbar2(n)*gam(2)*te(n)/ti(n)
	2		-0.333*(qsafe**2)*(f(2)/ep)*HNUION(2)*XNE/XNI(2)*gam(2)
 	C(3,4) =0.5*EP*VRHAT(2) - 0.5*EP*XNUSTAR(2,1)*FP*SQRT(XMAS(2)/
     1       	XMAS(1))*VTHETA(1) + 0.5*EP*fp*XNUATomSTAR(2)*VTHETA(2) +
     2		(QSAFE**2)*(F(2)*FP/EP)*((VTHETA(2)/3.)+(presshat(2)/2.) -
	3		(erhat2/2.) - (xlphim(n)*te(n)/(2.*vth(2)*bthet)))
	4		+ 0.5*qsafe*fp*vtheta(2)*HNUION(2)*XNE/XNI(2)*gam(2)
	C(3,1)=(-1./3.)*(QSAFE**2)*(F(2)/EP)*HNUION(2)*XNE/XNI(2)*gam(1) -
	2		0.25*qsafe*(te(n)/ti(n))*zbar2(n)*gam(1)
	C(3,2) =  0.5*EP*XNUSTAR(2,1)*FP*VTHETA(2) +
	1			0.5*QSAFE*FP*HNUION(2)*vtheta(2)*XNE/XNI(2)*gam(1)	-
	2	0.5*(qsafe**2)*f(2)*(fp/ep)*xlphim(n)*te(n)*atnum(1)/
     3	(vth(2)*bthet)
	 
	S(3) =   - 0.5*QSAFE*(VPHIHAT(2)**2) 
     1	- 0.5*qsafe*(FP**2)*(VTHETA(2)**2)
     3	+ (QSAFE*HNUION(2)*XNE/XNI(2))*
     4	(-0.5*(FP**2)*SINDENO(2)*vtheta(2) + 
     5	0.333*(QSAFE*F(2)/EP)*COSDENO(2))
C	********COS MOM, SPECIES 1 ***************************************
	C(2,1) = C(1,2)
	C(2,2) = -1.0*C(1,1) 
	C(2,3) = C(1,4)	 
	C(2,4) = -1.0*C(1,3)

	S(2)=-0.5*(QSAFE**2)*(F(1)*FP/EP)*(vtheta(1)+presshat(1)-erhat1)  
     4	-(QSAFE*HNUION(1)*XNE/XNI(1))*
	5	(0.5*(FP**2)*vtheta(1)*COSDENO(1)+
     6	0.33*(QSAFE*F(1)/EP)*SINDENO(1))
	7	-0.25*qsafe*(ep**2)*vphisin(1,n)*vphicos(1,n)
C	***********COS MOM SPECIES 2	**************************************
	C(4,3) =  C(3,4)
 	C(4,4) = -1.0*C(3,3) 
	C(4,1) =   C(3,2)
	C(4,2) = -1.0*C(3,1)

	S(4)= -.5*(QSAFE**2)*(F(2)*FP/EP)*(vtheta(2)+presshat(2)-erhat2)  
     4	-(QSAFE*HNUION(2)*XNE/XNI(2))*
     5	(0.5*(FP**2)*vtheta(2)*COSDENO(2)+
     6	0.33*(QSAFE*F(2)/EP)*SINDENO(2))
	7	-0.25*qsafe*(ep**2)*vphisin(2,n)*vphicos(2,n)

 
	xlphimm = xlphim(n)*te(n)/(vth(1)*bthet)
575	COS1 = COSDEN(1)
	COS2 = COSDEN(2)
	SIN1 = SINDEN(1)
	SIN2 = SINDEN(2)		
	CALL LSLRG(4,C,4,S,1,A)
c	limit poloidal variation to epsilon
c	if(ABS(a(1)).ge.1.0) THEN
c		A(1) = 1.0
c		IF(A(1).LT.0.0) a(1) = -1.0
c	ENDIF
c	if(ABS(a(2)).ge.1.0) THEN
c		A(2) = 1.0
c		IF(A(2).LT.0.0) a(2) = -1.0
c	ENDIF
c	if(ABS(a(3)).ge.1.0/ep) THEN
c		A(3) = 1.0/ep
c		IF(A(3).LT.0.0) a(3) = -1.0/ep
c	ENDIF
c	if(ABS(a(4)).ge.1.0) THEN
c		A(4) = 1.0
c		IF(A(4).LT.0.0) a(4) = -1.0
c	ENDIF


	COSDEN(1) = A(1)
	SINDEN(1) = A(2)
	COSDEN(2) = A(3)
	SINDEN(2) = A(4)
c	zg = cosden(2)
c	if(zg.ge.1.0) zg = 1.0 
c	cospot = (gam(1)*cosden(1)+gam(2)*zg)/phite
	
	 
	
25	continue	 
	COSDEN(1) = (1.-epd)*COSDEN(1)+epd*COS1
	COSDEN(2) = (1.-epd)*COSDEN(2)+epd*COS2
	SINDEN(1) = (1.-epd)*SINDEN(1)+epd*SIN1
	SINDEN(2) = (1.-epd)*SINDEN(2)+epd*SIN2
	cospot = (gam(1)*cosden(1)+gam(2)*cosden(2))/phite
	sinpot = (gam(1)*sinden(1)+gam(2)*sinden(2))/phite
 
c	electron sin & cos
	cosdene =(xni(1)*atnum(1)*cosden(1)+xni(2)*atnum(2)*cosden(2))/xne 
	sindene =(xni(1)*atnum(1)*sinden(1)+xni(2)*atnum(2)*sinden(2))/xne 

c	vtheta asymmetries
	vthetcos(1) = -1.*(1. + cosden(1))
  	vthetcos(2) = -1.*(1. + cosden(2))
 	vthetsin(1) = -1.*sinden(1)
      vthetsin(2) = -1.*sinden(2)

c	****will explicitly substitute in for vethcos vthetsin to avoid iteration difficulty
	if(ioptedge.eq.1) then
	  vthetcos(1) = vthetcos(1) - (xne/xni(1))*hnuion(1)*
     1			  (sindene+sindeno(1))
	  vthetcos(2) = vthetcos(2) - (xne/xni(2))*hnuion(2)*
     1			  (sindene+sindeno(2))
	  vthetsin(1) = vthetsin(1) + (xne/xni(1))*hnuion(1)*
     1			  (cosdene+cosdeno(1))
	  vthetsin(2) = vthetsin(2) + (xne/xni(2))*hnuion(2)*
     1			  (cosdene+cosdeno(2))
	endif

	

c	IF(IOPTPOL.NE.0) GO TO 396

C	CONSTRUCT THE TOROIDAL VELOCITY ASYMMETRY FROM RADIAL MOM BAL, CONSTRUCT
C	THE POTENTIAL ASYMMETRY FROM MAXWELL-BOLTZMAN = ELECTRON POL MOM BAL

C		POTENTIAL POLOIDAL VARIATION
	cospot = (gam(1)*cosden(1)+gam(2)*cosden(2))/phite
	sinpot = (gam(1)*sinden(1)+gam(2)*sinden(2))/phite
c	if(abs(cospot).gt.0.9) cospot = 0.9*cospot/abs(cospot)
c	if(abs(sinpot).gt.0.9) sinpot = 0.9*sinpot/abs(sinpot)
	spot(n) = sinpot
	cpot(n) = cospot
	 
C		velocity asymmetry terms(product [vphi/vth]*[vphicos/ep], etc.)

c	vphicos(1,n)=xlphim(n)*te(n)*cosdene/(vth(1)*bthet) - presshat(1) 
c	1	-1.*(1.+cosden(1))*vtheta(1)
c	2	-hnuion(1)*(xne/xni(1))*(sindene+sindeno(1))/fp
c
c     	vphicos(2,n)=xlphim(n)*te(n)*cosdene/(vth(2)*bthet) - presshat(2) 
c	1	-1.*(1.+cosden(2))*vtheta(2)
c	2	-hnuion(2)*(xne/xni(2))*(sindene+sindeno(2))/fp
c
c	chk(1) = xlphim(n)*te(n)*cosdene/(vth(2)*bthet)
c     	chk(2) =-1.*(1.+cosden(2))*vtheta(2)
c
c	vphisin(1,n)=xlphim(n)*te(n)*sindene/(vth(1)*bthet) - 
c     2	sinden(1)*vtheta(1) +
c	1	(xne/xni(1))*hnuion(1)*(cosdene+cosdeno(1))/fp
c
c	vphisin(2,n)= xlphim(n)*te(n)*sindene/vth(2) - 
c     2	sinden(2)*vtheta(2) +
c	1	(xne/xni(2))*hnuion(2)*(cosdene+cosdeno(2))/fp

	vphicos(1,n) = xlphim(n)*te(n)*cosdene/(vth(1)*bthet)
     1		-hnuion(1)*(xne/xni(1))*(sindene+sindeno(1))/fp
     2		+(1.+cosden(1))*((2.*erex(n)/(vth(1)*bthet))
     3		-2.*presshat(1)-vtor1(n)/vth(1))
     	vphicos(2,n) = xlphim(n)*te(n)*cosdene/(vth(2)*bthet)
     1		-hnuion(2)*(xne/xni(2))*(sindene+sindeno(2))/fp
     2		+(1.+cosden(2))*((2.*erex(n)/(vth(2)*bthet))
     3		-2.*presshat(2)-vtor2(n)/vth(2))
	vphisin(1,n) = xlphim(n)*te(n)*sindene/(vth(1)*bthet)
     1	+	(xne/xni(1))*hnuion(1)*(cosdene+cosdeno(1))/fp
	2	-	sinden(1)*(2.*presshat(1) + vtor1(n)/vth(1)
	3	-	2.*erex(n)/(vth(1)*bthet))
	vphisin(2,n) = xlphim(n)*te(n)*sindene/(vth(2)*bthet)
     1	+	(xne/xni(2))*hnuion(2)*(cosdene+cosdeno(2))/fp
     2	-	sinden(2)*(2.*presshat(2) + vtor2(n)/vth(2)
     3	-	2.*erex(n)/(vth(2)*bthet))


	VPHIS1 = VPHISIN(1,N)
 	VPHIS2 = VPHISIN(2,N)
	VPHIC1 = VPHICOS(1,N)
 	VPHIC2 = VPHICOS(2,N)

	GOTO 499
396	CONTINUE

	GOTO 70

C	NEW CALC OF PHI_C,S FROM RAD MOM BAL & VPHI_C,S FROM TOROIDAL MOM BAL
c	toroidal asymmetry solution	
	gyrohat1 =0.5*(ti(n)*rminor/(bphi*rmajor))*(xlpm(n)+xlvm(n))
c	*********for interpretation are using measured nudtot at xnuanom
	gyrohat1 = 0.
	gyrohat2 = gyrohat1/zbar2(n)
	
	hnucol(1) = rminor*xnuc12(n)/vth(1)
	hnucol(2) = rminor*xnuc21(n)/vth(2)

	
c	*************using xnuamom to represent measured drag less atomic************
	xnuanom1 = xnudtot1(n) - coldno(n)*svata(n) - xnuioni(n) 
     1										     - xnuionb(n)
      xnuanom2 = xnudtot2(n)  
	
	hnuanom1 = rminor*xnuanom1/vth(1)
	hnuanom2 = rminor*xnuanom2/vth(2)

	zz1 = xmompol(n,1)/(yni(n,1)*xmas(1)*(vth(1)**2))
	zz2 = xmompol(n,2)/(yni(n,2)*xmas(2)*(vth(2)**2))


	VPHIS1 = VPHISIN(1,N)
	VPHIS2 = VPHISIN(2,N)
	VPHIC1 = VPHICOS(1,N)
	VPHIC2 = VPHICOS(2,N)
c	cos mom species 1
	c(1,1) = rmajor*gyrohat1/vth(1) + 0.5*fp*vtheta(1)
	c(1,2) = 1.5*(fp**2)*f(1)*qsafe/ep - 0.5*(vrad1(n)/vth(1))*
	1		(rminor*xlvm(n)) + 0.5*(hnuatom(1)+hnucol(1)+hnuanom1)
	c(1,3) = 0.0
	c(1,4) =-0.5*hnucol(1)*sqrt(xmas(1)/xmas(2))
	s(1) = 0.5*vrhat(1)*(ep*fp/qsafe)*(1. + cosden(1))
	s(1) = s(1) - 0.5*vphihat(1)*(hnuatom(1)+hnuamom1)*(2.0+cosden(1))
	s(1) = s(1) + rminor*zz1
	s(1)=s(1)-0.5*hnucol(1)*(vphihat(1) - sqrt(xmas(1)/xmas(2))*
     1	vphiHat(2))*(2.+ cosden(1) + cosden(2))
	s(1) = s(1) + 0.5*(vrad1(n)/vth(1))*vphihat(1)*
	1			(rminor*xlvm(n)*(2.+ cosden(1)) - 1.0)
      s(1) = s(1) + 0.5*(qsafe*f(1)*(fp**2)/ep)*
	1		(-1.*(cosden(1)+3.)*vtheta(1)-fp*hnuion(1)*(xne/xni(1))*
     2		(sinden(2)+sindeno(1)) + 3.*vphihat(1)) 
c	cos mom species 2
     	c(3,3) = rmajor*gyrohat2/vth(2) + 0.5*fp*vtheta(2)
 	c(3,4) = 1.5*(fp**2)*f(2)*qsafe/ep - 0.5*(vrad2(n)/vth(2))*
	2		rminor*xlvm(n)
	1		+0.5*(hnuatom(2)+hnucol(2)+hnuanom2)
	c(3,1) = 0.0
	c(3,2) =-0.5*hnucol(2)*sqrt(xmas(2)/xmas(1))
	s(3) = 0.5*vrhat(2)*(ep*fp/qsafe)*(1. + cosden(2))
	s(3) = s(3)-0.5*vphihat(2)*(hnuatom(2)+hnuamom2)*(2.0 + cosden(2))
	s(3) = s(3) + rminor*zz2
	s(3)=s(3)-0.5*hnucol(2)*(vphihat(2)-sqrt(xmas(2)/xmas(1))*
     1	vphiHat(1))*(2.+ cosden(1) + cosden(2))
	s(3) = s(3) + 0.5*(vrad2(n)/vth(2))*vphihat(2)*
	1			(rminor*xlvm(n)*(2.+ cosden(2)) - 1.0)
      s(3) = s(3) + 0.5*(qsafe*f(2)*(fp**2)/ep)*
	1		(-1.*(cosden(2)+3.)*vtheta(2)-fp*hnuion(2)*(xne/xni(2))*
     2		(sinden(1)+sindeno(2)) + 3.*vphihat(2)) 
 

c	sin mom species 1
	c(2,2) = -1.0*c(1,1)
	c(2,1) =      c(1,2)
	c(2,4) =  0.0
      c(2,3) =      c(1,4)
	s(2) = 0.5*vrhat(1)*(ep*fp/qsafe)*sinden(1)
 	s(2) = s(2) - 0.5*vphihat(1)*(hnuatom(1)+hnuamom1)*sinden(1)
	s(2) = s(2) - (rmajor*gyrohat1/vth(1))*vphihat(1)
	s(2)=s(2) - 0.5*hnucol(1)*(vphihat(1)-sqrt(xmas(1)/xmas(2))*
     1	vphiHat(2))*(sinden(1) + sinden(2))
	s(2) = s(2) + 0.5*(vrad1(n)/vth(1))*vphihat(1)*
	1			(rminor*xlvm(n)*sinden(1)) 
	s(2) = s(2) + 0.5*(qsafe*f(1)*(fp**2)/ep)*(-1.*sinden(1)*vtheta(1)
	1		+ fp*hnuion(1)*(xne/xni(1))*(cosden(2)+cosdeno(1)))
	s(2) = s(2) +0.5*fp*vtheta(1)*vphihat(1)
c	sin mom species 2
	c(4,4) = -1.0*c(3,3)
 	c(4,3) =      c(3,4)
	c(4,2) =  0.0
      c(4,1) =      c(3,2)
	s(4) = 0.5*vrhat(2)*(ep*fp/qsafe)*sinden(2)
	s(4) = s(4) - 0.5*vphihat(2)*(hnuatom(2)+hnuamom2)*sinden(2)
 	s(4) = s(4) - 0.5*hnucol(2)*(vphihat(2)-sqrt(xmas(2)/xmas(1))*
     1	vphiHat(1))*(sinden(1) + sinden(2))
	s(4) = s(4) + 0.5*(vrad2(n)/vth(2))*vphihat(2)*
	1			(rminor*xlvm(n)*sinden(2)) 
	s(4) = s(4) + 0.5*(qsafe*f(2)*(fp**2)/ep)*(-1.*sinden(2)*vtheta(2)
 	1		+ fp*hnuion(2)*(xne/xni(2))*(cosden(2)+cosdeno(2)))
 	s(4) = s(4) +0.5*fp*vtheta(2)*vphihat(2)

	CALL LSLRG(4,C,4,S,1,A)
c	these are the products (vphi/vth)*(vphi_s,c/ep) ************************
	vphisin(1,n) = A(1)
	VPHICOS(1,N) = A(2)
	VPHISIN(2,N) = A(3)
	VPHICOS(2,N) = A(4)

	
	VPHICOS(1,N) = (1.-EPD)*VPHICOS(1,N) + EPD*VPHIC1
	VPHICOS(2,N) = (1.-EPD)*VPHICOS(2,N) + EPD*VPHIC2				 
	VPHISIN(1,N) = (1.-EPD)*VPHISIN(1,N) + EPD*VPHIS1
	VPHISIN(2,N) = (1.-EPD)*VPHISIN(2,N) + EPD*VPHIS2

	VPHIS1 = VPHISIN(1,N)
 	VPHIS2 = VPHISIN(2,N)
	VPHIC1 = VPHICOS(1,N)
	VPHIC2 = VPHICOS(2,N)


C	POTENTIAL ASYMMETRY FROM RAD MOM BAL
	VT1 = VTHETA(1)*VTH(1)*FP
	VT2 = VTHETA(2)*VTH(2)*FP  

	COSPHI(1) =-1.0*(1.+COSDEN(1)) + 
	1	(BTHET*(VPHICOS(1,N)+VPHIA(1)*COSDEN(1)) +
     2 BPHI*(VT1+HNUION(1)*VTH(1)*(XNE/XNI(1))*(SINDEN(2)+SINDENO(1))) +
	3	PRESS(1)*BTHET*(1.+COSDEN(1)))/RADE
	COSPHI(2) =-1.0*(1.+COSDEN(2)) + 
 	1	(BTHET*(VPHICOS(2,N)+VPHIA(2)*COSDEN(2)) -
     2 BPHI*(VT2+HNUION(2)*VTH(2)*(XNE/XNI(2))*(SINDEN(1)+SINDENO(2))) +
	3	PRESS(2)*BTHET*(1.+COSDEN(2)))/RADE
	SINPHI(1) = -1.0*SINDEN(1) +
	1	(BTHET*(VPHISIN(1,N)+VPHIA(1)*SINDEN(1)) -
     2	BPHI*HNUION(1)*VTH(1)*(XNE/XNI(1))*(COSDEN(2)+COSDENO(1)) +
 	3	PRESS(1)*BTHET*(SINDEN(1)))/RADE
	SINPHI(2) = -1.0*SINDEN(2) +
	1	(BTHET*(VPHISIN(2,N)+VPHIA(2)*SINDEN(2)) -
     2	BPHI*HNUION(2)*VTH(2)*(XNE/XNI(2))*(COSDEN(1)+COSDENO(2)) +
 	3	PRESS(2)*BTHET*(SINDEN(2)))/RADE

C	sinphi(1) = sinden(1)*(-1. + vth(1)*bthet*vphihat(1)/rade 
C	1	+ press(1)/rade) - 
C     2	cosden(2)*(vth(1)*bphi*hnuion(1)*(xne/xni(1))/rade) +
C	3(vth(1)/rade)*(bthet*vphis1-bphi*hnuion(1)*(xne/xni(1))*
C     4	cosdeno(1))
C     	sinphi(2) = sinden(2)*(-1.+vth(2)*bthet*vphihat(2)/rade 
C     1	+ press(2)/rade) - 
C     2	cosden(1)*(vth(2)*bphi*hnuion(2)*(xne/xni(2))/rade) +
C	3(vth(2)/rade)*(bthet*vphis2-bphi*hnuion(2)*(xne/xni(2))*
C     4	cosdeno(2))
C	cosphi(1) = cosden(1)*(-1.+vth(1)*bthet*vphihat(1)/rade 
C	1	+ press(1)/rade) + 
C     2	sinden(2)*(vth(1)*bphi*hnuion(1)*(xne/xni(1))/rade) +
C	3	(-1. + (vth(1)/rade)*(bthet*vphic1+fp*bphi*vtheta(1)+
C	4	bphi*hnuion(1)*(xne/xni(1))*sindeno(1)) + press(1)/rade)
C 	cosphi(2) = cosden(2)*(-1. + vth(2)*bthet*vphihat(2)/rade 
C	1	+ press(2)/rade) + 
C     2	sinden(1)*(vth(2)*bphi*hnuion(2)*(xne/xni(2))/rade) +
C	3	(-1. + (vth(2)/rade)*(bthet*vphic2+fp*bphi*vtheta(2)+
C	4	bphi*hnuion(2)*(xne/xni(2))*sindeno(2)) + press(2)/rade)


      
	COSPOT = COSPHI(1)
	SINPOT = SINPHI(1) 
70	CONTINUE
C		BOLTZMAN POTENTIAL POLOIDAL VARIATION
 	cospot = (gam(1)*cosden(1)+gam(2)*cosden(2))/phite
	sinpot = (gam(1)*sinden(1)+gam(2)*sinden(2))/phite

 


499	CONTINUE 
	
	
c	10/7/05
	xrad1 = ep*((1+cosden(1)+vthetcos(1))-2.*xlvm(n)) 
	xrad2 = ep*((1+cosden(2)+vthetcos(2))-2.*xlvm(n))
	xtht1=ep*(vphisin(1,n)*(1.+vthetcos(1))-vthetsin(1)*
	1	(1.+vphicos(1,n)))
	xtht2=ep*(vphisin(2,n)*(1.+vthetcos(2))-vthetsin(2)*
	1	(1.+vphicos(2,n)))
	xnuinert1(n) = 0.5*(vrad1(n)*xrad1+vtheta(1)*xtht1)/rmajor 
	xnuinert2(n) = 0.5*(vrad2(n)*xrad2+vtheta(2)*xtht2)/rmajor      
C	FROM VTHETAHAT S/C TO VTHETA S/C 	 
	VTHTCOS(1,N) = VTHETCOS(1)*EP*VTHETA(1)*FP*VTH(1)
	VTHTCOS(2,N) = VTHETCOS(2)*EP*VTHETA(2)*FP*VTH(2)
	VTHTSIN(1,N) = VTHETSIN(1)*EP*VTHETA(1)*FP*VTH(1)
	VTHTSIN(2,N) = VTHETSIN(2)*EP*VTHETA(2)*FP*VTH(2)
	VTHTCOS(1,N) = VTHETCOS(1)
	VTHTCOS(2,N) = VTHETCOS(2)
	VTHTSIN(1,N) = VTHETSIN(1)
	VTHTSIN(2,N) = VTHETSIN(2)
c	inertial pinch term
	pinert(n) = (xmas(1)*vphia(1)*vrad1(n)/rmajor)*
	1			((rmajor*xlvm(n))-1.)/(eq*atnum(1)*bthet)
	pinert(n) = pinert(1) + 0.5*(xmas(1)*vphia(1)*vth(1)*fp*ep/
	1	rmajor)*(vthetsin(1)*(3.+cosden(1))-vthetcos(1)*sinden(1) +
     2	vtheta(1)*sinden(1))/(eq*atnum(1)*bthet)
	pinert(n) = pinert(1) + 0.5*(xmas(1)*vtheta(1)*(vth(1)**2)*
     1		fp*ep/rmajor)*vphisin(1,n)/(eq*atnum(1)*bthet)

2001	xv(1) = 0.25*(ep**2)*(vphisin(1,n)*vphicos(1,n) + 
     1  (cosden(1)*vphisin(1,n) + sinden(1)*vphicos(1,n))*
     2	(vphia(1)/VTH(1)))
     	xv(2) = 0.25*(ep**2)*(vphisin(2,n)*vphicos(2,n) + 
     1  (cosden(2)*vphisin(2,n) + sinden(2)*vphicos(2,n))*
     2	(vphia(2)/VTH(2)))

C		del**2 term in inertial source term in sin moment     
      yv(1)=0.5*(ep**2)*(0.5*(vphicos(1,n)**2)+1.5*(vphisin(1,n)**2)+
     1 (cosden(1)*vphicos(1,n) + 3.*sinden(1)*vphisin(1,n))*
     2	vphia(1)/VTH(1))			 
	yv(2)=0.5*(ep**2)*(0.5*(vphicos(2,n)**2)+1.5*(vphisin(2,n)**2)+
     1 (cosden(2)*vphicos(2,n) + 3.*sinden(2)*vphisin(2,n))*
     2	vphia(2)/VTH(2)) 

C		ep**2 term in friction in unity moment
 	fric(1) = -0.5*(EP**2)*(1.+COSDEN(1)+COSDEN(1)**2+SINDEN(1)**2)
	fric(2) = -0.5*(EP**2)*(1.+COSDEN(2)+COSDEN(2)**2+SINDEN(2)**2)

c	viscous quadratic terms in sin & cos eqs
 	do 45 j=1,2 
	dela(j) = -1.*(ep**2)*sinden(j)*vphisin(j,n)/4.
 
c	cc1(j)=0.125*(ep**2)*((sinden(j)**2)/3. - cosden(j)/3. - 
c	2		(cosden(j)**2)/3. -5.)
c	cc2(j)=0.25*(ep**2)*((1.-cosden(j))*((vphi(j)/cs(j))-vphicos(j,n)) 
c	1	+ 	sinden(j)*vphisin(j))
c	sc1(j)= 0.125*(ep**2)*sinden(j)*(3.+ 2.*cosden(j))
c	sc2(j)= 0.75*(ep**2)*((vphicos(j,n)-(vphi(j)/cs(j)))*sinden(j) -
c     1	 vphisin(j,n)*(1.+cosden(j)))
	vc1(j) = (ep**2)*(cosden(j)+1.)/24.
	vc2(j) = -1.*(qsafe**2)*f(j)*fp*(ep**2)*
	1		((vphia(j)/VTH(j))-vphicos(j,n))/4. 
	vc3(j) = (qsafe**2)*f(j)*fp*(ep**2)*
	1		(vtheta(j)*sinden(j)/24.+vphisin(j,n)/4.)
      vc4(j) = (qsafe**2)*f(j)*fp*(ep**2)*
     1		(5.*vtheta(j)/8. - ((vphia(j)/VTH(j))-vphicos(j,n))/4.)
      vs1(j) = (ep**2)*(3./8. + cosden(j)/4.)
 	vs2(j) = -1.*(qsafe**2)*f(j)*fp*(ep**2)*
	1		((vphia(j)/VTH(j))-vphicos(j,n))*3./4. 
 	vs3(j) = -1.*(qsafe**2)*f(j)*fp*(ep**2)*vphisin(j,n)*3./4.
      vs4(j) = (qsafe**2)*f(j)*fp*(ep**2)*vphisin(j,n)*3./4.	 
45	continue 
			  
75	continue
c	*********converge the asymmetries on each vtheta iterate***************
	esym =0.001

c	XX = ABS(PS1/VPHIS1 - 1.)
c 	IF(XX.GT.esym) GOTO 91
c	XX = ABS(PS2/VPHIS2 - 1.)
c	IF(XX.GT.esym) GOTO 91
c	XX = ABS(PC1/VPHIC1 - 1.)
c 	IF(XX.GT.esym) GOTO 91
c 	XX = ABS(PC2/VPHIC2 - 1.)
c 	IF(XX.GT.esym) GOTO 91

	GOTO 91
90	continue
	WRITE (6,900) 'POLOIDAL ASYMMETRY SOLUTION NOT CONVERGED',nmesh 
 
91	continue	
c	write(6,301) vphia(1),vtheta(1),vtheta(2),sinden(1),sinden(2),
c     1			cosden(1),cosden(2)	
	if(ioptedgereal.eq.3) ioptedge = 3
c	if(nn.gt.100) zc = 1.
	vt(nn) = vtheta(1) 
c	do not iterate the nmesh = 25 solution (which diverges)
c	if(nmesh.eq.25) goto 150
C		CHECK CONVERGENCE OF ROTATION/ASYMMETRY CALCULATION
	if(ioptvthet.eq.1) goto 150 
	EPOL = 0.01
	IF(NN.LT.1) GOTO 100
	XX = ABS(VTHET1/VTHETA(1) - 1.)
 	IF(XX.GT.epol) GOTO 100
	if(ioptedge.ne.3) then
	XX = ABS(VTHET2/VTHETA(2) - 1.)
 	IF(XX.GT.epol) GOTO 100
	endif
c	XX = ABS(SIN1/SINDEN(1) - 1.)
c	IF(XX.GT.esym) GOTO 100
c	XX = ABS(SIN2/SINDEN(2) - 1.)
c	IF(XX.GT.esym) GOTO 100
c	XX = ABS(COS1/COSDEN(1) - 1.)
c 	IF(XX.GT.esym) GOTO 100
c 	XX = ABS(COS2/COSDEN(2) - 1.)
c 	IF(XX.GT.esym) GOTO 100
c	XX = ABS(CPOTC/COSPOT - 1.)
c  	IF(XX.GT.esym) GOTO 100
c 	XX = ABS(SPOTS/SINPOT - 1.)
c 	IF(XX.GT.esym) GOTO 100

c	if(zc.ge.1.0) then
	
      goto 150 
c	endif
c	zcc = 1.
c	c9 = 1.
100	continue	
900	format((1X,A40,i4))
	WRITE (6,900) 'POLOIDAL ROTATION SOLUTION NOT CONVERGED',nmesh 
150	CONTINUE
	
	
	 
c	if(ioptpert.eq.0) goto 155 
c	if(ioptpert.eq.2) goto 155 
c	if(niter.gt.1) goto 155 
c 	if(ioptpert2.eq.2) then
c  	ioptpert = 2
c	endif
c	if(ioptpert.eq.1) goto 155 
c	goto 5
155	continue 
C	UPDATE DRAG FREQUENCY
	par = (1.-((RMINOR/AMINOR)**2))	 
	dxn0 = xn0-xnped
	dxt0 = t0-tped
	dxv0 = vphexp-vped
	DO 750 J = 1,2
c	xlvm(n) = -1.*(torv(n+1)-torv(n))/(delma*(torv(n+1)+torv(n))) 


c	**************************************
c	if(it.gt.1) xlna(j) = 1./xlnm(nmesh)
c	xlva(j) = 1./xlvm(nmesh)
c	xlva(j) = xlv1
c	****************************************
c	XLINV=AMINOR*((1./XLNa(J))+(1./XLTa(J))+(1./XLVA(J)))
	
c	xlinv=aminor*(xlvm2(nmesh)+xlpm(nmesh))
c	if(j.eq.1) xlinv=aminor*(xlvm(nmesh)+xlpm(nmesh))

c	GRAD(J) = XLINV
c	xlvm(n) = 1./ylvbarx
	grad(j) = rhor(n)*AMINOR*(xlpm(n)+xlvm(n))

	epd = 0.0
	zz = vphia(j)/vth(j)
c	thetinert(j) = -1.*vthetcos(j)*(vphisin(j,n)+sinden(j)) +
c	1	vthetsin(j)*(2.+vphicos(j,n)+cosden(j)) +
c     2	vthetsin(j)+vphisin(j)+sinden(j)	
	
  
      THETw(J) = EPD*THETw(J)+(1.-EPD)*((4. + COSDEN(J))*vphisin(j,n)/zz 
	1			+ SINDEN(J)*(1.-vphicos(j,n)/zz))
c	thetw(j) = (4.+cosden(j))*(sinpot-(vtheta(j)*vth(j)/vphia(j))*
c	1	(sinpot+sinden(j))) + sinden(j)*((vtheta(j)*vth(j)/vphia(j))*
c	2		(2.+cospot+cosden(j))-cospot)
750	continue
	do 751 j=1,2
	k=2
	zz = ATNUM(J)
 	if(j.eq.2) THEN
	k = 1
	zz = ZBAR2(N)
	ENDIF

	XNUDRAG(J)=GRAD(J)*THETw(J)*TEMP(J)/
	2			(2.*(RMAJOR**2)*zz*(BPHI))
	if(xnudrag(j).le.0.0.and.ioptzerodrag.eq.0) xnudrag(j) = 0.0
c	add anomalous drag
c	if(ioptdrag.eq.1) xnudrag(j) = xnudrag(j) + anomdrag(j)   
751	continue
c*********************change***********************
c	if(ioptedge.eq.100) ioptedge = 1
	Gy(n) = grad(1)
	cosion(n) = cosden(1)
	cosimp(n) = cosden(2)
	sinion(n) = sinden(1)
	sinimp(n) = sinden(2) 
	thetwid1(n) = thetw(1)
	thetwid2(n) = thetw(2)
c	frequencies	for toroidal velocity integration
	m = n 
 	zz = 1.
	XNUDRAGp(m,1)=xlpm(m)*THETwid1(m)*TI(m)/
	1		(2.*(RMAJOR**2)*zz*(BPHI))
	zz = zbar2(m)
	XNUDRAGp(m,2)=xlpm(m)*THETwid2(m)*TI(m)/
 	1		(2.*(RMAJOR**2)*zz*(BPHI))
	xx = vtheta(1)*vth(1)*fp
	yy = vtheta(2)*vth(2)*fp
 
	ynuinert(m,1) = 0.5*ep*((vrad1(m)/rmajor)*
	1				(1.+cosion(m)+vphicos(1,m))+
	2	(xx/rmajor)*(vphisin(1,m)*(1.+vthetcos(1)+cosion(m))-
	3				 vthetsin(1)*(1.+vphicos(1,m))-
     4				 vphicos(1,m)*sinion(m)))	
      ynuinert(m,2) = 0.5*ep*((vrad2(m)/rmajor)*
     1				(1.+cosimp(m)+vphicos(2,m))+
	2	(xx/rmajor)*(vphisin(2,m)*(1.+vthetcos(2)+cosimp(m))-
	3				 vthetsin(2)*(1.+vphicos(2,m))-
     4				 vphicos(2,m)*sinimp(m)))	
	xnudragatomic(j) = xnuioni(j) + xnuati(j) + xnuionb(j) 
	zz = 1. 
	chiphi(m,1)=0.5*thetwid1(m)*rhor(m)*ti(m)/
	1	((RMAJOR**2)*zz*(BPHI))	- vrad1(m)
      zz = zbar2(m)
	chiphi(m,2)=0.5*thetwid2(m)*rhor(m)*ti(m)/
	1			((RMAJOR**2)*zz*(BPHI))	- vrad2(m)
	coefv(m,1) = (chiphi(m,1)/delna)+xnudragp(m,1)+ynuinert(m,1)+
     1				xnuc12(m)+ xnudragatomic(j) +
     2				xnudraganom1(n)	 
	coefv(m,2) = (chiphi(m,2)/delna)+xnudragp(m,2)+ynuinert(m,2)+
     1				xnuc21(m) + xnudraganom2(n)


 				  
800	return
	END 	