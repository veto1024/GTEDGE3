	SUBROUTINE PARAM(NMESH)
	INCLUDE 'SOLDIV.FI'
	parameter (jq=2)
	dimension xlp(jq),G(Jq),xmu(jq),
     1			GRAD(JQ)			
c	sign toroidal rotation
298	FORMAT(1X,'OK') 
	atnum2 = atnum(2)
	atnum(2) = zbar2(nmesh) 
	signv = 1.
	if(pbeam.lt.0.0) signv = -1. 
C  CALCULATES ROTATION & IMPURITY TRANSPORT PARAMETERS
C	BASIC PARAMETERS
	se = sqrt(0.5*(1.+elong**2)) 
	EP = AMINOR*se/RMAJOR
	rminor = aminor*se 
	XMPROT = 1.673E-27 
	XMASSELECTRON = 9.11E-31
	XMAS(1) = AION*XMPROT
	XMAS(2) = AIMP*XMPROT
	ATNUM(1) = ZION
	ATNUM(2) = ZIMP
 	EQ = 1.6E-19
	XK = 1.6E-19
	EP0 = 8.854E-12
	
	IF(IZEFF.EQ.1) XNI(2) = (ZEFF-1.)*XNI(1)/(ZIMP*(ZIMP-ZEFF))
	ALPHA = XNI(2)*(ZIMP**2)/(XNI(1)*(zion**2))
	
C		CYCLOTRON FREQUENCIES & GYRORADII
	DO 10 J = 1,2
	OMII(J) = ATNUM(J)*EQ*abs(BPHI)/XMAS(J)
c	if(j.eq.2) omii(j) = zbar2(nmesh)*eq*abs(bphi)/xmas(j) 
	VTH(J) = SQRT(2.*XK*TEMP(J)/XMAS(J))
	RHOTi(J) = VTH(J)/OMII(J)
10	CONTINUE
	rhot1(nmesh) = rhoti(1) 
	CSE = SQRT(2.*XK*TEP/XMASSELECTRON)
c	effective central values consistent with profiles
c	if(igscl.eq.1) then 
c	t0 = (temp(1)-tped)/((0.75)**alpht) + tped   
c	xne0 = ((atnum(1)*xn(1)+atnum(2)*xn(2))-xnped)/((0.75)**alphn) + 
c    2		xnped 
c	endif

 
C		POLOIDAL FIELD
		QSAFE = Q95
		BTHET = EP*abs(BPHI)/QSAFE
		FP = BTHET/BPHI
	
	DO 20 J = 1,2
C		COLLISION FREQUENCIES
	DO 20 K = 1,2 
C			COULOMB LOGARITHM
	Y = SQRT(XNI(K))*(ATNUM(K)**2)*ATNUM(J) 
	X = (EP0/EQ)**1.5	
	COULOG(J,K) = LOG(12.*3.1416*(TEMP(J)**1.5)*X/Y)
C			REDUCED MASS
	XMR = XMAS(J)*XMAS(K)/(XMAS(J)+XMAS(K))
		
C			FREQUENCIES
					 
 	XNUC(J,K)=XNI(K)*(ATNUM(K)**2)*(ATNUM(J)**2)*(EQ**2.5)*COULOG(J,K)
	2		  / (47.25*SQRT(XMR)*(EP0**2)*(TEMP(J)**1.5))
c	IF(K.EQ.J) XNU(J,K) = XNU(J,K)/1.414
c	CONST = 1.E-13*(0.7*SQRT(XMASS(J)*XMPROT)/XMASS(K) +
c    2		0.9*SQRT(XMPROT/XMASS(J)))	 
c	XNU(J,K)=CONST*XN(K)*((ATNUM(J)*ATNUM(K))**2)*COULOG(J,K)/TEMP(J) 
20	CONTINUE
	
	IF(INU.NE.1) GOTO 21

c		BRAGINSKI COLLISION FREQUENCIES
	XMR11 = XMAS(1)*(1.+XMAS(1)/XMAS(1))
	XMR12 = XMAS(1)*(1.+XMAS(1)/XMAS(2))
 	XMR21 = XMAS(2)*(1.+XMAS(2)/XMAS(1))
	XMR22 = XMAS(2)*(1.+XMAS(2)/XMAS(2))

	C1 = 1./((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5)) 
	XNUC(1,1) = 3.34*(COULOG(1,1)*(ATNUM(1)**4)*1.E-6*XNI(1))/
	2			(C1*SQRT(XMR11*1E3)*(TEMP(1)**1.5))
	XNUC(1,2)=3.34*(COULOG(1,2)*((ATNUM(1)*ATNUM(2))**2)*1.E-6*XNI(2))
	2			/(C1*SQRT(XMR12*1E3)*(TEMP(1)**1.5))
 	XNUC(2,1)=3.34*(COULOG(2,1)*((ATNUM(1)*ATNUM(2))**2)*1.E-6*XNI(1))
     2			/(C1*SQRT(XMR21*1E3)*(TEMP(1)**1.5))
	XNUC(2,2) = 3.34*(COULOG(2,2)*(ATNUM(2)**4)*1.E-6*XNI(2))/
	2			(C1*SQRT(XMR22*1E3)*(TEMP(2)**1.5))
21	CONTINUE
 	

	XNELECTRON = XNI(1)*(ATNUM(1)**2) + XNI(2)*(zbar2(nmesh)**2) 
	
	XNUEI = XNELECTRON/(6.4E14*((1.E-3*TEP)**1.5)) 
	  	
C	write(6,150) xnu(1,1),xnu(1,2),xnu(2,1),xnu(2,2),XNUEI
 	 	
125	FORMAT(1X,'CHIEXP-I=',E9.3,1X,'CHI CHANG-HINTON=',E9.3)	 

C		NORMALIZED COLLISION FREQUENCIES 
C		generalized to multiple ions
	DO 30 J = 1,2
	XNUSTAR(J,J) = 0.0
	DO 25 K = 1,2
	IF(IVISC.NE.0) GOTO 22 
	XNUSTAR(J,J) = XNUSTAR(J,J) + XNUC(J,K)*ABS(QSAFE)*RMAJOR/VTH(J)
	IF(K.NE.J) XNUSTAR(J,K) =  XNUC(J,K)*ABS(QSAFE)*RMAJOR/VTH(J)
	GOTO 25
22	XNUSTAR(J,K) = XNUC(J,K)*ABS(QSAFE)*RMAJOR/VTH(J)
25	CONTINUE 
30	CONTINUE
	

	XNUEISTAR = XNUEI*ABS(QSAFE)*RMAJOR/CSE
	 
C		NORMALIZED 'IONIZATION FREQUENCIES'
	DO 35 J = 1,2 
	XNUIONSTAR(J) = SION(J)/((XNI(J)*VTH(J))/(QSAFE*RMAJOR))
	if(ioptedge.eq.1) xnuionstar(j)=xnuioni(j)/(VTH(J)/(QSAFE*RMAJOR))
35	CONTINUE	 		 
	

C		VISCOUS NORMALIZED FREQUENCIES

C			VISCOUS
	XLL = XMAS(2)/(ZIMP*XMAS(1)) 
	XV1 = 1.+ ALPHA*(1.-XLL)*((VPHIa(1)/VTH(1))**2)/(1.+ALPHA)
	XV2 = 1.+ (XLL-1.)*((VPHIa(2)/VTH(2))**2)/(XLL*(1.+ALPHA))
	xv1 = 1.0
	XV2 = 1.0  
	
	G(1) = (XNUC(1,1)/XNUC(1,2))*XV1
	G(2) = (XNUC(2,2)/XNUC(2,1))*XV2
	cvisc = 1.0
	DO 40 J = 1,2
c	if(j.eq.2) cvisc = 1.5 
	XMU(J) = 1.5*SQRT(EP)*G(J)/
	2		 ((1.+XNUSTAR(J,J))*(1.+(XNUSTAR(J,J)*(EP**(1.5)))))
      F(J) = cvisc*XNUSTAR(J,J)/
	2		((1.+XNUSTAR(J,J))*((EP**(1.5)) + XNUSTAR(J,J)))

	
	ETA01HAT(NMESH) = F(1)*VTH(1)*QSAFE*rmajor
 	ETA02HAT(NMESH) = F(2)*VTH(2)*QSAFE*rmajor

 
40	CONTINUE
	XMU10 = XMU(1)/XV1
	XMU20 = XMU(2)/XV2
	XMU1 = 1.5*(EP**2)*G(1)*F(1)/XNUSTAR(1,2)
	XMU2 = 1.5*(EP**2)*G(2)*F(2)/XNUSTAR(2,1) 
	

	  
c	correct ion density to electron density
	xnelectron = xnI(1)*atnum(1) + xnI(2)*atnum(2)
	xnn = xnelectron 
c	XN0 = XNPED + (XNN-XNPED)/((1.-((RMINOR/AMINOR)**2))**ALPHN)
	xn0 = xne0/(atnum(1) - atnum(2)*xnI(2)/xnI(1))	
c 	V0 = VPED + (ABS(VPHI(1))-VPED)/((1.-((RMINOR/AMINOR)**2))**ALPHV)
	v0 = vphexp
C	T0 = TPED + (TEMP(1)-TPED)/((1.-((RMINOR/AMINOR)**2))**ALPHT)
c	pedestal-to-central ratios
C	ratnpedctr = xnped/xn0
C	rattpedctr = tped/t0
C	ratvpedctr = vped/v0 
150	format(6e10.3)	
202	format(6e10.3)
	
 
 
	 
C	PRESSURE GRADIENT	
C		DIFFERENCE	P'ION-P'IMP
c45	DELP = (TEMP(1)/BTHET)*(((1./(ZIMP*XLNa(2))) - (1./(XLNa(1))))
c     2		+ ((1./(ZIMP*XLTa(2))) - (1./(XLTa(1)))))	 
	
C		INDIVIDUAL
c	DO 50 J = 1,2 
c	xlp(j) = 1./((1./XLNa(J))+(1./XLTa(J)))

c	PRESS(J)=-1.*((TEMP(J)/(ATNUM(J)*BTHET)))/xlp(j)
C	press(j) = 0.
50	CONTINUE
	atnum(2) = atnum2 

	RETURN
	END