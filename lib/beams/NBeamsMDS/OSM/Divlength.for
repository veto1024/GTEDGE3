	SUBROUTINE DIVLENGTH
	INCLUDE 'SOLDIV.FI' 
	real zds
C	CALCULATES AVERAGE FIELD LINE LENGTH IN DIVERTOR
C	ZP=PLASMA CNTRLN HGHT, ZM=DIV COIL CNTRLN HGHT, ZX=XPT HGHT
c	ZDS=DIVERTOR STRIKE PT HGHT

	xdsep =  rsep1 - rx
	zhght =  abs(zx-zsep1)
	
     
	goto 250
c				
	ALPHADIV = XDSEP/ZHGHT 
	DELN = EPDIV*DELN 
	ZP = 0.0
	ZXOLD = ZX
	ZX = -1.*ELONG*AMINOR
	IF(IOPTDIV.EQ.1.AND.IOPTSN.EQ.2) ZX = ELONG*AMINOR
	ZM = 1.5*ZX 
	DZ = ZHGHT/100.
	IF(IOPTDIV.EQ.2.OR.IOPTDIV.EQ.1.AND.IOPTSN.EQ.1) THEN 
	ZDS = ZX - ZHGHT
	Z = ZDS + 0.5*DZ
	ENDIF 
	IF(IOPTDIV.EQ.1.AND.IOPTSN.EQ.2) THEN
	ZDS = ZX + ZHGHT	
	Z = ZDS - 0.5*DZ
	ENDIF
	DX = DELN/25.
	ZINT = 0.0 
	DO 200 I =1,100
	XINT = 0.0 
	X0 = 0.5*DX
	DO 100 J =1,25
	X1 = X0 +ALPHADIV*(ZX-Z)
	H0 = (Z-ZP)**2 + X1**2
	HM = (Z-ZM)**2 + X1**2
	HINT = SQRT(H0*HM)*DX/((ZX-ZM)*SQRT(H0)-(ZP-ZX)*SQRT(HM))
	XINT = XINT + HINT/DELN	
	X0 = X0 + DX
100	CONTINUE     
     	ZINT = ZINT + XINT*DZ/ZHGHT
	IF(IOPTDIV.EQ.2.OR.IOPTDIV.EQ.1.AND.IOPTSN.EQ.1) Z = Z + DZ
	IF(IOPTDIV.EQ.1.AND.IOPTSN.EQ.2) Z = Z - DZ
200	CONTINUE
	XENHANCE1 = (ELONG/SQRT(0.5*(1.+ELONG**2)))*ZINT 
c	calculate Lphi at each x0 and average
250	alphadiv = xdsep/zhght
      deln = deln*epdiv
	zp = 0.0
	sign = 1.
	if(ioptdiv.eq.2.or.ioptsn.eq.1)   sign = -1.
c		the next statement is a trick used to calc. the usn case 
 	if(ioptdiv.eq.1.and.ioptsn.eq.2)  sign = -1.
	zx = sign*elong*aminor
	zm = 1.5*zx
	dx = deln/25.
	dz = zhght/100.
 	xint = 0.0
	x0 = 0.5*dx 
	do 400 i=1,25
	zint = 0.0
	z = zx + sign*0.5*dz
	do 300 j=1,100
	x1 = x0 + sign*alphadiv*(z-zx)
	hpz = sqrt((z-zp)**2 + x1**2)
	hpx = sqrt((zx-zp)**2 + x1**2)
	hmz = sqrt((z-zm)**2 + x1**2)
	hmx = sqrt((zx-zm)**2 + x1**2)	 
	hxz = sqrt((z-zx)**2 + x1**2)
	denom = ((hxz*hmx/(hpx*hmz)) - hxz/hpz) 
	hint = (alphadiv*x0+(1.+alphadiv)*(z-zx))/denom
	zint = zint + sign*hint*dz/(aminor*sqrt(0.5*(1.+elong**2)))
	z = z + sign*dz
300	continue
	xint = xint + zint*dx/deln
	x0 = x0 + dx
400	continue
	xenhance2 = xint
	
      
	
C	correct beta based on plasma Btheta to Btheta at divertor strike pt
	x0 = 0.0
	x1 = x0 + sign*alphadiv*zhght
	zd = sign*zhght
	hmd = sqrt((zd-zm)**2 + x1**2)
	hpd = sqrt((zd-zp)**2 + x1**2)
	betacorrect = aminor*sqrt(0.5*(1.+elong**2))*((hmx/(hpx*hmd))-
	2											  (1./hpd))	
	betacorrect = abs(betacorrect)

	xenhance = xenhance2*betacorrect

     	DELN = DELN/EPDIV
     


      RETURN 
	END
