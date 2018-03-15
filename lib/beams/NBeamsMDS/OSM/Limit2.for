	SUBROUTINE LIMIT2
	INCLUDE 'SOLDIV.FI'
	PARAMETER (IINT = 4)

C		CALCULATES DENSITY LIMIT OF EQ.(16)


C	EFFECTIVE CIRCULAR RADIUS
	SE = SQRT((1.+ELONG**2)/2.)
C	RADIAL HEAT DIFFUSIVITIES FOR USE IN DISRUPTION ANALYSIS
C		INPUT CHIRC
	CHIRCORE(1) = CHIRC
	CHIRCORE(1) = 1.0
C		CONSISTENT HEAT REMOVAL CHIR USING EXP LT^-1 IN PEDESTAL
	CHIRCORE(2) = FCOND*FLUXHEAT/(XNPED*XK*TPED*GRADTBAR(3))
	CHIRCORE(2) = 0.5
C		CONSISTENT HEAT REMOVAL CHIR USING EXP LT^-1 IN TRANS BARRIER
	CHIRCORE(3) = FCOND*FLUXHEAT/(XNBAR*XK*TBAR*GRADTBAR(4))
	
      CHIRCORE(3) =0.25
C		CONSISTENT HEAT REMOVAL CHIR USING CORE NET HEAT REMOVAL

	 
	CHIRCORE(4) = FCOND*FLUXHEAT*AMINOR/
	2			  (XNAV*(T0-TPED)*XK) 
	CHIRCORE(4) = 0.1	 
	DO 60 J=1,4
	CHIR = CHIRCORE(J)
	DRCORE = CHIR
C	SOLVE QUADRATIC EQUATIONS	 
	  BB = CHIR*((5.5/(SE*AMINOR))**2)
	  X = 4.*(-1.*HEAT)*(-1.*RAD)/(BB**2)
	IF(X.LE.-1.) GOTO 25
	  Y = SQRT(1. + X)
	  DEN1(J) = -1.0*BB*(1 + Y)/(2.*RAD)
C		EVALUATE N-T INSTABILITY DENSITY LIMIT 
C		EQS. 39 & 40
25	A = 3.*TJ0*(SN - DRCORE*((5.5/(SE*AMINOR))**2)) + 2.*RADN
	B=3.*XNJ0*(SN-DRCORE*((5.5/(SE*AMINOR))**2))+3.*TJ0*ST-HEAT+RADN2
	C = 3.*XNJ0*ST
	YP = -1.*B*(1. + SQRT(1.-4.*A*C/(B**2)))/(2.*A)	  
	YM = -1.*B*(1. - SQRT(1.-4.*A*C/(B**2)))/(2.*A)
C		EQ. 44
C	XN0	= (1.+ALPHAN)*XNAV - ALPHAN*XNPED
      G = XNJ0/XN0
	FF = 1. + ALPHAN -ALPHAN*XNPED/XNAV
	AAP = (CHIR*((5.5/(SE*AMINOR))**2)*G + 2.*YP*RADN/XN0)/
	2       (-2.*RADN2/(XN0**2))
	AAM = (CHIR*((5.5/(SE*AMINOR))**2)*G + 2.*YM*RADN/XN0)/
	2       (-2.*RADN2/(XN0**2))
	BBP = (4.*HEAT*RADN2/(XN0**2))/((CHIR*((5.5/(SE*AMINOR))**2)*G
     2	+   2.*YP*RADN/XN0)**2)		 
	BBM = (4.*HEAT*RADN2/(XN0**2))/((CHIR*((5.5/(SE*AMINOR))**2)*G
     2	 +  2.*YM*RADN/XN0)**2)

	IF(BBP.LE.-1.) GOTO 50	
      DEN2PP(J) = AAP*(1. + SQRT(1. + BBP))/FF
	DEN2PM(J) = AAP*(1. - SQRT(1. + BBP))/FF
50	IF(BBM.LE.-1.) GOTO 60 
	DEN2MP(J) = AAM*(1. + SQRT(1. + BBM))/FF
	DEN2MM(J) = AAM*(1. - SQRT(1. + BBM))/FF
60	CONTINUE 
100	FORMAT (4E12.4) 
	RETURN
	END 	