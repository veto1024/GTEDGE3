	SUBROUTINE DISTR
	INCLUDE 'SOLDIV.FI'
C		DISTRIBUTES N & T

	AP = 39.44*RMAJOR*AMINOR*SQRT((1.+ELONG**2)/2.)
c	ap = 39.44*rmajor*aminor*sqrt(elong) 
 	VP = 19.72*RMAJOR*(AMINOR**2)*ELONG
	PSEP = XK*XNSEP*TSEP
 	CBAL = CBALLOON
	FAC = CBAL*((BFIELD**2)/(2.52E-6))*DELTB/((Q95**2)*RMAJOR*PSEP)
c	XNAVPED = (XNAV/XNEL)/APED
c      XNCTRPED = (XNAVPED*(1.+ALPHAN)-ALPHAN)
	
C 	XN0 = (1.+ALPHAN)*XNAV/(1.0 + ALPHAN/XNCTRPED)
C	ZT = 1.0/TCTRPED
C	ZN = 1.0/XNCTRPED
C	A1 = ((1.0+ALPHAN*ZN)/(1.0+ALPHAN))
C	A2 = ((1.0-ZT)*(1.0-ZN)/(ALPHAT2+ALPHAN+1.0))
C	A3 = ZT*(1.0-ZN)/(ALPHAN+1.0)
C	A4 = ZN*(1.0-ZT)/(ALPHAT2+1.0)
C	B1 = A1 + A2 + A3 + ZN*ZT
C	T0 = A1*TAV/B1
C      X = (1.+ALPHAN+ALPHAT2)*XNAV*TAV/(XNSEP*TSEP*(1.+FAC))
C      Y = (ALPHAT2/(1.+ALPHAN))*(XNCTRPED+ALPHAN*(2.+ALPHAN+ALPHAT2)/
C	2    (1.+ALPHAT2))
C	Z = XNCTRPED + ALPHAN/(1.+ALPHAT2)
	
c	RATNLINEAV = XNEL/XNAV
c	this expression valid for -1<alphn<1, including limits 1
	if(jjoptped.eq.10) xnctrped = xn0/xnped  
	if(jjoptped.eq.10) tctrped = t0/tped 
 	y = xnctrped - 1.
	
	ratnlineav = (1.+y*(1.-alphan/3. + alphan*(alphan-1.)/10. - 
     2		alphan*(alphan-1.)*(alphan-2.)/42.))/
	3		(1. + y/(1.+alphan))
	IF(JJOPTPED.NE.9) GOTO 10
 	XNPED = APED*XNAV*RATNLINEAV
c	ALPHAN = (XNCTRPED*(APED*RATNLINEAV)-1.)/(1.-APED*RATNLINEAV)
	
	XN0 = (1.+ALPHAN)*XNAV - ALPHAN*XNPED
 	TPED = XNSOL*TSOL*(1.+FAC)/XNPED
	alphat = alphat2
	X1 = (TAV/(APED*RATNLINEAV)) - ALPHAT*TPED*(1./(1.+ALPHAT) +
	2	 (XNCTRPED-1.)/((1.+ALPHAN)*(1.+ALPHAN+ALPHAT)))
     	X2 = 1./(1.+ALPHAT) + (XNCTRPED-1.)/(1.+ALPHAN+ALPHAT)
 	T0 = X1/X2 
     	GOTO 15
10	IF(JJOPTPED.NE.10) GOTO 15
C	XNPED = sqrt(1.+fac)*XNSOL
c	ALPHAN = (XNCTRPED*(APED*RATNLINEAV)-1.)/(1.-APED*RATNLINEAV)
C	IF(ITERN.GT.15)
C     2		ALPHAN = (XNCTRPED*(XNPED/XNAV)-1.)/(1.-XNPED/XNAV)
	XN0 = (1.+ALPHAN)*XNAV - ALPHAN*XNPED
C 	TPED = sqrt(1.+fac)*TSOL
	alphat = alphat2
c	X1 = (TAV/(APED*RATNLINEAV)) - ALPHAT*TPED*(1./(1.+ALPHAT) +
c 	2	 (XNCTRPED-1.)/((1.+ALPHAN)*(1.+ALPHAN+ALPHAT)))
c     	X2 = 1./(1.+ALPHAT) + (XNCTRPED-1.)/(1.+ALPHAN+ALPHAT)
	x1 = xnav*tav - (xn0-xnped)*tped*((1./(1.+alphan))-
	2	(1./(1.+alphan+alphat))) - xnped*tped*(1.- (1./(1.+alphat)))
      x2 = ((xn0-xnped)/(1.+alphan+alphat)) + xnped/(1.+alphat) 
 	T0 = X1/X2 
 
15	CONTINUE 
C	T0 = TPED*(X - Y)/Z 
	DELC = 1./56.
 
C	SE = SQRT(0.5*(1.+ELONG**2))
C	DELC = (1.0-DELTB/(AMINOR*SE))/56.0
	RHO(1) = 0.5*DELC
	TE(1) = TPED + A2*(T0-TPED)*((1.0 - (RHO(1)**2))**ALPHAT2) +
	2			   A1*(T0-TPED)*((1.0 - (RHO(1)))**ALPHAT1)
	XNC(1) = XNPED + (XN0-XNPED)*((1.0 - (RHO(1)**2))**ALPHAN)
	XNEL = XNC(1)*DELC
	DO 35 I = 2,55
	RHO(I) = RHO(I-1) + DELC
	TE(I) = TPED + A2*(T0-TPED)*((1.0 - (RHO(I)**2))**ALPHAT2) +
	2			   A1*(T0-TPED)*((1.0 - (RHO(I)))**ALPHAT1)
	XNC(I) = XNPED + (XN0-XNPED)*((1.0 - (RHO(I)**2))**ALPHAN)
	XNEL = XNEL + XNC(I)*DELC
35	CONTINUE
	TE(56) = TPED 
	XNC(56) = XNPED
	TE(57) = TBAR
	XNC(57) = XNBAR
C	AVG TEMP IN SOL = TSEP(1. - EXP(-DELN/DELT))
	TE(58) = TSOL*(1. - EXP(-1.*DELRATNT))
C	AVG DENSITY IN SOL =NSEP( 1. - EXP(-1.))= 0.632
	XNC(58) = 0.632*XNSOL
	XNEL = XNEL + XNC(56)*DELC
c	RELATING XNEL TO XNAV
c	XNEL = (1.333/(1.0+1.0/XNCTRPED))*XNAV
c	XNEL = XNAV 


	xnel = ratNlineav*xnav 

	JIT = 0
	DO 50 I = 1,55
	IF(TE(I).LE.0.) TE(I) = TPED
	IF(XNC(I).LE.0.) JIT = 1
	IF(XNC(I).LE.0.) XNC(I) = 1.E17
50	CONTINUE

	RETURN
	END