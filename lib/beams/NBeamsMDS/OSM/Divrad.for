      SUBROUTINE DIVRAD
	INCLUDE 'SOLDIV.FI'
	PARAMETER (MM=5)
	DIMENSION TSR(MM),TDR(MM),XNR(MM),RADD(MM)
C		CALCULATES IMPURITY RADIATION COOLING IN SOL & DIVERTOR
C		USES LACKNER/POST MODEL
      DOUBLE PRECISION TEV, CELZ, DLZZDT
C	CALCULATES LACKNER INTEGRAL AT 5 RADIAL LOCATIONS & INTEGRATES 
C	RADIALLY.  EXP T, CHOPPED EXP N DISTRIBUTIONS
	DELT = DELN/DELRATNT 
	TSR(1) = TSEP
	TSR(2) = TSEP*EXP(-0.25*DELN/DELT)
	TSR(3) = TSEP*EXP(-0.50*DELN/DELT)
	TSR(4) = TSEP*EXP(-0.75*DELN/DELT)
	TSR(5) = TSEP*EXP(-1.00*DELN/DELT) 
	TDR(1) = TD
 	TDR(2) = TD*EXP(-0.25*DELN/DELT)
	TDR(3) = TD*EXP(-0.50*DELN/DELT)
	TDR(4) = TD*EXP(-0.75*DELN/DELT)
	TDR(5) = TD*EXP(-1.00*DELN/DELT) 
	XNR(1) = XNSEP
	XNR(2) = XNSEP*EXP(-0.25)
	XNR(3) = XNSEP*EXP(-0.50)
	XNR(4) = XNSEP*EXP(-0.75)
	XNR(5) = XNSEP*EXP(-1.00)
 
	QZMULT = 1.0
	RADDIV = 0.0
	RAD1 = 0.
	
	dx = xlpar1/50.
	xloc = 0.
	DELTAN = (XNSEP-XND)/50.
	DO 50 J=1,5 
	TDD = TDR(J)
	IF(TDR(J).GE.TSR(J)) TDD = TSR(J) - 50.
	DELTA = (TSR(J) - TDD)/50.
	T = TDD + 0.5*DELTA
 	XNA = XND + 0.5*DELTAN
	
	DQRADD = 0.
	DO 25 N = 1, 49
	xloc = xloc + dx
	if(xloc.le.dellt1) fon = xnod/xnd
	if(xloc.gt.dellt1.and.xloc.le.(delx1+dellt1)) fon = xnodiv/xndiv
	if(xloc.gt.(delx1+dellt1)) fon = xnosol/xnsol
	TEV = T
	stev = t
	IZ1 = IZINJECT
	IZ2 = IZINTRIN

	IF(IOPTRAD.EQ.1) GOTO 10
	CELZ = 0.
	if(iz1.ne.4.and.iz1.ne.6.and.iz1.ne.74) goto 4
c*************temp avoid fo affect on radiation*******
c	 goto 4
	
	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz1,sTEv,fon,scelz,sdlzzdt,ZAV)
	DQRAD = 1.E-13*sCELZ*FZINJECT*qzmultdiv
	IF(J.EQ.1) 
  	2  RADDIV = RADDIV + 1.E-13*(3.5*sCELZ/(T)-sDLZZDT)*qzmultdiv*
     3  FZINJECT*XNA/XK
	IF(J.EQ.1) 
     2  RAD1 = RAD1 + 1.E-13*(3.5*sCELZ/(T))*qzmultdiv*FZINJECT*XNA/XK
	goto 5 

4	CALL cefits (IZ1, TEV, celz, 1, dlzZdt)
	DQRAD = 1.E-13*CELZ*FZINJECT*qzmultdiv
	IF(J.EQ.1) 
	2  RADDIV = RADDIV + 1.E-13*(3.5*CELZ/(T)-DLZZDT)*qzmultdiv*
     3  FZINJECT*XNA/XK
	IF(J.EQ.1) 
     2  RAD1 = RAD1 + 1.E-13*(3.5*CELZ/(T))*qzmultdiv*FZINJECT*XNA/XK
 
	CELZ = 0.
	
	
5	if(iz2.ne.4.and.iz2.ne.6.and.iz2.ne.74) goto 6
c*************temp avoid fo affect on radiation*******
c	 goto 6 
	
 	if(fon.lt.1.e-5) fon = 1.e-5
	CALL CXRCEFITS(iz2,sTEv,fon,scelz,sdlzzdt,ZAV) 
		DQRAD = DQRAD + 1.E-13*sCELZ*FZINTRIN*QZMULTDIV
	IF(J.EQ.1) 
	2  RADDIV = RADDIV + 1.E-13*(3.5*sCELZ/(T)-sDLZZDT)*qzmultdiv*
     3  FZINTRIN*XNA/XK
	IF(J.EQ.1) 
     2  RAD1 = RAD1 + 1.E-13*(3.5*sCELZ/(T))*qzmultdiv*FZINTRIN*XNA/XK 
	goto 7

6	CALL CEFITS (IZ2, TEV, celz, 1, dlZzdt)
	DQRAD = DQRAD + 1.E-13*CELZ*FZINTRIN*QZMULTDIV
	IF(J.EQ.1) 
	2  RADDIV = RADDIV + 1.E-13*(3.5*CELZ/(T)-DLZZDT)*qzmultdiv*
     3  FZINTRIN*XNA/XK
	IF(J.EQ.1) 
     2  RAD1 = RAD1 + 1.E-13*(3.5*CELZ/(T))*qzmultdiv*FZINTRIN*XNA/XK
 
 7	continue 
	GOTO 15

10	DQRAD = 0.
	SCELZ = 0. 
	CALL NCEFITS (IZ1,T,TAURES,SCELZ,ZAV) 
	DQRAD = 1.E-13*SCELZ*FZINJECT
	SCELZ = 0.
	CALL NCEFITS (IZ2,T,TAURES,SCELZ,ZAV)
C	CALL CXRCEFITS(IZ2,T,FON,SCELZ,DLZZDT1,ZAV)
	DQRAD = DQRAD + 1.E-13*SCELZ*FZINTRIN*QZMULTDIV 
15	DQRADD = DQRADD + SQRT(T)*DQRAD*DELTA
	T = T + DELTA
	XNA = XNA + DELTAN
25	CONTINUE
	IF(J.EQ.1) RADDIV = RADDIV*XLPAR/50. 
	IF(J.EQ.1) RAD1 = RAD1*XLPAR/50. 

	Y = 25.3 - 1.15*LOG10(1.E-6*XNR(J)) + 2.3*LOG10(TSR(J))
	IF(TSR(J).LT.50.) Y=23.4-1.15*LOG10(1.E-6*XNR(J))+
	2					3.45*LOG10(TSR(J))
	ZEFF=(1.+(IZINJECT**2)*FZINJECT+16.*FBE+4.*FHE+(IZINTRIN**2)*
     2  FZINTRIN)/
     3	(1.+ IZINJECT*FZINJECT + 4.*FBE + 2.*FHE + IZINTRIN*FZINTRIN) 
	REDKAPPA2 = 2.0
      C5 = Y*ZEFF/(3.07E3*REDKAPPA2) 
	RADD(J) = XNR(J)*TSR(J)*SQRT(DQRADD/C5)
50	CONTINUE

	DQRAD = (DELN/8.)*RADD(1)
	DO 100 J=2,4
	DQRAD = DQRAD + (DELN/4.)*RADD(J)
100	CONTINUE
	DQRAD = DQRAD + (DELN/8.)*RADD(5)
	
	goto 1000
c	simple divertor radiation model (12/7/01)
	dqradz = 0.
	dqradbrem = 0. 
	IZ1 = IZINJECT
	IZ2 = IZINTRIN
 
	do 200 j=1,3
	if(j.eq.1) then
	stev = tsep
	tev = tsep
	xnm = xnsep
	fon = xnosol/xnsep
	delz = xlperp
	dep = epdiv
	endif
	if(j.eq.2) then
	stev = tdiv
	tev = tdiv
	xnm = xndiv
	fon = xnodiv/xndiv
	delz = delx
	dep = 2.*epdiv
	endif
	if(j.eq.3) then
	stev = td
	tev = td
	xnm = xnd
	fon = xnod/xnd
	delz = dellt
	dep = 2.*epdiv
	endif
	if(fon.lt.1.e-5) fon = 1.e-5
c	injected impurity
	if(iz1.ne.4.and.iz1.ne.6.and.iz1.ne.74) goto 125
	CALL CXRCEFITS(iz1,sTEv,fon,scelz,sdlzzdt,ZAV)
	DQRADz = dqradz + 
	2	1.E-13*sCELZ*FZINJECT*qzmultdiv*delz*dep*(xnm**2)
	
	goto 150 

125	CALL cefits (IZ1, TEV, celz, 1, dlzZdt)
	DQRADz = dqradz + 
	2	1.E-13*CELZ*FZINJECT*qzmultdiv*delz*dep*(xnm**2)
	
c	intrinsic impurity	
150	if(iz2.ne.4.and.iz2.ne.6.and.iz2.ne.74) goto 175
	CALL CXRCEFITS(iz2,sTEv,fon,scelz,sdlzzdt,ZAV) 
		DQRADz = DQRADz + 
     2		1.E-13*sCELZ*FZINTRIN*QZMULTDIV*delz*dep*(xnm**2)
	
      goto 190
	   
175	CALL CEFITS (IZ2, TEV, celz, 1, dlZzdt)
		DQRADz = DQRADz + 
     2		1.E-13*CELZ*FZINTRIN*QZMULTDIV*delz*dep*(xnm**2)

190	continue 
c	bremsstrahlung
	dqradbrem = dqradbrem +	(xnm**2)*sqrt(stev)*delz*dep
200	continue

	dqradz = 1.e-13*deln*dqradz
	dqradbrem = 4.8e-37*deln*zeff*dqradbrem 
	dqrad = dqradbrem + dqradz 
1000  RETURN
      END

	