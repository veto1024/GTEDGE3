	SUBROUTINE DIVSOL2
c	INCLUDE 'SOLDIV.FI'

c	  Solves 1D Sol-Div density, temp., flow & potential eqs
c	  5/6/08
	  parameter (NCAP=51,jq=2,jp=6)
      dimension AA(NCAP),dpsi(NCAP),AALPHA(NCAP),ACOEF(NCAP,NCAP),
     1		SORC(NCAP),QHEAT(NCAP),SRAD(NCAP),SAT(NCAP),SOL(NCAP),
     2		XT(NCAP),YN(NCAP),DVDPSI(NCAP),gpart(ncap),difp(ncap),
     4		xnupart(ncap),xnumom(ncap),DPRESS(NCAP) 


      double precision TEDP,XLZ2DP,DLDT2DP,Tdbl,xlzdbl,dlzdbl

c	  GEOMETRY
c		xlpar1 and xlpar2 are the parallel distances from top of plasma to divertor
c	    plate on the outboard and inboard, respectively.
c	   	NCAP-1 equal mesh intervals	of length dpsi(n)
	do 5, n = 1, ncap-1
      dpsi(n) = (xlpar1+xlpar2)/(ncap-1.)
5	continue
C		cells with heat & part fluxes from plasma
	xl1 = (divdepth1/betag)
	xl2 = (divdepth2/betag)
	x1 = 0.
	x2 = 0.
	do 10 n = 1 , ncap-1
	x1 = x1 + dpsi(n)
	x2 = x2 + dpsi (ncap-n)
	if(x1.lt.xl1) m1 = n
	if(x2.lt.xl2) m2 = ncap-n
10	continue 
	delheat = deln
	delpart = deln 
	do 11 n = 1, ncap-1
	qheat(n) = fluxheat/delheat
	gpart(n) = fluxpat/delpart
	if(n.lt.m1.or.n.gt.m2) then 
		qheat(n) = 0.0 
		gpart(n) = 0.0
	endif
11	continue
c		atomic physics frequencies & initial densities
	xnupart(1) = xnod*sviond - xnd*svrecd
	xnumom(1) = xnod*svatd + xnd*svrecd
	yn(1) = xnd
	do 15 n = 2,m1-1
	xnupart(n) = xnodiv*sviondiv - xndiv*svrecdiv
 	xnumom(n) = xnodiv*svatdiv + xndiv*svrecdiv
	yn(n) = xndiv
15	continue
	xnupart(m1) = xnosolxpt*sviondiv 
	xnumom(m1) = xnosolxpt*svatdiv 
	xnupart(m1+1) = xnosolxpt*svionsol 
	xnumom(m1+1) = xnosolxpt*svatsol 
	yn(m1) = xndiv
	yn(m1+1) = xnsol
	do 16 n = m1+1,ncap-1-m2+1
	xnupart(n) = xnosol*svionsol
	xnumom(n) = xnosol*svatsol
	yn(n) = xnsol
16	continue
	xnupart(ncap-1-m2-1) = xnosolxpt*svionsol 
 	xnumom(ncap-1-m2-1) = xnosolxpt*svatsol 
	xnupart(ncap-1-m2) = xnosolxpt*sviondiv 
 	xnumom(ncap-1-m2) = xnosolxpt*svatd
	yn(ncap-1-m2-1) = xnsol
	yn(ncap-1-ms) = xndiv
	do 17 n = ncap-m2,ncap-2
	xnupart(n) = xnodiv*sviondiv - xndiv*svrecdiv
  	xnumom(n) = xnodiv*svatdiv + xndiv*svrecdiv
	yn(n) = xndiv
17	continue
	xnupart(ncap-1) = xnod*sviond - xnd*svrecd
 	xnumom(ncap-1) = xnod*svatd + xnd*svrecd   
	yn(ncap-1) = xnd
C			COULOMB LOGARITHM & PARALLEL CONDUCTIVITY
	XNN = 0.5*(XNSEP+XND)
	XTT = 0.5*(TSEP+TD)	
 	Y = 25.3 - 1.15*LOG10(1.E-6*XNN) + 2.3*LOG10(XTT)
 	IF(XTT.LT.50.) Y = 23.4-1.15*LOG10(1.E-6*XNN)+3.45*LOG10(XTT)
	ZEFF=(1.+(IZINJECT**2)*FZINJECT + 4.*FHE +(IZINTRIN**2)*FZINTRIN)/
	2      (1.+ IZINJECT*FZINJECT + 2.*FHE + IZINTRIN*FZINTRIN) 
	XKAPPA = (3.07E4/(ZEFF*Y))
 
C	  FEBS SOLUTION OF HEAT BALANCE EQ ***********************

C		EQ 13.13 FUSION PLASMA PHYSICS
	XK = 1.6E-19 
	X = SECEL
	XPI = 3.14159
	XMASEL = 9.108E-31
C		assuming TI=TE at plate
	GAMSHEATHIN = 2. + 2./(1.-X) +
     1	 0.5*LOG((((1-X)**2)*(XMAS(1)/XMASEL))/(4.*XPI))  
	GAMSHEATHOUT = GAMSHEATHIN 
C		SETUP
	ACOEF(2,2) = (2./7.)*XKAPPA/DPSI(2)
	ACOEF(2,3) = -1.*ACOEF(2,2)
	SORC(2) = 0.5*(QHEAT(2)-SRAD(2)-SAT(2))*(DPSI(1)+DPSI(2))
	SORC(2) = SORC(2) + 0.5*(QHEAT(1)-SRAD(1)-SAT(1))*DPSI(1) -
     1		  GAMSHEATHIN*XND1*SQRT(2./XMAS(1))*((XK*TD1)**1.5)
      DO 25 N =3,NCAP-2
	ACOEF(N,N-1) = (-2./7.)*XKAPPA/DPSI(N-1)
	ACOEF(N,N+1) = (-2./7.)*XKAPPA/DPSI(N)
	ACOEF(N,N)   = (2./7.)*XKAPPA*(1./DPSI(N-1) + 1./DPSI(N))
	SORC(N) = 0.5*(QHEAT(N)-SRAD(N)-SAT(N))*(DPSI(N)+DPSI(N-1))
25	CONTINUE
 	ACOEF(N-1,N-1) = (2./7.)*XKAPPA/DPSI(N-2)
	ACOEF(N-1,N-2) = -1.*ACOEF(N-1,N-1)
	SORC(N-1) = 0.5*(QHEAT(N-1)-SRAD(N-1)-SAT(N-1))*
	1								(DPSI(N-1)+DPSI(N-2))
	SORC(N-1) = SORC(N-1) + 0.5*(QHEAT(N)-SRAD(N)-SAT(N))*DPSI(N-1) -
     1		  GAMSHEATHOUT*XND2*SQRT(2./XMAS(1))*((SK*TD2)**1.5)
C		FORWARD ELIM
	AA(2) = ACOEF(2,3)/ACOEF(2,2)
	AALPHA(2) = SORC(2)/ACOEF(2,2)
	DO 35 N = 3,NCAP-1
	AA(N) = ACOEF(N,N+1)/(ACOEF(N,N)-ACOEF(N,N-1)*AA(N-I))
	AALPHA(N) = (SORC(N) - ACOEF(N,N-1)*AALPHA(N-1))/
	1			(ACOEF(N,N)-ACOEF(N,N-1)*AA(N-I))
35	CONTINUE
C		BACKWARD SUB
	SOL(NCAP-1) = AALPHA(NCAP-1)
	DO 45 M = 1, NCAP-3
	N = (NCAP-1) - M
	SOL(N) = AALPHA(N) -AA(N)*SOL(N+1) 
45	CONTINUE
C		CALC TEMPERATURE FROM SOL = T**7/2
	DO 55 N = 2,NCAP-1
	XT(N) = (SOL(N))**(2./7.)
55	CONTINUE			     			

C	SOLVE FOR PRESSURE DISTRIBUTION****************************

C		SET UP DIFFUSION COEF
	DO 65 N = 1, NCAP
	XTRAN = XMAS(1)*(xnupart(n)+xnumom(n)+dvdpsi(n)+gpart(n)/yn(n))
      DIFP(N) = 2./XTRAN
65	CONTINUE		 
C		SET UP EQUATIONS
	ACOEF(2,2) = DIFP(1)/DPSI(1) + 0.5*(DIFP(2)+DIFP(3))/DPSI(2) -
     1			 0.5*(DPSI(1)+DPSI(2))*XNUPART(2)/XT(2)
      ACOEF(2,3) = -0.5*(DIFP(2)+DIFP(3))/DPSI(2)
	SORC(2) = 0.5*(XNUPART(1)*YN(1)+GPART(1))*DPSI(1) +
     1			   GPART(2)*0.5*(DPSI(1)+DPSI(2)) +
     2			   DIFP(1)*YN(1)*XK*XT(1)/DPSI(1)	 				
	DO 75 N = 3,NCAP-2
	ACOEF(N,N-1) = -0.5*(DIFP(N)+DIFP(N-1))/DPSI(N-1)
	ACOEF(N,N) = 0.5*((DIFP(N)+DIFP(N+1))/DPSI(N) + 
	1			(DIFP(N)+DIFP(N-1))/DPSI(N-1) -
     2			0.5*(XNUPART(N)/XT(N))*(DPSI(N)+DPSI(N-1)))
      ACOEF(N,N+1) =  -0.5*(DIFP(N)+DIFP(N+1))/DPSI(N+1)
	SORC(N) = 0.5*GPART(N)*(DPSI(N)+DPSI(N-1))
75	CONTINUE
	N = NCAP 
	ACOEF(N-1,N-2) = -0.5*(DIFP(N-1)+DIFP(N-2))/DPSI(N-2)
	ACOEF(N-1,N-1) =  0.5*(DIFP(N-1)+DIFP(N-2))/DPSI(N-2) +
	1				  DIFP(N)/DPSI(N-1) -
     2				  0.5*(DPSI(N-1)+DPSI(N-2))*XNUPART(N-1)/XT(N-1)
      SORC(N-1) = 0.5*(XNUPART(N)*YN(N)+GPART(N))*DPSI(N-1) +
     1			0.5*GPART(N-1)*(DPSI(N-1)+DPSI(N-2)) +
     2			DIFP(N)*YN(N)*XK*XT(N)/DPSI(N-1)			
C		FORWARD ELIM
	AA(2) = ACOEF(2,3)/ACOEF(2,2)
	AALPHA(2) = SORC(2)/ACOEF(2,2)
	DO 85 N = 3,NCAP-1
	AA(N) = ACOEF(N,N+1)/(ACOEF(N,N)-ACOEF(N,N-1)*AA(N-I))
	AALPHA(N) = (SORC(N) - ACOEF(N,N-1)*AALPHA(N-1))/
	1			(ACOEF(N,N)-ACOEF(N,N-1)*AA(N-I))
85	CONTINUE
C		BACKWARD SUB
	SOL(NCAP-1) = AALPHA(NCAP-1)
	DO 95 M = 1, NCAP-3
	N = (NCAP-1) - M
	SOL(N) = AALPHA(N) -AA(N)*SOL(N+1) 
95	CONTINUE
C		CALCULATE DENSITY FROM PRESSURE SOLUTION
	DO 100 N = 2,NCAP-1
	DPRESS(N) = SOL(N)
	YN(N) = SOL(N)/(XK*XT(N))
100	CONTINUE
 

	RETURN 
	END