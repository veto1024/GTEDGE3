      SUBROUTINE DIVSOL(soln,solv,sole,Hmiller)
      include 'soldiv.fi'
		
c	  Solves 1D Sol-Div density, temp., flow & potential eqs
c	  5/6/08
	parameter (NCAP=51,jq=2,jp=6)
      dimension AA(NCAP),dpsi(NCAP),AALPHA(NCAP),ACOEF1(NCAP,NCAP),
     1		SORC(NCAP),QHEAT(NCAP),SRAD(NCAP),SAT(NCAP),SOL(NCAP),
     2		XT(NCAP),YN(NCAP),DVDPSI(NCAP),gpart(ncap),difp(ncap),
     3		xnupart(ncap),xnumom(ncap),DPRESS(NCAP) 
	real xmas1
      double precision TEDP,XLZ2DP,DLDT2DP,Tdbl,xlzdbl,dlzdbl
	xmas1 = 3.346e-27
c	  GEOMETRY
c		xlpar1 and xlpar2 are the parallel distances from top of plasma to divertor
c	    plate on the outboard and inboard, respectively.
c	   	NCAP-1 equal mesh intervals	of length dpsi(n) 
	open(500,file='bugging.txt',status='old')
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
	if(x1.lt.xl1) mone = n
	if(x2.lt.xl2) mtwo = (ncap-1)-n
10	continue
      write(*,*) mtwo
      delheat = deln
	delpart = deln
c	Particle and heat flux only after x-point
c      dQ/dr ~ -Qperp/delE
c      dGamma/dr~-Gammaperp/delN
	do 11 n = 1, ncap-1
	qheat(n) = fluxheat/delheat
	gpart(n) = fluxpart/delpart
	if(n.le.mone.or.n.ge.mtwo) then 
		qheat(n) = 0.0 
		gpart(n) = 0.0
	endif
11	continue
c		atomic physics frequencies & initial densities
	xnupart(1) = xnod*sviond - xnd*svrecd
	xnumom(1) = xnod*svatd + xnd*svrecd
	yn(1) = xnd
	do 15 n = 2,mone-1
	xnupart(n) = xnodiv*sviondiv - xndiv*svrecdiv
 	xnumom(n) = xnodiv*svatdiv + xndiv*svrecdiv
	yn(n) = xndiv
15	continue
	xnupart(mone) = xnosolxpt*sviondiv 
	xnumom(mone) = xnosolxpt*svatdiv 
	xnupart(mone+1) = xnosolxpt*svionsol 
	xnumom(mone+1) = xnosolxpt*svatsol 
	yn(mone) = xndiv
	yn(mone+1) = xnsol
	do 16 n = mone+2,mtwo-2
	xnupart(n) = xnosol*svionsol
	xnumom(n) = xnosol*svatsol
	yn(n) = xnsol
16	continue
	xnupart(mtwo-1) = xnosolxpt*svionsol 
 	xnumom(mtwo-1) = xnosolxpt*svatsol 
	xnupart(mtwo) = xnosolxpt*sviondiv 
 	xnumom(mtwo) = xnosolxpt*svatdiv
	yn(mtwo-1) = xnsol
	yn(mtwo) = xndiv
	do 17 n = mtwo+1,ncap-1
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
     1	 0.5*LOG((((1-X)**2)*(XMAS1/XMASEL))/(4.*XPI))  
	GAMSHEATHOUT = GAMSHEATHIN 
C		SETUP
      ACOEF1(2,2) = (2./7.)*XKAPPA/DPSI(2)
      ACOEF1(2,3) = -1.*ACOEF1(2,2)
	SORC(2) = 0.5*(QHEAT(2)-SRAD(2)-SAT(2))*(DPSI(1)+DPSI(2))
c	td used to be td1 here - solves eq(6) in PoP 16 042502
      SORC(2) = SORC(2) + 0.5*(QHEAT(1)-SRAD(1)-SAT(1))*DPSI(1) -
     1		  GAMSHEATHIN*XND1*SQRT(2./XMAS1)*((XK*td)**1.5)
      DO 25 N =3,NCAP-2
	ACOEF1(N,N-1) = (-2./7.)*XKAPPA/DPSI(N-1)
	ACOEF1(N,N+1) = (-2./7.)*XKAPPA/DPSI(N)
	ACOEF1(N,N)   = (2./7.)*XKAPPA*(1./DPSI(N-1) + 1./DPSI(N))
	SORC(N) = 0.5*(QHEAT(N)-SRAD(N)-SAT(N))*(DPSI(N)+DPSI(N-1))
25	CONTINUE
 	ACOEF1(ncap-1,ncap-1) = (2./7.)*XKAPPA/DPSI(ncap-2)
	ACOEF1(ncap-1,ncap-2) = -1.*ACOEF1(ncap-1,ncap-1)
	SORC(ncap-1) = 0.5*(QHEAT(ncap-1)-SRAD(ncap-1)-SAT(ncap-1))*
	1								(DPSI(ncap-1)+DPSI(ncap-2))
c	 td used to be td2 here
      SORC(ncap-1) = SORC(ncap-1) + 0.5*(QHEAT(ncap)-SRAD(ncap)-
	1 SAT(ncap))*DPSI(ncap-1)-GAMSHEATHOUT*XND2*SQRT(2./XMAS1)*
     2 ((xk*td)**1.5)
C		FORWARD ELIM		 
	AA(2) = ACOEF1(2,3)/ACOEF1(2,2)
	AALPHA(2) = SORC(2)/ACOEF1(2,2)
	DO 35 N = 3,NCAP-2
      AA(N) = ACOEF1(N,N+1)/(ACOEF1(N,N)-ACOEF1(N,N-1)*aa(N-1))
	AALPHA(N) = (SORC(N) - ACOEF1(N,N-1)*AALPHA(N-1))/
	1			(ACOEF1(N,N)-ACOEF1(N,N-1)*AA(N-1))
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
	DO 65 N = 1, NCAP-1
      XTRAN = XMAS1*(xnupart(n)+xnumom(n)+dvdpsi(n)+gpart(n)/yn(n))
      DIFP(N) = 2./XTRAN
65	CONTINUE		 
C		SET UP EQUATIONS
	ACOEF1(2,2) = DIFP(1)/DPSI(1) + 0.5*(DIFP(2)+DIFP(3))/DPSI(2) -
     1			 0.5*(DPSI(1)+DPSI(2))*XNUPART(2)/XT(2)
      ACOEF1(2,3) = -0.5*(DIFP(2)+DIFP(3))/DPSI(2)
	SORC(2) = 0.5*(XNUPART(1)*YN(1)+GPART(1))*DPSI(1) +
     1			   GPART(2)*0.5*(DPSI(1)+DPSI(2)) +
     2			   DIFP(1)*YN(1)*XK*XT(1)/DPSI(1)	 				
	DO 75 N = 3,NCAP-2
	ACOEF1(N,N-1) = -0.5*(DIFP(N)+DIFP(N-1))/DPSI(N-1)
	ACOEF1(N,N) = 0.5*((DIFP(N)+DIFP(N+1))/DPSI(N) + 
	1			(DIFP(N)+DIFP(N-1))/DPSI(N-1) -
     2			0.5*(XNUPART(N)/XT(N))*(DPSI(N)+DPSI(N-1)))
      ACOEF1(N,N+1) =  -0.5*(DIFP(N)+DIFP(N+1))/DPSI(N+1)
	SORC(N) = 0.5*GPART(N)*(DPSI(N)+DPSI(N-1))
75	CONTINUE
	N = NCAP 
	ACOEF1(N-1,N-2) = -0.5*(DIFP(N-1)+DIFP(N-2))/DPSI(N-2)
	ACOEF1(N-1,N-1) =  0.5*(DIFP(N-1)+DIFP(N-2))/DPSI(N-2) +
	1				  DIFP(N)/DPSI(N-1) -
     2				  0.5*(DPSI(N-1)+DPSI(N-2))*XNUPART(N-1)/XT(N-2)
	SORC(N-1) = 0.5*(XNUPART(N)*YN(N)+GPART(N))*DPSI(N-1) +
     1			0.5*GPART(N-1)*(DPSI(N-1)+DPSI(N-2)) +
     2			DIFP(N)*YN(N)*XK*XT(N)/DPSI(N-1)
      			
C		FORWARD ELIM
	AA(2) = ACOEF1(2,3)/ACOEF1(2,2)
	AALPHA(2) = SORC(2)/ACOEF1(2,2)
	DO 85 N = 3,NCAP-2
	AA(N) = ACOEF1(N,N+1)/(ACOEF1(N,N)-ACOEF1(N,N-1)*AA(N-1))
	AALPHA(N) = (SORC(N) - ACOEF1(N,N-1)*AALPHA(N-1))/
	1			(ACOEF1(N,N)-ACOEF1(N,N-1)*AA(N-1))
85	CONTINUE
C		BACKWARD SUB
	SOL(NCAP-1) = AALPHA(NCAP-1)
	DO 95 M = 1, NCAP-3
	N = (NCAP-1) - M
	SOL(N) = AALPHA(N) -AA(N)*SOL(N+1) 
95	CONTINUE
C		CALCULATE DENSITY FROM PRESSURE SOLUTION
	DO 100 N = 2,NCAP-2
      DPRESS(N) = SOL(N)
	YN(N) = SOL(N)/(XK*XT(N))
100	CONTINUE
		 			
      RETURN
	END
