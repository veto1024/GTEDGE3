      SUBROUTINE TOROTATE(NMESH)
	INCLUDE 'SOLDIV.FI'

C	INTEGRATES TOROIDAL VELOCITY EQS (3/SPECIES) INWARD; ITERATES SOLUTION

C	SET UP EQUATIONS
	N = NMESH 
	AM = (aminor*SQRT(0.5*(1.+ELONG**2)))
	ep = am*rhor(n)/rmajor 
	sb = bphi/abs(bphi)
	xlanom01(n) = xlpm(n)-1./(am*rhor(n))
	xlanom02(n) = xlanom01(n)
 	

	
 
			
c	bfield = abs(bphi)
 
	
	eta401(n) = TI(N)/(ATNUM(1)*BFIELD)
	eta402(n) = TI(N)/(zbar2(n)*BFIELD)
	eta001(n) = eta01hat(n)
	eta002(n) = eta02hat(n)
	etanom01(n)= etanomhat1
	etanom02(n)= etanomhat1/zbar2(n)
	
     	AV(1,1) = (RMAJOR*VRAD1(N)+etanom01(n))
	AV(1,2) = 0.5*eta401(N)*sinion(n) 
	AV(1,3) = -0.5*eta401(n)*(COSION(N)+4.) 		
	AV(1,4) = 0.0
	AV(1,5) = 0.0
	AV(1,6) = 0.0
	AV(2,1) = eta401(n)*SINION(N)  +  2.*ep*etanom01(n) +
     1	AM*RHOR(N)*VRAD1(N)*(COSION(N)+3.)
	AV(2,2) = RMAJOR*VRAD1(N) + etanom01(n)
	AV(2,3) = 0.0
	AV(2,4) = 0.0
 	AV(2,5) = 0.0
	AV(2,6) = 0.0
	AV(3,1) = am*RHOR(N)*VRAD1(N)*SINION(N) -
     1		  eta401(n)*(COSION(N)+4.)
	AV(3,2) = 0.0
	AV(3,3) = RMAJOR*VRAD1(N) + etanom01(n)
	AV(3,4) = 0.0
	AV(3,5) = 0.0
	AV(3,6) = 0.0 
	AV(4,1) = 0.0
	AV(4,2) = 0.0
	AV(4,3) = 0.0
	AV(4,4) = (RMAJOR*VRAD2(N)+etanom02(n))
	AV(4,5) = 0.5*eta402(N)*sinimp(n) 
	AV(4,6) = -0.5*eta402(n)*(COSImp(N)+4.) 
	AV(5,1) = 0.0
	AV(5,2) = 0.0
	AV(5,3) = 0.0
	AV(5,4) = eta402(n)*SINImp(N)  +  2.*ep*etanom02(n) +
     1	AM*RHOR(N)*VRAD2(N)*(COSImp(N)+3.)
	AV(5,5) = RMAJOR*VRAD2(N) + etanom02(n)
    	AV(5,6) = 0.0
	AV(6,1) = 0.0
	AV(6,2) = 0.0
	AV(6,3) = 0.0
	AV(6,4) = am*RHOR(N)*VRAD2(N)*SINImp(N) -
     1		  eta402(n)*(COSImp(N)+4.)
      AV(6,5) = 0.0
	AV(6,6) = RMAJOR*VRAD2(N) + etanom02(n)

				
	BV(1,1) = VRAD1(N)*ep*(2.+cosion(n)) + 
     1  RMAJOR*(XNUC12(N)+XNUATI(N)) + AM*RHOR(N)*SORCOS(1,N)/YNI(N,1)-
     2	velthet1(n)*ep*(sinion(n) + vthtsin(1,n)) - 
     3	etanom01(n)*xlanom01(n)
	DNSDR =	-1.0*ep*(yni(n+1,1)*SINION(N+1)-yni(n,1)*SINION(N))/DELNA
	DNCDR =	-1.0*ep*(yni(n+1,1)*COSION(N+1)-yni(n,1)*COSION(N))/DELNA
	dscmult = 0.0
	dnsdr = dscmult*dnsdr
	dncdr = dscmult*dncdr
 	xlncm1(n) = DNCDR/(EP*YNI(N,1))
	xlnsm1(n) = DNsDR/(EP*YNI(N,1))
	BV(1,2) = -0.5*eta401(n)*(SINION(N)*(XLTIM(N)+xlnm(n))+
     1	xlnsm1(n)) + 0.5*RMAJOR*SORCOS(1,N)/YNI(N,1) -
	2	0.5*velthet1(n)*(vthtsin(1,n) + sinion(n))
	BV(1,3) = 0.5*eta401(n)*((XLTIM(N)+xlnm(n))*(COSION(N)+4.) -
	1	4./(AM*RHOR(N)) + xlncm1(n))	+
	2	0.5*velthet1(n)*(vthtcos(1,n) + cosion(n))
      BV(1,4) = -1.*RMAJOR*XNUC12(N)
	BV(1,5) = 0.0
	BV(1,6) = 0.0
	BV11(N) = BV(1,1)
	BV12(N) = BV(1,2)
	BV13(N) = BV(1,3)
	BV14(N) = BV(1,4)
	BV(2,1) = AM*RHOR(N)*XNUC12(N)*(2.+COSION(N)+COSIMP(N)) +
     2	AM*RHOR(N)*XNUATI(N)*(2.+COSION(N)+COSDENO(N))+
     3	RMAJOR*SORCOS(1,N)/YNI(N,1)+2.*vrad1(n) -
	4	2.*ep*etanom01(n)*xlanom01(n)
      BV(2,2) = RMAJOR*(XNUC12(N)+XNUATI(N)) 
      BV(2,2) = BV(2,2)+3.0*ETA01HAT(N)/(rmajor*(qsafe**2))
	bv(2,2) = bv(2,2) - etanom01(n)*xlanom01(n) 
      BV(2,3) = eta401(n)*(XLTIM(N)+XLNM(N))/ep +
	1			velthet1(n)/ep
	BV(2,4) = -1.*AM*RHOR(N)*XNUC12(N)*(2.+COSIMP(N)+COSION(N))
	BV(2,5) = -1.*RMAJOR*XNUC12(N)
	BV(2,6) = 0.0
	BV21(N) = BV(2,1)
	BV22(N) = BV(2,2)
	BV23(N) = BV(2,3)
	BV24(N) = BV(2,4)
	BV25(N) = BV(2,5)
	bv33(n) = bv(3,3)
	BV11(N) = BV(1,1)
	BV12(N) = BV(1,2)
	BV13(N) = BV(1,3)
	BV14(N) = BV(1,4)
	BV32(N) = BV(3,2)
	BV31(N) = BV(3,1)
	AV11(N) = AV(1,1)
	AV12(N) = AV(1,2)
	AV13(N) = AV(1,3)


	BV(3,1) = -2.*velthet1(n) +	  	
     1	AM*RHOR(N)*XNUC12(N)*(SINIMP(N)+SINION(N)) +
     2	AM*RHOR(N)*XNUATI(N)*(SINION(N)+SINDENO(N))+
     3	RMAJOR*SORSIN(1,N)/YNI(N,1)
	BV(3,2) = -1.0*eta401(n)*(XLTIM(N)+XLNM(N))/ep -
	1		velthet1(n)/ep
	BV(3,3) = RMAJOR*(XNUC12(N)+XNUATI(N))
    	BV(3,3) = BV(3,3) + 3.0*ETA01HAT(N)/(rmajor*(qsafe**2))	+
	1	vrad1(n) - etanom01(n)*xlanom01(n)
	BV(3,4) = -1.*AM*RHOR(N)*XNUC12(N)*(SINIMP(N)+SINION(N))
	BV(3,5) = 0.0
	BV(3,6) = -1.*RMAJOR*XNUC12(N)
	BV(4,1) = -1.*RMAJOR*XNUC21(N)
 	BV(4,2) = 0.0
	BV(4,3) = 0.0
	BV(4,4) = VRAD2(N)*ep*(2.+cosimp(n)) + 
     1  RMAJOR*(XNUC21(N)) + AM*RHOR(N)*SORCOS(2,N)/YNI(N,2)-
     2	velthet2(n)*ep*(sinimp(n) + vthtsin(2,n)) - 
     3	etanom02(n)*xlanom02(n)
	bv44(n) = bv(4,4)
 	DNSDR =	-1.*ep*(yni(n+1,2)*SINIMP(N+1)-yni(n,2)*SINIMP(N))/DELNA 
	DNCDR =	-1.*ep*(yni(n+1,2)*COSIMP(N+1)-yni(n,2)*COSIMP(N))/DELNA
	dnsdr = dscmult*dnsdr
	dncdr = dscmult*dncdr

	xlncm2(n) = DNCDR/(EP*YNI(N,2))
 	xlnsm2(n) = DNsDR/(EP*YNI(N,2))
	BV(4,5) = -0.5*eta402(n)*(SINImp(N)*(XLTIM(N)+xlnm(n))+
     1	xlnsm2(n)) + 0.5*RMAJOR*SORCOS(2,N)/YNI(N,2) -
     2	0.5*velthet2(n)*(vthtsin(2,n) + sinimp(n))
	BV(4,6) = 0.5*eta402(n)*((XLTIM(N)+xlnm(n))*(COSImp(N)+4.) -
	1	4./(AM*RHOR(N)) + xlncm2(n))	+
	2	0.5*velthet2(n)*(vthtcos(2,n) + cosimp(n))
 
	
	BV(5,1) = -1.*AM*RHOR(N)*XNUC21(N)*(2.+COSIMP(N)+COSION(N))
 	BV(5,2) = -1.*RMAJOR*XNUC21(N)
	BV(5,3) = 0.0
	BV(5,4) = AM*RHOR(N)*XNUC21(N)*(2.+COSION(N)+COSIMP(N)) +
     3	2.*vrad2(n) - 	2.*ep*etanom02(n)*xlanom02(n)
      BV(5,5) = RMAJOR*(XNUC21(N)) 
      BV(5,5) = BV(5,5)+3.0*ETA02HAT(N)/(rmajor*(qsafe**2))
	bv(5,5) = bv(5,5) - etanom02(n)*xlanom02(n) 
      BV(5,6) = eta402(n)*(XLTIM(N)+XLNM(N))/ep +
	1			velthet2(n)/ep

	BV(6,1) = -1.*AM*RHOR(N)*XNUC21(N)*(SINIMP(N)+SINION(N))
	BV(6,2) = 0.0
	BV(6,3) = -1.*RMAJOR*XNUC21(N)
	BV(6,4) = -2.*velthet2(n) +	  	
     1	AM*RHOR(N)*XNUC21(N)*(SINIMP(N)+SINION(N)) 
     	BV(6,5) = -1.0*eta402(n)*(XLTIM(N)+XLNM(N))/ep -
	1		velthet2(n)/ep
	BV(6,6) = RMAJOR*(XNUC21(N))
    	BV(6,6) = BV(6,6) + 3.0*ETA02HAT(N)/(rmajor*(qsafe**2))	+
	1	vrad2(n) - etanom02(n)*xlanom02(n)

	
 	SV(1) = (ATNUM(1)*EQ/XMAS(1))*(EPHIA+VRAD1(N)*BTHET) +
     1 	XMOMTOR1(N)/(YNI(N,1)*XMAS(1))
	SV(2) = EP*((ATNUM(1)*EQ/XMAS(1))*(EPHIA*(1.+COSION(N))
	1	+VRAD1(N)*BTHET*COSION(N))+ XMOMTOR1(N)/(YNI(N,1)*XMAS(1)))     
     	SV(2) = SV(2) + ETA01HAT(N)*(VTHTCOS(1,N)
	1		 - 2.)*velthet1(n)*(fp/(am*rhor(n)))/(rmajor)      
      SV(3) = EP*SINION(N)*(ATNUM(1)*EQ/XMAS(1))*
	1	(EPHIA+VRAD1(N)*BTHET) 
      SV(3) = SV(3) - ETA01HAT(N)*VTHTSIN(1,N)*velthet1(n)*
	1		(fp/(am*rhor(n)))/rmajor
	SV(4) = (ATNUM(2)*EQ/XMAS(2))*(EPHIA+VRAD2(N)*BTHET) +
     1 	XMOMTOR2(N)/(YNI(N,2)*XMAS(2))
	SV(5) = EP*((ATNUM(2)*EQ/XMAS(2))*(EPHIA*(1.+COSImp(N))
	1	+VRAD2(N)*BTHET*COSImp(N))+ XMOMTOR2(N)/(YNI(N,2)*XMAS(2)))     
     	SV(5) = SV(5) + ETA02HAT(N)*(VTHTCOS(2,N)
	1		 - 2.)*velthet2(n)*(fp/(am*rhor(n)))/(rmajor)      
      SV(6) = EP*SINImp(N)*(ATNUM(2)*EQ/XMAS(2))*
	1	(EPHIA+VRAD2(N)*BTHET) 
      SV(6) = SV(6) - ETA02HAT(N)*VTHTSIN(2,N)*velthet2(n)*
	1		(fp/(am*rhor(n)))/rmajor



     	
	SV1(N) = SV(1)
	SV2(N) = SV(2)
	SV3(N) = SV(3)
	SV4(N) = SV(4)
	SV5(N) = SV(5)
	SV6(N) = SV(6) 

c	do 150 j=1,6
c	sv(j) = -1.*sv(j)
c150	continue			 






	RETURN
	END