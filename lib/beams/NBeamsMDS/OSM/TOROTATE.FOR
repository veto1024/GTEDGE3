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
 	
	etanomhat1 = 3.0e3
c	if(n.gt.20) etanomhat1 = 0.0

	etanomhat2 = etanomhat1/zbar2(n)
	
			
c	bfield = bphi
 
	
	eta401(n) = TI(N)/(ATNUM(1)*BFIELD)

	eta402(n) = TI(N)/(zbar2(n)*BFIELD)

	
	eta001(n) = eta01hat(n)
	eta002(n) = eta02hat(n)
	gq=1.0

	etanom01(n)= etanomhat1
	etanom02(n)= etanomhat2

c	include anomalous flux with radial flux
	vradanom(n) = 0.0 
	vrad1(n) = vrad1(n) + vradanom(n)
	vrad2(n) = vrad2(n) + vradanom(n) 	

	etanom01(n)= xnudraganomal(n)
	etanom02(n)= etanom01(n)/ZBAR2(N)
	
     	AV(1,1,n) = (RMAJOR*VRAD1(N))
	AV(1,2,n) = 0.5*eta401(N)*sinion(n) 
	AV(1,3,n) = -0.5*eta401(n)*(COSION(N)+3.) 		
	AV(1,4,n) = 0.0
	AV(1,5,n) = 0.0
	AV(1,6,n) = 0.0
	AV(2,1,n) = eta401(n)*SINION(N)  +  
     1	AM*RHOR(N)*VRAD1(N)*(COSION(N)+3.)
	AV(2,2,n) = RMAJOR*VRAD1(N) 
	AV(2,3,n) = 0.0
	AV(2,4,n) = 0.0
 	AV(2,5,n) = 0.0
	AV(2,6,n) = 0.0
	AV(3,1,n) = am*RHOR(N)*VRAD1(N)*SINION(N) -
     1		  eta401(n)*(COSION(N)+3.)
	AV(3,2,n) = 0.0
	AV(3,3,n) = RMAJOR*VRAD1(N) 
	AV(3,4,n) = 0.0
	AV(3,5,n) = 0.0
	AV(3,6,n) = 0.0 
	AV(4,1,n) = 0.0
	AV(4,2,n) = 0.0
	AV(4,3,n) = 0.0
	AV(4,4,n) = (RMAJOR*VRAD2(N))
	AV(4,5,n) = 0.5*eta402(N)*sinimp(n) 
	AV(4,6,n) = -0.5*eta402(n)*(COSImp(N)+3.) 
	AV(5,1,n) = 0.0
	AV(5,2,n) = 0.0
	AV(5,3,n) = 0.0
	AV(5,4,n) = eta402(n)*SINImp(N)  +  
     1	AM*RHOR(N)*VRAD2(N)*(COSImp(N)+3.)
	AV(5,5,n) = RMAJOR*VRAD2(N) 
    	AV(5,6,n) = 0.0
	AV(6,1,n) = 0.0
	AV(6,2,n) = 0.0
	AV(6,3,n) = 0.0
	AV(6,4,n) = am*RHOR(N)*VRAD2(N)*SINImp(N) -
     1		  eta402(n)*(COSImp(N)+3.)
      AV(6,5,n) = 0.0
	AV(6,6,n) = RMAJOR*VRAD2(N) 

				
	BV(1,1,n) = VRAD1(N)*ep*(2.+cosion(n)) + 
     1  RMAJOR*(XNUC12(N)+XNUATI(N)) + AM*RHOR(N)*SORCOS(1,N)/YNI(N,1)-
     2	velthet1(n)*ep*(sinion(n) + vthtsin(1,n)) 
	DNSDR =	-1.0*ep*(yni(n+1,1)*SINION(N+1)-yni(n,1)*SINION(N))/DELNA
	DNCDR =	-1.0*ep*(yni(n+1,1)*COSION(N+1)-yni(n,1)*COSION(N))/DELNA
	dscmult = 0.0
	dnsdr = dscmult*dnsdr
	dncdr = dscmult*dncdr
 	xlncm1(n) = DNCDR/(EP*YNI(N,1))
	xlnsm1(n) = DNsDR/(EP*YNI(N,1))
	BV(1,2,n) = -0.5*eta401(n)*(SINION(N)*(XLTIM(N)+xlnm(n))+
     1	xlnsm1(n)) + 0.5*RMAJOR*SORCOS(1,N)/YNI(N,1) -
	2	velthet1(n)*(vthtsin(1,n) + sinion(n))+vrad1(n)
	BV(1,3,n) = 0.5*eta401(n)*((XLTIM(N)+xlnm(n))*(COSION(N)+3.) -
	1	4./(AM*RHOR(N)) + xlncm1(n))	+
	2	velthet1(n)*(2.+vthtcos(1,n) + cosion(n))
      BV(1,4,n) = -1.*RMAJOR*XNUC12(N)
	BV(1,5,n) = 0.0
	BV(1,6,n) = 0.0
	
	BV(2,1,n) = AM*RHOR(N)*XNUC12(N)*(2.+COSION(N)+COSIMP(N)) +
     2	AM*RHOR(N)*XNUATI(N)*(2.+COSION(N)+COSDENO(N))+
     3	RMAJOR*SORCOS(1,N)/YNI(N,1)+2.*vrad1(n) 
      BV(2,2,n) = RMAJOR*(XNUC12(N)+XNUATI(N)) 
      BV(2,2,n) = BV(2,2,n)+3.0*gq*ETA01HAT(N)/(rmajor*(qsafe**2))
	
      BV(2,3,n) = eta401(n)*(XLTIM(N)+XLNM(N))/ep +
	1			velthet1(n)/ep
	BV(2,4,n) = -1.*AM*RHOR(N)*XNUC12(N)*(2.+COSIMP(N)+COSION(N))
	BV(2,5,n) = -1.*RMAJOR*XNUC12(N)
	BV(2,6,n) = 0.0
	


	BV(3,1,n) = -1.*velthet1(n) +	  	
     1	AM*RHOR(N)*XNUC12(N)*(SINIMP(N)+SINION(N)) +
     2	AM*RHOR(N)*XNUATI(N)*(SINION(N)+SINDENO(N))+
     3	RMAJOR*SORSIN(1,N)/YNI(N,1)
	BV(3,2,n) = -1.0*eta401(n)*(XLTIM(N)+XLNM(N))/ep -
	1		velthet1(n)/ep
	BV(3,3,n) = RMAJOR*(XNUC12(N)+XNUATI(N))
    	BV(3,3,n) = BV(3,3,n) + 3.0*gq*ETA01HAT(N)/(rmajor*(qsafe**2))	+
	1	vrad1(n)
	BV(3,4,n) = -1.*AM*RHOR(N)*XNUC12(N)*(SINIMP(N)+SINION(N))
	BV(3,5,n) = 0.0
	BV(3,6,n) = -1.*RMAJOR*XNUC12(N)
	BV(4,1,n) = -1.*RMAJOR*XNUC21(N)
 	BV(4,2,n) = 0.0
	BV(4,3,n) = 0.0
	BV(4,4,n) = VRAD2(N)*ep*(2.+cosimp(n)) + 
     1  RMAJOR*(XNUC21(N)) + AM*RHOR(N)*SORCOS(2,N)/YNI(N,2)-
     2	velthet2(n)*ep*(sinimp(n) + vthtsin(2,n)) 
	
 	DNSDR =	-1.*ep*(yni(n+1,2)*SINIMP(N+1)-yni(n,2)*SINIMP(N))/DELNA 
	DNCDR =	-1.*ep*(yni(n+1,2)*COSIMP(N+1)-yni(n,2)*COSIMP(N))/DELNA
	dnsdr = dscmult*dnsdr
	dncdr = dscmult*dncdr

	xlncm2(n) = DNCDR/(EP*YNI(N,2))
 	xlnsm2(n) = DNsDR/(EP*YNI(N,2))
	BV(4,5,n) = -0.5*eta402(n)*(SINImp(N)*(XLTIM(N)+xlnm(n))+
     1	xlnsm2(n)) + 0.5*RMAJOR*SORCOS(2,N)/YNI(N,2) -
     2	velthet2(n)*(vthtsin(2,n) + sinimp(n))+Vrad2(n)
	BV(4,6,n) = 0.5*eta402(n)*((XLTIM(N)+xlnm(n))*(COSImp(N)+3.) -
	1	4./(AM*RHOR(N)) + xlncm2(n))	+
	2	velthet2(n)*(2.+vthtcos(2,n) + cosimp(n))
 
	
	BV(5,1,n) = -1.*AM*RHOR(N)*XNUC21(N)*(2.+COSIMP(N)+COSION(N))
 	BV(5,2,n) = -1.*RMAJOR*XNUC21(N)
	BV(5,3,n) = 0.0
	BV(5,4,n) = AM*RHOR(N)*XNUC21(N)*(2.+COSION(N)+COSIMP(N)) +
     3	2.*vrad2(n) 
      BV(5,5,n) = RMAJOR*(XNUC21(N)) 
      BV(5,5,n) = BV(5,5,n)+3.0*gq*ETA02HAT(N)/(rmajor*(qsafe**2))
	 
      BV(5,6,n) = eta402(n)*(XLTIM(N)+XLNM(N))/ep +
	1			velthet2(n)/ep

	BV(6,1,n) = -1.*AM*RHOR(N)*XNUC21(N)*(SINIMP(N)+SINION(N))
	BV(6,2,n) = 0.0
	BV(6,3,n) = -1.*RMAJOR*XNUC21(N)
	BV(6,4,n) = -1.*velthet2(n) +	RMAJOR*SORSIN(2,N)/YNI(N,2) + 	
     1	AM*RHOR(N)*XNUC21(N)*(SINIMP(N)+SINION(N)) 
     	BV(6,5,n) = -1.0*eta402(n)*(XLTIM(N)+XLNM(N))/ep -
	1		velthet2(n)/ep
	BV(6,6,n) = RMAJOR*(XNUC21(N))
    	BV(6,6,n) = BV(6,6,n) + 3.0*gq*ETA02HAT(N)/(rmajor*(qsafe**2))	+
	1	vrad2(n) 

	 
 	SV(1,n) = (ATNUM(1)*EQ/XMAS(1))*(EPHIA+VRAD1(N)*BTHET) +
     1 	XMOMTOR1(N)/(YNI(N,1)*XMAS(1))
	SV(2,n) = EP*((ATNUM(1)*EQ/XMAS(1))*(EPHIA*(1.+COSION(N))
	1	+VRAD1(N)*BTHET*COSION(N))+ XMOMTOR1(N)/(YNI(N,1)*XMAS(1)))     
     	SV(2,n) = SV(2,n) + gq*ETA01HAT(N)*(VTHTCOS(1,N)
	1		 - 2.)*velthet1(n)*(fp/(am*rhor(n)))/(rmajor)      
      SV(3,n) = EP*SINION(N)*(ATNUM(1)*EQ/XMAS(1))*
	1	(EPHIA+VRAD1(N)*BTHET) 
      SV(3,n) = SV(3,n) + gq*ETA01HAT(N)*VTHTSIN(1,N)*velthet1(n)*
	1		(fp/(am*rhor(n)))/rmajor
	SV(4,n) = (ATNUM(2)*EQ/XMAS(2))*(EPHIA+VRAD2(N)*BTHET) +
     1 	XMOMTOR2(N)/(YNI(N,2)*XMAS(2))
	SV(5,n) = EP*((ATNUM(2)*EQ/XMAS(2))*(EPHIA*(1.+COSImp(N))
	1	+VRAD2(N)*BTHET*COSImp(N))+ XMOMTOR2(N)/(YNI(N,2)*XMAS(2)))     
     	SV(5,n) = SV(5,n) + gq*ETA02HAT(N)*(VTHTCOS(2,N)
	1		 - 2.)*velthet2(n)*(fp/(am*rhor(n)))/(rmajor)      
      SV(6,n) = EP*SINImp(N)*(ATNUM(2)*EQ/XMAS(2))*
	1	(EPHIA+VRAD2(N)*BTHET) 
      SV(6,n) = SV(6,n) + gq*ETA02HAT(N)*VTHTSIN(2,N)*velthet2(n)*
	1		(fp/(am*rhor(n)))/rmajor

	

c	do 150 j=1,6
c	sv(j) = -1.*sv(j)
c150	continue			 

	vrad1(n) = vrad1(n) - vradanom(n)
	vrad2(n) = vrad2(n) - vradanom(n) 




	RETURN
	END