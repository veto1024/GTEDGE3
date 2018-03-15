	subroutine febs
	include "soldiv.fi"
	parameter (lm = 26)
	dimension alit(lm,lm),abig(lm),alph(lm),sor(lm),sornb(lm),
     1			denom(lm)	 
c	forward elimination, backwards substitution solution
c	of diffusion equation for ion density in edge pedestal
C	ION FLUX IN PEDESTAL DUE TO NEUTRAL BEAM PARTICLES
	areaped = 39.44*RMAJOR*AMINOR*SQRT((1.+ELONG**2)/2.)
	delma = delna
 	gamnb =  ((pbeam*1.e6)/(xk*eb*1.e3))/areaped
	gamin = gamion(1,1)
	
c	"LITTLE A" COEFFICIENTS

c		first mp (use values at mp 1 for mp 0)
	ff = atnum(1) + fracz*atnum(2) 
	alit(1,1) = (0.5*(diffA(1)+diffA(2))/(delma**2))+0.5*coefp(2)/delma
	1			+ 2.*diffA(1)/(delma**2)
	     				 
	alit(1,2) = -0.5*(diffA(1)+diffA(2))/(delma**2) +0.5*coefp(1)/delma

	sor(1) = (gamin/delma) + 0.5*(xnuioni(1)+xnuioni(2))*ff*yni(1,1)+
     1		 dens(1)*xnuionb(1)	
	sornb(1) = dens(1)*xnuionb(1) 
c		interior mp
	do 100 n = 2,23
	alit(n,n) = (diffA(n)+0.5*(diffA(n-1)+diffA(n+1)))/(delma**2) +
	1			0.5*(coefp(n+1)-coefp(n-1))/delma
    
      alit(n,n-1)=-0.5*((diffA(n-1)+diffA(n))/(delma**2)-coefp(n)/delma)
	alit(n,n+1)=-0.5*((diffA(n+1)+diffA(n))/(delma**2)+coefp(n)/delma)
     
	sor(n) = 0.5*(dens(n)*xnuionb(n)+dens(n-1)*xnuionb(n-1))   +
	1		 0.5*(xnuioni(n)+xnuioni(n-1))*ff*yni(n,1)
	sornb(n) = 0.5*(dens(n)*xnuionb(n)+dens(n-1)*xnuionb(n-1)) 
100	continue
c		last mp
	alit(24,24) = (0.5*(diffA(23)+diffA(24))+2.*diffA(25))/(delma**2) 
c     1			(coefp(25)-coefp(24)-0.5*coefp(23))/delma
     1		+	 0.5*coefp(24)/delma
      alit(24,23) =-0.5*((diffA(23)+diffA(24))/(delma**2)+coefp(24)/
     1	delma)
	sor(24) = (2.*diffA(25)-coefp(24)*delma)*yni(25,1)/(delma**2) + 
	1			0.5*(dens(23)*xnuionb(23)+dens(24)*xnuionb(24))	+
	2			0.5*(xnuioni(23)+xnuioni(24))*ff*yni(24,1)

	sornb(24) = 0.5*(dens(23)*xnuionb(23)+dens(24)*xnuionb(24))
C      FORWARD ELIMINATION		
	abig(1) = alit(1,2)/alit(1,1)
	alph(1) = sor(1)/alit(1,1)
	do 200 n = 2, 24
	denom(n) = (alit(n,n)-alit(n,n-1)*abig(n-1))
	abig(n) = alit(n,n+1)/denom(n)
	alph(n) = (sor(n)-alit(n,n-1)*alph(n-1))/denom(n)
200	continue
C	 BACKWARDS SUBSTITUTION
	yni(24,1) = alph(24)
	do 300 nn = 1,23
	n = 24-nn 
	yni(n,1) = alph(n) - abig(n)*yni(n+1,1)
300	continue
	y = vpinchi(1)
	return
	end