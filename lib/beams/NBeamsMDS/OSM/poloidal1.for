	SUBROUTINE POLOIDAL(nmesh)
	INCLUDE 'soldiv.FI'
C  SOLVES THE 6 COUPLED EQUATIONS FOR ION & IMPURITY NORMALIZED
C  POLOIDAL VELOCITIES AND SIN & COS COMPONENTS OF DENSITY ASYM.
	PARAMETER (Jq=4,jk=2,jv=201,nm=26) 
	DIMENSION A(Jq), B(Jq), S(Jq),C(Jq,Jq),gam(jq),vrhat(jk),
	2	vthetsin(jk),vthetcos(jk),hnuat(jk),hnuion(jk),sc1(jk),
	3	sc2(jk),cc1(jk),cc2(jk),vc1(jk),vc2(jk),vc3(jk),vc4(jk),
	4	vs1(jk),vs2(jk),vs3(jk),vs4(jk),dela(jk),dels(jk),
	5	VPHICOS(JK),VPHISIN(JK),XV(JK),YV(JK),FRIC(JK),grad(jk),
	6	vt(jv),ss(jk)
298	format(1x,'OK')
301	format(7e11.4)
	
	nit = 0
	n = nmesh	
	iopteq = 2
	C9 = 0.0
	c8 = 0.0
	c7 = 0.0
	epol = 1.e-3
		
	xne = xni(1)*atnum(1)+xni(2)*atnum(2)  
	GAM(1) = XNi(1)*ATNUM(1)/(XNi(1)*ATNUM(1)+XNi(2)*ATNUM(2))
 	GAM(2) = XNi(2)*ATNUM(2)/(XNi(1)*ATNUM(1)+XNi(2)*ATNUM(2))
	if(ioptedge.ne.0) xne = xni(1)*atnum(1) + xni(2)*atnum(2)	 
		
	SIGNB = fp/abs(fp)
 	 
C	POLOIDAL ROTATION EQS
			
c	couple to toroidal rotation
c	if(it.gt.4) then 
c	vphia(1) = rmajor*omegt(n,1)	
c	vphia(2) = rmajor*omegt(n,4) 
c	endif
c	5/03 formulation
	
c	goto 1003
	if(ioptedge.eq.1) 	then
  	bfield = bphi
	x1 = atnum(1)*eq*bfield/(xmas(1)*VTH(1))
	y1 = qsafe*rmajor/vth(1)
	vr1 = gamion(nmesh,1)/yni(nmesh,1)
	vrad1(nmesh) = vr1

	if(vr1.gt.10.) vr1 = 10.
	vrhat(1) = x1*y1*vr1
	x2 = atnum(2)*eq*bfield/(xmas(2)*VTH(2))
 	y2 = qsafe*rmajor/vth(2)

	R2 = 0.01
c	sputtering yield for C in the 100-1000eV range is about .02,
c	allow for half othis to be trapped in divertor , R2 = 0.01
	vr2 = gamion(nmesh,2)/(yni(nmesh,2))

	vrad2(nmesh) = vr2
	vrhat(2) = x2*y2*vr2
	
 	endif
	 
	DO 100 NN=1,201
 	epd =0.0
	vt(nn) =0.0 
 	VTHET1 = VTHETa(1)
	VTHET2 = VTHETa(2)

c	************** 
	zc =1.
c	**************	
1003	nlimit = 2
	if(ioptedge.eq.3) nlimit = 1
c	Rutherford friction factor C1		
	cf = 0.35 
	cf = 1.0
1004	DO 1005 J=1,2
	 	K=2
	IF(J.EQ.2) K=1
	zcc = zc	

	A(J)= -1.0*FP*(QSAFE**2)*F(J)*(1.0+(zcc*(2./3.)*cosden(j)+COSPOT+
	1		c7*dela(j))) 
     2   +	cf*FP*XNUSTAR(J,K)*(1.+c9*fric(j)) - 
     3	QSAFE*ep*(VPHIa(J)/VTH(J))*(SINDEN(J) + SINPOT)
      B(J)=-1.*cf*FP*XNUSTAR(J,K)*SQRT(XMAS(J)/
     1	XMAS(K))*(1.+c9*fric(k))
	s(j) = -1.*vrhat(j) 
  	S(J) = s(j) - 1.0*QSAFE*ep*(((VPHIa(J)/VTH(J))**2)*SINPOT + 
	1	(VPHIa(J)/VTH(J))*((PRESS(J)/VTH(J))*SINPOT + 
     2	0.5*(VPHIa(J)/VTH(J))*SINDEN(J))) -
     3(QSAFE**2)*F(J)*FP*(((VPHIa(J)+PRESS(J))/VTH(J))*
     4	COSPOT+c7*dels(j))
	S(J) = S(J) - 0.25*QSAFE*EP*sinden(j) 
	s(j) = s(j)	- 0.25*QSAFE*EP*atnum(j)*(teP/temp(j))*
     2	((1.+cosden(j))*(gam(j)*sinden(j)+gam(k)*sinden(k))-
     3		 sinden(j)*(gam(j)*cosden(j)+gam(k)*cosden(k))	)

	nit = nit + 1

c	******************************************
	if(ioptedge.eq.1) then

		hnuat(1) = rminor*coldno(n)*svata(n)/VTH(1)
		avion = 0.5*(xnuioni(n)+ xnuioni(n+1))
		if(n.gt.1) avion = 0.5*(xnuioni(n)+xnuioni(n-1))		
		hnuion(1) =rminor*avion/VTH(1)
		hnuat(2) = 0.
		hnuion(2)=0.0

     	 	
c	*******************************************
	  a(j) = a(j) + qsafe*fp*hnuat(j)/ep
	  xy = (xne/xni(j))*(cosdene+cosdeno(j))-(cosden(j)+cosdeno(j)) 
       xz = (xne/xni(j))*(sindene+sindeno(j))-(sinden(j)+sindeno(j)) 
        yb = (1.+cosden(j))*xy+sinden(j)*xz
	  a(j) = a(j) + 0.5*fp*ep*hnuion(j)*yb
	  yx = (vphia(j)/VTH(j))*ep*(cosdene+cosdeno(j))
	  yz = 0.667*qsafe*f(j)*(sindene+sindeno(j))
	  s(j) = s(j) - (xne/xni(j))*hnuion(j)*ep*(yx-yz)
	endif
	
 	
c	write(6,301) a(j),b(j),s(j),f(j),xnustar(j,k),cs(1),cs(2)

1005	CONTINUE 

1006	continue
	  	
	
      VTHETa(1) = (A(2)*S(1)-B(1)*S(2))/(A(2)*A(1) - B(2)*B(1))
	if(ioptedge.ne.3)
     1	VTHETa(2) = (A(1)*S(2)-B(2)*S(1))/(A(1)*A(2) - B(1)*B(2))
	
1015	CONTINUE
	epd = 0.5
c	if(nmesh.eq.25.and.nn.gt.50) epd = 0.1
	IVTHET = IVTHETREAL
	if(nn.gt.1) then
 	VTHETa(1) = (1.-epd)*VTHETa(1)+epd*VTHET1
	VTHETa(2) = (1.-epd)*VTHETa(2)+epd*VTHET2
	endif
c	***************************************************
	xx = vtheta(1)*vth(1)*fp
c	if(abs(xx).gt.1.e4) vtheta(1) = 1.e4*vtheta(1)/abs(xx)
	yy = vtheta(2)*vth(2)*fp
c	if(abs(xx).gt.1.e4) vtheta(2) = 1.e4*vtheta(2)/abs(xx)
c	****************************************************
C	DENSITY ASYMMETRIES	
c	cosine moment species 1		
	C(1,1)=F(1)*(QSAFE**2)*FP*VTHETa(1)*((1./3.)-c7*vc1(1))+c7*vc2(1) 	
	1 + 0.5*(EP**2)*cf*FP*XNUSTAR(1,2)*VTHETa(2)*SQRT(XMAS(1)/XMAS(2))
	C(1,2) = -0.5*QSAFE*EP*((FP**2)*(VTHETa(1)**2)-0.5) +
	2  0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(1) + c7*vc3(1)
	C(1,3) = 0.5*FP*(EP**2)*cf*XNUSTAR(1,2)*VTHETa(1)
	C(1,4) = 0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(2)
	S(1) = -0.5*F(1)*(QSAFE**2)*FP*(VTHETa(1)*(1.+COSPOT+c8*cc1(1))-
	1	((VPHIa(1)+PRESS(1))/VTH(1))*COSPOT+c8*cc2(1))	
	s(1) = s(1) - QSAFE*EP*xv(1)*c9	+ c7*vc4(1)
	if(ioptedge.eq.1) then
	  c(1,1) = c(1,1)+0.5*ep*qsafe*hnuion(1)*fp*vtheta(1)
	  c(1,2) = c(1,2)-0.5*ep*qsafe*hnuat(1)*hnuion(1)
	  xy = fp*vtheta(1)*cosdeno(1)-hnuion(1)*sindeno(1)
	  xz = cosdeno(1)*(1.+xne/xni(1))+(xne/xni(1))*cosdene
	  yz = (xne/xni(1))*(sindene+sindeno(1))
	  s(1) = s(1)-0.5*ep*qsafe*hnuat(1)*xy-qsafe*hnuion(1)*fp*
     1         (0.5*vtheta(1)*xz +0.333*(qsafe/ep)*f(1)*yz)
      endif	   
c	sine moment species 1
	C(2,1) = 0.5*QSAFE*EP*((FP**2)*(VTHETa(1)**2)-0.5) -
	1	0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(1) +c7*vs3(1)	
	C(2,2) = (QSAFE**2)*F(1)*FP*VTHETa(1)*((1./3.)+vs1(1))+vs2(1) -
	1  0.5*(EP**2)*cf*FP*XNUSTAR(1,2)*VTHETa(2)*SQRT(XMAS(1)/XMAS(2))
	C(2,4) = 0.5*(EP**2)*cf*FP*XNUSTAR(1,2)*VTHETa(1)
 	C(2,3) = -0.25*QSAFE*ep*ATNUM(1)*(TEP/TEMP(1))*GAM(2)
	S(2) = -0.5*(QSAFE**2)*F(1)*FP*(VTHETa(1)*(sinpot+c8*sc1(1))-
	1	((VPHIa(1)+PRESS(1))/VTH(1))*SINPOT+c8*sc2(1))
	S(2) = S(2) - 0.5*QSAFE*EP*(((VPHIa(1)/VTH(1))**2)+c9*yv(1))
	2  - 0.5*QSAFE*(FP**2)*EP*(VTHETa(1)**2) + c7*vs4(1)
	if(ioptedge.eq.1) then
	  c(2,1) = c(2,1) + 0.5*ep*qsafe*hnuion(1)*hnuat(1)
	  c(2,2) = c(2,2) + 0.5*ep*qsafe*hnuion(1)*fp*vtheta(1)
	  xy = fp*vtheta(1)*sindeno(1)+hnuion(1)*cosdeno(1)
	  xz = sindeno(1)*(1.+xne/xni(1))+(xne/xni(1))*sindene
	  yz = (xne/xni(1))*(cosdene+cosdeno(1))
	  s(2) = s(2)-0.5*ep*qsafe*hnuat(1)*xy-qsafe*hnuion(1)*fp*
     1         (0.5*vtheta(1)*xz +0.333*(qsafe/ep)*f(1)*yz)
	endif
c	cosine moment species 2
	C(3,3) = F(2)*FP*(QSAFE**2)*VTHETa(2)*((1./3.)-vc1(2))+c7*vc2(2) -
	1  0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(1)*SQRT(XMAS(2)/XMAS(1))
	C(3,4) = -0.5*QSAFE*EP*((FP**2)*(VTHETa(2)**2)-0.5) +
	2	0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(2) +c7*vc3(2)			
	C(3,1) = 0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(2)
	C(3,2) = 0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(1)
	S(3) = -0.5*F(2)*(QSAFE**2)*FP*(VTHETa(2)*(1.+COSPOT+c8*cc1(2))-
	1	((VPHIa(2)+PRESS(2))/VTH(2))*COSPOT+c8*cc2(2))	
	s(3) = s(3) - QSAFE*EP*xv(2)*c9	 +c7*vc4(2)
	if(ioptedge.eq.1) then
	  c(3,3) = c(3,3)+0.5*ep*qsafe*hnuion(2)*fp*vtheta(2)
	  c(3,4) = c(3,4)-0.5*ep*qsafe*hnuat(2)*hnuion(2)
	  xy = fp*vtheta(2)*cosdeno(2)-hnuion(2)*sindeno(2)
	  xz = cosdeno(2)*(1.+xne/xni(2))+(xne/xni(2))*cosdene
	  yz = (xne/xni(2))*(sindene+sindeno(2))
	  s(3) = s(3)-0.5*ep*qsafe*hnuat(2)*xy-qsafe*hnuion(2)*fp*
     1         (0.5*vtheta(2)*xz +0.333*(qsafe/ep)*f(2)*yz)
      endif	   

c	sine moment species 2
	C(4,3) = 0.5*QSAFE*EP*((FP**2)*(VTHETa(2)**2)-0.5) -
	1	0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(2) +c7*vs3(2)
	C(4,4) = (QSAFE**2)*F(2)*FP*VTHETa(2)*((1./3.)+vs1(2))+c7*vs2(2) -
	1  0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(1)*SQRT(XMAS(2)/XMAS(1))
	C(4,2) = 0.5*(EP**2)*cf*FP*XNUSTAR(2,1)*VTHETa(2)
 	C(4,1) = -0.25*QSAFE*ep*ATNUM(2)*(TEP/TEMP(2))*GAM(1)
	S(4) = -0.5*(QSAFE**2)*F(2)*FP*(VTHETa(2)*(sinpot+c8*sc1(2))-
	1	(((VPHIa(2)+PRESS(2))/VTH(2))*SINPOT+c8*sc2(2))) 
	S(4) = S(4) -0.5*QSAFE*EP*(((VPHIa(2)/vth(2))**2)+C9*YV(2)) - 
     2	0.5*QSAFE*(FP**2)*EP*(VTHETa(2)**2) + c7*vs4(2)
	if(ioptedge.eq.1) then
	  c(4,3) = c(4,3) + 0.5*ep*qsafe*hnuion(2)*hnuat(2)
	  c(2,2) = c(2,2) + 0.5*ep*qsafe*hnuion(2)*fp*vtheta(2)
	  xy = fp*vtheta(2)*sindeno(2)+hnuion(2)*cosdeno(2)
	  xz = sindeno(2)*(1.+xne/xni(2))+(xne/xni(2))*sindene
	  yz = (xne/xni(2))*(cosdene+cosdeno(2))
	  s(4) = s(4)-0.5*ep*qsafe*hnuat(2)*xy-qsafe*hnuion(2)*fp*
     1         (0.5*vtheta(2)*xz +0.333*(qsafe/ep)*f(2)*yz)
 	endif
 
	xx = (0.5*qsafe*ep)

	do 3098 il = 1,4
	s(il) = s(il)/xx
	do 3098	jm = 1,4
	c(il,jm) = c(il,jm)/xx
3098	continue 

	COS1 = COSDEN(1)
	COS2 = COSDEN(2)
	SIN1 = SINDEN(1)
	SIN2 = SINDEN(2)		
	CALL LSLRG(4,C,4,S,1,A)
c	limit poloidal variation to epsilon
c	if(ABS(a(1)).ge.1.0) THEN
c		A(1) = 1.0
c		IF(A(1).LT.0.0) a(1) = -1.0
c	ENDIF
c	if(ABS(a(2)).ge.1.0) THEN
c		A(2) = 1.0
c		IF(A(2).LT.0.0) a(2) = -1.0
c	ENDIF
c	if(ABS(a(3)).ge.1.0/ep) THEN
c		A(3) = 1.0/ep
c		IF(A(3).LT.0.0) a(3) = -1.0/ep
c	ENDIF
c	if(ABS(a(4)).ge.1.0) THEN
c		A(4) = 1.0
c		IF(A(4).LT.0.0) a(4) = -1.0
c	ENDIF


	COSDEN(1) = A(1)
	SINDEN(1) = A(2)
	COSDEN(2) = A(3)
	SINDEN(2) = A(4)
c	zg = cosden(2)
c	if(zg.ge.1.0) zg = 1.0 
c	cospot = (gam(1)*cosden(1)+gam(2)*zg)/phite
	
	
	
25	continue	 
	COSDEN(1) = (1.-epd)*COSDEN(1)+epd*COS1
	COSDEN(2) = (1.-epd)*COSDEN(2)+epd*COS2
	SINDEN(1) = (1.-epd)*SINDEN(1)+epd*SIN1
	SINDEN(2) = (1.-epd)*SINDEN(2)+epd*SIN2
	cospot = (gam(1)*cosden(1)+gam(2)*cosden(2))/phite
	sinpot = (gam(1)*sinden(1)+gam(2)*sinden(2))/phite
 
c	electron sin & cos
	cosdene =(xni(1)*atnum(1)*cosden(1)+xni(2)*atnum(2)*cosden(2))/xne 
	sindene =(xni(1)*atnum(1)*sinden(1)+xni(2)*atnum(2)*sinden(2))/xne 
C		POTENTIAL POLOIDAL VARIATION
	cospot = (gam(1)*cosden(1)+gam(2)*cosden(2))/phite
	sinpot = (gam(1)*sinden(1)+gam(2)*sinden(2))/phite

		 
C		velocity asymmetry terms(product [vphi/vth]*[vphicos/ep], etc.)

	vphicos(1)=-1.0*vtheta(1)*(2.+cosden(1)+cospot)  
	1		+ (cospot*press(1) + (1.+cospot)*vphia(1))/VTH(1)
     	vphicos(2)=-1.0*vtheta(2)*(2.+cosden(2)+cospot)  
     1		+ (cospot*press(2) + (1.+cospot)*vphia(2))/VTH(2)
     	
	vphisin(1)=-1.0*vtheta(1)*(sinden(1)+sinpot)
	2		+	sinpot*(press(1)+vphia(1))/VTH(1)
	vphisin(2)=-1.0*vtheta(2)*(sinden(2)+sinpot)
	2		+	sinpot*(press(2)+vphia(2))/VTH(2)
	

		 
 	if(ioptedge.eq.1) then
	  vphicos(1) = vphicos(1) - (xne/xni(1))*hnuion(1)*
     1			  (sindene+sindeno(1))/fp
	  vphicos(2) = vphicos(2) - (xne/xni(2))*hnuion(2)*
     1			  (sindene+sindeno(2))/fp
	  vphisin(1) = vphisin(1) + (xne/xni(1))*hnuion(1)*
	1			   (cosdene+cosdeno(1))/fp
      vphisin(2) = vphisin(2) + (xne/xni(2))*hnuion(2)*
	1			   (cosdene+cosdeno(2))/fp	
     	endif 
c		poloidal variation in poloidal velocity
c		(product [vtheta/fp*vth]*[vethetacos/ep] etc.)
	vthetcos(1) = -1.*(1. + cosden(1))
 	vthetcos(2) = -1.*(1. + cosden(2))
 	vthetsin(1) = -1.*sinden(1)
      vthetsin(2) = -1.*sinden(2)
	if(ioptedge.eq.1) then
	  vthetcos(1) = vthetcos(1) - (xne/xni(1))*hnuion(1)*
     1			  (sindene+sindeno(1))/fp
	  vthetcos(2) = vthetcos(2) - (xne/xni(2))*hnuion(2)*
     1			  (sindene+sindeno(2))/fp
	  vthetsin(1) = vthetsin(1) + (xne/xni(1))*xnuioni(1)*
     1			  (cosdene+cosdeno(1))/fp
	  vthetsin(2) = vthetsin(2) + (xne/xni(2))*xnuioni(2)*
     1			  (cosdene+cosdeno(2))/fp
	endif
C	FROM VTHETAHAT S/C TO VTHETA S/C 	 
	VTHTCOS(1,N) = VTHETCOS(1)*EP*VTHETA(1)*FP*VTH(1)
	VTHTCOS(2,N) = VTHETCOS(2)*EP*VTHETA(2)*FP*VTH(2)
	VTHTSIN(1,N) = VTHETSIN(1)*EP*VTHETA(1)*FP*VTH(1)
	VTHTSIN(2,N) = VTHETSIN(2)*EP*VTHETA(2)*FP*VTH(2)
	VTHTCOS(1,N) = VTHETCOS(1)
	VTHTCOS(2,N) = VTHETCOS(2)
	VTHTSIN(1,N) = VTHETSIN(1)
	VTHTSIN(2,N) = VTHETSIN(2)
c	inertial pinch term
	pinert(n) = (xmas(1)*vphia(1)*vrad1(n)/rmajor)*
	1			((rmajor*xlvm(n))-1.)/(eq*atnum(1)*bthet)
	pinert(n) = pinert(1) + 0.5*(xmas(1)*vphia(1)*vth(1)*fp*ep/
	1	rmajor)*(vthetsin(1)*(3.+cosden(1))-vthetcos(1)*sinden(1) +
     2	vtheta(1)*sinden(1))/(eq*atnum(1)*bthet)
	pinert(n) = pinert(1) + 0.5*(xmas(1)*vtheta(1)*(vth(1)**2)*
     1		fp*ep/rmajor)*vphisin(1)/(eq*atnum(1)*bthet)
      
2001	xv(1) = 0.25*(ep**2)*(vphisin(1)*vphicos(1) + 
     1  (cosden(1)*vphisin(1) + sinden(1)*vphicos(1))*(vphia(1)/VTH(1)))
     	xv(2) = 0.25*(ep**2)*(vphisin(2)*vphicos(2) + 
     1  (cosden(2)*vphisin(2) + sinden(2)*vphicos(2))*(vphia(2)/VTH(2)))

C		del**2 term in inertial source term in sin moment     
      yv(1)=0.5*(ep**2)*(0.5*(vphicos(1)**2)+1.5*(vphisin(1)**2)+
     1 (cosden(1)*vphicos(1) + 3.*sinden(1)*vphisin(1))*vphia(1)/VTH(1))			 
	yv(2)=0.5*(ep**2)*(0.5*(vphicos(2)**2)+1.5*(vphisin(2)**2)+
     1 (cosden(2)*vphicos(2) + 3.*sinden(2)*vphisin(2))*vphia(2)/VTH(2)) 

C		ep**2 term in friction in unity moment
 	fric(1) = -0.5*(EP**2)*(1.+COSDEN(1)+COSDEN(1)**2+SINDEN(1)**2)
	fric(2) = -0.5*(EP**2)*(1.+COSDEN(2)+COSDEN(2)**2+SINDEN(2)**2)

c	viscous quadratic terms in sin & cos eqs
 	do 45 j=1,2 
	dela(j) = -1.*(ep**2)*sinden(j)*vphisin(j)/4.
 
c	cc1(j)=0.125*(ep**2)*((sinden(j)**2)/3. - cosden(j)/3. - 
c	2		(cosden(j)**2)/3. -5.)
c	cc2(j)=0.25*(ep**2)*((1.-cosden(j))*((vphi(j)/cs(j))-vphicos(j)) +
c	1		sinden(j)*vphisin(j))
c	sc1(j)= 0.125*(ep**2)*sinden(j)*(3.+ 2.*cosden(j))
c	sc2(j)= 0.75*(ep**2)*((vphicos(j)-(vphi(j)/cs(j)))*sinden(j) -
c     1	 vphisin(j)*(1.+cosden(j)))
	vc1(j) = (ep**2)*(cosden(j)+1.)/24.
	vc2(j) = -1.*(qsafe**2)*f(j)*fp*(ep**2)*
	1		((vphia(j)/VTH(j))-vphicos(j))/4. 
	vc3(j) = (qsafe**2)*f(j)*fp*(ep**2)*
	1		(vtheta(j)*sinden(j)/24.+vphisin(j)/4.)
      vc4(j) = (qsafe**2)*f(j)*fp*(ep**2)*
     1		(5.*vtheta(j)/8. - ((vphia(j)/VTH(j))-vphicos(j))/4.)
      vs1(j) = (ep**2)*(3./8. + cosden(j)/4.)
 	vs2(j) = -1.*(qsafe**2)*f(j)*fp*(ep**2)*
	1		((vphia(j)/VTH(j))-vphicos(j))*3./4. 
 	vs3(j) = -1.*(qsafe**2)*f(j)*fp*(ep**2)*vphisin(j)*3./4.
      vs4(j) = (qsafe**2)*f(j)*fp*(ep**2)*vphisin(j)*3./4.	 
45	continue 




75	continue
	
c	write(6,301) vphia(1),vtheta(1),vtheta(2),sinden(1),sinden(2),
c     1			cosden(1),cosden(2)	
	if(ioptedgereal.eq.3) ioptedge = 3
c	if(nn.gt.100) zc = 1.
	vt(nn) = vtheta(1) 
c	do not iterate the nmesh = 25 solution (which diverges)
c	if(nmesh.eq.25) goto 150
C		CHECK CONVERGENCE OF ROTATION/ASYMMETRY CALCULATION
	IF(NN.EQ.1) GOTO 100
	XX = ABS(VTHET1/VTHETA(1) - 1.)
 	IF(XX.GT.epol) GOTO 100
	if(ioptedge.ne.3) then
	XX = ABS(VTHET2/VTHETA(2) - 1.)
 	IF(XX.GT.epol) GOTO 100
	endif
	XX = ABS(SIN1/SINDEN(1) - 1.)
	IF(XX.GT.epol) GOTO 100
	XX = ABS(SIN2/SINDEN(2) - 1.)
 	IF(XX.GT.epol) GOTO 100
	XX = ABS(COS1/COSDEN(1) - 1.)
 	IF(XX.GT.epol) GOTO 100
 	XX = ABS(COS2/COSDEN(2) - 1.)
 	IF(XX.GT.epol) GOTO 100
c	if(zc.ge.1.0) then
	
      goto 150 
c	endif
c	zcc = 1.
c	c9 = 1.
100	continue	
900	format((1X,A40,i4))
	WRITE (6,900) 'POLOIDAL ROTATION SOLUTION NOT CONVERGED',nmesh 
150	CONTINUE
	
	
	 
c	if(ioptpert.eq.0) goto 155 
c	if(ioptpert.eq.2) goto 155 
c	if(niter.gt.1) goto 155 
c 	if(ioptpert2.eq.2) then
c  	ioptpert = 2
c	endif
c	if(ioptpert.eq.1) goto 155 
c	goto 5
155	continue 
C	UPDATE DRAG FREQUENCY
	par = (1.-((RMINOR/AMINOR)**2))	 
	dxn0 = xn0-xnped
	dxt0 = t0-tped
	dxv0 = vphexp-vped
	DO 750 J = 1,2
c	xlvm(n) = -1.*(torv(n+1)-torv(n))/(delma*(torv(n+1)+torv(n))) 


c	**************************************
c	if(it.gt.1) xlna(j) = 1./xlnm(nmesh)
c	xlva(j) = 1./xlvm(nmesh)
c	xlva(j) = xlv1
c	****************************************
c	XLINV=AMINOR*((1./XLNa(J))+(1./XLTa(J))+(1./XLVA(J)))
	
c	xlinv=aminor*(xlvm2(nmesh)+xlpm(nmesh))
c	if(j.eq.1) xlinv=aminor*(xlvm(nmesh)+xlpm(nmesh))

c	GRAD(J) = XLINV
c	xlvm(n) = 1./ylvbarx
	grad(j) = rhor(n)*AMINOR*(xlpm(n)+xlvm(n))

	epd = 0.0
	zz = vphia(j)/vth(j)
	

      THETw(J) = EPD*THETw(J)+(1.-EPD)*((4. + COSDEN(J))*vphisin(j)/zz 
	1			+ SINDEN(J)*(1.-vphicos(j)/zz))
c	thetw(j) = (4.+cosden(j))*(sinpot-(vtheta(j)*vth(j)/vphia(j))*
c	1	(sinpot+sinden(j))) + sinden(j)*((vtheta(j)*vth(j)/vphia(j))*
c	2		(2.+cospot+cosden(j))-cospot)
750	continue
	do 751 j=1,2
	k=2
	zz = ATNUM(J)
 	if(j.eq.2) THEN
	k = 1
	zz = ZBAR2(N)
	ENDIF

	XNUDRAG(J)=GRAD(J)*THETw(J)*TEMP(J)/
	2			(2.*(RMAJOR**2)*zz*(BPHI))
	if(xnudrag(j).le.0.0.and.ioptzerodrag.eq.0) xnudrag(j) = 0.0
c	add anomalous drag
	if(ioptdrag.eq.1) xnudrag(j) = xnudrag(j) + anomdrag(j)   
751	continue

	Gy(n) = grad(1)
	cosion(n) = cosden(1)
	cosimp(n) = cosden(2)
	sinion(n) = sinden(1)
	sinimp(n) = sinden(2) 
	thetwid1(n) = thetw(1)
	thetwid2(n) = thetw(2)
	return
	END 	