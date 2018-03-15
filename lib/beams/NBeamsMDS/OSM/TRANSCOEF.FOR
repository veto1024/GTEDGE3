	SUBROUTINE TRANSCOEF
	INCLUDE 'SOLDIV.FI'
	DOUBLE PRECISION TEDP,XLZ1DP,DLDT1DP,XLZ2DP,DLDT2DP 
c		evaluates exp. transport coef from input gsl, n & T
	tisolx = tsepexi
	tesolx = tsepexe
	tsolx = 0.5*(tisolx+tesolx)
	tepedx = tpedexe
	tipedx = tpedexi
	tpedx = 0.5*(tepedx+tipedx)
	tebarx = 0.5*(tesolx+tepedx)
	tibarx = 0.5*(tisolx+tipedx)
	xnbarx = 0.5*(xnpedex+xnsepex) 
	xnsolx = xnsepex
	xnpedx = xnpedex
	xnbarx = 0.5*(xnsolx+xnpedx)
	xlnbarx = 1./gradnbar(3)
	xltebarx = 1./gradtebar(3)
	xltibarx = 1./gradtbar(3)
c	***********use grad scale lengths from data18 to evaluate exp chi's *******
	xltebarx = yltebarx
	xltibarx = yltibarx
c	***************************************************************************	 


	PSOLx = 2.*XK*XNSOL*TSOL
 	PPEDx = 2.*XK*XNPED*TPED
	xv = 1.
	OMEGAPOL = (1.6E-19*BFIELD/XMASS)*BETAG
	CSBAR = SQRT(2.*XK*TISOLx/XMASS)
	RHOBAR = CSBAR/OMEGAPOL
	RHOBARI = RHOBAR*SQRT((TISOLx+TIPEDx)/(2.*TISOLx))
C	evaluate atomic physics cross sections
	CNEUT = 1.
 	IF(IOPTELN.EQ.0) CNEUT = 0.
 	TBARi = TiBARx
	TBARe = TeBARx		
  	IF(TBARi.LT.1.E-1) TBARi = 1.05E-1 
	IF(TBARi.GT.1E3) TBARi = .95E3
	IF(TBARe.LT.1.E-1) TBARe = 1.05E-1 
 	IF(TBARe.GT.1E3) TBARe = .95E3 
     	XNBARR = XNBARx
	IF(XNBARx.GT.1E22) XNBARR = 0.95E22
  	IF(XNBARx.LT.1E16) XNBARR = 1.1E16
	TNBAR = TBARi
	IF(TNBAR.GE.1000) TNBAR	= 995.
	CALL INTERP(TBARe,TBARi,TNBAR,XNBARR)
	SVELBAR = SEL(1)
	SVELBARN= SELN(1)
	SVCXBAR = SCX(1)
	SVATBAR = SEL(1) + SCX(1)
	SVIONBAR = SION(1)
	SVTOTBAR = SION(1) + SEL(1) + SCX(1) + CNEUT*SELN(1)*XNOBAR/XNBARx 
	Y = SVCXBAR + SVELBAR + CNEUT*SVELBARN*XNOBAR/XNBARx
	ASCAT = Y/(Y + SVIONBAR)
	RSCATTB = ASCAT
	VOBAR = SQRT(XK*TnBAR/XMASS)
	EIONBAR = 17.5
	IF(XNBARx.LE.1.E21) 
     2     EIONBAR = 17.5 + (5.+37.5/TBARe)*LOG10(1.E21/XNBARx)
    	IF(XNBARx.GT.1.E21)
     2    EIONBAR = (30.6 - 16.4*EXP(-5.E19/XNBARx))*
     3              EXP(5.45/(TBARe*EXP((XNBARx/1.37E20)**0.26)))
C
	SIGION = SVIONBAR/VOBAR
  	SIGTRAN = (SVIONBAR + SVATBAR - 0.667*SVELBAR)/VOBAR
      
	FLUXHEATE = HEATFRACE*FLUXHEAT
	FLUXHEATI = (1.-HEATFRACE)*FLUXHEAT   
	XNOBARAV = XNOBAR*(1.-WTXPT) + XNOBARXPT*WTXPT
	XNOCOLDBARAV = XNOCOLDBAR*(1.-WTXPT) + XNOCOLDBARXPT*WTXPT
 	EIONM = EIONBAR
	FCTR = 1.
	FCTR2 = 1.
c	 evaluate impurity radiation
	C2 = 1.E-13
	TEDP = TEBARx
	NZIMP1 = IZINJECT
	NZIMP2 = IZINTRIN
	CALL CEFITS(NZIMP1,TEDP,XLZ1DP,1,DLDT1DP)
 	CALL CEFITS(NZIMP2,TEDP,XLZ2DP,1,DLDT2DP)
	XFO = XNOBARAV/XNBARx
c	  IF(NZIMP1.EQ.4.or.nzimp1.eq.6.or.nzimp1.eq.74)  then 
c		CALL CXRCEFITS(NZIMP1,TEDP,XFO,XLZ1DP,DLDT1DP,ZAV1)
c	endif
c
c	IF(NZIMP2.EQ.4.or.nzimp2.eq.6.or.nzimp2.eq.74)  then 
		CALL CXRCEFITS(NZIMP2,TEDP,XFO,XLZ2DP,DLDT2DP,ZAV2)
	zimpmean= zav2
c	endif
	XNZLZBAR =XNBARx*(CFZINJTB*XLZ1DP + CFZINTTB*XLZ2DP)*QZMULTCORE*C2  
	DLZDTBAR =XNBARx*(CFZINJTB*DLDT1DP+CFZINTTB*DLDT2DP)*QZMULTCORE*C2
	ylz = xlz2dp
c	*********************************************************************
c	evaluate transport coefficients	
c		assume density width governs
c	widthnx = 0.0 
	widthx = widthnx
	DPART = (FLUXPART/XNBARx)-0.5*XNOBARAV*SVIONBAR*FCTR*WIDTHX
	1	-VPINCH 
	DEXP = XLNBARx*DPART 
	DE = (((FLUXHEATE)/(XK*TEBARX))-cc*FLUXPART*FCTR2)/XNBARX +	
     2 0.5*WIDTHX*((XNZLZBAR/(XK*TEBARX))+XNOBARAV*SVIONBAR*
     3			(EIONM/TEBARX+cc))
	if(de.le.0.0) DE = ((FLUXHEATE)/(XK*TEBARX*xnbarx)) 
	DI =  (((FLUXHEATI)/(XK*TIBARX))-cc*FLUXPART*FCTR)/XNBARX +
     2  	   0.5*WIDTHX*(1.5*XNOCOLDBARAV*SVATBAR +
     3       cc*XNOBARAV*SVIONBAR)
	if(di.le.0.0) DI = ((FLUXHEATI)/(XK*TIBARX*xnbarx)) 
	CHIEXPE = XLTEBARX*DE 
	CHIEXPI = XLTIBARX*DI
	DQE1= 0.5*WIDTHX*XNBARX*(XNZLZBAR+XNOBARAV*SVIONBAR*EIONM*XK)
     	DQI1= 0.5*WIDTHX*XNBARX*1.5*XNOCOLDBARAV*SVATBAR*XK*TIBARX 
	DG1 = 0.5*WIDTHX*XNBARX*XNOBARAV*SVIONBAR	
c	******if ioptranscoef = 0, evaluate chi's from input gradient scale lengths
c	******				 = 1, use input chi's

	
c	*********************************************************************
c	evaluate chi's for top of pedestal--atomic physics pedestal only
	titop = 0.5*(tiintop + tpedexi)
	tetop = 0.5*(teintop + tpedexe)
	xntop = 0.5*(xnintop + xnpedex) 
c	***********temperature grad scale lengths on top*****************
c	dx = pedrhote - rhotein
c	xltetop = aminor*sqrt(0.5*(1.+elong**2))*dx/log(teintop/tpedexe) 
c	dx = pedrhoti - rhotiin
c	xltitop = aminor*sqrt(0.5*(1.+elong**2))*dx/log(tiintop/tpedexi)
c***********************************************
	drho = pedrhon - rhonin
	xlntop = aminor*sqrt(0.5*(1.+elong**2))*drho/log(xnintop/xnpedex)
c	xltitop = aminor*sqrt(0.5*(1.+elong**2))*drho/log(tiintop/tpedexi) 
c	xltetop = aminor*sqrt(0.5*(1.+elong**2))*drho/log(teintop/tpedexe) 

	drho = 1.0 - pedrhon
	yln	= aminor*sqrt(0.5*(1.+elong**2))*drho/log(xnpedex/xnsepex)
c	ylti = aminor*sqrt(0.5*(1.+elong**2))*drho/log(tpedexi/tsepexi)
c	ylte = aminor*sqrt(0.5*(1.+elong**2))*drho/log(tpedexe/tsepexe)



	fluxparttop = fluxpart - xnbarx*svionbar*xnobarav*widthx
	fluxheatetop =  fluxheate +	
     1	widthx*xnbarx*(xnzlzbar + xnobarav*svionbar*eionm*xk)
	detop = fluxheatetop/(xntop*xk*tetop)-2.5*fluxparttop/xntop
	if(detop.le.0.0) detop =  fluxheatetop/(xntop*xk*tetop)
   	fluxheatitop=fluxheati + widthx*1.5*xnocoldbarav*svatbar*tibarx*xk
	ditop = fluxheatitop/(xntop*xk*titop)-2.5*fluxparttop/xntop
	if(ditop.le.0.0) ditop = fluxheatitop/(xntop*xk*titop)

	if(ioptranscoef.eq.1) goto 101
c	input xltitop & xltetop, compute corresponding chi's

	chetop = detop*xltetop
	chitop = ditop*xltitop 
c	xltitop = chitop/ditop
	
c	use experimental chi to get the heat out with exp Lt
c	chixpe = chiexpe 
c	chixpi = chiexpi 
c	ylti = chixpi/di
c	ylte = chixpe/de
c	input ylti and ylte, compute corresponding chi's
	chixpi = ylti*di
	chixpe = ylte*de
	if(chixpi.lt.0.1) chixpi = 0.1
	if(chixpe.lt.0.1) chixpe = 0.1

	goto 105
101	continue 
c	use input chi's and compute gradient scale lengths
	ylti = chixpi/di
	ylte = chixpe/de 
	xltitop = chitop/ditop
	xltetop = chetop/detop
105	 continue

c	***********************************************************************      

c	evaluation of pressure gradient scale lengths
c		experimental

	game = tebarx/(tebarx+tibarx)
	gamiA = tibarx/(tebarx+tibarx) 
	xlpexp = 1./((1/xlnbarx)+(gamiA/xltibarx)+(game/xltebarx))
c		theoretical
	XV = (0.5*gamiA)**0.333
	fshape = 1.
	DPDRCRITnom = (1./((Q95**2)*RMAJOR))*((BFIELD**2)/2.52E-6) 
	PCRIT = 2.667*(BFIELD**2/2.514E-6)*((RHOBARI/RMAJOR)**.667)*
     2	XV*FSHAPE/(Q95**2)
	pressbar = xk*xnbarx*(tibarx+tebarx)
	ccritex3 = pressbar/pcrit
	cmhdex = pressbar/(xlpexp*dpdrcritnom)
	xlpth1 = pressbar/dpdrcritnom
	xlpth2 = pcrit/dpdrcritnom
c	evaluate theoretical widths
c		neutral penetration	
	deltran = 1./(XNBARx*SQRT(3.*SIGION*SIGTRAN))
 	cmfpex = widthnx/deltran
c		critical pressure grad + critical pressure R&D
	deltap2 = ccritex3*PCRIT*LOG(2.*(ccritex3*PCRIT/PSEP)-1.)/
     2	(cmhdex*DPDRCRITnom)
c		critical pressure R&D 
c		use 'exp' trans coef determined above
	BFAC = ((FLUXPART/XNBARx)-VPINCH)/(dexp) +
     2   GAMIA*((FLUXHEATI/(XK*TIBARx))-cc*FLUXPART)/
     3	(XNBARx*chiexpi) + 
     4    GAME*((FLUXHEATE/(XK*TEBARx))-cc*FLUXPART)/
     5	(XNBARx*chiexpe)
     	
      AFAC = 0.5*(GAMIA*(1.5*XNOCOLDBARAV*SVATBAR
     3       + cc*XNOBARAV*SVIONBAR)/(CHIexpI) +
     4	  GAME*((XNZLZBAR/(XK*TEBARx))+XNOBARAV*SVIONBAR*
     3 (EIONM/TEBArx) +cc)/(CHIexpE) -
	5  XNOBARAV*SVIONBAR/(dexp))
	G = LOG(2.*(ccritex3*PCRIT/PSEP) - 1.)*FSHAPE
 
C		FIX 6/21/02
	IF((1.+4.*AFAC*G/(BFAC**2)).LE.0.0) AFAC = 0.
	  IF(AFAC.EQ.0.0) DELTap1 = G/BFAC
	  IF(AFAC.NE.0.0)  THEN
      	DELTap11 = (BFAC/(2.*AFAC))*(SQRT(1.+4.*AFAC*G/(BFAC**2))-1.)
	  endif
c	with deltan separately determined
		AFAC = 0.5*(GAMIA*(1.5*XNOCOLDBARAV*SVATBAR
     3       + cc*XNOBARAV*SVIONBAR)/(chiexpi) +
     4	  GAME*((XNZLZBAR/(XK*TEBARx))+XNOBARAV*SVIONBAR*
     3 ((EIONM/TEBARx) +cc))/(chiexpe))
    	BFAC = ((FLUXPART/XNBARx)-VPINCH-0.5*XNOBARAV*SVIONBAR*WIDTHNx)/
	1	(dexp) +
     2   GAMIA*((FLUXHEATI/(XK*TIBARx))-cc*FLUXPART)/
     3	(XNBARx*chiexpi) + 
     4    GAME*((FLUXHEATE/(XK*TEBARx))-cc*FLUXPART)/
     5	(XNBARx*chiexpe)
 
	  IF(AFAC.EQ.0.0) DELTap12 = G/BFAC
	  IF(AFAC.NE.0.0)  THEN
		DELTap12 = (BFAC/(2.*AFAC))*
	2						(SQRT(1.+4.*AFAC*G/(BFAC**2))-1.)
        ENDIF
	
C	EVALUATE PRESSURE WIDTH USING EXP Lp
c		pcrit = 1/2(pped+psep)
	psolx = xk*xnsepex*(tisolx+tesolx) 
	deltap3 = xlpexp*log(2.*(pcrit/psolx) - 1.)
c		pcrit = psep
	deltap4 = xlpexp*log(pcrit/psolx) 
	presspedx = xk*xnpedex*(tepedx+tipedx)
	ccritex4 = presspedx/pcrit      
c	***********************************************************
C		evaluate of l-h threshold eq 7 pop, 8/02
c	xkredge = 1./widthnx 
	xkredge = 1.0 
	chihmi = chiexpi
	chihme = chiexpe
	chihmi = chixpi
	chihme = chixpe

c	chihme = 1.0
c	chihmi = 1.0

	alphae =  ((xnu*XNZLZBAR/(xk*TEBARX))-DLZDTBAR/xk)*qzmultcore +
	1XNOBARAV*SVIONBAR*(xnu*EIONM/TEBARX+2.5*(xnu-1.))
	alphai =  (1.5*(xnu-1.)*XNOCOLDBARAV*SVATBAR +
     3       cc*(xnu-1.)*XNOBARAV*SVIONBAR)
c	correct for time-dependent W
	alphae = alphae + 1.5*xnu*dlnW_dt*fheate
	alphai = alphai + 1.5*xnu*dlnW_dt*(1.-fheate) 
	fluxpartav = fluxpart - 0.5*xnbarx*XNOBARav*SVIONBAR*WIDTHNX
c	fluxheattbe = fluxheate + 0.5*WIDTHX*((XNZLZBAR/(XK*TEBARX))+
c    1	XNOBARAV*SVIONBAR*1.5*xnbarx)
c	fluxheattbi = fluxheati + 0.5*WIDTHX*(1.5*XNOCOLDBARAV*SVATBAR)*
c	1							xnbarx
c	goto 124	
	zz = 1.25*fluxpartav/xnbarx
	qeav = zz*(sqrt(1.+chihme*(alphae-chihme*(xkredge**2))/
	1		(xnu*(zz**2)))	+ 1.)*xnbarx*xk*tebarx	 
 	qiav = zz*(sqrt(1.+chihmi*(alphai-chihmi*(xkredge**2))/
	1		(xnu*(zz**2)))	+ 1.)*xnbarx*xk*tibarx
	qesep = qeav-0.5*xnbarx*(XNZLZBAR+XNOBARav*SVIONBAR*EIONM*xk)
	1												*widthnx
	qisep = qiav-0.5*xnbarx*XNOCOLDBARav*SVATBAR*1.5*xk*widthnx
c	****note that xpt neutral densities are used here *******
	surfB = 39.479*aminor*sqrt((1.+elong**2)/2.)*rmajor 
      pthreshe = qesep*surfB*1e-6
	pthreshi = qisep*surfB*1e-6
	xxx = chihme*(alphae-chihme*(xkredge**2))/
	1		(xnu*(zz**2))
	yyy = chihmi*(alphai-chihmi*(xkredge**2))/
	1		(xnu*(zz**2))
	zzz = XNOBARAV*SVIONBAR*(xnu*EIONM/TEBARX+2.5*(xnu-1.))
	zzy	= ((XNZLZBAR/(XK*TEBARX))-DLZDTBAR/xk)*qzmultcore
	zzq = ((XNZLZBAR/(XK*TEBARX)))
c	2nd formulation 3/24/06
c124	pthreshe = AP*xnbarx*xk*tebarx*(alphae-chihme*(xkredge**2))*1e-6 
c	pthreshi = AP*xnbarx*xk*tibarx*(alphai-chihmi*(xkredge**2))*1e-6
c	pthreshe = pthreshe*gradtebar(4)/xnu
c	pthreshi = pthreshi*gradtbar(4)/xnu
 
	psepe = AP*fluxheat*fheate*1e-6
	psepi = AP*fluxheat*(1.-fheate)*1e-6
	psepthe = psepe/pthreshe
	psepthi = psepi/pthreshi
	ylamdae=chihme/(alphae-xnu*psepe/(xnbarx*tebarx*gradtebar(4)*ap))
	ylamdai=chihmi/(alphai-xnu*psepi/(xnbarx*tibarx*gradtbar(4)*ap))

	RETURN
	END