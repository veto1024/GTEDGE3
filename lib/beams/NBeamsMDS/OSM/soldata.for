
      SUBROUTINE SOLDATA(nbdep1,nbdep2,nbdep3,erextot,fpsi0)
	INCLUDE 'SOLDIV.FI'

	real ibal,imod,ioptsep,vplas,ioptgeom,dchirat,changeouter,tw,
	1     gammid,xlmid,czd(25),cmid,delsol,rwall,ioptdel,ioptpump,
	2     hxpt
	real nbdep1(51),nbdep2(51),nbdep3(51), erextot(51), fpsi0(51)

	
	NAMELIST/GEOM/IOPTSN,ZX,RX,RSEP1,ZSEP1,RSEP2,ZSEP2,RDOME,PHIDOME,
	2		XRING,ZRING,DELLT1,DELLT2,YPUMP1,YPUMP2,
     3		RPUMP1,RPUMP2,R1,R2,PHITRAN,RINROUT,XLWPL1,XLWPL2,XLUP1,
     4		XLWPF1,XLWPF2,XLWSPL,elong,AMINOR,RMAJOR,TRIANG,DELXreal
     5		,DELXPT,zdome,xdome

	
	NAMELIST/CORE/IPFLUX,IBAL,IMOD,ISHEATH,FLUXPART,FLUXHEAT,FRACRAD,
     2    	TSEP,TD,TSOL,TBAR,TPED,TPLAS,TAV,XND,XNOD,XNPLAS,XNPED,
     3        XNBAR,XNSEP,XNSOL,XNDIV,XNOSOL,XNODIV,IOPTSEP,DELRATNT,
     4		DELRATNV,IOPTDEL,PFUSION,FLUXNEUTIN,FLUXIONIN,IOPTPROF,
	5		IOPTOHM,POHMIN
	

	NAMELIST/PARAM/B,IOPTQ95,Q95,GAMSHEATH,IOPTRAD,TAURES,IZINJECT,
	1        FZINJECT,AZINJECT,AZTRACE,IZTRACE,REDKAPPA,SOLNMULT,
     2 SOLTMULT,IZINTRIN,FZINTRIN,FHE,FBE,VPLAS,DELSOL,IOPTSOL,XNATOM,
     3    XLZRAD,FS,ALB,IOPTGEOM,CHANGE,CHANGEOUTER,DELNMAX,DELNMIN,
     4        DCHIRAT,ENHANCE,QZMULTDIV,CHIRSOL,FUDGE,REDKAPPA2,
	5  XDSEP,ASYMPART,ASYMHEAT,ioptdiv,deld,divdepth,xenhance,ZHGHT
	NAMELIST/GASNEUT/REFLECT,TW,GAMDIV,GAMMID,XLMID,CZD,CMID,EIOND,
     2                 EIONSOL,RWALL,EEX
	NAMELIST/XSECT1/IOPTLYMAN,ELAST,ELASTN
	NAMELIST/XSECT2/CX,EION
	NAMELIST/XSECT3/REC,TINT,ZNINT,FREC
	NAMELIST/GASNEUT2/BETAG,FPF,FPD,THETA,FMOL,FUELPF,FUELPLOUT,
     2	FUELPLIN,FUELMP,SPELLET,RWPF,RWPL,HXPT,
     3	RWSOL,XMSOL,HNDIV,HVSOL,HVDIV,HTSOL,HTDIV,DELTB,DELPED,
     4	RWDP,IOPTPUMP,DELPUMP,RPUMP,PUMPSPEED,HNSOL,wtxpt,
	5	EPDIV,FRACXPT,IMAT,IMOL,FEX,TSOR,AMU,wallmult,iopttheta,
	6	REABSORBD,REABSORBDIV,REABSORBSOL,REABSORBCOR,IOPTELN,IOPTSOLN
 	NAMELIST/PEDBAR/NBAR,NPED,GRADTPEDIS,XNU,ZEFF,Z0,CFZINJ,
     2	CFZINT,cfzinjtb,cfzinttb,ALPHAT1,ALPHAT2,IOPTPEDWIDTH,CPED,
     3	A1,A2,CHIREDGE,PLASMACUR,PBEAM,ENBI,PHIIN,RADMULT,TAVMULT,
	2    IOPTPART,HRAT,RHOGAUSS,SIGMA,ZEFFC,DPERP,TAUE,TAURATIO,HCONF,
	3	c89,CHEI,ALPHAN,RCOEF,XNCTRPED,TCTRPED,NNEUT,IOPTPED,FION,
	4	RXLTIXLTE,ERAD,VPHI,VTHET,VRAD,XNUBAR,XLQ,XLE,XLVTHET,ISHAPE,
	5	C89N,CC,XLVPHI,XNAV,AMASS,XNRHOPED,RSYN,IPROF,JJOPTPED,
	6 IOPTGRAD,FCOND,CBALLOON,HN,QZMULTCORE,PEDFAC,CHIRTBI,CHIRTBE,
	7	cPINCH,VPINCH,DTBI,IOPTLN,CIE,CNE,HEATFRACE,IOPTPINCH,CPED2,
     8  CPED1,JOPTPED2,fluxheatianom,fluxheateanom,IGRADPED,FI,IOPTCONF,
	9	IDRAKE,SHEARQ,CBALL,FISEP,CHITBLI,CHITBLE,SHEAR0,CSHEAR,ASS,
	1	IOPTTAUE,IOPTTRAN,MODPED,FSHAPE,CMHD,ccrit,cmfp

	NAMELIST/EXPDATA/XNPEDEX,XNSEPEX,TPEDEXI,TSEPEXI,GRADNBAR,GRADTBAR
     2,GRADNPED,GRADTPED,GRADTEBAR,CHIRC,CHIBAR,CHIRMARF,APED,TPEDAV, 
	3DREDGE,TPEDEXE,TSEPEXE,WIDTHNX,WIDTHTEX,WIDTHTIX, 
	4  dln_dt,dlnw_dt,wmhd,pradcorex,ssi95,chii0,chie0,xkrtherm,
	5  dlnnped_dt,dlnpped_dt
	NAMELIST/EDGE/IOPTEDGE,PHITE,IGRAD,XLT11,XLV1,FRACZ,ZION,ZIMP,
     1		AION,AIMP,BPHI,YLTIbarx,YLNBARX,YLVBARX,eb,abeam,rtan,
     2		alphain,VPHIAPEDX,DELNA,ephia,enh,sindeno,yltebarx,
     3		joptedped,SHEARM,dedr,sheare,ioptvphi,ioptpinchi,
     4	ioptvdif,ioptxlvm,sheareconst,fheate,ioptran,
	5	chixpi,chixpe,chitop,chetop,teintop,tiintop,xnintop,xltitop,
     5		xltetop,pedrhoti,pedrhote, rhotein,rhotiin,ioptdrag,ylti,
	6		ylte,ioptzerodrag,ioptapproach,iopterad,ioptvthet,
	7		ioptranscoef,erex,  ioptexpdata,xnesepreal,
	8		pedrhon,rhonin,OMEGTSEP1,OMEGTSEP2,OMEGTSEP1C,OMEGTSEP2C,
	9		OMEGTSEP1S,OMEGTSEP2S,anomdrag
	1		,delbsep,ylntop,VTHEXP,TORV,ioptpol


	NAMELIST/EXPROFILE1/EXNE
	NAMELIST/EXPROFILE2/XLNE,exlv
	NAMELIST/EXPROFILE3/XTE
	NAMELIST/EXPROFILE4/EXLTE
	NAMELIST/EXPROFILE5/XTI
	NAMELIST/EXPROFILE6/EXLTI

	NAMELIST/TIMEDEP/dlnn_dt,dlnwe_dt,dlnwi_dt
		
	OPEN(UNIT=10,FILE='DATA10',STATUS='OLD')
      READ (10,GEOM)
      CLOSE (UNIT=10,STATUS='KEEP') 
	OPEN(UNIT=11,FILE='DATA11',STATUS='OLD')
      READ (11,CORE)
 	CLOSE (UNIT=11,STATUS='KEEP')

     
	OPEN(UNIT=12,FILE='DATA12',STATUS='OLD')
      READ (12,PARAM)
      CLOSE (UNIT=12,STATUS='KEEP')
	OPEN(UNIT=13,FILE='DATA13',STATUS='OLD')
      READ (13,GASNEUT)
      CLOSE (UNIT=13,STATUS='KEEP')
C	OPEN(UNIT=14,FILE='DATA14',STATUS='OLD')
C      READ(14,XSECT1)
C      READ(14,XSECT2)
C	READ(14,XSECT3)
C	CLOSE(UNIT=14,STATUS='KEEP')
		IOPTLYMAN = 0
	ELAST(1,1) = 2.7E-14
	ELAST(1,2) = 2.7E-14
	ELAST(1,3) = 4.6E-14
	ELAST(1,4) = 9.2E-14
	ELAST(1,5) = 9.2E-14
	ELAST(1,6) = 9.2E-14 
	ELAST(2,1) = 2.7E-14
	ELAST(2,2) = 2.7E-14
	ELAST(2,3) = 4.6E-14
	ELAST(2,4) = 9.2E-14
	ELAST(2,5) = 9.2E-14
	ELAST(2,6) = 9.2E-14 
	ELAST(3,1) = 4.6E-14
	ELAST(3,2) = 4.6E-14
	ELAST(3,3) = 6.8E-14
	ELAST(3,4) = 9.0E-14
	ELAST(3,5) = 9.0E-14
	ELAST(3,6) = 9.0E-14
	ELAST(4,1) = 9.2E-14
	ELAST(4,2) = 9.2E-14
	ELAST(4,3) = 9.0E-14
	ELAST(4,4) = 16.0E-14
	ELAST(4,5) = 16.0E-14
	ELAST(4,6) = 16.0E-14
	ELAST(5,1) = 9.2E-14
	ELAST(5,2) = 9.2E-14
	ELAST(5,3) = 9.0E-14
	ELAST(5,4) = 16.0E-14
	ELAST(5,5) = 16.0E-14
	ELAST(5,6) = 16.0E-14
	ELAST(6,1) = 9.2E-14
	ELAST(6,2) = 9.2E-14
	ELAST(6,3) = 9.0E-14
	ELAST(6,4) = 16.0E-14
	ELAST(6,5) = 16.0E-14
	ELAST(6,6) = 16.0E-14 
	ELASTN(1) = 2.7E-15
	ELASTN(2) = 6.8E-15
	ELASTN(3) = 16.0E-15
	ELASTN(4) = 16.0E-15
	ELASTN(5) = 16.0E-15
	ELASTN(6) = 16.0E-15 

	CX(1,1)    = 8.0E-15
	CX(1,2)    = 8.0E-15
	CX(1,3)    = 1.2E-14
	CX(1,4)    = 2.8E-14
	CX(1,5)    = 5.6E-14
	CX(1,6)    = 1.2E-13
	CX(2,1)    = 8.0E-15
	CX(2,2)    = 8.0E-15
	CX(2,3)    = 1.2E-14
	CX(2,4)    = 2.8E-14
	CX(2,5)    = 5.6E-14
	CX(2,6)    = 1.2E-13 
	CX(3,1)    = 1.2E-14
	CX(3,2)    = 1.2E-14
	CX(3,3)    = 1.5E-14
	CX(3,4)    = 2.9E-14
	CX(3,5)    = 5.8E-14
	CX(3,6)    = 1.2E-13
	CX(4,1)    = 2.8E-14
	CX(4,2)    = 2.8E-14
	CX(4,3)    = 2.9E-14
	CX(4,4)    = 3.7E-14
	CX(4,5)    = 7.4E-14
	CX(4,6)    = 1.5E-13
	CX(5,1)    = 2.8E-14
	CX(5,2)    = 2.8E-14
	CX(5,3)    = 2.9E-14
	CX(5,4)    = 3.7E-14
	CX(5,5)    = 7.4E-14
	CX(5,6)    = 1.5E-13
	CX(6,1)    = 2.8E-14
	CX(6,2)    = 2.8E-14
	CX(6,3)    = 2.9E-14
	CX(6,4)    = 3.7E-14
	CX(6,5)    = 7.4E-14
	CX(6,6)    = 1.5E-13 
	eion(1,1)  = 3.0E-21
	eion(1,2)  = 7.6E-21
	eion(1,3)  = 5.3E-15
	eion(1,4)  = 3.1E-14
	eion(1,5)  = 1.5E-14
	EION(1,6) =	 7.5E-15
	eion(2,1)  = 3.0E-21
	eion(2,2)  = 7.6E-21
	eion(2,3)  = 8.0E-15
	eion(2,4)  = 3.7E-14
	eion(2,5)  = 1.9E-14
	EION(2,6) =	 9 5E-15
	eion(3,1)  = 3.0E-21
	eion(3,2)  = 7.6E-21
	eion(3,3)  = 1.2E-14
	eion(3,4)  = 4.0E-14
	eion(3,5)  = 2.0E-14
	EION(3,6) =	 1.0E-14
	eion(4,1)  = 3.0E-21
	eion(4,2)  = 7.6E-21
	eion(4,3)  = 2.2E-14
	eion(4,4)  = 6.0E-14
	eion(4,5)  = 3.0E-14
	EION(4,6) =	 1.5E-14
	eion(5,1)  = 3.0E-21
	eion(5,2)  = 7.6E-21
	eion(5,3)  = 2.4E-14
	eion(5,4)  = 8.0E-14
	eion(5,5)  = 4.0E-14
	EION(5,6) =	 2.0E-14
	
	REC(1,1)   = 3.0E-16
	REC(1,2)   = 7.0E-17
	REC(1,3)   = 8.0E-20
	REC(1,4)   = 7.0E-21
	REC(1,5)   = 4.0E-22
	REC(2,1)   = 2.5E-15
	REC(2,2)   = 1.1E-18
	REC(2,3)   = 8.0E-20
	REC(2,4)   = 7.0E-21
	REC(2,5)   = 4.0E-22
	REC(3,1)   = 2.5E-14
	REC(3,2)   = 2.5E-18
	REC(3,3)   = 8.0E-20
	REC(3,4)   = 7.0E-21
	REC(3,5)   = 4.0E-22
	REC(4,1)   = 2.0E-13
	REC(4,2)   = 7.0E-18
	REC(4,3)   = 8.0E-20
	REC(4,4)   = 7.0E-21
	REC(4,5)   = 4.0E-22
	REC(5,1)   = 1.0E-12
	REC(5,2)   = 2.0E-17
	REC(5,3)   = 8.0E-20
	REC(5,4)   = 7.0E-21
	REC(5,5)   = 4.0E-22
	TINT(1) = 0.1
	TINT(2) = 1.0
	TINT(3) = 10.0
	TINT(4) = 100.0
	TINT(5) = 1000.0
	TINT(6) = 10000.0
	ZNINT(1) = 1.0E16
	ZNINT(2) = 1.0E18
	ZNINT(3) = 1.0E20
	ZNINT(4) = 1.0E21
	ZNINT(5) = 1.0E22
	FREC = 1.
	
C	ALL DATA ARE SIGMA-V FROM E.W. THOMAS EVALUATION 9/94.
C	SCHULTZ ELASTIC, JANEV C-X,IONIZATION ****not****LYMAN SUPPRESSED
C	RECOMBINATION FROM POST ALB=JANEV W/COLLISION-RADIATIVE, POST ET AL 10/9/98
C	FIRST INDEX Tn= 1, 10, 100 eV; 2ND INDEX Ti=.1,1,10,100,1000 eV: EL & CX
C    	1ST INDEX LOG10Ne=16,18,20,21,22 MKS; 2ND INDEX Te=.1,1,10,100,1000 eV: ION & RECOM

     
      OPEN(UNIT=15,FILE='DATA15',STATUS='OLD')
      READ(15,GASNEUT2)
	CLOSE(UNIT=15,STATUS='KEEP')
	OPEN(UNIT=16,FILE='DATA16',STATUS='OLD')
      READ(16,PEDBAR)
	CLOSE(UNIT=16,STATUS='KEEP')
	OPEN(UNIT=17,FILE='DATA17',STATUS='OLD')
      READ (17,EXPDATA)
      CLOSE (UNIT=17,STATUS='KEEP')
	OPEN(UNIT=18,FILE='DATA18',STATUS='OLD')
      READ (18,EDGE) 	
      CLOSE (UNIT=18,STATUS='KEEP')	
 
 
	heatfrace = fheate

	ixx = ioptq95
	thetain = theta
	q95in = q95
	CHIRCORE(1) = CHIRC
	SPELLETREAL = SPELLET
	FUELMPREAL = FUELMP  


10	PFUSIONIN = PFUSION 
	BFIELD = B
  	XMASS = XNATOM*1.67E-27
	XK = 1.6E-19 

	FESEP = 1. - FISEP

C		SHEATH BOUNDARY CONDITION
	FSHEATH = 1.0
	IF(ISHEATH.EQ.1) FSHEATH = 1.0/BETAG
C		ATOMS ONLY RECYCLING
	IF(IMOL.NE.0) GOTO 100 
C	RWDP = 1.0	 
	FMOL = 0.0
C		ATOMS + GROUND STATE MOLECULES RECYCLING
100	IF(IMOL.EQ.1) FEX = 0.0
C		CALCULATE RATIO EXP TI & TE GRAD SCALE LENGTHS IN TB
	IF(JJOPTPED.EQ.9) FION = GRADTBAR(4)/GRADTEBAR(4)
	
C		INVERT EXPERIMENTAL GRAD SCALE LENGTHS
	DO 125 I=1,5
	IF(GRADNBAR(I).GT.0.0) GRADNBAR(I) = 1./GRADNBAR(I) 
	IF(GRADTBAR(I).GT.0.0) GRADTBAR(I) = 1./GRADTBAR(I)
	IF(GRADTEBAR(I).GT.0.0) GRADTEBAR(I) = 1./GRADTEBAR(I)
125	CONTINUE
	
C	CONSTRUCT EXPERIMENTAL TRANS BARRIER PARAMETERS
	XNTBEX = 0.5*(XNPEDEX + XNSEPEX)
	TTBEX  = 0.25*(TPEDEXE+TPEDEXI + TSEPEXE+TSEPEXI)
	TTBEXE = 0.5*(TPEDEXE + TSEPEXE)
	TTBEXI = 0.5*(TPEDEXI + TSEPEXI)  

C	SET CHIR IN TB IN CGS FOR GROWTH RATE CALCULATION
	
C	SET CHIR IN TB FOR PREDICTIVE MODE
							 
	CHITBI = 1.E4*CHIREDGE/(1.+CHEI)
	CHITBE = 1.E4*CHIREDGE/(1.+(1./CHEI))
	
c      CHIRTBI = CHITBI 	 
c 	CHIRTBE = CHITBE
	
	CHIRTBI = 1.E4*CHIRTBI
	CHIRTBE = 1.E4*CHIRTBE
	
 
C	CONFINEMENT ADJUSTMENT FACTOR
C		TRIANGULARITY DEPENDENCE OF H89, OSBORNE ET AL, PPCF 42,A175,2000
C	IF(JJOPTPED.EQ.10) HCONF = HCONF + 0.4*TRIANG	
	IF(IOPTTAUE.EQ.89) H89 = HCONF	
	IF(IOPTTAUE.EQ.98) H98 = HCONF	
	
     
c	set certain signals to used input exp data when code is used
c	to interpret experimental data
c		for interpretting edge D and chi, use exp nsol & tsol
	if(jjoptped.eq.11) then
	ioptsoln = 1



	ioptpedn = 1
	endif 
c	123302 @ 2600 H-mode   80-99% inter-ELM interval  NEW DATA 8/4/10
	 IF(JJOPTPED.EQ.9) IOPTSOLN = 1  
	


 
	EXNE(1) = 3.43e19 
	EXNE(2) = 3.43E19
	EXNE(3) = 3.43E19
	EXNE(4) = 3.43E19
	EXNE(5) = 3.42E19
	EXNE(6) = 3.41E19
	EXNE(7) = 3.40E19
	EXNE(8) = 3.39E19
	EXNE(9) = 3.38E19
	EXNE(10) =3.37E19
	EXNE(11) =3.36E19
	EXNE(12) =3.35E19
	EXNE(13) =3.34E19
	EXNE(14) =3.32E19
	EXNE(15) = 3.31E19
	EXNE(16) = 3.30E19
	EXNE(17) = 3.29E19
	EXNE(18) = 3.24E19 
	EXNE(19) = 3.15E19
	EXNE(20) = 2.95E19
	EXNE(21) = 2.70E19
	EXNE(22) = 2.40E19
	EXNE(23) = 1.87E19
	EXNE(24) = 1.50E19
	EXNE(25) = 0.75E19
	EXNE(26) = 0.50E19

	XLNE(1) =  2.0
	XLNE(2) =  2.0
	XLNE(3) =  2.0
	XLNE(4) =  2.0
	XLNE(5) =  2.0
	XLNE(6) =  2.0
	XLNE(7) =  2.0
	XLNE(8) =  2.0
	XLNE(9) =  2.0
	XLNE(10) = 2.0
	XLNE(11) = 2.0
	XLNE(12) = 1.0
	XLNE(13) = .530
	XLNE(14) = .268
	XLNE(15) = .156
	XLNE(16) = .094
	XLNE(17) = .069
	XLNE(18) = .062
	XLNE(19) = .050
	XLNE(20) = .044
	XLNE(21) = .037
	XLNE(22) = .031
	XLNE(23) = .019
	XLNE(24) = .015
	XLNE(25) = .012
	XLNE(26) = .010
	
	exlv(1) = .146
	exlv(2) = .143
	exlv(3) = .139
	exlv(4) = .133
	exlv(5) = .126
	exlv(6) = .122
	exlv(7) = .117
	exlv(8) = .112
	exlv(9) = .106
	exlv(10) = .101
	exlv(11) = .095
	exlv(12) = .090
	exlv(13) = .086
	exlv(14) = .082
	exlv(15) = .078
	exlv(16) = .073
	exlv(17) = .069
	exlv(18) = .066
	exlv(19) = .060
	exlv(20) = .054
	exlv(21) = .049
	exlv(22) = .044
	exlv(23) = .039
	exlv(24) = .036
	exlv(25) = .032

		EXLTE(1) = .362
		EXLTE(2) = .342
		EXLTE(3) = .307
		EXLTE(4) = .277
		EXLTE(5) = .252
		EXLTE(6) = .217
		EXLTE(7) = .195
		EXLTE(8) = .175
		EXLTE(9) = .140
		EXLTE(10) = .122
		EXLTE(11) = .107
		EXLTE(12) = .092
		EXLTE(13) = .080
		EXLTE(14) = .067
		EXLTE(15) = .062
		EXLTE(16) = .052
		EXLTE(17) = .047
		EXLTE(18) = .042
		EXLTE(19) = .032
		EXLTE(20) = .027
		EXLTE(21) = .025
		EXLTE(22) = .020
		EXLTE(23) = .015
		EXLTE(24) = .009
		EXLTE(25) = .007
		exlte(26) = .007
	XTI(1) = 1780.
	XTI(2) = 1740.
	XTI(3) = 1720.
	XTI(4) = 1682.
	XTI(5) = 1640.
	XTI(6) = 1608.
	XTI(7) = 1575.
	XTI(8) = 1525.
	XTI(9) = 1497.
	XTI(10) = 1465.
	XTI(11) = 1420.
	XTI(12) = 1370.
	XTI(13) = 1340.
	XTI(14) = 1285.
	XTI(15) = 1242.
	XTI(16) = 1250.
	XTI(17) = 1138.
	XTI(18) = 1073.
	XTI(19) = 950.
	XTI(20) = 850.
	XTI(21) = 720.
	XTI(22) = 590.
	XTI(23) = 475.
	XTI(24) = 425.
	XTI(25) = 400.
	XTI(26) = 250.

	EXLTI(1) = .156
	EXLTI(2) = .152
	exLTI(3) = .148
	EXLTI(4) = .144
	EXLTI(5) = .140
	EXLTI(6) = .136
	EXLTI(7) = .132
	EXLTI(8) = .128
	EXLTI(9) = .124
	EXLTI(10) = .120
	EXLTI(11) = .116
	EXLTI(12) = .111
	EXLTI(13) = .107
	EXLTI(14) = .102
	EXLTI(15) = .097
	EXLTI(16) = .092
	EXLTI(17) = .087
	EXLTI(18) = .082
	EXLTI(19) = .078
	EXLTI(20) = .073
	EXLTI(21) = .069
	EXLTI(22) = .064
	EXLTI(23) = .059
	EXLTI(24) = .055
	EXLTI(25) = .050
	EXLTI(26) =	.045
		
	XTE(1) = 1350.
	XTE(2) = 1350.
	XTE(3) = 1311.
	XTE(4) = 1285.
	XTE(5) = 1260.
	XTE(6) = 1238.
	XTE(7) = 1223.
	XTE(8) = 1200.
	XTE(9) = 1171.
	XTE(10) =1151.
	XTE(11) =1120.
	XTE(12) =1090.
	XTE(13) =1075.
	XTE(14) = 1050.
	XTE(15) = 1025.
	XTE(16) = 992.
	XTE(17) = 988.
	XTE(18) = 920.
	XTE(19) = 865.
	XTE(20) = 780.
	XTE(21) = 670.
	XTE(22) = 400.
	XTE(23) = 320.
	XTE(24) = 200.
	XTE(25) = 091.
	XTE(26) = 070.

		qedge(1) = 2.91
		qedge(2) = 2.97
		qedge(3) = 3.02
		qedge(4) = 3.06
		qedge(5) = 3.11
		qedge(6) = 3.17
		qedge(7) = 3.23
		qedge(8) = 3.25
		qedge(9) = 3.39
		qedge(10) = 3.45
		qedge(11) = 3.55
		qedge(12) = 3.61
		qedge(13) = 3.66
		qedge(14) = 3.75
		qedge(15) = 3.84
		qedge(16) = 3.92
		qedge(17) = 4.00
		qedge(18) = 4.08
		qedge(19) = 4.16
		qedge(20) = 4.27
		qedge(21) = 4.32
		qedge(22) = 4.42
		qedge(23) = 4.42
		qedge(24) = 4.42
		qedge(25) = 4.42
		qedge(26) = 4.42



		
		
		
		
		
		


 

	

	
c	necessary to copy following files from data18
c	123302 @ 2600  H-mode 80-99% between ELMs  NEW DATA 8-4-10
	erex(1) = 1.00e3
	erex(2) = 1.00e3
	erex(3) = 23.4e3
	erex(4) = 21.6e3
	erex(5) = 19.5e3
	erex(6) = 17.8e3
	erex(7) = 15.0e3
	erex(8) = 12.2e3
	erex(9) =  9.46e3
	erex(10) = 7.01e3
	erex(11) = 3.70e3
	erex(12) = 1.15e3
	erex(13) = -2.40e3
	erex(14) = -6.50e3
	erex(15) = -10.5e3
	erex(16) = -14.5e3
	erex(17) = -19.7e3
	erex(18) = -23.7e3
	erex(19) = -27.5e3
	erex(20) = -29.0e3
	erex(21) = -28.0e3
	erex(22) = -21.0e3
	erex(23) = -13.2e3
	erex(24) = -5.0e3
	erex(25) = 1.0e3

c     Radial electric field for rho [0,1], deltarho = 0.02
      erextot(1) =  12.0214E3
      erextot(2) =   9.6798E3
      erextot(3) =  12.4080E3
      erextot(4) =  16.2712E3
      erextot(5) =  20.4616E3
      erextot(6) =  25.4903E3
      erextot(7) =  31.4465E3
      erextot(8) =  36.2645E3
      erextot(9) =  40.8774E3
      erextot(10) =  45.2762E3
      erextot(11) =  49.2980E3
      erextot(12) =  52.9664E3
      erextot(13) =  56.2796E3
      erextot(14) =  59.2020E3
      erextot(15) =  61.7277E3
      erextot(16) =  63.8663E3
      erextot(17) =  65.6122E3
      erextot(18) =  66.9733E3
      erextot(19) =  67.9598E3
      erextot(20) =  68.5916E3
      erextot(21) =  68.8917E3
      erextot(22) =  68.8837E3
      erextot(23) =  68.5997E3
      erextot(24) =  68.0707E3
      erextot(25) =  67.3394E3
      erextot(26) =  66.4367E3
      erextot(27) =  65.4149E3
      erextot(28) =  64.3005E3
      erextot(29) =  63.1456E3
      erextot(30) =  61.9937E3
      erextot(31) =  60.8426E3
      erextot(32) =  59.7278E3
      erextot(33) =  58.6216E3
      erextot(34) =  57.5529E3
      erextot(35) =  56.4092E3
      erextot(36) =  55.1953E3
      erextot(37) =  54.0807E3
      erextot(38) =  52.7824E3
      erextot(39) =  51.1412E3
      erextot(40) =  48.9597E3
      erextot(41) =  46.1144E3
      erextot(42) =  42.1548E3
      erextot(43) =  36.7260E3
      erextot(44) =  29.5493E3
      erextot(45) =  20.5774E3
      erextot(46) =  10.6846E3
      erextot(47) =  -0.3236E3
      erextot(48) = -16.3703E3
      erextot(49) = -35.9385E3
      erextot(50) = -21.8190E3
      erextot(51) = -12.8250E3
c	rh sign convention + is down at outer midplane
	vthexp(1) = 1.0e3
	vthexp(2) = 1.0e3
	vthexp(3) = -1.53e3
	vthexp(4) = -1.67e3
	vthexp(5) = -1.85e3
	vthexp(6) = -2.08e3
	vthexp(7) = -2.25e3
	vthexp(8) = -2.46e3
	vthexp(9) = -2.66e3
	vthexp(10) =-2.83e3
	vthexp(11) =-3.08e3
	vthexp(12) =-3.28e3
	vthexp(13) =-3.47e3
	vthexp(14) =-3.71e3
	vthexp(15) =-3.95e3
	vthexp(16) =-4.16e3
	vthexp(17) =-4.41e3
	vthexp(18) =-4.47e3
	vthexp(19) =-4.36e3
	vthexp(20) =-3.65e3
	vthexp(21) =-3.15e3
	vthexp(22) =-2.15e3
	vthexp(23) =-0.70e3
	vthexp(24) = 0.50e3
	vthexp(25) = 1.0e3

	torv(1) =1.0e3
	torv(2) =1.0e3
	torv(3)	 =68.5e3 
	torv(4)  =66.0e3 
	torv(5)  =62.0e3 
	torv(6)  =58.5e3 
	torv(7)  =54.5e3 
	torv(8)  =50.5e3 
	torv(9)  =46.2e3 
	torv(10) =42.2e3 
	torv(11) =37.5e3 
	torv(12) =33.0e3 
	torv(13) =28.5e3 
	torv(14) =22.1e3 
	torv(15) =18.0e3 
	torv(16) =14.1e3 
	torv(17) =11.4e3 
	torv(18) =9.13e3 
	torv(19) =7.30e3
	torv(20) =6.35e3
	torv(21) =6.16e3
	torv(22) =6.70e3
	torv(23) =7.80e3
	torv(24) =9.20e3 
	torv(25) = 1.0e3

c     Neutral beam deposition profile (hofr1) from NBeams (Dr. John Mandrekas)
c     units [/s] - need to convert to GTEDGE rho coordinates using MATLAB script
c     nbeams2gtedge.m	

      NBdep1(1) =   5.0500
      NBdep1(2) =  10.5360
      NBdep1(3) =   6.9943
      NBdep1(4) =   6.2466
      NBdep1(5) =   5.8823
      NBdep1(6) =   5.6315
      NBdep1(7) =   5.4262
      NBdep1(8) =   5.2331
      NBdep1(9) =   5.0324
      NBdep1(10) =   4.8226
      NBdep1(11) =   4.5918
      NBdep1(12) =   4.3389
      NBdep1(13) =   4.0471
      NBdep1(14) =   3.6857
      NBdep1(15) =   3.2818
      NBdep1(16) =   3.0248
      NBdep1(17) =   2.8124
      NBdep1(18) =   2.6286
      NBdep1(19) =   2.4618
      NBdep1(20) =   2.3025
      NBdep1(21) =   2.1523
      NBdep1(22) =   2.0077
      NBdep1(23) =   1.8694
      NBdep1(24) =   1.7324
      NBdep1(25) =   1.6027
      NBdep1(26) =   1.4759
      NBdep1(27) =   1.3525
      NBdep1(28) =   1.2345
      NBdep1(29) =   1.1215
      NBdep1(30) =   1.0140
      NBdep1(31) =   0.9128
      NBdep1(32) =   0.8180
      NBdep1(33) =   0.7299
      NBdep1(34) =   0.6476
      NBdep1(35) =   0.5731
      NBdep1(36) =   0.5051
      NBdep1(37) =   0.4444
      NBdep1(38) =   0.3904
      NBdep1(39) =   0.3429
      NBdep1(40) =   0.3015
      NBdep1(41) =   0.2662
      NBdep1(42) =   0.2363
      NBdep1(43) =   0.2115
      NBdep1(44) =   0.1913
      NBdep1(45) =   0.1750
      NBdep1(46) =   0.1622
      NBdep1(47) =   0.1523
      NBdep1(48) =   0.1447
      NBdep1(49) =   0.1383
      NBdep1(50) =   0.1329
      NBdep1(51) =   0.0000

      NBdep2(1) =   4.0115
      NBdep2(2) =   8.3782
      NBdep2(3) =   5.5770
      NBdep2(4) =   5.0029
      NBdep2(5) =   4.7402
      NBdep2(6) =   4.5729
      NBdep2(7) =   4.4468
      NBdep2(8) =   4.3346
      NBdep2(9) =   4.2174
      NBdep2(10) =   4.0942
      NBdep2(11) =   3.9532
      NBdep2(12) =   3.7921
      NBdep2(13) =   3.5939
      NBdep2(14) =   3.3327
      NBdep2(15) =   3.0271
      NBdep2(16) =   2.8384
      NBdep2(17) =   2.6813
      NBdep2(18) =   2.5442
      NBdep2(19) =   2.4174
      NBdep2(20) =   2.2915
      NBdep2(21) =   2.1697
      NBdep2(22) =   2.0487
      NBdep2(23) =   1.9301
      NBdep2(24) =   1.8081
      NBdep2(25) =   1.6907
      NBdep2(26) =   1.5726
      NBdep2(27) =   1.4548
      NBdep2(28) =   1.3399
      NBdep2(29) =   1.2277
      NBdep2(30) =   1.1190
      NBdep2(31) =   1.0151
      NBdep2(32) =   0.9165
      NBdep2(33) =   0.8237
      NBdep2(34) =   0.7361
      NBdep2(35) =   0.6556
      NBdep2(36) =   0.5819
      NBdep2(37) =   0.5152
      NBdep2(38) =   0.4555
      NBdep2(39) =   0.4027
      NBdep2(40) =   0.3565
      NBdep2(41) =   0.3170
      NBdep2(42) =   0.2835
      NBdep2(43) =   0.2558
      NBdep2(44) =   0.2334
      NBdep2(45) =   0.2156
      NBdep2(46) =   0.2020
      NBdep2(47) =   0.1916
      NBdep2(48) =   0.1836
      NBdep2(49) =   0.1768
      NBdep2(50) =   0.1708
      NBdep2(51) =   0.0000

      NBdep3(1) =   3.6012
      NBdep3(2) =   7.5246
      NBdep3(3) =   5.0148
      NBdep3(4) =   4.5071
      NBdep3(5) =   4.2818
      NBdep3(6) =   4.1445
      NBdep3(7) =   4.0461
      NBdep3(8) =   3.9622
      NBdep3(9) =   3.8744
      NBdep3(10) =   3.7820
      NBdep3(11) =   3.6736
      NBdep3(12) =   3.5465
      NBdep3(13) =   3.3839
      NBdep3(14) =   3.1620
      NBdep3(15) =   2.8964
      NBdep3(16) =   2.7356
      NBdep3(17) =   2.6015
      NBdep3(18) =   2.4843
      NBdep3(19) =   2.3749
      NBdep3(20) =   2.2641
      NBdep3(21) =   2.1556
      NBdep3(22) =   2.0462
      NBdep3(23) =   1.9377
      NBdep3(24) =   1.8241
      NBdep3(25) =   1.7138
      NBdep3(26) =   1.6014
      NBdep3(27) =   1.4881
      NBdep3(28) =   1.3765
      NBdep3(29) =   1.2665
      NBdep3(30) =   1.1591
      NBdep3(31) =   1.0558
      NBdep3(32) =   0.9570
      NBdep3(33) =   0.8636
      NBdep3(34) =   0.7749
      NBdep3(35) =   0.6929
      NBdep3(36) =   0.6176
      NBdep3(37) =   0.5492
      NBdep3(38) =   0.4877
      NBdep3(39) =   0.4333
      NBdep3(40) =   0.3855
      NBdep3(41) =   0.3446
      NBdep3(42) =   0.3101
      NBdep3(43) =   0.2816
      NBdep3(44) =   0.2587
      NBdep3(45) =   0.2409
      NBdep3(46) =   0.2275
      NBdep3(47) =   0.2176
      NBdep3(48) =   0.2099
      NBdep3(49) =   0.2030
      NBdep3(50) =   0.1966
      NBdep3(51) =   0.0000

c      cosine of angle between injected NB particle and Bphi
      fpsi0(1) =  -0.4459
      fpsi0(2) =  -0.4429
      fpsi0(3) =  -0.4399
      fpsi0(4) =  -0.4370
      fpsi0(5) =  -0.4341
      fpsi0(6) =  -0.4311
      fpsi0(7) =  -0.4285
      fpsi0(8) =  -0.4256
      fpsi0(9) =  -0.4229
      fpsi0(10) =  -0.4201
      fpsi0(11) =  -0.4175
      fpsi0(12) =  -0.4149
      fpsi0(13) =  -0.4123
      fpsi0(14) =  -0.4096
      fpsi0(15) =  -0.4071
      fpsi0(16) =  -0.4047
      fpsi0(17) =  -0.4021
      fpsi0(18) =  -0.3997
      fpsi0(19) =  -0.3973
      fpsi0(20) =  -0.3949
      fpsi0(21) =  -0.3925
      fpsi0(22) =  -0.3901
      fpsi0(23) =  -0.3878
      fpsi0(24) =  -0.3856
      fpsi0(25) =  -0.3833
      fpsi0(26) =  -0.3811
      fpsi0(27) =  -0.3788
      fpsi0(28) =  -0.3767
      fpsi0(29) =  -0.3746
      fpsi0(30) =  -0.3723
      fpsi0(31) =  -0.3702
      fpsi0(32) =  -0.3683
      fpsi0(33) =  -0.3662
      fpsi0(34) =  -0.3641
      fpsi0(35) =  -0.3621
      fpsi0(36) =  -0.3602
      fpsi0(37) =  -0.3582
      fpsi0(38) =  -0.3562
      fpsi0(39) =  -0.3543
      fpsi0(40) =  -0.3523
      fpsi0(41) =  -0.3505
      fpsi0(42) =  -0.3486
      fpsi0(43) =  -0.3468
      fpsi0(44) =  -0.3450
      fpsi0(45) =  -0.3432
      fpsi0(46) =  -0.3414
      fpsi0(47) =  -0.3396
      fpsi0(48) =  -0.3379
      fpsi0(49) =  -0.3361
      fpsi0(50) =  -0.3345
      fpsi0(51) =  -0.3327


c	123302 @ 2600ms	 H-mode  80-99%	  NEW DATA 8-4-10
	rmajor = 1.746
	aminor = .599
	elong = 1.836
	triang = 0.374
	plasmacur = 1.50
   	B = 1.98
	bphi = -1.98
	q95 = 3.86
	pbeam = 8.66
  	betag = 0.25

	Rx = 1.447
	zx = - 1.223
	Rsep1 = 1.589
	rsep2 = 1.215
	zsep1 =-1.367
	zsep2 =-1.367
	ssi95 = 4.49
	pohmin = 0.303
  	

	fzintrin = 0.03
	cfzint =   0.03
	cfzinttb = 0.03


 	fuelpf = 0.0e21
 	fuelplout = 0.0e19
	fuelplin = 0.0
	fuelmp = 0.0


	cballoon = 2.2
	hconf = 1.05
	tauratio =0.5
	hrat = 1.0
	alphan =3.0
	alphat2 = 3.5
	delxpt = 0.1
	delxreal = 0.05


	
   	fheate = 0.4
      xnpedex = 3.84e19
	xnsepex = 0.75e19
	tpedexi = 1075.
	tsepexi =  498.
	tpedexe = 1000.
	tsepexe = 91.
	widthnx = .066
	widthtex = .056
	widthtix = .076
	gradnbar(3) = .049
	gradnbar(4) = .049
	gradTbar (3) = .076
 	gradTbar (4) = .076
	gradTebar(3) = .034
	gradTebar(4) = .034
	aped = 0.86
	xnctrped = 2.3
	tctrped	 = 5.0
	 

c	*************reverse vthet sign to make consistent with rhs convention used in code*******
c	do 408 n = 1,25
c	vthexp(n) = -1.*vthexp(n)
c408	continue	 




c     *************necessary to copy all 8 files from data19 above*****
c     ****necessary to copy 3 files from data20 below*******************
c	118897 @ 2140
			dlnn_dt(1) =       6.178e-1
			dlnn_dt(2) =       6.243e-1
			dlnn_dt(3) =       6.313e-1
			dlnn_dt(4) =       6.352e-1
			dlnn_dt(5) =       6.341e-1
			dlnn_dt(6) =       6.309e-1
			dlnn_dt(7) =       6.275e-1
			dlnn_dt(8) =       6.228e-1
			dlnn_dt(9) =       6.173e-1
			dlnn_dt(10) =      6.108e-1
			 dlnn_dt(11) =     6.169e-1
			dlnn_dt(12) =      6.210e-1
			dlnn_dt(13) =      5.994e-1
			dlnn_dt(14) =      5.449e-1
			dlnn_dt(15) =      4.538e-1
			dlnn_dt(16) =      3.103e-1
			dlnn_dt(17) =      1.115e-1
			 dlnn_dt(18) =    -8.607e-2
			 dlnn_dt(19) =    -4.165e-1
			dlnn_dt(20) =     -7.858e-1
			dlnn_dt(21) =     -9.755e-1
			dlnn_dt(22) =     -1.042e+0
			dlnn_dt(23) =     -1.059e+0
			dlnn_dt(24) =     -1.069e+0
			dlnn_dt(25) =     -1.079e+0


			 dlnwe_dt(1) =     6.876e-1
			 dlnwe_dt(2) =     6.794e-1
			 dlnwe_dt(3) =     6.670e-1
			 dlnwe_dt(4) =     6.481e-1
			 dlnwe_dt(5) =     6.208e-1
			 dlnwe_dt(6) =     5.872e-1
			 dlnwe_dt(7) =     5.544e-1
			 dlnwe_dt(8) =     5.080e-1
			  dlnwe_dt(9) =    4.522e-1
			dlnwe_dt(10) =     3.857e-1
			 dlnwe_dt(11) =    3.188e-1
			 dlnwe_dt(12) =    2.515e-1
			dlnwe_dt(13) =     1.343e-1
			dlnwe_dt(14) =    -2.072e-2
			dlnwe_dt(15) =    -2.004e-1
			dlnwe_dt(16) =    -3.920e-1
			 dlnwe_dt(17) =   -5.642e-1
			dlnwe_dt(18) =    -6.693e-1
			dlnwe_dt(19) =    -8.212e-1
			dlnwe_dt(20) =    -1.038e+0
			 dlnwe_dt(21) =   -1.215e+0
			 dlnwe_dt(22) =   -1.488e+0
			 dlnwe_dt(23) =   -1.866e+0
			 dlnwe_dt(24) =   -1.463e+0
			 dlnwe_dt(25) =   -3.041e+0

			dlnwi_dt(1) =         3.169e-1
			dlnwi_dt(2) =         3.247e-1
			dlnwi_dt(3) =         3.352e-1
			dlnwi_dt(4) =         3.444e-1
			dlnwi_dt(5) =         3.508e-1
			 dlnwi_dt(6) =        3.577e-1
			dlnwi_dt(7) =         3.649e-1
			dlnwi_dt(8) =         3.757e-1
			dlnwi_dt(9) =         3.892e-1
			dlnwi_dt(10) =        4.058e-1
			dlnwi_dt(11) =        4.401e-1
			dlnwi_dt(12) =        4.712e-1
			dlnwi_dt(13) =        4.861e-1
			dlnwi_dt(14) =        4.733e-1
			dlnwi_dt(15) =        4.308e-1
			 dlnwi_dt(16) =       3.438e-1
			dlnwi_dt(17) =        2.114e-1
			 dlnwi_dt(18) =       7.712e-2
			dlnwi_dt(19) =       -1.813e-1
			dlnwi_dt(20) =       -4.780e-1
			dlnwi_dt(21) =       -5.708e-1
			dlnwi_dt(22) =       -5.229e-1
			dlnwi_dt(23) =       -4.343e-1
			dlnwi_dt(24) =       -3.002e-1
			dlnwi_dt(25) =       -1.475e-1


			
c	*******************temp*****************
	do 1044 n = 1,25
	dlnwi_dt(n)=0.
	dlnwe_dt(n)=0.
	dlnn_dt(n) = 0.
1044	continue		     
c	***********************
		 goto 225
	write(*,*) delna
		 rhor(25) = 1.0
      	 dO 105, NN=1,24
		  n= 25-NN
		 rhor(n) = rhor(n+1) - delna/(aminor*SQRT(0.5*(1.+ELONG**2)))
 105		 continue 
c	     do 200 n=1,25
c	     xlne(n) = 0.5*(exne(n+1)+exne(n))*delna/(exne(n)-exne(n+1)) 
c	     exlte(n) = 0.5*(xte(n+1)+xte(n))*delna/(xte(n)-xte(n+1)) 
c	     exlti(n) = 0.5*(xti(n+1)+xti(n))*delna/(xti(n)-xti(n+1))
c	     if(n.eq.25) then
c	     xlne(n)=xlne(n-1)
c	     exlte(n)=exlte(n-1)
c	     exlti(n)=exlti(n-1)
c  	endif
c200	continue 
c	normalize Rich's gsl to FSA geometry
225	se = sqrt(0.5*(1.+(elong**2)))
	do 250 n=1,26
	xlne(n) = se*xlne(n)
	exlte(n) = se*exlte(n)
	exlti(n) = se*exlti(n)	    	
c	qedge(n) = q95
250	continue
c	save experimental impurity poloidal velocity
	rhor(25) = 1.0
      	 do 255, NN=1,24
		  n= 25-NN
		 rhor(n) = rhor(n+1) - delna/(aminor*SQRT(0.5*(1.+ELONG**2)))
255		 continue 
 
	do 275 n = 1, 25
	vpol_imp(n) = vthexp(n)
275	continue
	do 1150 n = 1,24
	shearho(n) =  (qedge(n+1)-qedge(n))/(rhor(n+1)-rhor(n))
	dEdredge(n) = (erex(n+1)-erex(n))/(rhor(n+1)-rhor(n))
1150	CONTINUE
	shearho(25) = shearho(24) 
	dEdredge(25) = dEdredge(24)
c	smooth derivatives
	do 1155 n =2,24
	shearhonew(n) = (shearho(n-1)+shearho(n)+shearho(n+1))/3.
	dEdrnew(n)= (dEdredge(n-1)+dEdredge(n)+dEdredge(n+1))/3.
1155	continue
 	shearhonew(1)=shearho(1)
	shearhonew(25) = shearho(25)
	dEdrnew(1) = dEdredge(1)
	dEdrnew(25) = dEdredge(25)
	do 1160 n=1,25
	shearho(n) = shearhonew(n)
	dEdredge(n) = Dedrnew(n)
1160	continue 
	goto 300 
c	input stochastic magnetic diffusion coefficient
	dmag(1) = .47e-7
	dmag(2) = .32e-7
	dmag(3) = .30e-7
	dmag(4) = .43e-7
	dmag(5) = .50e-7
	dmag(6) = .48e-7
	dmag(7) = .46e-7
	dmag(8) = .46e-7
	dmag(9) = .46e-7
	dmag(10)= .61e-7
	dmag(11)= .72e-7
	dmag(12)= 1.5e-7
	dmag(13)= 2.0e-7
	dmag(14)= 1.9e-7
	dmag(15)= 1.7e-7
	dmag(16)= 2.4e-7
	dmag(17)= 2.8e-7
	dmag(18)= 2.6e-7
	dmag(19)= 2.4e-7
	dmag(20)= 2.2e-7
	dmag(21)= 2.1e-7
	dmag(22)= 2.0e-7
	dmag(23)= 1.8e-7
	dmag(24)= 0.0
	dmag(25)= 0.0
	 


300   return
	end   
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	