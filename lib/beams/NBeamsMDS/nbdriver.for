c/    This program is a driver to test the NB module nbeams. It is not
c/    part of the nbeams package. Its purpose is to create a background
c/    plasma (geometry, MHD and plasma profiles) and then call the NB 
c/    routines to calculate various neutral beam heating and current drive 
c/    parameters.
c/
c/    The program reads the input file 'inbeams.dat' which contains a 
c/    namelist with various input variables, calls the main routine of
c/    the NB module (calcBeams) and prints out some information to the 
c/    file 'outbeams.dat'
c/
c/    Written by John Mandrekas, GIT, for the NTCC 08/29/00
c/
c/    Description of some of the input variables in the namelist
c/    (Only the local variables are included here. For a description
c/    of the variables in the calcBeams list see the extensive comments
c/    in the beginning of the calcbeams.f routine)
c/
c/    e0         - elongation at center
c/    ea         - elongation at edge									 
c/    shft0      - shift of magnetic axis (normalized to minor radius)  
c/    teavg      - Density-weighted volume average electron temperature (keV)
c/    tea        - Electron temperature at the edge (separatrix) (keV)
c/    tiavg      - Density-weighted volume average ion temperature (keV)
c/    tia        - Ion temperature at the edge (separatrix) (keV)
c/    alfat      - Temperature profile exponent
c/    denav      - Average electron density (10^20 /m^3)
c/    edgavr     - edge-to-average ratio, ne(a)/<ne>
c/    alfan      - Density profile exponent
c/    fion       - Concentration of each ion species (ni/ne)
         
c/
      implicit none
      include 'nbparams.inc'
      include 'BeamsLoc.inc'

c/    Variables in the calcBeams variable list:
c/    ----------------------------------------
      
	integer iflag, ie, ib, nbeams, inbfus, n, nion, maxiter
      integer nbshape(maxBeams), nbptype(maxBeams)

      real amb, zbeam, r0, a, b0, volp, nbcurTot, etanbTot, beamBeta,
     .     pNBAbsorbTot, pNBLossTot, beamFusTot, beamFusChTot,
     .     snDDTotal, snDTTotal, taus

      real ebeam(maxBeams), pbeam(maxBeams), rtang(maxBeams), 
     .     bwidth(maxBeams), bheigh(maxBeams), bgaussR(maxBeams), 
     .     bgaussZ(maxBeams), bzpos(maxBeams), pwrfrac(3,maxBeams), 
     .     aion(maxIons), zion(maxIons), ne20(mxrho), 
     .     ni20(mxrho,maxIons), tekev(mxrho), tikev(mxrho), zeff(mxrho),
     .     rnorm(mxrho), vprime(mxrho), dvol(mxrho), darea(mxrho),
     .     kappa(mxrho), dkappa(mxrho), shafr(mxrho), dshafr(mxrho), 
     .     hofr(mxrho,3,maxBeams), shinethru(3,maxBeams), jnbTot(mxrho),
     .     pnbe(mxrho), pnbi(mxrho), beamDens(mxrho), beamPress(mxrho),
     .     pbfuse(mxrho), pbfusi(mxrho), snBeamDD(mxrho),beamFus(mxrho),
     .     snBeamDT(mxrho), nbcur(maxBeams), etanb(maxBeams), 
     .     gammanb(maxBeams), pNBAbsorb(maxBeams), pNBLoss(maxBeams)
	 
	

c/    Local variables:
      
      integer i, j, ip1, nm1, nbinp, nbout
      real e0, te0, ti0, ea, shft0, dlr, tt, pi,
	.   den0, L_tor, sumzef, rn2, ssum

      real r(mxrho), rmid(mxrho)

c      namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape, 
c     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
c     .    bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0,
c     .    e0, ea, shft0, teavg, tiavg, tea, tia, alfat, denav, edgavr, 
c     .    alfan, fion

      namelist /nbin/ nbeams, inbfus, n, nion, maxiter, nbshape, 
     .    nbptype, amb, zbeam, ebeam, pbeam, rtang, bwidth, bheigh,
     .    bgaussR, bgaussZ, bzpos, pwrfrac, aion, zion, r0, a, b0,
     .    e0, ea, shft0

      pi = acos(-1.0)

c/    Assign I/O unit numbers:
c/    -----------------------
      nbinp = 10
      nbout = 20

c/    Open files:
c/    ----------
      open (unit=nbinp,file ='inbeams.dat',status ='old')
      read (nbinp,nbin)
	close(unit=nbinp,status='keep')

	open (nbout,file ='outbeams.dat',status ='unknown')
      
c/    BACKGROUND PLASMA:

c/    Grid and other geometric quantities:
c/    -----------------------------------
      nm1 = n - 1
      dlr = a / float(nm1)
      r(1) = 0.0
	r(n) = a
	open(unit=30,file='bugging.txt',status='old')
      do i = 2, nm1
          r(i) = r(i-1) + dlr
	    write(30,*) r(i)
	enddo
	close(30,status='keep')
	
c/    Interface grid:
c/    --------------
c/    Notice that i+1/2 corresponds to i+1 in the interface grid!
      rmid(1) = 0.0
      do i = 2, n
        rmid(i) = 0.5 * (r(i-1) + r(i))
      enddo

c/    Neutral Beam normalized grid:

      do i = 1, n
         rnorm(i) = r(i) / r(n)
         rn2 = rnorm(i)**2
         tt = 1.0 - rn2
         kappa(i)  = e0 + (ea - e0)*rn2                                         
         dkappa(i) = 2.0*(ea - e0) * rnorm(i)                                
         shafr(i)  = shft0 * tt
         dshafr(i) = -2.0 * shft0 * rnorm(i)
      enddo
         
c/    Cylindrical limit of  important metric  quantities:
c/    --------------------------------------------------

      L_tor = 2.0 * pi * r0
      vprime(1) = 0.0
      do i = 1, nm1
	ip1 = i + 1
        dvol(i) = L_tor * pi * kappa(i) * (rmid(ip1)**2 - rmid(i)**2)
	vprime(ip1) = L_tor * 2.0 * pi * kappa(i) * rmid(ip1)
	darea(i) = dvol(i) / L_tor 
      enddo

      dvol(n) = 0.0
      darea(n) = 0.0
	
c/    Plasma volumes:
c/    --------------
      volp = ssum(nm1, dvol, 1)
c/    Parabolic-like plasma profiles:
c/    ------------------------------
c      te0 = ((1.0 + alfat+alfan) / (1.0 + alfan - alfan*alfat*edgavr /
c     .   (1.0 + alfat))) * teavg
c      ti0 = te0 * tiavg / teavg
c      den0 = (1.0 + alfan - alfan*edgavr) * denav
	   tekev(1)   =  1.8774444
	   tekev(2)   =  1.8745041
	   tekev(3)   =  1.8710290
	   tekev(4)   =  1.8664843
	   tekev(5)   =  1.8603349
	   tekev(6)   =  1.8520472
	   tekev(7)   =  1.8412869
	   tekev(8)   =  1.8281659
	   tekev(9)   =  1.8130698
	   tekev(10)  =  1.7965096
	   tekev(11)  =  1.7778015
	   tekev(12)  =  1.7576991
	   tekev(13)  =  1.7356207
	   tekev(14)  =  1.7121780
	   tekev(15)  =  1.6870150
	   tekev(16)  =  1.6603658
	   tekev(17)  =  1.6323369
	   tekev(18)  =  1.6028997
	   tekev(19)  =  1.5721576
	   tekev(20)  =  1.5402097
	   tekev(21)  =  1.5071616
	   tekev(22)  =  1.4730589
	   tekev(23)  =  1.4379949
	   tekev(24)  =  1.4020641
	   tekev(25)  =  1.3653637
	   tekev(26)  =  1.3279943
	   tekev(27)  =  1.2900660
 	   tekev(28)  =  1.2516791
	   tekev(29)  =  1.2129402
	   tekev(30)  =  1.1739659
	   tekev(31)  =  1.1348767
	   tekev(32)  =  1.0957901
	   tekev(33)  =  1.0568345
	   tekev(34)  =  1.0181310
	   tekev(35)  = 0.97982222
	   tekev(36)  = 0.94202911
	   tekev(37)  = 0.90490506
	   tekev(38)  = 0.86858323
	   tekev(39)  = 0.83320932
	   tekev(40)  = 0.79894456
	   tekev(41)  = 0.76592445
	   tekev(42)  = 0.73431212
	   tekev(43)  = 0.70428307
	   tekev(44)  = 0.67595631
	   tekev(45)  = 0.64947884
	   tekev(46)  = 0.62468917
	   tekev(47)  = 0.59987883
	   tekev(48)  = 0.56434535
	   tekev(49)  = 0.47177339
	   tekev(50)  = 0.26508076
	   tekev(51)  = 0.10241394
	   tikev(1)   =   2.4168516
	   tikev(2)   =   2.4108240
	   tikev(3)   =   2.4028383
	   tikev(4)   =   2.3909365
	   tikev(5)   =   2.3732991
	   tikev(6)   =   2.3498187
	   tikev(7)   =   2.3222266
	   tikev(8)   =   2.2907560
	   tikev(9)   =   2.2560479
	   tikev(10)  =   2.2180910
	   tikev(11)  =   2.1773615
	   tikev(12)  =   2.1343294
	   tikev(13)  =   2.0889563
	   tikev(14)  =   2.0413928
	   tikev(15)  =   1.9917738
	   tikev(16)  =   1.9402262
	   tikev(17)  =   1.8869721
	   tikev(18)  =   1.8324005
	   tikev(19)  =   1.7769511
	   tikev(20)  =   1.7210120
	   tikev(21)  =   1.6650125
	   tikev(22)  =   1.6093757
	   tikev(23)  =   1.5545051
	   tikev(24)  =   1.5008487
	   tikev(25)  =   1.4487864
	   tikev(26)  =   1.3987655
	   tikev(27)  =   1.3511485
	   tikev(28)  =   1.3057930
	   tikev(29)  =   1.2625814
	   tikev(30)  =   1.2213038
	   tikev(31)  =   1.1817080
	   tikev(32)  =   1.1436440
	   tikev(33)  =   1.1069128
	   tikev(34)  =   1.0712950
	   tikev(35)  =   1.0365940
	   tikev(36)  =   1.0026270
	   tikev(37)  =  0.96918549
	   tikev(38)  =  0.93606693
	   tikev(39)  =  0.90308041
	   tikev(40)  =  0.87002404
	   tikev(41)  =  0.83669700
	   tikev(42)  =  0.80290909
	   tikev(43)  =  0.76845057
	   tikev(44)  =  0.73317842
	   tikev(45)  =  0.69708380
	   tikev(46)  =  0.66017821
	   tikev(47)  =  0.62247072
	   tikev(48)  =  0.58397345
	   tikev(49)  =  0.54470150
	   tikev(50)  =  0.50466495
	   tikev(51)  =  0.46387507
	   ne20(1)    =  0.91754780
	   ne20(2)    =  0.91447238
	   ne20(3)    =  0.90994650
	   ne20(4)    =  0.90253702
	   ne20(5)    =  0.89173409
	   ne20(6)    =  0.87870479
	   ne20(7)    =  0.86373714
	   ne20(8)    =  0.84725101
	   ne20(9)    =  0.82953817
	   ne20(10)   =  0.81084922
	   ne20(11)   =  0.79146844
	   ne20(12)   =  0.77165633
	   ne20(13)   =  0.75165426
	   ne20(14)   =  0.73168674
	   ne20(15)   =  0.71196090
	   ne20(16)   =  0.69266090
	   ne20(17)   =  0.67395802
	   ne20(18)   =  0.65601118
	   ne20(19)   =  0.63895886
	   ne20(20)   =  0.62292214
	   ne20(21)   =  0.60801058
	   ne20(22)   =  0.59431113
	   ne20(23)   =  0.58183836
	   ne20(24)   =  0.57068970
	   ne20(25)   =  0.56098043
	   ne20(26)   =  0.55250098
	   ne20(27)   =  0.54566670
	   ne20(28)   =  0.53982182
	   ne20(29)   =  0.53546828
	   ne20(30)   =  0.53248644
	   ne20(31)   =  0.53075668
	   ne20(32)   =  0.53015938
	   ne20(33)   =  0.53057493
	   ne20(34)   =  0.53188370
	   ne20(35)   =  0.53378859
	   ne20(36)   =  0.53656768
	   ne20(37)   =  0.53984002
	   ne20(38)   =  0.54311047
	   ne20(39)   =  0.54657133
	   ne20(40)   =  0.55017650
	   ne20(41)   =  0.55342035
	   ne20(42)   =  0.55586575
	   ne20(43)   =  0.55707560
	   ne20(44)   =  0.55661280
	   ne20(45)   =  0.55446579
	   ne20(46)   =  0.55074124
	   ne20(47)   =  0.54161816
	   ne20(48)   =  0.51716057
	   ne20(49)   =  0.43905650
	   ne20(50)   =  0.26259143
	   ne20(51)   =  0.10760898
	   ni20(1,1)  =  0.65799831
	   ni20(2,1)  =  0.65668160
	   ni20(3,1)  =  0.65468712
	   ni20(4,1)  =  0.65134295
	   ni20(5,1)  =  0.64629127
	   ni20(6,1)  =  0.63987830
	   ni20(7,1)  =  0.63207199
	   ni20(8,1)  =  0.62293308
	   ni20(9,1)  =  0.61247583
	   ni20(10,1) =  0.60066936
	   ni20(11,1) =  0.58757461
	   ni20(12,1) =  0.57329082
	   ni20(13,1) =  0.55790858
	   ni20(14,1) =  0.54151229
	   ni20(15,1) =  0.52421353
	   ni20(16,1) =  0.50609627
	   ni20(17,1) =  0.48721995
	   ni20(18,1) =  0.46772024
	   ni20(19,1) =  0.44812463
	   ni20(20,1) =  0.42907810
	   ni20(21,1) =  0.41116441
	   ni20(22,1) =  0.39489219
	   ni20(23,1) =  0.38060189
	   ni20(24,1) =  0.36886853
	   ni20(25,1) =  0.36005038
	   ni20(26,1) =  0.35395032
	   ni20(27,1) =  0.35022447
	   ni20(28,1) =  0.34723882
	   ni20(29,1) =  0.34577815
	   ni20(30,1) =  0.34555507
	   ni20(31,1) =  0.34628219
	   ni20(32,1) =  0.34767209
	   ni20(33,1) =  0.34943741
	   ni20(34,1) =  0.35129074
	   ni20(35,1) =  0.35210076
	   ni20(36,1) =  0.35293319
	   ni20(37,1) =  0.35436238
	   ni20(38,1) =  0.35666430
	   ni20(39,1) =  0.35936946
	   ni20(40,1) =  0.36276587
	   ni20(41,1) =  0.36712907
	   ni20(42,1) =  0.37259996
	   ni20(43,1) =  0.37931946
	   ni20(44,1) =  0.38742854
	   ni20(45,1) =  0.39670221
	   ni20(46,1) =  0.40704308
	   ni20(47,1) =  0.42212705
	   ni20(48,1) =  0.43871499
	   ni20(49,1) =  0.40665652
	   ni20(50,1) =  0.25244577
	   ni20(51,1) =  0.10371023
	  te0 = tekev(1)
	  ti0 = tikev(1)
	  den0 = ne20(1)

      do  i = 1, n
	  tt = 1.0 - (r(i) / a)**2
c	  tekev(i) = (te0 - tea) * tt**alfat + tea
c        tikev(i) = (ti0 - tia) * tt**alfat + tia
c        ne20(i) = (den0 - edgavr*denav) * tt**alfan + 
c     .     edgavr*denav
	do j = 2, nion
	   ni20(i,j) = (ne20(i) - ni20(i,j-1)*zion(j-1))/zion(j)
	enddo
      enddo
	
c/    Calculation of the plasma Zeff:
c/    ------------------------------
      
      do i = 1, n
         sumzef = 0.0
         do j = 1, nion
            sumzef = sumzef + ni20(i,j) * zion(j)**2 / ne20(i)
         enddo
         zeff(i) = sumzef
      enddo

c/    Call the beams package:

      call calcBeams(nbeams, amb, zbeam, ebeam, pbeam, inbfus,      
     .   rtang, nbshape, bwidth, bheigh, nbptype, bgaussR, bgaussZ,       
     .   bzpos, pwrfrac, maxiter, nion, aion, zion, ne20, ni20, tekev,    
     .   tikev, zeff, r0, a, b0, volp, n, rnorm, vprime, dvol, darea,     
     .   kappa, dkappa, shafr, dshafr, hofr, shinethru, jnbTot, pnbe,     
     .   pnbi, beamDens, beamPress, beamFus, pbfuse, pbfusi, snBeamDD,    
     .   snBeamDT, nbcur, etanb, gammanb, pNBAbsorb, pNBLoss, nbcurTot,   
     .   etanbTot, beamBeta, pNBAbsorbTot, pNBLossTot, beamFusTot,        
     .   beamFusChTot, snDTTotal, snDDTotal, iflag,taus)                       
                                                                          
c/    Write output quantities to file nbout:
	
 
c/    Global parameters
c/    -----------------
      write (nbout, 2000)
      write (nbout, 2100)
      write (nbout, 2200) pNBAbsorbTot, pNBLossTot, nbcurTot,
     .   etanbTot, beamBeta, taus, volp
      write (nbout, '(1x)')
      
      write (nbout, 2300)
      write (nbout, 2400)
      write (nbout, 2410) (pNBAbsorb(ib), ib = 1, nbeams)
      write (nbout, 2415) (pNBLoss(ib), ib = 1, nbeams)
      write (nbout, 2420) (nbcur(ib), ib = 1, nbeams)
      write (nbout, 2430) (etanb(ib), ib = 1, nbeams)
      write (nbout, 2440) (gammanb(ib), ib = 1, nbeams)
      write (nbout, '(1x)')
      write (nbout, 2450)
      write (nbout, 2460) (shinethru(1,ib), ib = 1, nbeams)
      write (nbout, 2470) (shinethru(2,ib), ib = 1, nbeams)
      write (nbout, 2480) (shinethru(3,ib), ib = 1, nbeams)
      write (nbout, '(1x)')
      
      if (inbfus.NE.0) then
         write (nbout, 2500)
         write (nbout, 2600) beamFusTot, beamFusChTot, snDTTotal, 
     .      snDDTotal
      endif
 
      write (nbout, '(1x)')

c/    Write the deposition profile for each beamline and energy group:
c/    ---------------------------------------------------------------
      write (nbout, *) 'Neutral Beam Deposition Profiles'
      write (nbout, *) '--------------------------------'
      do ib = 1, nbeams
         write (nbout, 1000) ib
         do i = 1, n
            write (nbout,1100) rnorm(i), (hofr(i,ie,ib),ie=1,3)
         enddo
         write (nbout, '(1x)')
      enddo

c/    Write the pitch angle profile for each beamline and energy group:
c/    ---------------------------------------------------------------
      write (nbout, *) 'Pitch Angle at Birth'
      write (nbout, *) '--------------------------------'
      do ib = 1, nbeams
         write (nbout, 1200) ib
         do i = 1, n
            write (nbout,1300) rnorm(i), (pitchangl(i,ie,ib),ie=1,3)
         enddo
         write (nbout, '(1x)')
      enddo

c/    Write heating and current drive profile information:
c/    ----------------------------------------------------
      write (nbout, 3000)
      do i = 1, n
         write (nbout, 3100) rnorm(i), jnbTot(i), pnbe(i), pnbi(i), 
     .	    beamDens(i), beamPress(i), beamFus(i), dvol(i),darea(i)
      enddo
 
 
 1000 format (1x, 'Beam ', i2/
     .   5x, 'rho', 6x, 'hofr_1', 7x, 'hofr_2', 7x, 'hofr_3')
 1100 format (1x, f8.5, 3(1x,e12.5))
 1200 format (1x, 'Angle ', i2/
     .   5x, 'rho', 6x, 'zeta_1', 7x, 'zeta_2', 7x, 'zeta_3')
 1300 format (1x, f8.5, 3(1x,e12.5))
 
 2000 format (1x, 'Global NB Heating and Current Drive Parameters',/
     .        1x, '----------------------------------------------')
 2100 format (1x, 'TOTALS:')
 2200 format (1x, 'Total Absorbed Power      = ', f9.4, ' MW'/
     .        1x, 'Total Lost Power          = ', f9.4, ' MW'/
     .        1x, 'Total NB Driven Current   = ', e12.5,' A'/
     .        1x, 'Total NBCD Efficiency     = ', e12.5,' A/W'/
     .        1x, 'Total Beam Beta           = ', e12.5, /
     .        1x, 'Taus                      = ',  e12.5, ' s'/
     .        1x, 'Volp                      = ',  e12.5, ' m3' )
     
 2300 format (1x, 'Information per Beamline and Energy Group:'/
     .        1x, '-----------------------------------------')
 2400 format (23x,'Beam 1', 10x, 'Beam 2', 10x, 'Beam 3'/
     .        23x,'------', 10x, '------', 10x, '------')
 2410 format (1x, 'Absorbed Power   ', 1x, e12.5, 2(4x, e12.5))
 2415 format (1x, 'Lost Power       ', 1x, e12.5, 2(4x, e12.5))
 2420 format (1x, 'NB driven current', 1x, e12.5, 2(4x, e12.5))
 2430 format (1x, 'NBCD efficiency  ', 1x, e12.5, 2(4x, e12.5))
 2440 format (1x, 'NBCD gamma       ', 1x, e12.5, 2(4x, e12.5))
 2450 format (1x, 'Beam Shinethrough')
 2460 format (1x, '   energy group 1', 1x, e12.5, 2(4x, e12.5))
 2470 format (1x, '   energy group 2', 1x, e12.5, 2(4x, e12.5))
 2480 format (1x, '   energy group 3', 1x, e12.5, 2(4x, e12.5))
 
 2500 format (1x, 'Beam-Target Parameters'/
     .        1x, '----------------------')
 2600 format (1x, 'Total Beam-Target Fusion Power   = ', e12.5, ' MW'/
     .        1x, 'Total Power to Charged Particles = ', e12.5, ' MW'/
     .        1x, 'Total DT Neutron Rate            = ', e12.5, ' n/s'/
     .        1x, 'Total DD Neutron Rate            = ', e12.5, ' n/s')
     
 3000 format (3x, 'rho', 5x, 'jnbtot', 6x, 'pNBe', 8x, 'pNBi', 8x,
     .   'nbfast', 6x, 'pressb', 6x, 'pfusb', 6x, 'dvol',6x,'dA')
 3100 format (1x, f6.4, 8(1x, e11.4))          
     
	 close(nbout,status='unknown')
      stop 
      end
