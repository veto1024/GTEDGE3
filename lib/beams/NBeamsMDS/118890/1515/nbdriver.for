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
	   tekev(1) = 1.6269347
	   tekev(2) = 1.6219242
	   tekev(3) = 1.6146232
	   tekev(4) = 1.6027477
	   tekev(5) = 1.5850658
	   tekev(6) = 1.5630100
	   tekev(7) = 1.5371928
	   tekev(8) = 1.5082804
	   tekev(9) = 1.4763314
	   tekev(10) = 1.4417496
	   tekev(11) = 1.4048861
	   tekev(12) = 1.3660636
	   tekev(13) = 1.3255970
	   tekev(14) = 1.2837871
	   tekev(15) = 1.2408738
	   tekev(16) = 1.1971330
	   tekev(17) = 1.1528110
	   tekev(18) = 1.1081286
	   tekev(19) = 1.0633001
	   tekev(20) = 1.0185209
	   tekev(21) = 0.97397197
	   tekev(22) = 0.92981856
	   tekev(23) = 0.88620827
	   tekev(24) = 0.84327387
	   tekev(25) = 0.80113223
	   tekev(26) = 0.75988371
	   tekev(27) = 0.71961691
	   tekev(28) = 0.68040146
	   tekev(29) = 0.64229459
	   tekev(30) = 0.60533739
	   tekev(31) = 0.56953858
	   tekev(32) = 0.53492237
	   tekev(33) = 0.50147078
	   tekev(34) = 0.46916669
	   tekev(35) = 0.43796125
	   tekev(36) = 0.40780784
	   tekev(37) = 0.37863215
	   tekev(38) = 0.35034684
	   tekev(39) = 0.32284965
	   tekev(40) = 0.29602332
	   tekev(41) = 0.26973290
	   tekev(42) = 0.24382956
	   tekev(43) = 0.21814762
	   tekev(44) = 0.19250601
	   tekev(45) = 0.16670741
	   tekev(46) = 0.14053433
	   tekev(47) = 0.11371689
	   tekev(48) = 0.085729880
	   tekev(49) = 0.055230791
	   tekev(50) = 0.027278985
	   tekev(51) = 0.020189745
	   tikev(1) = 1.5664581
	   tikev(2) = 1.5609616
	   tikev(3) = 1.5522934
	   tikev(4) = 1.5373974
	   tikev(5) = 1.5155886
	   tikev(6) = 1.4892396
	   tikev(7) = 1.4583461
	   tikev(8) = 1.4235376
	   tikev(9) = 1.3853541
	   tikev(10) = 1.3443145
	   tikev(11) = 1.3010110
	   tikev(12) = 1.2558779
	   tikev(13) = 1.2094811
	   tikev(14) = 1.1623414
	   tikev(15) = 1.1149875
	   tikev(16) = 1.0679473
	   tikev(17) = 1.0217462
	   tikev(18) = 0.97689962
	   tikev(19) = 0.93395110
	   tikev(20) = 0.89343524
	   tikev(21) = 0.85579552
	   tikev(22) = 0.82091785
	   tikev(23) = 0.78831912
	   tikev(24) = 0.75742397
	   tikev(25) = 0.72777393
	   tikev(26) = 0.69882471
	   tikev(27) = 0.67007042
	   tikev(28) = 0.64099650
	   tikev(29) = 0.61138517
	   tikev(30) = 0.58148592
	   tikev(31) = 0.55158339
	   tikev(32) = 0.52196408
	   tikev(33) = 0.49291580
	   tikev(34) = 0.46472896
	   tikev(35) = 0.43768249
	   tikev(36) = 0.41206698
	   tikev(37) = 0.38817820
	   tikev(38) = 0.36628474
	   tikev(39) = 0.34668272
	   tikev(40) = 0.32966037
	   tikev(41) = 0.31512980
	   tikev(42) = 0.30236471
	   tikev(43) = 0.29061517
	   tikev(44) = 0.27908789
	   tikev(45) = 0.26701518
	   tikev(46) = 0.25362235
	   tikev(47) = 0.23812641
	   tikev(48) = 0.21977871
	   tikev(49) = 0.19779257
	   tikev(50) = 0.17139301
	   tikev(51) = 0.13981006
	   ne20(1) = 0.50662315
	   ne20(2) = 0.50499985
	   ne20(3) = 0.50269506
	   ne20(4) = 0.49902882
	   ne20(5) = 0.49364908
	   ne20(6) = 0.48704434
	   ne20(7) = 0.47951396
	   ne20(8) = 0.47123422
	   ne20(9) = 0.46229730
	   ne20(10) = 0.45286351
	   ne20(11) = 0.44306835
	   ne20(12) = 0.43303577
	   ne20(13) = 0.42287703
	   ne20(14) = 0.41269443
	   ne20(15) = 0.40258073
	   ne20(16) = 0.39261879
	   ne20(17) = 0.38288224
	   ne20(18) = 0.37343392
	   ne20(19) = 0.36432841
	   ne20(20) = 0.35560888
	   ne20(21) = 0.34730883
	   ne20(22) = 0.33945374
	   ne20(23) = 0.33205994
	   ne20(24) = 0.32513439
	   ne20(25) = 0.31867636
	   ne20(26) = 0.31265529
	   ne20(27) = 0.30704482
	   ne20(28) = 0.30183568
	   ne20(29) = 0.29696609
	   ne20(30) = 0.29237391
	   ne20(31) = 0.28802163
	   ne20(32) = 0.28380413
	   ne20(33) = 0.27965490
	   ne20(34) = 0.27548112
	   ne20(35) = 0.27117433
	   ne20(36) = 0.26661236
	   ne20(37) = 0.26167891
	   ne20(38) = 0.25625690
	   ne20(39) = 0.25018309
	   ne20(40) = 0.24330966
	   ne20(41) = 0.23548946
	   ne20(42) = 0.22651745
	   ne20(43) = 0.21624227
	   ne20(44) = 0.20446430
	   ne20(45) = 0.19097041
	   ne20(46) = 0.17555136
	   ne20(47) = 0.15789714
	   ne20(48) = 0.13713438
	   ne20(49) = 0.10875013
	   ne20(50) = 0.062792786
	   ne20(51) = 0.032463599
	   ni20(1,1) = 0.47875309
	   ni20(2,1) = 0.47736464
	   ni20(3,1) = 0.47531205
	   ni20(4,1) = 0.47193273
	   ni20(5,1) = 0.46689171
	   ni20(6,1) = 0.46068335
	   ni20(7,1) = 0.45361793
	   ni20(8,1) = 0.44591301
	   ni20(9,1) = 0.43767155
	   ni20(10,1) = 0.42906500
	   ni20(11,1) = 0.42020999
	   ni20(12,1) = 0.41112702
	   ni20(13,1) = 0.40183740
	   ni20(14,1) = 0.39246355
	   ni20(15,1) = 0.38310768
	   ni20(16,1) = 0.37385956
	   ni20(17,1) = 0.36481397
	   ni20(18,1) = 0.35604828
	   ni20(19,1) = 0.34762133
	   ni20(20,1) = 0.33960652
	   ni20(21,1) = 0.33203939
	   ni20(22,1) = 0.32483497
	   ni20(23,1) = 0.31794228
	   ni20(24,1) = 0.31142243
	   ni20(25,1) = 0.30536395
	   ni20(26,1) = 0.29978299
	   ni20(27,1) = 0.29458828
	   ni20(28,1) = 0.28974746
	   ni20(29,1) = 0.28519452
	   ni20(30,1) = 0.28090404
	   ni20(31,1) = 0.27684300
	   ni20(32,1) = 0.27291232
	   ni20(33,1) = 0.26904510
	   ni20(34,1) = 0.26514996
	   ni20(35,1) = 0.26111214
	   ni20(36,1) = 0.25681371
	   ni20(37,1) = 0.25214660
	   ni20(38,1) = 0.24700487
	   ni20(39,1) = 0.24122889
	   ni20(40,1) = 0.23466847
	   ni20(41,1) = 0.22718073
	   ni20(42,1) = 0.21857329
	   ni20(43,1) = 0.20870564
	   ni20(44,1) = 0.19738348
	   ni20(45,1) = 0.18439192
	   ni20(46,1) = 0.16952788
	   ni20(47,1) = 0.15249520
	   ni20(48,1) = 0.13245345
	   ni20(49,1) = 0.10504354
	   ni20(50,1) = 0.060652772
	   ni20(51,1) = 0.031353233
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
