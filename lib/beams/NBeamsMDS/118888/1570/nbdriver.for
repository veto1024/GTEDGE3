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
	   tekev(1) = 2.0269216
	   tekev(2) = 2.0218314
	   tekev(3) = 2.0151545
	   tekev(4) = 2.0053041
	   tekev(5) = 1.9907954
	   tekev(6) = 1.9715262
	   tekev(7) = 1.9489061
	   tekev(8) = 1.9232749
	   tekev(9) = 1.8949943
	   tekev(10) = 1.8641062
	   tekev(11) = 1.8310656
	   tekev(12) = 1.7960952
	   tekev(13) = 1.7593622
	   tekev(14) = 1.7211111
	   tekev(15) = 1.6815624
	   tekev(16) = 1.6409155
	   tekev(17) = 1.5993561
	   tekev(18) = 1.5570548
	   tekev(19) = 1.5141720
	   tekev(20) = 1.4708550 
	   tekev(21) = 1.4272393
	   tekev(22) = 1.3834502
	   tekev(23) = 1.3395896
	   tekev(24) = 1.2957568
	   tekev(25) = 1.2520340
	   tekev(26) = 1.2084914
	   tekev(27) = 1.1651899
	   tekev(28) = 1.1221668
	   tekev(29) = 1.0794628
	   tekev(30) = 1.0370893
	   tekev(31) = 0.99505654
	   tekev(32) = 0.95335915
	   tekev(33) = 0.91197301
	   tekev(34) = 0.87086815
	   tekev(35) = 0.83000020
	   tekev(36) = 0.78931042
	   tekev(37) = 0.74872768
	   tekev(38) = 0.70816835
	   tekev(39) = 0.66753555
	   tekev(40) = 0.62672107
	   tekev(41) = 0.58560057
	   tekev(42) = 0.54404016
	   tekev(43) = 0.50189328
	   tekev(44) = 0.45899709
	   tekev(45) = 0.41517897
	   tekev(46) = 0.37025013
	   tekev(47) = 0.32399898
	   tekev(48) = 0.27546950
	   tekev(49) = 0.197195916
	   tekev(50) = 0.0479290792
	   tekev(51) = 0.0245520347
	   tikev(1) = 2.2341931
	   tikev(2) = 2.2296422
	   tikev(3) = 2.2237494
	   tikev(4) = 2.2151726
	   tikev(5) = 2.2025741
	   tikev(6) = 2.1851483
	   tikev(7) = 2.1631229
	   tikev(8) = 2.1380800
	   tikev(9) = 2.1088528
	   tikev(10) = 2.0761979
	   tikev(11) = 2.0399002
	   tikev(12) = 1.9999662
	   tikev(13) = 1.9565054
	   tikev(14) = 1.9095745
	   tikev(15) = 1.8592156
	   tikev(16) = 1.8054531
	   tikev(17) = 1.7483177
	   tikev(18) = 1.6879326
	   tikev(19) = 1.6242562
	   tikev(20) = 1.5573986 
	   tikev(21) = 1.4874884
	   tikev(22) = 1.4155703
	   tikev(23) = 1.3430421
	   tikev(24) = 1.2713004
	   tikev(25) = 1.2015457
	   tikev(26) = 1.1344570
	   tikev(27) = 1.0705873
	   tikev(28) = 1.0105488
	   tikev(29) = 0.95492739
	   tikev(30) = 0.90430455
	   tikev(31) = 0.85926050
	   tikev(32) = 0.82026079
	   tikev(33) = 0.78684088
	   tikev(34) = 0.75782151
	   tikev(35) = 0.73233322
	   tikev(36) = 0.70929328
	   tikev(37) = 0.68763530
	   tikev(38) = 0.66637921
	   tikev(39) = 0.64447144
	   tikev(40) = 0.62090456
	   tikev(41) = 0.59512321
	   tikev(42) = 0.56717039
	   tikev(43) = 0.53708800
	   tikev(44) = 0.50490813
	   tikev(45) = 0.47067933
	   tikev(46) = 0.43444728
	   tikev(47) = 0.39625654
	   tikev(48) = 0.35615780
	   tikev(49) = 0.31419037
	   tikev(50) = 0.27040326
	   tikev(51) = 0.23845654
	   ne20(1) = 0.44743947
	   ne20(2) = 0.44638813
	   ne20(3) = 0.44509904
	   ne20(4) = 0.44333447
	   ne20(5) = 0.44085749
	   ne20(6) = 0.43756060
	   ne20(7) = 0.43360558
	   ne20(8) = 0.42929454
	   ne20(9) = 0.42458385
	   ne20(10) = 0.41959089
	   ne20(11) = 0.41435227
	   ne20(12) = 0.40895092
	   ne20(13) = 0.40342120
	   ne20(14) = 0.39783731
	   ne20(15) = 0.39223543
	   ne20(16) = 0.38666270
	   ne20(17) = 0.38117014
	   ne20(18) = 0.37578279
	   ne20(19) = 0.37054838
	   ne20(20) = 0.36549187
	   ne20(21) = 0.36063769
	   ne20(22) = 0.35602046
	   ne20(23) = 0.35164563
	   ne20(24) = 0.34753713
	   ne20(25) = 0.34372948
	   ne20(26) = 0.34015902
	   ne20(27) = 0.33693456
	   ne20(28) = 0.33399916
	   ne20(29) = 0.33130249
	   ne20(30) = 0.32892843
	   ne20(31) = 0.32685744
	   ne20(32) = 0.32507003
	   ne20(33) = 0.32350335
	   ne20(34) = 0.32211946
	   ne20(35) = 0.32098853
	   ne20(36) = 0.32004860
	   ne20(37) = 0.31923773
	   ne20(38) = 0.31849389
	   ne20(39) = 0.31784067
	   ne20(40) = 0.31725160
	   ne20(41) = 0.31663441
	   ne20(42) = 0.31589683
	   ne20(43) = 0.31499049
	   ne20(44) = 0.31428912
	   ne20(45) = 0.31317645
	   ne20(46) = 0.31124185
	   ne20(47) = 0.30897490
	   ne20(48) = 0.30271805
	   ne20(49) = 0.27247675
	   ne20(50) = 0.15097012
	   ne20(51) = 0.031064493
	   ni20(1,1) = 0.39764799
	   ni20(2,1) = 0.39536762
	   ni20(3,1) = 0.39290777
	   ni20(4,1) = 0.39008897
	   ni20(5,1) = 0.38673301
	   ni20(6,1) = 0.38286224
	   ni20(7,1) = 0.37891629
	   ni20(8,1) = 0.37538713
	   ni20(9,1) = 0.37218215
	   ni20(10,1) = 0.36909194
	   ni20(11,1) = 0.36609075
	   ni20(12,1) = 0.36316815
	   ni20(13,1) = 0.36026544
	   ni20(14,1) = 0.35730383
	   ni20(15,1) = 0.35414845
	   ni20(16,1) = 0.35070280
	   ni20(17,1) = 0.34694175
	   ni20(18,1) = 0.34291529
	   ni20(19,1) = 0.33890290
	   ni20(20,1) = 0.33503047
	   ni20(21,1) = 0.33141873
	   ni20(22,1) = 0.32804501
	   ni20(23,1) = 0.32470012
	   ni20(24,1) = 0.32138286
	   ni20(25,1) = 0.31817812
	   ni20(26,1) = 0.31506742
	   ni20(27,1) = 0.31214918
	   ni20(28,1) = 0.30944787
	   ni20(29,1) = 0.30694420
	   ni20(30,1) = 0.30468146
	   ni20(31,1) = 0.30265089
	   ni20(32,1) = 0.30084378
	   ni20(33,1) = 0.29922899
	   ni20(34,1) = 0.29779774
	   ni20(35,1) = 0.29656544
	   ni20(36,1) = 0.29546582
	   ni20(37,1) = 0.29443263
	   ni20(38,1) = 0.29336739
	   ni20(39,1) = 0.29215876
	   ni20(40,1) = 0.29096893
	   ni20(41,1) = 0.28982029
	   ni20(42,1) = 0.28873525
	   ni20(43,1) = 0.28775373
	   ni20(44,1) = 0.28704302
	   ni20(45,1) = 0.28631868
	   ni20(46,1) = 0.28538427
	   ni20(47,1) = 0.28477490
	   ni20(48,1) = 0.28117824
	   ni20(49,1) = 0.25581078
	   ni20(50,1) = 0.14372019
	   ni20(51,1) = 0.030039120
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
