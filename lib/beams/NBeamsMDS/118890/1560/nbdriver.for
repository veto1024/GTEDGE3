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
	   tekev(1) = 1.8212879
	   tekev(2) = 1.8167179
	   tekev(3) = 1.8107274
	   tekev(4) = 1.8018957
	   tekev(5) = 1.7888941
	   tekev(6) = 1.7716376
	   tekev(7) = 1.7513962
	   tekev(8) = 1.7284783
	   tekev(9) = 1.7032122
	   tekev(10) = 1.6756415
	   tekev(11) = 1.6461759
	   tekev(12) = 1.6150212
	   tekev(13) = 1.5823252
	   tekev(14) = 1.5483159
	   tekev(15) = 1.5131914
	   tekev(16) = 1.4771343
	   tekev(17) = 1.4403138
	   tekev(18) = 1.4028860
	   tekev(19) = 1.3649952
	   tekev(20) = 1.3267752
	   tekev(21) = 1.2883486
	   tekev(22) = 1.2498272
	   tekev(23) = 1.2113085
	   tekev(24) = 1.1728782
	   tekev(25) = 1.1346144
	   tekev(26) = 1.0965778
	   tekev(27) = 1.0588191
	   tekev(28) = 1.0213855
	   tekev(29) = 0.98429195
	   tekev(30) = 0.94757061
	   tekev(31) = 0.91121052
	   tekev(32) = 0.87521755
	   tekev(33) = 0.83956588
	   tekev(34) = 0.80422623
	   tekev(35) = 0.76915834
	   tekev(36) = 0.73430654
	   tekev(37) = 0.69960413
	   tekev(38) = 0.66497461
	   tekev(39) = 0.63032857
	   tekev(40) = 0.59556566
	   tekev(41) = 0.56057195
	   tekev(42) = 0.52522292
	   tekev(43) = 0.48938004
	   tekev(44) = 0.45289204
	   tekev(45) = 0.41556780
	   tekev(46) = 0.37705472
	   tekev(47) = 0.33624613
	   tekev(48) = 0.28829169
	   tekev(49) = 0.21532502
	   tekev(50) = 0.10512895
	   tekev(51) = 0.036034448
	   tikev(1) = 2.1994154
	   tikev(2) = 2.1935968
	   tikev(3) = 2.1854520
	   tikev(4) = 2.1726548
	   tikev(5) = 2.1533149
	   tikev(6) = 2.1277266
	   tikev(7) = 2.0980797
	   tikev(8) = 2.0636053
	   tikev(9) = 2.0248720
	   tikev(10) = 1.9822216
	   tikev(11) = 1.9359294
	   tikev(12) = 1.8862574
	   tikev(13) = 1.8334689
	   tikev(14) = 1.7778242
	   tikev(15) = 1.7195939
	   tikev(16) = 1.6591052
	   tikev(17) = 1.5966074
	   tikev(18) = 1.5323978
	   tikev(19) = 1.4667368
	   tikev(20) = 1.3999256
	   tikev(21) = 1.3325083
	   tikev(22) = 1.2653819
	   tikev(23) = 1.1994659
	   tikev(24) = 1.1356878
	   tikev(25) = 1.0749744
	   tikev(26) = 1.0182519
	   tikev(27) = 0.96644733
	   tikev(28) = 0.92009411
	   tikev(29) = 0.87873533
	   tikev(30) = 0.84168700
	   tikev(31) = 0.80837800
	   tikev(32) = 0.77816557
	   tikev(33) = 0.75042922
	   tikev(34) = 0.72451820
	   tikev(35) = 0.69983928
	   tikev(36) = 0.67573642
	   tikev(37) = 0.65159953
	   tikev(38) = 0.62679356
	   tikev(39) = 0.60068338
	   tikev(40) = 0.57279480
	   tikev(41) = 0.54323050
	   tikev(42) = 0.51229388
	   tikev(43) = 0.48027357
	   tikev(44) = 0.44747759
	   tikev(45) = 0.41419875
	   tikev(46) = 0.38073655
	   tikev(47) = 0.34738802
	   tikev(48) = 0.31445218
	   tikev(49) = 0.28222642
	   tikev(50) = 0.25100995
	   tikev(51) = 0.22109747
	   ne20(1) = 0.45848258
	   ne20(2) = 0.45728153
	   ne20(3) = 0.45575262
	   ne20(4) = 0.45356799
	   ne20(5) = 0.45041868
	   ne20(6) = 0.44630857
	   ne20(7) = 0.44155991
	   ne20(8) = 0.43637508
	   ne20(9) = 0.43075959
	   ne20(10) = 0.42484090
	   ne20(11) = 0.41869578
	   ne20(12) = 0.41238569
	   ne20(13) = 0.40600630
	   ne20(14) = 0.39960702
	   ne20(15) = 0.39325205
	   ne20(16) = 0.38699155
	   ne20(17) = 0.38088759
	   ne20(18) = 0.37496375
	   ne20(19) = 0.36926961
	   ne20(20) = 0.36384657
	   ne20(21) = 0.35868228
	   ne20(22) = 0.35386080
	   ne20(23) = 0.34930592
	   ne20(24) = 0.34512653
	   ne20(25) = 0.34123658
	   ne20(26) = 0.33768800
	   ne20(27) = 0.33450735
	   ne20(28) = 0.33153983
	   ne20(29) = 0.32890213
	   ne20(30) = 0.32653883
	   ne20(31) = 0.32438654
	   ne20(32) = 0.32238663
	   ne20(33) = 0.32057802
	   ne20(34) = 0.31886609
	   ne20(35) = 0.31715622
	   ne20(36) = 0.31539866
	   ne20(37) = 0.31360831
	   ne20(38) = 0.31164163
	   ne20(39) = 0.30937909
	   ne20(40) = 0.30686877
	   ne20(41) = 0.30387472
	   ne20(42) = 0.30042312
	   ne20(43) = 0.29630240
	   ne20(44) = 0.29150913
	   ne20(45) = 0.28587409
	   ne20(46) = 0.27928009
	   ne20(47) = 0.27160592
	   ne20(48) = 0.26273440
	   ne20(49) = 0.25251362
	   ne20(50) = 0.23368536
	   ne20(51) = 0.079020123
	   ni20(1,1) = 0.40424907
	   ni20(2,1) = 0.40302681
	   ni20(3,1) = 0.40163307
	   ni20(4,1) = 0.39989635
	   ni20(5,1) = 0.39765541
	   ni20(6,1) = 0.39491818
	   ni20(7,1) = 0.39188093
	   ni20(8,1) = 0.38866001
	   ni20(9,1) = 0.38504994
	   ni20(10,1) = 0.38105181
	   ni20(11,1) = 0.37670702
	   ni20(12,1) = 0.37213346
	   ni20(13,1) = 0.36763592
	   ni20(14,1) = 0.36296889
	   ni20(15,1) = 0.35805808
	   ni20(16,1) = 0.35286908
	   ni20(17,1) = 0.34781887
	   ni20(18,1) = 0.34326227
	   ni20(19,1) = 0.33900369
	   ni20(20,1) = 0.33496327
	   ni20(21,1) = 0.33096369
	   ni20(22,1) = 0.32706898
	   ni20(23,1) = 0.32324952
	   ni20(24,1) = 0.31962362
	   ni20(25,1) = 0.31621140
	   ni20(26,1) = 0.31306613
	   ni20(27,1) = 0.31025915
	   ni20(28,1) = 0.30772990
	   ni20(29,1) = 0.30553082
	   ni20(30,1) = 0.30358573
	   ni20(31,1) = 0.30180725
	   ni20(32,1) = 0.30009455
	   ni20(33,1) = 0.29853347
	   ni20(34,1) = 0.29705993
	   ni20(35,1) = 0.29560985
	   ni20(36,1) = 0.29416227
	   ni20(37,1) = 0.29272582
	   ni20(38,1) = 0.29115694
	   ni20(39,1) = 0.28933866
	   ni20(40,1) = 0.28730570
	   ni20(41,1) = 0.28485169
	   ni20(42,1) = 0.28197334
	   ni20(43,1) = 0.27854354
	   ni20(44,1) = 0.27460367
	   ni20(45,1) = 0.27005682
	   ni20(46,1) = 0.26484213
	   ni20(47,1) = 0.25889852
	   ni20(48,1) = 0.25214460
	   ni20(49,1) = 0.24443383
	   ni20(50,1) = 0.22863371
	   ni20(51,1) = 0.077950762
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
