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
	   tekev(1) = 1.8832686
	   tekev(2) = 1.8781308
	   tekev(3) = 1.8711249
	   tekev(4) = 1.8603828
	   tekev(5) = 1.8443631
	   tekev(6) = 1.8234325
	   tekev(7) = 1.7994619
	   tekev(8) = 1.7720188
	   tekev(9) = 1.7417413
	   tekev(10) = 1.7089868
	   tekev(11) = 1.6739145
	   tekev(12) = 1.6368429
	   tekev(13) = 1.5980545
	   tekev(14) = 1.5578087
	   tekev(15) = 1.5163506
	   tekev(16) = 1.4739126
	   tekev(17) = 1.4307132
	   tekev(18) = 1.3869693
	   tekev(19) = 1.3428603
	   tekev(20) = 1.2985664
	   tekev(21) = 1.2542571
	   tekev(22) = 1.2100817
	   tekev(23) = 1.1661840
	   tekev(24) = 1.1226815
	   tekev(25) = 1.0796968
	   tekev(26) = 1.0373196
	   tekev(27) = 0.99563962
	   tekev(28) = 0.95473551
	   tekev(29) = 0.91464968
	   tekev(30) = 0.87543992
	   tekev(31) = 0.83713986
	   tekev(32) = 0.79976544
	   tekev(33) = 0.76331487
	   tekev(34) = 0.72778646
	   tekev(35) = 0.69315887
	   tekev(36) = 0.65939791
	   tekev(37) = 0.62645392
	   tekev(38) = 0.59426653
	   tekev(39) = 0.56275986
	   tekev(40) = 0.53184556
	   tekev(41) = 0.50141949
	   tekev(42) = 0.47136988
	   tekev(43) = 0.44156152
	   tekev(44) = 0.41183322
	   tekev(45) = 0.38194337
	   tekev(46) = 0.35135339
	   tekev(47) = 0.31828974
	   tekev(48) = 0.27635894
	   tekev(49) = 0.20758859
	   tekev(50) = 0.10540371
	   tekev(51) = 0.034165660
	   tikev(1) = 2.0151334
	   tikev(2) = 2.0098278
	   tikev(3) = 2.0025204
	   tikev(4) = 1.9912093
	   tikev(5) = 1.9742264
	   tikev(6) = 1.9518517
	   tikev(7) = 1.9260094
	   tikev(8) = 1.8961970
	   tikev(9) = 1.8630252
	   tikev(10) = 1.8269261
	   tikev(11) = 1.7880279
	   tikev(12) = 1.7466806
	   tikev(13) = 1.7032090
	   tikev(14) = 1.6579320
	   tikev(15) = 1.6111641
	   tikev(16) = 1.5632380
	   tikev(17) = 1.5144553
	   tikev(18) = 1.4651231
	   tikev(19) = 1.4155681
	   tikev(20) = 1.3660970
	   tikev(21) = 1.3169075
	   tikev(22) = 1.2680582
	   tikev(23) = 1.2196272
	   tikev(24) = 1.1716612
	   tikev(25) = 1.1242399
	   tikev(26) = 1.0774224
	   tikev(27) = 1.0312660
	   tikev(28) = 0.98584382
	   tikev(29) = 0.94121875
	   tikev(30) = 0.89745395
	   tikev(31) = 0.85461239
	   tikev(32) = 0.81275647
	   tikev(33) = 0.77195237
	   tikev(34) = 0.73226165
	   tikev(35) = 0.69375158
	   tikev(36) = 0.65648622
	   tikev(37) = 0.62053057
	   tikev(38) = 0.58595054
	   tikev(39) = 0.55280725
	   tikev(40) = 0.52115705
	   tikev(41) = 0.49107868
	   tikev(42) = 0.46263028
	   tikev(43) = 0.43586838
	   tikev(44) = 0.41087029
	   tikev(45) = 0.38769874
	   tikev(46) = 0.36636561
	   tikev(47) = 0.34602569
	   tikev(48) = 0.32487134
	   tikev(49) = 0.30109149
	   tikev(50) = 0.27284686
	   tikev(51) = 0.23831833
	   ne20(1) = 0.39462965
	   ne20(2) = 0.39374243
	   ne20(3) = 0.39266801
	   ne20(4) = 0.39121921
	   ne20(5) = 0.38920886
	   ne20(6) = 0.38651952
	   ne20(7) = 0.38325676
	   ne20(8) = 0.37963890
	   ne20(9) = 0.37571938
	   ne20(10) = 0.37150145
	   ne20(11) = 0.36708983
	   ne20(12) = 0.36248409
	   ne20(13) = 0.35775797
	   ne20(14) = 0.35293736
	   ne20(15) = 0.34806422
	   ne20(16) = 0.34317690
	   ne20(17) = 0.33830555
	   ne20(18) = 0.33347733
	   ne20(19) = 0.32872933
	   ne20(20) = 0.32406833
	   ne20(21) = 0.31953537
	   ne20(22) = 0.31513262
	   ne20(23) = 0.31088674
	   ne20(24) = 0.30680482
	   ne20(25) = 0.30289923
	   ne20(26) = 0.29917411
	   ne20(27) = 0.29564306
	   ne20(28) = 0.29228556
	   ne20(29) = 0.28914862
	   ne20(30) = 0.28614566
	   ne20(31) = 0.28334905
	   ne20(32) = 0.28072142
	   ne20(33) = 0.27822363
	   ne20(34) = 0.27588514
	   ne20(35) = 0.27366469
	   ne20(36) = 0.27152423
	   ne20(37) = 0.26947425
	   ne20(38) = 0.26746242
	   ne20(39) = 0.26546177
	   ne20(40) = 0.26345758
	   ne20(41) = 0.26138170
	   ne20(42) = 0.25922481
	   ne20(43) = 0.25693960
	   ne20(44) = 0.25445078
	   ne20(45) = 0.25177409
	   ne20(46) = 0.24878922
	   ne20(47) = 0.24552982
	   ne20(48) = 0.24196718
	   ne20(49) = 0.23706265
	   ne20(50) = 0.21659900
	   ne20(51) = 0.099052042
	   ni20(1,1) = 0.34668707
	   ni20(2,1) = 0.34571723
	   ni20(3,1) = 0.34463962
	   ni20(4,1) = 0.34334647
	   ni20(5,1) = 0.34173003
	   ni20(6,1) = 0.33974238
	   ni20(7,1) = 0.33752689
	   ni20(8,1) = 0.33536026
	   ni20(9,1) = 0.33324640
	   ni20(10,1) = 0.33069424
	   ni20(11,1) = 0.32754876
	   ni20(12,1) = 0.32419951
	   ni20(13,1) = 0.32086820
	   ni20(14,1) = 0.31752450
	   ni20(15,1) = 0.31418578
	   ni20(16,1) = 0.31078693
	   ni20(17,1) = 0.30731775
	   ni20(18,1) = 0.30379375
	   ni20(19,1) = 0.30021074
	   ni20(20,1) = 0.29658561
	   ni20(21,1) = 0.29291872
	   ni20(22,1) = 0.28922550
	   ni20(23,1) = 0.28554451
	   ni20(24,1) = 0.28193387
	   ni20(25,1) = 0.27861451
	   ni20(26,1) = 0.27560501
	   ni20(27,1) = 0.27272814
	   ni20(28,1) = 0.26991166
	   ni20(29,1) = 0.26724762
	   ni20(30,1) = 0.26467520
	   ni20(31,1) = 0.26227363
	   ni20(32,1) = 0.26004790
	   ni20(33,1) = 0.25797496
	   ni20(34,1) = 0.25607889
	   ni20(35,1) = 0.25432334
	   ni20(36,1) = 0.25266099
	   ni20(37,1) = 0.25111883
	   ni20(38,1) = 0.24965217
	   ni20(39,1) = 0.24822202
	   ni20(40,1) = 0.24682752
	   ni20(41,1) = 0.24541096
	   ni20(42,1) = 0.24396796
	   ni20(43,1) = 0.24244076
	   ni20(44,1) = 0.24071638
	   ni20(45,1) = 0.23879017
	   ni20(46,1) = 0.23651660
	   ni20(47,1) = 0.23385014
	   ni20(48,1) = 0.23157334
	   ni20(49,1) = 0.22890768
	   ni20(50,1) = 0.21175512
	   ni20(51,1) = 0.098029233
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
	   ni20(i,j) = (ne20(i) - ni20(i,1))/zion(j)
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
