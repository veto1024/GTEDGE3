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
	   tekev(1) = 1.4893080
	   tekev(2) = 1.4856134
	   tekev(3) = 1.4807199
	   tekev(4) = 1.4734286
	   tekev(5) = 1.4626239
	   tekev(6) = 1.4482215
	   tekev(7) = 1.4312639
	   tekev(8) = 1.4118733
	   tekev(9) = 1.3904122
	   tekev(10) = 1.3668351
	   tekev(11) = 1.3413980
	   tekev(12) = 1.3143476
	   tekev(13) = 1.2857830
	   tekev(14) = 1.2558559
	   tekev(15) = 1.2247263
	   tekev(16) = 1.1925417
	   tekev(17) = 1.1594423
	   tekev(18) = 1.1255767
	   tekev(19) = 1.0910572
	   tekev(20) = 1.0560019
	   tekev(21) = 1.0205316
	   tekev(22) = 0.98475845
	   tekev(23) = 0.94877725
	   tekev(24) = 0.91269012
	   tekev(25) = 0.87658320
	   tekev(26) = 0.84054008
	   tekev(27) = 0.80463688
	   tekev(28) = 0.76894460
	   tekev(29) = 0.73352623
	   tekev(30) = 0.69843766
	   tekev(31) = 0.66372716
	   tekev(32) = 0.62943633
	   tekev(33) = 0.59559939
	   tekev(34) = 0.56224179
	   tekev(35) = 0.52937895
	   tekev(36) = 0.49701165
	   tekev(37) = 0.46513772
	   tekev(38) = 0.43373222
	   tekev(39) = 0.40275835
	   tekev(40) = 0.37216178
	   tekev(41) = 0.34186732
	   tekev(42) = 0.31178172
	   tekev(43) = 0.28179957
	   tekev(44) = 0.25181613
	   tekev(45) = 0.22175970
	   tekev(46) = 0.19164594
	   tekev(47) = 0.16166320
	   tekev(48) = 0.13227787
	   tekev(49) = 0.10430768
	   tekev(50) = 0.078888166
	   tekev(51) = 0.057228175
	   tikev(1) = 1.8392977
	   tikev(2) = 1.8345398
	   tikev(3) = 1.8279815
	   tikev(4) = 1.8178222
	   tikev(5) = 1.8025279
	   tikev(6) = 1.7822470
	   tikev(7) = 1.7586707
	   tikev(8) = 1.7313187
	   tikev(9) = 1.7007193
	   tikev(10) = 1.6671989
	   tikev(11) = 1.6308335
	   tikev(12) = 1.5918911
	   tikev(13) = 1.5506247
	   tikev(14) = 1.5072704
	   tikev(15) = 1.4621043
	   tikev(16) = 1.4153378
	   tikev(17) = 1.3671953
	   tikev(18) = 1.3179503
	   tikev(19) = 1.2678210
	   tikev(20) = 1.2170639
	   tikev(21) = 1.1659090
	   tikev(22) = 1.1146029
	   tikev(23) = 1.0633847
	   tikev(24) = 1.0124956
	   tikev(25) = 0.96217721
	   tikev(26) = 0.91267182
	   tikev(27) = 0.86422019
	   tikev(28) = 0.81705420
	   tikev(29) = 0.77142453
	   tikev(30) = 0.72757507
	   tikev(31) = 0.68574633
	   tikev(32) = 0.64616117
	   tikev(33) = 0.60908316
	   tikev(34) = 0.57475278
	   tikev(35) = 0.54339176
	   tikev(36) = 0.51524813
	   tikev(37) = 0.49020244
	   tikev(38) = 0.46750093
	   tikev(39) = 0.44633765
	   tikev(40) = 0.42591570
	   tikev(41) = 0.40542592
	   tikev(42) = 0.38405862
	   tikev(43) = 0.36148012
	   tikev(44) = 0.33839511
	   tikev(45) = 0.31564658
	   tikev(46) = 0.29408589
	   tikev(47) = 0.27453709
	   tikev(48) = 0.25787740
	   tikev(49) = 0.24492295
	   tikev(50) = 0.23649867
	   tikev(51) = 0.23355717
	   ne20(1) = 0.47005719
	   ne20(2) = 0.46876424
	   ne20(3) = 0.46706418
	   ne20(4) = 0.46454991
	   ne20(5) = 0.46086949
	   ne20(6) = 0.45611315
	   ne20(7) = 0.45071474
	   ne20(8) = 0.44468209
	   ne20(9) = 0.43820966
	   ne20(10) = 0.43130922
	   ne20(11) = 0.42410866
	   ne20(12) = 0.41669268
	   ne20(13) = 0.40912262
	   ne20(14) = 0.40147428
	   ne20(15) = 0.39381335
	   ne20(16) = 0.38619464
	   ne20(17) = 0.37867023
	   ne20(18) = 0.37128615
	   ne20(19) = 0.36408113
	   ne20(20) = 0.35707257
	   ne20(21) = 0.35030007
	   ne20(22) = 0.34378409
	   ne20(23) = 0.33753267
	   ne20(24) = 0.33153543
	   ne20(25) = 0.32581523
	   ne20(26) = 0.32036312
	   ne20(27) = 0.31513762
	   ne20(28) = 0.31015111
	   ne20(29) = 0.30536227
	   ne20(30) = 0.30073347
	   ne20(31) = 0.29623898
	   ne20(32) = 0.29181663
	   ne20(33) = 0.28742176
	   ne20(34) = 0.28299502
	   ne20(35) = 0.27847294
	   ne20(36) = 0.27377243
	   ne20(37) = 0.26881676
	   ne20(38) = 0.26354263
	   ne20(39) = 0.25783646
	   ne20(40) = 0.25160205
	   ne20(41) = 0.24473419
	   ne20(42) = 0.23713512
	   ne20(43) = 0.22867685
	   ne20(44) = 0.21923700
	   ne20(45) = 0.20868287
	   ne20(46) = 0.19687750
	   ne20(47) = 0.18367659
	   ne20(48) = 0.16889993
	   ne20(49) = 0.15230167
	   ne20(50) = 0.13334427
	   ne20(51) = 0.11060181
	   ni20(1,1) = 0.43763440
	   ni20(2,1) = 0.43723506
	   ni20(3,1) = 0.43632669
	   ni20(4,1) = 0.43440024
	   ni20(5,1) = 0.43103354
	   ni20(6,1) = 0.42650055
	   ni20(7,1) = 0.42150641
	   ni20(8,1) = 0.41610206
	   ni20(9,1) = 0.41048050
	   ni20(10,1) = 0.40461958
	   ni20(11,1) = 0.39855980
	   ni20(12,1) = 0.39222604
	   ni20(13,1) = 0.38560924
	   ni20(14,1) = 0.37877866
	   ni20(15,1) = 0.37186571
	   ni20(16,1) = 0.36494060
	   ni20(17,1) = 0.35808465
	   ni20(18,1) = 0.35136725
	   ni20(19,1) = 0.34482197
	   ni20(20,1) = 0.33840066
	   ni20(21,1) = 0.33211803
	   ni20(22,1) = 0.32603368
	   ni20(23,1) = 0.32020420
	   ni20(24,1) = 0.31461199
	   ni20(25,1) = 0.30923437
	   ni20(26,1) = 0.30405153
	   ni20(27,1) = 0.29905865
	   ni20(28,1) = 0.29428838
	   ni20(29,1) = 0.28971867
	   ni20(30,1) = 0.28531158
	   ni20(31,1) = 0.28104152
	   ni20(32,1) = 0.27684733
	   ni20(33,1) = 0.27267271
	   ni20(34,1) = 0.26845273
	   ni20(35,1) = 0.26412170
	   ni20(36,1) = 0.25960366
	   ni20(37,1) = 0.25483834
	   ni20(38,1) = 0.24979418
	   ni20(39,1) = 0.24436133
	   ni20(40,1) = 0.23843544
	   ni20(41,1) = 0.23191306
	   ni20(42,1) = 0.22469587
	   ni20(43,1) = 0.21665287
	   ni20(44,1) = 0.20766981
	   ni20(45,1) = 0.19761755
	   ni20(46,1) = 0.18637724
	   ni20(47,1) = 0.17390483
	   ni20(48,1) = 0.16013349
	   ni20(49,1) = 0.14480887
	   ni20(50,1) = 0.12705557
	   ni20(51,1) = 0.10541320
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
