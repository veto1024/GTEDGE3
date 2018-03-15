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
     tekev(1)  =1.8768246591
     tekev(2)  =1.8547741588
     tekev(3)  =1.83190556824
     tekev(4)  =1.80823851509
     tekev(5)  =1.783792627
     tekev(6)  =1.75858753165
     tekev(7)  =1.73264285671
     tekev(8)  =1.70597822984
     tekev(9)  =1.67861327872
     tekev(10)  =1.65056763101
     tekev(11)  =1.62186091437
     tekev(12)  =1.59251275648
     tekev(13)  =1.562542785
     tekev(14)  =1.53197062761
     tekev(15)  =1.50081591197
     tekev(16)  =1.46909826574
     tekev(17)  =1.4368373166
     tekev(18)  =1.40405269222
     tekev(19)  =1.37076402026
     tekev(20)  =1.33699092838
     tekev(21)  =1.30275304427
     tekev(22)  =1.26806999558
     tekev(23)  =1.23296140999
     tekev(24)  =1.19744691515
     tekev(25)  =1.16154613875
     tekev(26)  =1.12527870845
     tekev(27)  =1.08866425191
     tekev(28)  =1.0517223968
     tekev(29)  =1.0144727708
     tekev(30)  =0.976935001561
     tekev(31)  =0.939128716765
     tekev(32)  =0.901073544073
     tekev(33)  =0.862789111155
     tekev(34)  =0.824295045679
     tekev(35)  =0.785610975312
     tekev(36)  =0.746756527722
     tekev(37)  =0.707751330578
     tekev(38)  =0.668615011547
     tekev(39)  =0.629367198298
     tekev(40)  =0.590027518497
     tekev(41)  =0.550615599814
     tekev(42)  =0.511151069916
     tekev(43)  =0.471653556471
     tekev(44)  =0.432142687148
     tekev(45)  =0.392638089613
     tekev(46)  =0.353159391535
     tekev(47)  =0.313726220583
     tekev(48)  =0.274358204423
     tekev(49)  =0.235074970724
     tekev(50)  =0.195896147154
     tekev(51)  =0.156841361381
     tikev(1)  =2.30373318812
     tikev(2)  =2.26838660192
     tikev(3)  =2.23113145578
     tikev(4)  =2.19207627903
     tikev(5)  =2.15132960101
     tikev(6)  =2.10899995107
     tikev(7)  =2.06519585855
     tikev(8)  =2.02002585279
     tikev(9)  =1.97359846313
     tikev(10)  =1.92602221891
     tikev(11)  =1.87740564948
     tikev(12)  =1.82785728416
     tikev(13)  =1.77748565231
     tikev(14)  =1.72639928327
     tikev(15)  =1.67470670637
     tikev(16)  =1.62251645097
     tikev(17)  =1.56993704638
     tikev(18)  =1.51707702197
     tikev(19)  =1.46404490708
     tikev(20)  =1.41094923103
     tikev(21)  =1.35789852317
     tikev(22)  =1.30500131285
     tikev(23)  =1.25236612941
     tikev(24)  =1.20010150218
     tikev(25)  =1.14831596051
     tikev(26)  =1.09711803374
     tikev(27)  =1.04661625121
     tikev(28)  =0.996919142255
     tikev(29)  =0.948135236226
     tikev(30)  =0.900373062459
     tikev(31)  =0.853741150296
     tikev(32)  =0.808348029078
     tikev(33)  =0.764302228147
     tikev(34)  =0.721712276843
     tikev(35)  =0.680686704507
     tikev(36)  =0.641334040481
     tikev(37)  =0.603762814105
     tikev(38)  =0.568081554721
     tikev(39)  =0.534398791669
     tikev(40)  =0.502823054291
     tikev(41)  =0.473462871929
     tikev(42)  =0.446426773922
     tikev(43)  =0.421823289612
     tikev(44)  =0.39976094834
     tikev(45)  =0.380348279447
     tikev(46)  =0.363693812275
     tikev(47)  =0.349906076164
     tikev(48)  =0.339093600456
     tikev(49)  =0.331364914491
     tikev(50)  =0.326828547611
     tikev(51)  =0.325593029156
     ne20(1)  =0.46317615728
     ne20(2)  =0.454782791851
     ne20(3)  =0.446776477149
     ne20(4)  =0.439142531331
     ne20(5)  =0.431866272554
     ne20(6)  =0.424933018975
     ne20(7)  =0.418328088751
     ne20(8)  =0.412036800037
     ne20(9)  =0.406044470993
     ne20(10)  =0.400336419773
     ne20(11)  =0.394897964535
     ne20(12)  =0.389714423437
     ne20(13)  =0.384771114634
     ne20(14)  =0.380053356283
     ne20(15)  =0.375546466543
     ne20(16)  =0.371235763568
     ne20(17)  =0.367106565516
     ne20(18)  =0.363144190545
     ne20(19)  =0.35933395681
     ne20(20)  =0.355661182469
     ne20(21)  =0.352111185679
     ne20(22)  =0.348669284596
     ne20(23)  =0.345320797378
     ne20(24)  =0.34205104218
     ne20(25)  =0.338845337161
     ne20(26)  =0.335689000476
     ne20(27)  =0.332567350283
     ne20(28)  =0.329465704739
     ne20(29)  =0.326369382
     ne20(30)  =0.323263700223
     ne20(31)  =0.320133977566
     ne20(32)  =0.316965532184
     ne20(33)  =0.313743682236
     ne20(34)  =0.310453745877
     ne20(35)  =0.307081041264
     ne20(36)  =0.303610886555
     ne20(37)  =0.300028599906
     ne20(38)  =0.296319499475
     ne20(39)  =0.292468903417
     ne20(40)  =0.28846212989
     ne20(41)  =0.284284497051
     ne20(42)  =0.279921323056
     ne20(43)  =0.275357926063
     ne20(44)  =0.270579624227
     ne20(45)  =0.265571735707
     ne20(46)  =0.260319578659
     ne20(47)  =0.25480847124
     ne20(48)  =0.249023731607
     ne20(49)  =0.242950677916
     ne20(50)  =0.236574628324
     ne20(51)  =0.229880900988
     ni20(1,1)  =0.403437636806
     ni20(2,1)  =0.397400523606
     ni20(3,1)  =0.391631923216
     ni20(4,1)  =0.386121078024
     ni20(5,1)  =0.380857230419
     ni20(6,1)  =0.37582962279
     ni20(7,1)  =0.371027497524
     ni20(8,1)  =0.366440097012
     ni20(9,1)  =0.362056663641
     ni20(10,1)  =0.3578664398
     ni20(11,1)  =0.353858667878
     ni20(12,1)  =0.350022590264
     ni20(13,1)  =0.346347449346
     ni20(14,1)  =0.342822487512
     ni20(15,1)  =0.339436947153
     ni20(16,1)  =0.336180070656
     ni20(17,1)  =0.33304110041
     ni20(18,1)  =0.330009278803
     ni20(19,1)  =0.327073848225
     ni20(20,1)  =0.324224051063
     ni20(21,1)  =0.321449129708
     ni20(22,1)  =0.318738326547
     ni20(23,1)  =0.316080883969
     ni20(24,1)  =0.313466044362
     ni20(25,1)  =0.310883050117
     ni20(26,1)  =0.30832114362
     ni20(27,1)  =0.305769567261
     ni20(28,1)  =0.303217563428
     ni20(29,1)  =0.300654374511
     ni20(30,1)  =0.298069242897
     ni20(31,1)  =0.295451410976
     ni20(32,1)  =0.292790121136
     ni20(33,1)  =0.290074615766
     ni20(34,1)  =0.287294137255
     ni20(35,1)  =0.28443792799
     ni20(36,1)  =0.281495230362
     ni20(37,1)  =0.278455286758
     ni20(38,1)  =0.275307339568
     ni20(39,1)  =0.272040631179
     ni20(40,1)  =0.268644403981
     ni20(41,1)  =0.265107900363
     ni20(42,1)  =0.261420362712
     ni20(43,1)  =0.257571033418
     ni20(44,1)  =0.253549154869
     ni20(45,1)  =0.249343969455
     ni20(46,1)  =0.244944719563
     ni20(47,1)  =0.240340647583
     ni20(48,1)  =0.235520995902
     ni20(49,1)  =0.230475006911
     ni20(50,1)  =0.225191922996
     ni20(51,1)  =0.219660986548
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