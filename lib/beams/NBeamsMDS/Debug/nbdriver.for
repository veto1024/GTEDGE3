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
     tekev(1)  =1.67623728527
     tekev(2)  =1.65318034514
     tekev(3)  =1.62856585431
     tekev(4)  =1.60246896785
     tekev(5)  =1.57496484084
     tekev(6)  =1.54612862834
     tekev(7)  =1.51603548542
     tekev(8)  =1.48476056717
     tekev(9)  =1.45237902864
     tekev(10)  =1.41896602491
     tekev(11)  =1.38459671104
     tekev(12)  =1.34934624212
     tekev(13)  =1.31328977322
     tekev(14)  =1.27650245939
     tekev(15)  =1.23905945573
     tekev(16)  =1.20103591729
     tekev(17)  =1.16250699914
     tekev(18)  =1.12354785637
     tekev(19)  =1.08423364403
     tekev(20)  =1.04463951721
     tekev(21)  =1.00484063096
     tekev(22)  =0.964912140375
     tekev(23)  =0.924929200512
     tekev(24)  =0.884966966444
     tekev(25)  =0.845100593243
     tekev(26)  =0.805405235981
     tekev(27)  =0.765956049728
     tekev(28)  =0.726828189556
     tekev(29)  =0.688096810536
     tekev(30)  =0.649837067739
     tekev(31)  =0.612124116237
     tekev(32)  =0.5750331111
     tekev(33)  =0.5386392074
     tekev(34)  =0.503017560208
     tekev(35)  =0.468243324595
     tekev(36)  =0.434391655632
     tekev(37)  =0.401537708392
     tekev(38)  =0.369756637944
     tekev(39)  =0.33912359936
     tekev(40)  =0.309713747711
     tekev(41)  =0.281602238069
     tekev(42)  =0.254864225505
     tekev(43)  =0.229574865089
     tekev(44)  =0.205809311893
     tekev(45)  =0.183642720989
     tekev(46)  =0.163150247447
     tekev(47)  =0.144407046339
     tekev(48)  =0.127488272736
     tekev(49)  =0.112469081709
     tekev(50)  =0.0994246283296
     tekev(51)  =0.0884300676685
     tikev(1)  =1.71170373287
     tikev(2)  =1.6615404356
     tikev(3)  =1.61216319797
     tikev(4)  =1.56357389457
     tikev(5)  =1.51577439997
     tikev(6)  =1.46876658876
     tikev(7)  =1.42255233552
     tikev(8)  =1.37713351484
     tikev(9)  =1.3325120013
     tikev(10)  =1.28868966948
     tikev(11)  =1.24566839396
     tikev(12)  =1.20345004932
     tikev(13)  =1.16203651015
     tikev(14)  =1.12142965104
     tikev(15)  =1.08163134656
     tikev(16)  =1.04264347129
     tikev(17)  =1.00446789982
     tikev(18)  =0.967106506733
     tikev(19)  =0.93056116661
     tikev(20)  =0.894833754032
     tikev(21)  =0.859926143584
     tikev(22)  =0.825840209848
     tikev(23)  =0.792577827405
     tikev(24)  =0.76014087084
     tikev(25)  =0.728531214734
     tikev(26)  =0.69775073367
     tikev(27)  =0.667801302231
     tikev(28)  =0.638684794999
     tikev(29)  =0.610403086557
     tikev(30)  =0.582958051487
     tikev(31)  =0.556351564373
     tikev(32)  =0.530585499797
     tikev(33)  =0.50566173234
     tikev(34)  =0.481582136587
     tikev(35)  =0.45834858712
     tikev(36)  =0.435962958521
     tikev(37)  =0.414427125373
     tikev(38)  =0.393742962258
     tikev(39)  =0.373912343759
     tikev(40)  =0.35493714446
     tikev(41)  =0.336819238941
     tikev(42)  =0.319560501787
     tikev(43)  =0.303162807579
     tikev(44)  =0.2876280309
     tikev(45)  =0.272958046333
     tikev(46)  =0.259154728461
     tikev(47)  =0.246219951865
     tikev(48)  =0.23415559113
     tikev(49)  =0.222963520836
     tikev(50)  =0.212645615568
     tikev(51)  =0.203203749907
     ne20(1)  =0.507331172645
     ne20(2)  =0.495268198044
     ne20(3)  =0.483715525787
     ne20(4)  =0.472652476628
     ne20(5)  =0.462058371321
     ne20(6)  =0.451912530621
     ne20(7)  =0.442194275281
     ne20(8)  =0.432882926057
     ne20(9)  =0.423957803701
     ne20(10)  =0.415398228969
     ne20(11)  =0.407183522614
     ne20(12)  =0.399293005391
     ne20(13)  =0.391705998054
     ne20(14)  =0.384401821358
     ne20(15)  =0.377359796055
     ne20(16)  =0.370559242902
     ne20(17)  =0.363979482651
     ne20(18)  =0.357599836057
     ne20(19)  =0.351399623875
     ne20(20)  =0.345358166858
     ne20(21)  =0.339454785761
     ne20(22)  =0.333668801338
     ne20(23)  =0.327979534343
     ne20(24)  =0.32236630553
     ne20(25)  =0.316808435654
     ne20(26)  =0.311285245469
     ne20(27)  =0.305776055728
     ne20(28)  =0.300260187187
     ne20(29)  =0.2947169606
     ne20(30)  =0.28912569672
     ne20(31)  =0.283465716302
     ne20(32)  =0.2777163401
     ne20(33)  =0.271856888868
     ne20(34)  =0.265866683361
     ne20(35)  =0.259725044333
     ne20(36)  =0.253411292538
     ne20(37)  =0.246904748729
     ne20(38)  =0.240184733663
     ne20(39)  =0.233230568091
     ne20(40)  =0.22602157277
     ne20(41)  =0.218537068453
     ne20(42)  =0.210756375894
     ne20(43)  =0.202658815847
     ne20(44)  =0.194223709067
     ne20(45)  =0.185430376309
     ne20(46)  =0.176258138325
     ne20(47)  =0.16668631587
     ne20(48)  =0.156694229699
     ne20(49)  =0.146261200566
     ne20(50)  =0.135366549225
     ne20(51)  =0.12398959643
     ni20(1,1)  =0.492679625148
     ni20(2,1)  =0.478326944695
     ni20(3,1)  =0.464750101088
     ni20(4,1)  =0.451917671779
     ni20(5,1)  =0.439798234222
     ni20(6,1)  =0.428360365868
     ni20(7,1)  =0.417572644171
     ni20(8,1)  =0.407403646582
     ni20(9,1)  =0.397821950556
     ni20(10,1)  =0.388796133544
     ni20(11,1)  =0.380294772998
     ni20(12,1)  =0.372286446373
     ni20(13,1)  =0.364739731119
     ni20(14,1)  =0.357623204691
     ni20(15,1)  =0.350905444541
     ni20(16,1)  =0.34455502812
     ni20(17,1)  =0.338540532883
     ni20(18,1)  =0.332830536281
     ni20(19,1)  =0.327393615768
     ni20(20,1)  =0.322198348795
     ni20(21,1)  =0.317213312816
     ni20(22,1)  =0.312407085283
     ni20(23,1)  =0.30774824365
     ni20(24,1)  =0.303205365367
     ni20(25,1)  =0.298747027889
     ni20(26,1)  =0.294341808668
     ni20(27,1)  =0.289958285157
     ni20(28,1)  =0.285565034807
     ni20(29,1)  =0.281130635073
     ni20(30,1)  =0.276623663406
     ni20(31,1)  =0.272012697259
     ni20(32,1)  =0.267266314085
     ni20(33,1)  =0.262353091337
     ni20(34,1)  =0.257241606466
     ni20(35,1)  =0.251900436927
     ni20(36,1)  =0.24629816017
     ni20(37,1)  =0.24040335365
     ni20(38,1)  =0.234184594819
     ni20(39,1)  =0.227610461129
     ni20(40,1)  =0.220649530034
     ni20(41,1)  =0.213270378985
     ni20(42,1)  =0.205441585435
     ni20(43,1)  =0.197131726838
     ni20(44,1)  =0.188309380645
     ni20(45,1)  =0.17894312431
     ni20(46,1)  =0.169001535284
     ni20(47,1)  =0.158453191021
     ni20(48,1)  =0.147266668974
     ni20(49,1)  =0.135410546595
     ni20(50,1)  =0.122853401336
     ni20(51,1)  =0.109563810651
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