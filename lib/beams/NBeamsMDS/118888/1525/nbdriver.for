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
	   tekev(1) =   1.4035515
	   tekev(2) =   1.4015922
	   tekev(3) =   1.3992248
	   tekev(4) =   1.3960413
	   tekev(5) =   1.3916335
	   tekev(6) =   1.3855938
	   tekev(7) =   1.3776222
	   tekev(8) =   1.3676971
	   tekev(9) =   1.3560577
	   tekev(10) =  1.3429997
	   tekev(11) =  1.3280147
	   tekev(12) =  1.3116049
	   tekev(13) =  1.2933756
	   tekev(14) =  1.2735166
	   tekev(15) =  1.2520678
	   tekev(16) =  1.2289563
	   tekev(17) =  1.2042476
	   tekev(18) =  1.1779946
	   tekev(19) =  1.1502484
	   tekev(20) =  1.1210557 
	   tekev(21) =  1.0904669
	   tekev(22) =  1.0585413
	   tekev(23) =  1.0253449
	   tekev(24) = 0.99096826
	   tekev(25) = 0.95550245
	   tekev(26) = 0.91901803
	   tekev(27) = 0.88163342
	   tekev(28) = 0.84344325
	   tekev(29) = 0.80455955
	   tekev(30) = 0.76510234
	   tekev(31) = 0.72519692
	   tekev(32) = 0.68497466
	   tekev(33) = 0.64457497
	   tekev(34) = 0.60414240
	   tekev(35) = 0.56382931
	   tekev(36) = 0.52379434
	   tekev(37) = 0.48420148
	   tekev(38) = 0.44522227
	   tekev(39) = 0.40703418
	   tekev(40) = 0.36980948
	   tekev(41) = 0.33372934
	   tekev(42) = 0.29895573
	   tekev(43) = 0.26561571
	   tekev(44) = 0.23376134
	   tekev(45) = 0.20328085
	   tekev(46) = 0.17376321
	   tekev(47) = 0.14433722
	   tekev(48) = 0.11384351
	   tekev(49) = 0.082274866
	   tekev(50) = 0.053358432
	   tekev(51) = 0.032876627
	   tikev(1) =  1.4744065
	   tikev(2) =  1.4694183
	   tikev(3) =  1.4617613
	   tikev(4) =  1.4488300
	   tikev(5) =  1.4298830
	   tikev(6) =  1.4069735
	   tikev(7) =  1.3802957
	   tikev(8) =  1.3504733
	   tikev(9) =  1.3180368
	   tikev(10) = 1.2834918
	   tikev(11) = 1.2473506
	   tikev(12) = 1.2101304
	   tikev(13) = 1.1723612
	   tikev(14) = 1.1345533
	   tikev(15) = 1.0972252
	   tikev(16) = 1.0609191
	   tikev(17) = 1.0260957
	   tikev(18) = 0.99331895
	   tikev(19) = 0.96311476
	   tikev(20) = 0.93599228 
	   tikev(21) = 0.91247589
	   tikev(22) = 0.89241746
	   tikev(23) = 0.87498469
	   tikev(24) = 0.85910217
	   tikev(25) = 0.84364546
	   tikev(26) = 0.82762628
	   tikev(27) = 0.80998197
	   tikev(28) = 0.78965333
	   tikev(29) = 0.76584791
	   tikev(30) = 0.73878168
	   tikev(31) = 0.70890233
	   tikev(32) = 0.67676272
	   tikev(33) = 0.64283641
	   tikev(34) = 0.60763114
	   tikev(35) = 0.57164221
	   tikev(36) = 0.53537195
	   tikev(37) = 0.49932178
	   tikev(38) = 0.46399096
	   tikev(39) = 0.42988234
	   tikev(40) = 0.39749639
	   tikev(41) = 0.36723137
	   tikev(42) = 0.33912149
	   tikev(43) = 0.31311932
	   tikev(44) = 0.28917450
	   tikev(45) = 0.26724357
	   tikev(46) = 0.24727921
	   tikev(47) = 0.22923775
	   tikev(48) = 0.21306653
	   tikev(49) = 0.19871620
	   tikev(50) = 0.18615186
	   tikev(51) = 0.17531932
	   ne20(1) =  0.47999990
	   ne20(2) =  0.47858268
	   ne20(3) =  0.47663438
	   ne20(4) =  0.47362393
	   ne20(5) =  0.46918200
	   ne20(6) =  0.46355579
	   ne20(7) =  0.45724135
	   ne20(8) =  0.45017154
	   ne20(9) =  0.44250569
	   ne20(10) = 0.43440458
	   ne20(11) = 0.42592987
	   ne20(12) = 0.41718507
	   ne20(13) = 0.40826456
	   ne20(14) = 0.39925334
	   ne20(15) = 0.39022778
	   ne20(16) = 0.38125659
	   ne20(17) = 0.37239961
	   ne20(18) = 0.36371072
	   ne20(19) = 0.35523762
	   ne20(20) = 0.34701966 
	   ne20(21) = 0.33908966
	   ne20(22) = 0.33147215
	   ne20(23) = 0.32418599
	   ne20(24) = 0.31724209
	   ne20(25) = 0.31064515
	   ne20(26) = 0.30437410
	   ne20(27) = 0.29842925
	   ne20(28) = 0.29279699
	   ne20(29) = 0.28745173
	   ne20(30) = 0.28233921
	   ne20(31) = 0.27743405
	   ne20(32) = 0.27268971
	   ne20(33) = 0.26804170
	   ne20(34) = 0.26342902
	   ne20(35) = 0.25877935
	   ne20(36) = 0.25401157
	   ne20(37) = 0.24905419
	   ne20(38) = 0.24380050
	   ne20(39) = 0.23815152
	   ne20(40) = 0.23199912
	   ne20(41) = 0.22522464
	   ne20(42) = 0.21772676
	   ne20(43) = 0.20934472
	   ne20(44) = 0.19996840
	   ne20(45) = 0.18943805
	   ne20(46) = 0.17760940
	   ne20(47) = 0.16431446
	   ne20(48) = 0.14937444
	   ne20(49) = 0.13240343
	   ne20(50) = 0.11091435
	   ne20(51) = 0.071682116
	   ni20(1,1) =  0.44365256
	   ni20(2,1) =  0.44296567
	   ni20(3,1) =  0.44167145
	   ni20(4,1) =  0.43916255
	   ni20(5,1) =  0.43505601
	   ni20(6,1) =  0.42981577
	   ni20(7,1) =  0.42407326
	   ni20(8,1) =  0.41791733
	   ni20(9,1) =  0.41148520
	   ni20(10,1) = 0.40469599
	   ni20(11,1) = 0.39760019
	   ni20(12,1) = 0.39007643
	   ni20(13,1) = 0.38208848
	   ni20(14,1) = 0.37384468
	   ni20(15,1) = 0.36567716
	   ni20(16,1) = 0.35768364
	   ni20(17,1) = 0.34982911
	   ni20(18,1) = 0.34209535
	   ni20(19,1) = 0.33453920
	   ni20(20,1) = 0.32722218 
	   ni20(21,1) = 0.32018163
	   ni20(22,1) = 0.31340703
	   ni20(23,1) = 0.30689533
	   ni20(24,1) = 0.30066197
	   ni20(25,1) = 0.29472054
	   ni20(26,1) = 0.28905528
	   ni20(27,1) = 0.28363610
	   ni20(28,1) = 0.27843272
	   ni20(29,1) = 0.27343297
	   ni20(30,1) = 0.26863396
	   ni20(31,1) = 0.26400975
	   ni20(32,1) = 0.25949432
	   ni20(33,1) = 0.25500222
	   ni20(34,1) = 0.25048216
	   ni20(35,1) = 0.24587890
	   ni20(36,1) = 0.24112351
	   ni20(37,1) = 0.23613007
	   ni20(38,1) = 0.23077744
	   ni20(39,1) = 0.22496995
	   ni20(40,1) = 0.21860704
	   ni20(41,1) = 0.21159919
	   ni20(42,1) = 0.20392509
	   ni20(43,1) = 0.19551640
	   ni20(44,1) = 0.18632554
	   ni20(45,1) = 0.17623515
	   ni20(46,1) = 0.16511768
	   ni20(47,1) = 0.15277370
	   ni20(48,1) = 0.13898608
	   ni20(49,1) = 0.12334245
	   ni20(50,1) = 0.10347334
	   ni20(51,1) = 0.066970565
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
