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