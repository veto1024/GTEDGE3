      subroutine transm (igeom, inon, lmfp, Li, Lj, thetaij, 
     .    alphaj, Lperp, delta_n, Tij)

c/    This subroutine calculates the transmission coefficient
c/    used in the TEP theory for rectangular regions.

c/    New version using direct integration over the angle -phi-
c/    01/13/98, John Mandrekas, GIT
c/    05/25/99, jm, new version using Bickley functions
      
c/    VARIABLES:
c/    ----------
c/    igeom    : flag determining the geometry
c/               igeom = 1 : intersecting sides
c/               igeom = 2 : non-intersecting sides
c/    Li       : length of the -i- side (m)
c/    Lj       : length of the -j- side (m)
c/    thetaij  : angle between the i and j sides (rad)
c/    alphaj   : angle needed in non-intersecting sides (rad)
c/    Lperp    : vertical height (for igeom = 2)
c/    lmfp     : total mean free path (m)

c/    inon     : dummy variable for compatibility with older version
c/    delta_n  : dummy variable for compatibility with older version

c/    Tij      : the value of the transmission coefficient

      implicit none
      include 'geometry.fi'

      data pi /3.1415926/, nph /21/

      integer igeom, nph, inon
      real Li, Lj, Lperp, lmfp, thetaij, alphaj, Tij, ss, pi, delta_n
      real t_ij

      external t_ij

c/    Assign variables in the common block:
c/    ------------------------------------
      nph_pnts = nph
      i_geom = igeom
      L_i = Li
      L_j = Lj
      L_perp = Lperp
      l_mfp = lmfp

      theta_ij = thetaij
      alpha_j = alphaj

      call qgauss_10 (t_ij, 0.0, L_i, ss)

      Tij = 2.0 * ss / (pi*Li)

      return
      end

c///////////////////////////////////////////////////////////////////////

      real function t_ij (x)

c/    This function returns the integral of the transmission coefficient
c/    over the angle phi for a given value of ksi_i = x (see the paper in 
c/    Nuclear Fusion 34 (1994) 1385)

c/    John Mandrekas, GIT, 12/11/97
c/    04/07/99 : Introduced Bickley function 

      implicit none

      include 'geometry.fi'

      integer maxph, i
      parameter (maxph = 51)
      real x, tvec(maxph), sina, cosa, phi_min, phi_max, d_phi, phi, 
     .   ksi_j, el_of_phi, pi, tanth, sinth, reslt
      real bickley

      data pi /3.141592654/

      save pi

c/    Calculate maximum and minimum angles for integration:
c/    ----------------------------------------------------

      tanth = tan (theta_ij)
      sinth = sin (theta_ij)

      if (i_geom.EQ.1) then
	 phi_min = atan (tanth / (1.0 - x * tanth / (L_j * sinth)))
	 phi_max = pi

      else if (i_geom.EQ.2) then

         sina = sin (alpha_j)
         cosa = cos (alpha_j)
         phi_min = atan ((L_perp + L_j * sina) /
     .      (L_perp / tanth  + L_j * cosa - x))
         phi_max = atan (L_perp / (L_perp / tanth - x))

      endif

      if (phi_min.LT.0.0) phi_min = pi + phi_min
      if (phi_max.LT.0.0) phi_max = pi + phi_max

      d_phi = (phi_max - phi_min) / float (nph_pnts - 1)
	
      do i = 1, nph_pnts
	 phi = phi_min + (i-1) * d_phi
	 if (i_geom.EQ.1) then
	    el_of_phi = x * sinth / sin (phi - theta_ij)
	 else
	    ksi_j = ((L_perp / tan(theta_ij) - x) * tan (phi)-L_perp) /
     .         (sina - cosa * tan(phi))
	    el_of_phi = (L_perp + ksi_j*sina)/ sin (phi)
	 endif

	 tvec(i) = sin (phi) * bickley (el_of_phi/l_mfp)

      enddo

c/    Perform integration wrt the angle phi:
c/    -------------------------------------
      call simpson (tvec, d_phi, nph_pnts, reslt)

      t_ij = reslt

      return
      end

c///////////////////////////////////////////////////////////////////////

      subroutine qgauss_10 (func, a, b, ss)
     
c/    This subroutine returns as ss the integral of the function -func-
c/    between -a- and -b-, by 20 point Gauss-Legendre integration. The
c/    function is evaluated exactly twenty times at interior points in
c/    the range of integration.
c/    Written by John Mandrekas, GIT, 08/12/93
     
c/    References:
c/    1) Numerical Recipes, Chapter 4.5
c/    2) Abramowitz and Stegun, Handbook of Mathematical Functions.

      implicit none

      real a, b, ss, func
      external func
      integer i, npoints
      real dx, xm, xr

      parameter (npoints = 10)
      real x_i(npoints), w_i(npoints)
      save x_i, w_i

      data x_i/0.0765265211, 0.2277858511, 0.3737060887, 0.5108670019,
     .         0.6360536807, 0.7463319064, 0.8391169718, 0.9122344282,
     .         0.9639719272, 0.9931285991/

      data w_i/0.1527533871, 0.1491729864, 0.1420961093, 0.1316886384,
     .         0.1181945319, 0.1019301198, 0.0832767415, 0.0626720483,
     .         0.0406014298, 0.0176140071/

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      ss = 0.0
      do i = 1, npoints
         dx = xr * x_i(i)
         ss = ss + w_i(i) * (func(xm+dx) + func(xm-dx))
      enddo
     
      ss = xr * ss

      return
      end                                                              

c///////////////////////////////////////////////////////////////////////

      real function bickley(x)

c/    This is an approximation for the Ki3 Bickley-Naylor
c/    function. 
c/    Ref: F.G. Lether, J. Quant. Spectrosc. Radiat. Transfer,
c/         43 (1990) 187.
c/    Created by John Mandrekas, on 4/7/99 for the GTNEUT code.

      real x, numer, denom, Ki3

      numer = 1.88571 + 1.56855*x + 0.202569*x**2
      denom = 2.40135 + 1.43711*x + 0.161532*x**2

      Ki3 = exp(-x) * numer / (Sqrt(x+1.0) * denom)

      bickley = Ki3

      return
      end

c///////////////////////////////////////////////////////////////////////

      subroutine simpson(f, h, n, s)

c///////////////////////////////////////////////////////////////////////
c/                                                                     /
c/    This subroutine calculates the integral of the function -f-      /
c/    which is defined at -n- discrete points with equal step -h-      /
c/    using  Simpson's rule. If the number of points -n- is even,      /
c/    the trapezoidal rule is used for the last two points.            /
c/    Note: use this routine for cases, where the function is given    /
c/    in discrete points. If an analytic form of the function is       /
c/    available, it's better to use an adaptive scheme.                /
c/    Created by John Mandrekas, GIT, 1/14/92                          /
c/    Cosmetic changes, 06/21/95                                       /
c/                                                                     /
c///////////////////////////////////////////////////////////////////////

      implicit none
      integer i, n
      real s, h, f(n)

      s = 0.0
      if (n.LE.1) return

      do i = 3, n, 2
        s = s + f(i-2) + 4.0 * f(i-1) + f(i)
      end do
      s = s / 3.0

      if (mod(n,2).EQ.0) then
        s = s + 0.5 * (f(n-1) + f(n))
      endif

      s = s * h

      return
      end
