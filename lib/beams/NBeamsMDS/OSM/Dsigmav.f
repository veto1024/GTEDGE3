      subroutine Dsigmav(t, svDT, dsvDT)
      
c/    This function calculates the reactivity for the T(d,n)4He reaction
c/    as a function of the ion temperature, using the formalism by
c/    Bosch and Hale (Nuclear Fusion, 32, 611, 1992).
c/    The reactivity is valid for 0.2 keV < T_i < 100 keV
c/  
c/    11/28/2001, jm: Converted to a subroutine and calculation of the
c/                    derivative d(sv) / dT was added.

c/    VARIABLES:
c/    ---------
c/    t      : ion temperature in keV
c/    svDT   : <sv>_DT in m^3 / s  
c/    dsvDT  : d(<sv>_DT) / dT in m^3 / (s * keV)
      
      implicit none
      real t, svDT, dsvDT
      real BG, mc2, c1, c2, c3, c4, c5, c6, c7, theta, ksi, term1,
     .   term2, numer, denom, dnum, dden, dthtdt, dksidt, dfksdt      
      
      data BG/34.3827/, mc2/1124656.00/, c1/1.17302e-09/, 
     .     c2/1.51361e-02/, c3/7.51886e-02/, c4/4.60643e-03/, 
     .     c5/1.35000e-02/, c6/-1.0675e-04/, c7/1.36600e-05/ 
     
      if (t.le.0.0) then
        svDT = 0.0                        
        dsvDT = 0.0  
        return
      endif 
      
      term1 = t * (c2 + t*(c4 + t*c6))
      term2 = 1.0 + t * (c3 + t*(c5 + t*c7))
      theta = t / (1.0 - term1 / term2)
      ksi = (0.25 * BG**2 / theta)** (0.3333334)
      
      svDT = 1.0e-6*c1 * theta * sqrt(ksi/(mc2*t**3)) * exp(-3.0*ksi)

c/    Calculation of the derivative d(sv_DT) / dT
      
c/    d Theta / dT:

      numer = t*term2
      denom = term2 - term1

      dnum = 1. + 2.*c3*t + 3.*c5*t**2 + 4.*c7*t**3
      dden = c3-c2 + 2.*(c5-c4)*t + 3.*(c7-c6)*t**2

      dthtdt = (dnum*denom - numer*dden) / denom**2

c/    d(ksi) / dT:
      
      dksidt = -0.21 * (BG**(0.66667) / theta**(1.33334)) * dthtdt
      
c/    d(Sqrt(ksi / mc2*T**3)) / dT:

      dfksdt = (1./(t**3 * sqrt(mc2))) * 
     .         (0.5*t**1.5 * dksidt / sqrt(ksi) - 1.5*sqrt(ksi*t))

      dsvDT = 1.0e-6 * c1 * exp(-3.0*ksi) * 
     .         (-3.0*theta*sqrt(ksi/(mc2*t**3))*dksidt +
     .          sqrt(ksi/(mc2*t**3)) * dthtdt +
     .          theta * dfksdt)

      return
      end
