      real function sigmav(t)
      
c/    This function calculates the reactivity for the T(d,n)4He reaction
c/    as a function of the ion temperature, using the formalism by
c/    Bosch and Hale (Nuclear Fusion, 32, 611, 1992).
c/    The reactivity is valid for 0.2 keV < T_i < 100 keV

c/    t      : ion temperature in keV
c/    sigmav : <sv>_DT in m^3 / s  
      
      implicit none
      real t, BG, mc2, c1, c2, c3, c4, c5, c6, c7, theta, ksi, term1,
     .   term2 
      
      data BG/34.3827/, mc2/1124656.00/, c1/1.17302e-09/, 
     .     c2/1.51361e-02/, c3/7.51886e-02/, c4/4.60643e-03/, 
     .     c5/1.35000e-02/, c6/-1.0675e-04/, c7/1.36600e-05/ 
     
      if (t.le.0.0) then
        sigmav = 0.0                        
        
        return
      endif 
      
      term1 = t * (c2 + t*(c4 + t*c6))
      term2 = 1.0 + t * (c3 + t*(c5 + t*c7))
      theta = t / (1.0 - term1 / term2)
      ksi = (0.25 * BG**2 / theta)** (0.3333334)
      
      sigmav = 1.0e-6*c1 * theta * sqrt(ksi/(mc2*t**3)) * exp(-3.0*ksi)
      
      return
      end
