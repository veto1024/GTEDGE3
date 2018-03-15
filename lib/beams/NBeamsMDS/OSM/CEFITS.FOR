      subroutine cefits (iz, te, celz, ider, dlzdt)

c/    This subroutine calculates the corona equilibrium radiative
c/    rates for some elements of interest using fits to ADPAK 
c/    calculations. The fits have been made using JANDEL's TableCurve
c/    program
c/    John Mandrekas, 10-22-95, GIT
c/    10/08/96, jm: Added data for Carbon
c/    04/16/97, jm: Extended Carbon range to lower temperatures
c/    04/16/97, jm: Added data for Oxygen (5 regions)
c/    04/17/97, jm: Added data for Fe (3 regions)
c/    04/17/97, jm: Added option not to calculate derivatives
c/    04/17/97, jm: Extended Ar range to 5 eV

c/    PARAMETERS
c/    ----------
c/    iz      : atomic number of impurity species
c/    te      : electron temperature (eV)
c/    celz    : corona equilibrium radiative rate (erg-cm^3 /s)
c/    ider    : if equal to 1, calculate the derivative
c/    dlzdt   : the derivative dLz/dTe

c/    Notice that all  parameters should be double precision on
c/    32-bit machines.

c/    The subroutine currently contains data for the following
c/    species:

c/    Z        Element
c/    -        -------
c/     6       Carbon
c/     8       Oxygen
c/    10       Neon
c/    18       Argon
c/    26       Iron
c/    36       Krypton

c/    If the element requested is not included, the routine exits
c/    with an error message.
c/    If te_ev is less than the minimum temperature supported for
c/    species iz, it is set equal to the minimum temperature.
c/    If te_ev is larger than the maximum  temperature supported for
c/    this species, the Lz is scaled as  sqrt(Te/Te_max) i.e., pure
c/    bremmstrahlung.


      integer iz, ibrem, ider
      double precision te, te_ev, celz, dlzdt, dte, te_min, te_max,
     .   tefix
      double precision neon1, neon2, neon3, argon1, argon2, argon3,
     .   krypton1, krypton2, krypton3, carbon0, carbon1, carbon2,
     .   oxy1, oxy2, oxy3, oxy4, oxy5, iron1, iron2, iron3
      external neon1, neon2, neon3, argon1, argon2, argon3, krypton1,
     .   krypton2, krypton3, carbon0, carbon1, carbon2, oxy1, oxy2,
     .   oxy3, oxy4, oxy5, iron1, iron2, iron3

      te_ev = te
      ibrem = 0


c/    CARBON
c/    ======
      if (iz.eq.6) then
	if (te_ev.LT.3.0) te_ev = 3.0
	if (te_ev.GT.10000.0) then
	   te_true = te_ev
	   te_ev = 10000.00
	   tefix = te_ev
           ibrem = 1
        endif

c/    Carbon range 0: 3 eV to 10 eV:
c/    --------------------------------

        if ((te_ev.GE.3.0).AND.(te_ev.LE.10.00)) then
	   te_min = 3.0
	   te_max = 10.0
	   dte = 0.1
	   celz = carbon0 (te_ev)

c/    Calculate derivative dL_z / dT_e
           if (ider.NE.0) 
     .      call calcderiv (carbon0, te_ev, te_max, te_min, dte, dlzdt)

c/    Carbon range 1: 10 eV to 100 eV:
c/    --------------------------------

        else if ((te_ev.GE.10.0).AND.(te_ev.LE.100.00)) then
	   te_min = 10.0
	   te_max = 100.0
	   dte = 0.5
	   celz = carbon1 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (carbon1, te_ev, te_max, te_min, dte, dlzdt)

c/    Carbon range 2: 100 eV to 10,000 eV:
c/    -----------------------------------

        else if ((te_ev.GT.100.0).AND.(te_ev.LE.1.0e4)) then
	   te_min = 100.0
	   te_max = 10000.00
	   celz = carbon2 (te_ev)
	   if (te_ev.GT.100.0.AND.te_ev.LE.120.0) dte = 0.5
	   if (te_ev.GT.120.0.AND.te_ev.LE.1000.0) dte = 5.0
	   if (te_ev.GT.1000.0.AND.te_ev.LE.10000.0) dte = 10.00

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (carbon2, te_ev, te_max, te_min, dte, dlzdt)

        endif

c/    OXYGEN 
c/    =======
      else if (iz.eq.8) then
	if (te_ev.LT.5.0) te_ev = 5.0
	if (te_ev.GT.50000.0) then
	   te_true = te_ev
	   te_ev = 50000.00
	   tefix = te_ev
           ibrem = 1
        endif

c/    Oxygen range 1: 5 eV to 25 eV:
c/    --------------------------------

        if ((te_ev.GE.5.0).AND.(te_ev.LT.25.00)) then
	   te_min = 5.0
	   te_max = 25.0
	   dte = 1.0
	   celz = oxy1 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (oxy1, te_ev, te_max, te_min, dte, dlzdt)

c/    Oxygen range 2: 25 eV to 100 eV:
c/    --------------------------------

        else if ((te_ev.GE.25.0).AND.(te_ev.LT.100.0)) then
	   te_min = 25.0
	   te_max = 100.0
	   dte = 5.0
	   celz = oxy2 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (oxy2, te_ev, te_max, te_min, dte, dlzdt)

c/    Oxygen range 3: 100 eV to 3 keV:
c/    --------------------------------

        else if ((te_ev.GE.100.0).AND.(te_ev.LT.3.0e3)) then
	   te_min = 100.0
	   te_max = 3000.0
	   dte = 10.0
	   celz = oxy3 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (oxy3, te_ev, te_max, te_min, dte, dlzdt)

c/    Oxygen range 4: 3 keV to 10  keV:
c/    --------------------------------

        else if ((te_ev.GE.3.0e3).AND.(te_ev.LT.1.0e4)) then
	   te_min = 3.0e3
	   te_max = 1.0e4
	   dte = 100.0
	   celz = oxy4 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .       call calcderiv (oxy4, te_ev, te_max, te_min, dte, dlzdt)

c/    Oxygen range 5: 10 keV to 50  keV:
c/    --------------------------------

        else if ((te_ev.GE.1.0e4).AND.(te_ev.LT.5.0e4)) then
	   te_min = 1.0e4
	   te_max = 5.0e4
	   dte = 1000.0
	   celz = oxy5 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (oxy5, te_ev, te_max, te_min, dte, dlzdt)

        endif

c/    NEON
c/    ====
      else if (iz.eq.10) then
	if (te_ev.LT.10.0) te_ev = 10.0
	if (te_ev.GT.10000.0) then
	   te_true = te_ev
	   te_ev = 10000.0
	   tefix = te_ev
           ibrem = 1
        endif

c/    Neon range 1: 10 eV to 100 eV:
c     ==============================

        if ((te_ev.GE.10.0).AND.(te_ev.LE.100.00)) then
	   te_min = 10.0
	   te_max = 100.0
	   dte = 10.0
	   celz = neon1 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (neon1, te_ev, te_max, te_min, dte, dlzdt)

c/    Neon range 2: 100 eV to 1000 eV:
c     ================================

        else if ((te_ev.GT.100.0).AND.(te_ev.LE.1000.00)) then
	   te_min = 100.0
	   te_max = 1000.0
	   dte = 10.0
	   celz = neon2 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (neon2, te_ev, te_max, te_min, dte, dlzdt)

c/    Neon range 3: 1000 eV to 10000 eV:
c     ==================================

        else if ((te_ev.GT.1000.0).AND.(te_ev.LE.10000.00)) then
	   te_min = 1000.0
	   te_max = 10000.0
	   dte = 100.0
	   celz = neon3 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (neon3, te_ev, te_max, te_min, dte, dlzdt)

        endif

c/    ARGON
c/    =====
      else if (iz.eq.18) then
	if (te_ev.LT.5.0) te_ev = 5.0
	if (te_ev.GT.1.0e5) then
	   te_true = te_ev
	   te_ev = 1.0e5
	   tefix = te_ev
           ibrem = 1
        endif

c/    Argon range 1: 5.0 eV to 100 eV:
c     ===============================

        if ((te_ev.GE.5.0).AND.(te_ev.LE.100.00)) then
	   te_min = 5.0
	   te_max = 100.0
	   dte = 1.0
	   celz = argon1 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (argon1, te_ev, te_max, te_min, dte, dlzdt)

c/    Argon range 2: 100 eV to 1000 eV:
c     =================================

        else if ((te_ev.GT.100.0).AND.(te_ev.LE.1000.00)) then
	   te_min = 100.0
	   te_max = 1000.0
	   dte = 10.0
	   celz = argon2 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (argon2, te_ev, te_max, te_min, dte, dlzdt)

c/    Argon range 3: 1000 eV to 100000 eV:
c     ===================================

        else if ((te_ev.GT.1000.0).AND.(te_ev.LE.1.0e5)) then
	   te_min = 1000.0
	   te_max = 1.0e5
	   dte = 50.0
	   celz = argon3 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (argon3, te_ev, te_max, te_min, dte, dlzdt)

        endif

c/    IRON 
c/    ====
      else if (iz.eq.26) then
	if (te_ev.LT.10.0) te_ev = 10.0
	if (te_ev.GT.1.0e5) then
	   te_true = te_ev
	   te_ev = 1.0e5
	   tefix = te_ev
           ibrem = 1
        endif

c/    Iron range 1: 10 eV to 1000 eV:
c     ==============================

        if ((te_ev.GE.10.0).AND.(te_ev.LE.1000.00)) then
	   te_min = 10.0
	   te_max = 1000.0
	   dte = 1.0
	   celz = iron1 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .       call calcderiv (iron1, te_ev, te_max, te_min, dte, dlzdt)

c/    Iron range 2: 1 keV to 10 keV:
c/    ==============================
        else if ((te_ev.GE.1.0e3).AND.(te_ev.LE.1.0e4)) then
	   te_min = 1.0e3
	   te_max = 1.0e4
	   dte = 100.0
	   celz = iron2 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .       call calcderiv (iron2, te_ev, te_max, te_min, dte, dlzdt)

c/    Iron range 3: 10 keV to 100 keV:
c/    ================================
        else if ((te_ev.GE.1.0e4).AND.(te_ev.LE.1.0e5)) then
	   te_min = 1.0e4
	   te_max = 1.0e5
	   dte = 1000.0
	   celz = iron3 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .       call calcderiv (iron3, te_ev, te_max, te_min, dte, dlzdt)

	endif

c/    KRYPTON
c/    =======
      else if (iz.eq.36) then
	if (te_ev.LT.50.0) te_ev = 50.0
	if (te_ev.GT.15000.0) then
	   te_true = te_ev
	   te_ev = 15000.0
	   tefix = te_ev
           ibrem = 1
        endif

c/    Krypton range 1: 50 eV to 1000 eV:
c     ==================================

        if ((te_ev.GE.50.0).AND.(te_ev.LE.1000.00)) then
	   te_min = 50.0
	   te_max = 1000.0
	   dte = 10.0
	   celz = krypton1 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (krypton1, te_ev, te_max, te_min, dte, dlzdt)

c/    Krypton range 2: 1000 eV to 5000 eV:
c     ====================================

        else if ((te_ev.GT.1000.0).AND.(te_ev.LE.5000.00)) then
	   te_min = 1000.0
	   te_max = 5000.0
	   dte = 100.0
	   celz = krypton2 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (krypton2, te_ev, te_max, te_min, dte, dlzdt)

c/    Krypton range 3: 5000 eV to 15000 eV:
c     =====================================

        else if ((te_ev.GT.5000.0).AND.(te_ev.LE.15000.00)) then
	   te_min = 5000.0
	   te_max = 15000.0
	   dte = 100.0
	   celz = krypton3 (te_ev)

c/    Calculate derivative dL_z / dT_e

           if (ider.NE.0) 
     .      call calcderiv (krypton3, te_ev, te_max, te_min, dte, dlzdt)

        endif

      endif

c/    If te_ev > Te_max, extrapolate as bremsstrahlung:
        
      if (ibrem.eq.1) then 
	 dlzdt = 0.5 * celz / sqrt(te_true * tefix)
	 celz = celz * sqrt (te_true / tefix)
      endif
      return
      end

c///////////////////////////////////////////////////////////////////////

      subroutine calcderiv (func, te_ev, te_max, te_min, dte, dlzdt)

c/    This subroutine calculates the derivative dL_z/dT_e using
c/    the appropriate function that is passed through the calling
c/    statement.

c/    func    : the fitting function for the particular impurity
c/              species and temperature range.
c/    te_ev   : the temperature at which the derivative is desired
c/    te_min, : the minimum temperature of this range
c/    te_max  : the maximum temperature of this range
c/    dte     : the delta_Te for this region
c/    dlzdt   : the derivative dLz / dTe

      double precision func, te_ev, te_max, te_min, dte, dlzdt
      external func

      if ((te_ev + dte).GT.te_max) then
         dlzdt = (func(te_ev) - func(te_ev - dte)) / dte
      else if ((te_ev - dte).LT.te_min) then
         dlzdt = (func(te_ev + dte) - func(te_ev)) / dte
      else
         dlzdt = (func(te_ev+dte) - func(te_ev-dte)) /
     .       (2.0 * dte)
      endif
        
      return
      end

c//////////////////// FUNCTION NEON1 ///////////////////////////////////

      double precision function neon1 (x)

C**** NEON (10 - 100 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 6113  y^(-1)=a+bx+cx^2+dx^3+ex^4+fx^5
C**** a=  1.231858377047056E+19
C**** b= -1.167122564478018E+18
C**** c=  5.061989206828907E+16
C**** d= -1115566359746292
C**** e=  11958368146564.65
C**** f= -45752515431.83521

      DOUBLE PRECISION X, Y
      Y = 1.231858377047056E+19 + X * (-1.167122564478018E+18+X*(
     .    5.061989206828907E+16+X*(-1115566359746292.+X*(
     .    11958368146564.65+X*(-45752515431.83521)))))
      Y = 1 / Y

      neon1 = y

      RETURN
      END

c//////////////////// FUNCTION NEON2 ///////////////////////////////////

      double precision function neon2 (x)

C**** NEON (100 - 1000 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 6206  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4+jx^5+k/x^5
C**** a=  6.125085947534055E-18
C**** b= -1.371766505945192E-20
C**** c= -1.679205173830267E-15
C**** d=  1.908988537733591E-23
C**** e=  2.807508219884995E-13
C**** f= -1.611537514888576E-26
C**** g= -2.775589396662275E-11
C**** h=  7.568155225631293E-30
C**** i=  1.493103633158381E-09
C**** j= -1.518756042183371E-33
C**** k= -3.353847345183818E-08

      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (-1.371766505945192E-20+X*(1.908988537733591E-23+X*(
     .    -1.611537514888576E-26+X*(7.568155225631293E-30+X*(-
     .    1.518756042183371E-33))))) + 6.125085947534055E-18 + N * (
     .    -1.679205173830267E-15+N*(2.807508219884995E-13+N*(-
     .    2.775589396662275E-11+N*(1.493103633158381E-09+N*(-
     .    3.353847345183818E-08)))))

      neon2 = y

      RETURN
      END
 
c//////////////////// FUNCTION NEON3 ///////////////////////////////////

      double precision function neon3 (x)

C**** NEON (1 - 10 keV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 6204  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4
C**** a=  3.026448719747039E-20
C**** b= -6.628024598014404E-24
C**** c= -6.909982558373779E-17
C**** d=  8.889358233529658E-28
C**** e=  1.029616271505005E-13
C**** f= -6.151825334164403E-32
C**** g= -7.203273172627677E-11
C**** h=  1.711992089539565E-36
C**** i=  2.031978973186038E-08

      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (-6.628024598014404E-24+X*(8.889358233529658E-28+X*(
     .    -6.151825334164403E-32+X*(1.711992089539565E-36)))) + 
     .    3.026448719747039E-20 + N * (-6.909982558373779E-17+N*(
     .    1.029616271505005E-13+N*(-7.203273172627677E-11+N*(
     .    2.031978973186038E-08))))

      neon3 = y

      RETURN
      END
 
c//////////////////// FUNCTION ARGON1 //////////////////////////////////


      DOUBLE PRECISION FUNCTION ARGON1(X)
C----------------------------------------------------------*
C**** Argon (5 - 100 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 7104  lny=(a+cx+ex^2)/(1+bx+dx^2+fx^3)
C**** r2=0.9906022626695329
C**** r2adj=0.9900024070952477
C**** StdErr=5.737292113813911E-20
C**** Fval=2002.763253416617
C**** a= -42.65475553798339
C**** b= -0.04406746013666619
C**** c= 1.976572632183478
C**** d= 0.0007937234001294226
C**** e= -0.03710404474982667
C**** f= 4.823660058271981E-07
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-42.65475553798339+X*(1.976572632183478+X*(-
     .    0.03710404474982667))) / (1.0+X*(-0.04406746013666619+X*(
     .    0.0007937234001294226+X*(4.823660058271981E-07))))
      Y = DEXP(Y)
      ARGON1 = Y
      RETURN
      END
 
c//////////////////// FUNCTION ARGON2 //////////////////////////////////

      DOUBLE PRECISION FUNCTION ARGON2(X)

C**** ARGON (100 - 1000 eV)

C**** TableCurve D:\WORKSPC\ARGONFIT\ARGON2.FOR Oct 16, 1996 9:03:13 AM
C**** Argon 100 - 1000 eV
C**** X= Te
C**** Y= Lz
C**** Eqn# 7104  lny=(a+cx+ex^2)/(1+bx+dx^2+fx^3)
C**** r2=0.9996946261752012
C**** r2adj=0.9996909171813777
C**** StdErr=3.274435515970801E-21
C**** Fval=324093.8153640019
C**** a= -43.72013282403189
C**** b= -0.003635592110764488
C**** c= 0.1665181119513338
C**** d= 7.084158574514305E-06
C**** e= -0.0003261764129687459
C**** f= 7.810575603077244E-11
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-43.72013282403189+X*(0.1665181119513338+X*(-
     .    0.0003261764129687459))) / (1.0+X*(-0.003635592110764488+X*(
     .    7.084158574514305E-06+X*(7.810575603077244E-11))))
      Y = DEXP(Y)
      ARGON2 = Y

      RETURN
      END
 
c//////////////////// FUNCTION ARGON3 //////////////////////////////////

      DOUBLE PRECISION FUNCTION ARGON3 (X)

C**** Argon 1000 to 100000 eV
C**** X= Te
C**** Y= Lz
C**** Eqn# 6206  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4+jx^5+k/x^5
C**** r2=0.9997293838807432
C**** r2adj=0.9997232964015779
C**** StdErr=5.527830175323746E-23
C**** Fval=181019.2975373837
C**** a= 3.150351699141569E-21
C**** b= 8.913074800293119E-25
C**** c= 3.129874709400417E-17
C**** d= -3.00095716911636E-29
C**** e= 4.846429512967016E-13
C**** f= 5.252628921061997E-34
C**** g= -1.432308868861869E-09
C**** h= -4.308069609855959E-39
C**** i= 1.414931981625638E-06
C**** j= 1.338048174038565E-44
C**** k= -0.0004718581892746704
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (8.913074800293119E-25+X*(-3.000957169116360E-29+X*(
     .    5.252628921061997E-34+X*(-4.308069609855959E-39+X*(
     .    1.338048174038565E-44))))) + 3.150351699141569E-21 + N * (
     .    3.129874709400417E-17+N*(4.846429512967016E-13+N*(-
     .    1.432308868861869E-09+N*(1.414931981625638E-06+N*(-
     .    0.0004718581892746704)))))
      ARGON3 = Y

      RETURN
      END
 
c//////////////////// FUNCTION KRYPTON1 ////////////////////////////////

      DOUBLE PRECISION FUNCTION KRYPTON1 (X)


c**** Kr (50 -1000 eV)
c**** X= Te
c**** Y= Lz
c**** Eqn# 6007  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8+jx^9+kx^(10)
c**** a= -4.247627744434805E-20
c**** b= -3.96855504681203E-21
c**** c=  2.802242360320274E-22
c**** d= -2.019876543038057E-24
c**** e=  7.743995823715699E-27
c**** f= -1.911563758126309E-29
c**** g=  3.151930534530092E-32
c**** h= -3.430943976743534E-35
c**** i=  2.358777152872231E-38
c**** j= -9.26504248442168E-42
c**** k=  1.583783095228654E-45

      Double Precision X, Y
      Y = -4.247627744434805E-20 + X * (-3.968555046812030E-21+X*(
     .    2.802242360320274E-22+X*(-2.019876543038057E-24+X*(
     .    7.743995823715699E-27+X*(-1.911563758126309E-29+X*(
     .    3.151930534530092E-32+X*(-3.430943976743534E-35+X*(
     .    2.358777152872231E-38+X*(-9.265042484421680E-42+X*(
     .    1.583783095228654E-45))))))))))
     
      krypton1 = y
      
      Return
      End
 
c//////////////////// FUNCTION KRYPTON2 ////////////////////////////////

      DOUBLE PRECISION FUNCTION KRYPTON2 (X)
      
c**** Kr (1000 - 5000 eV)
c**** X= Te
c**** Y= Lz
c**** Eqn# 6124  y^(0.5)=a+bx+cx^2+dx^3+ex^4+fx^5
c**** a=  2.527410904977941E-09
c**** b= -2.910755002646449E-12
c**** c=  2.006288706198077E-15
c**** d= -6.917713545445094E-19
c**** e=  1.13005805441513E-22
c**** f= -7.017717884910174E-27
      
      Double Precision X, Y
      Y = 2.527410904977941E-09 + X * (-2.910755002646449E-12+X*(
     .    2.006288706198077E-15+X*(-6.917713545445094E-19+X*(
     .    1.130058054415130E-22+X*(-7.017717884910174E-27)))))
      Y = Y * Y
      krypton2 = y
      
      Return
      End
 
c//////////////////// FUNCTION KRYPTON3 ////////////////////////////////

      DOUBLE PRECISION FUNCTION KRYPTON3 (X)
      
c**** Kr (5.0 to 15.0 keV)
c**** X=Te
c**** Y=Lz
c**** Eqn# 6113  y^(-1)=a+bx+cx^2+dx^3+ex^4+fx^5
c**** a= -1.07321284691047E+19
c**** b=  5087138924288157
c**** c= -219107668186.0897
c**** d= -26042499.97607108
c**** e=  2630.950795171409
c**** f= -0.06568445378271925
      
      Double Precision X, Y
      Y = -1.073212846910470E+19 + X * (5087138924288157.+X*(-
     .    219107668186.0897+X*(-26042499.97607108+X*(2630.950795171409+X
     .    *(-0.06568445378271925)))))
      Y = 1 / Y
      krypton3 = y
      
      Return
      End

c///////CARBON DATA ////////////////////////////////////////////////////

c/    Range from 3 to 10 eV:

      DOUBLE PRECISION FUNCTION CARBON0(X)
C----------------------------------------------------------*
C**** Carbon (3-10 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 1552  y^(-1)=a+be^x+ce^(-x)
C**** r2=0.9886010612518466
C**** r2adj=0.9865894838257018
C**** StdErr=3.491976853553102E-20
C**** Fval=780.547184948068
C**** a= 7.771844618017265E+17
C**** b= 59122888943831.57
C**** c= 2.426539413730379E+20
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION X1, X2
      X1 = DEXP(X)
      X2 = DEXP(-X)
      Y = 7.771844618017265E+17 + 59122888943831.57 * X1 + 
     .    2.426539413730379E+20 * X2
      Y = 1 / Y
      CARBON0 = Y
      RETURN
      END


      DOUBLE PRECISION FUNCTION CARBON1(X)

C**** Range 10 - 100 eV
C**** X= Te
C**** Y= Lz
C**** Eqn# 6206  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4+jx^5+k/x^5
C**** r2=0.9999952235181305
C**** r2adj=0.9999945584383765
C**** StdErr=1.56766897330871E-22
C**** Fval=1674864.891501524
C**** a= 1.639722967154732E-18
C**** b= -2.654182242083568E-20
C**** c= -6.37423721093842E-17
C**** d= 2.458509606418251E-22
C**** e= 1.590289149116386E-15
C**** f= -9.488347829679447E-25
C**** g= -2.435345151607877E-14
C**** h= -1.316645200326443E-27
C**** i= 2.038708041891294E-13
C**** j= 1.459022440553818E-29
C**** k= -6.468093918696965E-13
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (-2.654182242083568E-20+X*(2.458509606418251E-22+X*(
     .    -9.488347829679447E-25+X*(-1.316645200326443E-27+X*(
     .    1.459022440553818E-29))))) + 1.639722967154732E-18 + N * (
     .    -6.374237210938420E-17+N*(1.590289149116386E-15+N*(-
     .    2.435345151607877E-14+N*(2.038708041891294E-13+N*(-
     .    6.468093918696965E-13)))))
      CARBON1 = Y
      RETURN
      END

      DOUBLE PRECISION FUNCTION CARBON2(X)

C**** Range: 100 - 10000 eV

C**** X= Te
C**** Y= Lz
C**** Eqn# 6206  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4+jx^5+k/x^5
C**** r2=0.9998905133338273
C**** r2adj=0.9998880504435862
C**** StdErr=7.031561417044288E-24
C**** Fval=447494.0818463843
C**** a= -5.770913821900679E-22
C**** b= 5.333211787469455E-25
C**** c= 1.042550985018631E-18
C**** d= -1.404553291126376E-28
C**** e= -3.134926132115641E-16
C**** f= 2.164576220631324E-32
C**** g= 8.18472272168479E-14
C**** h= -1.691926173232754E-36
C**** i= -6.157523360890277E-12
C**** j= 5.202445752080436E-41
C**** k= 9.84843322304437E-11
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (5.333211787469455E-25+X*(-1.404553291126376E-28+X*(
     .    2.164576220631324E-32+X*(-1.691926173232754E-36+X*(
     .    5.202445752080436E-41)))))  -5.770913821900679E-22 + N * (
     .    1.042550985018631E-18+N*(-3.134926132115641E-16+N*(
     .    8.184722721684790E-14+N*(-6.157523360890277E-12+N*(
     .    9.848433223044370E-11)))))
      CARBON2 = Y

      RETURN
      END

c///////OXYGEN DATA ////////////////////////////////////////////////////

      DOUBLE PRECISION FUNCTION OXY1(X)
C----------------------------------------------------------*
C**** Oxygen (5-25 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 1003  y=a+bx+cx^2
C**** r2=0.999678973099738
C**** r2adj=0.9996690444327195
C**** StdErr=5.755430985585936E-21
C**** Fval=152586.1840297616
C**** a= -2.526944504853684E-19
C**** b= 7.032493188960815E-20
C**** c= -5.380401893701647E-22
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = -2.526944504853684E-19 + X * (7.032493188960815E-20+X*(-
     .    5.380401893701647E-22))
      OXY1 = Y
      RETURN
      END

      DOUBLE PRECISION FUNCTION OXY2(X)
C----------------------------------------------------------*
C**** Oxygen (25 - 100 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 1417  lny=a+b/x^(1.5)+ce^(-x)
C**** r2=0.952931476897593
C**** r2adj=0.9514757493789618
C**** StdErr=6.17708858998423E-20
C**** Fval=992.0354260190122
C**** a= -47.79687557784926
C**** b= 951.2480473186996
C**** c= -86758921098.56986
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION X1, X2
      X1 = 1.0 / (X*DSQRT(X))
      X2 = DEXP(-X)
      Y = -47.79687557784926 + 951.2480473186996 * X1 - 
     .    86758921098.56986 * X2
      Y = DEXP(Y)
      OXY2 = Y
      RETURN
      END

      DOUBLE PRECISION FUNCTION OXY3(X)
C----------------------------------------------------------*
C**** Oxygen (100 - 3000 eV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 1610  y^(-1)=a+b/x^(0.5)+clnx/x
C**** r2=0.9984397985282939
C**** r2adj=0.9984303808131729
C**** StdErr=1.105846129663765E-22
C**** Fval=159345.7731850993
C**** a= 1.569330703918676E+21
C**** b= -5.720727773028839E+22
C**** c= 9.620830252321734E+22
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION X1, X2
      X1 = 1.0 / DSQRT(X)
      X2 = DLOG(X) / X
      Y = 1.569330703918676E+21 - 5.720727773028839E+22 * X1 + 
     .    9.620830252321734E+22 * X2
      Y = 1 / Y
      OXY3 = Y
      RETURN
      END

      DOUBLE PRECISION FUNCTION OXY4(X)
C----------------------------------------------------------*
C**** Oxygen (3 to 10 keV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 6205  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4+jx^5
C**** r2=0.9999014286675548
C**** r2adj=0.9998994170077089
C**** StdErr=5.941108447011281E-25
C**** Fval=553408.14198754
C**** a= -2.975292054062271E-18
C**** b= 5.361343468740014E-22
C**** c= 1.072449986847168E-14
C**** d= -6.27069278557776E-26
C**** e= -2.418587102472807E-11
C**** f= 4.59799217188739E-30
C**** g= 3.100397720978912E-08
C**** h= -1.920293070180654E-34
C**** i= -1.723016446709003E-05
C**** j= 3.485046841440938E-39
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (5.361343468740014E-22+X*(-6.270692785577760E-26+X*(
     .    4.597992171887390E-30+X*(-1.920293070180654E-34+X*(
     .    3.485046841440938E-39))))) -2.975292054062271E-18 + N * (
     .    1.072449986847168E-14+N*(-2.418587102472807E-11+N*(
     .    3.100397720978912E-08+N*(-1.723016446709003E-05))))
C     
      OXY4 = Y
C
      RETURN
      END

      DOUBLE PRECISION FUNCTION OXY5(X)
C----------------------------------------------------------*
C**** Oxygen (10 - 50 keV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 4104  y=a+bx+cxlnx+dx^(0.5)lnx+e/x^2
C**** r2=0.9999999931488669
C**** r2adj=0.9999999930796635
C**** StdErr=3.043945500636726E-26
C**** Fval=18099195766.41754
C**** a= 3.830559189516505E-22
C**** b= -7.620779923225267E-26
C**** c= 4.802543659133072E-27
C**** d= 1.47015764315235E-24
C**** e= 1.63358460430098E-15
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION X1, X2, X3, X4
      X1 = X
      X2 = X * DLOG(X)
      X3 = DSQRT(X) * DLOG(X)
      X4 = 1.0 / (X*X)
      Y = 3.830559189516505E-22 - 7.620779923225267E-26 * X1 + 
     .    4.802543659133072E-27 * X2 + 1.470157643152350E-24 * X3 + 
     .    1.633584604300980E-15 * X4
      OXY5 = Y
      RETURN
      END

c/////////// Fe fits //////////////////////////////////////////////////

c/    Range 1: 10 to 1000 eV:

      DOUBLE PRECISION FUNCTION IRON1(X)
C----------------------------------------------------------*
C**** Iron (10 eV to 1 keV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 6308  y=a+blnx+c(lnx)^2+d(lnx)^3+e(lnx)^4+f(lnx)^5+
C**** g(lnx)^6+h(lnx)^7+i(lnx)^8+j(lnx)^9+k(lnx)^(10)
C**** r2=0.9995914402253218
C**** r2adj=0.9995822497191429
C**** StdErr=1.60011099702924E-20
C**** Fval=119884.4908547695
C**** a= -3.767068451905311E-16
C**** b= 1.167049534140273E-15
C**** c= -1.564027819892049E-15
C**** d= 1.201070742337944E-15
C**** e= -5.863672972354394E-16
C**** f= 1.900981581024757E-16
C**** g= -4.140487571695907E-17
C**** h= 5.981149466428279E-18
C**** i= -5.488822966452407E-19
C**** j= 2.894537848560817E-20
C**** k= -6.675935765344652E-22
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      X = DLOG(X)
      Y = -3.767068451905311E-16 + X * (1.167049534140273E-15+X*(-
     .    1.564027819892049E-15+X*(1.201070742337944E-15+X*(-
     .    5.863672972354394E-16+X*(1.900981581024757E-16+X*(-
     .    4.140487571695907E-17+X*(5.981149466428279E-18+X*(-
     .    5.488822966452407E-19+X*(2.894537848560817E-20+X*(-
     .    6.675935765344652E-22))))))))))
      IRON1 = Y
      RETURN
      END

c/    Range 2: 1 to 10 keV:

      DOUBLE PRECISION FUNCTION IRON2(X)
C----------------------------------------------------------*
C**** Iron (1 - 10 keV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 7104  lny=(a+cx+ex^2)/(1+bx+dx^2+fx^3)
C**** r2=0.9997039263428143
C**** r2adj=0.9997003303064922
C**** StdErr=1.585382369505486E-21
C**** Fval=334277.2526562476
C**** a= -42.19483219321717
C**** b= -0.0007557892111234372
C**** c= 0.0324452894801835
C**** d= 3.080297295076421E-07
C**** e= -1.346534855784905E-05
C**** f= -2.046361594047959E-13
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-42.19483219321717+X*(0.03244528948018350+X*(-
     .    1.346534855784905E-05))) / (1.0+X*(-0.0007557892111234372+X*(
     .    3.080297295076421E-07+X*(-2.046361594047959E-13))))
      Y = DEXP(Y)
      IRON2 = Y
      RETURN
      END

c/    Range 3: 10 to 100 keV:

      DOUBLE PRECISION FUNCTION IRON3(X)
C----------------------------------------------------------*
C**** Iron (10 to 100 keV)
C**** X= Te
C**** Y= Lz
C**** Eqn# 6205  y=a+bx+c/x+dx^2+e/x^2+fx^3+g/x^3+hx^4+i/x^4+jx^5
C**** r2=0.9989779697808669
C**** r2adj=0.9989571120212927
C**** StdErr=1.068723015716669E-22
C**** Fval=53325.03590293295
C**** a= 1.008718401006994E-18
C**** b= -2.568991921587864E-23
C**** c= -2.120110890240461E-14
C**** d= 3.922740389077567E-28
C**** e= 2.745022743025846E-10
C**** f= -3.293025031887385E-33
C**** g= -1.87650535765289E-06
C**** h= 1.37444717239993E-38
C**** i= 0.005231668669431576
C**** j= -2.082388289650421E-44
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      DOUBLE PRECISION N
      N = 1.0 / X
      Y = X * (-2.568991921587864E-23+X*(3.922740389077567E-28+X*(
     .    -3.293025031887385E-33+X*(1.374447172399930E-38+X*(-
     .    2.082388289650421E-44))))) + 1.008718401006994E-18 + N * (
     .    -2.120110890240461E-14+N*(2.745022743025846E-10+N*(-
     .    1.876505357652890E-06+N*(0.005231668669431576))))
      IRON3 = Y
      RETURN
      END
