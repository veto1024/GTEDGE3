      subroutine ncefits (iz, elecTemp, tau_res, Lz_eff, zav)

c///////////////////////////////////////////////////////////////////////
c/                                                                     /
c/    This subroutine calculates the effective cooling rate (Lz) and   /
c/    the average atomic number <Z> for selected impurities, under     /
c/    non-coronal equilibrium assumptions, taking into account the     /
c/    effect of finite residence time.                                 /
c/                                                                     /
c/    The fits have been produced with the TableCurve 2D package using /
c/    data from numerical simulations.                                 /
c/                                                                     /
c/    Written by John Mandrekas, GIT, 07/31/97                         /
c/                                                                     /
c/    Input and output variables:                                      /
c/    --------------------------                                       /
c/    iz       : atomic number of impurity species                     /
c/    elecTemp : electron temperature (eV)                             /
c/    tau_res  : impurity residence time (s)                           /
c/    Lz_eff   : effective cooling rate (erg-cm^3/s)                   /
c/    Zav      : local value of <Z>                                    /
c/                                                                     /
c///////////////////////////////////////////////////////////////////////

      implicit none
      logical first
      integer iz, it_ar
      parameter (it_ar = 6)
      real elecTemp, tau_res, tau, Lz_eff, Zav, teev, gf, ltau
      real tau_ar(it_ar), ltauar(it_ar)
      double precision logte, loglz, loglz1, loglz2, avz, avz1, avz2, 
     .    dtau

      integer i, itp, itm

      data tau_ar /1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 1.0/
      data first /.true./
      save first
      
      if (first) then
	 do i = 1, it_ar
	    ltauar(i) = alog(tau_ar(i))
	 enddo
	 first = .false.
      endif

c/    Impose limits in electron temperature and tau:
c/    The fits are valid in the range 10 eV < T_e < 1 keV and
c/    for tau < 1 s. If the input variables are outside this
c/    range, they are forced to be inside:
c/    For Kr the low temperature level is 30 eV and the smallest
c/    tau is 1.0e-6
     
c/    Define local values of temperature and tau:
c/    ------------------------------------------
      teev = elecTemp
      tau  = tau_res

      if (teev.GT.1000.0) teev = 1000.0

      if ((iz.EQ.10).OR.(iz.EQ.18)) then
         if (tau.GT.1.0) tau = 1.0
         if (teev.LT.10.0) teev = 10.0
         if (tau.LT.1.0e-4) tau = 1.0e-4
      else if (iz.EQ.36) then
         if (tau.GT.0.1) tau = 0.1
         if (tau.LT.1.0e-6) tau = 1.0e-6
         if (teev.LT.30.0) teev = 30.0
      endif
      
c/    Neon fits (3D fits from Table Curve 3D)
c/    ---------------------------------------
      if (iz.EQ.10) then
	 dtau = dble(alog(tau))
	 logte = dble(alog(teev))
	 call neLzeff(logte, dtau, loglz)
	 call nezavg (logte, dtau, avz)
	 Zav = avz
	 Lz_eff = DEXP(loglz)

c/    Argon fits (logarithmic interpolation in tau):
c/    ---------------------------------------------
      else if (iz.EQ.18) then
	 ltau = alog(tau)
	 logte = dble(alog(teev))
	 do i = 1, it_ar
	   itp = i
	   itm = i-1
	   if (ltauar(i).GT.ltau) go to 10
	 enddo

 10      gf = (ltau-ltauar(itm)) / (ltauar(itp) - ltauar(itm)) 

	 if (itp.EQ.2) then
	    call argon4Lz(logte, logLz2)
	    call argon4z(logte, avz2)
	    call argon5Lz(logte, logLz1)
	    call argon5z(logte, avz1)
	 else if (itp.EQ.3) then
	    call argon3Lz(logte, logLz2)
	    call argon3z(logte, avz2)
	    call argon4Lz(logte, logLz1)
	    call argon4z(logte, avz1)
	 else if (itp.EQ.4) then
	    call argon2Lz(logte, logLz2)
	    call argon2z(logte, avz2)
	    call argon3Lz(logte, logLz1)
	    call argon3z(logte, avz1)
	 else if (itp.EQ.5) then
	    call argon1Lz(logte, logLz2)
	    call argon1z(logte, avz2)
	    call argon2Lz(logte, logLz1)
	    call argon2z(logte, avz1)
	 else if (itp.EQ.6) then
	    call argon0Lz(logte, logLz2)
	    call argon0z(logte, avz2)
	    call argon1Lz(logte, logLz1)
	    call argon1z(logte, avz1)
         endif

c/    Interpolate:
c/    ------------
	 loglz = logLz1 + gf * (logLz2 - logLz1)
	 avz = avz1 + gf * (avz2 - avz1)

	 Zav = avz
	 Lz_eff = DEXP(loglz)

c/    Kr 3D fit from TableCurve:
c/    -------------------------
      else if (iz.EQ.36) then
	 dtau = dble(alog(tau))
	 logte = dble(alog(teev))
	 call krLzeff(logte, dtau, loglz)
	 call krzavg (logte, dtau, avz)
	 Zav = avz
	 Lz_eff = DEXP(loglz)

      else 
	 write (6, 1000) iz
	 return
      endif

 1000 format (1x, 'Sorry! No data for element with Z = ', i3)
      return
      end

c/    Fitting routines for Ne:
*----------------------------------------------------------*
      SUBROUTINE EVALCSI(order, x, y, c, z)
*----------------------------------------------------------*
      INTEGER order
      DOUBLE PRECISION x,y,c(*),z
      INTEGER nc,j,m,iv
      DOUBLE PRECISION cx(12),cy(12),v(70)
      x=(x-(2.302600000000000D0))/(1.465880687853593D0)
      y=(y-(-9.210300000000000D0))/(2.931729544718568D0)
      IF(x.LT.0.0) x=0.0
      IF(x.GT.3.14159265358979323846) x=3.14159265358979323846
      IF(y.LT.0.0) y=0.0
      IF(y.GT.3.14159265358979323846) y=3.14159265358979323846

      if (order.EQ.5) then
        nc=2
      else if (order.EQ.9) then
        nc=3
      else if (order.EQ.14) then
        nc=4
      else if (order.EQ.20) then
        nc=5
      else if (order.EQ.27) then
        nc=6
      else if (order.EQ.35) then
        nc=7
      else if (order.EQ.44) then
        nc=8
      else if (order.EQ.54) then
        nc=9
      else if (order.EQ.65) then
        nc=10
      else
        z=0.0
        RETURN
      endif
      cx(1)=1.D0
      cy(1)=1.D0
      DO 10 j=1,nc
        cx(j+1)=DCOS(j*x)
        cy(j+1)=DCOS(j*y)
   10 CONTINUE
      iv=1
      DO 30 j=1,nc+1
        DO 20 m=j,1,-1
          v(iv)=cx(m)*cy(j-m+1)
          iv=iv+1
   20   CONTINUE
   30 CONTINUE
      z=0.0
      DO 40 j=1,order+1
        z = z + c(j)*v(j)
   40 CONTINUE
      RETURN
      END

*----------------------------------------------------------*
      SUBROUTINE neLzeff(x0,y0,z)
*----------------------------------------------------------*
* TableCurve 3D
* File Source= g:\imprty\corona\lzeff3d.txt
* Date= Jun 17, 1997
* Time= 2:52:39 PM
* Data Source= Neon
* X0 = Te
* Y0 = tau
* Z= Lzeff
* Eqn#= 534
* Eqn= Cosine Series Bivariate Order 5
* r2=0.993072122950754
* r2adj=0.9924582604274032
* StdErr=0.1318781610323391
* Fstat=1705.797920360059
* a= -43.23285956703503
* b= 1.32086416610929
* c= 1.123574264933135
* d= -0.252308955763543
* e= -1.166658676771189
* f= 0.2664309465339077
* g= -0.6588326411000396
* h= -0.3336800390752378
* i= -0.2616573468518904
* j= 0.009165909790695556
* k= -0.3254341930855056
* l= 0.2651892300424894
* m= -0.04982093654780001
* n= -0.08582600961212271
* o= 0.03785096451610444
* p= 0.04979145490959691
* q= 0.2729579802318522
* r= 0.1245359892100098
* s= 0.04062715801615898
* t= -0.0633707599336731
* u= 0.07057534834177396
*----------------------------------------------------------*
      DOUBLE PRECISION x0, y0, x,y,z
      DOUBLE PRECISION c(20+1)
      DATA c(1)/-43.23285956703503D0/
      DATA c(2)/1.320864166109290D0/
      DATA c(3)/1.123574264933135D0/
      DATA c(4)/-0.2523089557635430D0/
      DATA c(5)/-1.166658676771189D0/
      DATA c(6)/0.2664309465339077D0/
      DATA c(7)/-0.6588326411000396D0/
      DATA c(8)/-0.3336800390752378D0/
      DATA c(9)/-0.2616573468518904D0/
      DATA c(10)/0.009165909790695556D0/
      DATA c(11)/-0.3254341930855056D0/
      DATA c(12)/0.2651892300424894D0/
      DATA c(13)/-0.04982093654780001D0/
      DATA c(14)/-0.08582600961212271D0/
      DATA c(15)/0.03785096451610444D0/
      DATA c(16)/0.04979145490959691D0/
      DATA c(17)/0.2729579802318522D0/
      DATA c(18)/0.1245359892100098D0/
      DATA c(19)/0.04062715801615898D0/
      DATA c(20)/-0.06337075993367310D0/
      DATA c(21)/0.07057534834177396D0/

      x = x0
      y = y0

      CALL EVALCSI(20,x,y,c,z)
      RETURN
      END

c//////////////////////////////////////////////////////////////////////

      SUBROUTINE EVALSIG(order, x, y, c, z)

c/    Modified by John Mandrekas, to replace SELECT CASE with
c/    an if-the-else loop

      INTEGER order
      DOUBLE PRECISION x,y,c(*),z
      INTEGER pcnt,j,m,iv
      DOUBLE PRECISION px(12),py(12),v(70),ctr
      x=(x-(4.605200000000000D0))/(2.302600000000000D0)
      y=(y-(-4.605150000000000D0))/(4.605150000000000D0)

      if (order.EQ.5) then
        pcnt=3
      else if (order.EQ.9) then
        pcnt=4
      else if (order.EQ.14) then
        pcnt=5
      else if (order.EQ.20) then
        pcnt=6
      else if (order.EQ.27) then
        pcnt=7
      else if (order.EQ.35) then
        pcnt=8
      else if (order.EQ.44) then
        pcnt=9
      else if (order.EQ.54) then
        pcnt=10
      else if (order.EQ.65) then
        pcnt=11
      else
        z=0.0
        RETURN
      endif
      
      px(1)=1.D0
      py(1)=1.D0
      px(2)=x
      py(2)=y

      DO 10 j=3,pcnt
        ctr=-1.D0+DFLOAT(j-2)*(2.D0/(DFLOAT(pcnt)-1.D0))
        px(j)=-1.D0+2.D0/(1.D0+DEXP(-(x-ctr)/0.12D0))
        py(j)=-1.D0+2.D0/(1.D0+DEXP(-(y-ctr)/0.12D0))
   10 CONTINUE

      iv=1
      DO 30 j=1,pcnt
        DO 20 m=j,1,-1
          v(iv)=px(m)*py(j-m+1)
          iv=iv+1
   20   CONTINUE
   30 CONTINUE
      z=0.0
      DO 40 j=1,order+1
        z = z+c(j)*v(j)
   40 CONTINUE
      RETURN
      END

c//////////////////////////////////////////////////////////////////////

      SUBROUTINE nezavg(x0,y0,z)
c----------------------------------------------------------*
c TableCurve 3D
c File Source= e:\workspace\cefits\zavg3d.txt
c Date= Aug 4, 1997
c Time= 7:53:08 PM
c Data Source= Ne <Z>
c X0 = Te
c Y0 = tau
c Z = Zavg
c Eqn#= 605
c Eqn= Sigmoid Series Bivariate Order 6
c r2=0.9978553738075728
c r2adj=0.997594288879799
c StdErr=0.09570631055367352
c Fstat=3980.743645286862
c a= 5.333348461724801
c b= 1.318117953145604
c c= 4.087593551215054
c d= 0.548281721811961
c e= 0.6448150513972668
c f= -0.427414399152828
c g= 1.06338057911186
c h= -0.4103597066378987
c i= 0.1729422922028607
c j= -0.3352609521093104
c k= -0.03251174076433288
c l= 0.5217084687870015
c m= -0.06514097975553723
c n= -0.06342133454217161
c o= -0.5657419369343312
c p= 0.01803337942511036
c q= -0.8137682878792351
c r= 0.4142674506964295
c s= 0.03237150538628928
c t= 0.2813270961922968
c u= -0.4714531051966177
c v= 0.1963800853271892
c aa= 0.2686586947467848
c ab= -0.2924975020343322
c ac= -0.07344565198430797
c ad= 0.006074729268764003
c ae= -0.1443344083173938
c af= -1.178480962396732
c----------------------------------------------------------*
      DOUBLE PRECISION x0, y0, x,y,z
      DOUBLE PRECISION c(27+1)
      DATA c(1)/5.333348461724801D0/
      DATA c(2)/1.318117953145604D0/
      DATA c(3)/4.087593551215054D0/
      DATA c(4)/0.5482817218119610D0/
      DATA c(5)/0.6448150513972668D0/
      DATA c(6)/-0.4274143991528280D0/
      DATA c(7)/1.063380579111860D0/
      DATA c(8)/-0.4103597066378987D0/
      DATA c(9)/0.1729422922028607D0/
      DATA c(10)/-0.3352609521093104D0/
      DATA c(11)/-0.03251174076433288D0/
      DATA c(12)/0.5217084687870015D0/
      DATA c(13)/-0.06514097975553723D0/
      DATA c(14)/-0.06342133454217161D0/
      DATA c(15)/-0.5657419369343312D0/
      DATA c(16)/0.01803337942511036D0/
      DATA c(17)/-0.8137682878792351D0/
      DATA c(18)/0.4142674506964295D0/
      DATA c(19)/0.03237150538628928D0/
      DATA c(20)/0.2813270961922968D0/
      DATA c(21)/-0.4714531051966177D0/
      DATA c(22)/0.1963800853271892D0/
      DATA c(23)/0.2686586947467848D0/
      DATA c(24)/-0.2924975020343322D0/
      DATA c(25)/-0.07344565198430797D0/
      DATA c(26)/0.006074729268764003D0/
      DATA c(27)/-0.1443344083173938D0/
      DATA c(28)/-1.178480962396732D0/

      x = x0
      y = y0

      CALL EVALSIG(27,x,y,c,z)

      RETURN
      END

c///////////////////////////////////////////////////////////////////////
c/    Fitting routines for Ar:

C----------------------------------------------------------*
      SUBROUTINE ARGON0LZ(X,Y)
C----------------------------------------------------------*
C**** Aug 4, 1997 12:25:33 PM
C**** Argon CE
C**** X= log_TeV
C**** Y= log_Lzeff
C**** Eqn# 7007  y=(a+cx+ex^2+gx^3+ix^4)/(1+bx+dx^2+fx^3+hx^4)
C**** r2=0.99900605356637
C**** r2adj=0.9986747380884934
C**** StdErr=0.04265572961102707
C**** Fval=3517.816523283592
C**** a= -42.38250194528416
C**** b= -0.8588803576324078
C**** c= 36.60854939055573
C**** d= 0.2734680693001029
C**** e= -11.72334296489076
C**** f= -0.03813083157368667
C**** g= 1.644005554139061
C**** h= 0.001970262473937455
C**** i= -0.08542432585174578
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-42.38250194528416+X*(36.60854939055573+X*(-11.72334296489076
     .    +X*(1.644005554139061+X*(-0.08542432585174578))))) / (1.0+X*(
     .    -0.8588803576324078+X*(0.2734680693001029+X*(-
     .    0.03813083157368667+X*(0.001970262473937455)))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON0Z(X,Y)
C----------------------------------------------------------*
C**** Argon CE
C**** X= log_TeV
C**** Y= Zav
C**** Eqn# 7009  y=(a+cx+ex^2+gx^3+ix^4+kx^5)/(1+bx+dx^2+fx^3+hx^4+jx^5)
C**** r2=0.9998252992074379
C**** r2adj=0.9997484308587105
C**** StdErr=0.06289092426617098
C**** Fval=14879.98846378972
C**** a= -0.3511837450506403
C**** b= -0.9088390355186895
C**** c= 1.888941118900576
C**** d= 0.3125041530929955
C**** e= -1.471366233966557
C**** f= -0.049693230186344
C**** g= 0.4529369689190328
C**** h= 0.0035032498024078
C**** i= -0.06196150942351487
C**** j= -7.699023882919249E-05
C**** k= 0.003154312559101326
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-0.3511837450506403+X*(1.888941118900576+X*(-
     .    1.471366233966557+X*(0.4529369689190328+X*(-
     .    0.06196150942351487+X*(0.003154312559101326)))))) / (1.0+X*(
     .    -0.9088390355186895+X*(0.3125041530929955+X*(-
     .    0.04969323018634400+X*(0.003503249802407800+X*(-
     .    7.699023882919249E-05))))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON1LZ(X,Y)
C----------------------------------------------------------*
C**** Argon tau = 1e-2 s
C**** X= log_TeV
C**** Y= log_Lzeff
C**** Eqn# 7007  y=(a+cx+ex^2+gx^3+ix^4)/(1+bx+dx^2+fx^3+hx^4)
C**** r2=0.9992528521065143
C**** r2adj=0.9990038028086857
C**** StdErr=0.02538639788363597
C**** Fval=4680.980851135572
C**** a= -42.33260224498212
C**** b= -0.8583198487312616
C**** c= 36.52666797761916
C**** d= 0.2745909921635287
C**** e= -11.74455835358558
C**** f= -0.03867360423515321
C**** g= 1.661884533631507
C**** h= 0.002029849593447703
C**** i= -0.08759429431375242
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-42.33260224498212+X*(36.52666797761916+X*(-11.74455835358558
     .    +X*(1.661884533631507+X*(-0.08759429431375242))))) / (1.0+X*(
     .    -0.8583198487312616+X*(0.2745909921635287+X*(-
     .    0.03867360423515321+X*(0.002029849593447703)))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON1Z(X,Y)
C----------------------------------------------------------*
C**** Argon tau = 1.e-2
C**** X= log_TeV
C**** Y= Zav
C**** Eqn# 7009  y=(a+cx+ex^2+gx^3+ix^4+kx^5)/(1+bx+dx^2+fx^3+hx^4+jx^5)
C**** r2=0.9999279134449143
C**** r2adj=0.9998961953606766
C**** StdErr=0.03842323473303066
C**** Fval=36065.15211979798
C**** a= -0.139509793280015
C**** b= -0.8964192004588974
C**** c= 1.518414650535741
C**** d= 0.3068254345643882
C**** e= -1.211887510222264
C**** f= -0.04929388599198845
C**** g= 0.369377813535042
C**** h= 0.003622011670727258
C**** i= -0.04957705277144021
C**** j= -9.107562539881708E-05
C**** k= 0.002469536628254192
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-0.1395097932800150+X*(1.518414650535741+X*(-
     .    1.211887510222264+X*(0.3693778135350420+X*(-
     .    0.04957705277144021+X*(0.002469536628254192)))))) / (1.0+X*(
     .    -0.8964192004588974+X*(0.3068254345643882+X*(-
     .    0.04929388599198845+X*(0.003622011670727258+X*(-
     .    9.107562539881708E-05))))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON2LZ(X0,Y)
C----------------------------------------------------------*
C**** Argon tau = 3e-3
C**** X= log_TeV
C**** Y= log_Lzeff
C**** Eqn# 7208  y=(a+clnx+e(lnx)^2+g(lnx)^3+i(lnx)^4)/(1+blnx+d(lnx)^2+f(lnx)^3+h(lnx)^4+j(lnx)^5)
C**** r2=0.999533888556116
C**** r2adj=0.9993546149238529
C**** StdErr=0.01525037507763896
C**** Fval=6433.229016395087
C**** a= -41.84620861544481
C**** b= -2.617823991463277
C**** c= 109.5846694632419
C**** d= 2.579903197139291
C**** e= -107.7978302162973
C**** f= -1.135889764645552
C**** g= 47.13388761889036
C**** h= 0.1911134601625584
C**** i= -7.728012780968235
C**** j= -0.001384512956359468
C----------------------------------------------------------*
      DOUBLE PRECISION X0, X, Y
      X = DLOG(X0)
      Y = (-41.84620861544481+X*(109.5846694632419+X*(-107.7978302162973
     .    +X*(47.13388761889036+X*(-7.728012780968235))))) / (1.0+X*(
     .    -2.617823991463277+X*(2.579903197139291+X*(-1.135889764645552+
     .    X*(0.1911134601625584+X*(-0.001384512956359468))))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON2Z(X0,Y)
C----------------------------------------------------------*
C**** Argon tau = 3e-3
C**** X= log_TeV
C**** Y= Zav
C**** Eqn# 7206  y=(a+clnx+e(lnx)^2+g(lnx)^3)/(1+blnx+d(lnx)^2+
C****              f(lnx)^3+h(lnx)^4)
C**** r2=0.999900498938209
C**** r2adj=0.9998720700634116
C**** StdErr=0.0393724493636286
C**** Fval=41632.16803529677
C**** a= 1.640031769795291
C**** b= -2.372829972278729
C**** c= -2.776559297187706
C**** d= 2.107061481731835
C**** e= 1.510100187086858
C**** f= -0.8293449565278169
C**** g= -0.2520125972891318
C**** h= 0.1224907282424437
C----------------------------------------------------------*
      DOUBLE PRECISION X, X0, Y
      X = DLOG(X0)
      Y = (1.640031769795291+X*(-2.776559297187706+X*(1.510100187086858+
     .    X*(-0.2520125972891318)))) / (1.0+X*(-2.372829972278729+X*(
     .    2.107061481731835+X*(-0.8293449565278169+X*(0.1224907282424437
     .    )))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON3LZ(X,Y)
C----------------------------------------------------------*
C**** Argon tau=1e-3
C**** X= log_TeV
C**** Y= log_Lzeff
C**** Eqn# 7007  y=(a+cx+ex^2+gx^3+ix^4)/(1+bx+dx^2+fx^3+hx^4)
C**** r2=0.999562142674707
C**** r2adj=0.9994161902329427
C**** StdErr=0.01066071253323152
C**** Fval=7989.971383990876
C**** a= -42.26898903433604
C**** b= -0.8363189334942396
C**** c= 35.53500081562899
C**** d= 0.2672865221399071
C**** e= -11.40082756143058
C**** f= -0.03846349674702897
C**** g= 1.644548360564156
C**** h= 0.002113832612117444
C**** i= -0.09045393631686603
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-42.26898903433604+X*(35.53500081562899+X*(-11.40082756143058
     .    +X*(1.644548360564156+X*(-0.09045393631686603))))) / (1.0+X*(
     .    -0.8363189334942396+X*(0.2672865221399071+X*(-
     .    0.03846349674702897+X*(0.002113832612117444)))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON3Z(X,Y)
C----------------------------------------------------------*
C**** Argon tau=1e-3
C**** X= log_TeV
C**** Y= Zav
C**** Eqn# 7006  y=(a+cx+ex^2+gx^3)/(1+bx+dx^2+fx^3+hx^4)
C**** r2=0.9998915982126643
C**** r2adj=0.9998606262734255
C**** StdErr=0.03564051825558597
C**** Fval=38213.46632329514
C**** a= 0.5254247509502632
C**** b= -0.687790425938857
C**** c= 0.3096395702204568
C**** d= 0.1786323158483258
C**** e= -0.2141568673737633
C**** f= -0.02099817612625054
C**** g= 0.02815006542766565
C**** h= 0.0009994541898981667
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (0.5254247509502632+X*(0.3096395702204568+X*(-
     .    0.2141568673737633+X*(0.02815006542766565)))) / (1.0+X*(
     .    -0.6877904259388570+X*(0.1786323158483258+X*(-
     .    0.02099817612625054+X*(0.0009994541898981667)))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON4LZ(X,Y)
C----------------------------------------------------------*
C**** Argon tau=3e-4
C**** X= log_TeV
C**** Y= log_Lzeff
C**** Eqn# 7007  y=(a+cx+ex^2+gx^3+ix^4)/(1+bx+dx^2+fx^3+hx^4)
C**** r2=0.9996558748190148
C**** r2adj=0.999541166425353
C**** StdErr=0.006934252628265474
C**** Fval=10167.21749873011
C**** a= -42.31651002963976
C**** b= -0.790191093899583
C**** c= 33.66646051305002
C**** d= 0.2512030149717891
C**** e= -10.73444259008938
C**** f= -0.0375625862195634
C**** g= 1.604036134375151
C**** h= 0.002234533551695805
C**** i= -0.09507372296293638
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = (-42.31651002963976+X*(33.66646051305002+X*(-10.73444259008938
     .    +X*(1.604036134375151+X*(-0.09507372296293638))))) / (1.0+X*(
     .    -0.7901910938995830+X*(0.2512030149717891+X*(-
     .    0.03756258621956340+X*(0.002234533551695805)))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON4Z(X,Y)
C----------------------------------------------------------*
C**** Argon tau=3e-4
C**** X= log_TeV
C**** Y= Zav
C**** Eqn# 6004  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7
C**** r2=0.9998510061391379
C**** r2adj=0.9998084364646059
C**** StdErr=0.03304130536631604
C**** Fval=27801.41314956605
C**** a= -87.81697977367694
C**** b= 182.6506540142857
C**** c= -153.8510368873545
C**** d= 69.07805703363328
C**** e= -17.7019511277191
C**** f= 2.599437743769436
C**** g= -0.2034928673019496
C**** h= 0.006581547515335523
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = -87.81697977367694 + X * (182.6506540142857+X*(-
     .    153.8510368873545+X*(69.07805703363328+X*(-17.70195112771910+X
     .    *(2.599437743769436+X*(-0.2034928673019496+X*(
     .    0.006581547515335523)))))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON5LZ(X,Y)
C----------------------------------------------------------*
C**** Argon tau = 1e-4
C**** X= log_TeV
C**** Y= log_Lzeff
C**** Eqn# 6005  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8
C**** r2=0.9991174702032939
C**** r2adj=0.9988232936043919
C**** StdErr=0.008584561525650195
C**** Fval=3962.371762135777
C**** a= -212.1702327595376
C**** b= 340.4516109247146
C**** c= -289.5190403864003
C**** d= 136.8124046392919
C**** e= -39.14077575841633
C**** f= 6.938095515774705
C**** g= -0.7453056613901851
C**** h= 0.0444712334280873
C**** i= -0.001131573325904769
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = -212.1702327595376 + X * (340.4516109247146+X*(-
     .    289.5190403864003+X*(136.8124046392919+X*(-39.14077575841633+X
     .    *(6.938095515774705+X*(-0.7453056613901851+X*(
     .    0.04447123342808730+X*(-0.001131573325904769))))))))
      RETURN
      END
C
C----------------------------------------------------------*
      SUBROUTINE ARGON5Z(X,Y)
C----------------------------------------------------------*
C**** TableCurve D:\WORKSPC\NOCEFITS\ARGON\ARGON5Z.FOR Aug 4, 1997 12:59:35 PM
C**** Argon tau = 1e-4
C**** X= log_TeV
C**** Y= Zav
C**** Eqn# 6004  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7
C**** r2=0.9998768961068195
C**** r2adj=0.9998417235659107
C**** StdErr=0.02306965161613673
C**** Fval=33649.19690183208
C**** a= -11.66049767121149
C**** b= 37.70124410303969
C**** c= -38.99223950136468
C**** d= 19.91634626271103
C**** e= -5.470925393879433
C**** f= 0.832268921389279
C**** g= -0.06615414281417286
C**** h= 0.002145733931788694
C----------------------------------------------------------*
      DOUBLE PRECISION X, Y
      Y = -11.66049767121149 + X * (37.70124410303969+X*(-
     .    38.99223950136468+X*(19.91634626271103+X*(-5.470925393879433+X
     .    *(0.8322689213892790+X*(-0.06615414281417286+X*(
     .    0.002145733931788694)))))))
      RETURN
      END

c///////////////////////////////////////////////////////////////////////
c/    Fitting routines for Kr:

      SUBROUTINE EVALFSI(order, x, y, c, z)
*----------------------------------------------------------*
      INTEGER order
      DOUBLE PRECISION x,y,c(*),z
      INTEGER i,sicnt,j,k,m,n1,n2
      DOUBLE PRECISION vx(10),vy(10),v(70)
      x=(x-(3.401200000000000D0))/(1.116185446892080D0)
      y=(y-(-13.81600000000000D0))/(3.664829043588456D0)
      IF(x.LT.0.0) x=0.0
      IF(x.GT.3.14159265358979323846) x=3.14159265358979323846
      IF(y.LT.0.0) y=0.0
      IF(y.GT.3.14159265358979323846) y=3.14159265358979323846

      if (order.EQ.12) then
        sicnt=2
      else if (order.EQ.24) then
        sicnt=3
      else if (order.EQ.40) then
        sicnt=4
      else if (order.EQ.60) then
        sicnt=5
      else
        z=0.0
        RETURN
      endif

      j=1
      DO 10 i=1,sicnt
        vx(j)=DCOS(i*x)
        vx(j+1)=DSIN(i*x)
        vy(j)=DCOS(i*y)
        vy(j+1)=DSIN(i*y)
        j=j+2
   10 CONTINUE

      if (sicnt.EQ.2) then
        m=4
        n1=2
        n2=2
      else if (sicnt.EQ.3) then
        m=6
        n1=2
        n2=4
      else if (sicnt.EQ.4) then
        m=8
        n1=4
        n2=6
      else if (sicnt.EQ.5) then
        m=10
        n1=4
        n2=8
      endif

      v(1)=1.0
      k=2
      DO 20 i=1,m
        v(k)=vx(i)
        k=k+1
        v(k)=vy(i)
        k=k+1
   20 CONTINUE
      DO 40, i=1,n1
        IF(i.EQ.3) n2=n2-2
        DO 30 j=i,n2
          v(k)=vx(i)*vy(j)
          k=k+1
          IF(j.NE.i) THEN
            v(k)=vx(j)*vy(i)
            k=k+1
          END IF
   30   CONTINUE
   40 CONTINUE
      z=0.0
      DO 50 j=1,order+1
        z = z+ c(j)*v(j)
   50 CONTINUE
      RETURN
      END

      SUBROUTINE EVALCPOLY(order, logx, logy, x, y, c, z)
c----------------------------------------------------------*
      INTEGER order,logx,logy
      DOUBLE PRECISION x,y,c(*),z
      INTEGER tcnt,j,m,iv
      DOUBLE PRECISION tx(12),ty(12),v(70)
      IF(logx.NE.1) THEN
        x=(x-(5.154500000000000D0))/(1.753300000000000D0)
      ELSE
        x=(DLOG(x)-(1.578389759227518D0))/(0.3542614486980182D0)
      END IF
      IF(logy.NE.1) THEN
        y=(y-(-8.059300000000000D0))/(5.756700000000000D0)
      ELSE
        y=(DLOG(y)-(0.000000000000000D0))/(0.000000000000000D0)
      END IF

      if (order.EQ.5) then
        tcnt=3
      else if (order.EQ.9) then
        tcnt=4
      else if (order.EQ.14) then
        tcnt=5
      else if (order.EQ.20) then
        tcnt=6
      else if (order.EQ.27) then
        tcnt=7
      else if (order.EQ.35) then
        tcnt=8
      else if (order.EQ.44) then
        tcnt=9
      else if (order.EQ.54) then
        tcnt=10
      else if (order.EQ.65) then
        tcnt=11
      else
        z=0.0
        RETURN
      endif

      IF(tcnt.GT.6) THEN
        IF(x.LT.-1.D0) x=-1.D0
        IF(x.GT. 1.D0) x= 1.D0
        IF(y.LT.-1.D0) y=-1.D0
        IF(y.GT. 1.D0) y= 1.D0
      END IF
      tx(1)=1.D0
      ty(1)=1.D0
      tx(2)=x
      ty(2)=y
      DO 10 j=3,tcnt
        tx(j)=2*x*tx(j-1)-tx(j-2)
        ty(j)=2*y*ty(j-1)-ty(j-2)
   10 CONTINUE
      iv=1
      DO 30 j=1,tcnt
        DO 20 m=j,1,-1
          v(iv)=tx(m)*ty(j-m+1)
          iv=iv+1
   20   CONTINUE
   30 CONTINUE
      z=0.0
      DO 40, j=1,order+1
        z = z + c(j)*v(j)
   40 CONTINUE
      RETURN
      END

c----------------------------------------------------------*
      SUBROUTINE krLzeff(x0, y0, z)
c----------------------------------------------------------*
c TableCurve 3D
c File Source= e:\workspace\cefits\lzeff3d.txt
c Date= Aug 11, 1997
c Time= 9:44:48 PM
c Data Source= Kr Lz
c X= Te
c Y= tau
c Z= Lzeff
c Eqn#= 423
c Eqn= Chebyshev LnX,Y Bivariate Polynomial Order 4
c r2=0.9962594736486986
c r2adj=0.9958829105931984
c StdErr=0.07228458817985718
c Fstat=2853.664870592914
c a= -40.59384464075425
c b= 0.5969059166511411
c c= -1.626506396551109
c d= -0.1423725415256514
c e= 1.032148666200518
c f= 0.3177447018374038
c g= -0.1960408018253158
c h= -0.3701880930881759
c i= -0.02860883025889682
c j= 0.1196929740165592
c k= 0.01435607377767004
c l= -0.1594257613091802
c m= -0.2040822117091972
c n= -0.2220904134615838
c o= -0.04912957109645291
c----------------------------------------------------------*
      DOUBLE PRECISION x0,y0,x,y,z
      DOUBLE PRECISION c(14+1)
      DATA c(1)/-40.59384464075425D0/
      DATA c(2)/0.5969059166511411D0/
      DATA c(3)/-1.626506396551109D0/
      DATA c(4)/-0.1423725415256514D0/
      DATA c(5)/1.032148666200518D0/
      DATA c(6)/0.3177447018374038D0/
      DATA c(7)/-0.1960408018253158D0/
      DATA c(8)/-0.3701880930881759D0/
      DATA c(9)/-0.02860883025889682D0/
      DATA c(10)/0.1196929740165592D0/
      DATA c(11)/0.01435607377767004D0/
      DATA c(12)/-0.1594257613091802D0/
      DATA c(13)/-0.2040822117091972D0/
      DATA c(14)/-0.2220904134615838D0/
      DATA c(15)/-0.04912957109645291D0/

      x = x0
      y = y0

      CALL EVALCPOLY(14,1,0,x,y,c,z)
      RETURN
      END
c/    Fitting routines for Kr <Z>

c----------------------------------------------------------*
      SUBROUTINE krzavg(x0,y0,z)
c----------------------------------------------------------*
c TableCurve 3D
c Date= Aug 8, 1997
c Data Source= Kr <Z>
c X= Te
c Y= tau
c Z= Zavg
c Eqn#= 521
c Eqn= Fourier Series Bivariate Order 2x2
c r2=0.9987996415029612
c r2adj=0.9986962993806996
c StdErr=0.2126429365928602
c Fstat=10539.7363845172
c a= 8.605909451537082
c b= -4.046148428074152
c c= -7.221022865532085
c d= -1.274257108080557
c e= 3.52123472623367
c f= 0.2349160972383984
c g= 0.2323702058926043
c h= 0.4230616969117336
c i= 0.4516609631946268
c j= 3.658486287866877
c k= -0.7634369664971292
c l= 2.340301574911032
c m= 0.1991483939827037
c----------------------------------------------------------*
      DOUBLE PRECISION x0, y0, x,y,z
      DOUBLE PRECISION c(12+1)
      DATA c(1)/8.605909451537082D0/
      DATA c(2)/-4.046148428074152D0/
      DATA c(3)/-7.221022865532085D0/
      DATA c(4)/-1.274257108080557D0/
      DATA c(5)/3.521234726233670D0/
      DATA c(6)/0.2349160972383984D0/
      DATA c(7)/0.2323702058926043D0/
      DATA c(8)/0.4230616969117336D0/
      DATA c(9)/0.4516609631946268D0/
      DATA c(10)/3.658486287866877D0/
      DATA c(11)/-0.7634369664971292D0/
      DATA c(12)/2.340301574911032D0/
      DATA c(13)/0.1991483939827037D0/

      x = x0
      y = y0

      CALL EVALFSI(12,x,y,c,z)
      RETURN
      END
