integer nit,npsi0,nmesh,ioptionorb,ioptxtran,n0,nd,nn,n,nx,m,npsi				 
parameter (nit=9,npsi0=22,nmesh=25)
real Icur,Btor,radminor,rmajor,xm,eq,radwall,dr,rm,xmu0,zmass,
	1      delna,rminor0,rminord,thetax,deltathetx,elong,deltarx,
	2      rminorx,xpi,xx,xy,ersum,xk,xxi,fphid,fphi0,fphix,rxloss,
	3      tion,psi0,wtrap,theta0,thetad,h0,hd,delphi,xpot,xflux,a,b,
	4      c,abc,x,value,value1,fvalue,enbi,fvalue1,valuev,fvvalue,
	5      fvvalue1,evalue,evalue1,fevalue,fevalue1,avtfloss,avtvloss,
	6      avteloss,avffloss,avfvloss,avfeloss,tflosav,tvlosav,telosav,
	7      fflosav,fvlosav,felosav

	real      v0m(nit,nit),v0p(nit,nit),xke0m(nit,nit),xke0p(nit,nit),
	 1   	  angle0(nit),angleD(nit),tfloss(npsi0,nit),
     2		  tvloss(npsi0,nit),teloss(npsi0,nit),ffloss(npsi0,nit),
     3		  fvloss(npsi0,nit),feloss(npsi0,nit),epmin(npsi0,nit)
     4		  ephi(nmesh),xkeDm(nit,nit),xkeDp(nit,nit),xti(nmesh),
     5		  tflos(nit),tvlos(nit),telos(nit),fflos(nit),fvlos(nit),
     6          felos(nit),erav(nit),rhor(nmesh),emin(npsi0,nit),
     7          cos0(npsi0),erex(nmesh),ephi(nmesh)