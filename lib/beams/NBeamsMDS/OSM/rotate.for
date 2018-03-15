	Subroutine rotate(m)
c	calculates toroidal velocity integration at mp=m
	include 'soldiv.fi'
	parameter (jq=25)


c	form 1 both vphi's eval from toroidal mom bal
	yy1(m) = atnum(1)*eq*(ephia+bthet*vrad1(m))/xmas(1) +
	1		 xmtor(1)/(yni(m,1)*xmas(1))
	yy2(m) = zbar2(m)*eq*(ephia+bthet*vrad2(m))/xmas(2) +
	1		 xmtor(2)/(yni(m,2)*xmas(2))
	ss1(m) = (chiphi(m,1)/delna)*torvel(m+1,1)+yy1(m)
	ss2(m) = (chiphi(m,2)/delna)*torvel(m+1,2)+yy2(m)
	torvel(m,1) = (xnuc12(m)*ss2(m)+coefv(m,2)*ss1(m))/
	1			  (coefv(m,1)*coefv(m,2)-xnuc12(m)*xnuc21(m)) 
	torvel(m,2) = (xnuc21(m)*ss1(m)+coefv(m,1)*ss2(m))/
	1			  (coefv(m,1)*coefv(m,2)-xnuc12(m)*xnuc21(m))
c	form 2, impurity vphi eval from rad mom bal
	denom = coefv(m,1) - xnuc12(n)
	radterm = ti(m)*xlpm(m)/(bthet*(1/atnum(1)-1/zbar2(m)))-
	1		  (velthet1(m)-velthet2(m))/fp
	ss1(m) = (chiphi(m,1)/delna)*veltor(m+1,1)+yy1(m)
	veltor(m,1) = (ss1(m)-xnuc12(m)*radterm)/denom
	veltor(m,2) = veltor(m,1) - radterm

	return
	end

	
	