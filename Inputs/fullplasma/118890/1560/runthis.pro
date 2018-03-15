PRO runthis,g
g=readg(118890,1560)
;g=readg(118890,1515)
;g=readg(166606,1950)
;g=readg(149468,1905)

calculate_bfield,bp,brad,bt,bvert,g

openw, 2,'btor.txt'
printf, 2, bt
close, 2

openw, 2,'bpol.txt'
printf, 2, bp
close, 2

openw, 2,'R.txt'
printf, 2, g.R
close, 2

openw, 2,'Z.txt'
printf, 2, g.Z
close, 2

openw, 2,'BDRY.txt'
printf, 2, g.BDRY
close, 2

openw, 2,'lim.txt'
printf, 2, g.LIM
close, 2

openw, 2,'psirz.txt'
printf, 2, g.PSIRZ
close, 2

openw, 2,'EPOTEN.txt'
printf, 2, g.EPOTEN
close, 2
end

