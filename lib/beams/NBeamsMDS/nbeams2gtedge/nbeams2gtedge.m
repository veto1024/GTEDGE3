close all
clear all
clc
%converts NBeams deposition profile ranging from rho = [0,1] to rho =
%[0.864,1] for GTEDGE.
%Uses a quadratic fit for data up to n-1 points and a line for the last 2
%points
%change 'd' to hofr_1 from NBeams
%output nbdep goes directly into soldata input file
%output nbsr is neutral beam source rate to core before IOL is calculated
%--------------------------------------------------------------------------
%INPUTS
%define deposition profile from NBeams
d = [.33884E+01
.69976E+01
.46200E+01
.40910E+01
.38046E+01
.35969E+01
.34063E+01
.32280E+01
.30513E+01
.28753E+01
.27034E+01
.25336E+01
.23686E+01
.22099E+01
.20599E+01
.19195E+01
.16895E+01
.15316E+01
.14151E+01
.13214E+01
.12432E+01
.11768E+01
.11192E+01
.10699E+01
.10253E+01
.98532E+00
.95078E+00
.92115E+00
.89538E+00
.87380E+00
.85755E+00
.84251E+00
.83036E+00
.82131E+00
.81511E+00
.81192E+00
.81005E+00
.81127E+00
.81006E+00
.80902E+00
.80842E+00
.80700E+00
.80216E+00
.79338E+00
.78387E+00
.76855E+00
.74234E+00
.68778E+00
.56207E+00
.32558E+00
.00000E+00];

rho = [ .00000
 .02000
 .04000
 .06000
 .08000
 .10000
 .12000
 .14000
 .16000
 .18000
 .20000
 .22000
 .24000
 .26000
 .28000
 .30000
 .32000
 .34000
 .36000
 .38000
 .40000
 .42000
 .44000
 .46000
 .48000
 .50000
 .52000
 .54000
 .56000
 .58000
 .60000
 .62000
 .64000
 .66000
 .68000
 .70000
 .72000
 .74000
 .76000
 .78000
 .80000
 .82000
 .84000
 .86000
 .88000
 .90000
 .92000
 .94000
 .96000
 .98000
1.00000];

Pbeam = 2.186; %NBI Power [MW]
Enbi = 76910; %Beam energy [eV]
%--------------------------------------------------------------------------
%Define different rho vectors
rhoGT = [0.85563499
0.86164999
0.86766499
0.87368000
0.87969601
0.88571101
0.89172602
0.89774102
0.90375602
0.90977198
0.91578698
0.92180198
0.92781699
0.93383300
0.93984801
0.94586301
0.95187801
0.95789301
0.96390897
0.96992397
0.97593898
0.98195398
0.98796999
0.99398500
1.0000000
]';
rhoGT=rhoGT.'
rhoNB = [0:0.04:1.0]';

%cut NB output down to coarse mesh in edge region
cd = d(end-4:end);
crho = rhoNB(end-4:end);

%fit parabola to n-1 coarse deposition points
coefs1 = polyfit(crho(1:end-1),cd(1:end-1),2);
parabola = polyval(coefs1,rhoGT(1:end-7));

%fit last 2 points of coarse deposition to line
coefs2 = polyfit(crho(end-1:end),cd(end-1:end),1);
line = polyval(coefs2,rhoGT(end-7:end));

%plots
nbdep = [parabola(1:end-1);line];
figure(1)
hold on
plot(rhoNB(end-4:end),d(end-4:end),'b','LineWidth',2)
plot(rhoGT,nbdep,'r','LineWidth',2)
title('Neutral Beam Deposition Profile H(\rho) in the Edge','FontSize',30)
xlabel('Normalized Radius \rho','FontSize',30)
ylabel('H(\rho)','FontSize',30)

figure(2)
plot(rho,d,'k','LineWidth',2)
title('Total Neutral Beam Deposition Profile H(\rho)','FontSize',30)
xlabel('Normalized Radius \rho','FontSize',30)
ylabel('H(\rho)','FontSize',30)

%--------------------------------------------------------------------------
%print output
num = 1:length(d);
% for k = 1:length(num)
%     fprintf(1,'nbdep(%0.0f) = %8.4f\n',num(k),nbdep(k))
% end

for k = 1:length(num)
    fprintf(1,'      NBrho(%0.0f) = %8.4f\n',num(k),rho(k))
end

for k = 1:length(num)
    fprintf(1,'      nbdep(%0.0f) = %8.4f\n',num(k),d(k))
end

%nbsr = 0.624e25*Pbeam/Enbi*trapz(rhoNB,d);
%fprintf(1,'nbsr = %8.3e',nbsr)
    