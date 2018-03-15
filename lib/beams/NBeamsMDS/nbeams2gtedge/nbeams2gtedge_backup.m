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
d = [.59976E+01
.12143E+02
.80229E+01
.71124E+01
.66318E+01
.62760E+01
.59631E+01
.56653E+01
.53718E+01
.50801E+01
.38927E+01
.33362E+01
.29553E+01
.26654E+01
.24326E+01
.22403E+01
.20772E+01
.19364E+01
.18135E+01
.17042E+01
.16060E+01
.15172E+01
.14361E+01
.13612E+01
.12916E+01
.12268E+01
.11660E+01
.11083E+01
.10534E+01
.10011E+01
.95065E+00
.90222E+00
.85559E+00
.81030E+00
.76611E+00
.72330E+00
.68151E+00
.64051E+00
.60043E+00
.56119E+00
.52282E+00
.48530E+00
.44878E+00
.41370E+00
.38213E+00
.35185E+00
.31212E+00
.27004E+00
.22135E+00
.17225E+00
.00000E+00
];

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
1.00000
];

Pbeam = 8.66; %NBI Power [MW]
Enbi = 74000; %Beam energy [eV]
%--------------------------------------------------------------------------
%Define different rho vectors
rhoGT = [0.864,0.870,0.876,0.881,0.887,0.893,0.898,0.904,0.910,0.915,0.921,0.927,0.932,0.938,0.944,0.9490,...
        0.955,0.960,0.966,0.972,0.977,0.983,0.989,0.994,1.000]';
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
    