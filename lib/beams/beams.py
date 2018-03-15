#!/usr/bin/python

##########################################################
#
#   NBeams Function
#    
#   Generates the input files and calls NBeams
#
#   NBeamsMDS is located in the beams.py folder
#    
##########################################################

import os
from lib.funcs.dataGen import dataGrab,newScatPlot,smooth
from lib.funcs.physFuncs import zeffCalc
from lib.graphs.popups import popup
from subprocess import Popen, PIPE
import linecache
from math import pi,sqrt,exp,cos

class beamDataClass():
    
    # Defines a class to hold profiles
    
    def __init__(self,data):
        for k, v in data.items():
            setattr(self,k,v)


def dirCheck(beamPath):
    
    #    Verifies that NBeamsMDS.exe is where it should be
    
    if os.path.isfile(os.path.join(beamPath,"NBeamsMDS.exe")) == False:
        raise Exception("NBeamsMDS.exe not found")
    
    
def inbeamsGen(data):
    
    ###########################################################################
    #
    #   inbeams.dat generation modules
    #
    #   Writes an inbeams.dat input file based on user input from input files  
    #
    ###########################################################################    
    path=os.path.dirname(__file__)
    
    f=open(path+"/NBeamsMDS/Debug/inbeams.dat","w+")
    f.write("$nbin\n")
    f.write("nbeams = 1\n")
    f.write("inbfus = 1\n")
    f.write("amb = 2.0\n")
    f.write("zbeam = 1.0\n")
    f.write("ebeam = "+str(data.ebeam)+"\n")
    f.write("pbeam = "+str(data.pbeam)+"\n")
    f.write("rtang = "+str(data.rtang)+"\n")
    f.write("nbshape = 1\n")
    f.write("bwidth = 0.12\n")                                                   # Is this default?
    f.write("bheigh = 0.48\n")
    f.write("bgaussR = 0.066\n")
    f.write("bgaussZ = 0.18\n")
    f.write("bzpos = 0.0\n")
    f.write("nbptype = 1\n")
    f.write("maxiter = 2\n")
    f.write("pwrfrac(1,1) = "+str(data.pwrfrac1)+"   "+str(data.pwrfrac2)+"   "+str(data.pwrfrac3)+"\n")
    f.write("a = "+str(data.aminor)+"\n")
    f.write("r0 = "+str(data.rmajor)+"\n")
    f.write("b0 = "+str(data.bknot)+"\n")
    f.write("n = 51\n")
    f.write("e0 = "+str(data.epsknot)+"\n")
    f.write("ea = "+str(data.epssep)+"\n")
    f.write("shft0 = "+str(data.shftknot)+"\n")
    f.write("nion = 2\n")
    f.write("aion = 2.0 12.0\n")
    f.write("zion = 1.0 6.0\n")
    f.write("$end\n")
    f.close()
    
def nbdriverGen(te,ti,ne,ni):
    
    ###########################################################################
    #
    #   nbdriver.for generation module
    #
    #   Writes an nbdriver.fo input file based on user input from input files  
    #   Utilizes templates in the beams/ folder
    #
    ###########################################################################        

    path=os.path.dirname(__file__)

    f=open(path+"/NBeamsMDS/Debug/nbdriver.for","w+")
    a=open(path+"/nbdriverTemplate1.txt","r")
    b=open(path+"/nbdriverTemplate2.txt","r")
    x=a.readlines()
    f.writelines(x)
    for x in range(len(te)):
        lineList=["     ",str("tekev("+str(x+1)+")"),"=",str(te[x])]
        line=lineList[0]+lineList[1].rjust(5)+lineList[2].rjust(3)+lineList[3].rjust(8)+"\n"
        f.write(line)
    for x in range(len(ti)):
        lineList=["     ",str("tikev("+str(x+1)+")"),"=",str(ti[x])]
        line=lineList[0]+lineList[1].rjust(5)+lineList[2].rjust(3)+lineList[3].rjust(8)+"\n"
        f.write(line)        
    for x in range(len(ne)):
        lineList=["     ",str("ne20("+str(x+1)+")"),"=",str(ne[x])]
        line=lineList[0]+lineList[1].rjust(5)+lineList[2].rjust(3)+lineList[3].rjust(8)+"\n"
        f.write(line)   
    for x in range(len(ni)):
        lineList=["     ",str("ni20("+str(x+1)+","+"1)"),"=",str(ni[x])]
        line=lineList[0]+lineList[1].rjust(5)+lineList[2].rjust(3)+lineList[3].rjust(8)+"\n"
        f.write(line)
    y=b.readlines()
    f.writelines(y)
    f.close()
    a.close()
    b.close()
    
def depGrab(filePath):
    
    def dataPull(f,counter):
        depProf1=[]
        depProf2=[]
        depProf3=[]
        for a in range(51):
            line=linecache.getline(f+"outbeams.dat",counter+a+1)
            line=line.split()
            depProf1.append(float(line[1]))
            depProf2.append(float(line[2]))        
            depProf3.append(float(line[3])) 
        
        
        return depProf1,depProf2,depProf3
    
    
    f=open(filePath+"/outbeams.dat","r+")
    
    
    count=0
    
    for line in f:
        try: 
            if (line.split()[0] != "rho"):
                count=count+1
            else:
                count=count+1
                f.close()
                return dataPull(filePath,count)
        except:
            count=count+1
            pass
        
    
def run(nbiFile,te,ti,ne,ni,nc,zbar2,rmajor,elong,fastiol,run=True):


    ###########################################################################
    #
    #   NBeams main run routine
    #
    #   Pulls data from nbi file, normalizes profiles to keV/E20
    #   Interpolates to the 51 points necessary for NBeams
    #   Runs NBeams and outputs deposition profiles 
    #
    ###########################################################################  
  
    dicData=dataGrab(const=nbiFile)
    
    
    #   Turn dictionary-based data into class-based data and normalize
    
    beamData=beamDataClass(dicData)
    beamData.rhor=[0. + x*(1./(len(te)-1)) for x in range(len(te))]
    beamData.tekev=[x/1000. for x in te]
    beamData.tikev=[x/1000. for x in ti]
    beamData.ne20=[x/1.e20 for x in ne]
    beamData.ni20=[x/1.e20 for x in ni] 
    
    #   Establish rho vectors for 201 and 51 point datasets
    
    rhor=[0.+x*(1./(len(te)-1)) for x in range(len(te))]
    nbrho=[0.+x*(1./50) for x in range(51)]
    
    #   Build interpolation functions
    
    tekevinterp=smooth(rhor,beamData.tekev)
    tikevinterp=smooth(rhor,beamData.tikev)
    ne20interp=smooth(rhor,beamData.ne20)
    ni20interp=smooth(rhor,beamData.ni20)
    
    
    #   Create 51-point profiles for NBeams
    
    tekev,tikev,ne20,ni20=[],[],[],[]
    
    for a in range(51):
        
        tekev.append(float(tekevinterp(nbrho[a])))
        tikev.append(float(tikevinterp(nbrho[a])))
        ne20.append(float(ne20interp(nbrho[a])))
        ni20.append(float(ni20interp(nbrho[a])))
    
    #   Checks to see if NBeams is alive
    
    path=os.path.dirname(__file__)+"/NBeamsMDS/Debug/" 
    dirCheck(path)
    
    #   Generate inbeams and nbdriver
    
    inbeamsGen(beamData)
    nbdriverGen(tekev,tikev,ne20,ni20)
        
    
    ###########################################################################
    #
    #    Does NBeams need to be compiled here? Apparently not!
    #
    ###########################################################################
    
    tempcwd=os.getcwd()
    os.chdir(path)
    if run==True:    
        ps=Popen("NBeamsMDS.exe",stdin=PIPE)
    else:
        pass
        
    os.chdir(tempcwd)
        
    depProfs=depGrab(path)
    
    ###########################################################################
    #
    #   Interpolate deposition profiles to 201 points
    #
    ###########################################################################
    
    depProfInterp1=smooth(nbrho,depProfs[0])
    depProfInterp2=smooth(nbrho,depProfs[1])
    depProfInterp3=smooth(nbrho,depProfs[2]) 
    
    depProf1Full=[float(depProfInterp1(a)) for a in rhor]
    depProf2Full=[float(depProfInterp2(a)) for a in rhor]
    depProf3Full=[float(depProfInterp3(a)) for a in rhor]    
#    


    ###########################################################################
    #
    #   Calculate particle and energy sources
    #
    ###########################################################################
    
    beamData.vp = 19.72*rmajor*(beamData.aminor**2)*.5*(1+elong**2)/2.
    se=sqrt(.5*(1+elong**2))
    delma=1./(len(rhor)-1)
    
    volm = 4.*(pi**2)*rmajor*beamData.aminor*delma*se
    nbi,ebeam,alphain,abeam=75000.,beamData.ebeam,0.6475,2.                                   # Huh?
    
    Snbi=[((0.624E25*beamData.pbeam)/(beamData.vp*nbi))*beamData.pwrfrac1*nbd1+
    2.*beamData.pwrfrac2*nbd2+
    3.*beamData.pwrfrac3*nbd3 for nbd1,nbd2,nbd3 in zip(depProf1Full,depProf2Full,depProf3Full)]
    
    zeff=zeffCalc(ni,nc,zbar2)

    
    xlam1 = [5.5E17*ebeam/(abeam*0.5*a*(b**0.5)) for a,b in zip(ni,zeff)]  # Using ion density instead of dens() function here
    atten1 = [1-exp(-1*(delma/cos(alphain))/a) for a in xlam1]
    atten1.reverse()
    atten1.pop()
    atten1.insert(0,0)

    xlam2 = [5.5E17*(ebeam/2.)/(abeam*0.5*a*(b**0.5)) for a,b in zip(ni,zeff)]  # Using ion density instead of dens() function here
    atten2 = [1-exp(-1*(delma/cos(alphain))/a) for a in xlam2]
    atten2.reverse()
    atten2.pop()
    atten2.insert(0,0)
    
    xlam3 = [5.5E17*(ebeam/3.)/(abeam*0.5*a*(b**0.5)) for a,b in zip(ni,zeff)]  # Using ion density instead of dens() function here
    atten3 = [1-exp(-1*(delma/cos(alphain))/a) for a in xlam3]
    atten3.reverse()
    atten3.pop()
    atten3.insert(0,0)
        
    unatten1=[1.]
    unatten2=[1.]
    unatten3=[1.]  
    
    for a in range(200):
        unatten1.append(unatten1[-1]*(1-atten1[a]))
        unatten2.append(unatten2[-1]*(1-atten2[a]))    
        unatten3.append(unatten3[-1]*(1-atten3[a]))    
        
    xpartdot1=[a*b*(((beamData.pwrfrac1*beamData.pbeam*1.E6)/(1.6E-19*(ebeam/1.)*1.E3)))/volm for a,b in zip(unatten1,atten1)]    
    xpartdot2=[a*b*(((beamData.pwrfrac2*beamData.pbeam*1.E6)/(1.6E-19*(ebeam/2.)*1.E3)))/volm for a,b in zip(unatten2,atten2)]    
    xpartdot3=[a*b*(((beamData.pwrfrac3*beamData.pbeam*1.E6)/(1.6E-19*(ebeam/3.)*1.E3)))/volm for a,b in zip(unatten3,atten3)]    
    
    # qnb calculation. Note that NBIeloss = 0, reducing this formula significantly
    
    qnb=[(a+b/2.+c/3.)*1.6E-19*1.E3*ebeam for a,b,c in zip(xpartdot1,xpartdot2,xpartdot3)]  # What is this 80?                                   
    
    #   TODO: Figure this division out
    
    qnbi=qnb
    qnbe=0.
    
    #######################################################################
    #   Calculate toroidal momentum source rates
    #######################################################################
 
    NBIspin=1.
    abeam = 2.
    xmbeam= abeam*1.673E-27 
    
    torque1=[sqrt(2.*xmbeam/(1.6E-19*ebeam*1.E3))*beamData.rtang*beamData.pwrfrac1*beamData.pbeam*1.E6*(1.-a)*NBIspin for a in fastiol[0]]
    torque2=[sqrt(2.*xmbeam/(1.6E-19*ebeam*1.E3/2.))*beamData.rtang*beamData.pwrfrac2*beamData.pbeam*1.E6*(1.-a)*NBIspin for a in fastiol[1]]
    torque3=[sqrt(2.*xmbeam/(1.6E-19*ebeam*1.E3/3.))*beamData.rtang*beamData.pwrfrac3*beamData.pbeam*1.E6*(1.-a)*NBIspin for a in fastiol[2]]
    
    xmphi1=[a*b*c/rmajor for a,b,c in zip(unatten1,atten1,torque1)]
    xmphi2=[a*b*c/rmajor for a,b,c in zip(unatten2,atten2,torque2)]
    xmphi3=[a*b*c/rmajor for a,b,c in zip(unatten3,atten3,torque3)]

    totmomphi = [(a+b+c)/volm for a,b,c in zip(xmphi1,xmphi2,xmphi3)]
    
    #   Determine momentumd istribution between deterium and carbon
    #   Distribute equally

    frac=[b/(a+b) for a,b in zip(ni,nc)]
    
    # Calculate momentum source terms, leaving room for anomolous torque
    
    # TODO: Anomolous torque
    anomtorq=[0.]*len(ne)
    
    xmomtor1=[(1-a)*(b+c)  for a,b,c in zip(frac,totmomphi,anomtorq)]
    xmomtor2=[a*(b+c) for a,b,c in zip(frac,totmomphi,anomtorq)]
    
    return [depProf1Full,depProf2Full,depProf3Full],Snbi,qnbi,qnbe,xmomtor1,xmomtor2
    
    
    
    
    
    