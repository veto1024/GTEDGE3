#!/usr/bin/python

###########################################################
#
#   Miscellaneous functions for generating profiles needed
#   for various calculations in GTEDGE
#
#   Includes the following calculations:
#
#   cxcoolCalc  - Calculates cooling via charge exchange
#   zeffCalc    - Calculates zeffective
#   qieCalc     - Calculates heating from ions to electrons
#   xniCalc     - Calculates dueterium density from zbar2/experimental nz1
#   xnCcalc     - Calculates carbon density from zbar2/nz1, 
#                 possibly needs to be revised for multi-impurity plasmas
#   coulCalc    - Calculates coulomb logarithm
#   xnucCalc    - Calculates scattering between 2 species
#   xnuatiCalc  - Calculates atomic scattering
#   zbar2Calc   - Calculates average charge state for carbon impurity
#
###########################################################

from math import *
import numpy as np

def cxcoolCalc(xni,xti,xnuati):
    
    
    
    xk=1.6E-19
    cxcool=[]
    
    #  Set inner-most flux surface to have no charge exchange
    #  This allows the average temeprature to go to previous 
    #  point for tiav calculation and thus give a full profile
    
    cxcool.append(0.)
    dens=[xni[0]]
    tiav=[xti[0]]
    for n in range(1,len(xni)):
        dens.append((xni[n]+xni[n-1])/2.)
        tiav.append((xti[n]+xti[n-1])/2.)
    cxcool=[1.5*a*b*xk*c for a,b,c in zip(dens,tiav,xnuati)]
    
    return cxcool
    
def zeffCalc(xni,xnC,zbar2):
    
    zeff=[]
    atnum=1.0
    for xni,xnC,zbar2 in zip(xni,xnC,zbar2):
        zefftemp=(xni*(atnum**2)+xnC*(zbar2**2))/(xni*atnum+xnC*zbar2)
        zeff.append(zefftemp)
    
    return zeff
     
def qieCalc(xni,xnC,xne,xmas1,xte,xti,zbar2):
    
    ###########################################################################
    #   
    #    Rate of heat transfer from ions to electrons
    #    zbar2 in GTEDGE used 6.0
    #    This zbar2 is experimentally inferred zbar2    
    #
    ###########################################################################
    
    #  Atomic number of ion and impurities
    
    atnum=1.0
    eq1=1.6E-19
    ep0 = 8.854E-12         # What is this?  
    qie=[0.]*len(xni)
    
    zeff=zeffCalc(xni,xnC,zbar2)
    
    for n in range(len(xni)-1):
        
        # TODO: Should verify this calculation
        
        teav = 0.5*(xte[n]+xte[n+1])
        tiav = 0.5*(xti[n]+xti[n+1])
        a = 0.5*(xni[n]+xni[n+1])*atnum + 0.5*(xnC[n]+xnC[n+1])*zbar2[n]   
        b = (ep0/eq1)**1.5   
        couloge = log(12.*pi*(teav**1.5)*b/sqrt(a))
        cequil = 7.9e-42*couloge*zeff[n]/xmas1

        qie[n] = cequil*xne[n]*(tiav-teav)/(teav**1.5)
     
    return qie
     
def  xniCalc(xne,zbar2,fracz):
    
    ###########################################################################
    #
    #   This might be garbage
    #
    ###########################################################################
    
    xni = [a/(1.+b*c) for a,b,c in zip(xne,fracz,zbar2)]
    return xni

def xnCcalc(fz1,zbar2):
    xnC = [a/b for a,b in zip(fz1,zbar2)]
    return xnC

def coulCalc(atnum1,atnum2,zbar2,xti,xnC):                                     # To be verified
    
    ###########################################################################
    #
    #    Coulomb logarithm calculation
    #
    #    atnums:
    #    1. = Deuterium
    #    6. = Carbon
    #
    #    Returns a list for the given logarithm
    ###########################################################################

    if float(atnum1)==1.: z1=[1.]*len(zbar2)
    else: z1 = [6.]*len(zbar2)
    if float(atnum2)==1: z2 = [1.]*len(zbar2)
    else: z2 = [6.]*len(zbar2)
    
    ep0=8.854E-12
    eq=1.6E-19
    for x in xnC:
        if x<=1E9: 
            raise ValueError("Carbon density too low in Coulomb Logarithm")
            
    a = [(ep0/eq)**1.5] * len(zbar2)
    b = [sqrt(xnC)*(y**2)*z for xnC,y,z in zip(xnC,z2,z1)]
    
    return [log(12.*pi*(xti**1.5)*m/n) for xti,m,n in zip(xti,a,b)]                  
    
    
    
def xnucCalc(data):
    
    ###########################################################################
    #
    #    Calculation of scattering cross sections
    #
    ###########################################################################
    atnum=[1.,6.]

    def xmrCalc(data,a,b):
        if a==1: mass1=data.xmas1
        else: mass1=data.xmas2
        
        if b==1: mass2=data.xmas1
        else: mass2=data.xmas2   
        
        return mass1*(1+(mass1/mass2))
    
    zbar2=data.zbar2
    xti=data.xti
    xni=data.xni
    xnC=data.xnC
    
    xmrArray=np.array([[xmrCalc(data,1,1),xmrCalc(data,1,2)],
                       [xmrCalc(data,2,1),xmrCalc(data,2,2)]])
    coulArray=np.array([[coulCalc(1.,1.,zbar2,xti,xnC),coulCalc(1.,6.,zbar2,xti,xnC)],
                        [coulCalc(6.,1.,zbar2,xti,xnC),coulCalc(6.,6.,zbar2,xti,xnC)]])
    
    C1 = 1./((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5)) 
     
    xnucArray=np.zeros((2,2),dtype=np.ndarray)

    xnucArray[0,0]=[3.34*(c*(atnum[0]**4)*1.E-6*b)/(C1*sqrt(xmrArray[0,0]*1.E3)*(a**1.5)) for a,b,c in zip(xti,xni,coulArray[0,0])]
    xnucArray[0,1]=[3.34*(d*(atnum[0]*c)**2)*1.E-6*a/(C1*sqrt(xmrArray[0,1]*1.E3)*(b**1.5)) for a,b,c,d in zip(xnC,xti,zbar2,coulArray[0,1])]
    xnucArray[1,0]=[3.34*(d*(atnum[0]*c)**2)*1.E-6*a/(C1*sqrt(xmrArray[1,0]*1.E3)*(b**1.5)) for a,b,c,d in zip(xni,xti,zbar2,coulArray[1,0])]
    xnucArray[1,1]=[3.34*(d*(a**4)*1.E-6*c)/(C1*sqrt(xmrArray[1,1]*1.E3)*(b**1.5)) for a,b,c,d in zip(zbar2,xti,xnC,coulArray[1,1])]
    
    return xnucArray
    
#TODO def xnuatiCalc(args**):
#    ####################################
#    # Give me code, max
#    ####################################
#
#    
#
#
#
#    return xnuati
    
    def zbar2Calc(data):
        
        pass
        
        return zbar2
     