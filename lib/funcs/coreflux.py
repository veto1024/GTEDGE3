#!/usr/bin/python

##########################################################
#
#    Calculation of the particle flux
#    for the full plasma for GTEDGE V3
#
#    
##########################################################

from math import exp
from lib.funcs.physFuncs import qieCalc,cxcoolCalc,xniCalc
from lib.funcs.dataGen import dictAppend


def snbiIOL(snbi,fiolProf):
    
    #   Reduce NBI particle source term by IOL
    #   Check to see if this is done right because IOL gives dF/drho, which
    #   is updated from (1-Fiol)
    
    snbi = [a*(1-b) for a,b in zip(snbi,fiolProf)]
    return snbi    
    

def fluxCore(coreData,iolFlag=False):

    ##########################################################
    #
    #    Main calculation of the radial particle flux
    #    for the full plasma
    #
    #    iolFlag turns fast and thermal IOL on/off
    #
    ##########################################################
    

    xne=coreData.xne
    xni=coreData.xni
    xnC=coreData.xnC
    zbar2=coreData.zbar2
    xti=coreData.xti
    xte=coreData.xte
    fracz=coreData.fracz
    snbi=coreData.Snbi
    qnbi=coreData.qnbi
    rhor=coreData.rhor
    forbl=coreData.fiol
    eorbl=coreData.eiol
    xnuati=coreData.xnuati
    xnuioni=coreData.xnuioni
    Fast_iol=[0.]*len(coreData.rhor)
        
    
    mesh=len(coreData.rhor)
    
    if iolFlag==True:
        snbi=snbiIOL(snbi,Fast_iol)

        

    ###########################################################################
    
    
    cxcool=cxcoolCalc(xni,xti,xnuati)
    cmulteq=1.0                                                                # Why?
    qie=qieCalc(xni,xnC,xne,coreData.xmas1,xte,xti,zbar2)
    
    delma=1./(mesh-1)
    
    gamma=[0.]*mesh
    qheati=[0.]*mesh
    qheate=[0.]*mesh
    elong=coreData.elong
    aminor=coreData.aminor
    
#    rhor=[0.+a*1./mesh for a in range(mesh+1)]
#    for n in range(mesh-2,0,-1):
#        se=sqrt(0.5*(1.+elong**2))
#        
#        rhor[n] = rhor[n+1] - delma/(aminor*se)

    
    for n in range(2,mesh):
        
        #####################################################
        #
        #   Particle flux calculation w/ flags to deactive
        #   IOL if necessary
        #   
        #
        #####################################################
        dens = (xni[n]+xni[n-1])/2.
        srprim=snbi[n] + dens*xnuioni[n]*(1.+fracz[n]*zbar2[n])
        xpon=exp(-2.*(forbl[n]-forbl[n-1]))
        if iolFlag==False: xpon=1.
        gamma[n]=(rhor[n-1]/rhor[n])*gamma[n-1]*xpon+srprim*delma                  # What is delma? Is this correct here?
        
        xponq=exp(-1.*(eorbl[n]-eorbl[n-1]))
        if iolFlag==False: xponq=1.
        srprimq=qnbi[n]-cxcool[n]-1.*qie[n]
        qheati[n]=(rhor[n-1]/rhor[n])*qheati[n-1]*xponq+srprimq*delma
        
        # TODO: electron heating
        
    return gamma,qheati,qheate
        
        
        
        
        
        
    