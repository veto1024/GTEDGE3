#!/usr/bin/python

###########################################################
#
#   Calculation of radial particle flux
#   Uses 25-point GTEDGE model
#
#   iolFlag = 1 w/ IOL, 0 w/out IOL
#
#   6/12/2017
#
#   Particle flux and heat flux for deuterium implemented
#   
#
###########################################################

###########################################################
#
#    WARNING: DEPRECATED
#
###########################################################
from math import *
from physFuncs import cxcoolCalc,qieCalc

def fluxCalc(data):
    
    gamma=[0.]*25
    gammahat=[0.]*25
    qheati=[0.]
    qhatheati=[0.]
    qheate=[0.]
    gammaC=[0.]*25
    gammaHeati=[0.]*25
    gammahatHeati=[0.]*25
    enh=1.0                                                                    # What are you?
    
    aminor,elong,fracz,delma,xk,BCnbi,fluxheat,fheate,xmas1,fluxpart=data['aminor'],data['elong'],data['fracz'],data['delma'],data['xk'],data['BCnbi'],data['fluxheat'],data['fheate'],data['xmas1'],data['fluxpart']
    rhor,sNBI,xnuioni,zbar2,xnuioni,xni,xnC,xne,xti,xte,dlnn_dt,dlnwi_dt,xnuati,fionb,qnb,eorbl=data['rhor'],data['Snbi'].values(),data['xnuioni'].values(),data['zbar2'].values(),data['xnuioni'].values(),data['xni'].values(),data['xnC'].values(),data['xne'].values(),data['xti'].values(),data['xte'].values(),data['dlnn_dt'].values(),data['dlnwi_dt'].values(),data['xnuati'].values(),data['fionb'].values(),data['qnb'].values(),data['eorbl'].values()
    
    cxcool=cxcoolCalc(xni,xti,xnuati)
    cmulteq=1.0                                                                # Why?
    qie=qieCalc(xni,xnC,xne,xmas1,xte,xti,zbar2)
    forbl=data['forbl'].values()

    forbl[24]=forbl[23]
    radminor = aminor*sqrt(0.5*(1+elong**2))
    gamma[24]=enh*fluxpart
    
    
    for n in range(23,-1,-1):
        
        radminor=aminor*sqrt(.5*(1+elong**2))
        #########################################################
        #
        #  Particle flux calculation 
        #  
        #  sNBI includes IOL corrections, so unsure of validity,
        #  but the magnitude is so small that may be insignificant
        #           
        #########################################################
        
        srprim=sNBI[n]+((xni[n]+xni[n+1])/2)*xnuioni[n]*(1+fracz*zbar2[n])
        
        gamma[n]=(rhor[n+1]/rhor[n])*(gamma[n+1]*1.-srprim*delma)
        
 
    # Time-dependant corrections
        dnped_dt = 0.5*(dlnn_dt[n]*xni[n] +	dlnn_dt[n+1]*xni[n+1])
        gamma[n]=gamma[n]+dnped_dt*delma
        
    # Change in BCs
      
    yn=0.
    for n in range(24):
        yn=yn + 0.5*(dlnn_dt[n]+dlnn_dt[n+1])*xni[n]*delma       
    gamma=[x - yn for x in gamma]

    
    ########################################################
    #   Ion heat flux calculation
    #   Differential slab geometry
    #
    #   gammaHeati does not include NBI or IOL losses
    #   Don't question it, run by Stacey
    ########################################################
        
#    gammaHeate[24]=fluxheat*fheate
    gammaHeati[24]=(1.-fheate)*fluxheat                                    # fheate simply given as 0.4 in GTEDGE
    for n in range(23,-1,-1):
        gammaHeati[n] = gammaHeati[n+1] + delma*(cxcool[n] + cmulteq*qie[n])

    # Time-dependent corrections
        uni = 0.25*(xni[n]+xni[n+1])*xk*(xti[n]+xti[n+1])
        dwiped_dt= 0.5*(dlnwi_dt[n]+dlnwi_dt[n+1])*1.5*uni 
        
        gammaHeati[n] = gammaHeati[n]+dwiped_dt*delma
    
    # Change in BCs
    uni,yi=0.,0.
    for n in range(24):
        uni = 0.25*(xni[n]+xni[n+1])*(xti[n]+xti[n+1])*xk
        yi = yi + 0.5*(dlnwi_dt[n]+dlnwi_dt[n+1])*1.5*uni*delma
    
    gammaHeati=[x - yi for x in gammaHeati]

    #############################################################
    #
    #  Flux calculations with IOL
    #
    ################################################################
    
    radminor=aminor*sqrt(.5*(1+elong**2))
    
    #  Using BCnbi as the SNbi from n=0
    
    gammahat[0]=BCnbi + ((xni[0]+xni[1])/2.)*xnuioni[0]*(1.+fracz*zbar2[0])*0.02*radminor


    xponDB=[]
    srprimDB=[]
    for n in range(1,25):

#            radiusInt=rhor[n]*radminor
        #########################################################
        #  Particle flux calculation            
        #########################################################
        xpon = exp(-2.*(forbl[n]-forbl[n-1]))           
        srprim=sNBI[n]+((xni[n]+xni[n-1])/2.)*xnuioni[n]*(1+fracz*zbar2[n])
 #       xponDB.append(xpon)
#        srprimDB.append((xni[n]+xni[n-1]/2.))
#        srprimDB.append(srprim)
        gammahat[n]=(rhor[n-1]/rhor[n])*gammahat[n-1]*xpon+srprim*delma
        
        
    # Time-dependant corrections
        
    uni,yn,yi=0.,0.,0.
    
    for n in range(24):
        uni=0.25*(xni[n]+xni[n+1])*(xti[n]+xti[n+1])*xk #Trace me
        yn = yn + 0.5 *(dlnn_dt[n]+dlnn_dt[n+1])*xni[n]*delma
    
    gammahat=[x-yn for x in gammahat]
    
    ########################################################
    #   Ion heat flux calculation
    #   Differential slab geometry
    # 
    #   This is my guess for how to do this
    ########################################################
    gammaDB=[]    
    eorbl[24]=eorbl[23]
    BCenbi=0.
    gammahatHeati[0]=gammaHeati[0]-BCenbi
    for n in range(23,-1,-1):
        xpon=exp(-1.*(eorbl[n]-eorbl[n+1]))                                    # So... for n+1 = 25, what do?

        srprim=qnb[n]-(qie[n]+cmulteq*cxcool[n])                                # Has qnb been reduced by FIOL? 
        gammahatHeati[n] = (rhor[n+1]/rhor[n])*gammahatHeati[n+1]*xpon - delma*srprim
        srprimDB.append(srprim)
        gammaDB.append((rhor[n+1]/rhor[n]))
    
    # Time-dependent corrections
        uni = 0.25*(xni[n]+xni[n+1])*xk*(xti[n]+xti[n+1])
        dwiped_dt= 0.5*(dlnwi_dt[n]+dlnwi_dt[n+1])*1.5*uni 
        
        gammahatHeati[n] = gammahatHeati[n]+dwiped_dt*delma
    
    for n in range(23,-1,-1):
        	uni = 0.25*(xni[n]+xni[n+1])*(xti[n]+xti[n+1])*xk
         	yi = yi + 0.5*(dlnwi_dt[n]+dlnwi_dt[n+1])*1.5*uni*delma
    
    gammahatHeati=[x - yi for x in gammahatHeati]  
    
    return gamma,gammahat,gammaHeati,gammahatHeati
