#!/usr/bin/python

############################################3
#
#
#    Experimetnal Chi calculation functions
#
#
############################################

from lib.funcs.dataGen import newScatPlot
from math import sqrt

class chiClass():
    
    def __init__(self,data):
        self.atnum1=1.
        #self.chie=self.chieCalc(data)   
        self.qtotal=data.qHeati
        self.heatin=self.heatinCalc(data)
        self.heatvisc=self.heatviscCalc(data)
        self.conv15=self.conv15calc(data)
        self.conv25=self.conv25calc(data)


        
    def chieCalc(self,data):
        gameltemp=[self.atnum1*b+c*d for b,c,d in zip(data.gamma,data.zbar2,data.gamC)]
        return [a*((b/(c*data.xk*f))-1.5*g/h) for a,b,c,f,g,h in zip(data.exlte,data.gamheate,data.xne,data.xte,gameltemp,data.xne)]


###############################################################################
#
#   Main Chi calculation
#   conv15 = True = subtract away convective heat flux with 3/2 coefficient in front
#   conv25 = True = substract away convective heat flux with 1/2 coefficient in front
#   visc = True = subtract away viscous heat flux
#   pressure = True = subtract away work done by the plasma on the pressure tensor
#
#   conv15 and conv25 can't be simultaneously true
###############################################################################

    def chiiCalc(self,data,conv15=False,conv25=False,visc=False,pressure=False):
        qcond=data.qHeati
        

        if ((conv15==True) and (conv25==True)):
            raise("Can't use Conv15 and Conv25 simultaneously")
            
        if conv25==True:
            qcond=[a - b for a,b in zip(qcond,self.conv25)]
            
        if conv15==True:
            qcond=[a - b for a,b in zip(qcond,self.conv15)] 
            
        if pressure==True:
            qcond=[a - b for a,b in zip(qcond,self.heatin)]
        
        if visc==True:
            qcond=[a-b for a,b in zip(qcond,self.heatvisc)]
        
        return [d*(a)/(b*c*data.xk) for a,b,c,d in zip(data.exlti,data.xni,data.xti,qcond)]
        
    def conv15calc(self,data):
        return [1.5*data.xk*a*b for a,b in zip(data.gamma,data.xti)]
        
    def conv25calc(self,data):
        return [2.5*data.xk*a*b for a,b in zip(data.gamma,data.xti)]
        
    def heatviscCalc(self,data):
        se=sqrt(0.5*(1+data.elong**2))
        ep=data.aminor*se/data.rmajor
        
        fp=[a/data.bphi for a in data.bthet]

        xnustar11=0.

        # QUESTIONABLE CALCULATION COMMENTED OUT
#        for a,b in zip(data.xnuc[0,1],data.vpolD):
#            xnustar11=xnustar11+a*abs(data.q95)*data.rmajor/b
#
#       Also, why the fuck is it xnuc12 in GTEDGE?! WHAT ARE YOU DOING WITH YOUR LIFE, STACEY?!

        xnustar11=[a*abs(data.q95)*data.rmajor/b for a,b in zip(data.xnuc[0,0],data.vpolD)]    
        
        eff=[a/((1+a)*((ep**1.5)+a)) for a in xnustar11]
        vrad1=[a/b for a,b in zip(data.gamma,data.xni)]
          
        eta0=[a*data.xmas1*b*c*data.rmajor*d for a,b,c,d in zip(data.xni,data.vpolD,data.q,eff)]
        eta4=[(a*data.xmas1*b*data.xk)/(data.eq1*abs(data.bphi)) for a,b in zip(data.xni,data.xti)]
   
       
#       Calculate viscous heating:  a=vtord, b=fp, c = eta0, d=vrad1, f = eta4, g= vpold
#   TODO: THIs does not match POP17-052504

        return [a*(b*c*d-0.5*f*(4.*a+g))-0.5*g*(c*d+f*(a+0.5*g)) for a,b,c,d,f,g in zip(data.vtorD,fp,eta0,vrad1,eta4,data.vpolD)]
             
    def heatinCalc(self,data):
        return [a*0.5*data.xmas1*(b**2+c**2) for a,b,c in zip(data.gamma,data.vtorD,data.vpolD)]
        
 
class chiorbClass():

    
    def __init__(self,data):
        self.atnum1=1.
        #self.chie=self.chieCalc(data)   
        self.qtotal=data.qhatHeati
        self.heatin=self.heatinCalc(data)
        self.heatvisc=self.heatviscCalc(data)
        self.conv15=self.conv15calc(data)
        self.conv25=self.conv25calc(data)


        
    def chieCalc(self,data):
        gameltemp=[self.atnum1*b+c*d for b,c,d in zip(data.gammahat,data.zbar2,data.gamhatC)]
        return [a*((b/(c*data.xk*f))-1.5*g/h) for a,b,c,f,g,h in zip(data.exlte,data.gamheate,data.xne,data.xte,gameltemp,data.xne)]


###############################################################################
#
#   Main Chi-hat calculation
#   conv15 = True = subtract away convective heat flux with 3/2 coefficient in front
#   conv25 = True = substract away convective heat flux with 1/2 coefficient in front
#   visc = True = subtract away viscous heat flux
#   pressure = True = subtract away work done by the plasma on the pressure tensor
#
#   conv15 and conv25 can't be simultaneously true
###############################################################################

    def chiiCalc(self,data,conv15=False,conv25=False,visc=False,pressure=False):
        qcond=data.qhatHeati
        

        if ((conv15==True) and (conv25==True)):
            raise("Can't use Conv15 and Conv25 simultaneously")
            
        if conv25==True:
            qcond=[a - b for a,b in zip(qcond,self.conv25)]
            
        if conv15==True:
            qcond=[a - b for a,b in zip(qcond,self.conv15)] 
            
        if pressure==True:
            qcond=[a - b for a,b in zip(qcond,self.heatin)]
        
        if visc==True:
            qcond=[a-b for a,b in zip(qcond,self.heatvisc)]
        
        return [d*(a)/(b*c*data.xk) for a,b,c,d in zip(data.exlti,data.xni,data.xti,qcond)]
        
    def conv15calc(self,data):
        return [1.5*data.xk*a*b for a,b in zip(data.gammahat,data.xti)]
        
    def conv25calc(self,data):
        return [2.5*data.xk*a*b for a,b in zip(data.gammahat,data.xti)]
        
    def heatviscCalc(self,data):
        se=sqrt(0.5*(1+data.elong**2))
        ep=data.aminor*se/data.rmajor
        
        fp=[a/data.bphi for a in data.bthet]
        xnustar11=0.

        # QUESTIONABLE CALCULATION COMMENTED OUT
#        for a,b in zip(data.xnuc[0,1],data.vpolD):
#            xnustar11=xnustar11+a*abs(data.q95)*data.rmajor/b
#
#       Also, why the fuck is it xnuc12 in GTEDGE?! WHAT ARE YOU DOING WITH YOUR LIFE, STACEY?!

        xnustar11=[a*abs(data.q95)*data.rmajor/b for a,b in zip(data.xnuc[0,0],data.vpolD)]    
        
        eff=[a/((1+a)*((ep**1.5)+a)) for a in xnustar11]
        vrad1=[a/b for a,b in zip(data.gammahat,data.xni)]
          
        eta0=[a*data.xmas1*b*c*data.rmajor*d for a,b,c,d in zip(data.xni,data.vpolD,data.q,eff)]
        eta4=[(a*data.xmas1*b*data.xk)/(data.eq1*abs(data.bphi)) for a,b in zip(data.xni,data.xti)]
   
       
#       Calculate viscous heating:  a=vtord, b=fp, c = eta0, d=vrad1, f = eta4, g= vpold
#   TODO: THIs does not match POP17-052504
        asymR=0.1/data.rmajor
        return [asymR*a*(b*c*d-0.5*f*(4.*a+g))-0.5*g*(c*d+f*(a+0.5*g)) for a,b,c,d,f,g in zip(data.vtorDhat,fp,eta0,vrad1,eta4,data.vpolD)]
             
    def heatinCalc(self,data):
        return [a*0.5*data.xmas1*(b**2+c**2) for a,b,c in zip(data.gammahat,data.vtorDhat,data.vpolD)]
     
