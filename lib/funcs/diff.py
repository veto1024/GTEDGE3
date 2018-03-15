#!/usr/bin/python

############################################3
#
#
#    Experimetnal Diffusion coefficient calculation
#    functions
#
#
############################################

def diffCalcs(data):
    
    diffCoNF=[]
    
     
    for xnuc12,xti,nudrag,bthet,gammahat,xni,exlti in zip(data.xnuc[0,1],data.xti,
                                        data.nudrag,data.bthet,
                                        data.gammahat,
                                        data.xni,
                                        data.exlti):
                                            
        diffCoNFtemp=data.xmas1*xnuc12*xti*(1+(nudrag/xnuc12)-(1./6.))/(data.eq1*(bthet**2))
        diffCoNF.append(diffCoNFtemp)

    diffCoLame=[-1.*a*c/b for a,b,c in zip(data.gamma,data.xni,data.exlni)]
    
    
    return diffCoNF,diffCoLame
        