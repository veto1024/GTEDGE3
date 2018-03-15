#!/usr/bin/python

##############################################
#
#   Calculations of nudrag using old perturbation 
#   theory
#
#   DEPRECATED
#
##############################################

# Nudrag perturbation calculation with IR correction
def nuDragPert(switch,data):
    delv1=[]
    delv0=[]
    # Initial zeroth calculation
    xnudrageff=[]
    xnudtot2=[]
    xnudrageff0=[]
    xnudrageff=[]
    mIOL=[]
    delv1New=[]
    debugDic={'term1':[],'term2':[],'term3':[],'term4':[]}
    for y11,y22,xni,xnC,torvhat,xnuc12,xnuc21,yy2,yy1 in zip(data['y11'].values(), 
                         data['y22'].values(), 
                         data['xni'].values(), 
                         data['xnC'].values(),
                         data['torvhat'].values(),
                         data['xnuc12'].values(),
                         data['xnuc21'].values(),
                         data['yy2'].values(),
                         data['yy1'].values()):
                             
        if switch==1:
            tempxnudrageff0=((y11+y22)/((xni*data['xmas1']+xnC*data['xmas2'])*torvhat))
            
            
            xnudrageff0.append(tempxnudrageff0)
            debugDic['term1'].append(xni)
            debugDic['term2'].append(y22)
            debugDic['term3'].append(xni*data['xmas1'])
            debugDic['term4'].append(torvhat)        
                                     
            tempDelv0=((y11-(xni*data['xmas1']*tempxnudrageff0*torvhat))/(xni*data['xmas1']*(xnuc12+tempxnudrageff0)))
            # What is Delv0 expanded out? Determine ordering
            delv0.append(tempDelv0)
            
            tempxnudtot2=(y22+(xnC*data['xmas2']*xnuc21*tempDelv0))/(xnC*data['xmas2']*torvhat)
            xnudtot2.append(tempxnudtot2)
        
            for n in range(1):
                tempxnudrageff=((((xni*data['xmas1'])+(xnC*data['xmas2']))*tempxnudrageff0)-
                     (xnC*data['xmas2']*tempxnudtot2))/(xni*data['xmas1'])
                delv1=(y11-xni*data['xmas1']*tempxnudrageff*torvhat)/(xni*data['xmas1']*(xnuc12+tempxnudrageff))
                 
                xnudrageff.append(tempxnudrageff)
                delv1New.append(delv1)
        
        # New method of calculating nudrag with Miol term and perturbation theory      
        elif switch==2:
            pass        
             
    torvDPert=[]    
    for delv1,torvhat in zip(delv1New,data['torvhat'].values()):
        temp=torvhat+delv1
        torvDPert.append(temp)
             
    return xnudrageff,torvDPert
    
if __name__=="__main__":
    pass