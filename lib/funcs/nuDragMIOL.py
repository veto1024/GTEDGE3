##############################################
#
#   Calculations of nudrag using MIOL 
#   
#
##############################################

    
# Nudrag calculation with new perturbation theory method (8/30/2016)
    

def nuDragMIOL(switch,data,y11=False,y22=False):
    if y11==False: y11=data['y11'].values()
    if y22==False: y22=data['y22'].values()

    xnudrag1,xnudrag10,delVpert1,pertTerm=[],[],[],[]
    debugDic={'term1':[],'term2':[],'term3':[],'term4':[]}
    if switch==1:
        for y11,y22,xni,xnC,torvhat,vtordhat,xnuc12,xnuc21,yy2,yy1 in zip(y11,y22, 
                             data['xni'].values(), 
                             data['xnC'].values(),
                             data['torvhat'].values(),
                             data['vtordhat'].values(),
                             data['xnuc12'].values(),
                             data['xnuc21'].values(),
                             data['yy2'].values(),
                             data['yy1'].values()):
            # temp remove IR
            torvhat=torvhat+yy2
            vtordhat=vtordhat+yy1
                    
            
            xnudrag10temp=y11+xni*data['xmas1']*xnuc12*(yy1-yy2)
            xnudrag10temp=xnudrag10temp/(xni*data['xmas1']*(torvhat-yy2))
            xnudrag10.append(xnudrag10temp)
            
            pertTermtemp=y11+xni*data['xmas1']*xnuc12*(yy1-yy2)-xni*data['xmas1']*xnudrag10temp*(torvhat-yy2)
            pertTermtemp=pertTermtemp/(xni*data['xmas1']*(xnudrag10temp+xnuc12))
            pertTerm.append(pertTermtemp)
            
    
            debugDic['term1'].append(y11)
            debugDic['term2'].append(xni*data['xmas1']*xnuc12*(yy1-yy2))
            debugDic['term3'].append(xni*data['xmas1']*(torvhat-yy2))        
            
            delVpert1temp=y11+xni*data['xmas1']*xnuc12*(yy1-yy2)-xni*data['xmas1']*xnudrag10temp*(torvhat-yy2)
            delVpert1temp=delVpert1temp/(xni*data['xmas1']*(xnudrag10temp+xnuc12))
            
            delVpert1.append(delVpert1temp)
            
            xnudrag1temp=y11-xni*data['xmas1']*((delVpert1temp)-(yy1-yy2))
            xnudrag1temp=xnudrag1temp/(xni*data['xmas1']*(torvhat+delVpert1temp-yy2))
            
            xnudrag1.append(xnudrag1temp)            
            debugDic['term4'].append(xnudrag10temp)
        
        
    elif switch==2:
        vtordhatpert=[]
        for y11,y22,xni,xnC,torvhat,vtordhat,xnuc12,xnuc21,yy2,yy1 in zip(y11,y22,
                             data['xni'].values(), 
                             data['xnC'].values(),
                             data['torvhat'].values(),
                             data['vtordhat'].values(),
                             data['xnuc12'].values(),
                             data['xnuc21'].values(),
                             data['yy2'].values(),
                             data['yy1'].values()):
        # temp remove IR
            torvhat=torvhat+yy2
            vtordhat=vtordhat+yy1
            
            xnudrag10temp=y11/(xni*data['xmas1']*(torvhat-yy2))
            xnudrag10.append(xnudrag10temp)
            
            pertTermtemp=y11-(data['xmas1']*xni*xnudrag10temp*(torvhat-yy2))
            pertTermtemp=pertTermtemp/((data['xmas1']*xni)*(xnuc12))
            pertTerm.append(pertTermtemp)
            
            vtordhatperttemp=pertTermtemp+yy1-yy2+torvhat
            vtordhatpert.append(vtordhatperttemp)
            
            xnudrag1temp=y11-(xni*data['xmas1']*xnuc12*pertTermtemp)
            xnudrag1temp=xnudrag1temp/(vtordhatperttemp-yy1)
            xnudrag1.append(xnudrag1temp)
            
    
        return xnudrag10,xnudrag1,vtordhatpert
        
    elif switch==3:
        vtordhatpert=[]
        xnudrageff0=[]
        delv0=[]
        xnudragk=[]
        xnudragj=[]
        db=[]
        for y11,y22,xni,xnC,torvhat,vtordhat,xnuc12,xnuc21,yy2,yy1 in zip(data['y11'].values(), 
                             data['y22'].values(), 
                             data['xni'].values(), 
                             data['xnC'].values(),
                             data['torvhat'].values(),
                             data['vtordhat'].values(),
                             data['xnuc12'].values(),
                             data['xnuc21'].values(),
                             data['yy2'].values(),
                             data['yy1'].values()):
        # temp remove IR
            torvhat=torvhat+yy2
            vtordhat=vtordhat+yy1        
            
            xnudrageff0temp=(y11+y22)/((torvhat-yy2)*(xni*data['xmas1']+xnC*(data['xmas2'])))
            xnudrageff0.append(xnudrageff0temp)
            dbtemp=yy2
            db.append(yy2)
            
            delv0temp=(y11-xni*data['xmas1']*xnudrageff0temp*(torvhat-yy2))/((xni*data['xmas1']*(xnudrageff0temp+xnuc12)))
            delv0.append(delv0temp)
            
            xnudragktemp=(y22+xnC*data['xmas2']*delv0temp)/(xnC*data['xmas2']*(torvhat-yy2))
            xnudragk.append(xnudragktemp)
    
            xnudragjtemp=(xnudrageff0temp*(xni*data['xmas1']+xnC*data['xmas2']))-xnC*data['xmas2']*xnudragktemp
            xnudragjtemp=xnudragjtemp/(xni*data['xmas1'])
            xnudragj.append(xnudragjtemp)
            
    
        return xnudrageff0,xnudragj
    if switch==4:     
        
        
        # 10/3/16 derivation with full velocities
        # DelV = yy1 - yy2
        vtordPert=[]
        xnudrag0=[]
        delv0=[]
        xnudragj1=[]
        vtordhatPert=[]
        for y11,y22,xni,xnC,torvhat,vtordhat,xnuc12,xnuc21,yy2,yy1 in zip(y11,y22,
                             data['xni'].values(), 
                             data['xnC'].values(),
                             data['torvhat'].values(),
                             data['vtordhat'].values(),
                             data['xnuc12'].values(),
                             data['xnuc21'].values(),
                             data['yy2'].values(),
                             data['yy1'].values()):
    
            torvhat=torvhat+yy2
            vtordhat=vtordhat+yy1  
    
    #            
            xnudrag0temp=y11+y22
            xnudrag0temp=xnudrag0temp/((xni*data['xmas1']+xnC*data['xmas2'])*torvhat+(xni*data['xmas1']*(yy1-yy2)))
            xnudrag0.append(xnudrag0temp)
    
            delv0temp=y11-xni*data['xmas1']*xnudrag0temp*torvhat
            delv0temp=delv0temp/(xni*data['xmas1']*(xnuc12+xnudrag0temp))
            delv0.append(delv0temp)
    
            # Deuterium perturbation velocity w/out IR correction
            vtordPerttemp=torvhat+delv0temp
            vtordhatPerttemp=torvhat+delv0temp-yy1
            
            vtordPert.append(vtordPerttemp)
            vtordhatPert.append(vtordhatPerttemp)
            
            xnudrag1temp=y11+y22
            xnudrag1temp=xnudrag1temp/((xni*data['xmas1']+xnC*data['xmas2'])*torvhat+(xni*data['xmas1']*delv0temp))  
            xnudrag1.append(xnudrag1temp)
        return   xnudrag0,xnudrag1,delv0,vtordPert,vtordhatPert
    
    if switch==5:     
        
        
        # 10/3/16 derivation with IOL corrected velocities
        # DelV = yy1 - yy2
        
        vtordPert=[]
        delv0=[]
        xnudragk=[]
        xnudragj1=[]
        vtordhatPert=[]
        xnudrag0=[]
        for y11,y22,xni,xnC,torvhat,vtordhat,xnuc12,xnuc21,yy2,yy1,vtord in zip(y11,y22,
                             data['xni'].values(), 
                             data['xnC'].values(),
                             data['torvhat'].values(),
                             data['vtordhat'].values(),
                             data['xnuc12'].values(),
                             data['xnuc21'].values(),
                             data['yy2'].values(),
                             data['yy1'].values(),
                             data['vtord'].values()):
      
    
    #            
            xnudrag0temp=y11+y22
            xnudrag0temp=xnudrag0temp/((xni*data['xmas1']+xnC*data['xmas2'])*torvhat+(xni*data['xmas1']*(yy1-yy2)))
            xnudrag0.append(xnudrag0temp)
    
            delv0temp=y11-xni*data['xmas1']*xnudrag0temp*torvhat
            delv0temp=delv0temp/(xni*data['xmas1']*(xnuc12+xnudrag0temp))
            delv0.append(delv0temp)
    
            # Deuterium perturbation velocity w/ IOL correction
            vtordPerttemp=torvhat+delv0temp+yy1
            vtordhatPerttemp=torvhat+delv0temp
            
            vtordPert.append(vtordPerttemp)
            vtordhatPert.append(vtordhatPerttemp)
            
            xnudrag1temp=y11+y22
            xnudrag1temp=xnudrag1temp/((xni*data['xmas1']+xnC*data['xmas2'])*torvhat+(xni*data['xmas1']*delv0temp))  
            xnudrag1.append(xnudrag1temp)
            
        
            
    if switch==6:     
        
        
        # Nick 12/2016 perturbation theory from new edgecalc
        # Equivalent to ndrag=8 I think
        
        vtordPert=[]
        xnudrageff1=[]
        delv1=[]
        xnudrageff2=[]
        xnudtot1=[]
        delv2=[]
        xnudrag1=[]
        xnudrag2=[]
        delv0=[]
        
        for y11,y22,xni,xnC,torvhat,xnuc12,xnuc21,yy2,yy1 in zip(y11,y22,
                             data.xni, 
                             data.xnC,
                             data.vtorChat,
                             data.xnuc[0,1],
                             data.xnuc[1,0],
                             data.intrinC,
                             data.intrin):
      
            torv=torvhat+yy2
            
            
            delv0temp=yy1-yy2
            delv0.append(delv0temp)
            xnudrageff1temp=(y11+y22)/((xni*data.xmas1+xnC*data.xmas2)*torv+xni*data.xmas1*delv0temp)
            xnudrageff1.append(xnudrageff1temp)

     
            delv1temp=(y11-xni*data.xmas1*xnudrageff1temp*torv)/(xni*data.xmas1*(xnuc12+xnudrageff1temp))
            delv1.append(delv1temp)
#            print xnuc12,xnudrageff1temp

            # PROBLEM: When nudrag is negative, it can come close to xnuc12 in magnitude
            #          and blow up the pertrubation theory.

            xnudrageff2temp=(y22+xnC*data.xmas2*xnuc21*delv1temp)/(xnC*data.xmas2*torv)
            xnudrageff2.append(xnudrageff2temp)
       
            xnudtot1temp=(y11+y22-xni*data.xmas1*xnudrageff1temp*delv1temp)/((xni*data.xmas1+xnC*data.xmas2)*torv)
            xnudtot1.append(xnudtot1temp)
            
            #delv1 in ndrag=8
            
            delv2temp=(y11-xni*data.xmas1*xnudtot1temp*torv)/(xni*data.xmas1*(xnuc12*xnudtot1temp))
            delv2.append(delv2temp)
            

            xnudrag1temp=(y11-xni*data.xmas1*xnuc12*delv2temp)/(xni*data.xmas1*(torv+delv2temp))
            xnudrag1.append(xnudrag1temp)
 
            
            xnudrag2temp=(y22+xnC*data.xmas2*xnuc21*delv2temp)/(xnC*data.xmas2*torv)
            xnudrag2.append(xnudrag2temp)
            
            vtordPerttemp=torv+delv1temp
            vtordPert.append(vtordPerttemp)
    
    
        return   xnudrag1,delv0,vtordPert  
        
    if switch==7:
        
        # Nudrag=9 from Nick's GTEDGE        
        
        delV1,xnudtot1,xnudtot2=[],[],[]
        xmas1,xmas2=data['xmas1'],data['xmas2']
        
        for vtord,torv,y11,xni,y22,xnC,xnuc12,xnuc21 in zip(data['vtord'].values(),
                                              data['torv'].values(),y11,
                                              data['xni'].values(),y22,
                                              data['xnC'].values(),
                                              data['xnuc12'].values(),
                                              data['xnuc21'].values()):
                                                  
            delV1temp=vtord-torv
            delV1.append(delV1temp)
            
            xnudtot1temp=(y11-xni*xmas1*xnuc12*delV1temp)/(xni*xmas1*vtord)
            xnudtot1.append(xnudtot1temp)
            
            xnudtot2temp=(y22+xni*xmas2*xnuc21*delV1temp)/(xnC*xmas2*torv)
            xnudtot2.append(xnudtot2temp)
            
        return delV1,xnudtot1,xnudtot2
            
            
            
        
if __name__=="__main__":
    pass