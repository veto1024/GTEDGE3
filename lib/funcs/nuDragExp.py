##############################################
#
#   Calculations of nudrag using experimental data
#   
#
##############################################

def nuDragExp(switch,data,y11=False,y22=False):
    if y11==False: y11=data.y11
    if y22==False: y22=data.y22
    # Nudrag calculation with experimental data corrected for IR
    if switch==1:
        delV1,xnudtot1,xnudtot2=[],[],[]
        xmas1,xmas2=data.xmas1,data.xmas2
        
        for vtordhat,torvhat,y11,xni,y22,xnC,xnuc12,xnuc21 in zip(data.vtordhat,
                                              data.torvhat,y11,
                                              data.xni,y22,
                                              data.xnC,
                                              data.xnuc12,
                                              data.xnuc21):
                                                  
            delV1temp=vtordhat-torvhat
            delV1.append(delV1temp)
            
            xnudtot1temp=(y11-xni*xmas1*xnuc12*delV1temp)/(xni*xmas1*vtordhat)
            xnudtot1.append(xnudtot1temp)
            
            xnudtot2temp=(y22+xni*xmas2*xnuc21*delV1temp)/(xnC*xmas2*torvhat)
            xnudtot2.append(xnudtot2temp)
            
        return xnudtot1
          
    # Nudrag calculation with experimental data uncorrected for IR in nudrag
    if switch==2:  
        delV1,xnudtot1,xnudtot2=[],[],[]
        xmas1,xmas2=data.xmas1,data.xmas2
        
        for vtord,torv,y11,xni,y22,xnC,xnuc12,xnuc21 in zip(data.vtord,
                                              data.torv,y11,
                                              data.xni,y22,
                                              data.xnC,
                                              data.xnuc12,
                                              data.xnuc21):
                                                  
            delV1temp=vtord-torv
            delV1.append(delV1temp)
            
            xnudtot1temp=(y11-xni*xmas1*xnuc12*delV1temp)/(xni*xmas1*vtord)
            xnudtot1.append(xnudtot1temp)
            
            xnudtot2temp=(y22+xni*xmas2*xnuc21*delV1temp)/(xnC*xmas2*torv)
            xnudtot2.append(xnudtot2temp)
            
        return xnudtot1
        
if __name__=="__main__":
    pass