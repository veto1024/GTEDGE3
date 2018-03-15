#!/usr/bin/python

##################################################
#
#   Calculation of total deuterium momentum terms
#
#
###################################################

   
def y11(data):
    y11=[]
    ephia,xmas1,xmas2=data.ephia,data.xmas1,data.xmas2
    for xmomtor1,xni,bthet,gammahat in zip(data.xmomtor1,
                                                 data.xni,
                                                 data.bthet,
                                                 data.gammahat):
                                                 
         y11temp=xmomtor1+data.xk*(xni*ephia+bthet*gammahat)
         y11.append(y11temp)
    return y11
    
if __name__=="__main__":
    pass