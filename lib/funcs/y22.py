#!/usr/bin/python

##################################################
#
#   Calculation of total carbon momentum terms
#
#
###################################################

   
def y22(data):
    y22=[]
    ephia,xmas1,xmas2=data.ephia,data.xmas1,data.xmas2
    for xmomtor2,xnC,bthet,gamC,zbar2 in zip(data.xmomtor2,
                                                 data.xnC,
                                                 data.bthet,
                                                 data.gamC,
                                                 data.zbar2):
         y22temp=xmomtor2+zbar2*data.xk*(xnC*ephia+bthet*gamC)
         y22.append(y22temp)
    return y22
    
if __name__=="__main__":
    pass