#!/usr/bin/python

##########################################################
#
#    Graphical presentation of momentum exchange frequency
#    analysis
#    
# 
##########################################################

from lib.graphs.graphs import prettyCompare


class s000000():
    shotid=000000
    def __init__(self,data,tid):

        self.data=data
        self.tid=tid
        self.runDic={0000:self.runid0000}
        self.runDic[tid]() 
    def runid0000(self):
        data=self.data
        prettyID=[r'$\nu_{d,j}$']
        title=r"Inferred $\nu_{d,j}$"
                  
        adjustments={}
        caption=[r"Inferred $\nu_{d,j}$ for DIII-D shot "+str(self.shotid)+"."+str(self.tid), 
                 "using pertrubation theory"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],data.nudrag)],
#                             yrange=[0.,2.E4],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$\nu_{d,j}$",r'$\left[s^{-1}\right]$'],
                             caption=caption,
                             adjust=adjustments,
#                             textSpace=.025,
#                             marginBottom=.18,
#                             yLabelAdj=-1.2,
#                             marginLeft=.2,
                             size=(6,9),
                             toSN=4)


class s118888():
    shotid=118888
    def __init__(self,data,tid):

        self.data=data
        self.runDic={1525:self.runid1525,}
        self.tid=tid
        self.runDic[tid]() 
    def runid1525(self):
        data=self.data
        prettyID=[r'$\nu_{d,j}$']
        title=r"Inferred $\nu_{d,j}$"
                  
        adjustments={}
        caption=[r"Inferred $\nu_{d,j}$ for DIII-D shot "+str(self.shotid)+"."+str(self.tid), 
                 "using pertrubation theory"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],data.nudrag)],
                             yrange=[0.,2.E4],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$\nu_{d,j}$",r'$\left[s^{-1}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.18,
                             yLabelAdj=-1.2,
                             marginLeft=.2,
                             size=(6,9),
                             toSN=4)
    def runid1570(self):
        data=self.data
        prettyID=[r'$\nu_{d,j}$']
        title=r"Inferred $\nu_{d,j}$"
                  
        adjustments={}
        caption=[r"Inferred $\nu_{d,j}$ for DIII-D shot "+str(self.shotid)+"."+str(self.tid), 
                 "using pertrubation theory"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],data.nudrag)],
                             yrange=[0.,2.E4],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$\nu_{d,j}$",r'$\left[s^{-1}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.18,
                             yLabelAdj=-1.2,
                             marginLeft=.2,
                             size=(6,9),
                             toSN=4)

class s118890():
    shotid=118890
    def __init__(self,data,tid):

        self.data=data
        self.runDic={1515:self.runid1515,}
        self.tid=tid
        self.runDic[tid]() 
    def runid1515(self):
        data=self.data
        prettyID=[r'$\nu_{d,j}$']
        title=r"Inferred $\nu_{d,j}$"
                  
        adjustments={}
        caption=[r"Inferred $\nu_{d,j}$ for DIII-D shot "+str(self.shotid)+"."+str(self.tid), 
                 "using pertrubation theory"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],data.nudrag)],
                             yrange=[0.,2.E4],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$\nu_{d,j}$",r'$\left[s^{-1}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.18,
                             yLabelAdj=-1.2,
                             marginLeft=.2,
                             size=(6,9),
                             toSN=4)    
class s144977():
    shotid=144977
    def __init__(self,data,tid):

        self.data=data
        self.runDic={3000:self.runid3000}
        self.runDic[tid]() 
    def runid3000(self):
        data=self.data
        prettyID=[r'$\nu_{d,j}$']
        title=r"Inferred $\nu_{d,j}$"
                  
        adjustments={}
        caption=[r"Inferred $\nu_{d,j}$ for DIII-D shot 144977.3000", 
                 "using pertrubation theory"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],data.nudrag)],
                             yrange=[0.,2.E4],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$\nu_{d,j}$",r'$\left[s^{-1}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.18,
                             yLabelAdj=-1.2,
                             marginLeft=.2,
                             size=(6,9),
                             toSN=4)
                             
class s166606():
    shotid=166606
    def __init__(self,data,tid):

        self.data=data
        self.runDic={1950:self.runid1950}
        self.runDic[tid]() 
    def runid1950(self):
        data=self.data
        prettyID=[r'$\nu_{d,j}$']
        title=r"Inferred $\nu_{d,j}$"
                  
        adjustments={}
        caption=[r"Inferred $\nu_{d,j}$ for DIII-D shot 166606.1950",
                 "using pertrubation theory"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],data.nudrag)],
                             yrange=[0,4.E3],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$\nu_{d,j}$",r'$\left[s^{-1}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.18,
                             yLabelAdj=-.6,
                             marginLeft=.15,
                             size=(6,9),
                             toSN=3)