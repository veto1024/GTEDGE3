#!/usr/bin/python

##########################################################
#
#    Graphical presentation of chi analysis results
#    
# 
##########################################################

from lib.graphs.graphs import prettyCompare


class s000000():
    shotid=000000
#   Diffusion coefficient graph template    
    
    def __init__(self,data,tid):

        self.data=data
        self.tid=tid
        self.runDic={0000:self.runid0000}
        self.runDic[tid]() 
    def runid0000(self):
        diff=self.data.diff
        prettyID=[r'$D^{pinch}_j$',r'$D^{std}_{j}$']
        title=r"$\bf{D}_j$ Considering Non-diffusive Effects"
                  
        adjustments={}
        caption=[r"Main ion particle diffusion coefficient in the edge for",
                 r"DIII-D shot "+str(self.shotid)+"."+str(self.tid)+" considering IOL and pinch",
                 r"velocity, $D^{pinch}_{j}$, compared with standard",
                 r"definition, $D^{std}_j.$"]
        
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],diff[0]),(prettyID[1],diff[1])],
#                             yrange=[-.5,.5],
                             datalabels=[prettyID[0],prettyID[1]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
#                             textSpace=.025,
#                             marginBottom=.2,
#                             yLabelAdj=-.6,
#                             marginLeft=.15,
                             size=(6,9))
        print self.data.rhor          
        adjustments={}
        caption=["Main ion particle diffusion coefficient in the core for",
                 "DIII-D shot 000000.0000 considering IOL and pinch velocity"]
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],diff[0][:-30])],
                             corePlot=True,
#                             yrange=[-.5,.5],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
#                             textSpace=.02,
                             size=(6,9))

class s118888():
    shotid=118888
#   Diffusion coefficient for shot 118888  
    
    def __init__(self,data,tid):

        self.data=data
        self.diff=data.diff
        self.runDic={1525:self.runid1525,
                     1570:self.runid1570}
        self.runDic[tid]() 
    def runid1525(self):
        diff=self.diff
        prettyID=[r'$D^{pinch}_j$',r'$D^{std}_{j}$']
        title=r"$\bf{D}_j$ Considering Non-diffusive Effects"
                  
        adjustments={}
        caption=[r"Main ion particle diffusion coefficient in the edge for",
                 r"DIII-D shot 118888.1525 considering IOL and pinch",
                 r"velocity, $D^{pinch}_{j}$, compared with standard",
                 r"definition, $D^{std}_j.$"]
        
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],diff[0]),(prettyID[1],diff[1])],
                             yrange=[-.5,.5],
                             datalabels=[prettyID[0],prettyID[1]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.2,
                             yLabelAdj=-.7,
                             marginLeft=.15,
                             size=(6,9))
                  
        adjustments={}
        caption=["Main ion particle diffusion coefficient in the core for",
                 "DIII-D shot 118888.1525 considering IOL and pinch velocity"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],diff[0][:-30])],
                             corePlot=True,yrange=[-.5,.5],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02,
                             size=(6,9))
    def runid1570(self):
        diff=self.diff
        prettyID=[r'$D^{pinch}_j$',r'$D^{std}_{j}$']
        title=r"$\bf{D}_j$ Considering Non-diffusive Effects"
                  
        adjustments={}
        caption=[r"Main ion particle diffusion coefficient in the edge for",
                 r"DIII-D shot 118888.1570 considering IOL and pinch",
                 r"velocity, $D^{pinch}_{j}$, compared with standard",
                 r"definition, $D^{std}_j.$"]
        
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],diff[0]),(prettyID[1],diff[1])],
                             yrange=[-.5,.5],
                             datalabels=[prettyID[0],prettyID[1]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.2,
                             yLabelAdj=-.6,
                             marginLeft=.15,
                             size=(6,9))
                  
        adjustments={}
        caption=["Main ion particle diffusion coefficient in the core for",
                 "DIII-D shot 118888.1570 considering IOL and pinch velocity"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],diff[0][:-30])],
                             corePlot=True,yrange=[-.5,.5],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02,
                             size=(6,9))

class s118890():
    shotid=118890
#   Diffusion coefficient for shot 118890   
    
    def __init__(self,data,tid):

        self.data=data
        self.diff=data.diff
        self.runDic={1515:self.runid1515}
        self.runDic[tid]() 
    def runid1515(self):
        diff=self.diff
        prettyID=[r'$D^{pinch}_j$',r'$D^{std}_{j}$']
        title=r"$\bf{D}_j$ Considering Non-diffusive Effects"
                  
        adjustments={}
        caption=[r"Main ion particle diffusion coefficient in the edge for",
                 r"DIII-D shot 118890.1515 considering IOL and pinch",
                 r"velocity, $D^{pinch}_{j}$, compared with standard",
                 r"definition, $D^{std}_j.$"]
        
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],diff[0]),(prettyID[1],diff[1])],
                             yrange=[-.5,.5],
                             datalabels=[prettyID[0],prettyID[1]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.2,
                             yLabelAdj=-.6,
                             marginLeft=.15,
                             size=(6,9))
                  
        adjustments={}
        caption=["Main ion particle diffusion coefficient in the core for",
                 "DIII-D shot 118890.1515 considering IOL and pinch velocity"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],diff[0][:-30])],
                             corePlot=True,yrange=[-.5,.5],
                             datalabels=[prettyID[0]],
                             title=title,
                             ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02,
                             size=(6,9))



class s144977():
    shotid=144977
    def __init__(self,data,tid):

        self.data=data
        self.diff=data.diff
        self.runDic={3000:self.runid3000}
        self.runDic[tid]() 
    def runid3000(self):
        diff=self.diff
        prettyID=[r'$D^{pinch}_j$',r'$D^{std}_{j}$']
        title=r"$\bf{D}_j$ Considering Non-diffusive Effects"
                  
        adjustments={}
        caption=[r"Main ion particle diffusion coefficient in the edge for",
                 r"DIII-D shot 144977.3000 considering IOL and pinch",
                 r"velocity, $D^{pinch}_{j}$, compared with standard",
                 r"definition, $D^{std}_j.$"]
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],diff[0]),(prettyID[1],diff[1])],
                             yrange=[-.5,.5],
                             datalabels=[prettyID[0],prettyID[1]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.2,
                             yLabelAdj=-.6,
                             marginLeft=.15,
                             size=(6,9))
                  
        adjustments={}
        caption=["Main ion particle diffusion coefficient in the core for",
                 "DIII-D shot 144977.3000 considering IOL and pinch velocity"]
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],diff[0][:-30])],
                             corePlot=True,yrange=[-.5,.5],
                             datalabels=[prettyID[0]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],caption=caption,adjust=adjustments,textSpace=.02,size=(6,9))

     
                             
class s166606():
    shotid=166606
    def __init__(self,data,tid):

        self.data=data
        self.diff=data.diff
        self.runDic={1950:self.runid1950}
        self.runDic[tid]() 
    def runid1950(self):
        diff=self.diff
        prettyID=[r'$D^{pinch}_j$',r'$D^{std}_{j}$']
        title=r"$\bf{D}_j$ Considering Non-diffusive Effects"
                  
        adjustments={}
        caption=[r"Main ion particle diffusion coefficient in the edge for DIII-D shot",
                 r"166606.1950 considering IOL and pinch velocity, $D^{pinch}_{j}$,",
                 r"compared with standard definition, $D^{std}_j$."]
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],diff[0]),(prettyID[1],diff[1])],
                             yrange=[-.5,.5],
                             datalabels=[prettyID[0],prettyID[1]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.025,
                             marginBottom=.18,
                             yLabelAdj=-.6,
                             marginLeft=.15,
                             size=(6,9))


                  
        adjustments={}
        caption=["Main ion particle diffusion coefficient in the core for",
                 "DIII-D shot 166606.1950 considering IOL and pinch velocity"]
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],diff[0][:-30])],
                             corePlot=True,yrange=[-.5,.5],
                             datalabels=[prettyID[0]],
                             title=title,ylabel=[r"$D_j$",r"$\left[\frac{m^2}{s}\right]$"],caption=caption,adjust=adjustments,textSpace=.02,size=(6,9))

                             
