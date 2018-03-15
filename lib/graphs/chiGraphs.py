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
    #   Template for chigraph classes
    
    def __init__(self,data,chis,qs,tid):

        self.data=data
        self.tid=tid
        self.chis=chis
        self.qs=qs
        self.runDic={0000:self.runid0000}
        self.runDic[tid]() 
    def runid0000(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$']
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={}       
        caption=["Main ion heat diffusivities in the edge for DIII-D shot "+str(self.shotid)+"."+str(self.tid)+" correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]   
        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2]),(prettyID[3],chis[3]),(prettyID[4],chis[4])],
#                             yrange=[-4.,1.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
#                             textSpace=.095,
#                             marginBottom=.15,
#                             marginLeft=0.075,
                             size=(16,12))
        print self.data.rhor          
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 000000.0000 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,
#                             yrange=[-1.,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.22)
                             
                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']        
        caption=["Heat flux terms for DIII-D shot 000000.0000 for the edge plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={
                }
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
#                             yrange=[-4.E4,4.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 000000.0000 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,
#                             yrange=[-1.E5,1.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.1)      
        
class s118888():
    shotid=118888
    #   Chi graphs for shot 118888
    
    def __init__(self,data,chis,qs,tid):

        self.data=data
        self.chis=chis
        self.qs=qs
        rundic={1525:self.runid1525,
                1570:self.runid1570}
        rundic[tid]()
    def runid1525(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$']
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={1:1.4,2:2.8,3:-1.4}       
        caption=["Main ion heat diffusivities in the edge for DIII-D shot 118888.1525 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]   
        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2]),(prettyID[3],chis[3]),(prettyID[4],chis[4])],
                             yrange=[0,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'',r'$\left[\frac{m^2}{s}\right]$'],
                                 caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
                             textSpace=.095,
                             marginBottom=.15,
                             marginLeft=0.1,
                             yLabelAdj=-0.5,
                             size=(16,12))
                  
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 118888.1525 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,
#                             yrange=[-1.,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.22)

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118888.1525 for the edge plasma" ]
        title=r"$\bf{Q}$ term comparison"
        adjustments={
                }
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
                             yrange=[-.5E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118888.1525 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,
#                             yrange=[-1.E5,1.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.1)

    def runid1570(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$']
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={1:1.4,2:2.8,3:-1.4}      
        caption=["Main ion heat diffusivities in the edge for DIII-D shot 118888.1570 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]   
        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2]),(prettyID[3],chis[3]),(prettyID[4],chis[4])],
                             yrange=[0.,2.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'',r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
                             textSpace=.095,
                             marginBottom=.15,
                             marginLeft=0.1,
                             yLabelAdj=-0.5,
                             size=(16,12))
                  
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 118888.1570 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,
                             yrange=[0.,2.5],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.22)

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118888.1570 for the edge plasma" ]
        title=r"$\bf{Q}$ term comparison"
        adjustments={
                }
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
                             yrange=[-0.5E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118888.1570 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,
                             yrange=[-1.2E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.1)

        
class s118890():
    shotid=118890
    #   Chi graphs for shot 118890
    
    def __init__(self,data,chis,qs,tid):

        self.data=data
        self.chis=chis
        self.qs=qs
        rundic={1515:self.runid1515,
                1560:self.runid1560}
        rundic[tid]()
    def runid1515(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$']
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={1:1.4,2:2.8,3:-1.4}       
        caption=["Main ion heat diffusivities in the edge for DIII-D shot 118890.1515 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]   
        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2]),(prettyID[3],chis[3]),(prettyID[4],chis[4])],
                             yrange=[0,2.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                                 caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
                             textSpace=.095,
                             marginBottom=.15,
                             marginLeft=0.1,
                             size=(16,12))
                  
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 118890.1515 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,
#                             yrange=[-1.,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.22)

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118890.1515 for the edge plasma" ]
        title=r"$\bf{Q}$ term comparison"
        adjustments={
                }
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
                             yrange=[-.5E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118890.1515 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,
#                             yrange=[-1.E5,1.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.1)

    def runid1560(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$']#,
#                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$',
#                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$']
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={1:1.4}#,2:2.8,3:-1.4}      
        caption=["Main ion heat diffusivities in the edge for DIII-D shot 118890.1560 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]   
        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2])],#(prettyID[3],chis[3]),(prettyID[4],chis[4])],
                             yrange=[0.,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2]],#,prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'',r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
                             textSpace=.045,
                             marginBottom=.15,
                             marginLeft=0.1,
                             size=(16,12))
                  
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 118890.1560 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,
                             yrange=[0.,2.5],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.22)

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118890.1560 for the edge plasma" ]
        title=r"$\bf{Q}$ term comparison"
        adjustments={
                }
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
                             yrange=[-0.5E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        
        caption=["Heat flux terms for DIII-D shot 118890.1560 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,
                             yrange=[-1.2E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],
                             caption=caption,
                             adjust=adjustments,
                             textSpace=.1)

    
class s144977():
    shotid=144977
    def __init__(self,data,chis,qs,tid):

        self.data=data
        self.chis=chis
        self.qs=qs
        self.runDic={3000:self.runid3000}
        self.runDic[tid]()

    def runid3000(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$',
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$']
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={1:1.4,4:1.3,3:-1.3}
        caption=["Main ion heat diffusivities in the edge for DIII-D shot 144977.3000 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]        
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2]),(prettyID[3],chis[3]),(prettyID[4],chis[4])],
                             yrange=[-4.,1.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
                             textSpace=.095,
                             marginBottom=.15,
                             marginLeft=0.075,
                             size=(16,12))
                  
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 144977.3000 correcting for IOL, convective heat flux",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,yrange=[-1.,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],caption=caption,adjust=adjustments,textSpace=.22)

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        caption=["Heat flux terms for DIII-D shot 144977.3000 for the edge plasma" ]
        title=r"$\bf{Q}$ term comparison"
        adjustments={1:.95,3:.75}
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
                             yrange=[-4.E4,4.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],caption=caption,adjust=adjustments,textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        caption=["Heat flux terms for DIII-D shot 144977.3000 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,yrange=[-1.E5,1.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],caption=caption,adjust=adjustments,textSpace=.1)      
                             
class s166606():
    shotid=166606
    def __init__(self,data,chis,qs,tid):

        self.data=data
        self.chis=chis
        self.qs=qs
        self.runDic={1950:self.runid1950}
        self.runDic[tid]()        
    def runid1950(self):
        chis,qs=self.chis,self.qs
        prettyID=[r'$Q^{cond}_j=Q^{total}_{j}$',   # 0 
                  r'$Q^{cond}_j=Q^{total}_j$ w/ IOL',  # 1 
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j$',  # 2
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$', #3
                  r'$Q^{cond}_j=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$'] #4
        title=r"Comparison of Ion Heat Diffusivity Considering Various Non-diffusive Mechanisms"
                  
        adjustments={0:.6,2:-.60,3:-1.9,4:.25}
        caption=["Main ion heat diffusivities in the edge for DIII-D shot 166606.1950 correcting for IOL, convective heat flux,",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
#        caption='''
#            Main ion heat diffusivities in the edge
#            for DIII-D shot # 166606.1950 correcting
#            for IOL,convective heat flux, work done
#            on the pressure tensor, and viscous 
#            heating (assuming 10% asymmetry).'''
        prettyCompare(('rhor',self.data.rhor ),[(prettyID[0],chis[0]),(prettyID[1],chis[1]),(prettyID[2],chis[2]),(prettyID[3],chis[3]),(prettyID[4],chis[4])],
                             yrange=[-.4,.1],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,
                             ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],
                             caption=caption,
                             adjust=adjustments,
                             xTickremove=None,
                             textSpace=.095,
                             marginBottom=.15,
                             marginLeft=0.075,
                             size=(16,12))


                  
        adjustments={}
        caption=["Main ion heat diffusivities in the core for DIII-D shot 166606.1950 correcting for IOL, convective heat flux,",
                "work done on the pressure tensor, and viscous heating (assuming 10% asymmetry)"]
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],chis[0][:-30]),(prettyID[1],chis[1][:-30]),(prettyID[2],chis[2][:-30]),(prettyID[3],chis[3][:-30]),(prettyID[4],chis[4][:-30])],
                             corePlot=True,yrange=[-1.,5.],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,ylabel=[r"$\chi_j$",r'$\left[\frac{m^2}{s}\right]$'],caption=caption,adjust=adjustments,textSpace=.22)

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        caption=["Heat flux terms for DIII-D shot 166606.1950 for the edge plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={1:.95,3:.75}
        prettyCompare(('rhor',self.data.rhor),[(prettyID[0],qs[0]),(prettyID[1],qs[1]),(prettyID[2],qs[2]),(prettyID[3],qs[3]),(prettyID[4],qs[4])],
                             yrange=[-2.5E4,2.5E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],caption=caption,adjust=adjustments,textSpace=.02)     
                             

                             
        prettyID=[r'$Q^{total}_j$',
                  r'$Q^{total}_j$ w/ IOL',
                  r'$Q^{conv15}_j$',
                  r'$Q^{heatin}_j$',
                  r'$Q^{visc}_j$']
        caption=["Heat flux terms for DIII-D shot 166606.1950 for the core plasma"]
        title=r"$\bf{Q}$ term comparison"
        adjustments={}
        prettyCompare(('rhor',self.data.rhor[:-30]),[(prettyID[0],qs[0][:-30]),(prettyID[1],qs[1][:-30]),(prettyID[2],qs[2][:-30]),(prettyID[3],qs[3][:-30]),(prettyID[4],qs[4][:-30])],
                             corePlot=True,yrange=[-1.E5,1.E4],
                             datalabels=[prettyID[0],prettyID[1],prettyID[2],prettyID[3],prettyID[4]],
                             title=title,ylabel=[r"$q$",r"$\left[\frac{W}{m^2}\right]$"],caption=caption,adjust=adjustments,textSpace=.1)  