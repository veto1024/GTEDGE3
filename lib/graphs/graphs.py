#!/usr/bin/python

##########################################################
#
#    Graphical presentation of various analysis results
#    
# 
##########################################################

#from matplotlib import rc
import matplotlib.pyplot as plt
from math import ceil,log
import numpy as np
from lib.funcs.nuDragExp import *
from lib.funcs.nuDragMIOL import *
from lib.funcs.nuDragPert import *
from lib.funcs.chi import *
from lib.funcs.diff import *
from lib.funcs.y11 import *
from lib.funcs.y22 import *
import pandas as pd


#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

                  
def expDelv(data):
        # Experimental DelV
    expDelV=[]
    expDelVfluid=[]
    delViol=[]
    for vtordhat,torvhat,yy1,yy2 in zip(data['vtordhat'].values(),
                                        data['torvhat'].values(),
                                        data['yy1'].values(),
                                        data['yy2'].values()):
        expDelV.append((vtordhat+yy1)-(torvhat+yy2))
        expDelVfluid.append(vtordhat-torvhat)
        delViol.append(yy1-yy2)
#    caption=r'$\Delta V$ for various calculations using experimental Deuterium data'
    multiScatPlot(data['rhor'],[expDelV,expDelVfluid,delViol,nuDragMIOL(4,data)[2],nuDragMIOL(5,data)[2]],title=r"$V_D-V_C$",
                  legend=[r'$V^{tot,meas}_D-V^{tot,meas}_C$',
                          r'$V^{fluid,meas}_D-V^{fluid,meas}_C$',
                          r'$V^{iol}_D-V^{iol}_C$',
                          r'$V^{tot,pert}_D-V^{tot,meas}_C$ w/out IOL cor.',
                          r'$V^{fluid,pert}_D-V^{fluid,meas}_C w/ IOL corr.$'],
                  yLabel=r"$V_D-V_C$ (m/s)",
                  colors=['blue','red','black','green','orange'],fontsize=20)
                  
    
def gammaPlot(data,gamma,gammahat):
    
#   Plot gamma vs. gammahat
    
    plot=plt.figure()
    fig1=plot.add_subplot(111)
    fig1.set_xlabel(r'$\rho$',fontsize=20)
    fig1.set_ylabel(r'$\Gamma_{ion}$',fontsize=20)
    fig1.set_title('Radial particle flux comparison, edge')
    gamma=fig1.scatter(data.rhor,gamma,marker='o',color='red')
    gammahat=fig1.scatter(data.rhor,gammahat,marker='o',color='blue')
    fig1.legend((gamma,gammahat),(r'$\Gamma_{r}$',r'$\hat{\Gamma}_{r}$'),scatterpoints=1,loc='upper left',ncol=1,fontsize=18)
    
def gammaPlotFull(data):
#   Plot gamma vs. gammahat
    gamma=data.gamma
    gammahat=data.gammahat
    
    plot=plt.figure()
    fig1=plot.add_subplot(111)
    fig1.set_xlabel(r'$\rho$',fontsize=20)
    fig1.set_ylabel(r'$\Gamma_{ion}$',fontsize=20)
    fig1.set_title('Radial particle flux comparison, full plasma')
    gamma=fig1.scatter(data.rhor,gamma,marker='o',color='red')
    gammahat=fig1.scatter(data.rhor,gammahat,marker='o',color='blue')
    fig1.legend((gamma,gammahat),(r'$\Gamma_{r}$',r'$\hat{\Gamma}_{r}$'),scatterpoints=1,loc='upper left',ncol=1,fontsize=18)    
    plt.show()
    
    
def diffCoPlot(data,ndrag):
    
    plot=plt.figure()
    fig1=plot.add_subplot(111)
    fig1.set_xlabel(r'$\rho$',fontsize=20)
    fig1.set_ylabel(r'$D$',fontsize=20)
    fig1.set_title('Diffusion coefficient')
    NF=fig1.scatter(data['rhor'],diffCalcs(data,ndrag)[0],marker='o',color='red')
    VP=fig1.scatter(data['rhor'],diffCalcs(data,ndrag)[1],marker='o',color='blue')
#    fig1.set_ybound(lower=-1.,upper=2.)
    fig1.legend((NF,VP),(r'$D_{NF}$',r'$D_{VP}$'),scatterpoints=1,loc='upper right',ncol=1,fontsize=18)


def nudragCheck(data,ndrag):
    
    mag=[]
    for n,vphi,nudrag in zip(data['xni'].values(),
                             data['torvhat'].values(),
                             ndrag):
    
        magtemp=n*data['xmas1']*vphi*nudrag
        mag.append(magtemp)
                             
    plot=plt.figure()
    fig1=plot.add_subplot(111)
    fig1.set_xlabel(r'$\rho$',fontsize=20)
    fig1.set_ylabel(r'$n_jm_j\hat{V}_{\phi}\nu_{drag}$',fontsize=20)
    fig1.set_title("Total momentum input density")
    fig1.scatter(data['rhor'],mag,marker='o',color='blue')
    
    plot2=plt.figure()
    fig2=plot2.add_subplot(111)
    fig2.set_xlabel(r'$\rho$',fontsize=20)
    fig2.set_ylabel(r'$\nu_{drag}$',fontsize=20)
    fig2.set_title(r'$\nu_{drag}$ from perturbation theory')
    fig2.scatter(data['rhor'],ndrag,marker='o',color='red')

#    fig1.legend((nudrag,magplot,nuc12),(r'$\nu_{d}$',r'$\nu$',r'$\nu_{12}$'),scatterpoints=1,loc='upper right',ncol=1,fontsize=18)
#    fig1.legend((nudrag,magplot),(r'$\nu_{d}$',r'$\nu$'),scatterpoints=1,loc='upper right',ncol=1,fontsize=18)

def chiDebug(data,funcs):
    
    funcs.dataGen.newScatPlot(data.rhor,data.qhatHeati,ylabel="Qhat")
    funcs.dataGen.newScatPlot(data.rhor,data.chihat.heatin,ylabel="qhat_in")
    funcs.dataGen.newScatPlot(data.rhor,data.chihat.heatvisc,ylabel="qhat_visc")
    funcs.dataGen.newScatPlot(data.rhor,data.chihat.gamconvi,ylabel="qhat_conv")

    funcs.dataGen.newScatPlot(data.rhor,data.qHeati,ylabel="Q")
    funcs.dataGen.newScatPlot(data.rhor,data.chi.heatin,ylabel="q_in")
    funcs.dataGen.newScatPlot(data.rhor,data.chi.heatvisc,ylabel="q_visc")
    funcs.dataGen.newScatPlot(data.rhor,data.chi.gamconvi,ylabel="q_conv")
    
def prettyCompare(xdata,ydata,yrange=[0.,1.],ylabel="No Label",datalabels=False,
                  title="No Title",caption=False,floatCaption=None,
                  logSpace=False,corePlot=False,plotType="line",
                  adjust=None,textSpace=0.,size=(12,9),xTickremove=None,
                  yTickremove=None, marginLeft=None,marginRight=None,
                  marginBottom=None,marginTop=None,
                  xLabelAdj=0.,yLabelAdj=0., xlabel=r"$\rho$",
                  toSN=False):
        
    prettyData=[xdata]
    for a in ydata:
        prettyData.append(a)
    pandasFrame=pd.DataFrame.from_items(prettyData)

#    These are the "Tableau 20" colors as RGB.    

#    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
#                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
#                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
#                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
#                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  

    tableau20 = [(255,0,0), (0,0,255), (0,255,0), (255,0,255),    
                 (128,0,0), (128,0,128), (0,0,128), (0,128,128)]              
#    
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)      

    fig=plt.figure(figsize=size)
    
#   Remove the plot frame lines. They are unnecessary chartjunk.    
    ax = plt.subplot(111)    
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)    

    plt.subplots_adjust(bottom=marginBottom,left=marginLeft,top=marginTop,right=marginRight)
#   Ensure that the axis ticks only show up on the bottom and left of the plot.    
#   Ticks on the right and top of the plot are generally unnecessary chartjunk.    
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left()   
      
#   Limit the range of the plot to only where the data is.    
#   Avoid unnecessary whitespace.    
    plt.ylim(yrange[0],yrange[1])    
    if corePlot==True:
        plt.xlim(0,0.85+textSpace)     # core
        xrange=[0,.85]
    else:
        plt.xlim(.85,1+textSpace)    # edge
        xrange=[.85,1.]

#   Make sure your axis ticks are large enough to be easily read.    
#   You don't want your viewers squinting to read your plot. 
    if logSpace==False:
        if toSN!=False:
            if yTickremove!=None:
                plt.yticks(np.linspace(yrange[0],yrange[1],num=6)[:(-1*yTickremove)], [str(round(x,2)) for x in np.linspace(yrange[0]/(1.*10**toSN),yrange[1]/(1.*10**toSN),num=6)], fontsize=30)    
            else:
#                print np.linspace(yrange[0]/(10**toSN),yrange[1]/(1.*10**toSN),num=6)
                plt.yticks(np.linspace(yrange[0],yrange[1],num=6), [str(round(x,2)) for x in np.linspace(yrange[0]/(1.*10**toSN),yrange[1]/(1.*10**toSN),num=6)], fontsize=26)    
        elif yTickremove!=None:
            plt.yticks(np.linspace(yrange[0],yrange[1],num=6)[:(-1*yTickremove)], [str(round(x,2)) for x in np.linspace(yrange[0],yrange[1],num=6)], fontsize=30)    
        else:
            plt.yticks(np.linspace(yrange[0],yrange[1],num=6), [str(round(x,2)) for x in np.linspace(yrange[0],yrange[1],num=6)], fontsize=26)    
            
    else:
        if yTickremove!=None:
            plt.yticks([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)][:(-1*yTickremove)],[str(round(x,2)) for x in [log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]], fontsize=26)
        else:
            plt.yticks([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)],[str(round(x,2)) for x in [log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]], fontsize=26)
    
    plt.tick_params(axis='x',pad=18)



    if xTickremove!=None:
        plt.xticks(np.linspace(xrange[0],xrange[1],num=6)[:(-1*xTickremove)], [str(round(x,2)) for x in np.linspace(xrange[0],xrange[1],num=6)], fontsize=26)    
    else:
        plt.xticks(np.linspace(xrange[0],xrange[1],num=6), [str(round(x,2)) for x in np.linspace(xrange[0],xrange[1],num=6)], fontsize=26)    

        
#   Provide tick lines across the plot to help your viewers trace along    
#   the axis ticks. Make sure that the lines are light and small so they    
#   don't obscure the primary data lines. 
    if logSpace==False:        
        for y in np.linspace(yrange[0],yrange[1],num=6):    
            plt.plot(np.linspace(xrange[0],xrange[1],num=6), [y] * len(np.linspace(xrange[0],xrange[1],num=6)), "--", lw=0.5, color="black", alpha=0.3)     
    else:
        for y in [log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]:    
            plt.plot([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)], [y] * len([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]), "--", lw=0.5, color="black", alpha=0.3)     
        
#   Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
    plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                    labelbottom="on", left="off", right="off", labelleft="on")  
    ystep=(yrange[1]-yrange[0])/25.
    for rank, column in enumerate(datalabels):    
#       Plot each line separately with its own color, using the Tableau 20    
#       color set in order.    
        if plotType=="line":
            plt.plot(pandasFrame.rhor.values,
                     pandasFrame[column.replace("\n", " ")].values,    
                    lw=2.5, color=tableau20[rank])  
            y_pos=pandasFrame[column.replace("\n", " ")].values[-1]
        elif plotType=="scat":        
            plt.scatter(pandasFrame.rhor.values,
                    pandasFrame[column.replace("\n", " ")].values,    
                    color=tableau20[rank])                   
            y_pos=pandasFrame[column.replace("\n", " ")].values[-1]
        else:
            raise ("No line type requested")
        try:
            if rank in adjust.keys():
                y_pos += adjust[rank]*ystep

        except:
            pass
        
    # Again, make sure that all labels are large enough to be easily read    
    # by the viewer.    
        plt.text(xrange[1]+0.0025, y_pos, column, fontsize=20, color=tableau20[rank])  
    # matplotlib's title() call centers the title on the plot, but not the graph,    
    # so I used the text() call to customize where the title goes.    
      
    # Make the title big enough so it spans the entire plot, but don't make it    
    # so big that it requires two lines to show.    
      
    # Note that if the title is descriptive enough, it is unnecessary to include    
    # axis labels; they are self-evident, in this plot's case.    
#    plt.text((xrange[1]-xrange[0])/2.+xrange[0]+0.015, (yrange[1]-yrange[0])/20.+yrange[1], title, fontsize=26, ha="center")      
    plt.text((xrange[1]-xrange[0])/2.+xrange[0]+.5*textSpace, (yrange[1]-yrange[0])/20.+yrange[1], title, fontsize=26, ha="center")      

    if caption!=False:
        for num,cap in enumerate(caption):
            plt.text(xrange[0],yrange[0]-(yrange[1]-yrange[0])/7.5-(num)*(yrange[1]-yrange[0])/22.5,cap,fontsize=16,ha="left")
    if floatCaption!=None:
        plt.text(floatCaption[1][0],floatCaption[1][1],floatCaption[0],fontsize=22,ha="left",va="top")
        
    for num,label in enumerate(ylabel):
        plt.text(xrange[0]-(xrange[1]-xrange[0])/10.+yLabelAdj*(xrange[1]-xrange[0])/10.,(yrange[1]-yrange[0])/2.+yrange[0]-num*(yrange[1]-yrange[0])/10.,label,fontsize=38,ha="center")
            
    plt.text((xrange[1]-xrange[0])/2.+xrange[0],yrange[0]-(yrange[1]-yrange[0])/12.5+xLabelAdj*(yrange[1]-yrange[0]/10.),xlabel,fontsize=38,ha="center")
    