#!/usr/bin/python

#########################################
#
#   GTEDGE V2
#  
#  Main script
#
#########################################

import matplotlib.pyplot as plt
import lib.beams as beams
import lib.funcs as funcs
import lib.graphs as graphs
from MaxPlasma import plasma
#from lib.funcs.dataGen import newScatPlot,smooth,multiScatPlot
from scipy.interpolate import UnivariateSpline
from numpy import interp, pi
from math import sqrt
import inspect,datetime,warnings,sys


def maxDebug(coreData):
    
    maxiandic={}
    maxiandic['rhor']=coreData.rhor
    maxiandic['fiol']=coreData.fiol
    maxiandic['eiol']=coreData.eiol
    maxiandic['miol']=coreData.miol
    maxiandic['xti']=coreData.xti
    maxiandic['xni']=coreData.xni
    maxiandic['xte']=coreData.xne
    maxiandic['fiolC']=coreData.fiolC
    maxiandic['eiolC']=coreData.eiolC
    maxiandic['miolC']=coreData.miolC
    maxiandic['xnC']=coreData.xnC
    maxiandic['xre']=coreData.xer
    maxList=maxiandic.keys()
    
    funcs.dataGen.csvDump("fioldebug.csv",maxiandic)


def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stderr.write(warnings.formatwarning(message, category, filename, lineno))
    
###############################################################################
#
#
#    MAIN PROGRAM 
#
#
###############################################################################

if __name__== "__main__":
    
    ###########################################################################
    #
    #    New GTEDGE 3 with full core data and interpretation
    #
    #   Data provided as 201-point CSV files with the following naming conv.:
    #
    #   Input/GT3Profs_shotid_timeid.csv                Profiles w/ Core
    #   Input/GT3Consts_shotid_timeid.csv               Constants
    #   Input/nudraginputs_shotid_timeid_2v2.csv        GTEDGE profiles
    #   Input/GT3NBIConsts_144977_3000                  NBI-Relevent Constants
    #
    ############################################################################

#    shotid=118888
#    timeid=1525
#    runid="r88"
#    nbRun=True
#    IOL=True
    
#    shotid=118888
#    timeid=1570
#    runid="r88"
#    nbRun=True
#    IOL=True
	
#    shotid=118890
#    timeid=1515
#    runid="r90"
#    nbRun=True
#    IOL=True
    
#    shotid=118890
#    timeid=1560
#    runid="r90"
#    nbRun=True	
#    IOL=True
#    
#    shotid=166606
#    timeid=1950
#    runid="j8099"
#    nbRun=True
#    IOL=True
    
#    shotid=144977
#    timeid=3000
#    runid="j3000"
#    nbRun=True
#    IOL=True
    
    try:
        shotid,timeid,runid
    except:
        raise RuntimeError("Shotid/Timeid/Runid not found or formatted incorrectly")
    
    errLog=open("GTEDGE3warn.log","a+")
    errLog.write("\n")
    errLog.write("Time: %s \n" % str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    errLog.write("%s.%s.%s  nbRun=%s IOL=%s \n" % (str(shotid),str(timeid),str(runid),str(nbRun),str(IOL)))
    
    sys.stderr=errLog
    warnings.simplefilter('always', UserWarning)
    warnings.simplefilter('always', RuntimeWarning)
    warnings.showwarning=customwarn

    
    fileNames=funcs.dataCat.CatalogueCall(shotid,timeid,runid)

    
    GTedgeProfsFile=fileNames['GTEDGEsupp']
    constsFile=fileNames['GT3Consts']
    coreProfsFile=fileNames['GT3Profs']
    nbiFile=fileNames['GT3NBI']
    
#   GTEDGE Data
    
    GTEDGEdata=funcs.dataGen.dataGrab(profs=GTedgeProfsFile,const=constsFile)
    coreDataDic=funcs.dataGen.dataGrab(profs=coreProfsFile,const=constsFile)
    


    ###########################################################################
    #
    #   Augment coreData to allow things such as
    #   coreData.xni = []
    #   coreData.xti = []
    #   coreData.aminor = #.#
    #
    #   coreData will include 1) the coreDataDic dictionary, which includes
    #   all the mathematical functions, and 2) the above augmenting data
    #
    #   This removes the need to do things like coreData['key'].values()
    #
    ###########################################################################
    
    coreData=funcs.dataGen.dataAug(coreDataDic)
    coreData.rhor=[0.+x*(1./(len(coreData.xte)-1)) for x in range(len(coreData.xte))]    
    #   Calculate carbon density from experimental data
    #   Set zbar2 = 6. if data not provided
    try:
        coreData.zbar2
    except:
        coreData.zbar2=[6.]*201
    coreData.xnC=funcs.physFuncs.xnCcalc(coreData.fz1,coreData.zbar2)                    # Possibly not legit    
    #      Set q=4. if data not provided
    #   Note: This is "qedge" legacy from GTEDGE and should just be q
    #   in reality
    try:
        coreData.q
    except:
        coreData.q=[4.]*201
    
    coreData.xnuc=funcs.physFuncs.xnucCalc(coreData)

    
 
    

    ###########################################################################
    #
    #   Remnants of GTEDGE
    #   Here, we give GTEtoGT3 the GT3 rho vector, GTEDGE rho vector, and
    #   the GTEDGE 25-point profile to create a GT3 profile by interpolating
    #   and then setting to 0 everything up to the edge.
    #
    ###########################################################################

    coreData.xnuati=funcs.dataGen.GTEtoGT3(coreData.rhor,GTEDGEdata['rhor'],GTEDGEdata['xnuati'].values())
    coreData.xnuioni=funcs.dataGen.GTEtoGT3(coreData.rhor,GTEDGEdata['rhor'],GTEDGEdata['xnuioni'].values())
#    funcs.dataGen.newScatPlot(coreData.rhor,coreData.xnuati)
 #   funcs.dataGen.newScatPlot(coreData.rhor,coreData.xnuioni)
    
    ###################################################################
    #
    #   Run MaxPlasma here
    #    Make sure your ER file is in MaxPlasma/inputs
    #
    ###################################################################

    myPlasma=plasma(str(shotid)+"_"+str(runid)+"_"+str(timeid)+" Shot")
    myPlasma.solve(fileNames,data=coreData) 

 
    ###################################################################
    #
    #   Display summary of shot being run
    #   Will error if shot/run IDs do not match between Max's code
    #   and this code
    #   Also takes in flag to say if NBeams is running
    #
    ###################################################################
 
    funcs.dataGen.runSummary(funcs.dataGen.shotID(coreProfsFile,fileNames['Erspl']),nbeams=nbRun,IOLflag=IOL)
    
    ###################################################################
    #
    #   IOL Calculations
    #
    ###################################################################
    maxrhor=[0.+x*(1./(len(myPlasma.tiol[1].F_orb)-1)) for x in range(len(myPlasma.tiol[1].F_orb))]
        
#    coreData.bthet=myPlasma
#    coreData.bphi=myPlasma
    
    ep = [a*coreData.aminor/coreData.rmajor for a in coreData.rhor]
    coreData.bthet=[a * abs(coreData.bphi) / coreData.q95 for a in ep]
    
    coreData.fiol=interp(coreData.rhor,maxrhor,myPlasma.tiol[1].F_orb)
    coreData.fiol=coreData.fiol.tolist()
    coreData.miol=interp(coreData.rhor,maxrhor,myPlasma.tiol[1].M_orb)
    coreData.miol=coreData.miol.tolist()
    coreData.eiol=interp(coreData.rhor,maxrhor,myPlasma.tiol[1].E_orb)
    coreData.eiol=coreData.eiol.tolist()
    
    coreData.fiolC=interp(coreData.rhor,maxrhor,myPlasma.tiol[6].F_orb)
    coreData.fiolC=coreData.fiolC.tolist()
    coreData.miolC=interp(coreData.rhor,maxrhor,myPlasma.tiol[6].M_orb)
    coreData.miolC=coreData.miolC.tolist()
    coreData.eiolC=interp(coreData.rhor,maxrhor,myPlasma.tiol[6].E_orb)
    coreData.eiolC=coreData.eiolC.tolist()
        
    if IOL==False:
        coreData.fiol,coreData.miol,coreData.eiol=[0.]*len(coreData.fiol),[0.]*len(coreData.miol),[0.]*len(coreData.eiol)
        coreData.fiolC,coreData.miolC,coreData.eiolC=[0.]*len(coreData.fiolC),[0.]*len(coreData.miolC),[0.]*len(coreData.eiolC)
    ###################################################################
    #
    #    Plug in GTNeut here
    #
    ###################################################################

    ###################################################################
    #
    #   NBeams call, set data to nbiDep attribute in coreData
    #   Data obtained as NBIresults are as follows:
    #
    #   NBIresults[0]=nbiDep list (0-2)
    #   NBIresults[1]=snbi distribution
    #   NBIresults[2]=qnbi distribution
    #   NBIresults[3]=qnbe distribution 
    #   NBIresults[4]=beam input toroidal momentum deterium
    #   NBIresults[5]=beam input toroidal momentum carbon
    ###################################################################    
    
    #   Need FIOL data from Max's plasma
    
    fastiolzip = [[0.]*201,[0.]*201,[0.]*201]
    NBIresults=beams.run(nbiFile,coreData.xte,coreData.xti,coreData.xne,coreData.xni,coreData.xnC,coreData.zbar2,coreData.rmajor,coreData.elong,fastiolzip,run=nbRun)
   
    coreData.Snbi=NBIresults[1]
    coreData.qnbi=NBIresults[2]
    coreData.qnbe=NBIresults[3]
    coreData.xmomtor1=NBIresults[4]
    coreData.xmomtor2=NBIresults[5]
    
    ###################################################################
    #
    #   CXR data to be generated here. Currently set zbar2=6.0
    #   
    #
    ###################################################################
    

    coreData.gamC=[0.]*len(coreData.rhor)
    
    gamma,qHeati,qHeate=funcs.coreflux.fluxCore(coreData,iolFlag=False)
    gammahat,qhatHeati,qhatHeate=funcs.coreflux.fluxCore(coreData,iolFlag=True)
   
    
    funcs.dataGen.dictAppend(coreData,coreData.rhor,(('gamma',gamma),('gammahat',gammahat),('qHeati',qHeati),('qhatHeati',qhatHeati),('qHeate',qHeate),('qhatHeate',qhatHeate),('gamC',coreData.gamC)))

    coreData.y11=funcs.y11.y11(coreData)
    coreData.y22=funcs.y22.y22(coreData) 
    
    ###################################################################
    #
    #    Intrinsic Rotation
    #
    ###################################################################
    
    rlossiol=0.5

    coreData.intrin=[2./sqrt(pi)*morbl*sqrt(2*coreData.xk*Temp/coreData.xmas1) for morbl,Temp in zip(coreData.miol,coreData.xti)]
    coreData.intrinC=[2./sqrt(pi)*morbl*sqrt(2*coreData.xk*Temp/coreData.xmas2) for morbl,Temp in zip(coreData.miolC,coreData.xti)]
    
    coreData.vtorChat=[a-b for a,b in zip(coreData.vtorC,coreData.intrinC)]

 
#    coreData.vtorChat=[a-b for a,b in zip(coreData.vtorC,coreData.yy2)]    
    
    ##################################################################
    #
    #   Calculate deuterium velocities with Carbon data or exp. data
    #
    ##################################################################
    
    nudragResults=funcs.nuDragMIOL.nuDragMIOL(6,coreData,y11=coreData.y11,y22=coreData.y22)
  
    if ((hasattr(coreData,'vtorD')==False) and (hasattr(coreData,'vtord')==False)):
        coreData.vtorD=nudragResults[2]
        coreData.vtorDhat=[a-b for a,b in zip(coreData.vtorD,coreData.intrin)]
    else:
        coreData.vtorDhat=[a-b for a,b in zip(coreData.vtorD,coreData.intrin)]
    
    ##################################################################
    #
    #   Remove bad points from crazy perturbation theory calculation
    #
    ##################################################################    
    
    try:
        for a in fileNames['vtorID']:
            funcs.dataGen.pointRemove(coreData.vtorD,a) 
            funcs.dataGen.pointRemove(coreData.vtorDhat,a)
    except:
        print "Failed to remove points"
        pass
    coreData.nudrag=nudragResults[0]


#   TODO: Assuming Vpolc = 0.4*V_pold until better answer figured out
    
    coreData.vpolD=[a/.4 for a in coreData.vpolC]
    
    
    
    chi=funcs.chi.chiClass(coreData)
    chihat=funcs.chi.chiorbClass(coreData)  

    coreData.diff=funcs.diff.diffCalcs(coreData)
###############################################################################
#
#    Chi inference from expermental data
#
#   Chicalc(data,conv15=False,conv25=False,visc=False,pressure=False)
#
#   Setting an argument to True subtracts away that argument from qcond
#
###############################################################################

    chi1=chi.chiiCalc(coreData) 
    chi2=chihat.chiiCalc(coreData) 
    chi3=chihat.chiiCalc(coreData,conv25=True) 
    chi4=chihat.chiiCalc(coreData,conv25=True,pressure=True) 
    chi5=chihat.chiiCalc(coreData,conv25=True,pressure=True,visc=True) 
    chiList=[chi1,chi2,chi3,chi4,chi5]
    
    


#    q1=chi.qtotal 
#    q2=chihat.qtotal 
    q3=chihat.conv25 
    q4=chihat.heatin 
    q5=chihat.heatvisc 
    q2=[a-b-c-d for a,b,c,d in zip(chihat.qtotal,chihat.conv25,chihat.heatin,chihat.heatvisc)]
    q1=[a-b for a,b in zip(chi.qtotal,q2)]
    qList=[q1,q2,q3,q4,q5]
#    
    
    chiGraphsDic={}
    diffGraphsDic={}
    nuGraphsDic={}
    
    for name,obj in inspect.getmembers(graphs.chiGraphs):
        if inspect.isclass(obj):
            chiGraphsDic[obj.shotid]=obj
            
    for name,obj in inspect.getmembers(graphs.diffGraphs):
        if inspect.isclass(obj):
            diffGraphsDic[obj.shotid]=obj
            
    for name,obj in inspect.getmembers(graphs.nuGraphs):
        if inspect.isclass(obj):
            nuGraphsDic[obj.shotid]=obj
            
    
    try:
        chiGraphsDic[shotid](coreData,chiList,qList,timeid)
    except IndexError:
        print "fuckballs"
        raise IndexError("Index Error in ChiGraphDic")
    except:
        graphs.popups.popup("Chi graphs not defined for this shot.timeid, using default",title="Warning")
        chiGraphsDic[0000](coreData,chiList,qList,0000)
        
    try:
        diffGraphsDic[shotid](coreData,timeid)
    except IndexError:
        raise IndexError("Index Error in diffGraphsDic")
    except:
        graphs.popups.popup("Diffusion graphs not defined for this shot.timeid, using default",title="Warning")
        diffGraphsDic[0000](coreData,0000)
        
    try:    
        nuGraphsDic[shotid](coreData,timeid)
    except IndexError:
        raise IndexError("Index Error in nuGraphsDic")
    except:
        graphs.popups.popup("Nu drag graphs not defined for this shot.timeid, using default",title="Warning")
        nuGraphsDic[0000](coreData,0000)



    #   Dump FIOL debug data
    maxDebug(coreData)

    funcs.dataGen.variableFlushtoCSV("log.csv",coreData)
    plt.show(block=True)