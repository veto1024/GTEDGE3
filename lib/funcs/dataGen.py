#!/usr/bin/python
#
############################################
#
#
#   Data Generation/Plotting routines
#
#
############################################

import numpy as np
import scipy as sp
import os
import re
import sys
import csv
from matplotlib import rc
import glob
import brewery
from brewery import ds
import matplotlib.pyplot as plt
import lib.graphs.popups as pup
#matplotlib.use('Qt4Agg', warn=False)

class DataFileError(Exception):
    def __init__(self,value):
        self.value=value
    def __str__(self):
        return repr(self.value)

class profile:
    
    def __init__(self,x,y):
        # Defines internal smoothing function
        self.x=x
        self.y=y
            
        self.smooth=smooth(self.x,self.y)
        
    def __call__(self,x):
        # Allows one to call a specific value from GTEDGE 25 values        
        
        for a,b in zip(self.x,self.y):
            if x==a: return float(b)
    def lscale(self):
        result=[]

        for a,b in zip(self.y,self.derivValues):
            result.append(a/b)
        return result
        
    def lscaleCall(self,x):
        for a,b in zip(self.x,self.lscale()):
            if x==a: return b

        
    def scatPlot(self):
        # Prints a scatter plot of the current function
        plt.scatter(self.x,self.y,color='red')
        plt.show()
        
    def derivative(self,x,smoothCon=False,knots=False):
        
        # Returns the value of the derivative using the default smoothing function
        # for arbitrary values of rho.
        if smoothCon!=False:
            return float(smooth(self.x,self.y,smoothCon=smoothCon).derivative(1)(x))
        elif knots!=False:
            return float(smooth(self.x,self.y,knots=knots).derivative(1)(x))
        else:
            return float(self.smooth.derivative(1)(x))
    
    def values(self):
        # Returns a list of all the values of the function at the GTEDGE 25 points
        return self.y
        
    def interDerivs(self,rho,knotList):
        
        # After smoothing, use this function to define the 25 points for values of
        # the derivative
        
        def __call__(self,x):
            # Returns a value of the derivative from fit for  GTEDGE 25 points
            return self.derivValues(x)
        from scipy import interpolate
        y=smooth(rho,[self.smooth.derivative(1)(a) for a in rho])
        self.derivSpline=interpolate.LSQUnivariateSpline(rho,[y(x) for x in rho],bbox=[rho[0],rho[-1]],t=knotList)
        self.derivValues=[self.derivSpline(x) for x in self.x]
        
    def units(self,units):
        
        # TODO: Include units

        
        self.units=units
            
    def derivGrab(self,rho):
           for a,b in zip(self.x,self.derivValues):
               if a==rho: return b
                   
class dataAug:
    def __init__(self,inputDic):
        self.dic=inputDic
        for k, v in inputDic.items():
            try:
                setattr(self,k,v.y)
            except:
                pass
            if type(v)==float:
                setattr(self,k,v)
            if type(v)==int:
                setattr(self,k,v)

        
        
    
def gtedgePull(gtedgeInput,profs,consts):
    
    # Automatically generates .csv file of pyinput.txt files
    # to provide GTEDGE input values
    
    os.chdir(gtedgeInput)
    pyinputs=glob.glob('pyinput*')
    pyinputsConst=[]
    pyinputsProfs=[]
    sourcesConst=[]
    sourcesProfs=[]
    tempVarsConst=[]
    tempVarsProfs=[]
    varsConst=[]
    varsProfs=[]
    varNamesConst=[]
    varNamesProfs=[]
    constHeaders=[]
    profsHeaders=[]
    
    firstLine=True
    for filename in pyinputs:
        with open(filename) as openfile:
            num_lines=sum(1 for line in openfile)
            openfile.seek(0)
            if num_lines==2:
                for line in openfile:
                    if firstLine==True:
                        line=line.strip()
                        varNamesConst.append(line.split(','))
                        constHeaders.append(line.split(','))
                        firstLine=False
                    else:
                        tempVarsConst.append(line.split())
                        varsConst.append(line.split())
                sourcesConst.append({"file":filename+".csv","fields":varNamesConst[0]})
            else:
                for line in openfile:
                    if firstLine==True:
                        line=line.strip()
                        varNamesProfs.append(line.split(','))
                        profsHeaders.append(line.split(','))
                        firstLine=False
                    else:
                        tempVarsProfs.append(line.split())
                        varsProfs.append(line.split())
                sourcesProfs.append({"file":filename+".csv","fields":varNamesProfs[0]})
            firstLine=True

        with open(filename+'.csv','wb') as csvfile:
            if num_lines==2:
                csvWriter=csv.writer(csvfile,delimiter=',')
                csvWriter.writerow(varNamesConst[0])
                for a in tempVarsConst:
                    csvWriter.writerow(a)
                varNamesConst=[]
                tempVarsConst=[]
            else:
                csvWriter=csv.writer(csvfile,delimiter=',')
                csvWriter.writerow(varNamesProfs[0])
                for a in tempVarsProfs:
                    csvWriter.writerow(a)
                varNamesProfs=[]
                tempVarsProfs=[]
                         
    all_fields_Const=brewery.FieldList(["file"])
    for source in sourcesConst:
        for field in source["fields"]:
            pass
#            if field not in all_fields_Const.fields():
#                all_fields_Const.append(field)
    breweryConstOut=ds.CSVDataTarget(gtedgeInput+"_const.csv",)
    breweryConstOut.fields=brewery.FieldList(source["fields"])
    breweryConstOut.initialize()
    for source in sourcesConst:
        path=source["file"]
        src=ds.CSVDataSource(path,read_header=True,skip_rows=0)
        src.fields=ds.FieldList(source["fields"])
        src.initialize()
        breweryConstOut.field_names=(source["fields"])
        for record in src.records():
            record["file"]=path
            breweryConstOut.append(record)
        src.finalize()
    breweryConstOut.finalize()
    breweryConstOut.close_file
    
    all_fields_Profs=brewery.FieldList(["file"])
    for source in sourcesProfs:
        for field in source["fields"]:
            if field not in all_fields_Const:
                all_fields_Profs.append(field)
    breweryProfsOut=ds.CSVDataTarget(gtedgeInput+"_profs_temp.csv")
    breweryProfsOut.fields=brewery.FieldList(all_fields_Profs)
    breweryProfsOut.initialize()
    for source in sourcesProfs:
        path=source["file"]
        src=ds.CSVDataSource(path,read_header=False,skip_rows=1)
        src.fields=ds.FieldList(source["fields"])
        src.initialize()
        for record in src.records():
            record["file"]=path
            breweryProfsOut.append(record)
        src.finalize()
    breweryProfsOut.finalize()
    breweryProfsOut.close_file
    
    in_file=open(gtedgeInput+"_profs_temp.csv")
    out_file=open(gtedgeInput+"_profs.csv","wb+")
    
    for line in in_file:
        line=re.sub('\,+',',',line)
        out_file.write(line)
    out_file.close()
        
    
#    profsHeaders=sum(profsHeaders,[])        
#    all_fields=brewery.FieldList(profsHeaders)
#    breweryProfsOut=brewery.ds.CSVDataTarget(gtedgeInput+"_profs.csv")
#    breweryProfsOut.fields=all_fields
#    breweryProfsOut.initialize()
#    for source in pyinputsProfs:
#        src=brewery.ds.CSVDataSource(source,read_header=False,skip_rows=0)
#        src.fields=brewery.ds.FieldList(profsHeaders)
#        src.initialize()
#        for record in varsProfs:
#            breweryProfsOut.append(record)
#        src.finalize()
#    breweryProfsOut.finalize()
#    breweryProfsOut.close_file

#    newConst=open(gtedgeInput+"_consts.csv","wb+")
#    headers=[]
#    for file in pyinputsCons:
#        f=open(file)
#        headers.append(f.read().split(,))
        
    
    os.chdir(os.pardir)
    
    
    file_path_profs=os.path.relpath(profs)
    file_path_const=os.path.relpath(consts)
    
    f_profs=open(file_path_profs,'r')
    f_consts=open(file_path_const,'r')
    
def dataGrab(profs=False,const=False):
    
    # Pulls profile and constant data from files and puts them in a variable
    # dictionary. Includes spline fits for profiles.
    
    from scipy import interpolate
    
    # File input declaration
    
    inputDic={}
    if profs!=False:
        file_path_profs=os.path.relpath(profs)
        f_profs=open(file_path_profs,'r')
    if const!=False:
        file_path_const=os.path.relpath(const)
        f_consts=open(file_path_const,'r')
    
    # CSV module read in and append to inputDic dictionary
    if profs!=False:
        reader=csv.DictReader(f_profs)
    
        for row in reader:
            for column,value in row.iteritems():
                inputDic.setdefault(column,[]).append(float(value))
    
    # Reads in constants into the dictionary
    if const!=False:
        reader_consts=csv.DictReader(f_consts)  
        for row in reader_consts:
            for column,value in row.iteritems():
                inputDic.setdefault(column,[]).append(float(value)) 
            
    # Generate spline data for profiles != rho_r or constants
    for x in inputDic.keys():
        if((len(inputDic[x])>1) and (x!='rhor')):
#            spline=smooth(inputDic['rhor'],inputDic[x])
            inputDic[x]=profile(inputDic['rhor'],inputDic[x])
        elif (x!='rhor'):
            inputDic[x]=inputDic[x][0]
    
    # Attach units to profiles
    # TO BE IMPLEMENTED
    
#    m=open("Inputs/units.txt")
#    for line in m:
#        line=line.strip()
#        newLine=line.split("=")
#        (a,b)=(newLine[0],newLine[1])
#        if a!='rhor':
#            try:        
#                inputDic[a].units(str(b))
#            except(KeyError):
#                pass
    
    return inputDic

def scatPlot(x,y,color='red'):
    c=color
    plt.scatter(x,y,c)
    plt.show()

def linPlot(x,y):

    plt.plot(x,y,color='blue')

    return 

def smooth(inx,iny,smoothCon=None,knots=None):
    from scipy import interpolate
    x=[inx[r] for r in range(len(inx))]
    y=[iny[s] for s in range(len(iny))]
#    for a,b in zip(x,y):
#        if abs(b)<=1e-19:
#            print "value removed!"
#            x.remove(a)
#            y.remove(b)
    # Remove values greater than 2 std away
#    std=np.std(y)
#    avg=np.average(y)
#    for a,b in zip(x,y):
#        if abs(avg-b)>2.*std:
#            x.remove(a)
#            y.remove(b)

    try:
        tempArray=np.column_stack((x,y))
    except:
#        print x
 #       print y
        raise "Crash!"
    for n in range(len(tempArray)):
        try: 
            if tempArray[n,1]==0.:
                tempArray=np.delete(tempArray,n,0)
        except:
            pass
#            print y
    if smoothCon!=None:
        a=[tempArray[n][0] for n in range(len(tempArray))]
        b=[tempArray[n][1] for n in range(len(tempArray))]
        spline=interpolate.UnivariateSpline(a,b,bbox=[tempArray[0][0],tempArray[-1][0]],s=smoothCon)
    elif knots!=None:
        a=[tempArray[n][0] for n in range(len(tempArray))]
        b=[tempArray[n][1] for n in range(len(tempArray))]
        spline=interpolate.LSQUnivariateSpline(a,b,bbox=[tempArray[0][0],tempArray[-1][0]],t=knots)
    else:
        a=[tempArray[n][0] for n in range(len(tempArray))]
        b=[tempArray[n][1] for n in range(len(tempArray))]
        spline=interpolate.UnivariateSpline(a,b,bbox=[tempArray[0][0],tempArray[-1][0]])
    return spline

def newScatPlot(x,y,xminmax=False,yminmax=False,xlabel=False,ylabel=False):
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    if xminmax!=False:
        ax1.set_xlim(left=xminmax[0],right=xminmax[1])
    if yminmax!=False:
        ax1.set_ylim(bottom=yminmax[0],top=yminmax[1])
    if xlabel!=False:
        ax1.set_xlabel(xlabel)
    if ylabel!=False:
        ax1.set_ylabel(ylabel)        
    ax1.scatter(x,y)
    fig.show()
    
def multiScatPlot(x,yArray,yLabel,legend=None,title="No Title",colors=None,fontsize=10,caption=None,yminmax=False):
    fig=plt.figure()    
    ax=fig.add_subplot(1,1,1)
    ax.set_xlabel('rho')
    ax.set_ylabel(yLabel)
    ax.set_title(title)
    if yminmax!=False:
        ax.set_ylim(bottom=yminmax[0],top=yminmax[1])
    if caption!=None:
        fig.text(0.1,.1,caption)
    try:
        if legend!=None:
            for y,color,legend in zip(yArray,colors,legend):
               ax.scatter(x,y,color=color,label=legend,s=40)
        else:
            for y,color in zip(yArray,colors):
               ax.scatter(x,y,color=color,s=40)            
    except:
        raise Warning("Number of arrays, colors, and legends do not match")
    plt.legend(scatterpoints=len(yArray),loc='upper left',ncol=1,fontsize=fontsize)
    fig.show()

def multiLinePlot(x,yArray,yLabel,legend,title="No Title",colors=None,fontsize=10,caption=None):
    fig=plt.figure()    
    ax=fig.add_subplot(1,1,1)
    ax.set_xlabel('rho')
    ax.set_ylabel(yLabel)
    ax.set_title(title)
#    if caption!=None:
#        fig.text(0.1,.1,caption)
    try:
        for y,color,legend in zip(yArray,colors,legend):
           ax.plot(x,y,color=color,label=legend)
    except:
        raise Warning("Number of arrays, colors, and legends do not match")
    plt.legend(scatterpoints=len(yArray),loc='upper left',ncol=1,fontsize=fontsize)
    
def fitDerivs(x,name,knotList=None):
    plt.clf()
    plt.scatter(x,[data[name].derivative(a) for a in x],color='red')
    y=smooth(x,[data[name].derivative(a) for a in x],knots=knotList)
    plt.plot(x,y(x),color='blue')
        
    
def fitPlots(name,sC=None,knotList=None):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.clear()
    ax.set_xlabel('rho')
    ax.set_ylabel(name)
    ax.set_title(name)
    ax.scatter(rho,data[name].values(),color='red')
    if sC!=None:
        y=smooth(rho,data[name].values(),smoothCon=sC)
        ax.plot(rho,y(rho))   
    if knotList!=None:
#        knotList.insert(0,rho[0])
#        knotList.append(rho[-1])
        y=smooth(rho,data[name].values(),knots=knotList)
        ax.plot(rho,y(rho))                         


def dicReview(inputDic):
    for x in inputDic.keys():
        if x!='rhor':
            try:
                newScatPlot(inputDic['rhor'],inputDic[x].values(),xlabel='rhor',ylabel=x)
            except:
                continue
    
    
def csvDump(file,varDic=False,listArray=False):
    import numpy as np
    if varDic!=False:
        varDicHeaders=varDic.keys()
        
#        try:
#            varDicHeaders.remove('rhor')
#        except:
#            pass
        
        tempNameList=np.array([varDicHeaders])
#        tempNameList=varDicHeaders
        tempList=[]
        for n in range(len(varDicHeaders)):
#            if n!='rhor':
            tempList.append(varDic[tempNameList[0][n]])
        tempArray=np.array(tempList)
        tempArray=tempArray.transpose()
        with open(file,"wb") as f:
            writer=csv.writer(f)
            writer.writerow(varDicHeaders)
            for n in range(len(varDic[tempNameList[0][0]])):
                writer.writerow(tempArray[n])
    if listArray!=False:
        headers=listArray[0]
        values=listArray[1]
        tempArray=np.array(values)
        tempArray=tempArray.transpose()
        with open(file,"wb") as f:
            writer=csv.writer(f)
            writer.writerow(headers)
            for n in range(len(values[0])):
                writer.writerow(tempArray[n])
        
def variableFlushtoCSV(file,data):
    
    names=dir(data)
    dumpdic={}
#    print names
    for word in names[:]:
        if word.startswith("__"):
            names.remove(word)
    for x in names:
        y = getattr(data,x)
        flag = type(y) is list
        if flag:
            dumpdic[x]=y

    csvDump(file,dumpdic)
    return
                
    
            
def scatAndLine(xin1,yin1,xin2,yin2,xminmax=False,yminmax=False,xlabel=False,ylabel=False):
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    if xminmax!=False:
        ax1.set_xlim(left=xminmax[0],right=xminmax[1])
    if yminmax!=False:
        ax1.set_ylim(bottom=yminmax[0],top=yminmax[1])
    if xlabel!=False:
        ax1.set_xlabel(xlabel)
    if ylabel!=False:
        ax1.set_ylabel(ylabel)        
    ax1.scatter(xin1,yin1)
    ax1.plot(xin2,yin2)

def dictAppend(data,rhor,valueTuples):
    for x in valueTuples:
        data.dic[x[0]]=profile(rhor,x[1])
        setattr(data,x[0],x[1])
    return
    

def GTEtoGT3(corerho,rho,prof):
    
    smoother=smooth(rho,prof)
    fullprof=[]
    for a in rho:
        if a>0.85:
            fullprof.append(float(smoother(a)))
    remLen=len(corerho)-len(rho)
    fullprof.reverse()
    fullprof.extend([0.]*remLen)
    fullprof.reverse()
    
    return fullprof
    
def pointRemove(a,num):
    a[num]=(a[num+1]+a[num-1])/2.
    
def shotID(runID,maxrunID):
    temp=re.split(r'\_',runID)
    temp[-1]=temp[-1].split('.')[0]
    templist1=[]
    
    for a in temp:
        try:
            templist1.append(int(a))
        except:
            continue
        
    temp=re.split(r'\_',maxrunID)
    temp[0]=temp[0][1:]
    templist2=[]
    for a in temp:
        try:
            templist2.append(int(a))
        except:
            continue
    if len(templist2)!=2:
        raise DataFileError("Input file names not formatted properly")
    for a,b in zip(templist1,templist2):
        if a!=b:
            raise DataFileError("WARNING: Input file shots/runids do not match")
            
    return templist1

def runSummary(idList,nbeams=None,IOLflag=None):

    string="Shot #: "+str(idList[0])+"\n"+ \
           "Run id: "+str(idList[1])+"\n"
    if nbeams==True:
        string=string+"NBeams Loaded \n"
    elif nbeams==False:
        string=string+"NBeams NOT Loaded \n"
    else: raise Warning("Nbeams switch not set")
    if IOLflag==True:
        string=string+"IOL being corrected \n"
    elif IOLflag==False:
        string=string+"IOL NOT being corrected \n"
    else: raise Warning("IOL switch not set")    
    pup.popup(msg=string,title="Run Summary")

    