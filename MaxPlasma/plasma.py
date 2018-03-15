# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 11:56:27 2017

@author: max
"""


from readinput import readinput
from getexp import getexp
from background import background
from thermaliol import thermaliol
from neutrals import neutrals
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.dont_write_bytecode = True



        
####################################################################
####################################################################

class plasma():
    '''Represents a teacher.'''
    def __init__(self, shotlabel):
        #Include all specs attributes in plasma class
        #readinput.__init__(self)
        #Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel

    def solve(self,fileCat,data=False):
        if data==False:
            print('Working on {}.'.format(self.shotlabel))
            print('Getting input parameters...')
            self.param = readinput(fileCat)
            print('Getting experimental data...')
            self.exp = getexp(fileCat)
            print('Creating background plasma...')
            self.brnd = background(self.param,self.exp)
            print('Calculating Ion Orbit Loss...')
            self.tiol = thermaliol(self.brnd)
            print('Calculating neutrals...')
            self.ntrl = neutrals(self.brnd)
            print('Plasma Solved')
        else:
            self.tiol={}
            print('Working on {} using experimental profiles.'.format(self.shotlabel))
            print('Getting input parameters...')
            self.param = readinput(fileCat)
            print('Getting sexy experimental data...')
            self.exp = getexp(fileCat)
            print('Creating background plasma...')
            print "theta points = "+str(self.param.thetapts)
            self.brnd = background(self.param,self.exp,data=data)
            print('Calculating Ion Orbit Loss...')
            self.tiol[1]=thermaliol(self.brnd,1)
            self.tiol[6]=thermaliol(self.brnd,6)
            print('Calculating neutrals...')
            self.ntrl = neutrals(self.brnd)
            print('Plasma Solved')
        return 

####################################################################
####################################################################
if __name__=="__main__":
    
    myshot = plasma('ITER Base Case')
    myshot.solve()
    
    mill_plot = plt.figure(1,figsize=(8,8))
    ax1 = mill_plot.add_subplot(111)
    ax1.axis('equal')
    CS = ax1.contourf(myshot.brnd.R,myshot.brnd.Z,myshot.brnd.B_p,500) #plot something calculated by miller
    plt.colorbar(CS)
    
    iol_plot1 = plt.figure(2,figsize=(8,8))
    ax1 = iol_plot1.add_subplot(111)
    ax1.set_title(r'$F_{orb}$',fontsize=20)
    ax1.set_ylabel(r'cumulative particle loss fraction',fontsize=20)
    ax1.set_xlabel(r'minor radius',fontsize=20)
    ax1.grid(b=True,which='both',axis='both')
    ax1.set_xlim(0.75, 1)
    ax1.plot(np.linspace(0,1,myshot.param.rpts),myshot.tiol.F_orb,label='Forb',color='black',lw=3)
    
    iol_plot2 = plt.figure(3,figsize=(8,8))
    ax1 = iol_plot2.add_subplot(111)
    ax1.set_title(r'$M_{orb}$',fontsize=20)
    ax1.set_ylabel(r'cumulative momentum loss fraction',fontsize=20)
    ax1.set_xlabel(r'minor radius',fontsize=20)
    ax1.grid(b=True,which='both',axis='both')
    ax1.set_xlim(0.75, 1)
    ax1.plot(np.linspace(0,1,myshot.param.rpts),myshot.tiol.M_orb,label='Morb',color='black',lw=3)
    
    iol_plot3 = plt.figure(4,figsize=(8,8))
    ax1 = iol_plot3.add_subplot(111)
    ax1.set_title(r'$E_{orb}$',fontsize=20)
    ax1.set_ylabel(r'cumulative energy loss fraction',fontsize=20)
    ax1.set_xlabel(r'minor radius',fontsize=20)
    ax1.grid(b=True,which='both',axis='both')
    ax1.set_xlim(0.75, 1)
    ax1.plot(np.linspace(0,1,myshot.param.rpts),myshot.tiol.E_orb,label='Eorb',color='black',lw=3)