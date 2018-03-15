# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:54:54 2017

@author: max
"""

import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os

class getexp():
    '''background plasma stuff.'''
    def __init__(self,fileCat):
        
        self.readexpdata(fileCat)
    def readexpdata(self,fileCat):
        '''Read experimental data.'''
        ## Experimental data is in R,Z. R is in columns
        ## and Z is in rows so the data is arranged as it is physically
        #R matrix
  
        path=os.path.dirname(__file__)

        with open(path+"/inputs/"+fileCat['R']) as f:
            newlist = f.read().split()#.splitlines() 
            array1D = [float(i) for i in np.asarray(newlist)]
            
        self.R = np.zeros((65,65))
        for i in range(0,65):
            self.R[:,i]=array1D[i]
        
        #Z matrix
        with open(path+"/inputs/"+fileCat['Z']) as f:
            newlist = f.read().split()#.splitlines() 
            array1D = [float(i) for i in np.asarray(newlist)]
        self.Z = np.zeros((65,65))
        for i in range(0,65):
            self.Z[i,:]=array1D[i]
        #Z_exp = np.flipud(Z_exp)
        
        #Poloidal magnetic field
        with open(path+"/inputs/"+fileCat['bpol']) as f:
            newlist = f.read().split()#.splitlines() 
            array1D = [float(i) for i in np.asarray(newlist)]
        self.B_p = np.resize(array1D,(65,65))
        self.B_p_in = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
                                self.B_p.flatten(),
                                (self.R, self.Z),
                                method='cubic')
    
        #Toroidal magnetic Field
        with open(path+"/inputs/"+fileCat['btor']) as f:
            newlist = f.read().split()#.splitlines() 
            array1D = [float(i) for i in np.asarray(newlist)] 
        self.B_phi = np.resize(array1D,(65,65))
        
        #Seperatrix trace
        with open(path+"/inputs/"+fileCat['BDRY']) as f:
            newlist = f.read().split()#.splitlines() 
            array1D = [float(i) for i in np.asarray(newlist)]
        self.bdry = np.resize(array1D,(82,2)) #There are 82 sets of R,Z datapoints tracing out the LCFS
        

            
        #E_radial and electric potential
        E_r_dat = np.loadtxt(path+"/inputs/"+fileCat['Erspl'],skiprows=2) 
        E_r_x = np.around(E_r_dat[:,0],decimals=3)
        E_r_y = E_r_dat[:,2]
        self.E_r = np.column_stack((E_r_x,E_r_y))

            
    
        ####################################################################
        # IDENTIFY VERTICES IN LIM BASED ON THE ANGLES THEY CORRESPOND TO.
        # THIS PREVENTS UNNECESSARY POINTS FROM MESSING WITH THE MESH.
        ####################################################################
    
        #lim (first wall)
        with open(path+"/inputs/"+fileCat['lim']) as f:
            newlist = f.read().split()#.splitlines() 
            array1D = [float(i) for i in np.asarray(newlist)]
        lim = np.reshape(array1D, (-1, 2))    
    
        adotb = (lim[:,0]-np.roll(lim[:,0],1))*(lim[:,0]-np.roll(lim[:,0],-1)) + (lim[:,1]-np.roll(lim[:,1],1))*(lim[:,1]-np.roll(lim[:,1],-1))
        mag_a = np.sqrt((lim[:,0]-np.roll(lim[:,0],1))**2+(lim[:,1]-np.roll(lim[:,1],1))**2)
        mag_b = np.sqrt((lim[:,0]-np.roll(lim[:,0],-1))**2+(lim[:,1]-np.roll(lim[:,1],-1))**2)
        
        lim_angles = np.arccos(adotb/(mag_a*mag_b))/pi
        self.lim_vertex = np.zeros((lim.shape[0],2))
        
        for i in range(0,lim.shape[0]):
            if lim_angles[i] <= 0.9:
                self.lim_vertex[i,:]=lim[i,:]
            else:
                self.lim_vertex[i,:]=0
        self.lim_vertex = self.lim_vertex[np.all(self.lim_vertex != 0, axis=1)] #removing zeros from array
        #need to add in an additional criteria to also remove points that are extremely close to other points, even if they create a sufficiently large angle

            
        return