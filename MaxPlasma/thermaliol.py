# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:48:34 2017

@author: max
"""

import numpy as np
from math import pi, sqrt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.constants import e
import matplotlib as mpl
import matplotlib.ticker as mtick
from scipy.stats import maxwell
from scipy.special import gamma, gammaincc
from scipy.interpolate import InterpolatedUnivariateSpline, griddata
from scipy.integrate import quad, quadrature
import sys
import math

AU=1.660539040*(10**-27)
mass = {1:2.01410178*AU,
        2:4.002602*AU,
        3:6.94*AU,
        4:9.0121831*AU,
        5:10.81*AU,
        6:12.011*AU,
        7:14.007*AU,
        8:15.999*AU,
        9:18.99840316*AU,
        10:20.1797*AU}


class thermaliol():
    '''background plasma stuff.'''
    def __init__(self,brnd,species):
        self.ni={1:brnd.ni,
                 6:brnd.nC}
        self.calctiol(brnd,species)

#   TODO: Bring in charge states

    def calctiol(self,brnd,species):
                
        '''Calculate thermal Ion Orbit Loss.'''
        numcos = 100
        coslist = np.linspace(-1,1,num=numcos)
        Tprofile = brnd.Ti_kev.T[0]
        
        polpts = len(brnd.r[-1])
        radpts = len(brnd.r.T[-1])
        thetas = brnd.theta[-1]
        thetanums = np.arange(polpts)
        rs = brnd.r.T[-1]
        
        #Launch point values
        r0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.r)[-1]
        theta0      = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.theta)[-1]
        R0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.R)[-1]
        Z0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.Z)[-1]
        ni0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],self.ni[species])[-1]
#        ni0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.ni)[-1]
        ne0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.ne)[-1]
        Ti0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.Ti_kev)[-1]
        Te0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.Te_kev)[-1]
        Bp0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.B_p)[-1]
        Bt0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.B_phi)[-1]
        B0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.B_tot)[-1]
        f0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.f_phi)[-1]
        Psi0        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.Psi)[-1]
        Psi_norm0   = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.Psi_norm) [-1] 
        phi0        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.E_pot)[-1]#*1000 #now in volts
        xi0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.linspace(-1,1,num=numcos)[:,None,None],np.ones(brnd.R.shape))[1]
        
       
        #Destination Point Values
        theta1      = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.theta[-1][:,None,None,None])[-1]
        R1          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.R[-1][:,None,None,None])[-1]
        f1          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.f_phi[-1][:,None,None,None])[-1]
        B1          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.B_tot[-1][:,None,None,None])[-1]
        Psi1        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.Psi[-1][:,None,None,None])[-1]
        phi1        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.E_pot[-1][:,None,None,None])[-1]#*1000 #now in volts
        
        a = (np.abs(R0/R1)*f0/f1*xi0)**2 - 1 + (1 - xi0**2)*np.abs(B1/B0)
        b = 2*species*e*(Psi0-Psi1)/(R1*mass[species]*f1) * (np.abs(R0/R1)*f0/f1*xi0)
        c = (species*e*(Psi0-Psi1)/(R1*mass[species]*f1))**2 - species*e*(phi0-phi1)/mass[species]
        
        v_sep_1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
        v_sep_2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
        
        v_sep = np.zeros(r0.shape)
        v_sep = np.where(
                np.logical_and(
                        np.logical_or(v_sep_1<=0,np.isnan(v_sep_1)),
                        np.logical_or(v_sep_2<=0,np.isnan(v_sep_2))
                        ),
                np.nan,v_sep)
        v_sep = np.where(
                np.logical_or(
                        np.logical_and(v_sep_1>0,np.logical_or(v_sep_2<=0,np.isnan(v_sep_2))),
                        np.logical_and(np.logical_and(v_sep_1>0,v_sep_2>0),v_sep_1<v_sep_2)
                        )
                ,v_sep_1,v_sep)
        v_sep = np.where(
                np.logical_or(
                        np.logical_and(v_sep_2>0,np.logical_or(v_sep_1<=0,np.isnan(v_sep_1))),
                        np.logical_and(np.logical_and(v_sep_1>0,v_sep_2>0),v_sep_2<=v_sep_1)
                        )
                ,v_sep_2,v_sep)
        
        ## PREP FOR FLOSS, MLOSS, AND ELOSS CALCUALTIONS
        v_sep_min = np.nanmin(np.nanmin(v_sep,axis=0),axis=2).T
        v_sep_min[-1] = 0
        T_matrix = np.zeros(v_sep_min.shape)
        for indx,row in enumerate(T_matrix):
            T_matrix[indx,:] = Tprofile[indx]
        zeta_matrix = np.zeros(v_sep_min.shape)
        for indx,column in enumerate(zeta_matrix.T):
            zeta_matrix[:,indx] = coslist[indx]
        
        eps_min = mass[species] * v_sep_min**2 / (2*T_matrix*1E3*1.6021E-19)
        
        ## F_orb calculation
        integrand = gammaincc(3/2,eps_min)
        print ("np.sum(integrand,axis=1) = ",np.sum(integrand,axis=1))
        print ("2/(numcos-1) = ",2./(numcos-1.))
        print ("2*gamma(3/2) = ",2.*gamma(3./2.))
        self.F_orb = np.sum(integrand,axis=1)*(2./(numcos-1.))/(2.*gamma(3./2.))
        print ("self.F_orb = ",self.F_orb)
        
        ## M_orb calculation
        integrand = zeta_matrix*gammaincc(2.,eps_min)
        self.M_orb = np.sum(integrand,axis=1)*(2./(numcos-1.))/(2.*gamma(2.))
        
        ## E_orb calculation
        integrand = gammaincc(5/2,eps_min)
        self.E_orb = np.sum(integrand,axis=1)*(2./(numcos-1.))/(5./2.*gamma(2.))
            
        return