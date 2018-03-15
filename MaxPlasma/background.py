# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:05:08 2017

@author: max
"""

from math import pi,ceil,sin,acos
import numpy as np
from scipy.constants import mu_0, elementary_charge, k
from helpers import PolyArea
from scipy.interpolate import UnivariateSpline,interp1d

class background():
    '''background plasma stuff.'''
    def __init__(self,p,experimental,data=False):
        self.createbackground(p,experimental,data)

    def createbackground(self,p,experimental,data):
        '''Create background plasma using the miller model.'''
        ## CREATE r AND theta MATRICES
        thetapts =  int(4 * ceil(float(p.thetapts)/4))+1
        r1d = np.linspace(0,p.a,p.rpts)
        theta1d = np.linspace(0,2*pi,thetapts)
        self.theta, self.r = np.meshgrid(theta1d,r1d)
        if data==False:
        ## CREATE DENSITY, TEMPERATURE, PRESSURE, AND CURRENT DENSITY MATRICES
            self.ni     = np.where(self.r<0.9*p.a,
                             (p.ni0-p.ni9)*(1-(self.r/p.a)**2)**p.nu_ni + p.ni9,
                             (p.ni_sep-p.ni9)/(0.1*p.a)*(self.r-0.9*p.a)+p.ni9)
    
            self.ne     = np.where(self.r<0.9*p.a,
                             (p.ne0-p.ne9)*(1-(self.r/p.a)**2)**p.nu_ne + p.ne9,
                             (p.ne_sep-p.ne9)/(0.1*p.a)*(self.r-0.9*p.a)+p.ne9)
    
            self.Ti_kev = np.where(self.r<0.9*p.a,
                             (p.Ti0-p.Ti9)*(1-(self.r/p.a)**2)**p.nu_Ti + p.Ti9,
                             (p.Ti_sep-p.Ti9)/(0.1*p.a)*(self.r-0.9*p.a)+p.Ti9)
            self.Ti_K = self.Ti_kev * 1.159E7
            self.Ti_J = self.Ti_kev * 1000 * elementary_charge
            
            self.Te_kev = np.where(self.r<0.9*p.a,
                             (p.Te0-p.Te9)*(1-(self.r/p.a)**2)**p.nu_Te + p.Te9,
                             (p.Te_sep-p.Te9)/(0.1*p.a)*(self.r-0.9*p.a)+p.Te9) 
            self.Te_K = self.Te_kev * 1.159E7
            self.Te_J = self.Te_kev * 1000 * elementary_charge
            
        else:
            
            ###################################################################
            #
            #   Interpolation keeps giving me garbage. Fuck you interpolation
            #
            ###################################################################

            thetaones=[1.]*thetapts
            self.ni={}
            ni_list=data.xni[::2]
            ni_list.pop()
            self.ni[1]=np.array([ni_list,]*thetapts).transpose()
            
            ne_list=data.xne[::2]
            ne_list.pop()
            self.ne=np.array([ne_list,]*thetapts).transpose()
            
            ti_list=data.xti[::2]
            ti_list.pop()
            ti_list=[a/1000. for a in ti_list]
            self.Ti_kev=np.array([ti_list,]*thetapts).transpose()
            
            self.Ti_K = self.Ti_kev * 1.159E7
            self.Ti_J = self.Ti_kev * 1000 * elementary_charge
            
            te_list=data.xte[::2]
            te_list.pop()
            te_list=[a/1000. for a in te_list]
            self.Te_kev=np.array([te_list,]*thetapts).transpose()
            
            
            self.Te_K = self.Te_kev * 1.159E7
            self.Te_J = self.Te_kev * 1000 * elementary_charge            

            ###################################################################
            #
            #   Carbon density
            #
            ###################################################################
            
            nC_list=data.xnC[::2]
            nC_list.pop()
            self.nC=np.array([nC_list,]*thetapts).transpose()
            

        self.pressure = [(a+b) * k * c for a,b,c in zip(self.nC,self.ni,self.Ti_K)]

        
        self.j_r = p.j0*(1-(self.r/p.a)**2)**p.nu_j
       
        ## CREATE kappa, tri AND RELATED MATRICES
        upperhalf   = (self.theta>=0)&(self.theta<pi)
        kappa  = np.where(upperhalf, 
                         p.kappa_up / (p.a**p.s_k_up) * self.r**p.s_k_up,
                         p.kappa_lo / (p.a**p.s_k_lo) * self.r**p.s_k_lo)
        
        tri_lo = sin(3*pi/2 - acos((p.xpt[0]-p.R0_a)/p.a))
        tri    = np.where(upperhalf,
                         p.tri_up * self.r/p.a,
                         tri_lo * self.r/p.a)

        s_tri  = np.where(upperhalf,
                         self.r*p.tri_up/(p.a*np.sqrt(1-tri)),
                         self.r*tri_lo/(p.a*np.sqrt(1-tri)))
        
        ## MODIFY kappa USING THE X-MILLER MODEL
        
        
        ## CALCULATE INITIAL R,Z WITH NO SHAFRANOV SHIFT
        ## (NECESSARY TO GET ESTIMATES OF L_r WHEN CALCULATING SHAFRANOV SHIFT)
        R0 = np.ones(self.r.shape) * p.R0_a 
        self.R = R0 + self.r * np.cos(self.theta+np.arcsin(tri)*np.sin(self.theta))
        self.Z = kappa*self.r*np.sin(self.theta)
        
        # THIS CALCULATES A MATRIX OF THE LENGTHS OF EACH SECTION OF EACH FLUX
        # SURFACE AND THEN SUMS THEM TO GET THE PERIMETER IN 2D OF EACH FLUX
        # SURFACE (VALUE OF r).
        L_seg = np.sqrt((self.Z-np.roll(self.Z,-1,axis=1))**2 + (self.R-np.roll(self.R,-1,axis=1))**2)
        L_seg [:,-1] = 0        
        L_r = np.tile(np.sum(L_seg, axis=1), (thetapts, 1)).T
        
        #CALCULATE CROSS-SECTIONAL AREA CORRESPONDING TO EACH r AND ASSOCIATED
        #DIFFERENTIAL AREAS
        area = np.zeros(self.r.shape)
        for i in range(0,p.rpts):
            area[i,:] = PolyArea(self.R[i,:],self.Z[i,:])
    
        diff_area = area - np.roll(area,1,axis=0)
        diff_area[0,:]=0
        
        diff_vol = diff_area * 2*pi*p.R0_a #approx because it uses R0_a instead of shifted R0
        vol = np.cumsum(diff_vol,axis=0)
        
        #Calculate each differential I and sum to get cumulative I
        j_r_ave = np.roll((self.j_r + np.roll(self.j_r,-1, axis=0))/2,1,axis=0)
        j_r_ave[0,:]=0
        diff_I = diff_area * j_r_ave
        I = np.cumsum(diff_I, axis=0)
        self.IP = I[-1,0]  
        
        #Calculate B_p_bar
        B_p_bar = mu_0 * I / L_r
        B_p_bar[0,:]=0
        
        #Calculate li
        li = (np.cumsum(B_p_bar**2 * diff_vol, axis=0) / vol) / (2*B_p_bar**2)
        li[0,:]=0
        
        #Calculate beta_p
        beta_p = 2*mu_0*(np.cumsum(self.pressure*diff_vol,axis=0)/vol-self.pressure) / B_p_bar**2
    
        #Calculate dR0dr
        dR0dr = np.zeros(self.r.shape)
        R0 = np.zeros(self.r.shape)
    
        f = 2*(kappa**2+1)/(3*kappa**2+1)*(beta_p+li/2)+1/2*(kappa**2-1)/(3*kappa**2+1)
        f[0,:] = f[1,:] ############ NEED TO REVISIT, SHOULD EXTRAPOLATE SOMEHOW
        
        dR0dr[-1,:] = -p.a*f[-1,:]/p.R0_a
        R0[-1,:] = p.R0_a
        
        for i in range(p.rpts-2,-1,-1):
            R0[i,:] = dR0dr[i+1,:] * (self.r[i,:]-self.r[i+1,:]) + R0[i+1,:]
            dR0dr[i,:] = -self.r[i,:]*f[i,:]/R0[i,:]
        
        #NOW USE UPDATED R0 AND dR0dr to get new R,Z.
        self.R = R0 + self.r * np.cos(self.theta+np.arcsin(tri)*np.sin(self.theta))
        self.Z = kappa*self.r*np.sin(self.theta)
        
        ## RECALCULATE GRAD-r
        dkappa_dtheta   = np.gradient(kappa, edge_order=1)[1] * thetapts/(2*pi)
        dkappa_dr       = np.gradient(kappa, edge_order=1)[0] * p.rpts/p.a
    
        dkappa_dtheta[-1] = dkappa_dtheta[-2]
        dkappa_dr[-1] = dkappa_dr[-2]
    
        dZ_dtheta       = self.r*(kappa*np.cos(self.theta)+dkappa_dtheta*np.sin(self.theta))
        dZ_dr           = np.sin(self.theta)*(self.r*dkappa_dr + kappa)
        dR_dr           = dR0dr - np.sin(self.theta + np.sin(self.theta)*np.arcsin(tri))*(np.sin(self.theta)*s_tri) + np.cos(self.theta+np.sin(self.theta)*np.arcsin(tri))
        dR_dtheta       = -self.r*np.sin(self.theta+np.sin(self.theta)*np.arcsin(tri))*(1+np.cos(self.theta)*np.arcsin(tri))
    
        abs_grad_r = np.sqrt(dZ_dtheta**2 + dR_dtheta**2) / np.abs(dR_dr*dZ_dtheta - dR_dtheta*dZ_dr)
        
        ## WE WANT TO CALCULATE THE POLOIDAL FIELD STRENGTH EVERYWHERE
        ## THE PROBLEM IS THAT WE'VE GOT 2 EQUATIONS IN 3 UNKNOWNS. HOWEVER, IF WE ASSUME THAT THE POLOIDAL
        ## INTEGRAL OF THE FLUX SURFACE AVERAGE OF THE POLOIDAL MAGNETIC FIELD IS APPROX. THE SAME AS THE
        ## POLOIDAL INTEGRAL OF THE ACTUAL POLOIDAL MAGNETIC FIELD, THEN WE CAN CALCULATE THE Q PROFILE
        self.B_phi = p.B_phi_0 * self.R[0,0] / self.R
        
        #Calculate initial crappy guess on q
        q = p.B_phi_0*self.R[0,0] / (2*pi*B_p_bar) * np.tile(np.sum(L_seg/self.R**2,axis=1), (thetapts, 1)).T #Equation 16 in the miller paper. The last term is how I'm doing a flux surface average
        q[0,:]=q[1,:]
        
        dPsidr = (p.B_phi_0 * self.R[0,0]) / (2*pi*q)*np.tile(np.sum(L_seg/(self.R*abs_grad_r),axis=1), (thetapts, 1)).T
        
        self.Psi = np.zeros(self.r.shape)
        for index,row in enumerate(self.r):
            if index >= 1:
                self.Psi[index] = dPsidr[index]*(self.r[index,0]-self.r[index-1,0]) + self.Psi[index-1]
        self.Psi_norm = self.Psi / self.Psi[-1,0]
        
        self.B_p = dPsidr * 1/self.R * abs_grad_r
        self.B_p[0,:] = 0
        
        
        self.B_phi = p.B_phi_0 * self.R[0,0] / self.R
        self.B_tot = np.sqrt(self.B_p**2 + self.B_phi**2)
        self.f_phi = self.B_phi/self.B_tot
        #######################################################################
        ## CALCULATE ELECTRIC POTENTIAL FROM EXPERIMENTAL RADIAL ELECTRIC FIELD DATA
        
        E_r_fit = UnivariateSpline(experimental.E_r[:,0], experimental.E_r[:,1])
        self.E_pot = np.zeros(self.r.shape)
        for i,rval in np.ndenumerate(self.r):
            self.E_pot[i[0],i[1]] = E_r_fit.integral(rval/p.a, 1.0)*1000.

        def xmiller(self):
            return
            
        return