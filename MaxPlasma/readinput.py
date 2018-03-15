# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:10:25 2017

@author: max
"""
from math import ceil
import os
class readinput():
    '''Read input file toplasma.'''
    
    def __init__(self, fileCat):
        self.inputfile = 'MaxPlasma/'+fileCat['toplasma']
        self.readtoneut(self.inputfile)

    def readtoneut(self, inputfile):
        print 'self.inputfile = '+ self.inputfile
        #Treating these as class variables. Could also treat as instance variables.
        with open(self.inputfile, 'r') as infile:
            for line in infile:
                # Typical line: variable = value
                variable, value = line.split('=')
                variable = variable.strip()  # remove leading/traling blanks
                if variable == 'd3d_iter':
                    self.d3d_iter = int(value)
                if variable == 'a':
                    self.a = float(value)
                elif variable == 'B_phi_0':
                    self.B_phi_0 = float(value)
                elif variable == 'R0_a':
                    self.R0_a = float(value)
                elif variable == 'Z0':
                    self.Z0 = float(value)
                elif variable == 'xpt':
                    self.xpt = eval(value)
                elif variable == 'kappa_up':
                    self.kappa_up = float(value)
                elif variable == 'kappa_lo':
                    self.kappa_lo = float(value)
                elif variable == 'tri_up':
                    self.tri_up = float(value)
                elif variable == 'tri_lo':
                    self.tri_lo = float(value)
                elif variable == 'xmil':
                    self.xmil = int(value)
                elif variable == 'thetapts_approx':
                    thetapts_approx = int(value)
                    self.thetapts =  int(4 * ceil(float(thetapts_approx)/4))+1 
                elif variable == 'rmeshnum_p':
                    self.rmeshnum_p = int(value)
                elif variable == 'rpts':
                    self.rpts = int(value)                
                elif variable == 'd3d_input':
                    self.d3d_input = int(value)
                elif variable == 'ni0':
                    self.ni0 = float(value)
                elif variable == 'ni9':
                    self.ni9 = float(value)
                elif variable == 'ni_sep':
                    self.ni_sep = float(value)
                elif variable == 'nu_ni':
                    self.nu_ni = float(value)
                elif variable == 'ne0':
                    self.ne0 = float(value)
                elif variable == 'ne9':
                    self.ne9 = float(value)
                elif variable == 'ne_sep':
                    self.ne_sep = float(value)
                elif variable == 'nu_ne':
                    self.nu_ne = float(value)
                elif variable == 'Ti0':
                    self.Ti0 = float(value)
                elif variable == 'Ti9':
                    self.Ti9 = float(value)
                elif variable == 'Ti_sep':
                    self.Ti_sep = float(value)
                elif variable == 'nu_Ti':
                    self.nu_Ti = float(value)
                elif variable == 'Te0':
                    self.Te0 = float(value)
                elif variable == 'Te9':
                    self.Te9 = float(value)
                elif variable == 'Te_sep':
                    self.Te_sep = float(value)
                elif variable == 'nu_Te':
                    self.nu_Te = float(value)
                elif variable == 'j0':
                    self.j0 = float(value)
                elif variable == 'j_sep':
                    self.j_sep = float(value)
                elif variable == 'nu_j':
                    self.nu_j = float(value)
                elif variable == 's_k_up':
                    self.s_k_up = float(value)
                elif variable == 's_k_lo':
                    self.s_k_lo = float(value)
                elif variable == 'xtheta1':
                    self.xtheta1 = float(value)
                elif variable == 'xtheta2':
                    self.xtheta2 = float(value)
                elif variable == 'xtheta3':
                    self.xtheta3 = float(value)
                elif variable == 'xtheta4':
                    self.xtheta4 = float(value)
            infile.close()    

    def showparams(self):
        '''Tell my details.'''
        print('**PARAMETERS FOR SHOT \'{}\'.'.format(self.shotlabel))
        for key in vars(params).iteritems():
            if key[0][0]!='_' and key[0]!='line' and key[0]!='infile' and key[0]!='variable' and key[0]!='value':
                print ('{} = {}'.format(key[0],key[1]))
        return '**END OF PARAMETERS**'