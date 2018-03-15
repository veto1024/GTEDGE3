# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:54:44 2017

@author: max
"""

class neutrals():
    '''neutrals calculation.'''
    def __init__(self,brnd):
        self.makeneutmesh(brnd)

    def makeneutmesh(self,brnd):
        '''make the mesh to give to gtneut'''
