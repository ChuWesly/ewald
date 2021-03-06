# -*- coding: utf-8 -*-
"""
Calculate the electromagnetic potential in a periodic box. Fisrt we will show 
The limitation width the traditional approach and then we will implement a 
ewals summation.
Created on Sat Mar  6 16:35:36 2021

@author: Wesly
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import sys

dx = 0.5 # [nm] stepsize x of the periodical box
dy = 0.5 # [nm] stepsize y of the periodical box
dz = 0.5 # [nm] stepsize z of the periodical box
sizex = 10# [nm] size of the periodic box in x
sizey = 10# [nm] size of the periodic box in y
sizez = 10# [nm] size of the periodic box in z
e0 = 8.854e-12 # [F/m] vacuum permittivity
er = 1 #[-] relative epermittivity
qe = 1.602e-19 # [C] electron charge
lengthunit = 1e-9 #[-] scaling factor from meter to nm since distance are given in nm

def CoulombPotential(r, q):
    '''return electric potential (V) as function of distance (r) and charge (q)
    distance in m and charge in coulomb'''
    v = q/(4*math.pi*e0*er*r)
    return v

def plot1dpotential(xvalues, yvalues, name = '1dpot'):
    '''plot potential in 1d, xvalues in nm and yvalues in V'''
    fig = plt.figure(1, figsize=((4,3)))
    ax = fig.add_subplot(111)
    ax.plot(xvalues, yvalues)
    ax.set_ylabel('Potential (V)')
    ax.set_xlabel("Distance (nm)")
    fig.savefig(f'{name}.png', bbox_inches = 'tight', pad_inches = 0.3)
    plt.close(fig)
    
def plot2dheatmap(pot, name = '2dpot'):
    '''plot potential in 2d heatmap'''
    fig = plt.figure(1, figsize=((4,3)))
    ax = fig.add_subplot(111)
    heat = ax.imshow(pot, cmap='viridis')
    bar = plt.colorbar(heat)
    bar.set_label("Potential(V)")
    ax.set_ylabel('Nx')
    ax.set_xlabel("Ny")
    fig.savefig(f'{name}.png', bbox_inches = 'tight', pad_inches = 0.3)
    plt.show()
    plt.close(fig)

class charge():
    '''class for holding charges and its xyz coordinate. The xyz coordinates 
    are in system coordinates so q is float and xyz positive ints'''
    def __init__(self, q, x, y, z):
        if not (isinstance(x, int) and isinstance(y, int) and isinstance(z, int)):
            sys.exit("coordinates for charge are not ints")
        self.q = q
        self.x = x
        self.y = y
        self.z = z

class grid():
    '''object to store all the potential values'''
    def __init__(self):
        self.Nx = int(sizex/dx)
        self.Ny = int(sizey/dy)
        self.Nz = int(sizez/dz)
        self.Nsize = self.Nx*self.Ny*self.Nz
        self.pot = np.zeros((self.Nx, self.Ny, self.Nz))
        self.charges = []
        
    def add_charge(self, charge):
        '''adds one charge the system and calculates the potential'''
        self.charges.append(charge)
        chargepot = self.calcpot(charge)
        self.pot += chargepot

    def calcpot(self, charge):
        '''calculation of the potential of one charge xyz are in units of the 
        grid'''
        pot_charge = np.zeros((self.Nx, self.Ny, self.Nz))
        for x in range(self.Nx):
            distx = (x-charge.x)*dx*lengthunit #[m] distance from charge in x
            for y in range(self.Ny):
                disty = (y-charge.y)*dy*lengthunit #[m] distance from charge in y
                for z in range(self.Nz):
                    distz = (z-charge.z)*dz*lengthunit #[m] distance from charge in z
                    if x == charge.x and y == charge.y and z == charge.z:
                        '''avoid singularity'''
                        continue
                    dist = math.sqrt(distx**2+disty**2+distz**2)
                    pot_charge[x,y,z] = CoulombPotential(dist, charge.q)
        return pot_charge

testcharge = charge(qe,3,2,0)
mincharge = charge(-qe,3,3,0)
test_grid = grid()
test_grid.add_charge(testcharge)
test_grid.add_charge(mincharge)
heatmap = test_grid.pot[:,:,0]
plot2dheatmap(heatmap)
