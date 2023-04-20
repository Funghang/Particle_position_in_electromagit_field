# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 03:14:57 2021

@author: DELL
"""

import numpy as np
from scipy.integrate import odeint

#magnetic field
B = np.array((0,0,1)) 
#Electric field
E = np.array((0,0,0))   

def cyclo_radius(q, m, v0):
   
    B_sq = np.matmul(B,B) 
    v0_par = (np.matmul(v0 , B)) * B / B_sq
    v0_perp = v0 - v0_par
    v0_perp_abs = np.sqrt( np.matmul(v0_perp , v0_perp))
    return m * v0_perp_abs / (np.abs(q) * np.sqrt(B_sq))  #rho = m.v_perp / (|q|B) 

def lorentz_force(X, t, q_over_m):
        
        v = X[3:]   #X defines the particle's position and velocity at time, X=[x,y,z,vx,vy,vz] 
        drdt = v
        dvdt = q_over_m * (E + np.cross(v, B))  #F = ma = (q/m)[E + v×B]
        return np.hstack((drdt, dvdt))

def calc_trajectory(q, m, r0=None, v0 = np.array((1,0,1))):

    if r0 is None:
        rho = cyclo_radius(q, m, v0) 
        vp = np.array((v0[1],-v0[0],0))
        r0 = -np.sign(q) * vp * rho
    tf = 50 #Final time
    N = 10 * tf #number of time steps
    t = np.linspace(0, tf, N) #time grid
    X0 = np.hstack((r0, v0)) #Initial positon and velocity components
    X = odeint(lorentz_force, X0, t, args=(q/m,)) #numerical integration of the equation of motion
    return X


me, qe = 1, -1
v0 = np.array((0,0,0))
X = calc_trajectory(qe, me, r0=(0,0,0), v0=(0.1,0,0))
X2 = calc_trajectory(-2*qe, me, r0=(0,0,0), v0=(0.1,0,0))
X3 = calc_trajectory(-qe, me, r0=(0,0,0), v0=(0.1,0,0))
print (X)