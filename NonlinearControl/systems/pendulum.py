# -*- coding: utf-8 -*-
"""
Based on code from Mohammad Asif Zaman (https://github.com/zaman13/Double-Pendulum-Motion-Animation/blob/master/Python%20Code/Double_Pendulum_v1.py)

Double pendulum motion animation using FuncAnimation()
"""

from __future__ import print_function   
from scipy.integrate import odeint 

import time
import math
import numpy as np
import pylab as py


#import matplotlib.pyplot as plt

from matplotlib import animation, rc
from IPython.display import HTML
from matplotlib import pyplot as plt
x0 = [-np.pi/2.2, 0]    # initial conditions. 
# x[0] = angle 
# x[1] = angular velocity 

tfinal = 25.0       # Final time. Simulation time = 0 to tfinal.
Nt = 751

def get_params():
    L = 1.4              # length of pendulum 1 (in meter)
    g = 9.8               # gravitatioanl acceleration constant (m/s^2)
    return L, g

# Differential equations describing the system
def pendulum_system(x, t, controller):    
    L, g = get_params()
    dx = np.zeros(2)
    dx[0] = x[1]   # d(theta 1)
    u = controller(x)
    dx[1] = -(g/L)*np.sin(x[0]) + u
    return dx

def simulate_pendulum(x0, tfinal, Nt, controller):
    t = np.linspace(0, tfinal, Nt)
    sol = odeint(pendulum_system, x0, t, args=(controller,))
    sol0 = sol[:,0]     # theta_1 
    sol1 = sol[:,1]     # omega 1
    return sol0, sol1, t[2]-t[1]

def make_animation(dt, theta_vec):
    
    def get_mass_pos(theta):
        L, g = get_params()
        x = L*np.sin(theta)
        y = -L*np.cos(theta)
        return x, y
    
    # initialization function: plot the background of each frame
    def init():
        line1.set_data([], [])
        line4.set_data([], [])
        line5.set_data([], [])
        time_string.set_text('')
        return  line4, line5, line1, time_string
    
    # animation function.  This is called sequentially
    def animate(i):
        # Motion trail sizes. Defined in terms of indices. Length will vary with the time step, dt. E.g. 5 indices will span a lower distance if the time step is reduced.
        trail1 = 6              # length of motion trail of weight 1 
        line1.set_data(posx[i:max(1,i-trail1):-1], posy[i:max(1,i-trail1):-1])   # marker + line of first weight
        line4.set_data([posx[i], 0], [posy[i],0])                # line connecting origin to weight 1
        line5.set_data([0, 0], [0, 0])
        time_string.set_text(time_template % (i*dt))
        return  line4,line5,line1, time_string
    
    posx, posy = get_mass_pos(theta_vec)
    
    fig = plt.figure()
    L,g = get_params()
    ax = plt.axes(xlim=(-L-0.5, L+0.5), ylim=(-2.5, 1.5))
    line1, = ax.plot([], [], 'o-',color = '#d2eeff',markersize = 12, markerfacecolor = '#0077BE',lw=2, markevery=10000, markeredgecolor = 'k')   # line for Earth
    line4, = ax.plot([], [], color='k', linestyle='-', linewidth=2)
    line5, = ax.plot([], [], 'o', color='k', markersize = 10)
    time_template = 'Time = %.1f s'
    time_string = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    
    ax.get_xaxis().set_ticks([])    # enable this to hide x axis ticks
    ax.get_yaxis().set_ticks([])    # enable this to hide y axis ticks
    ax.set_aspect('equal', adjustable='box')
    
    # MAIN ANIMATION PLOTTING
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Nt, interval=1000*(dt)*0.8, blit=True)
    plt.close(fig)
    return anim

# Example controller function
def controller(x):
    return 0  # No control input

# thetavec, thetadotvec, dt = simulate_pendulum(x0, tfinal, Nt, controller)
# anim = make_animation(dt, thetavec)

# Comment out the following lines if you do not want to save the animation to file
# anim.save('pendulum.gif', fps=1.0/(dt), writer = 'imagemagick')
