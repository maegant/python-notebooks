
from __future__ import print_function   
from scipy.integrate import odeint 

import time
import math
import numpy as np
import pylab as py


#import matplotlib.pyplot as plt

from matplotlib import animation
import matplotlib.patches as patches

from IPython.display import HTML
from matplotlib import pyplot as plt
x0 = [0, 0, 5*np.pi/6, 0]    # initial conditions. 
# x[0] = cart position
# x[1] = pendulum angle 
# x[2] = cart velocity
# x[3] = angular velocity 

tfinal = 25.0       # Final time. Simulation time = 0 to tfinal.
Nt = 751

def get_params():
    L = 1.4              # length of pendulum 1 (in meter)
    g = 9.8               # gravitatioanl acceleration constant (m/s^2)
    mc = 5.0            # mass of cart
    mp = 1.0            # mass of pendulum
    rh = 0.5           # Rectangle Height
    rw = 0.5           # Rectangle Width
    return L, g, mc, mp, rh, rw

def get_model(x):
    L, g, mc, mp, rh, rw = get_params()
    M = np.array([[mc + mp, mp * L * np.cos(x[2])],
                  [mp * L * np.cos(x[2]), mp * L**2]])
    H = np.array([[- mp * L * x[3]**2 * np.sin(x[2])], [mp * g * L * np.sin(x[2])]])
    B = np.array([[1],[0]])
    return M, H, B

# Differential equations describing the system
def cartpole_system(x, t, controller):    
        
    dx = np.zeros(4)
    dx[0] = x[1]   # d(pos)
    dx[2] = x[3]   # d(theta)
    
    u = controller(x) # force on cart only applied to cart position
    M, H, B = get_model(x)
    qddot = np.linalg.inv(M) @ (-H + (np.dot(B,u)))
    dx[1] = qddot[0,0]
    dx[3] = qddot[1,0]
    return dx

def simulate_cartpole(x0, tfinal, Nt, controller):
    t = np.linspace(0, tfinal, Nt)
    sol = odeint(cartpole_system, x0, t, args=(controller,))
    sol0 = sol[:,0]     # pos 
    sol1 = sol[:,1]     # vel
    sol2 = sol[:,2]     # theta
    sol3 = sol[:,3]     # omega
    return sol0, sol1, sol2, sol3, t[2]-t[1]

def make_animation(dt, cartpos_vec,  theta_vec):
    
    def get_mass_pos(theta):
        L, g, mc, mp, rh, rw = get_params()
        x = L*np.sin(theta)
        y = -L*np.cos(theta)
        return x, y
        
    # initialization function: plot the background of each frame
    def init():
        line1.set_data([], [])
        line4.set_data([], [])
        line5.set_data([], [])
        time_string.set_text('')
        return  line4, line5, line1, cart, time_string
    
    # animation function.  This is called sequentially
    def animate(i):
        # Motion trail sizes. Defined in terms of indices. Length will vary with the time step, dt. E.g. 5 indices will span a lower distance if the time step is reduced.
        trail1 = 6              # length of motion trail of weight 1 
        line1.set_data(posx[i:max(1,i-trail1):-1]+cartpos_vec[i:max(1,i-trail1):-1], posy[i:max(1,i-trail1):-1])   # marker + line of pendulum mass
        line4.set_data([posx[i]+cartpos_vec[i], cartpos_vec[i]], [posy[i],0])                # line connecting pivot point to pendulum mass
        line5.set_data([cartpos_vec[i], cartpos_vec[i]], [0, 0])            # pivot point of pendulum
        # cart.set_data([cartpos_vec[i] - 0.1, cartpos_vec[i] + 0.1, cartpos_vec[i] + 0.1, cartpos_vec[i] - 0.1, cartpos_vec[i] - 0.1],
                    #   [-0.1, -0.1, 0.1, 0.1, -0.1]) 
        cart.set_xy((cartpos_vec[i]-rw/2, -rh/2))
        time_string.set_text(time_template % (i*dt))
        return  line4,line5,line1,cart, time_string
    
    posx, posy = get_mass_pos(theta_vec)
    
    fig, ax = plt.subplots()
    L, g, mc, mp, rh, rw = get_params()
    cart = patches.Rectangle((-rw/2, -rh/2), rw, rh, facecolor = 'blue')
    ax.add_patch(cart)
    
    line1, = ax.plot([], [], 'o-',color = '#d2eeff',markersize = 12, markerfacecolor = '#0077BE',lw=2, markevery=10000, markeredgecolor = 'k')   # line for Earth
    line4, = ax.plot([], [], color='k', linestyle='-', linewidth=2)
    line5, = ax.plot([], [], 'o', color='k', markersize = 10)
    time_template = 'Time = %.1f s'
    time_string = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    
    ax.set_xlim(min(cartpos_vec) - L, max(cartpos_vec) + L)
    ax.set_ylim(-L-0.5, L+0.5)
    ax.get_xaxis().set_ticks([])    # enable this to hide x axis ticks
    ax.get_yaxis().set_ticks([])    # enable this to hide y axis ticks
    ax.set_aspect('equal', adjustable='box') 
    
    # MAIN ANIMATION PLOTTING
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Nt, interval=1000*(dt)*0.8, blit=True)
    plt.close(fig)
    # plt.show()
    return anim

# Example controller function
def controller(x):
    return 0  # No control input

# Comment out the following lines if you do not want to save the animation to file
#anim.save('double_pendulum_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
# anim.save('double_pendulum_animation.gif', fps=1.0/(dt), writer = 'imagemagick')

# Run the simulation
posvec, velvec, thetavec, thetadotvec, dt = simulate_cartpole(x0, tfinal, Nt, controller)
anim = make_animation(dt, posvec, thetavec)