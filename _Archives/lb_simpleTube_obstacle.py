#!/usr/bin/python3
# Copyright (C) 2015 Universite de Geneve, Switzerland
# E-mail contact: jonas.latt@unige.ch
#
# 2D flow around a cylinder
#

from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import progressbar
import imageio

###### Flow definition #########################################################
maxIter = 20000  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 300, 100 # Number of lattice nodes.
ly = ny-1         # Height of the domain in lattice units.
cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
uLB     = 0.04                  # Velocity in lattice units.
nulb    = uLB*r/Re;             # Viscoscity in lattice units.
omega = 1 / (3*nulb+0.5);       # Relaxation parameter.

###### Lattice Constants #######################################################
v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
t = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

col1 = array([0, 1, 2])
col2 = array([3, 4, 5])
col3 = array([6, 7, 8])

###### Saving results #######################################################
giffing = False
ProgressVisualisation = True
plots = 500

if ProgressVisualisation:
    widgets = [' [',
            progressbar.Timer(format= 'elapsed time: %(elapsed)s'),
            '] ',
            progressbar.Bar('*'),' (',
            progressbar.ETA(), ') ',
            ]
    bar = progressbar.ProgressBar(max_value=maxIter,widgets=widgets).start()

######################## TESTING
    
visualisation = True
velocityTest = False
populationTest = False

###### Function Definitions ####################################################
def macroscopic(fin):
    rho = sum(fin, axis=0)
    u = zeros((2, nx, ny)) # quantité de mouvement
    for i in range(9):
        u[0,:,:] += v[i,0] * fin[i,:,:]
        u[1,:,:] += v[i,1] * fin[i,:,:]
    u /= rho #qté demouvement /densité = vitesse
    return rho, u

def equilibrium(rho, u):              # Equilibrium distribution function.
    usqr = 3/2 * (u[0]**2 + u[1]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,:,:] + v[i,1]*u[1,:,:])
        feq[i,:,:] = rho*t[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq


###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
# Creation of a mask with 1/0 values, defining the shape of the obstacle.

def obstacle_fun(x, y):
    y1 = y >= 35
    y2 = y <= 65
    x1 = x >= 135
    x2 = x <= 165

    return (y1 & y2 & x1 & x2)

obstacle = fromfunction(obstacle_fun, (nx,ny))

# openPath = invert(obstacle)

# Initial velocity profile: almost zero, with a slight perturbation to trigger
# the instability.
def inivel(d, x, y):
    return (1-d) * uLB * (1 + 1e-4*sin(y/ly*2*pi)) 

vel = fromfunction(inivel, (2,nx,ny))

vel2 = zeros((2,nx,ny))

# Initialization of the populations at equilibrium with the given velocity.
fin = equilibrium(1, vel2)
fout = equilibrium(1, vel2)

##################### TESTING

velTop = []
velMidLeft = []
velMidRight = []
velBot = []
velObstacle = []

popSumOut = []
popSumIn = []
popSumEq = []

###### Main time loop ##########################################################
for time in range(maxIter):

    ####### Visualize progress 

    # print("Time step : ", time,"/",maxIter, end="\r", flush=True)
    if ProgressVisualisation: bar.update(time)

    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(
                            roll(fout[i,:,:], v[i,0], axis=0),
                            v[i,1], axis=1 )

    # right wall: outflow condition. -> sortie gradient nul
    fin[col3,nx-1,:] = fin[col3,nx-2,:] 

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Left wall: inflow condition.  formule zhou-he condition de vitesse (velocity boundary conditions)
    u[:,0,:] = vel[:,0,:]
    rho[0,:] = 1/(1-u[0,0,:]) * ( sum(fin[col2,0,:], axis=0) +
                                  2*sum(fin[col3,0,:], axis=0))
    
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # Collision step.
    # fout[:,openPath] = fin[:,openPath] - omega * (fin[:,openPath] - feq[:,openPath])
    fout = fin - omega * (fin - feq)


 
    
    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(
                            roll(fout[i,:,:], v[i,0], axis=0),
                            v[i,1], axis=1 )
        
       # # Bounce-back condition for obstacle.
    for i in range(9):
        fout[i,obstacle] = fin[8-i,obstacle]
    #     # top
    #     fout[i,:,0] = fin[8-i,:,0]
    #     # bottom
    #     fout[i,:,ny-1] = fin[8-i,:,ny-1]
    #     # right
    #     fout[i,nx-1,1:ny-1] = fin[8-i,nx-1,1:ny-1]
    #     # left
    #     # fout[i,0,:] = fout[8-i,0,:]
    
    ##################### TESTING
    
    if velocityTest and time%100 == 0:
        velTop.append(abs(sum(u[:,75,1:50])))
        velMidLeft.append(abs(sum(u[:,151:201,100])))
        velMidRight.append(abs(sum(u[:,251:301,100])))
        velBot.append(abs(sum(u[:,75,151:201])))
        # velObstacle.append(abs(sum(u[:,obstacle])))
    
    if populationTest : 
        out = 0
        out = sum(fout[:,1:300,1:50])
        out += sum(fout[:,151:200,50:151])
        out += sum(fout[:,251:300,50:151])
        out += sum(fout[:,1:300,151:300])
        popSumOut.append(out)

        sumin = 0
        sumin = sum(fin[:,1:300,1:50])
        sumin += sum(fin[:,151:200,50:151])
        sumin += sum(fin[:,251:300,50:151])
        sumin += sum(fin[:,1:300,151:300])
        popSumIn.append(sumin)

        eq = 0
        eq = sum(feq[:,1:300,1:50])
        eq += sum(feq[:,151:200,50:151])
        eq += sum(feq[:,251:300,50:151])
        eq += sum(feq[:,1:300,151:300])
        popSumEq.append(eq)

    
    ##################################


    if visualisation: 
        if (time%10==0):
            plt.clf()
            plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
            plt.pause(.01)
            plt.cla()


############################ TESTING

# displaying velocities

if velocityTest:
    x = arange(0,len(velTop),1)

    # plt.plot(x,velObstacle)
    # plt.title('Velocity inside obstacle')
    # plt.show()

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    fig.suptitle('Velocities of the system')
    ax1.plot(x,velTop)
    ax1.set_title('Top cylinder')
    ax1.set_ylim([-0.5,4])
    ax2.plot(x,velMidLeft)
    ax2.set_title('Mid cylinder left')
    ax2.set_ylim([-0.5,4])
    ax3.plot(x,velMidRight)
    ax3.set_title('Mid cylinder right')
    ax3.set_ylim([-0.5,4])
    ax4.plot(x,velBot)
    ax4.set_title('Bot cylinder')
    ax4.set_ylim([-0.5,4])
    plt.show()

if populationTest:
    x = arange(0,len(popSumOut),1)
    plt.plot(x,popSumOut)
    plt.title('Population Sum Out')
    plt.show()

    plt.plot(x,popSumIn)
    plt.title('Population Sum In')
    plt.show()

    plt.plot(x,popSumEq)
    plt.title('Population Sum Eq')
    plt.show()


