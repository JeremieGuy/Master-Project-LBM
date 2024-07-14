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

###### Flow definition #########################################################
maxIter = 15000  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 301, 201 # Number of lattice nodes.
ly = ny-1         # Height of the domain in lattice units.
cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
uLB     = 0.04                  # Velocity in lattice units.
nulb    = uLB*r/Re;             # Viscoscity in lattice units.
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
velocity = 0.04             #inital velocity
cs2 = 1/3                # sound veocity adapted to lattice units     

###### Lattice Constants #######################################################
v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
t = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

col1 = array([0, 1, 2])
col2 = array([3, 4, 5])
col3 = array([6, 7, 8])

###### Function Definitions ####################################################
def macroscopic(fin):
    rho = sum(fin, axis=0)
    u = zeros((2, nx, ny))
    for i in range(9):
        u[0,:,:] += v[i,0] * fin[i,:,:]
        u[1,:,:] += v[i,1] * fin[i,:,:]
    u /= rho
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

# def obstacle_fun(x, y):
#     return (x-cx)**2+(y-cy)**2<r**2

# obstacle = fromfunction(obstacle_fun, (nx,ny))

# Initial velocity profile: almost zero, with a slight perturbation to trigger
# the instability.
# def inivel(d, x, y):
#     return (1-d) * uLB * (1 + 1e-4*sin(y/ly*2*pi))

# vel = fromfunction(inivel, (2,nx,ny))

def iniVel():
    vel = zeros((2,nx,ny))
    vel[0,0,:] = velocity
    return vel

vel = iniVel()

# Initialization of the populations at equilibrium with the given velocity.
fin = equilibrium(1, vel)

##################### DEFINING CONTROL VARIABLES #####################

systemCheck = True
velocityEvolution = False
StopInOut = True
visualize = False

latticePopulation = []

rho_u_left = [] # rho*u at x = 50
rho_u_center = [] # rho*u at x = 150
rho_u_right = [] # rho*u at x = 250

bilanIn = []
bilanOut = []

velocityEvo = []
plots = 100
plotTime = []

if velocityEvolution:
    new_dir = "VelocityProfile_" + str(maxIter) + "_it"
    if not os.path.exists("./"+new_dir):
        os.mkdir(new_dir)
        print("Made new directory : " + new_dir)

###### Main time loop ##########################################################
for time in range(maxIter):

    if StopInOut:
        if time<=(maxIter//2):
        # Right wall: outflow condition.
            fin[col3,-1,:] = fin[col3,-2,:] 
    else : 
        fin[col3,-1,:] = fin[col3,-2,:]

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    if StopInOut:
        if time<=(maxIter//2):
        # Left wall: inflow condition.
            u[:,0,1:ny-2] = vel[:,0,1:ny-2]
            rho[0,1:ny-2] = 1/(1-u[0,0,1:ny-2]) * ( sum(fin[col2,0,1:ny-2], axis=0) +
                                        2*sum(fin[col3,0,1:ny-2], axis=0) )
    else : 
        u[:,0,1:ny-2] = vel[:,0,1:ny-2]
        rho[0,1:ny-2] = 1/(1-u[0,0,1:ny-2]) * ( sum(fin[col2,0,1:ny-2], axis=0) +
                                    2*sum(fin[col3,0,1:ny-2], axis=0) )
        
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # bilan entrée-sortie
    bilanIn.append(sum(fin[:,0,:]))
    bilanOut.append(sum(fin[:,nx-1,:]))

    # Collision step.
    fout = fin - omega * (fin - feq)

    # Bounce-back condition for obstacle.
    for i in range(9):
        # top
        fout[i,:,0] = fin[8-i,:,0]
        # bottom
        fout[i,:,ny-1] = fin[8-i,:,ny-1]

    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(
                            roll(fout[i,:,:], v[i,0], axis=0),
                            v[i,1], axis=1 )
 
    # Visualization of the velocity.
    if (time%10==0) and visualize:
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        plt.title("iteration : %d/%d" % (time, maxIter))
        # plt.savefig("vel.{0:03d}.png".format(time//100))
        plt.pause(.01)
        plt.cla()
    else : 
        print("iteration : " + str(time) + "/" + str(maxIter), end="\r")
    
    # population control : 
    latticePopulation.append(sum(fin))

    rho_u_left.append(sum(fin[:,50,:]*u[0,50,:]))
    rho_u_center.append(sum(fin[:,150,:]*u[0,150,:]))
    rho_u_right.append(sum(fin[:,250,:]*u[0,250,:]))

    if(time%plots==0):
        velocityEvo.append(u[0,nx//2,:])
        plotTime.append(time)


#################### SYSTEM CHECKING ###################### 
# VELOCITY
if velocityEvolution:
    x = arange(0,len(velocityEvo[0]),1)
    plotnb = arange(0,len(velocityEvo),1)
    for u_evo,i in zip(velocityEvo,plotnb):
        print("making graphs ...", end='\r')
        plt.clf()
        plt.plot(x,u_evo)
        plt.axis([-10,len(x)+10,-0.005,0.05])
        plt.title("Velocity profile Evolution at " + str(plots*(i+1)) + " it")
        plt.xlabel("Lattice width coordinates")
        plt.ylabel("Velocity")
        plt.savefig("./" + new_dir + "/profile_" + str(i) + ".png")
    
    print("Velocity graphs done.")

if systemCheck :
    savefiles = True

    # velocity profile after system stabilisation at coordinates x = nx/2

    # final colucity at the center
    ux = u[0,nx//2,:]

    # theorical variables lattice dependant
    deltaRho = abs(mean(rho[0,:]) - mean(rho[nx-1,:]))
    deltaP = deltaRho*cs2
    print("Rho left : ", mean(rho[0,:]))
    print("Rho right : ", mean(rho[nx-1,:]))

    R = ny//2
    umax = u[0,nx//2,R]
    r = abs(arange((-ny//2)+1,(ny//2)+1,1))

    # expected and theorical velocities
    expectedU = [umax*(1-(i/R)**2) for i in r]
    uformula = [deltaP*(R**2-i**2)/(4*nulb*nx) for i in r]

    mse_expectedU = mean(((ux - expectedU)**2))
    mse_physicU = mean((ux - uformula)**2)

    # VELOCITY PROFILES

    plt.close()
    plt.plot(arange(0,ny,1), ux, label="Real profile at x=150")
    plt.plot(arange(0,ny,1), expectedU, label = "Expected profile")
    plt.plot(arange(0,ny,1), uformula,label="Theoretical profile",)
    plt.title("Velocity profiles")
    plt.xlabel("Lattice width coordinates")
    plt.ylabel("Velocity")
    plt.legend()
    name = "./Monitoring_Tube/" + "Velocity_Profiles_" + str(maxIter)
    if StopInOut: name += "_stopInOutAtHalf"
    if savefiles: plt.savefig(name, bbox_inches='tight')
    plt.show()

    # POPULATION

    # total population
    plt.figure()
    plt.plot(arange(0,len(latticePopulation),1),latticePopulation)
    plt.title("Total population over the lattice")
    plt.xlabel("Iterations")
    plt.ylabel("Population")
    name = "./Monitoring_Tube/" + "pop_sum_" + str(maxIter)
    if StopInOut: 
        name += "_stopInOutAtHalf"
        plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
    if savefiles: plt.savefig(name, bbox_inches='tight')
    plt.show()


    # MOMENTUM

    plt.figure()
    plt.plot(arange(0,len(rho_u_left),1),rho_u_left, label="x = 50")
    plt.plot(arange(0,len(rho_u_center),1),rho_u_center, label="x = 150")
    plt.plot(arange(0,len(rho_u_right),1),rho_u_right, label = "x = 250")
    plt.title("Momentum density across lattice")
    plt.xlabel("Iterations")
    plt.ylabel("Momentum [kg*s/m^2]")
    plt.legend(title = "x coordinates")
    name = "./Monitoring_Tube/" + "rho_u_" + str(maxIter)
    if StopInOut: 
        name += "_stopInOutAtHalf"
        plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
    if savefiles: plt.savefig(name, bbox_inches='tight')
    plt.show()

    # INPUT & OUTPUT

    plt.figure()
    plt.plot(arange(0,len(bilanIn),1),bilanIn, label="Inflow")
    plt.plot(arange(0,len(bilanOut),1),bilanOut, label="Outflow")
    plt.title("Inflow & Outflow assessment")
    plt.xlabel("Iterations")
    plt.ylabel("Population")
    plt.legend()
    name = "./Monitoring_Tube/" + "in_out_" + str(maxIter)
    if StopInOut: 
        name += "_stopInOutAtHalf"
        plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
    if savefiles: plt.savefig(name, bbox_inches='tight')
    plt.show()


####################### COMMENTS & QUESTIONS #################################

# how to pass coordinates range as function parameters ? ex : f(x):print(array[x]) with x = 10:40 ?