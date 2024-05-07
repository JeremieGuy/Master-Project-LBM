#!/usr/bin/python3
# Copyright (C) 2015 Universite de Geneve, Switzerland
# E-mail contact: jonas.latt@unige.ch
#
# 2D flow around a cylinder
#

from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm

###### Flow definition #########################################################
maxIter = 20000  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 300, 100 # Number of lattice nodes.
ly = ny-1         # Height of the domain in lattice units.
cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
uLB     = 0.04                  # Velocity in lattice units.
nulb    = uLB*r/Re;             # Viscoscity in lattice units.
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
velocity = 0.04             #inital velocity

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

latticePopulation = []

populationLeft = [] # rho*u at x = 50
populationCenter = [] # rho*u at x = 150
populationRight = [] # rho*u at x = 250

bilan = []

###### Main time loop ##########################################################
for time in range(maxIter):

    if time<=(maxIter//2):
    # Right wall: outflow condition.
        fin[col3,-1,:] = fin[col3,-2,:] 

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    if time<=(maxIter//2):
    # Left wall: inflow condition.
        u[:,0,1:ny-2] = vel[:,0,1:ny-2]
        rho[0,1:ny-2] = 1/(1-u[0,0,1:ny-2]) * ( sum(fin[col2,0,1:ny-2], axis=0) +
                                    2*sum(fin[col3,0,1:ny-2], axis=0) )
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # bilan entrée-sortie
    bilan.append(sum(fin[:,0,:]-fin[:,nx-1,:]))

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
    if (time%10==0):
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        plt.title("iteration : %d/%d" % (time, maxIter))
        # plt.savefig("vel.{0:03d}.png".format(time//100))
        plt.pause(.01)
        plt.cla()
    
    # population control : 
    latticePopulation.append(sum(fin))

    populationLeft.append(sum(fin[:,50,:]*u[0,50,:]))
    populationCenter.append(sum(fin[:,150,:]*u[0,150,:]))
    populationRight.append(sum(fin[:,250,:]*u[0,250,:]))


#################### SYSTEM CHECKING ###################### 
# VELOCITY

# velocity profile after system stabilisation at coordinates x = nx/2

ux = u[0,nx//2,:]

deltaRho = abs(mean(rho[0,:]) - mean(rho[nx-1,:]))

R = ny//2
umax = u[0,nx//2,R]
r = abs(arange(-ny//2,ny//2,1))
expectedU = [umax*(1-(i/R)**2) for i in r]
uformula = [deltaRho*(R**2-i**2)/(4*nulb*nx) for i in r]

mse_expectedU = mean(((ux - expectedU)**2))
mse_physicU = mean((ux - uformula)**2)

# plt.close()
# fig, (ax1, ax2) = plt.subplots(1, 2)
# fig.suptitle("MSE : " + str(mse))
# ax1.plot(arange(0,ny,1),ux)
# ax1.set_title("Velocity profile at x=150")

# ax2.plot(arange(0,ny,1),expectedU)
# ax2.set_title("expected velocity profile")

plt.close()
plt.plot(arange(0,ny,1),ux, label="velocity profile at x=150")
plt.plot(arange(0,ny,1),expectedU, label = "expected velocity profile")
plt.plot(arange(0,ny,1),uformula,label="expected u with phyisc units",)
plt.title("MSE expectedU = " + str(mse_expectedU) + ", MSE physic U = " + str(mse_physicU))
plt.legend()
plt.show()

# POPULATION

# total population
plt.figure()
plt.plot(arange(0,len(latticePopulation),1),latticePopulation)
plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
plt.title("Sum of total population of the lattice")
plt.show()

# population accross lattice
# fig,(popl,popc,popr) = plt.subplots(1,3)
# fig.suptitle("rho*u accross lattice")

# popl.plot(arange(0,len(populationLeft),1),populationLeft)
# popl.set_title("x = 50")

# popc.plot(arange(0,len(populationCenter),1),populationCenter)
# popc.set_title("x = 150")

# popr.plot(arange(0,len(populationRight),1),populationCenter)
# popr.set_title("x = 250")

plt.figure()
plt.plot(arange(0,len(populationLeft),1),populationLeft, label="x = 50")
plt.plot(arange(0,len(populationCenter),1),populationCenter, label="x = 150")
plt.plot(arange(0,len(populationRight),1),populationRight, label = "x = 250")
plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
plt.title("rho*u accross lattice")
plt.legend(title = "x coordinates")
plt.show()

plt.figure()
plt.plot(arange(0,len(bilan),1),bilan)
plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
plt.title("Bilan entrées - sorties")
plt.show()


####################### COMMENTS & QUESTIONS #################################

# how to pass coordinates range as function parameters ? ex : f(x):print(array[x]) with x = 10:40 ?