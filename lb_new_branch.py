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
maxIter = 50000 # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 301, 201 # Number of lattice nodes.
ly = ny-1         # Height of the domain in lattice units.
cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
uLB     = 0.04                  # Velocity in lattice units.
nulb    = uLB*r/Re;             # Viscoscity in lattice units.
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
velocity = 0.08             #inital velocity
cs = 1/3                # sound veocity adapted to lattice units     

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

def obstacle_fun(x, y):
    y1 = y > 50
    y2 = y < 150
    x1 = x < 150
    x2 = x > 200
    x3 = x < 250
    # x4 = x >= 0

    return (y1 & y2 & x1)  | (y1 & y2 & x2 & x3)

obstacle = fromfunction(obstacle_fun, (nx,ny))
openPath = invert(obstacle)

def iniVel():
    vel = zeros((2,nx,ny))
    vel[0,0,1:50] = velocity
    return vel

vel = iniVel()

# Initialization of the populations at equilibrium with the given velocity.
fin = equilibrium(1, vel)

##################### DEFINING CONTROL VARIABLES #####################

systemCheck = True
savefiles = False
velocityEvolution = False
StopInOut = False
visualize = False
saveVisual = True

latticePopulation = []
latticePopulationInsideObstacle = []

rho_u_top = [] # rho*u at x = 50
rho_u_bot = [] # rho*u at x = 150¨
rho_u_left = []
rho_u_right = [] # rho*u at x = 250

bilanIn = []
bilanOut = []

velocityEvo = []
plots = 20
plotTime = []

if velocityEvolution:
    new_dir_velocity = "VelocityProfile_" + str(maxIter) + "_it"
    if not os.path.exists("./"+new_dir_velocity):
        os.mkdir(new_dir_velocity)
        print("Made new velocitiy directory : " + new_dir_velocity)

flow = []

if saveVisual :
    new_dir_visual = "Flow_" + str(maxIter) + "_it"
    if not os.path.exists("./"+new_dir_visual):
        os.mkdir(new_dir_visual)
        print("Made new flow directory : " + new_dir_visual)

###### Main time loop ##########################################################
for time in range(maxIter):

    if StopInOut:
        if time<=(maxIter//2):
        # low left wall: outflow condition.
            fin[col3,0,150:ny-1] = fin[col3,1,150:ny-1] 
    else:
        fin[col3,0,150:ny-1] = fin[col3,1,150:ny-1] 

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Top left wall: inflow condition.
    if StopInOut:
        if time<=(maxIter//2):
            u[:,0,1:50] = vel[:,0,1:50]
            rho[0,1:50] = 1/(1-u[0,0,1:50]) * ( sum(fin[col2,0,1:50], axis=0) +
                                        2*sum(fin[col3,0,1:50], axis=0) )
    else : 
        u[:,0,1:50] = vel[:,0,1:50]
        rho[0,1:50] = 1/(1-u[0,0,1:50]) * ( sum(fin[col2,0,1:50], axis=0) +
                                    2*sum(fin[col3,0,1:50], axis=0) )
        
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # bilan entrée-sortie
    bilanIn.append(sum(fin[:,0,1:50]))
    bilanOut.append(sum(fin[:,0,150:ny-1]))

    # Collision step.
    fout = fin - omega * (fin - feq)

    # Bounce-back condition for obstacle.
    for i in range(9):
        # top
        fout[i,:,0] = fin[8-i,:,0]
        # bottom
        fout[i,:,ny-1] = fin[8-i,:,ny-1]
        # right
        fout[i,nx-1,:] = fin[8-i,nx-1,:]
        # obstacle
        fout[i, obstacle] = fin[8-i, obstacle]

    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(
                            roll(fout[i,:,:], v[i,0], axis=0),
                            v[i,1], axis=1 )
 
    # Visualization of the velocity.
    if (time%10==0) and visualize:
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        # plt.title("iteration : %d/%d" % (time, maxIter))
        # plt.savefig("vel.{0:03d}.png".format(time//100))
        plt.pause(.01)
        plt.cla()
        print("iteration : " + str(time) + "/" + str(maxIter), end="\r")
    else : 
        print("iteration : " + str(time) + "/" + str(maxIter), end="\r")
    
    # population control : 
    latticePopulation.append(sum(fin[:,openPath]))
    latticePopulationInsideObstacle.append(sum(fin[:,obstacle]))

    rho_u_top.append(sum(fin[:,75,1:50]*u[0,75,1:50:]))
    rho_u_bot.append(sum(fin[:,75,150:ny-1]*u[0,75,150:ny-1]))
    rho_u_left.append(sum(fin[:,150:200,100]*u[1,150:200,100]))
    rho_u_right.append(sum(fin[:,250:nx-1,100]*u[1,250:nx-1,100]))

    if(time%plots==0):
        if velocityEvo:
            velocityEvo.append(u[0,nx//2,:])
            plotTime.append(time)
        if saveVisual:
            plt.clf()
            plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
            plt.savefig("./" + new_dir_visual + "/fluid_{0:05d}.png".format(time//plots))


#################### SYSTEM CHECKING ###################### 
# VELOCITY
# if velocityEvolution:
#     x = arange(0,len(velocityEvo[0]),1)
#     plotnb = arange(0,len(velocityEvo),1)
#     for u_evo,i in zip(velocityEvo,plotnb):
#         print("making graphs ...", end='\r')
#         plt.clf()
#         plt.plot(x,u_evo)
#         plt.axis([-10,len(x)+10,-0.005,0.05])
#         plt.title("Velocity profile Evolution at " + str(plots*(i+1)) + " it")
#         plt.xlabel("Lattice width coordinates")
#         plt.ylabel("Velocity")
#         plt.savefig("./" + new_dir + "/profile_" + str(i) + ".png")
    
#     print("Velocity graphs done.")

if systemCheck :

    # velocity profile after system stabilisation at coordinates x = nx/2

    # final velocities in different branches
    ux_top = u[0,75,1:51]
    ux_bot = u[0,75,150:ny-1]
    uy_left = u[1,150:200,100]
    uy_right = u[1,250:nx-1,100]

    # theorical variables lattice dependant
    deltaRho = abs(mean(rho[0,:]) - mean(rho[nx-1,:]))
    deltaP = deltaRho/cs**2
    print("Rho left : ", mean(rho[0,:]))
    print("Rho right : ", mean(rho[nx-1,:]))

    R = len(ux_top)//2
    umax_top = u[0,75,25]
    umax_bot = u[0,75,175]
    umax_left = u[1,175,100]
    umax_right = u[1,275,100]
    r = abs(arange(-R,R,1))

    # expected and theorical velocities
    expectedU_top = [umax_top*(1-(i/R)**2) for i in r]
    expectedU_bot = [umax_bot*(1-(i/R)**2) for i in r]
    expectedU_left = [umax_left*(1-(i/R)**2) for i in r]
    expectedU_right = [umax_right*(1-(i/R)**2) for i in r]
    # uformula = [deltaP*(R**2-i**2)/(4*nulb*nx) for i in r]

    # VELOCITY PROFILES
    plt.cla()
    x_axes = arange(0,len(ux_top),1)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.suptitle('Velocity profiles in branches')

    ax1.plot(x_axes,ux_top, label="Real profile")
    ax1.plot(x_axes,expectedU_top, label="Expected profile")
    ax1.legend()
    ax1.set_title("Velocity Profile top tube")
    ax1.set_xlabel("Tube section width")
    ax1.set_ylabel("Velocity")

    ax2.plot(x_axes,ux_bot, label="Real profile")
    ax2.plot(x_axes,expectedU_bot, label="Expected profile")
    ax2.legend()
    ax2.set_title("Velocity Profile bottom tube")
    ax2.set_xlabel("Tube section width")
    ax2.set_ylabel("Velocity")

    ax3.plot(x_axes,uy_left, label="Real profile")
    ax1.plot(x_axes,expectedU_left, label="Expected profile")
    ax3.legend()
    ax3.set_title("Velocity Profile left tube")
    ax3.set_xlabel("Tube section width")
    ax3.set_ylabel("Velocity")

    ax4.plot(x_axes,uy_right, label="Real profile")
    ax4.plot(x_axes,expectedU_right, label="Expected profile")
    ax4.legend()
    ax4.set_title("Velocity Profile right tube")
    ax4.set_xlabel("Tube section width")
    ax4.set_ylabel("Velocity")

    name = "./Monitoring_Branching/" + "Velocity_Profiles_" + str(maxIter)
    if StopInOut: name += "_stopInOutAtHalf"
    if savefiles: plt.savefig(name, bbox_inches='tight')

    plt.show()


    # POPULATION

    # total population
    plt.figure()
    plt.plot(arange(0,len(latticePopulation),1),latticePopulation, label="Population in open path")
    plt.plot(arange(0,len(latticePopulationInsideObstacle),1),latticePopulationInsideObstacle,label="Population in obstacle")
    plt.legend()
    plt.title("Total population over the lattice")
    plt.xlabel("Iterations")
    plt.ylabel("Population")
    name = "./Monitoring_Branching/" + "pop_sum_" + str(maxIter)
    if StopInOut: 
        name += "_stopInOutAtHalf"
        plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
    if savefiles: plt.savefig(name, bbox_inches='tight')
    plt.show()


    # MOMENTUM

    plt.figure()
    plt.plot(arange(0,len(rho_u_top),1),rho_u_top, label="top, x = 75")
    plt.plot(arange(0,len(rho_u_bot),1),rho_u_bot, label="bot, x = 75")
    plt.plot(arange(0,len(rho_u_left),1),rho_u_left, label="left, y = 100")
    plt.plot(arange(0,len(rho_u_right),1),rho_u_right, label="right, y = 100")
    plt.title("Momentum density across lattice")
    plt.xlabel("Iterations")
    plt.ylabel("Momentum [kg*s/m^2]")
    plt.legend(title = "coordinates & tube branch")
    name = "./Monitoring_Branching/" + "rho_u_" + str(maxIter)
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
    name = "./Monitoring_Branching/" + "in_out_" + str(maxIter)
    if StopInOut: 
        name += "_stopInOutAtHalf"
        plt.axvline(x=maxIter//2, color ="r", linestyle = 'dashed')
    if savefiles: plt.savefig(name, bbox_inches='tight')
    plt.show()


####################### COMMENTS & QUESTIONS #################################

# how to pass coordinates range as function parameters ? ex : f(x):print(array[x]) with x = 10:40 ?