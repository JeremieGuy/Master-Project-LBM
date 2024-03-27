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
maxIter = 2000  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 300, 200 # Number of lattice nodes.
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
plots = 500

if giffing:
    new_dir = "Iteration_" + str(maxIter) + "v2"

    widgets = [' [',
         progressbar.Timer(format= 'elapsed time: %(elapsed)s'),
         '] ',
           progressbar.Bar('*'),' (',
           progressbar.ETA(), ') ',
          ]
    
    os.mkdir(new_dir)
    bar = progressbar.ProgressBar(max_value=maxIter,widgets=widgets).start()

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
    y1 = y > 50
    y2 = y < 150
    x1 = x < 150
    x2 = x > 200
    x3 = x < 250

    return (y1 & y2 & x1) | (y1 & y2 & x2 & x3)

obstacle = fromfunction(obstacle_fun, (nx,ny))

openPath = invert(obstacle)

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

###### Main time loop ##########################################################
for time in range(maxIter):
    # low left wall: outflow condition. -> sortie gradient nul
    fin[col3,0,151:200] = fin[col3,1,151:200] 

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Left wall: inflow condition.  formule zhou-he condition de vitesse (velocity boundary conditions)
    u[:,1,1:49] = vel[:,1,1:49]
    rho[1,1:49] = 1/(1-u[0,1,1:49]) * ( sum(fin[col2,1,1:49], axis=0) +
                                  2*sum(fin[col3,1,1:49], axis=0))
    
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # Collision step.
    fout[:,openPath] = fin[:,openPath] - omega * (fin[:,openPath] - feq[:,openPath])


    # Bounce-back condition for obstacle.
    for i in range(9):
        fout[i, obstacle] = fin[8-i, obstacle]
        # top
        fout[i,:,0] = fin[8-i,:,0]
        # bottom
        fout[i,:,ny-1] = fin[8-i,:,ny-1]
        # right
        fout[i,nx-1,1:ny-2] = fin[8-i,nx-1,1:ny-2]
        # # left
        # fout[i,0,:] = fout[8-i,0,:]
    
    # Streaming step.
    for i in range(9):
        fin[i,:,:] = roll(
                            roll(fout[i,:,:], v[i,0], axis=0),
                            v[i,1], axis=1 )
    
    ##################### TESTING
    
    if time%100 == 0:
        velTop.append(sum(u[:,75,1:49]))

    ##################################

    # Visualization of the velocity.
    if giffing:
        if (time%plots==0):
            plt.clf()
            plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
            plt.savefig("./" + new_dir + "/fluid.{0:03d}.png".format(time//plots))
        bar.update(time)
    else : 
         if (time%10==0):
            plt.clf()
            plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
            plt.pause(.01)
            plt.cla()

###### Make GIF ##########################################################
if giffing: 
    frames = []
    t = maxIter//plots

    print("\nMaking Gif ...")
    for i in range(t):
        num = "{0:0=3d}".format(i)
        image = imageio.v2.imread(f"./" + new_dir + "/fluid." + num + ".png")
        frames.append(image)

    imageio.mimsave("./" + new_dir + ".gif", frames,duration = 80)


############################ TESTING
x = linspace()

# site 0 et 50 -> bounceback, 1 à 40 inlet


# propagation - condition de bord - collision

# mesure grandeur physique => juste vant la collision

# mesure variable physique -> resultat de propagation = fin => jsute avant collision il faut
# condition de bord partout utilisable pour collision

# bounceback noeud vs lien
# bounceback lien -> modifier opérateur de propagation





# somme Fi
# rho*u mesure en haut, tubes et bas , somme ()

# profil de vitesse -> ligne verticale (entre 0 et 50) vers 50,70 en x => tracer amplitude vitesse en horizontale
# idem tubes verticaux