# Author : Guy Jérémie
# inspired by Jonas Latt

from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import time

###### Flow definition #########################################################
maxIter = 2000  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 201, 101 # Number of lattice nodes.
# ly = ny-1         # Height of the domain in lattice units.
# cx, cy, r = nx//4, ny//2, ny//2 # Coordinates of the cylinder.
# rayon = ny//2           # rayon of tube section
R = ny//2
uLB     = 0.04                 # Velocity in lattice units.
nulb    = uLB*ny/Re;             # Viscoscity in lattice units. velocity*characteristic length (= H)/Raynolds
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

def monitorVelocity():
    ux = u[0,nx//2,:]

    # theorical variables lattice dependant
    deltaRho = abs(mean(rho[49,:]) - mean(rho[nx-50,:]))
    deltaP = deltaRho*cs2
    # print("Rho left : ", mean(rho[0,:]))
    # print("Rho right : ", mean(rho[nx-1,:]))

    # R = ny//2
    umax = u[0,nx//2,R]
    rho_center = mean(rho[nx//2,:])

    r = abs(arange((-ny//2)+1,(ny//2)+1,1))
    y_values = arange(0,ny,1)
    
    L = (nx-50)-50
    H = ny-1

    # expected and theorical velocities
    expectedU = [umax*(1-(i/R)**2) for i in r]
    uformulaTube = [deltaP*(R**2-i**2)/(4*nulb*rho_center*L) for i in r]
    uformulaPlane = [deltaP*(y*(H-y))/(2*nulb*rho_center*L) for y in y_values]

    valueFile = open(new_dir_monitoring + "/VeloctiyProfileValues_" + str(execTime) + ".txt", 'w')

    txt = "\nTheoretical values : \n"
    txt +=  "deltaRho : " + str(deltaRho) + "\n"
    txt += "deltaP : " + str(deltaP) + "\n"
    txt += "R : " + str(R) + "\n"
    txt += "viscosity(nulb) : " + str(nulb) + "\n"
    txt += "Rho at the center : " + str(rho_center) + "\n"
    txt += "L : " + str(L) + "\n"
    txt += "Max theoretical velocity : " + str(deltaP*(R**2)/(4*nulb*rho_center*L)) + "\n"

    valueFile.write(txt)
    valueFile.close()

    # VELOCITY PROFILES

    plt.close()
    plt.plot(y_values, ux, label="Real profile at x=150")
    plt.plot(y_values, expectedU, label = "Expected profile")
    plt.plot(y_values, uformulaTube,label="Theoretical profile Tube",)
    plt.plot(y_values, uformulaPlane, label = "Theoretical profile Plane")
    plt.title("Velocity profiles")
    plt.xlabel("Lattice width coordinates")
    plt.ylabel("Velocity")
    plt.legend()
    name = new_dir_monitoring + "/" + "Velocity_Profiles_" + str(execTime)
    # if StopInOut: name += "_stopInOutAtHalf"
    if savefiles: plt.savefig(name, bbox_inches='tight')
    if lookAtGraphs : plt.show()
    plt.clf()

def monitorRho():
    #  DENSITY PROFILE  

    rhoRange = arange(49,nx-50,1)
    meanRho = [mean(rho[i,:]) for i in rhoRange]
    plt.figure()
    plt.plot(rhoRange,meanRho)
    plt.title("Density profile accross lattice")
    plt.xlabel("Coordinate")
    plt.ylabel("Density")
    name = new_dir_monitoring + "/" + "rho_profile_" + str(execTime)
    if savefiles: plt.savefig(name, bbox_inches='tight')
    if lookAtGraphs : plt.show()
    plt.clf()

###### Setup: square obstacle and velocity inlet with perturbation ########
# 0 = open, 1 = bounceback, 2 = obstacle

# open path = 0
flags = zeros((nx,ny))

# bounceback = 1
flags[:,0] = 1
flags[:,ny-1] = 1

# obstacle = (80:120, 210:250)
# bounceback around obstalce
# flags[210:251,80] = 1
# flags[210:251,120] = 1
# flags[210,80:121] = 1
# flags[250,80:121] = 1

# inside obstacle = 2 = no update
# flags[211:250,81:120] = 2

# plt.imshow(flags)
# plt.show()

def iniVel():
    vel = zeros((2,nx,ny))
    distance = abs(arange((-ny//2)+1,(ny//2)+1,1))
    velocity_curve = [velocity*(1-(i/R)**2) for i in distance]
    vel[0,0,:] = velocity_curve
    return vel

vel = iniVel()

# Initialization of the populations at equilibrium with the given velocity.
fin = equilibrium(1, vel)

##################### DEFINING CONTROL VARIABLES #####################

systemCheck = True
savefiles = True
visualize = False
saveVisual = False
lookAtGraphs = False

plots = 50
plotTime = []

if savefiles : 
    new_dir_monitoring = "./Monitoring/Monitoring_tube_changed_nulb" + str(maxIter) + "_it"
    if not os.path.exists(new_dir_monitoring):
        os.mkdir(new_dir_monitoring)
        print("Made new monitoring directory : " + new_dir_monitoring)

if saveVisual :
    new_dir_visual = "./Flow/Flow_no_roll_flags_corrected_pop_and_flag_" + str(maxIter) + "_it"
    if not os.path.exists(new_dir_visual):
        os.mkdir(new_dir_visual)
        print("Made new flow directory : " + new_dir_visual)

start_time = time.time()
###### Main time loop ##########################################################
for execTime in range(maxIter):

    # Right wall: outflow condition.
    fin[col3,-1,:] = fin[col3,-2,:]

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Left wall: inflow condition.
    u[:,0,1:ny-2] = vel[:,0,1:ny-2]
    rho[0,1:ny-2] = 1/(1-u[0,0,1:ny-2]) * ( sum(fin[col2,0,1:ny-2], axis=0) +
                                2*sum(fin[col3,0,1:ny-2], axis=0) )
        
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # bilan entrée-sortie
    # bilanIn.append(sum(fin[:,0,:]))
    # bilanOut.append(sum(fin[:,nx-1,:]))

    # Collision step.
    fout = fin - omega * (fin - feq)

    # Bounce-back condition

    # for i in range(9):
    #     # top
    #     fout[i,:,0] = fin[8-i,:,0]
    #     # bottom
    #     fout[i,:,ny-1] = fin[8-i,:,ny-1]

    for x in range(nx):
        for y in range(ny):
            if flags[x,y] == 1:
                for i in range(9):
                    fout[i,x,y] = fin[8-i,x,y]

    # Streaming step.
    for x in range(nx):
        for y in range(ny):
            for i in range(9):
                if flags[x,y] == 0 or flags[x,y]==1:
                    next_x = x + v[i,0]
                    next_y = y + v[i,1]

                    if next_x < 0:
                        next_x = nx-1
                    if next_x >= nx:
                        next_x = 0
                    
                    if next_y < 0:
                        next_y = ny-1
                    if next_y >= ny:
                        next_y = 0
                    
                    fin[i,next_x,next_y] = fout[i,x,y]
  
 
    # Visualization of the velocity.
    # if (execTime%10==0) and visualize:
    #     plt.clf()
    #     plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    #     plt.title("iteration : %d/%d" % (execTime, maxIter))
    #     # plt.savefig("vel.{0:03d}.png".format(time//100))
    #     plt.pause(.01)
    #     plt.cla()
    # else : 
    #     print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")
    
    # population control : 
    # latticePopulation.append(sum(fin))

    # rho_u_left.append(sum(fin[:,50,:]*u[0,50,:]))
    # rho_u_center.append(sum(fin[:,150,:]*u[0,150,:]))
    # rho_u_right.append(sum(fin[:,250,:]*u[0,250,:]))

    # if(execTime%plots==0):
    #     if velocityEvo:
    #         velocityEvo.append(u[0,nx//2,:])
    #         plotTime.append(execTime)
    #     if saveVisual:
    #         plt.clf()
    #         plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    #         plt.savefig(new_dir_visual + "/fluid_{0:05d}.png".format(execTime//plots))

    print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")
    if(execTime%500==0):
        monitorVelocity()
        monitorRho()

end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

#################### SYSTEM CHECKING ###################### 


monitorVelocity()
monitorRho()



####################### COMMENTS & QUESTIONS #################################

# how to pass coordinates range as function parameters ? ex : f(x):print(array[x]) with x = 10:40 ?