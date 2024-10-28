# Author : Guy Jérémie
# inspired by Jonas Latt

from numpy import *
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import os
import time

########################### Flow definition #############################################

maxIter = 30000  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 301, 151 # Number of lattice nodes.
# ly = ny-1         # Height of the domain in lattice units.
# cx, cy, r = nx//4, ny//2, ny//2 # Coordinates of the cylinder.
# rayon = ny//2           # rayon of tube section
R = ny//2
uLB     = 0.04                 # Velocity in lattice units.
nulb    = uLB*ny/Re;             # Viscoscity in lattice units. velocity*characteristic length (= H)/Raynolds
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
velocity = 0.04             #inital velocity
cs2 = 1/3                # sound veocity adapted to lattice units     
gamma = 1 # solid fraction


########################## Lattice Constants ###########################################

v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
t = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

col1 = array([0, 1, 2])
col2 = array([3, 4, 5])
col3 = array([6, 7, 8])

#################### Main Function Definitions ######################################

# macroscopic variable computation
def macroscopic(fin):
    rho = zeros((nx,ny))
    rho[openPath] = sum(fin[:,openPath], axis=0)
    
    u = zeros((2, nx, ny))
    for i in range(9):
        u[0,openPath] += v[i,0] * fin[i,openPath]
        u[1,openPath] += v[i,1] * fin[i,openPath]
    u[:,openPath] /= rho[openPath]
    return rho, u

# Equilibrium distribution function (rho = array)
def equilibrium(rho, u):          
    usqr = 3/2 * (u[0,openPath]**2 + u[1,openPath]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,openPath] + v[i,1]*u[1,openPath])
        feq[i,openPath] = rho[openPath]*t[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq

# set intial equilibrium state with a fixed int density rho
def equilibriumInitial(rho, u):           
    usqr = 3/2 * (u[0,openPath]**2 + u[1,openPath]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,openPath] + v[i,1]*u[1,openPath])
        feq[i,openPath] = rho*t[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq

# initial velocity for a simple pipe flowing from left to right
def iniVel():
    vel = zeros((2,nx,ny))
    distance = abs(arange((-ny//2)+1,(ny//2)+1,1))
    velocity_curve = [velocity*(1-(i/R)**2) for i in distance]
    vel[0,0,:] = velocity_curve
    return vel

# initlai velocity for a localised inlet on the left border
def iniVel2(inlety1,inlety2):
    vel = zeros((2,nx,ny))
    inletSize = inlety2+1-inlety1
    distance = abs(arange((-inletSize//2)+1,(inletSize//2)+1,1))
    velocity_curve = [velocity*(1-(i/R)**2) for i in distance]
    vel[0,0,inlety1:(inlety2+1)] = velocity_curve
    return vel

# setting relevant BB nodes to 0 on the 31x23 system
def setBBNodeToZero():
    # top border
    fin[:,:,ny-1] = 0
    # bottom border
    fin[:,:,0] = 0
    # right border
    fin[:,nx-1,:] = 0

    # obstacle
    fin[:,0:249,52:99] = 0

     # top border
    fout[:,:,ny-1] = 0
    # bottom border
    fout[:,:,0] = 0
    # right border
    fout[:,nx-1,:] = 0

    # obstacle
    fout[:,0:249,52:99] = 0

############################ Monitoring functions #################################

def plotSystem():

    # norm = plt.Normalize(vmin=0,vmax=0.7)
    flagsname = ["open path", "bounceback", "inlet", "outlet", "blood clot"]
    plt.figure(figsize=(7.9,4))
    plt.title("Flags")
    values = unique(flags_plot.ravel())
    
    im = plt.imshow(flags_plot.transpose())
    colors = [ im.cmap(im.norm(value)) for value in values]
    patches = [ mpatches.Patch(color=colors[i], label=flagsname[i] ) for i in range(len(values)) ]
    plt.title("Flags")
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    plt.plot([10,10],[100,100], color="red", linewidth=3)
    plt.savefig(new_dir_monitoring + "/system.png",bbox_inches='tight')
    plt.show()
    plt.close()

def plotPopulation():
    plt.figure()
    # y_formatter = mp.ticker.ScalarFormatter(useOffset=False)
    # plt.gca().yaxis.set_major_formatter(y_formatter)
    plt.plot(arange(0,len(latticePopulation),1),latticePopulation)
    plt.title("Total population over the lattice")
    plt.xlabel("Iterations")
    plt.ylabel("Population")
    name = new_dir_monitoring + "/" + "pop_sum_" + str(maxIter)
    plt.savefig(name, bbox_inches='tight')
    plt.show()

def plotVelocityProfiles():

    #### Final velocities
    # Top  
    uTop = abs(u[0,150,1:52])
    umaxTop = max(uTop)
    # Left cylinder
    # uMid = abs(u[1,249:300,75])
    # umaxMid = max(uMid)
    # Bottom
    uBot = abs(u[0,150,99:150])
    umaxBot = max(uBot)

    # velocity Plotting variables
    tubesSize = len(uTop)
    r = abs(arange((-tubesSize//2)+1,(tubesSize//2)+1,1))
    R = tubesSize//2
    
    x = arange(0,tubesSize,1)

    # expected velocities
    expectedUTop = [umaxTop*(1-(i/R)**2) for i in r]
    # expectedUMid = [umaxMid*(1-(i/R)**2) for i in r]
    expectedUBot = [umaxBot*(1-(i/R)**2) for i in r]

    # plot
    plt.clf()

    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig, (ax1, ax3) = plt.subplots(1, 2)

    fig.suptitle('Velocity Profiles')

    ax1.plot(x,uTop, label = "Real Profile")
    ax1.plot(x,expectedUTop, label = "Expected Profile")
    ax1.set_title('Top cylinder')
    ax1.set_xlabel("Tube width coordinates")
    ax1.set_ylabel("Velocity")
    ax1.legend()
    ax1.set_ylim([-0.01,umaxTop+0.01])

    # ax2.plot(x,uMid, label = "Real Profile")
    # ax2.plot(x,expectedUMid, label = "Expected Profile")
    # ax2.set_title('Middle cylinder')
    # ax2.set_xlabel("Tube width coordinates")
    # ax2.set_ylabel("Velocity")
    # ax2.legend()
    # ax2.set_ylim([-0.01,umaxTop+0.01])

    ax3.plot(x,uBot, label = "Real Profile")
    ax3.plot(x,expectedUBot, label = "Expected Profile")
    ax3.set_title('Bottom cylinder')
    ax3.set_xlabel("Tube width coordinates")
    ax3.set_ylabel("Velocity")
    ax3.legend()
    ax3.set_ylim([-0.01,umaxTop+0.01])
    
    name = new_dir_monitoring + "/" + "velocity_profiles_" + str(maxIter)
    plt.savefig(name, bbox_inches='tight')
    plt.show()

def saveLastFrame():
    plt.clf()
    plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    plt.title("Velocity at iteration = " +  str(maxIter))
    plt.savefig(new_dir_monitoring + "/last_frame_velocity.png")

##################### DEFINING CONTROL VARIABLES #####################

savefiles = True

plots = 10
plotTime = []
latticePopulation = []

if savefiles : 
    new_dir_monitoring = "./Monitoring/u_clot_v1_gamma=" +str(gamma) + "_it=" + str(maxIter)
    if not os.path.exists(new_dir_monitoring):
        os.mkdir(new_dir_monitoring)
        print("Made new monitoring directory : " + new_dir_monitoring)

########################################### PLOTTING VARIABLES & FLAGS ############################

#  Bounceback nodes mask
bounceback = full((nx,ny),False)
bounceback[:,0] = True
bounceback[:,ny-1] = True
bounceback[nx-1,:] = True
bounceback[0:249,52:99] = True

# open path flags
openPath = invert(bounceback)

# clot
clot = full((nx,ny),False)
clot[249:300,60:91] = True

#### draw system
# open path = 0
flags_plot = zeros((nx,ny))

# bounceback = 1
flags_plot[:,0] = 1
flags_plot[:,ny-1] = 1
flags_plot[nx-1,:] = 1
flags_plot[0:249,52:99] = 1

# inlet = 2 
flags_plot[0,1:52] = 2

# outlet = 3
flags_plot[0,99:150] = 3

# clot = 4
flags_plot[249:300,60:91] = 4

plotSystem()

##################################### Control Variables ############################################

# lattice total population
latticePopulation = []

################################### System Initliaization ##########################################

# initial velocity
vel = iniVel2(1,51)

# # Initialization of the populations at equilibrium with the given density & velocity.
fin = equilibriumInitial(1, vel)
fout = equilibriumInitial(1, vel)

#### set BB nodes PDFs to 0 
setBBNodeToZero()

#################################### Main time loop #################################################
start_time = time.time()

for execTime in range(maxIter):

    # lower left wall: outflow condition.
    fin[col3,0,99:150] = fin[col3,1,99:150]

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Left wall: inflow condition.
    u[:,0,1:52] = vel[:,0,1:52]
    rho[0,1:52] = 1/(1-u[0,0,1:52]) * ( sum(fin[col2,0,1:52], axis=0) + 2*sum(fin[col3,0,1:52], axis=0) )

    # Compute equilibrium.
    feq = equilibrium(rho, u)
    fin[[0,1,2],0,1:52] = feq[[0,1,2],0,1:52] + fin[[8,7,6],0,1:52] - feq[[8,7,6],0,1:52]

    # Collision step for open path
    fout[:,openPath] = fin[:,openPath] - omega * (fin[:,openPath] - feq[:,openPath])

    # Partial collision for clot
    for i in range(9):
        fout[i,clot] = (1 - gamma)*fout[i,clot] + gamma*fin[8-i,clot]

    # Bounce-back condition
    for i in range(9):
        fout[i, bounceback] = fin[8-i, bounceback]
    
    # streaming step
    # i = 0
    fin[0,:,:] = roll(roll(fout[0,:,:],1,axis=0),1,axis=1)
    # i = 1
    fin[1,:,:] = roll(fout[1,:,:],1,axis=0)
    # i = 2
    fin[2,:,:] = roll(roll(fout[2,:,:],1,axis=0),-1,axis=1)
    # i = 3
    fin[3,:,:] = roll(fout[3,:,:],1,axis=1)
    # i = 4
    fin[4,:,:] = fout[4,:,:]
    # i = 5
    fin[5,:,:] = roll(fout[5,:,:],-1,axis=1)
    # i = 6
    fin[6,:,:] = roll(roll(fout[6,:,:],-1,axis=0),1,axis=1)
    # i = 7
    fin[7,:,:] = roll(fout[7,:,:],-1,axis=0)
    # i = 8
    fin[8,:,:] = roll(roll(fout[8,:,:],-1,axis=0),-1,axis=1)


    if (execTime%10==0):
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        # plt.savefig("vel.{0:03d}.png".format(time//100))

        plt.pause(.01)
        plt.cla()
    

    print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")

    latticePopulation.append(sum(fin))


end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

#################### SYSTEM CHECKING ###################### 

saveLastFrame()

plotPopulation()

plotVelocityProfiles()

####################### COMMENTS & QUESTIONS #################################
# When changing dimensions, Change : flags, Mask, Bounceback nodes, inlet, outlet, velocity profiles