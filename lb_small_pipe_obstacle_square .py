# Author : Guy Jérémie
# inspired by Jonas Latt

from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# from matplotlib import cm
import os
import time

###### Flow definition #########################################################
maxIter = 100  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 7, 5 # Number of lattice nodes.
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

inletCorner = True
outletCorner = False

if inletCorner: inletY = [0,ny]
else : inletY = [1,ny-1]

if outletCorner : outletY = [0,ny]
else : outletY = [1,ny-1]

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

def writeFileMonitoring(x,y):
    lineiter = "iteration : " + str(execTime) + "\n\n"
    line1 = str('%.5f'%(fin[2,x,y])) + " " + str('%.5f'%(fin[5,x,y])) + " " + str('%.5f'%(fin[8,x,y])) 
    line2 = str('%.5f'%(fin[1,x,y])) + " " + str('%.5f'%(fin[4,x,y])) + " " + str('%.5f'%(fin[7,x,y])) 
    line3 = str('%.5f'%(fin[0,x,y])) + " " + str('%.5f'%(fin[3,x,y])) + " " + str('%.5f'%(fin[6,x,y])) 
    line1_2 = " | " + str('%.5f'%(fout[6,x,y])) + " " + str('%.5f'%(fout[3,x,y])) + " " + str('%.5f'%(fout[0,x,y]))
    line2_2 = " | " + str('%.5f'%(fout[7,x,y])) + " " + str('%.5f'%(fout[4,x,y])) + " " + str('%.5f'%(fout[1,x,y])) + " \n"
    line3_2 = " | " + str('%.5f'%(fout[8,x,y])) + " " + str('%.5f'%(fout[5,x,y])) + " " + str('%.5f'%(fout[2,x,y])) + " \n"
    line1_3 = "       " + str('%.5f'%(u[0,x,y])) + "      " + str('%.5f'%(u[1,x,y])) + "\n"
    line4 = "\n---------------------------------------------------------------------------\n"

    monitorFile.write(lineiter+line1+line1_2+line1_3+line2+line2_2+line3+line3_2+line4)

def saveNodeImage(x,y, time):
    tmpOut = [fout[i,x,y] for i in range(9)]
    rtmpOut = reshape(tmpOut,[3,3])

    tmpIn = [fin[i,x,y] for i in range(9)]
    rtmpIn = reshape(tmpIn,[3,3])

    plt.clf()
    fig, axes = plt.subplots(1, 2,figsize=(12, 6))

    fig.suptitle("Fin & Fout for node [" + str(x) + "," + str(y) +"] with flag = 1, it = " + str(time), fontsize=16)

    im1 = axes[0].imshow(rtmpIn, cmap='viridis', norm=norm)
    axes[0].set_title('Fin')
    cbar1 = plt.colorbar(im1, ax=axes[0])
    cbar1.ax.set_position([axes[0].get_position().x1 + 0.01,
                    axes[0].get_position().y0,
                    0.02,
                    axes[0].get_position().height])
    axes[0].axis('off')

    im2 = axes[1].imshow(rtmpOut, cmap='viridis', norm=norm)
    axes[1].set_title('Fout')
    cbar2 = plt.colorbar(im2, ax=axes[1])
    cbar2.ax.set_position([axes[1].get_position().x1 + 0.01,
                    axes[1].get_position().y0,
                    0.02,
                    axes[1].get_position().height])
    axes[1].axis('off')

    plt.savefig(new_dir_node + "/node_{0:05d}.png".format(execTime//plots))
    plt.close()

def drawoutput():
    populationFile.write("Iteration : " + str(execTime) + "\n")
    for y in range(ny):
        line1 = ""
        line2 = ""
        line3 = ""

        for x in range(nx):
            line1 += str('%.5f'%(fin[2,x,y])) + " " + str('%.5f'%(fin[5,x,y])) + " " + str('%.5f'%(fin[8,x,y]))  + "    "
            line2 += str('%.5f'%(fin[1,x,y])) + " " + str('%.5f'%(fin[4,x,y])) + " " + str('%.5f'%(fin[7,x,y])) + "    "
            line3 += str('%.5f'%(fin[0,x,y])) + " " + str('%.5f'%(fin[3,x,y])) + " " + str('%.5f'%(fin[6,x,y])) + "    "
        
        line1 += "|   "
        line2 += "|   "
        line3 += "|   "

        for x in range(nx):
            line1 += str('%.5f'%(fout[6,x,y])) + " " + str('%.5f'%(fout[3,x,y])) + " " + str('%.5f'%(fout[0,x,y])) + "    "
            line2 += str('%.5f'%(fout[7,x,y])) + " " + str('%.5f'%(fout[4,x,y])) + " " + str('%.5f'%(fout[1,x,y])) + "    "
            line3 += str('%.5f'%(fout[8,x,y])) + " " + str('%.5f'%(fout[5,x,y])) + " " + str('%.5f'%(fout[2,x,y])) + "    "
        
        populationFile.write("\n" + line1 + "\n" + line2 + "\n" + line3 + "\n\n")
    
    populationFile.write("Horizontal velocity                                                           | Vertical velocity\n")

    for y in range(ny):
        line = ""
        
        for x in range(nx): line += str('%.5f'%(u[0,x,y])) + "    "

        line += "|   "

        for x in range(nx): line += str('%.5f'%(u[1,x,y])) + "    "
        
        populationFile.write("\n" + line + "\n")


    populationFile.write(dotline + "\n")

###### Setup: square obstacle and velocity inlet with perturbation ########
# 0 = open, 1 = bounceback, 2 = obstacle

# open path = 0
flags = zeros((nx,ny))

# bounceback = 1

flags[:,0] = 1
flags[:,ny-1] = 1

if inletCorner: 
    flags[0,0] = 0
    flags[0,ny-1] = 0

if outletCorner:
    flags[nx-1,0] = 0
    flags[nx-1,ny-1] = 0

# obstacle = (80:120, 210:250)
# bounceback around obstalce
flags[4,2] = 1
# flags[7,2] = 1

# inside obstacle = 2 = no update
# flags[131:160,36:65] = 2

flags_plot = copy(flags)

#inlet 
flags_plot[0,inletY[0]:inletY[1]] = 2

# outlet
flags_plot[-1,outletY[0]:outletY[1]] = 3


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
visualize = True
saveVisual = False
lookAtGraphs = False
monitorNode = True

plots = 10
plotTime = []

if savefiles : 
    new_dir_monitoring = "./Monitoring/Small_tube_square_obstacle_" + str(maxIter) + "_it_corner=" + str(inletCorner)
    if not os.path.exists(new_dir_monitoring):
        os.mkdir(new_dir_monitoring)
        print("Made new monitoring directory : " + new_dir_monitoring)

if saveVisual :
    new_dir_visual = "./Flow/Flow_no_roll_flags_corrected_pop_and_flag_" + str(maxIter) + "_it"
    if not os.path.exists(new_dir_visual):
        os.mkdir(new_dir_visual)
        print("Made new flow directory : " + new_dir_visual)

if monitorNode:
    new_dir_node = new_dir_monitoring + "/FoutEvo_node[130,50]"
    if not os.path.exists(new_dir_node):
        os.mkdir(new_dir_node)
        print("Made new node visualisation directory : " + new_dir_node)



monitorFile = open(new_dir_monitoring + "/monitor_obstacle_node_maxiter=" + str(maxIter) + "_corner=" + str(inletCorner) + ".txt", 'w')
populationFile = open(new_dir_monitoring + "/monitor_full_system_maxiter=" + str(maxIter) + "_corner=" + str(inletCorner) + ".txt", 'w')
nodeX, nodeY = 6, 2


linetitle = "Node : [" +str(nodeX) + "," + str(nodeY) +"], maxiter = " + str(maxIter) + "\n\n"
line = "Fin | Fout    Horizontal Velocity  Vertical velocity\n\n"
line1 = "2" + " " + "5" + " " + "8" + " | " + "6" + " " + "3" + " " + "0" + "       " + "u[0]" + "      " + "u[1]\n"
line2 = "1" + " " + "4" + " " + "7" + " | " + "7" + " " + "4" + " " + "1" + " \n"
line3 = "0" + " " + "3" + " " + "6" + " | " + "8" + " " + "5" + " " + "2" + " \n"
line4 = "\n-----------------------------------------------------------------------------\n"

monitorFile.write(linetitle+line+line1+line2+line3+line4)

dotline = "\n"
for i in range(121):
    dotline+="-"
dotline += "\n"

populationFile.write("MaxIter : " + str(maxIter) + ", system size : [" + str(nx) + "x"+ str(ny) + "]\n\n")
populationFile.write("System : \n\n" + str(transpose(flags_plot)) +"\n\n0 = open path, 1 = BB, 2 = inlet, 3 = outlet\n\n")

populationFile.write("Directions for respectively fin & fout : \n")

for y in range(ny):
    line1 = ""
    line2 = ""
    line3 = ""

    for x in range(nx):
        line1 += "2" + " " + "5" + " " + "8" + "    "
        line2 += "1" + " " + "4" + " " + "7" + "    "
        line3 += "0" + " " + "3" + " " + "6" + "    "
    
    line1 += "|   "
    line2 += "|   "
    line3 += "|   "

    for x in range(nx):
        line1 += "6" + " " + "3" + " " + "0" + "    "
        line2 += "7" + " " + "4" + " " + "1" + "    "
        line3 += "8" + " " + "5" + " " + "2" + "    "
    
    populationFile.write("\n" + line1 + "\n" + line2 + "\n" + line3 + "\n")

populationFile.write("\nHorizontal Velocity                                                      | Vertical Velocity\n\n")

line1 = "u[0,0,0]  u[0,1,0]  u[0,2,0]  u[0,3,0]  u[0,4,0]  u[0,5,0]  u[0,6,0]     |   u[1,0,0]  u[1,1,0]  u[1,2,0]  u[1,3,0]  u[1,4,0]  u[1,5,0]  u[1,6,0]\n\n"
line2 = "u[0,0,1]  u[0,1,1]  u[0,2,1]  u[0,3,1]  u[0,4,1]  u[0,5,1]  u[0,6,1]     |   u[1,0,1]  u[1,1,1]  u[1,2,1]  u[1,3,1]  u[1,4,1]  u[1,5,1]  u[1,6,1]\n\n"
line3 = "u[0,0,2]  u[0,1,2]  u[0,2,2]  u[0,3,2]  u[0,4,2]  u[0,5,2]  u[0,6,2]     |   u[1,0,2]  u[1,1,2]  u[1,2,2]  u[1,3,2]  u[1,4,2]  u[1,5,2]  u[1,6,2]\n\n"
line4 = "u[0,0,3]  u[0,1,3]  u[0,2,3]  u[0,3,3]  u[0,4,3]  u[0,5,3]  u[0,6,3]     |   u[1,0,3]  u[1,1,3]  u[1,2,3]  u[1,3,3]  u[1,4,3]  u[1,5,3]  u[1,6,3]\n\n"
line5 = "u[0,0,4]  u[0,1,4]  u[0,2,4]  u[0,3,4]  u[0,4,4]  u[0,5,4]  u[0,6,4]     |   u[1,0,4]  u[1,1,4]  u[1,2,4]  u[1,3,4]  u[1,4,4]  u[1,5,4]  u[1,6,4]\n\n"

populationFile.write(line1 + line2 + line3 + line4 + line5)

populationFile.write(dotline)

########################################### PLOTTING VARIABLES & FLAGS ############################

norm = plt.Normalize(vmin=0,vmax=0.7)
flagsname = ["open path", "bounceback", "inlet", "outlet"]
plt.figure(figsize=(10,4))
plt.title("Flags")
values = unique(flags_plot.ravel())
im = plt.imshow(flags_plot.transpose())
colors = [ im.cmap(im.norm(value)) for value in values]
patches = [ mpatches.Patch(color=colors[i], label=flagsname[i] ) for i in range(len(values)) ]
plt.title("Flags")
plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
plt.savefig(new_dir_monitoring + "/system.png")
plt.close()

###### Main time loop ##########################################################
start_time = time.time()

for execTime in range(maxIter):

    # Right wall: outflow condition.
    # fin[col3,-1,:] = fin[col3,-2,:]
    fin[col3,-1,outletY[0]:outletY[1]] = fin[col3,-2,outletY[0]:outletY[1]]

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Left wall: inflow condition.
    u[:,0,1:ny-1] = vel[:,0,1:ny-1]
    rho[0,1:ny-1] = 1/(1-u[0,0,1:ny-1]) * ( sum(fin[col2,0,1:ny-1], axis=0) +
                                2*sum(fin[col3,0,1:ny-1], axis=0) )
    # u[:,0,inletY[0]:outletY[1]] = vel[:,0,inletY[0]:outletY[1]]
    # rho[0,inletY[0]:outletY[1]] = 1/(1-u[0,0,inletY[0]:outletY[1]]) * ( sum(fin[col2,0,inletY[0]:outletY[1]], axis=0) +
    #                             2*sum(fin[col3,0,inletY[0]:outletY[1]], axis=0) )
        
    # Compute equilibrium.
    feq = equilibrium(rho, u)
    # fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]
    fin[[0,1,2],0,1:ny-1] = feq[[0,1,2],0,1:ny-1] + fin[[8,7,6],0,1:ny-1] - feq[[8,7,6],0,1:ny-1]
    # fin[[0,1,2],0,inletY[0]:inletY[1]] = feq[[0,1,2],0,inletY[0]:inletY[1]] + fin[[8,7,6],0,inletY[0]:inletY[1]] - feq[[8,7,6],0,inletY[0]:inletY[1]]

    if inletCorner:
        # top corner
        fin[col1,0,0] = fin[col1,0,1]
        fin[5,0,0] = fin[1,0,1]
        fin[8,0,0] = fin[2,0,1]

        # bottom corner
        fin[col1,0,-1] = fin[col1,0,-2]
        fin[3,0,-1] = fin[1,0,-2]
        fin[6,0,-1] = fin[0,0,-2]

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
    
    if (execTime%plots==0):
        # print("write")
        saveNodeImage(nodeX,nodeY, execTime)
    
    writeFileMonitoring(nodeX,nodeY)
    drawoutput()


    # if saveVisual and execTime%plots==0:
    #     plt.clf()
    #     plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    #     plt.savefig(new_dir_visual + "/fluid_{0:05d}.png".format(execTime//plots))

    print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")
    # if(execTime%500==0):
    #     monitorVelocity()
    #     monitorRho()

end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

#################### SYSTEM CHECKING ###################### 

####################### COMMENTS & QUESTIONS #################################
