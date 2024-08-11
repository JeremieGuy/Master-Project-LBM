# Author : Guy Jérémie
# inspired by Jonas Latt

from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import time

###### Flow definition #########################################################
maxIter = 5  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 5,4 # Number of lattice nodes.
# ly = ny-1         # Height of the domain in lattice units.
# cx, cy, r = nx//4, ny//2, ny//2 # Coordinates of the cylinder.
# rayon = ny//2           # rayon of tube section
R = ny//2
uLB     = 0.04                 # Velocity in lattice units.
nulb    = uLB*R/Re;             # Viscoscity in lattice units. velocity*rayon/Raynolds
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


###### Setup: square obstacle and velocity inlet with perturbation ########
# 0 = open, 1 = bounceback, 2 = obstacle

# open path = 0
flags = zeros((nx,ny))

# bounceback = 1
flags[:,0] = 1
flags[:,ny-1] = 1

directions = [0, 1, 2, 3, 4, 5, 6, 7, 8]
coordinates = [2,2]

def drawoutput():
    for y in y_range:
        line1 = ""
        line2 = ""
        line3 = ""

        for x in range(nx):
            line1 += str(fin[2,x,y]) + " " + str(fin[5,x,y]) + " " + str(fin[8,x,y]) + "    "
            line2 += str(fin[1,x,y]) + " " + str(fin[4,x,y]) + " " + str(fin[7,x,y]) + "    "
            line3 += str(fin[0,x,y]) + " " + str(fin[3,x,y]) + " " + str(fin[6,x,y]) + "    "
        
        line1 += "|   "
        line2 += "|   "
        line3 += "|   "

        for x in range(nx):
            line1 += str(fout[6,x,y]) + " " + str(fout[3,x,y]) + " " + str(fout[0,x,y]) + "    "
            line2 += str(fout[7,x,y]) + " " + str(fout[4,x,y]) + " " + str(fout[1,x,y]) + "    "
            line3 += str(fout[8,x,y]) + " " + str(fout[5,x,y]) + " " + str(fout[2,x,y]) + "    "
        
        populationFile.write("\n" + line1 + "\n" + line2 + "\n" + line3 + "\n")




##################### DEFINING CONTROL VARIABLES #####################


new_dir_monitoring = "./Monitoring/Monitoring_simple_system_inverted_affichage_fin_" + str(coordinates) + "_" + str(maxIter) + "_it"
if not os.path.exists(new_dir_monitoring):
    os.mkdir(new_dir_monitoring)
    print("Made new monitoring directory : " + new_dir_monitoring)


 
for i_dir in directions : 




    populationFile = open(new_dir_monitoring + "/population_bounceback_dir_fin=" + str(i_dir) + ".txt", 'w')

    populationFile.write("Flags : \n\n" + str(flags) +"\n\n")

    populationFile.write("Directions for respectively fin & fout : \n")
    y_range = arange(ny-1,-1,-1)
    for y in y_range:
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

    dotline = "\n"
    for i in range(121):
        dotline+="-"
    dotline += "\n"


    fin = zeros([9,nx,ny])
    fout = zeros([9,nx,ny])

    fin[i_dir,coordinates[0],coordinates[1]] = 1

    ###### Main time loop ##########################################################
    for execTime in range(maxIter):
        
        it = "\nIteration : " + str(execTime) + "\n"
        populationFile.write(dotline + it)

        # Collision step.
        # fout = fin - omega * (fin - feq)
        
        # fout = fin
        
        # Bounce-back condition
        bb_changes = "\nBounceback occured on : \n"
        bb_changes_bool = False



        for x in range(nx):
            for y in range(ny):
                if flags[x,y] == 1:
                    for i in range(9):
                        fout[i,x,y] = fin[8-i,x,y]
                        
                        
                        if fin[8-i,x,y] != 0:
                            bb_changes_bool = True
                            bb_changes += "- fin[" + str(8-i) + "," + str(x)+ "," + str(y)+ "]"
                            bb_changes += " to fout[" + str(i) + "," + str(x)+ "," + str(y)+ "]\n"
            
        if not bb_changes_bool:
            fout = fin
                            

        fin = zeros([9,nx,ny])

        col = "\nAfter Collision, Before Streaming\n"
        populationFile.write(col)

        drawoutput()

        if bb_changes_bool:
            populationFile.write(bb_changes)
        else : 
            populationFile.write("\nNo bounceback\n")

        # Streaming step.

        stream_changes = "\nStreaming moved value on : \n"
        stream_changes_bool = False

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
                        
                        if fout[i,x,y] != 0:
                            # print("stream")
                            stream_changes_bool = True
                            stream_changes += "- fout[" + str(i) + "," + str(x) + "," + str(y)+ "]"
                            stream_changes += " to fin[" + str(i) + ","+str(next_x)+ "," + str(next_y)+ "]\n"
        
        fout = zeros([9,nx,ny])
        
        stream = "\nAfter Streaming, Before Collision\n"

        populationFile.write(stream)

        drawoutput()

        if stream_changes_bool:
            populationFile.write(stream_changes)
        else : 
            populationFile.write("\nValue did not stream\n")

        
        # print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")

# print("finished simulation - " + str(execTime+1) + "/" + str(maxIter))
  
