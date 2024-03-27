import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

 # plot
plottingFrames = 10

# animation
fig, ax = plt.subplots()
ims = []

# distance function : used to compute cylinder diameter
def distance(x1,y1,x2,y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def main():
    ''' Variables '''
    # lattice dimensions = nb cells
    Nx = 300
    Ny = 200 # for 2D
    # Timesteps (simulation duration = iterations)
    Nt = 2000
    # time scale = kinematic velocity = relaxation time
    tau = .53

    # time and space steps
    deltaT = 1
    deltaX = 1

    '''Lattice Speeds & Weights''' # each direction has a set speed of 1 to go to a node in any direction
    # Nb of velocities (= neighbours)
    # NL = 3 # 1D
    NL = 9 # 2D

    # discreete velocities of each neihbouring node
    # 1D - nodes order [center, left, right]
    # cxs = np.array([0,-1,1])  
    # 2D - nodes order : [center, top left, top center, top right, medium center, low right, low center, left center, medium left]
    cxs = np.array([0,-1,0,1,1, 1, 0,-1,-1])
    cys = np.array([0, 1,1,1,0,-1,-1,-1, 0])
    
    # weights
    # weights = np.array([2/9,5/9,2/9]) # 1D
    weights = np.array([4/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36,1/9]) # 2D

    '''Initial conditions'''
    # "phase space " = distribution function of 3 dimensions : 1 or2 geometrical dimension(s), 1 velocity dimension (= mesoscopic velocity)
    # use of np.ones because more stable than np.zeroes (source : dude trust me)
    # f = np.ones([Nx//deltaX,NL]) # 1D
    # f = np.ones([Nx//deltaX,NL]) + .01 * np.random.randn(Nx, NL) # 1D + randomness (adds 'swirly' effect to fluid)
    # f = np.ones([Ny//deltaX,Nx//deltaX,NL]) # 2D
    # f = np.ones((Ny//deltaX,Nx//deltaX,NL)) + .01 * np.random.randn(Ny, Nx, NL) # 2D
    f = np.full((Ny//deltaX,Nx//deltaX,NL),0.0001)  # 2D

    # Initial directionnal velocity (to the right in this case)
    # f[:,2] = 2.3 # 1D
    # f[:,:, 4] = 2.3 # 2D
    f[:,:, 4] = 0 # 2D

    '''Optionnal : obstacle (2D)'''
    # define cylinder as array (same size as phase space) with initial False value 
    # cylinder = np.full((Ny, Nx), False)
    # # if point (x,y) is in cylinder, set value to True
    # for y in range(0,Ny):
    #     for x in range(0,Nx):
    #         if (distance(Nx//4,Ny//2,x,y)<10 or distance(3*Nx//4,Ny//2,x,y)<10 or distance(Nx//2,Ny//2,x,y)<14) : # set center at half of y height and 1/4 of x distance ('//' = integer division)
    #             cylinder[y][x] = True
    

    obstacle = np.full((Ny, Nx), False)

    for y in range(0,Ny):
            for x in range(0,Nx):
                    if y>50 and y<150 and x<150:
                            obstacle[y][x] = True
                    if y>50 and y<150 and x>200 and x<250:
                            obstacle[y][x] = True
    
    openPath = np.invert(obstacle)
    
    '''Main Loop'''
    time1 = time.time()
    for it in range(Nt//deltaT):
        # print("iteration : ", it)

        ''' Boundary Condition'''
        # Absorbing the waves on the last section of the lattice by setting them as the value right next to it
        # -> absorbs waves that were reflected form the  other side 
        f[150:200, -1, [1, 7, 8]] = f[150:200, -2, [1, 7, 8]]
        # f[:, 0, [[3, 4, 5]]] = f[:, 1, [[3, 4, 5]]]

        # rebound X 
        f[:, -1, :] = f[:, -1, [0,5,6,7,8,1,2,3,4]]
        # f[:, 0, :] = f[:, 0, [0,5,6,7,8,1,2,3,4]]

        # rebound Y
        # f[-1, :, :] = f[-1, :, [0,5,6,7,8,1,2,3,4]]
        # f[0, :, :] = f[0, :, [0,5,6,7,8,1,2,3,4]]

        ''' Streaming step '''
        # for  i, cx in zip(range(NL),cxs): # 1D, zip = (it,cx,cy), total = NL
        #     f[:,i] = np.roll(f[:,i],cx,axis = 1)
        # for [x,y] in openPath:
        for  i, cx, cy in zip(range(NL),cxs,cys): # 2D, zip = (it,cx,cy), total = NL
            f[:,:,i] = np.roll(f[:,:,i], cx, axis = 1)
            f[:,:,i] = np.roll(f[:,:,i], cy, axis = 0)
            # f[x,y,i] = np.roll(f[x,y,i], cx, axis = 1)
            # f[x,y,i] = np.roll(f[x,y,i], cy, axis = 0)

        # set boundaries
        boundaryF = f[obstacle,:] # set the nodes that are in the cylinder ... 
        boundaryF = boundaryF[:,[0,5,6,7,8,1,2,3,4]] # ... and bounce them back = define their velocity as the one in the opposite direction

        ''' Fluid variables '''

        # defining density to compute equilibrium set (1D)
        # rho = np.sum(f, 1)
        # momentum
        # ux = np.sum(f * cxs, 1) / rho

        # defining density to compute equilibrium set (2D)
        rho = np.sum(f, 2)
        # momentum
        ux = np.sum(f * cxs, 2) / rho
        uy = np.sum(f * cys, 2) / rho

        # Boundary step (obstacle - 2D only)
        f[obstacle,:] = boundaryF
        ux[obstacle] = 0
        uy[obstacle] = 0

        ''' collision step '''
        feq = np.ones(f.shape) # initialize equilibrium state
        # compute equilibrium state : 
        # for i , cx, w in zip((range(NL), cxs, weights)): # 1D
        #     feq[:, i] = w * rho * (1 + 3 * (cx*ux) + 9 * (cx*ux)**2 / 2 - 3  * (ux**2) / 2)
        for i, cx,cy,w in zip(range(NL),cxs,cys,weights): # 2D
            feq[:,:,i] = w * rho * (1 + 3 * (cx*ux + cy*uy) + 9 * (cx*ux + cy*uy)**2 / 2 - 3  * (ux**2 + uy**2) / 2)

        # compute next step
        f = f + -(1/tau) * (f - feq)

        ''' Plot '''
        if (it%plottingFrames == 0):
            # plt.imshow(np.sqrt(ux**2)) # 1D

            plt.imshow(np.sqrt(ux**2+uy**2), cmap="gist_heat") # 2D - full fluid
            # compute fluid high swirling zones 
            # dfydx = ux[2:, 1:-1] - ux[0:-2, 1:-1]
            # dfxdy = uy[1:-1, 2:] - uy[1:-1, 0:-2]
            # swirl = dfydx - dfxdy
            # plt.imshow(swirl, cmap="seismic")

            plt.pause(.01)
            plt.cla()

        # ''' Animation '''
        # # im = ax.plot(np.sqrt(ux**2+uy**2)) # F*ucked up but nice
        # im = ax.imshow(np.sqrt(ux**2+uy**2))
        # ims.append(im)
    
    time2 = time.time()

    print("temps : ", time2-time1)



# def animate(i):
#     return ims[i]

main()

# ani = animation.ArtistAnimation(fig, animate, interval=50, blit=True, repeat_delay=1000)
# plt.show()



''' OPEN BF '''