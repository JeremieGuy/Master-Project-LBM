import numpy as np
import sys
import matplotlib.pyplot as plt
import math
import cv2
import matplotlib.patches as patches
import matplotlib

nx = 301
ny = 201

def testMeshGrid():
    nx = np.linspace(-5,5,20)
    ny = np.linspace(-5,5,20)
    X,Y = np.meshgrid(nx,ny)

    # print(X)

    plt.plot(X,Y, 'o')
    plt.show()

def testMatrix():
    x = np.linspace(0,300)
    
    y = np.linspace(0,100,100)
    # print(y)
    # plt.plot(x,y, 'o')
    # plt.show()
    print(np.ones((30,10)))

def testFull():
    x = np.full((10,2),False)
    print(x)

def testDivide():
    x = 100
    print(x/4)
    print(int(x/4))
    print(x//4)

def testRange():
    for i in range(0,14):
        print(i)
    
def rolling():
    x =  [  [[ 0, 1,  2,  3], 
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]],
            
            [[ 0, 1,  2,  3], 
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]]]
    
def testSum():
    x = np.ones([10,5,3])
    b = [0,-1,1]
    print(x)
    print(sum(x,2))
    print(sum(x*b,2))

def testLinesapce():
    # x = np.arange(0,100,1)
    # print(x[1:49])
    # print(len(x[1:49]))
    x = abs(np.arange((-ny//2)+1,(ny//2)+1,1))
    print(x)
    print(len(x))

def testVel(nx,ny):
    np.set_printoptions(threshold=sys.maxsize)
    lattice = np.zeros((2,nx,ny))
    lattice[0,0,:] = 0.1
    print(np.size(lattice))
    print(lattice[0,:,:])
    plt.imshow(np.transpose(lattice[0,:,:]))
    plt.show()

def testUmax():
    distance = abs(np.arange(-ny//2,ny//2,1))
    print(distance)
    print(len(distance))
    umax = 0.04
    expectedU = [umax*(1-(i/ny)**2) for i in distance]
    plt.plot(np.arange(0,ny,1),expectedU)
    plt.show()

def testcoord(x):
    distance = abs(np.arange(-ny//2,ny//2,1))
    print(distance[int(x)])


def testflag():
    f = np.zeros((9,10,10))
    flags = np.zeros((10,10))

    for x in range(10):
        for y in range(10):
            for i in range(9):
                if x == y :
                    flags[x,y] = 1
                    f[i,x,y] = 0
    
    print(flags)
    print(f)

def testflagobstacle():
    # 0 = open, 1 = bounceback, 2 = obstacle

    # open path = 0
    flags = np.zeros((nx,ny))

    # bounceback = 1
    flags[:,0] = 1 # top
    flags[:,ny-1] = 1 # bot
    flags[nx-1,:] = 1 # right
    flags[0,0:51] = 1 # behind inlet
    # obstacle1 = (80:120, 210:250)
    # bounceback around obstalce
    flags[0:151,50] = 1
    flags[0:151,150] = 1
    flags[0,50:151] = 1
    flags[150,50:151] = 1

    # inside obstacle = 2 = no update
    flags[1:150,51:150] = 2

    # obstacle 2 = (200:250,50:150)
    # bounceback around obstalce
    flags[200:251,50] = 1
    flags[200:251,150] = 1
    flags[200,50:151] = 1
    flags[250,50:151] = 1

    # inside obstacle = 2 = no update
    flags[201:250,51:150] = 2

    plt.imshow(np.transpose(flags))
    plt.show()

def test1():

    x = np.arange(0,10,1)
    print(x)
    print(x[3:7])

def testwrite():
    f = np.zeros([4,3,4])
    # print(f)
    f[3,2,3] = 0.123456789
    # a = f[3,2,3]
    f = open("test.txt", 'w')
    # x = {f[3,2,3]:4d}
    # a = "{:.3f}".format(f[3,2,3]) + " " + "{:.3f}".format(f[3,2,3])
    # print(a)
    templ = '{0:.2f} {1:.2f}\n'
    f.write(templ.format(f[3,2,3], f[3,2,3]))
    f.close()

def testinivel():
    nx = 201
    ny = 101
    velocity = 0.04
    R = ny//2
    vel = np.zeros((2,nx,ny))
    distance = abs(np.arange((-ny//2)+1,(ny//2)+1,1))
    velocity_curve = [velocity*(1-(i/R)**2) for i in distance]
    vel[0,0,:] = velocity_curve
    plt.imshow(np.transpose(vel[0,:,:]))
    plt.colorbar(location = "left", label = "horizontal velocity [m/s]")
    plt.show()

def testAffichage():
    ny = 5
    print(np.arange(ny-1,-1,-1))

def testflagsmallsystem():
    # cv2.rectangle(img,(384,0),(510,128),(0,255,0),3)
    flags = np.zeros((7,5))
    flags[:,0] = 1
    flags[:,4] = 1
    flags[4,2] = 1
    flags[0,1:4] = 2
    flags[6,:] = 3
    # flags[5,2] = 1
    
    fig, ax = plt.subplots()
    ax.imshow(np.transpose(flags))
    rect = patches.Rectangle((5.5, 1.5), 1, 1, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(rect)

    # plt.imshow()
    plt.show()

def rangetest():
    
    ny = 20
    test = np.arange(0,ny,1)
    

    print(test)
    print(test[0:11])

def rolltest():
    mat = np.zeros([5,5])
    mat[2:5,2:5] = 1
    # mat2 = np.roll(np.roll(mat,-1,axis=0),-1,axis=1)
    mat2 =np.roll(mat,1,axis=0)
    print(mat2)

def notMatTest():
    mat = ([
        [True, False, False],
        [False, True, False],
        [False, False, True]
    ])

    print(np.invert(mat))

def testexepctedprofile():
    tubesSize = 51
    r = abs(np.arange((-tubesSize//2)+1,(tubesSize//2)+1,1))
    R = tubesSize//2
    print("r = ",r)
    print("R = ", R)
    test = 4
    test2 = [test*(1-(i/R)**2) for i in r]
    print(test2)

def testClotVelocities():

    u = np.random.randint(0, 100, size=(2, 20, 10))
    print(u)

    fig, ( ax2, ax3) = plt.subplots(2, 1)
    fig.suptitle('Velocity Profiles test')

    # ax1.imshow(np.sqrt(u[0]**2+u[1]**2).transpose())
    # plt.show()

    # plt.clf()
    clotVel = np.sqrt(u[0,5:15,1:9]**2+u[1,5:15,1:9]**2)
    ax2.imshow(clotVel.transpose())
    

    vel = np.sum(clotVel,axis=0)
    x = np.arange(0,len(vel),1)
    ax3.plot(x,vel)
    plt.show()

# testMatrix()
# testMeshGrid()
# testFull()
# testDivide()
# testRange()
# rolling()
# testSum()
# testLinesapce()
# testVel(30,10)
# testUmax()
# testcoord("20:30")s
# testflag()
# testflagobstacle()        
# test1()
# testwrite()
# testinivel()
# testAffichage()
testflagsmallsystem()
# rangetest()
# rolltest()
# notMatTest()
# testexepctedprofile()

# testClotVelocities()